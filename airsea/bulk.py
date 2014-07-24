from warnings import warn

import numpy as np
import seawater as sw
from functools import wraps
from oceans.sw_extras import visc, tcond
from .atmosphere import viscair, qsat
from .windstress import cdn
from .constant import (Charnock_alpha, R_roughness, gas_const_R, Qsat_coeff,
                       emiss_lw, sigmaSB, eps_air, kappa, CtoK, cp, g)


def to_array(a):
    return np.atleast_1d(np.asanyarray(a))


def ensure_atleast_1d(f):
    """Decorator o make all *args in to atleast_1d arrays."""
    @wraps(f)
    def wrapper(*args, **kw):
        newargs = [to_array(a) for a in args]
        ret = f(*newargs, **kw)
        return ret
    wrapper.__wrapped__ = f
    return wrapper


# Default values used in Fairall et al, 1996

CVB_depth = 600  # Depth of convective boundary layer in atmosphere [m].
beta_conv = 1.25  # Scaling constant for gustiness.
min_gustiness = 0.5  # min. "gustiness".  i.e. Unresolved fluctuations) [m s-1]
# should keep this strictly >0, otherwise bad stuff might happen (divide by
# zero errors).


def LKB(Reu):
    """Computes roughness Reynolds numbers for temperature and humidity
    Liu, Katsaros and Businger (1979), J. Atmos. Sci., 36, 1722-1735.

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    08/28/1998: version 1.1
    08/05/1999: version 1.2
    11/26/2010: Filipe Fernandes, Python translation.
    """

    # Convert input to numpy array.
    Reu = np.asanyarray(Reu)
    Reu = np.atleast_1d(Reu)

    Ret = 0.177 * np.ones_like(Reu)
    Req = 0.292 * np.ones_like(Reu)

    j = np.logical_and(Reu > 0.11, Reu <= 0.825)
    Ret[j] = 1.376 * Reu[j]**0.929
    Req[j] = 1.808 * Reu[j]**0.826

    j = np.logical_and(Reu > 0.825, Reu <= 3)
    Ret[j] = 1.026 / Reu[j]**0.599
    Req[j] = 1.393 / Reu[j]**0.528

    j = np.logical_and(Reu > 3, Reu <= 10)
    Ret[j] = 1.625 / Reu[j]**1.018
    Req[j] = 1.956 / Reu[j]**0.870

    j = np.logical_and(Reu > 10, Reu <= 30)
    Ret[j] = 4.661 / Reu[j]**1.475
    Req[j] = 4.994 / Reu[j]**1.297

    j = Reu > 30
    Ret[j] = 34.904 / Reu[j]**2.067
    Req[j] = 30.790 / Reu[j]**1.845

    return Ret, Req


def psi(zet, utc=True):
    """Computes the turbulent velocity profile function following TOGA/COARE,
    given zet = (z/L), L the Monin-Obukoff length scale.  Edson and Fairall
    TOGA COARE code (version 2.0) as modified to include Rogers' weighting
    factor to combine the Dyer and free convection forms for unstable
    conditions.
    """
    zet = np.asanyarray(zet)
    zet = np.atleast_1d(zet)

    c13 = 1./3.
    sq3 = np.sqrt(3.0)

    # Stable conditions.
    y = -4.7 * zet

    # Unstable conditions.
    j = zet < 0
    zneg = zet[j]

    # Nearly stable (standard functions).
    x = (1 - 16.0 * zneg)**0.25
    if utc:
        y1 = (2.0 * np.log((1 + x) / 2) + np.log((1 + x**2) / 2) - 2 *
              np.arctan(x) + np.pi/2)
    else:
        y1 = 2.0 * np.log((1 + x**2) / 2)

    # Free convective limit.
    x = (1 - 12.87 * zneg)**c13
    y2 = (1.5 * np.log((x**2 + x + 1) / 3) - sq3 *
          np.arctan((2 * x + 1) / sq3) + np.pi / sq3)

    # Weighted sum of the two.
    F = 1.0 / (1 + zneg**2)
    y[j] = F * y1 + (1 - F) * y2

    return y


@ensure_atleast_1d
def hfbulktc(ur, zr, Ta, zt, rh, zq, Pa, Ts, **kw):
    """Computes sensible and latent heat fluxes and other variables.

            Hs      = sensible heat flux INTO ocean [W m**-2]
            Hl      = latent heat flux INTO ocean [W m**-2]
            Hl_webb = Webb correction to latent heat flux INTO ocean [W m**-2]
            stress  = wind stress [N/m^2]
            U_star  = velocity friction scale [m s-1]
            T_star  = temperature scale [deg C]
            Q_star  = humidity scale [kg/kg]
            L       = Monin-Obukhov length [m]
            zetu    = zr/L
            CD      = drag coefficient
            CT      = temperature transfer coefficient (Stanton number)
            CQ      = moisture transfer coefficient (Dalton number)
            RI      = bulk Richardson number
            Dter    = cool-skin temperature difference (optional output) [C]
                      positive if surface is cooler than bulk (presently no
                      warm skin permitted by model)

    Based on the following buoy input data:

              ur     = wind speed [m s-1] measured at height zr [m]
              Ta     = air temperature [C] measured at height zt [m]
              rh     = relative humidity [%] measured at height zq [m]
              Pa     = air pressure [mb]
              Ts     = sea surface temperature [C]
              sal    = salinity [psu (PSS-78)]
                       (optional - only needed for cool-skin)
              dlw    = downwelling (INTO water) longwave radiation [W m**-2]
                       (optional - only needed for cool-skin)
              dsw    = measured insolation [W m**-2]
                       (optional - only needed for cool-skin)
              nsw    = net shortwave radiation INTO the water [W m**-2]
                       (optional - only needed for cool-skin)

    where ur, Ta, rh, Pa, Ts, zr, zt, and zq (and optional sal, dlw,
    dsw, and nsw) may be either row or column vectors and rh, Pa,
    zr, zt, and zq (and optional sal) may also be fixed scalars.

    Output variables are given as column vectors in A:

    1) without cool-skin correction:

      A=[Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI]

    2) with cool-skin correction:

      A=[Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI Dter]

    Code follows Edson and Fairall TOGA COARE code (version 2.0), modified
    to include Rogers' weighting factor for unstable conditions.  Code does
    include gustiness, and assumes that the marine boundary layer height is
    known and constant over time for simplicity. zr/L is limited to
    be <=3.0 to ensure that the code converges to nonzero stress and heat
    flux values for strongly stable conditions.  The bulk Richardson number
    is computed between the sea surface and zr as a diagnostic about whether
    turbulent boundary layer theory is applicable.  Code does not include
    warm layer effects to modify Ts.  See Fairall et al (1996), J. Geophys.
    Res., 101, 3747-3764, for description of full TOGA COARE code and
    comparison with data.

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    08-19-1998: version 1.1 (rewritten by RP to remove inconsistencies in
             virtual and real temperatures, improve loop structure,
             correct gustiness component of stress computation)
    04-09-1999: version 1.2 (rewritten by AA to clarify some variable names
            and include cool-skin effect and Webb correction to latent
            heat flux added to output matrix)
    08-05-1999: version 2.0
    11-26-2010: Filipe Fernandes, Python translation.
    """

    # Optional cool-skin stuff
    sal = kw.pop('sal', None)
    dlw = kw.pop('dlw', None)
    dsw = kw.pop('dsw', None)
    nsw = kw.pop('nsw', None)

    # Initialize various constants.
    tol = 0.001  # Tolerance on Re changes to make sure soln has converged.
    onethird = 1/3.
    o61 = 1 / eps_air-1  # 0.61 Moisture correction for temperature)

    visc = viscair(Ta)  # Viscosity.
    Qsats = Qsat_coeff * qsat(Ts, Pa)  # Saturation specific humidity.
    Q = (0.01 * rh) * qsat(Ta, Pa)  # Specific humidity of air [kg/kg].
    T = Ta + CtoK
    Tv = T * (1 + o61 * Q)  # Air virtual temperature.
    rho = (100 * Pa) / (gas_const_R * Tv)  # Air density.
    Dt = (Ta + 0.0098 * zt) - Ts  # Adiabatic temperature difference.
    Dq = Q - Qsats  # Humidity difference.

    # Compute initial neutral scaling coefficients.
    S = np.sqrt(ur**2 + min_gustiness**2)
    # Smith's neutral cd as first guess.
    cdnhf = np.sqrt(cdn(S, zr, drag='smith', Ta=Ta))

    z0t = 7.5 * 1e-5
    ctnhf = kappa / np.log(zt/z0t)

    z0q = z0t.copy()
    cqnhf = kappa / np.log(zq/z0q)

    U_star = cdnhf * S  # Includes gustiness.
    T_star = ctnhf * Dt
    Q_star = cqnhf * Dq

    Dter = 0
    Dqer = 0
    if sal:
        # Initial cool-skin thickness guess.
        delta = 0.001 * np.ones_like(Ts)

    Reu = Ret = Req = 0
    # Begin iteration loop to compute best U_star, T_star, and Q_star.
    for iter in range(80):
        # Save old values.
        ReuO = Reu.copy()
        RetO = Ret.copy()
        ReqO = Req.copy()

        # Compute Monin-Obukov length (NB - definition given as eqn (7)
        # of Fairall et al (1996) probably wrong, following, e.g.
        # Godfrey and Bellars (1991), JGR, 96, 22043-22048 and original code).
        bs = g * (T_star * (1 + o61 * Q) + o61 * T * Q_star) / Tv
        L = (U_star**2) / (kappa * bs)

        # Set upper limit on zr/L = 3.0 to force convergence under.
        # very stable conditions. Assume that zr, zt and zq comparable.
        index_limit = np.logical_and(L < zr/3, L > 0)
        L[index_limit] = zr[index_limit] / 3
        # Non-dimensional heights.
        zetu = zr/L
        zett = zt/L
        zetq = zq/L

        # Surface roughness.
        z0 = (Charnock_alpha / g) * U_star**2 + R_roughness * visc/U_star

        # Compute U_star correction for non-neutral conditions.
        cdnhf = kappa / (np.log(zr/z0) - psi(zetu, utc=True))
        U_star = cdnhf * S

        Reu = z0*U_star/visc  # Roughness Reynolds.
        Ret, Req = LKB(Reu)  # Compute other roughness Reynolds.

        # Compute t and q roughness scales from roughness R#s.
        z0t = visc * Ret / U_star
        z0q = visc * Req / U_star

        # Compute new transfer coefficients at measurement heights.
        cthf = kappa / (np.log(zt / z0t) - psi(zett, utc=False))
        cqhf = kappa / (np.log(zq / z0q) - psi(zetq, utc=False))

        # Compute new values of T_star, Q_star.
        T_star = cthf * (Dt + Dter)
        Q_star = cqhf * (Dq + Dqer)

        # Estimate new gustiness.
        Ws = U_star * (-CVB_depth / (kappa*L))**onethird
        wg = min_gustiness * np.ones_like(ur)
        j = zetu < 0  # Convection in unstable conditions only.
        wg[j] = np.max(min_gustiness, beta_conv * Ws[j])  # Set min gustiness.
        S = np.sqrt(ur**2 + wg**2)

        if sal:
            # Compute cool-skin parameters.
            delta, Dter, Dqer = cool_skin(sal, Ts-Dter, rho, cp, Pa, U_star,
                                          T_star, Q_star, dlw, dsw, nsw, delta,
                                          g, gas_const_R, CtoK, Qsat_coeff)

    ii = np.logical_or(np.abs(Reu-ReuO) > tol, np.abs(Ret-RetO) > tol)
    ii = np.logical_or(ii, np.abs(Req-ReqO) > tol)
    if ii.any():
        print('Algorithm did not converge for %s values.'
              ' Indices are: %s' % (ii.sum()), np.where(ii))
        warn('Not converged!')

    # Compute latent heat.
    Le = (2.501-0.00237 * (Ts-Dter)) * 1e6

    # Compute fluxes into ocean.
    Hs = rho * cp * U_star * T_star
    Hl = rho * Le * U_star * Q_star

    # compute transfer coefficients at measurement heights
    CD = (U_star / S)**2
    CT = U_star * T_star / (S * (Dt + Dter))  # Stanton number.
    CQ = U_star * Q_star / (S * (Dq + Dqer))  # Dalton number.

    # To compute mean stress, we don't want to include the effects of
    # gustiness which average out (in a vector sense).
    stress = rho * CD * S * ur

    # Compute bulk Richardson number (as a diagnostic) - the "T" is probably
    # not quite right - assumes T \ approx. Ts (good enough though).
    RI = g * zr * ((Dt + Dter) + o61 * T * (Dq + Dqer)) / (Tv * S**2)

    # Compute Webb correction to latent heat flux into ocean.
    # eq. 21.
    W = 1.61 * U_star * Q_star + (1 + 1.61 * Q) * U_star * T_star / T
    # eq. 22, Fairall et al. (1996), JGR, 101, p3751.
    Hl_webb = rho * Le * W * Q

    # Output array.
    if sal:
        # Output additional cool-skin parameters.
        return (Hs, Hl, Hl_webb, stress, U_star, T_star,
                Q_star, L, zetu, CD, CT, CQ, RI, Dter)
    else:
        return (Hs, Hl, Hl_webb, stress, U_star, T_star,
                Q_star, L, zetu, CD, CT, CQ, RI)


@ensure_atleast_1d
def cool_skin(sal, Tw, rhoa, cpa, Pa, U_star, T_star, Q_star,  dlw, dsw, nsw,
              delta, g, Rgas, CtoK, Qsat_coeff):
    """Computes the cool-skin parameters.  This code follows the FORTRAN
    program bulk_v25b.f.  For more details, see the cool-skin and warm layer
    paper by Fairall et al (1996), JGR, 101, 1295-1308.  All input variables
    should be vectors (either row or column), except Rgas, CtoK, Qsat_coeff,
    and g, which can be scalars.

    Parameter
    ---------
    sal : salinity [psu (PSS-78)]
    Tw : water surface temperature [C]
    rhoa : air density [kg/m^3]
    cpa : specific heat capacity of air [J/kg/C]
    Pa : air pressure [mb]
    U_star : friction velocity including gustiness [m s-1]
    T_star : temperature scale [C]
    Q_star : humidity scale [kg/kg]
    dlw : downwelling (INTO water) longwave radiation [W m**-2]
    dsw : measured insolation [W m**-2]
    nsw : net shortwave radiation INTO water [W m**-2]
    delta : cool-skin layer thickness [m]
    g : gravitational constant [m s**-2]
    Rgas : gas constant for dry air [J/kg/K]
    CtoK : conversion factor for deg C to K
    Qsat_coeff : saturation specific humidity coefficient

    Returns
    -------
    delta : cool-skin layer thickness [m]
    Dter : cool-skin temperature difference [C]; positive if
           surface is cooler than bulk (presently no warm skin
           permitted by model)
    Dqer : cool-skin specific humidity difference [kg/kg]

    04-09-1999: version 1.2 (contributed by AA)
    08-05-1999: version 2.0
    """

    alpha = sw.alpha(sal, Tw, 0)  # Thermal expansion coeff [1/C].
    beta_sal = sw.beta(sal, Tw, 0)  # Saline contraction coeff [1/psu].
    cpw = sw.cp(sal, Tw, 0)  # Specific heat capacity  [J/kg/C].
    rhow = sw.dens0(sal, Tw)  # Density at atmospheric press [kg/m^3].
    viscw = visc(sal, Tw, 0)  # kinematic viscosity of sea-water [m^2/s].
    tcondw = tcond(sal, Tw, 0)  # thermal conductivity of sea-water [W/m/K].

    # The following are values used for COARE.
    # alpha    = 2.1e-5*(Tw+3.2)**0.79
    # beta_sal = 0.026 / sal
    # cpw      = 4000
    # rhow     = 1022
    # viscw    = 1e-6
    # tcondw   = 0.6

    # Latent heat of water.
    Le = (2.501 - 0.00237 * Tw) * 1e6

    # Saturation specific humidity;.
    Qs = Qsat_coeff * qsat(Tw, Pa)

    # A big constant.
    bigc = (16 * g * cpw * (rhow * viscw)**3) / (tcondw**2 * rhoa**2)

    # Constant for correction of dq; slope of sat. vap.
    wetc = 0.622 * Le * Qs / (Rgas * (Tw + CtoK)**2)

    # Compute fluxes out of the ocean (i.e., up = positive).
    hsb = -rhoa * cpa * U_star * T_star
    hlb = -rhoa * Le * U_star * Q_star

    # Net longwave (positive up).
    nlw = -lwhf(Tw, dlw, dsw)

    # Total heat flux out of the water surface.
    qout = nlw + hsb + hlb

    # Compute deltaSc = fc*Sns, see sec. 2.4 (p. 1297-1298) in cool-skin paper
    # 3 choices; comment out those that are not used!
    deltaSc = np.zeros_like(Tw)
    ipos_nsw = np.where(nsw > 0)[0]

    # Paulson and Simpson (1981)
    deltaSc[ipos_nsw] = f_c(delta[ipos_nsw], 1) * nsw[ipos_nsw]
    # COARE approx. to Paulson.
    # deltaSc[ipos_nsw] = f_c(delta[ipos_nsw], 2) * nsw[ipos_nsw]
    # Hasse (1971).
    # deltaSc[ipos_nsw] = f_c(delta[ipos_nsw], 3) * nsw[ipos_nsw]

    qcol = qout - deltaSc

    # Initialize.
    alphaQb = np.zeros_like(Tw)
    lamda = np.zeros_like(Tw)
    Dter = np.zeros_like(Tw)

    ipos_qcol = qcol > 0
    # eqn. 17 in cool-skin paper.
    alphaQb[ipos_qcol] = (alpha[ipos_qcol] * qcol[ipos_qcol] + sal[ipos_qcol] *
                          beta_sal[ipos_qcol] * hlb[ipos_qcol] *
                          cpw[ipos_qcol] / Le[ipos_qcol])

    # eqn. 14 in cool-skin paper.
    lamda[ipos_qcol] = (6 / (1 + (bigc[ipos_qcol] * alphaQb[ipos_qcol] /
                                  U_star[ipos_qcol]**4)**0.75)**(1/3))

    # eqn. 12 in cool-skin paper.
    delta[ipos_qcol] = (lamda[ipos_qcol] * viscw[ipos_qcol] /
                        (np.sqrt(rhoa[ipos_qcol] / rhow[ipos_qcol]) *
                         U_star[ipos_qcol]))

    # eqn. 13 in cool-skin paper.
    Dter[ipos_qcol] = qcol[ipos_qcol] * delta[ipos_qcol] / tcondw[ipos_qcol]
    Dqer = wetc * Dter

    return delta, Dter, Dqer


@ensure_atleast_1d
def f_c(delta, option=1):
    """
    F_C: Computes the absorption coefficient fc.
    fc = F_C(delta, option) computes the absorption coefficient fc.

    Parameter
    ---------
    delta : Thickness of cool-skin layer [m]
            option  -  1  Use Paulson and Simpson (1981) data for seawater;
                          See also p. 1298 of Fairall et al (1996) JGR, 101,
                          cool-skin and warm-layer paper.
                       2  Use approximation to Paulson as given in
                          p. 1298 of Fairall et al (1996) JGR, 101, cool-skin
                          and warm-layer paper.
                       3  Use fc = const. = 0.19, as suggested by Hasse (1971).
    08-05-1999: version 1.2 (contributed by AA)
    """

    if option == 1:
        # Use Paulson and Simpson data
        # Wavelength bands for the coefficients [um]
        # 1) 0.2-0.6
        # 2) 0.6-0.9
        # 3) 0.9-1.2
        # 4) 1.2-1.5
        # 5) 1.5-1.8
        # 6) 1.8-2.1
        # 7) 2.1-2.4
        # 8) 2.4-2.7
        # 9) 2.7-3.0
        # F_i is the amplitude
        F_i = np.array([0.237, 0.360, 0.179, 0.087,
                        0.080, 0.0246, 0.025, 0.007, 0.0004])
        # F_i1 = repmat(F_i, len(delta), 1)
        # Gam_i is the absorption length [m].
        Gam_i = np.array([34.8, 2.27, 3.15e-2, 5.48e-3, 8.32e-4,
                          1.26e-4, 3.13e-4, 7.82e-5, 1.44e-5])
        # Gam_i1 = repmat(Gam_i, len(delta), 1)
        # delta1 = repmat(delta, 1, len(Gam_i))
        # fc is the absorption in the cool-skin layer of thickness delta.
        # fc = np.sum(F_i1 * (1 - (Gam_i1 / delta1) *
        # (1 - np.exp(-delta1 / Gam_i1))), 2)
        fc = np.sum(F_i * (1 - (Gam_i / delta[..., None]) *
                           (1 - np.exp(-delta[..., None] / Gam_i))), axis=1)
    elif option == 2:
        # Use COARE approximation to Paulson and Simpson data.
        fc = 0.137 + 11 * delta - (6.6e-5 / delta) * (1-np.exp(-delta/8e-4))
    elif option == 3:
        # Use Hasse simple approximation.
        fc = 0.19
    return fc


@ensure_atleast_1d
def lwhf(Ts, dlw, dsw=None):
    """
    LWHF: computes net longwave heat flux following Dickey et al (1994).
    qlw=LWHF(Ts,dlw) computes the net longwave heat flux into the ocean.
    Following Dickey et al (1994), J. Atmos. Oceanic Tech., 11, 1057-1078,
    the incident longwave flux can be corrected for sensor heating due to
    insolation if you are using Epply or Kipp & Zonen CG1 pyrgeometers.
    In this case, use qlw=LWHF(Ts,dlw,dsw). Epply is the default
    pyrgeometer; change code for the Kipp & Zonen instrument.

    Parameters
    ----------
    Ts : sea surface temperature [C]
    dlw : (measured) downward longwave flux [W m**-2]
    dsw : (measured) insolation [W m**-2] (needed for Eppley
          or Kipp & Zonen pyrgeometers)

    Returns
    qlw : net longwave heat flux [W m**-2]

    03-08-1997: version 1.0
    08-19-1998: version 1.1 (revised for non-Epply pyrgeometers by RP)
    04-09-1999: version 1.2 (included Kipp & Zonen CG1 pyrgeometers by AA)
    08-05-1999: version 2.0
    """
    Ts, dlw = np.broadcast_arrays(Ts, dlw)

    # Convert degC to degK.
    ts = Ts + CtoK
    # Correct dlw for sensor heating by insolation.
    if dsw:
        # This line is for Epply pyrgeometers.
        dlwc = dlw-0.036 * dsw
        # This line is for Kipp & Zonen CG1 pyrgeometers
        # (the offset is specified as 25 W m**-2 at 1000 W m**-2).
        # dlwc = dlw-0.025 * dsw
    else:
        dlwc = dlw

    # Compute upward gray-body longwave flux.
    lwup = -emiss_lw * sigmaSB * (ts**4)
    # Compute net flux.
    qlw = lwup + emiss_lw * dlwc
    return qlw
