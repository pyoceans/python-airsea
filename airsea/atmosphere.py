# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np
from .constants import eps_air, gas_const_R, CtoK, P_default

"""Most of these functions came from the Matlab air sea toolbox."""


def vapor(Tw):
    """Calculates heat of evaporation for pure water.

    Parameters
    ----------
    Tw : array_like
        water temperature [:math:`^\\circ` C]

    Returns
    -------
    L : array_like
        heat of evaporation [J kg :sup:`-1`]

    Notes
    -----
    This can be used to compute the fresh water flux from latent heat flux.
    Range of validity: 0 <= t <= 100 :math:`^\\circ` C.  No formulas are known
    to be available for the change of the heat of evaporation as function of
    salinity.

    Examples
    --------
    >>> from airsea import atmosphere as asea
    >>> asea.vapor([0.0, 5., 21., 33., 100.])
    array([ 2500390.     ,  2488557.86125,  2450741.02909,  2422299.01153,
            2256560.     ])

    References
    ----------
    .. [1] Landolt-Bornstein, Numerical Data and Functional Relationships in
    Science and Technology. New Series, Sundermann, J. (editor), vol. 3,
    subvol. a, Springer-Verlag, p. 256.

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    01/22/1999: version 1.0 (contributed by RO)
    08/05/1999: version 2.0
    11/26/2010: Filipe Fernandes, Python translation.
    """

    Tw = np.asarray(Tw)

    a0 = 2.50039e6
    a1 = -2.3683e3
    a2 = 4.31e-1
    a3 = -1.131e-2

    return a0 + a1 * Tw + a2 * Tw ** 2 + a3 * Tw ** 3


def air_dens(Ta, rh, Pa=P_default):
    """Computes the density of moist air.

    Parameters
    ----------
    Ta : array_like
        air temperature [:math:`^\\circ` C]
    rh : array_like
        relative humidity [percent]
    Pa : array_like, optional
        air pressure [mb]

    Returns
    -------
    rhoa : array_like
        air density [kg m :sup:`-3`]

    See Also
    --------
    TODO: qsat

    Examples
    --------
    >>> from airsea import atmosphere as asea
    >>> asea.air_dens([5., 15., 23.1], 95.)
    array([ 1.27361105,  1.22578105,  1.18750778])
    >>> asea.air_dens([5., 15., 23.1], [95, 100, 50], 900)
    array([ 1.12331233,  1.08031123,  1.05203796])

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    04/07/1999: version 1.2 (contributed by AA)
    08/05/1999: version 2.0
    11/26/2010: Filipe Fernandes, Python translation.
    """

    Ta, rh, Pa = np.asarray(Ta), np.asarray(rh), np.asarray(Pa)

    o61 = 1 / eps_air - 1  # 0.61 -> Moisture correction for temperature
    Q = (0.01 * rh) * qsat(Ta, Pa)  # Specific humidity of air [kg/kg]
    T = Ta + CtoK
    Tv = T * (1 + o61 * Q)  # Air virtual temperature.
    rhoa = (100 * Pa) / (gas_const_R * Tv)  # Air density [kg/m**3]

    return rhoa


def qsat(Ta, Pa=P_default, sflag=True):
    """Computes the specific humidity [kg/kg] at saturation at air temperature
    Ta [deg C].  Dependence on air pressure, Pa, is small, but is included as
    an optional input.

    Parameters
    ----------
    Ta : array_like
        air temperature [:math:`^\\circ` C]
    Pa : array_like, optional
        air pressure [mb]
    sflag : bool
        Salt water flag, True for Salt Water and False for Fresh Water.

    Returns
    -------
    q : array_like
        saturation specific humidity  [kg/kg]

    Examples
    --------
    >>> from airsea import atmosphere as asea
    >>> asea.qsat([10., 11.3, 15.5, 21.])
    array([ 0.00755174,  0.00823827,  0.0108498 ,  0.01536531])
    >>> asea.qsat([10., 11.3, 15.5, 21.], 910.)
    array([ 0.00846606,  0.00923618,  0.01216637,  0.01723551])

    Notes
    -----
    Version 1.0 used Tetens' formula for saturation vapor pressure from Buck
    (1981), J. App. Meteor., 1527-1532.  This version follows the saturation
    specific humidity computation in the COARE Fortran code v2.5b.
    This results in an increase of ~5% in latent heat flux compared to the
    calculation with version 1.0.

    .. [1] Buck (1981), J. App. Meteor., 1527-1532

    1997/08/03: version 1.0
    1999/07/04: version 1.2 (revised as above by AA)
    1990/05/08: version 2.0
    2011/10/28: Python version
    """
    
    Ta, Pa = np.asarray(Ta), np.asarray(Pa)
    
    # Original code
    # a = (1.004 * 6.112 * 0.6220) / Pa
    # q = a * np.exp((17.502 * Ta) / (240.97 + Ta))

    # Fortran code for COARE v2.5b
    ew = (6.1121 * (1.0007 + 3.46e-6 * Pa) *
          np.exp((17.502 * Ta) / (240.97 + Ta)))  # [mb]
    qsat = 0.62197 * (ew / (Pa - 0.378 * ew))  # [mb] -> [g/kg]
    qsat = (1.0 - 0.02 * sflag) * qsat  # flag for fresh (0) or salt (1) water

    # When it is multiplied by 0,98, i.e. sflag = 1, and T is SST, it is eq.3
    # from Fairall et al. 1996.  When sflag = 0 and T is air, it is eq.4 from
    # Fairall et al. 1996.

    return qsat


def satvap(Ta, Pa=P_default):
    """Computes saturation vapor pressure.

    Parameters
    ----------
    Ta : array_like
        air temperature [:math:`^\\circ` C]
    Pa : float, optional
        air pressure [mb]

    Returns
    -------
    q : array_like
         saturation vapor pressure [mb]

    See Also
    --------
    TODO: blwhf, relhumid

    Examples
    --------
    >>> from airsea import atmosphere as asea
    >>> asea.satvap(10.)
    12.324209279140598
    >>> asea.satvap([15., 25., 40.])
    array([ 17.11594161,  31.823173  ,  74.32668994])
    >>> asea.satvap([15., 25., 40.], 900.)
    array([ 17.10646652,  31.80464856,  74.2782608 ])
    >>> asea.satvap([15., 25., 40.],[900.,800., 1030])
    array([ 17.10646652,  31.78921152,  74.3307257 ])

    References
    ----------
    .. [1] Gill (1982), Atmos-Ocean Dynamics, 606.

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    03/08/1997: version 1.0
    08/27/1998: version 1.1 (corrected by RP)
    08/05/1999: version 2.0
    11/26/2010: Filipe Fernandes, Python translation.
    """

    Ta, Pa = np.asarray(Ta), np.asarray(Pa)
    ew = np.power(10., ((0.7859 + 0.03477 * Ta) / (1 + 0.00412 * Ta)))
    fw = 1 + 1e-6 * Pa * (4.5 + 0.0006 * Ta ** 2)
    return fw * ew


def rhadj(rh, rhmax):
    """Rescales RH to have a maximum of 100%.

    Parameters
    ----------
    rh : array_like
         relative humidity [percent]
    rhmax : float
            maximum relative humidity value measured [percent]

    Returns
    -------
    rhadj : array_like
            adjusted relative humidity [percent]

    Notes
    -----
    Re-scaling rh is needed to avoid the maximum values to exceed 100%. Assumes
    values between 93% and rhmax should be re-scaled to 93 - 100%.

    The calibration curves of rh sensors usually become nonlinear above ~ 90%,
    and may peak above 100% in this range above ~ 90% where their calibration
    becomes unstable.

    Examples
    --------
    >>> from airsea import atmosphere as asea
    >>> asea.rhadj([98., 95., 94., 93., 10., 20.], 98.)
    array([ 100. ,   95.8,   94.4,   93. ,   10. ,   20. ])

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    04/10/1998: version 1.0
    08/05/1999: version 2.0
    11/26/2010: Filipe Fernandes, Python translation.
    """

    rh, rhmax = np.asarray(rh), np.asarray(rhmax)

    rhl = 93.
    rhn = rh
    a = (100. - rhl) / (rhmax - rhl)
    drh = rh - rhl
    rhn[drh > 0] = rhl + a * drh[drh > 0]

    return rhn


def ep(rfd, Qlat, dt=3600):
    """Computes precipitation and evaporation accumulation from rainfall rate
    (rfd) and latent heat (Qlat).

    Parameters
    ----------
    rfd : array_like
          precipitation rate [ m s :sup:`-1`]
    Qlat : array_like
           latent heat flux [W m :sup:`-2`]
    dt : float
         delta t for series in seconds, assumes hourly data

    Returns
    -------
    P : array_like
        precipitation accumulation [m]
    E : array_like
        evaporation accumulation [m]

    Notes
    -----
    Convert precipitation from [mm min :sup:`-1`] to [ m s :sup:`-1`]

    Examples
    --------
    >>> import numpy as np
    >>> Qlat = [550., 450., 350.]
    >>> rfd = [10., 15., 35.]
    >>> from airsea import atmosphere as asea
    >>> asea.ep(np.array(rfd)/1000/60, Qlat)[0] # Evaporation
    array([ 0.00077268,  0.00140488,  0.00189659])
    >>> asea.ep(np.array(rfd)/1000/60, Qlat)[1] # Precipitation
    array([ 0.6,  1.5,  3.6])

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    08/05/1999: version 2.0
    11/26/2010: Filipe Fernandes, Python translation.
    """

    rfd, Qlat = np.asarray(rfd), np.asarray(Qlat)

    P = np.cumsum(rfd)
    P = P * dt

    Le = 2.5e6  # heat of vaporization [W/m**2]
    pw = 1025.  # density of seawater [kg/m**3] at 32 psu, 10 degC, 0 db
    E = np.cumsum(Qlat) * dt / (Le * pw)

    return E, P


def relhumid(Td, Tw, Pa=1020, p_typ='screen'):
    """Finds relative humidity from wet/dry thermometer readings using the
    psychrometric equation.

    # TODO: http://www.ncl.ucar.edu/Document/Functions/Built-in/relhum.shtml

    Parameters
    ----------
    Td : array_like
         dry bulb thermometer [:math:`^\\circ` C]
    Tw : array_like
         Wet thermometer [:math:`^\\circ` C]
    Pa : array_like, optional
        air pressure [mb]
    p_typ : str or int, optional
            'assman' for Assman-type forced ventilation
            'screen' for standard screen (natural ventilation)
            Default is 'screen'.

    Returns
    -------
    rh : array_like
         relative humidity  [%]

    See Also
    --------
    TODO: satvap

    Examples
    --------
    >>> from airsea import atmosphere as asea
    >>> asea.relhumid([10., 15., 30.1], [9.2, 10.1, 33.0], 900)
    array([  90.09582024,   51.89065871,  122.80282659])
    >>> asea.relhumid([10., 15., 30.1], [9.2, 10.1, 33.0], p_typ='assman')
    array([  90.34986252,   53.01105734,  122.53946619])

    References
    ----------
    .. [1] Sargent (1980), Meteorol. Mag. 109, 238-246.

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    08/28/1998: version 1.1 (contributed by RP)
    08/05/1999: version 2.0
    11/26/2010: Filipe Fernandes, Python translation.
    """
    Td, Tw, Pa = np.asarray(Td), np.asarray(Tw), np.asarray(Pa)

    # Psychrometric coefficient.
    if p_typ is 'screen':
        A = 0.000799  # Natural screens.
    elif p_typ is 'assman':
        A = 0.000667  # Assmann-type with forced ventilation.
    else:
        print('unknown psychrometer type: %s' % p_typ)
        A = np.NaN  # FIXME: Add an proper python error

    # Compute saturation vapor pressure for both temperatures.
    ed = satvap(Td, Pa)
    ewp = satvap(Tw, Pa)

    # The psychrometric equation.
    e = ewp - A * Pa * (Td - Tw)  # Ambient vapor pressure.
    rh = e / ed * 100
    return rh


def cloudcor(C, optns, lat):
    """Computes cloud correction factor for bulk long-wave flux as a function of
    the cloud fraction "C" for bulk long-wave flux formulae.

    Parameters
    ----------
    C : array_like FIXME
        cloud fraction FIXME
    optns : str, int
            'clarke', Clarke (1974) corrections for np.abs(lat)<50.
            'bunker', Bunker (1976) corrections for North Atlantic.
            [a1 a2] -> use a correction factor of [1 - a1*C - a2*C**2 ].
    lat : array_like FIXME: or float
          latitude in decimal degrees north [-90..+90].

    Returns
    -------
    Fc : array_like
         correction factor used as input by "blwhf"

    See Also
    --------
    TODO: None

    Notes
    -----
    In general, these are functions of the form
                1 - a_n * C**n
    Since the coefficients and powers depend a lot on the dominant cloud type
    which may vary from region to region and season to season, it is not clear
    which parametrization is best.

    The particular parametrization used here depends on the second input
    variable, for which no default is given to emphasize the fact that you
    really need to understand what you are doing here!

    There are several "built-in" formulae (from Fung et al) that all have a
    latitude-dependence of some kind.

    Examples
    --------
    >>> from airsea import atmosphere as asea
    >>> asea.cloudcor([0.8, 0.1, 1], 'clarke', 47)
    array([ 0.5328,  0.9927,  0.27  ])
    >>> asea.cloudcor(0.5, [0.2, 0.4], 47)
    0.80000000000000004

    References
    ----------
    .. [1] Fung et al (1984), Rev. of Geophys. and Space Phys., 22, 177-193.
    .. [2] Clarke (1974) FIXME
    .. [3] Bunker (1976) FIXME

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    03/12/1998: version 1.1 (contributed by RP)
    08/05/1999: version 2.0
    11/26/2010: Filipe Fernandes, Python translation.
    """
    # convert input to numpy array
    C, lat = np.asarray(C), np.asarray(lat)

    if type(optns) == str:
        if optns == 'clarke':
            a1 = 0
            if np.abs(lat) > 55:
                a2 = np.NaN
            elif np.abs(lat) > 45:
                a2 = 0.73
            elif np.abs(lat) > 35:
                a2 = 0.69
            elif np.abs(lat) > 25:
                a2 = 0.64
            elif np.abs(lat) > 15:
                a2 = 0.60
            elif np.abs(lat) > 7:
                a2 = 0.56
            elif np.abs(lat) > 2:
                a2 = 0.53
            else:
                a2 = 0.51
        elif optns == 'bunker':
            a2 = 0
            if np.abs(lat) > 75:
                a1 = 0.84
            elif np.abs(lat) > 65:
                a1 = 0.80
            elif np.abs(lat) > 55:
                a1 = 0.76
            elif np.abs(lat) > 45:
                a1 = 0.72
            elif np.abs(lat) > 35:
                a1 = 0.68
            elif np.abs(lat) > 25:
                a1 = 0.63
            elif np.abs(lat) > 15:
                a1 = 0.59
            elif np.abs(lat) > 7:
                a1 = 0.52
            else:
                a1 = 0.50
        else:
            print('Unrecognized option')  # Add a proper python error.
    else:
        a1, a2 = np.asarray(optns[0]), np.asarray(optns[1])

    Fc = 1 - a1 * C - a2 * C ** 2
    return Fc


def visc_air(Ta):
    """Computes the kinematic viscosity of dry air as a function of air
    temperature

    Parameters
    ----------
    Ta : array_like
         air temperature [:math:`^\\circ` C]

    Returns
    -------
    visa : array_like
           [m :sup:`2` s :sup:`-1`]

    See Also
    --------
    hfbulktc, cdn
    sw.visc

    Notes
    -----
    sw.visc from python seawater package

    Examples
    --------
    >>> from airsea import atmosphere as asea
    >>> asea.visc_air([[0.1, 5., 15],[22.8, 28.9, 31.4]])
    array([[  1.32686758e-05,   1.36964784e-05,   1.45857532e-05],
           [  1.52942886e-05,   1.58573695e-05,   1.60903922e-05]])

    References
    ----------
    .. [1] Andreas (1989), CRREL Report 89-11.

    Modifications: Original from COARE 3.0
    11/26/2010: Filipe Fernandes, Python translation.
    """
    # convert input to numpy array
    Ta = np.asarray(Ta)

    visa = 1.326e-5 * (1 + 6.542e-3 * Ta + 8.301e-6 * Ta ** 2 - 4.84e-9 *
                       Ta ** 3)
    return visa
