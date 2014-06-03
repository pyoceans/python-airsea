# -*- coding: utf-8 -*-
#
# windstress.py
#
# purpose:
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.github.io/
# created:  21-Aug-2013
# modified: Thu 08 May 2014 07:59:14 AM BRT
#
# obs:
#

import numpy as np

from .constants import kappa
from .atmosphere import visc_air


def cdn(sp, z, drag='largepond', Ta=10):
    """Computes neutral drag coefficient.
    Methods available are: Large & Pond (1981),  Vera (1983) or Smith (1988)

    Parameters
    ----------
    sp : array_like
         wind speed [m s :sup:`-1`]
    z : float, array_like
        measurement height [m]
    drag : str
           neutral drag by:
           'largepond' <-- default
           'smith'
           'vera'
    Ta : array_like, optional for drag='smith'
         air temperature [:math:`^\\circ` C]

    Returns
    -------
    cd : float, array_like
         neutral drag coefficient at 10 m
    u10 : array_like
          wind speed at 10 m [m s :sup:`-1`]

    See Also
    --------
    stress, spshft, visc_air

    Notes
    -----
    Vera (1983): range of fit to data is 1 to 25 [m s :sup:`-1`].

    Examples
    --------
    >>> from airsea import windstress as ws
    >>> ws.cdn([10., 0.2, 12., 20., 30., 50.], 10)
    (array([ 0.00115,  0.00115,  0.00127,  0.00179,  0.00244,  0.00374]),
     array([ 10. ,   0.2,  12. ,  20. ,  30. ,  50. ]))
    >>> ws.cdn([10., 0.2, 12., 20., 30., 50.], 15, 'vera')
    (array([ 0.00116157,  0.01545237,  0.00126151,  0.00174946,  0.00242021,
            0.00379521]),
     array([  9.66606155,   0.17761896,  11.58297824, 19.18652915,
            28.5750255 ,  47.06117334]))
    >>> ws.cdn([10., 0.2, 12., 20., 30., 50.], 20, 'smith', 20.)
    (array([ 0.00126578,  0.00140818,  0.00136533,  0.00173801,  0.00217435,
            0.00304636]),
     array([  9.41928554,   0.18778865,  11.27787697,  18.65250005,
            27.75712916,  45.6352786 ]))

    References
    ----------
    .. [1] Large and Pond (1981), J. Phys. Oceanog., 11, 324-336.
    .. [2] Smith (1988), J. Geophys. Res., 93, 311-326.
    .. [3] E. Vera (1983) FIXME eqn. 8 in Large, Morzel, and Crawford (1995),
    J. Phys. Oceanog., 25, 2959-2971.

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    03-08-1997: version 1.0
    08-26-1998: version 1.1 (vectorized by RP)
    08-05-1999: version 2.0
    11-26-2010: Filipe Fernandes, Python translation.
    """
    # convert input to numpy array
    sp, z, Ta = np.asarray(sp), np.asarray(z), np.asarray(Ta)

    tol = 0.00001  # Iteration end point.

    if drag == 'largepond':
        a = np.log(z / 10.) / kappa  # Log-layer correction factor.
        u10o = np.zeros(sp.shape)
        cd = 1.15e-3 * np.ones(sp.shape)
        u10 = sp / (1 + a * np.sqrt(cd))
        ii = np.abs(u10 - u10o) > tol

        while np.any(ii):
            u10o = u10
            cd = (4.9e-4 + 6.5e-5 * u10o)  # Compute cd(u10).
            cd[u10o < 10.15385] = 1.15e-3
            u10 = sp / (1 + a * np.sqrt(cd))  # Next iteration.
            # Keep going until iteration converges.
            ii = np.abs(u10 - u10o) > tol

    elif drag == 'smith':
        visc = visc_air(Ta)

        # Remove any sp==0 to prevent division by zero
        i = np.nonzero(sp == 0)

        #sp[i] = 0.1 * np.ones(len(i)) FIXME

        # initial guess
        ustaro = np.zeros(sp.shape)
        ustarn = 0.036 * sp

        # iterate to find z0 and ustar
        ii = np.abs(ustarn - ustaro) > tol
        while np.any(ii):
            ustaro = ustarn
            z0 = Charnock_alpha * ustaro ** 2 / g + R_roughness * visc / ustaro
            ustarn = sp * (kappa / np.log(z / z0))
            ii = np.abs(ustarn - ustaro) > tol

        sqrcd = kappa / np.log(10. / z0)
        cd = sqrcd ** 2
        u10 = ustarn / sqrcd
    elif drag == 'vera':
        # constants in fit for drag coefficient
        A = 2.717e-3
        B = 0.142e-3
        C = 0.0764e-3

        a = np.log(z / 10.) / kappa  # Log-layer correction factor.
        # Don't start iteration at 0 to prevent blowups.
        u10o = np.zeros(sp.shape) + 0.1
        cd = A / u10o + B + C * u10o
        u10 = sp / (1 + a * np.sqrt(cd))

        ii = np.abs(u10 - u10o) > tol
        while np.any(ii):
            u10o = u10
            cd = A / u10o + B + C * u10o
            u10 = sp / (1 + a * np.sqrt(cd))  # Next iteration.
            # Keep going until iteration converges.
            ii = np.abs(u10 - u10o) > tol
    else:
        print('Unknown method')  # TODO: Add a proper python error.

    return cd, u10


def spshft(sp, z1, z2, drag='largepond', Ta=10.):
    """
    Adjusts wind speed from height z1 to z2. Methods available are: Large &
    Pond (1981),  Vera (1983) or Smith (1988).

    Parameters
    ----------
    sp : array_like
          wind speed [m s :sup:`-1`]
    z1 : float
         measurement height [m]
    z2 : float
         desired height [m]
    drag : str
           neutral drag by:
           'largepond' <-- default
           'smith'
           'vera'
    Ta : array_like
         air temperature [:math:`^\\circ` C]

    Returns
    -------
    sp_adj : array_like
          predicted wind speed [m s :sup:`-1`]
    ustar : array_like
            friction velocity [m s :sup:`-1`]

    See Also
    --------
    cdn

    Notes
    -----
    TODO

    Examples
    --------
    >>> from airsea import windstress as ws
    >>> ws.spshft([10., 0.2, 12., 20., 30., 50.], 10, 10)[0]
    array([ 10. ,   0.2,  12. ,  20. ,  30. ,  50. ])
    >>> from airsea import windstress as ws
    >>> ws.spshft([10., 0.2, 12., 20., 30., 50.], 10, 8, 'smith', 20)
    (array([  9.79908171,   0.19583568,  11.74922628,  19.52618419,
            29.20068179,  48.40456082]), array([ 0.3601597 ,  0.00746483,  0.44952896,  0.84934708,  1.43283229,
            2.8599333 ]))
    >>> ws.spshft([10., 0.2, 12., 20., 30., 50.], 15, 10, 'vera')
    (array([  9.66606155,   0.17761896,  11.58297824,  19.18652915,
            28.5750255 ,  47.06117334]), array([ 0.32943742,  0.02207938,  0.41140089,  0.80250639,  1.40576781,
            2.89921535]))

    References
    ----------
    .. [1] Large and Pond (1981), J. Phys. Oceanog., 11, 324-336.
    .. [2] Smith (1988), J. Geophys. Res., 93, 311-326.
    .. [3] E. Vera (1983) FIXME eqn. 8 in Large, Morzel, and Crawford (1995),
    J. Phys. Oceanog., 25, 2959-2971.

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    03-08-1997: version 1.0
    08-27-1998: version 1.1 (revised to use cdn* efficiently by RP)
    08-05-1999: version 2.0
    11-26-2010: Filipe Fernandes, Python translation.
    """

    z1, z2, sp, Ta = map(np.asarray, (z1, z2, sp, Ta))

    # find cd and ustar
    if drag == 'largepond':
        cd10, sp10 = cdn(sp, z1, 'largepond')
    elif drag == 'smith':
        cd10, sp10 = cdn(sp, z1, 'smith', Ta)
    elif drag == 'vera':
        cd10, sp10 = cdn(sp, z1, 'vera')
    else:
        print('Unknown method')  # TODO: add a proper python error

    ustar = np.sqrt(cd10) * sp10
    sp_adj = sp10 + ustar * np.log(z2 / 10.) / kappa
    return sp_adj, ustar


def stress(sp, z=10., drag='largepond', rho_air=1.22, Ta=10.):
    """
    Computes the neutral wind stress.

    Parameters
    ----------
    sp : array_like
         wind speed [m s :sup:`-1`]
    z : float, array_like, optional
        measurement height [m]
    rho_air : array_like, optional
           air density [kg m :sup:`-3`]
    drag : str
           neutral drag by:
           'largepond' <-- default
           'smith'
           'vera'
    Ta : array_like, optional
         air temperature [:math:`^\\circ` C]

    Returns
    -------
    tau : array_like
          wind stress  [N m :sup:`-2`]

    See Also
    --------
    cdn

    Notes
    -----
    TODO

    Examples
    --------
    >>> from airsea import windstress as ws
    >>> ws.stress([10., 0.2, 12., 20., 30., 50.], 10)
    array([  1.40300000e-01,   5.61200000e-05,   2.23113600e-01,
             8.73520000e-01,   2.67912000e+00,   1.14070000e+01])
    >>> ws.stress([10., 0.2, 12., 20., 30., 50.], 15, 'smith', rho_air=1.02, Ta=23.)
    array([  1.21440074e-01,   5.32531576e-05,   1.88322389e-01,
             6.62091968e-01,   1.85325310e+00,   7.15282267e+00])
    >>> ws.stress([10., 0.2, 12., 20., 30., 50.], 8, 'vera')
    array([  1.50603698e-01,   7.16568379e-04,   2.37758830e-01,
             9.42518454e-01,   3.01119044e+00,   1.36422742e+01])

    References
    ----------
    .. [1] Large and Pond (1981), J. Phys. Oceanog., 11, 324-336.
    .. [2] Smith (1988), J. Geophys. Res., 93, 311-326.
    .. [3] E. Vera (1983) FIXME eqn. 8 in Large, Morzel, and Crawford (1995),
    J. Phys. Oceanog., 25, 2959-2971.

    Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    03-08-1997: version 1.0
    08-26-1998: version 1.1 (revised by RP)
    04-02-1999: versin 1.2 (air density option added by AA)
    08-05-1999: version 2.0
    11-26-2010: Filipe Fernandes, Python translation.
    """
    z, sp, Ta, rho_air = map(np.asarray, (z, sp, Ta, rho_air))

    # find cd and ustar
    if drag == 'largepond':
        cd, sp = cdn(sp, z, 'largepond')
    elif drag == 'smith':
        cd, sp = cdn(sp, z, 'smith', Ta)
    elif drag == 'vera':
        cd, sp = cdn(sp, z, 'vera')
    else:
        print('Unknown method')  # TODO: add a proper python error

    tau = rho_air * (cd * sp ** 2)

    return tau

""" Wave effects on wind stress, higly experimental !!! """

#"""
#Wdnotes:
#-------
#Ua measured at height za for the effects of the wave boundary layer following the empirical model presented by Large, Morzel, and Crawford (1995), J. Phys. Oceanog., 25, 2959-2971. In particular, an analytic expression was found for the omega function (`omegalmc`) shown in their Fig. 9b, which allows the 'true' wind speed (Ut10) and stress at 10 m (assumed above the wave boundary layer height) to be computed using `wave10` and the true wind speed (Uta) at the measurement height za using `wave`. The Large et al model assumes neutral stability (reasonable for high winds and wave conditions) and uses a 10 m neutral drag law (`cdnve`) based on Vera (1983; unpublished manuscript). This drag law follows Large and Pond (1982) for winds above 10 m/s but increases at lower wind speeds like Smith (1987). The wave field is specified by the significant wave height Hw.

#To compute 'true' wind speed Uta at za given Hw, use:
#>>> Uta = wave(Ua, za, Hw)

#To compute 'true' wind speed Ut at 10 m given Hw, use:
#Ut10, (Ut10-U10) = wave10(Ua, za, Hw)

#FIXME: add this as a test?
#To plot the predicted effects of wave distortion on the wind Ua measured at the height za for a range of significant wave heights,
#>>> Hw = range(0,8,2) # [m]
#>>> y = waveplt(za)
#subroutines called:
#>>> y = omegalmc(x)
#>>> cd10 = cdnve(u10)
#"""

#def omegalmc(x):
    #"""
    #Computes the log profile correction function due to wind distortion associated with surface waves.

    #Parameters
    #----------
    #x : array_like
        #za / Hw, where za is the measurement height and Hw is the dominant surface wave height. [FIXME]
        #Assumes x is a vector with all elements greater than zr.

    #Returns
    #-------
    #y : array_like
        #log profile correction [FIXME]

    #See Also
    #--------
    #TODO: wave10, wave, wavdist

    #Notes
    #-----
    #Functional form is simplified (analytic) version of empirical omega curves shown in Fig. 9b of Large et al. 1995, with the wave-induced roughness length xr = 0.15.

    #Examples
    #--------
    #>>> omegalmc([TODO])
    #array([TODO])

    #References
    #----------
    #.. [1] Large, Morzel, and Crawford (1995), J. Phys. Oceanog., 25, 2959-2971

    #Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    #03/08/1997: version 1.0
    #08/05/1999: version 2.0
    #11/26/2010: Filipe Fernandes, Python translation.
    #"""
    ## convert input to numpy array
    #x = np.asarray(x)

    #xr     = 0.15
    #ylimit = -6
    #y      = ylimit * np.ones( x.shape )
    #i      = np.nonzero(x < 3.2967)
    ## polynomial fit
    #a   = -2.6
    #p1  = -0.0199
    #p2  =  0.0144
    #p3  =  0.7660
    #p4  =  0.0654

    #x2  = x[i]**2
    #x3  = x2 * x[i]

    #y[i] = a * np.log( x[i] / xr ) + p1 * x3 + p2 * x2 + p3 * x[i] + p4
    #return y

#def wave10(Ua, za, Hw):
    #"""
    #Computes true 10 m wind speed U10 using the wind speed measured Ua at the height za and measured wave height Hw and the neutral log profile corrected for the effects of low-level distortion of the wind profile by surface waves.

    #Parameters
    #----------
    #Ua : array_like
         #wind speed [m s :sup:`-1`]
    #za : float, array_like
         #wind measurement height [m]
    #Hw : float, array_like
         #wave height [m]

    #Returns
    #-------
    #U10 : array_like
          #true 10 m wind speed [m s :sup:`-1`]
    #delU : array_like
           #difference between true and uncorrected 10 m wind speed [m s :sup:`-1`]

    #Notes
    #-----
    #TODO

    #Examples
    #--------
    #TODO

    #References
    #----------
    #.. [1] Large, Morzel, and Crawford (1995), J. Phys. Ocean., 25, 2959-2971.

    #Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    #8/31/98: version 1.1
    #8/5/99: version 2.0
    #11/26/2010: Filipe Fernandes, Python translation.
    #"""
    ## convert input to numpy array
    #Ua = np.asarray(Ua) ; za = np.asarray(za) ; Hw = np.asarray(Hw)

    #tol = 0.001 # change in u10 to stop iteration
    #zs  = 10 # reference height
    #Xia = za / Hw
    #Xis = zs / Hw

    ## compute uncorrected 10 m wind speed and ustar (as initial guess in iteration)
    #cd10, u10 = cdnve(Ua, za) # Vera (1983)
    #Ustar = np.sqrt(cd10) * u10

    ## compute corrected 10 m wind speed
    #U10  = u10
    #U10o = 0
    #k    = 0

    #while np.abs( U10 - U10o ).max() > tol and  k < 15:
        #U10o = U10
        #k    = k + 1
        #Ustar = np.sqrt( cdnve(U10, 10) * U10**2)
        #U10   = Ua + Ustar * (np.log( zs / za ) - omegalmc(Xis) + omegalmc(Xia) ) / kappa


    #if k == 15:
        #print 'Iteration may not have converged'

    #delU = U10 - u10

    #return U10, delU

#def wave(Ua, za, Hw):
    #"""
    #Computes the true wind speed Ut at the measurement height za using the wind speed Ua measured at za and measured wave height Hw.

    #Parameters
    #----------
    #Ua : array_like
         #wind speed [m s :sup:`-1`]
    #za : float, array_like
         #wind measurement height [m]
    #Hw : float, array_like
         #wave height [m]

    #Returns
    #-------
    #Ut : array_like
          #true wind speed [m s :sup:`-1`]

    #Notes
    #-----
    #TODO

    #Examples
    #--------
    #TODO

    #References
    #----------
    #.. [1] Large, Morzel, and Crawford (1995), J. Phys. Oceanog., 25, 2959-2971.

    #Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    #05/05/1997: version 1.0
    #07/28/1999: version 1.1
    #08/05/1999: version 2.0
    #11/26/2010: Filipe Fernandes, Python translation.
    #"""
    ## convert input to numpy array
    #Ua = np.asarray(Ua) ; za = np.asarray(za) ; Hw = np.asarray(Hw)

    #k   = 0.4 #FIXME: kappa
    #z10 = 10
    #A = np.log( z10 / za ) / k

    ## eliminate any Ua==0
    #jj = np.nonzero(Ua==0)
    #Ua[jj] = 0.01 * np.ones( Ua[jj].shape )

    ## compute uncorrected 10 m wind speed
    #u10 = Ua # initial guess
    #for n in range(1, 11, 1):
        #ustar = np.sqrt( cdnve(u10) * u10**2 )
        #u10   = Ua + ustar * A

    ## compute corrected 10 m wind speed
    #Ustar = ustar
    #U10   = u10 #initial guesses
    #Za    = za / Hw
    #Z10   = z10 / Hw

    #for n in range(1, 11, 1):
        #Ustar = np.sqrt( cdnve(U10) * U10**2, za ) # FIXME: za was missing from the original, which was broken...
        #U10   = Ua + Ustar * ( np.log( z10 / za ) - omegalmc(Z10) + omegalmc(Za) ) / k

    ## compute 'true' wind speed at za using U10, Ustar
    #Ut = U10 - Ustar * A

    #return Ut

#def wavdist(za):
    #"""
    #FIXME: ?
    #Compute the wave distortion effects on wind at za for the following significant wave heights Hw = range(0,10,2) in meters.

    #Parameters
    #----------
    #za : array_like
         #wind measurement height [m]

    #Returns
    #-------
    #TODO

    #See Also
    #--------
    #TODO

    #Examples
    #--------
    #TODO

    #References
    #----------
    #TODO

    #Modifications: Original from AIR_SEA TOOLBOX, Version 2.0
    #05/05/1997: version 1.0
    #04/10/1998: version 1.1
    #08/05/1999: version 2.0
    #11/26/2010: Filipe Fernandes, Python translation.
    #"""
    ## convert input to numpy array
    #za = np.asarray(za)

    #Hw = np.arange(0,10,2)
    #Ua = np.arange(0.01, 20.01, 0.01)
    #N  = len(Hw)
    #M  = len(Ua)
    ##Ut = np.zeros((M, N))
    #Ut = []

    #for n in range(0, N+1, 1):
        ##Ut = wave( Ua, za, Hw[n] )
        #Ut.append( wave( Ua, za, Hw[n] ) )

    ##plot(Ua,Ut)
    ##title(['Predicted effects of wave distortion on wind speed at height ',num2str(za),' m'])
    ##xlabel('Measured wind speed Ua (m/s)')
    ##ylabel('Predicted wind speed Ut (m/s)')
    ##text(10,2,'Wave height increment = 2 m')

    #return Ua, Ut
