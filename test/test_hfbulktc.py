"""
Variables:
  ur     = wind speed [m/s] measured at height zr [m]
  Ta     = air temperature [C] measured at height zt [m]
  rh     = relative humidity [%] measured at height zq [m]
  Pa     = air pressure [mb]
  Ts     = sea surface temperature [C]
  sal    = salinity [psu (PSS-78)]
  dlw    = downwelling (INTO water) longwave radiation [W/m^2]
  dsw    = measured insolation [W/m^2]
  nsw    = net shortwave radiation INTO the water [W/m^2]
  rain   = rain rate  [mm/hr]

04/07/1999: version 1.2 (contributed by AA)
08/05/1999: version 2.0
"""

import os
import unittest
import numpy as np
from oct2py import octave

from airsea.atmosphere import qsat

rootpath = os.path.dirname(__file__)
path = os.path.join(rootpath, 'air_sea')
_ = octave.addpath(octave.genpath(path))

test2_5b = np.loadtxt(os.path.join(rootpath, 'air_sea', 'test2_5b.dat'), comments='%')

ur = test2_5b[:, 1]
zr = 15.
Ta = test2_5b[:, 3]
zt = 15.
Pa = 1008 * np.ones_like(ur)
q = test2_5b[:, 4]
rh = q / qsat(Ta, Pa) / 10.
zq = 15.
Ts = test2_5b[:, 13]
sal = 30 * np.ones_like(ur)
dsw = test2_5b[:, 8]
nsw = 0.945 * dsw
dlw = test2_5b[:, 9]
rain = test2_5b[:, 10]


# No cool-skin; compare to COARE output in file no_skin.out.
A1 = hfbulktc(ur, zr, Ta, zt, rh, zq, Pa, Ts)


# YES cool-skin; compare to COARE output file yes_skin.out.
A2 = hfbulktc(ur, zr, Ta, zt, rh, zq, Pa, Ts, sal, dlw, dsw, nsw)
