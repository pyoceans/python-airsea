# -*- coding: utf-8 -*-

import numpy as np
from airsea.atmosphere import air_dens, qsat, satvap, vapor, visc_air
from airsea.windstress import stress

# paths relative to package root for running pytest
test_out = np.genfromtxt('tests/matlab/test_data_out.csv', delimiter=',')
test2_5b = np.genfromtxt('tests/matlab/test2_5b.dat', skip_header=1)

ur = test2_5b[:, 1]
zr = 15
Ta = test2_5b[:, 3]
zt = 15
Pa = 1008*np.ones(np.shape(ur))
q = test2_5b[:, 4]
zq = 15
Ts = test2_5b[:, 13]
sal = 30*np.ones(np.shape(ur))
dsw = test2_5b[:, 8]
nsw = 0.945*dsw
dlw = test2_5b[:, 9]
rain = test2_5b[:, 10]

rh  = q / qsat(Ta, Pa) / 10

# Test atmosphere functions


def test_airdens():
    test0 = air_dens(Ta, rh)
    test1 = air_dens(Ta, rh, Pa)
    assert (np.allclose(test0, test_out[:, 0]) &
            np.allclose(test1, test_out[:, 1]))


def test_qsat():
    test2 = qsat(Ta, Pa)
    assert np.allclose(test2, test_out[:, 2])


def test_satvap():
    test3 = satvap(Ta)
    test4 = satvap(Ta, Pa)
    assert (np.allclose(test3, test_out[:, 3]) &
            np.allclose(test4, test_out[:, 4]))


def test_vapor():
    test5 = vapor(Ts)
    assert np.allclose(test5, test_out[:, 5])


def test_viscair():
    test6 = visc_air(Ta)
    assert np.allclose(test6, test_out[:, 6])


# Test windstress functions

# Note that Matlab and Python versions have different tolerance thresholds 
# for testing convergence of u10


def test_stresslp():
    test7 = stress(ur, zr)
    test8 = stress(ur, zr, rho_air=test_out[:, 0])
    assert (np.allclose(test7, test_out[:, 7], atol=1e-6) &
            np.allclose(test8, test_out[:, 8], atol=1e-6))


def test_stresstc():
    test9 = stress(ur, zr, drag='smith')
    test10 = stress(ur, zr, Ta=Ta, drag='smith')
    test11 = stress(ur, zr, Ta=Ta, rho_air=np.mean(test_out[:, 0]), drag='smith')
    assert (np.allclose(test9, test_out[:, 9], atol=1e-6) &
            np.allclose(test10, test_out[:, 10], atol=1e-6) &
            np.allclose(test11, test_out[:, 11], atol=1e-6))


def test_stressve():
    test12 = stress(ur, zr, drag='vera')
    test13 = stress(ur, zr, rho_air=np.mean(test_out[:, 0]), drag='vera')
    assert (np.allclose(test12, test_out[:, 12], atol=1e-6) &
            np.allclose(test13, test_out[:, 13], atol=1e-6))
