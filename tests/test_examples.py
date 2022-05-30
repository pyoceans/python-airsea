# -*- coding: utf-8 -*-

import numpy as np
from airsea.atmosphere import *
from airsea.windstress import *

# paths relative to package root for running pytest
ex_out = np.genfromtxt('tests/matlab/test_examples_out.csv', delimiter=',')

# atmosphere
def test_vapor():
    ex0 = vapor([0.0, 5., 21., 33., 100.])
    ii = np.isfinite(ex_out[:, 0])
    assert np.allclose(ex0, ex_out[ii, 0])


def test_airdens():
    ex1 = air_dens([5., 15., 23.1], 95.)
    ex2 = air_dens([5., 15., 23.1], [95, 100, 50], 1000.)
    ii1 = np.isfinite(ex_out[:, 1])
    ii2 = np.isfinite(ex_out[:, 2])
    assert (np.allclose(ex1, ex_out[ii1, 1]) &
            np.allclose(ex2, ex_out[ii2, 2]))


def test_qsat():
    ex3 = qsat([10., 11.3, 15.5, 21.])
    ex4 = qsat([10., 11.3, 15.5, 21.], 910.)
    ii3 = np.isfinite(ex_out[:, 3])
    ii4 = np.isfinite(ex_out[:, 4])
    assert (np.allclose(ex3, ex_out[ii3, 3]) &
            np.allclose(ex4, ex_out[ii4, 4]))


def test_satvap():
    ex5 = satvap(10.)
    ex6 = satvap([15., 25., 40.])
    ex7 = satvap([15., 25., 40.], 900.)
    ex8 = satvap([15., 25., 40.], [900., 1000., 1030])
    ii5 = np.isfinite(ex_out[:, 5])
    ii6 = np.isfinite(ex_out[:, 6])
    ii7 = np.isfinite(ex_out[:, 7])
    ii8 = np.isfinite(ex_out[:, 8])
    assert (np.allclose(ex5, ex_out[ii5, 5]) &
            np.allclose(ex6, ex_out[ii6, 6]) &
            np.allclose(ex7, ex_out[ii7, 7]) &
            np.allclose(ex8, ex_out[ii8, 8]))


def test_rhadj():
    ex9 = rhadj([98., 95., 94., 93., 10., 20.], 98.)
    ii9 = np.isfinite(ex_out[:, 9])
    assert np.allclose(ex9, ex_out[ii9, 9])


def test_rhd():
    Qlat = np.array([550., 450., 350.]) # latent heat flux
    rfd = np.array([10., 15., 35.]) # precip rate (mm/min)

    # note difference in rfd input units: Matlab: mm/min, Python: m/s
    ex10,ex11 = ep(rfd/1000/60, Qlat) 
    ii10 = np.isfinite(ex_out[:, 10])
    ii11 = np.isfinite(ex_out[:, 11])
    assert (np.allclose(ex10, ex_out[ii10, 10]) &
            np.allclose(ex11, ex_out[ii11, 11]))


def test_relhumid():
    ex12 = relhumid([10., 15., 30.1], [9.2, 10.1, 33.0], 900);
    ex13 = relhumid([10., 15., 30.1], [9.2, 10.1, 33.0], 1020, 'assman');
    ii12 = np.isfinite(ex_out[:, 12])
    ii13 = np.isfinite(ex_out[:, 13])
    assert (np.allclose(ex12, ex_out[ii12, 12]) &
            np.allclose(ex13, ex_out[ii13, 13]))


def test_cloudcor():
    ex14 = cloudcor([0.8, 0.1, 1], 'clarke', 47);
    ex15 = cloudcor(0.5, [0.2, 0.4], 47);
    ii14 = np.isfinite(ex_out[:, 14])
    ii15 = np.isfinite(ex_out[:, 15])
    assert (np.allclose(ex14, ex_out[ii14, 14]) &
            np.allclose(ex15, ex_out[ii15, 15]))


def test_viscair():
    visc = visc_air([[0.1, 5., 15],[22.8, 28.9, 31.4]])
    ex16 = visc[0,:]
    ex17 = visc[1,:]
    ii = np.isfinite(ex_out[:, 16])
    assert (np.allclose(ex16, ex_out[ii, 16]) &
            np.allclose(ex17, ex_out[ii, 17]))


# windstress
def test_cdn():
    ex18,u10 = cdn([10., 0.2, 12., 20., 30., 50.], 10)
    ex19,u10 = cdn([10., 0.2, 12., 20., 30., 50.], 15, drag='vera')
    ex20,u10 = cdn([10., 0.2, 12., 20., 30., 50.], 20, 'smith', 20.)

    assert (np.allclose(ex18, ex_out[:, 18]) &
            np.allclose(ex19, ex_out[:, 19]) &
            np.allclose(ex20, ex_out[:, 20]))


def test_spshft():
    sp1, ustar1 = spshft([10., 0.2, 12., 20., 30., 50.], 10, 10)
    ex21 = sp1
    ex22 = ustar1

    sp2, ustar2 = spshft([10., 0.2, 12., 20., 30., 50.], 10, 8, 'smith', 20)
    ex23 = sp2
    ex24 = ustar2

    sp3, ustar3 = spshft([10., 0.2, 12., 20., 30., 50.], 15, 10, 'vera')
    ex25 = sp3
    ex26 = ustar3

    assert (np.allclose(ex21, ex_out[:, 21]) &
            np.allclose(ex22, ex_out[:, 22]) &
            np.allclose(ex23, ex_out[:, 23]) &
            np.allclose(ex24, ex_out[:, 24]) &
            np.allclose(ex25, ex_out[:, 25]) &
            np.allclose(ex26, ex_out[:, 26]))


# Note that Matlab and Python versions have different tolerance thresholds 
# for testing convergence of u10

def test_stress():
    ex27 = stress([10., 0.2, 12., 20., 30., 50.], 10)
    kw = dict(rho_air=1.02, Ta=23.)
    ex28 =  stress([10., 0.2, 12., 20., 30., 50.], 15, 'smith', **kw)
    ex29 = stress([10., 0.2, 12., 20., 30., 50.], 8, 'vera')

    print(ex29-ex_out[:, 29])

    assert (np.allclose(ex27, ex_out[:, 27]) &
            np.allclose(ex28, ex_out[:, 28]) &
            np.allclose(ex29, ex_out[:, 29], atol=1e-4, rtol=1e-4))


