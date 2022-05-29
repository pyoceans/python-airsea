import numpy as np
import airsea

test_out = np.genfromtxt('matlab/test_data_out.csv',delimiter=',')
test2_5b = np.genfromtxt('matlab/test2_5b.dat',skip_header=1)

ur  = test2_5b[:,1]
zr  = 15
Ta  = test2_5b[:,3]
zt  = 15
Pa  = 1008*np.ones(np.shape(ur))
q   = test2_5b[:,4]
zq  = 15
Ts  = test2_5b[:,13]
sal = 30*np.ones(np.shape(ur))
dsw = test2_5b[:,8]
nsw = 0.945*dsw
dlw = test2_5b[:,9]
rain= test2_5b[:,10]


rh  = q/airsea.atmosphere.qsat(Ta,Pa)/10

test = airsea.atmosphere.air_dens(Ta,rh)

