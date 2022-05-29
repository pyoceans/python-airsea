
% generate_test_data.m
% --------------------
%
% Generate test data for Python airsea package
%
% Uses test2_5b.dat file included in Matlab airsea v2.0. Makes calculations and stores
% the results in test_data_out.csv

% VARIBLES:
%
%   ur     = wind speed [m/s] measured at height zr [m] 
%   Ta     = air temperature [C] measured at height zt [m]
%   rh     = relative humidity [%] measured at height zq [m]
%   Pa     = air pressure [mb]
%   Ts     = sea surface temperature [C]
%   sal    = salinity [psu (PSS-78)]
%   dlw    = downwelling (INTO water) longwave radiation [W/m^2]
%   dsw    = measured insolation [W/m^2]
%   nsw    = net shortwave radiation INTO the water [W/m^2]
%   rain   = rain rate  [mm/hr]


% directory containing Matlab airsea toolbox
dir_airsea = '../../../../airsea-matlab/';

addpath(dir_airsea)

test_data_path = [dir_airsea 'test2_5b.dat'];
load(test_data_path);

% parse data (following t_hfbulktc.m)
ur  = test2_5b(:,2);
zr  = 15;
Ta  = test2_5b(:,4);
zt  = 15;
Pa  = 1008*ones(size(ur));
q   = test2_5b(:,5);
rh  = q./qsat(Ta,Pa)/10;
zq  = 15;
Ts  = test2_5b(:,14);
sal = 30*ones(size(ur));
dsw = test2_5b(:,9);
nsw = 0.945*dsw;
dlw = test2_5b(:,10);
rain= test2_5b(:,11);

ntests = 14;
out = NaN([length(ur),ntests]);

% *** air_dens ***
out(:,1) = air_dens(Ta,rh);
out(:,2) = air_dens(Ta,rh,Pa);

% *** qsat ***
out(:,3) = qsat(Ta,Pa);

% *** satvap ***
out(:,4) = satvap(Ta);
out(:,5) = satvap(Ta,Pa);

% *** vapor ***
out(:,6) = vapor(Ts);

% *** viscair ***
out(:,7) = viscair(Ta);

% *** stress - Large and Pond ***
out(:,8) = stresslp(ur,zr);
out(:,9) = stresslp(ur,zr,out(:,1));

% *** stress - Smith (TC) ***
out(:,10) = stresstc(ur,zr);
out(:,11) = stresstc(ur,zr,Ta);
out(:,12) = stresstc(ur,zr,Ta,mean(out(:,1)));

% *** stress - Vera ***
out(:,13) = stressve(ur,zr);
out(:,14) = stressve(ur,zr,mean(out(:,1)));

writematrix(out,'test_data_out.csv');
