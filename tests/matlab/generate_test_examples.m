% generate_test_examples.m
% --------------------
%
% Generate example test output for Python airsea package.

% These examples come from the Python docstrings

ex = NaN(6,30);

% atmosphere
ex(1:5,1) = vapor([0.0, 5., 21., 33., 100.]);

ex(1:3,2) = air_dens([5., 15., 23.1], 95.);
ex(1:3,3) = air_dens([5., 15., 23.1], [95, 100, 50], 1000.);

ex(1:4,4) = qsat([10., 11.3, 15.5, 21.]);
ex(1:4,5) = qsat([10., 11.3, 15.5, 21.], 910.);

ex(1,6) = satvap(10.);
ex(1:3,7) = satvap([15., 25., 40.]);
ex(1:3,8) = satvap([15., 25., 40.], 900.);
ex(1:3,9) = satvap([15., 25., 40.],[900.,1000., 1030]);

ex(1:6,10) = rhadj([98., 95., 94., 93., 10., 20.], 98.);

Qlat = [550., 450., 350.];
rfd = [10., 15., 35.]; % precip rate (mm/min)

% note difference in rfd input units: Matlab: mm/min, Python: m/s
[E,P] = ep(rfd, Qlat);
ex(1:3,11) = E;
ex(1:3,12) = P;

ex(1:3,13) = relhumid([10., 15., 30.1], [9.2, 10.1, 33.0], 900);
ex(1:3,14) = relhumid([10., 15., 30.1], [9.2, 10.1, 33.0], 1020, 'assman');

ex(1:3,15) = cloudcor([0.8, 0.1, 1], 'clarke', 47);
ex(1,16) = cloudcor(0.5, [0.2, 0.4], 47);

ex(1:3,17:18) = viscair([[0.1, 5., 15];[22.8, 28.9, 31.4]])';

% windstress
ex(1:6,19) = cdnlp([10., 0.2, 12., 20., 30., 50.], 10);
ex(1:6,20) = cdnve([10., 0.2, 12., 20., 30., 50.], 15);
ex(1:6,21) = cdntc([10., 0.2, 12., 20., 30., 50.], 20, 20.);

[sp1,ustar1] = spshftlp([10., 0.2, 12., 20., 30., 50.], 10, 10);
ex(1:6,22) = sp1;
ex(1:6,23) = ustar1;

[sp2,ustar2] = spshfttc([10., 0.2, 12., 20., 30., 50.], 10, 8, 20);
ex(1:6,24) = sp2;
ex(1:6,25) = ustar2;

[sp3,ustar3] = spshftve([10., 0.2, 12., 20., 30., 50.], 15, 10);
ex(1:6,26) = sp3;
ex(1:6,27) = ustar3;

ex(1:6,28) = stresslp([10., 0.2, 12., 20., 30., 50.], 10);
ex(1:6,29) = stresstc([10., 0.2, 12., 20., 30., 50.], 15, 23., 1.02);
ex(1:6,30) = stressve([10., 0.2, 12., 20., 30., 50.], 8);

writematrix(ex,'test_examples_out.csv');
