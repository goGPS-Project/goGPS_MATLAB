function jdtdb = utc2tdb (jdutc)

% convert UTC julian date to TDB julian date

% input

%  jdutc   = UTC julian date
%  tai_utc = TAI-UTC (seconds)

% output

%  jdtdb = TDB julian date 

% note: requires readleap in main script

% Reference Frames in Astronomy and Geophysics
% J. Kovalevsky et al., 1989, pp. 439-442

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtr = pi / 180.0;

% compute TDT julian date
 
leapsecond = findleap(jdutc);

corr = (32.184 + leapsecond) / 86400.0;
    
jdtdt = jdutc + corr;
    
% time argument for correction

t = (jdtdt - 2451545.0) / 36525.0;
    
% compute correction in microseconds

corr = 1656.675     * sin(dtr * (35999.3729 * t + 357.5287))...
       + 22.418     * sin(dtr * (32964.467  * t + 246.199))...
       + 13.84      * sin(dtr * (71998.746  * t + 355.057))...
       +  4.77      * sin(dtr * ( 3034.906  * t +  25.463))...
       +  4.677     * sin(dtr * (34777.259  * t + 230.394))...
       + 10.216 * t * sin(dtr * (35999.373  * t + 243.451))...
       +  0.171 * t * sin(dtr * (71998.746  * t + 240.98 ))...
       +  0.027 * t * sin(dtr * ( 1222.114  * t + 194.661))...
       +  0.027 * t * sin(dtr * ( 3034.906  * t + 336.061))...
       +  0.026 * t * sin(dtr * (  -20.186  * t +   9.382))...
       +  0.007 * t * sin(dtr * (29929.562  * t + 264.911))...
       +  0.006 * t * sin(dtr * (  150.678  * t +  59.775))...
       +  0.005 * t * sin(dtr * ( 9037.513  * t + 256.025))...
       +  0.043 * t * sin(dtr * (35999.373  * t + 151.121));

% convert corrections to days

corr = 0.000001 * corr / 86400.0;

% compute TDB julian date

jdtdb = jdtdt + corr;


