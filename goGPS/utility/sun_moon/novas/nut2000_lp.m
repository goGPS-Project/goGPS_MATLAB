function [dpsi, deps] = nut2000_lp (t)

% low precison nutation based on iau 2000a

% this function evaluates a short nutation series and returns approximate
% values for nutation in longitude and nutation in obliquity for a given
% tdb julian date. in this mode, only the largest 13 terms of the iau 2000a
% nutation series are evaluated.

% input

%  t = tdb time in julian centuries since j2000.0

% output

%  dpsi = nutation in longitude in arcseconds

%  deps = nutation in obliquity in arcseconds

% note: in low-accuracy mode, max error in dpsi < 0.05 arcsec,
% max error in deps < 0.02 arcsec, average error about 1/4 of max.

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

% largest 13 terms of iau 2000a nutation series, with precision
% of coefficients truncated

x = [0.0, 0.0, 0.0, 0.0, 1.0, -17.2064,-0.01747,  9.2052,  0.00091; ...
     0.0, 0.0, 2.0,-2.0, 2.0, -1.3171, -0.00017,  0.5730, -0.00030; ...
     0.0, 0.0, 2.0, 0.0, 2.0, -0.2276, -0.00002,  0.0978, -0.00005; ...
     0.0, 0.0, 0.0, 0.0, 2.0,  0.2075,  0.00002, -0.0897,  0.00005; ...
     0.0, 1.0, 0.0, 0.0, 0.0,  0.1476, -0.00036,  0.0074, -0.00002; ...
     0.0, 1.0, 2.0,-2.0, 2.0, -0.0517,  0.00012,  0.0224, -0.00007; ...
     1.0, 0.0, 0.0, 0.0, 0.0,  0.0711,  0.00001, -0.0007,  0.00000; ...
     0.0, 0.0, 2.0, 0.0, 1.0, -0.0387, -0.00004,  0.0201,  0.00000; ...
     1.0, 0.0, 2.0, 0.0, 2.0, -0.0301,  0.00000,  0.0129, -0.00001; ...
     0.0,-1.0, 2.0,-2.0, 2.0,  0.0216, -0.00005, -0.0096,  0.00003; ...
     0.0, 0.0, 2.0,-2.0, 1.0,  0.0128,  0.00001, -0.0069, -0.00000; ...
    -1.0, 0.0, 2.0, 0.0, 2.0,  0.0123,  0.00000, -0.0053,  0.00000; ...
    -1.0, 0.0, 0.0, 2.0, 0.0,  0.0157,  0.00000, -0.0001,  0.00000];

% transpose x matrix

x = x';

% ----------------------------------------------------
% remaining terms all have amplitudes < 0.01 arcsecond
% ----------------------------------------------------

% computation of fundamental arguments

[el, elp, f, d, om] = funarg (t);

dpsi = 0.0d0;

deps = 0.0d0;

% sum nutation series terms

for i = 13:-1:1
    
    arg = x(1, i) * el ...
        + x(2, i) * elp ...
        + x(3, i) * f ...
        + x(4, i) * d ...
        + x(5, i) * om;
    
    dpsi = (x(6, i) + x(7, i) * t) * sin(arg) + dpsi;
    
    deps = (x(8, i) + x(9, i) * t) * cos(arg) + deps;
    
end

% add in out-of-phase component of principal (18.6-year) term
% (to avoid small but long-term bias in results)

dpsi = dpsi + 0.0033d0 * cos(om);

deps = deps + 0.0015d0 * sin(om);


