function [ttjd, secdif] = times_novas (tdbjd)

% this function computes the terrestrial time (tt) julian date
% corresponding to a barycentric dynamical time (tdb) julian date.
% the expression used in this version is a truncated form of a
% longer and more precise series given by fairhead & bretagnon
% (1990) a&a 229, 240. the result is good to about 10 microseconds.

% input

%  tdbjd = tdb julian date

% output

%  ttjd = tt julian date

%  secdif = difference tdbjd - ttjd, in seconds

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

t = (tdbjd - t0) / 36525.0d0;

% expression given in usno circular 179, eq. 2.6

secdif = 0.001657d0 * sin(628.3076d0 * t + 6.2401d0) ...
    + 0.000022d0 * sin(575.3385d0 * t + 4.2970d0) ...
    + 0.000014d0 * sin(1256.6152d0 * t + 6.1969d0) ...
    + 0.000005d0 * sin(606.9777d0 * t + 4.0212d0) ...
    + 0.000005d0 * sin(52.9691d0 * t + 0.4444d0) ...
    + 0.000002d0 * sin(21.3299d0 * t + 5.5431d0) ...
    + 0.000010d0 * t * sin(628.3076d0 * t + 4.2490d0);

ttjd = tdbjd - secdif / 86400.0d0;
