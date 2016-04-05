function [dpsi, deps] = nod (t)

% this function returns the values for nutation in longitude and
% nutation in obliquity for a given tdb julian date.

%  t     = tdb time in julian centuries since j2000.0 (in)

%  dpsi  = nutation in longitude in arcseconds (out)

%  deps  = nutation in obliquity in arcseconds (out)

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.0d0 * 3600.0d0 / pi;

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

% get method/accuracy mode

mode = getmod;

t1 = t * 36525.0d0;

% evaluate nutation series

% resulting nutation in longitude and obliquity in arc seconds

if (mod (mode, 2) == 0)
    
    % high accuracy mode -- iers 2000a
    
    [dp, de] = nut2000a (t0, t1);
    
else
    
    % low accuracy mode -- iau 2000k
    
    [dp, de] = nut2000k (t0, t1);
    
end

dpsi = dp * seccon;

deps = de * seccon;


