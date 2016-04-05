function raeq = eqxra (tjd, k)

% this function computes the intermediate right ascension
% of the equinox at julian date tjd, using an analytical expression
% for the accumulated precession in right ascension.  for the
% true equinox the result is the equation of the origins.

% input

%  tjd = tdb julian date

%  k = equinox selection code

%      set k = 0 for mean equinox
%      set k = 1 for true equinox (equation of the origins)

% output

%  raeq = intermediate right ascension of the equinox, in hours (+ or -)

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

t = (tjd - t0) / 36525.0d0;

% for the true equinox, obtain the equation of the equinoxes in time seconds

if (k == 1)
    
    [a, a, ee, a, a] = etilt (tjd);
    
    eqeq = ee;
    
else
    
    eqeq = 0.0d0;
    
end

% precession in ra in arcseconds taken from capitaine et al. (2003),
% astronomy and astrophysics 412, 567-586, eq. (42)

precra = 0.014506d0 + ...
    (((( -0.0000000368d0 * t ...
    - 0.000029956d0) * t ...
    - 0.00000044d0) * t ...
    + 1.3915817d0) * t ...
    + 4612.156534d0) * t;

raeq = -(precra / 15.0d0 + eqeq) / 3600.0d0;


