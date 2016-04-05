function [oblm, oblt, eqeq, dpsi, deps] = etilt (tjd)

% this function computes quantities related to the orientation
% of the earth's rotation axis at julian date tjd.

%   tjd    = tdb julian date for orientation parameters (in)

%   oblm   = mean obliquity of the ecliptic in degrees at date tjd (out)

%   oblt   = true obliquity of the ecliptic in degrees at date tjd (out)

%   eqeq   = equation of the equinoxes in time seconds at date tjd (out)

%   dpsi   = nutation in longitude in arcseconds at date tjd (out)

%   deps   = nutation in obliquity in arcseconds at date tjd (out)

% note:  the equation of the equinoxes includes the complementary terms.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

global psicor epscor

seccon = 180.0d0 * 3600.0d0 / pi;

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

% get method/accuracy mode

mode = getmod;

t = (tjd - t0) / 36525.0d0;

% obtain nutation parameters in arcseconds

[psi_novas, eps] = nod (t);

% obtain complementary terms for equation of the equinoxes in arc seconds

if (mod (mode, 2) == 0)
    
    % high-accuracy mode
    
    cterms = eect2000 (tjd, 0.0d0) * seccon;
    
else
    
    % low-accuracy mode
    
    [el, elp, f, d, omega] = funarg (t);
    
    % series from iers conventions (2003), chapter 5,
    % table 5.2c, with some adjustments to coefficient values
    % copied from iers function eect2000, which has a more
    % complete series
    
    cterms = 2640.96d-6 * sin (omega) ...
        +   63.52d-6 * sin (2.0d0 * omega) ...
        +   11.75d-6 * sin (2.0d0 * f - 2.0d0 * d + 3.0d0 * omega) ...
        +   11.21d-6 * sin (2.0d0 * f - 2.0d0 * d +         omega) ...
        -    4.55d-6 * sin (2.0d0 * f - 2.0d0 * d + 2.0d0 * omega) ...
        +    2.02d-6 * sin (2.0d0 * f             + 3.0d0 * omega) ...
        +    1.98d-6 * sin (2.0d0 * f             +         omega) ...
        -    1.72d-6 * sin (3.0d0 * omega) ...
        -    0.87d-6 * t * sin (omega);
    
    % (terms smaller than 2 microarcseconds omitted)
    
end

delpsi = psi_novas + psicor;

deleps = eps + epscor;

% compute mean obliquity of the ecliptic in arcseconds

obm = obliq(t);

% compute true obliquity of the ecliptic in arcseconds

obt = obm + deleps;

% compute equation of the equinoxes in arcseconds

ee = delpsi * cos (obm / seccon) + cterms;

% convert to output units

oblm = obm / 3600.0d0;

oblt = obt / 3600.0d0;

eqeq = ee  / 15.0d0;

dpsi = delpsi;

deps = deleps;

