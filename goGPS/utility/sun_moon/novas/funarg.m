function [el, elp, f, d, omega] = funarg (t)

% this function computes fundamental arguments (mean elements)
% of the sun and moon.  see simon et al. (1994) astronomy and
% astrophysics 282, 663-683, especially sections 3.4-3.5.

% input

%  t = tdb time in julian centuries since j2000.0

% output

%  el = mean anomaly of the moon in radians at date tjd

%  elp = mean anomaly of the sun in radians at date tjd

%  f = mean longitude of the moon minus mean longitude of
%      the moon's ascending node in radians at date tjd

%  d = mean elongation of the moon from the sun in radians at date tjd

%  omega = mean longitude of the moon's ascending node in radians at date tjd

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.0d0 * 3600.0d0 / pi;

rev = 360.0d0 * 3600.0d0;

% fundamental (delaunay) arguments from simon et al. (1994)

% mean anomaly of the moon

el = rem (485868.249036d0 + ...
    t * (1717915923.2178d0 + ...
    t * (31.8792d0 + ...
    t * (0.051635d0 + ...
    t * (-0.00024470d0)))), rev) / seccon;

% mean anomaly of the sun

elp = rem (1287104.79305d0 + ...
    t * (129596581.0481d0 + ...
    t * (-0.5532d0 + ...
    t * (0.000136d0 + ...
    t * (- 0.00001149d0)))), rev) / seccon;

% mean argument of the latitude of the moon

f = rem (335779.526232d0 + ...
    t * (1739527262.8478d0 + ...
    t * (-12.7512d0 + ...
    t * (-0.001037d0 + ...
    t * (0.00000417d0)))), rev) / seccon;

% mean elongation of the moon from the sun

d = rem (1072260.70369d0 + ...
    t * (1602961601.2090d0 + ...
    t * (- 6.3706d0 + ...
    t * (0.006593d0 + ...
    t * (- 0.00003169d0)))), rev) / seccon;

% mean longitude of the ascending node of the moon (from simon
% section 3.4(b.3), precession = 5028.8200 arcsec/cy)

omega = rem (450160.398036d0 + ...
    t * (- 6962890.5431d0 + ...
    t * (7.4722d0 + ...
    t * (0.007702d0 + ...
    t * (- 0.00005939d0)))), rev) / seccon;


