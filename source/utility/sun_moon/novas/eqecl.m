function [elon, elat] = eqecl (tjd, icoord, ra, dec)

% this function converts right ascension and declination
% to ecliptic longitude and latitude

% input

%  tjd    = tt julian date of equator, equinox, and ecliptic
%           used for coordinates

%  icoord = coordinate system selection

%               set icoord = 0 for mean equator and equinox of date
%               set icoord = 1 for true equator and equinox of date

%               (ecliptic is always the mean plane)

%  ra     = right ascension in hours, referred to specified
%           equator and equinox of date

%  dec    = declination in degrees, referred to specified
%           equator and equinox of date

% output

%  elon = ecliptic longitude in degrees, referred to specified
%         ecliptic and equinox of date

%  elat = ecliptic latitude in degrees, referred to specified
%         ecliptic and equinox of date

% note:  to convert icrs ra and dec to ecliptic coordinates (mean
% ecliptic and equinox of j2000.0), set tjd = 0.d0 and icoord = 0.
% except for the input to this case, all coordinates are dynamical.

% ported from NOVAS3.0

%%%%%%%%%%%%%%%%%%%%%%

radcon = pi / 180.0d0;

% form position vector in equatorial system from input coordinates

r = ra * 15.0d0 * radcon;

d = dec * radcon;

pos1(1) = cos(d) * cos(r);

pos1(2) = cos(d) * sin(r);

pos1(3) = sin(d);

% convert the vector from equatorial to ecliptic system

pos2 = eqec (tjd, icoord, pos1);

% decompose ecliptic vector into ecliptic longitude and latitude

xyproj = sqrt(pos2(1)^2 + pos2(2)^2);

e = 0.0d0;

if (xyproj > 0.0d0)
    e = atan2 (pos2(2), pos2(1));
end

elon = e / radcon;

if (elon < 0.0d0)
    elon = elon + 360.0d0;
end

e = atan2 (pos2(3), xyproj);

elat = e / radcon;

