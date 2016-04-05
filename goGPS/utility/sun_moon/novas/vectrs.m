function [pos, vel] = vectrs (ra, dec, pmra, pmdec, parllx, rv)

% this function converts angular quantities related to a star's
% position and motion to vectors.

% input

%  ra     = right ascension in hours

%  dec    = declination in degrees

%  pmra   = proper motion in ra in milliarcseconds per year
               
%  pmdec  = proper motion in dec in milliarcseconds per year

%  parllx = parallax in milliarcseconds

%  rv     = radial velocity in kilometers/second

% output

%  pos = position vector, equatorial rectangular coordinates,
%           with respect to solar system barycenter, components in au

%  vel = velocity vector, equatorial rectangular coordinates,
%        with respect to solar system barycenter, components in au/day

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.0d0 * 3600.0d0 / pi;

% speed of light in kilometers/second

c = 1.0d-3 * 2997924580.0d0;

% au in kilometers

aukm = 1.0d-3 * 499.0047838061d0 * c;

% if parallax is unknown, undetermined, or zero, set it to 1e-6
% milliarcsecond, corresponding to a distance of 1 gigaparsec

paralx = parllx;

if ( paralx <= 0.0d0 )
    paralx = 1.0d-6;
end

% convert right ascension, declination, and parallax to position
% vector in equatorial system with units of au

dist = 1.0d0 / sin (paralx * 1.0d-3 / seccon);

r = ra * 54000.0d0 / seccon;

d = dec * 3600.0d0 / seccon;

cra = cos(r);

sra = sin(r);

cdc = cos(d);

sdc = sin(d);

pos(1) = dist * cdc * cra;

pos(2) = dist * cdc * sra;

pos(3) = dist * sdc;

% compute doppler factor, which accounts for change in
% light travel time to star

k = 1.d0 / (1.0d0 - rv / c);

% convert proper motion and radial velocity to orthogonal components
% of motion with units of au/day

pmr = pmra  / (paralx * 365.25d0) * k;

pmd = pmdec / (paralx * 365.25d0) * k;

rvl = rv * 86400.0d0 / aukm       * k;

% transform motion vector to equatorial system

vel(1) = - pmr * sra - pmd * sdc * cra + rvl * cdc * cra;

vel(2) =   pmr * cra - pmd * sdc * sra + rvl * cdc * sra;

vel(3) =               pmd * cdc       + rvl * sdc;

