function [ra, dec] = apstar (tjd, n, rai, deci, pmra, pmdec, parlax, radvel)

% apparent place of a star

%  tjd    = tt julian date for apparent place (in)

%  n      = body identification number for the earth (in) (no longer used)

%  rai    = icrs right ascension in hours (in)

%  deci   = icrs declination in degrees (in)

%  pmra   = icrs proper motion in ra in milliarcseconds/year (in)

%  pmdec  = icrs proper motion in dec in milliarcseconds/year (in)

%  parlax = parallax in milliarcseconds (in)

%  radvel = radial velocity in kilometers/second (in)

%  ra     = apparent right ascension in hours (out)

%  dec    = apparent declination in degrees (out)

% note: coordinate system for output ra and dec is equator and equinox of date

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

locatn = 0;

icoord = 1;

ttjd = tjd;

object = 'star';

star(1) = rai;

star(2) = deci;

star(3) = pmra;

star(4) = pmdec;

star(5) = parlax;

star(6) = radvel;

observ = zeros(6, 1);

skypos = place (ttjd, object, locatn, icoord, star, observ);

ra = skypos(4);

dec = skypos(5);