function [ra, dec, dis] = applan (tjd, l, n)

% geocentric place of the sun, moon and planets

%  tjd = tt julian date for apparent place (in)
%  l   = body identification number for desired object (in)
%  n   = body identification number for the earth (in) (no longer used)
%  ra  = apparent right ascension in hours (out)
%  dec = apparent declination in degrees (out)
%  dis = true distance from earth to planet in au (out)

% note: coordinate system for output ra and dec is equator and equinox of date

% note: 'planet' is used generically for any solar system body.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

locatn = 0;

icoord = 1;

ttjd = tjd;

object = horzcat('=', num2str(l));

star = zeros(6,1);

observ = zeros(6,1);

skypos = place (ttjd, object, locatn, icoord, star, observ);

ra = skypos(4);

dec = skypos(5);

dis = skypos(6);


