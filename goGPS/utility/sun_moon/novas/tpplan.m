function [ra, dec, dis] = tpplan (ujd, ttjd, object, glon, glat, ht)

% topocentric place of the sun, moon or planet

%  ujd    = ut1 julian date for topocentric place (in)
%  ttjd   = tt julian date (in)
%  object = target object id number (in)
%  glon   = geodetic (itrs) longitude (east +) of observer in degrees (in)
%  glat   = geodetic (itrs) latitude (north +) of observer in degrees (in)
%  ht     = height of observer in meters (in)
%  ra     = topocentric right ascension in hours (out)
%  dec    = topocentric declination in degrees (out)
%  dis    = true distance from observer to planet in au (out)

% note 1: coordinate system for output ra and dec is equator and equinox of date.

% note 2: ujd can also be greenwich apparent sidereal time in hours,
% equivalent to ut1 julian date, but this option will not be
% supported indefinitely. advise using ujd = ut1 julian date only.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

locatn = 1;

icoord = 1;

if (ujd > 100.0d0)
    
    deltat = (ttjd - ujd) * 86400.0d0;
    
else
    
    gast = mod (ujd, 24.0d0);
    
    if (gast < 0.0d0)
        
        gast = gast + 24.0d0;
        
    end
    
    deltat = 0.0d0;
    
end

observ(1) = glon;

observ(2) = glat;

observ(3) = ht;

observ(4) = deltat;

star = zeros(6, 1);

object = horzcat('=', num2str(object));

skypos = place (ttjd, object, locatn, icoord, star, observ);

ra = skypos(4);

dec = skypos(5);

dis = skypos(6);


