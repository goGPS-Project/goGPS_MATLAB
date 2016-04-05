function [pos, vel, ierr] = solsys (tjd, body, origin)

% purpose

% this is solsys version '2-da'. it is intended to provide
% an interface between the jpl direct-access solar system
% ephemerides and the 'novas' astrometric function library

% references

% standish, e.m. and newhall, x x (1988). "the jpl
% export planetary ephemeris"; jpl document dated 17 june 1988.

% kaplan, g.h. (1988). "novas: naval observatory vector astrometry
% subroutines"; usno document dated 20 october 1988

% input

%  tjd    = julian date of the desired time, on the tdb time scale

%  body   = body identification number for the
%           solar system object of interest;
%           mercury = 1,...,pluto = 9, sun = 10, moon = 11

%  origin = origin code; solar system barycenter = 0,
%           center of mass of the sun = 1

% output

%  pos  = position vector of 'body' at tjd;
%         equatorial rectangular coordinates in
%         au referred to the mean equator and
%         equinox of j2000.0

%  vel  = velocity vector of 'body' at tjd;
%         equatorial rectangular system referred
%         to the mean equator and equinox of
%         j2000.0, in au/day

%  ierr = 0 ... everything ok .
%       = 1 ... invalid value of 'body'
%       = 2 ... invalid value of 'origin'
%       = 3 ... error detected by jpl software;
%               tjd out of range

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

ierr = 0;

% perform sanity checks on the input body and origin

if (body < 1 || body > 11)
    
   ierr = 1;
   
   return;
   
elseif (origin < 0 || origin > 1)
    
   ierr = 2;
   
   return;
   
end

% select 'targ' according to value of 'body'

if (body == 10)
    
   targ = 11;
   
elseif (body == 11)
    
   targ = 10;
   
else
    
   targ = body;
   
end

% select 'cent' according to the value of 'origin'.

if (origin == 0)
    
   cent = 12;
   
elseif (origin == 1)
    
   cent = 11;
   
end

% call jpl routine to obtain position and velocity array 'posvel'

posvel = jplephem (tjd, targ, cent);

% decompose 'posvel' into position 'pos' and velocity 'vel'.

for i = 1:1:3
    
    pos(i) = posvel(i);
    
    vel(i) = posvel(i + 3);
    
end

