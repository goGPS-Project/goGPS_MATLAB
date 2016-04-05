function [pos, vel] = geopos (tjd, locatn, observ)

% this function computes the geocentric position and velocity
% of an observer on the surface of the earth or on a near-earth
% spacecraft. the final vectors are expressed in the gcrs.

%  tjd    = tt julian date (in)

%  locatn = integer code specifying location of observer (in)
%           set locatn=0 for observer at geocenter
%           set locatn=1 for observer on surface of earth
%           set locatn=2 for observer on near-earth spacecraft

%  observ = array of data specifying location of observer (in)
%           for locatn=0, this array not used

%           for locatn = 1,

%           observ(1) = geodetic longitude (wgs-84) of observer
%                       (east +) in degrees (in)
%           observ(2) = geodetic latitude (wgs-84) of observer
%                       (north +) in degrees (in)
%           observ(3) = height of observer above ellipsoid
%                       in meters (in)
%           observ(4) = value of delta-t in seconds (in)
%                       (delta-t=tt-ut1)
%           observ(5) = (not used, reserved for future use)
%           observ(6) = (not used, reserved for future use)

%           for locatn = 2,

%           observ(1) = geocentric x in kilometers
%           observ(2) = geocentric y in kilometers
%           observ(3) = geocentric z in kilometers
%           observ(4) = geocentric x-dot in kilometers/second
%           observ(5) = geocentric y-dot in kilometers/second
%           observ(6) = geocentric z-dot in kilometers/second
%           with respect to true equator and equinox of date

%  pos    = position vector of observer, with respect to origin
%           at geocenter, referred to gcrs axes, components in au (out)

%  vel    = velocity vector of observer, with respect to origin
%           at geocenter, referred to gcrs axes, components in au/day (out)

% note 1: if locatn = 1 and observ(4) = 0.0d0, the value of delta-t will
% be obtained from getdt, which provides the last value of delta-t
% defined by user via call to setdt.

% note 2: this subroutine calls subroutine terra for an observer
% on the surface of the earth. terra neglects polar motion, an
% approximation which may yield up to 15 meters error in position
% and several millimeters/sec error in velocity.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

% get au, the length of the astronomical unit in kilometers

au = astcon ('au', 1.0d-3);

if (locatn == 0)
    
    pos(1) = 0.0d0;
    pos(2) = 0.0d0;
    pos(3) = 0.0d0;
    vel(1) = 0.0d0;
    vel(2) = 0.0d0;
    vel(3) = 0.0d0;
    
    return
    
end

ttjd  = tjd;

% tdb is approximated by tt

tdbjd = tjd;

% get geocentric position and velocity vectors of observer wrt
% equator and equinox of date

if (locatn == 1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % observer on surface of earth
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get delta-t value
    
    if (observ(4) ~= 0.0d0)
        
        deltat = observ(4) / 86400.0d0;
        
    else
        
        deltat = getdt;
        
    end
    
    % using delta-t value, compute ut1 and sidereal time
    
    if (ttjd == 0.0d0)
        
        ut1jd = tdbjd - deltat;
        
    else
        
        ut1jd = ttjd - deltat;
        
    end
    
    eqinox;
    
    gmst = sidtim (ut1jd, 0.0d0, 0);
    
    [x, x, eqeq, x, x] = etilt (tdbjd);
    
    resume;
    
    gast = gmst + eqeq / 3600.0d0;
    
    % subroutine terra does the hard work, given sidereal time
    
    [pos1, vel1] = terra (observ(1), observ(2), observ(3), gast);
    
elseif (locatn == 2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % observer on near-earth spacecraft
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % convert units to au and au/day
    
    for j = 1:3
        
        pos1(j) = observ(j) / au;
        
        vel1(j) = observ(j + 3) / au * 86400.0d0;
        
    end
    
end

% transform geocentric position vector of observer to gcrs

pos2 = nutate (-tdbjd, pos1);

pos3 = preces (tdbjd, pos2, t0);

pos = frame (pos3, -1);

% transform geocentric velocity vector of observer to gcrs

vel2 = nutate (-tdbjd, vel1);

vel3 = preces (tdbjd, vel2, t0);

vel = frame (vel3, -1);



