function rv = radvl (pos, vel, velobs, star, dist)

% this function predicts the radial velocity of the observed
% object as it would be measured by spectroscopic means.  radial
% velocity is here defined as the radial velocity measure (z)
% times the speed of light.  for a solar system body, it applies
% to a fictitious emitter at the center of the observed object,
% assumed massless (no gravitational red shift), and does not
% in general apply to reflected light.  for stars, it includes
% all effects, such as gravitational red shift, contained
% in the catalog barycentric radial velocity measure, a scalar
% derived from spectroscopy.  nearby stars with a known kinematic
% velocity vector (obtained independently of spectroscopy) can be
% treated like solar system objects.  see lindegren & dravins
% (2003), astronomy & astrophysics 401, 1185-1201.
%
%      pos    = geometric position vector of object with respect to
%               observer, corrected for light-time, in au (in)
%      vel    = velocity vector of object with respect to solar
%               system barycenter, components in au/day (in)
%      velobs = velocity vector of observer with respect to solar
%               system barycenter, components in au/day (in)
%      star   = 3-element array of catalog data for a star, to be
%               non-zero if observed object is a star for which the
%               catalog radial velocity is consistent with
%               the iau definition of barycentric radial velocity
%               measure (otherwise all elements should be set to
%               0.d0 exactly) (in)
%               star(1) = catalog ra in hours
%               star(2) = catalog dec in degrees
%               star(3) = z*c, the catalog barycentric radial
%                         velocity measure times the speed of light,
%                         in kilometers/second
%               all three data elements must apply to the same
%               epoch (usually j2000.0 = jd 2451545.0 tt)
%      dist   = 3-element array of distances in au (in)
%               dist(1) = distance of observer from the geocenter
%               dist(2) = distance of observer from the sun
%               dist(3) = distance of object from the sun
%      rv     = the observed radial velocity measure times
%               the speed of light, in kilometers/second (out)

% note 1:  all the input arguments are bcrs quantities, expressed
% with respect to the icrs axes.  vel and velobs are kinematic
% velocities -- derived from geometry or dynamics, not spectroscopy.

% note 2:  if any element of array star is non-zero, the algorithm
% used will be consistent with the iau definition of stellar
% radial velocity, specifically, the barycentric radial velocity
% measure, which is derived from spectroscopy.  in that case,
% the vector vel can be very approximate -- or, for distant stars
% or galaxies, zero -- as it will be used only for a small geometric
% correction that is proportional to proper motion.

% note 3:  any of the distances in array dist can be set to zero
% (0.d0) if the corresponding general relativistic gravitational
% potential term is not to be evaluated.  these terms
% generally are important only at the meter/second level.  if
% the first two distances are both zero, an average value
% will be used for the relativistic term for the observer,
% appropriate for an observer on the surface of the earth.  the
% third distance, if given, is used only for solar system objects.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

radcon = pi / 180.0d0;

% get au, length of astronomical unit in meters

au = astcon ('au', 1.0d0);

% get c, the speed of light in meters/second

c = astcon ('c', 1.0d0);

% get gs, heliocentric gravitational constant

gs = astcon ('gs', 1.0d0);

% get ge, geocentric gravitational constant

ge = astcon ('ge', 1.0d0);

% (gs and ge are in meters^3/second^2)

c2 = c^2;

toms = au / 86400.0d0;

toms2 = toms^2;

rv = 0.0d0;

% compute length of position vector = distance to object in au

posmag = sqrt(pos(1)^2 + pos(2)^2 + pos(3)^2);

if (posmag < 1.0d-8)
    
    return
    
end

% determine how object is to be processed

dostar = star(1) ~= 0.0d0 || star(2) ~= 0.0d0 || star(3) ~= 0.0d0;

% compute unit vector toward object

for j = 1:3
    
    uk(j) = pos(j) / posmag;
    
end

% compute velocity-squared factors

v2 = (vel(1)^2 + vel(2)^2 + vel(3)^2) * toms2;

vo2 = (velobs(1)^2 + velobs(2)^2 + velobs(3)^2) * toms2;

% compute geopotential at observer, unless observer is geocentric

r = dist(1) * au;

phigeo = 0.0d0;

if (r > 1.0d6)
    
    phigeo = ge / r;
    
end

% compute solar potential at observer

r = dist(2) * au;

phisun = 0.0d0;

if (r > 1.0d8)
    
    phisun = gs / r;
    
end

% compute relativistic potential and velocity factor for observer

if (dist(1) ~= 0.0d0 || dist(2) ~= 0.0d0)
    
    % lindegren & dravins eq. (41), second factor in parentheses
    
    rel = 1.0d0 - (phigeo + phisun) / c2 - 0.5d0 * vo2 / c2;
    
else
    
    % lindegren & dravins eq. (42), inverse
    
    rel = 1.0d0 - 1.550d-8;
    
end

if (dostar == 1)
    
    % for stars, update barycentric radial velocity measure for change in view angle
    
    ra = star(1) * 15.0d0 * radcon;
    
    dc = star(2) * radcon;
    
    du(1) = uk(1) - (cos (dc) * cos (ra));
    
    du(2) = uk(2) - (cos (dc) * sin (ra));
    
    du(3) = uk(3) - (sin (dc));
    
    zc = star(3) * 1.0d3 + (vel(1) * du(1) + vel(2) * du(2) + vel(3) * du(3)) * toms;
    
    % compute observed radial velocity measure of a star (inverse of
    % lindegren & dravins eq. (41))
    
    zb1 = 1.d0 + zc / c;
    
    kvobs = (uk(1) * velobs(1) + uk(2) * velobs(2) + uk(3) * velobs(3)) * toms;
    
    zobs1 = zb1 * rel / (1.0d0 + kvobs / c);
    
else
    
    % compute solar potential at object, if within solar system
    
    r = dist(3) * au;
    
    phisun = 0.0d0;
    
    if (r > 1.0d8 && r < 1.0d16)
        
        phisun = gs / r;
        
    end
    
    % compute observed radial velocity measure of a planet or other
    % object -- including a nearby star -- where kinematic
    % barycentric velocity vector is known and gravitational
    % red shift is negligible (lindegren & dravins eq. (40),
    % applied as per s. klioner private communication (2006))
    
    kv = (uk(1) * vel(1) + uk(2) * vel(2) + uk(3) * vel(3)) * toms;
    
    zb1 = (1.0d0 + kv / c) / (10.0d0 - phisun / c2 - 0.5d0 * v2 / c2);
    
    kvobs = (uk(1) * velobs(1) + uk(2) * velobs(2) + uk(3) * velobs(3)) * toms;
    
    zobs1 = zb1 * rel / (1.0d0 + kvobs / c);
    
end

% convert observed radial velocity measure to kilometers/second

rv = (zobs1 - 1.0d0) * c / 1000.0d0;

