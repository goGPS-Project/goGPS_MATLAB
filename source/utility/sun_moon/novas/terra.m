function [pos, vel] = terra (glon, glat, ht, st)

% this function computes the position and velocity vectors of
% a terrestrial observer with respect to the geocenter.

% input

%  glon = longitude of observer with respect to reference
%           meridian (east +) in degrees

%  glat = geodetic latitude (north +) of observer in degrees

%  ht   = height of observer in meters

%  st   = local apparent sidereal time at reference meridian in hours

% output

%  pos = position vector of observer with respect to
%        geocenter, equatorial rectangular coordinates,
%        referred to true equator and equinox of date, components in au

%  vel = velocity vector of observer with respect to
%        geocenter, equatorial rectangular coordinates,
%        referred to true equator and equinox of date, components in au/day

% note 1:  if reference meridian is greenwich and st=0.d0, pos
% is effectively referred to equator and greenwich.

% note 2:  this function ignores polar motion, unless the
% observer's longitude and latitude have been corrected for it,
% and variation in the length of day (angular velocity of earth).
% neglect of polar motion may yield 15 meters error in position
% and of order 1 millimeter/sec error in velocity.  neglect of
% variations in length of day results in even smaller velocity
% errors.

% note 3:  the true equator and equinox of date do not form an
% inertial system.  therefore, with respect to an inertial system,
% the small velocity component, of order 0.1 millimeter/sec,
% due to the precession and nutation of the earth's axis, is not
% accounted for here.

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.0d0 * 3600.0d0 / pi;

% equatorial radius of earth in kilometers

erad = astcon('erad', 1.0d-3);

% flattening factor of earth ellipsoid

f = astcon('f', 1.0d0);

% nominal mean rotational angular velocity of earth in radians/second

omega = astcon('angvel', 1.0d0);

% au in kilometers

aukm = astcon('au', 1.0d-3);

% compute parameters relating to geodetic to geocentric conversion

df2 = (1.0d0 - f)^2;

phi = glat * 3600.0d0 / seccon;

sinphi = sin(phi);

cosphi = cos(phi);

c = 1.0d0 / sqrt (cosphi^2 + df2 * sinphi^2);

s = df2 * c;

ach = erad * c + ht / 1000.0d0;

ash = erad * s + ht / 1000.0d0;

% compute local sidereal time factors

stlocl = (st * 54000.0d0 + glon * 3600.0d0) / seccon;

sinst = sin(stlocl);

cosst = cos(stlocl);

% compute position vector components in km

pos(1) = ach * cosphi * cosst;

pos(2) = ach * cosphi * sinst;

pos(3) = ash * sinphi;

% compute velocity vector components in km/sec

vel(1) = -omega * ach * cosphi * sinst;

vel(2) =  omega * ach * cosphi * cosst;

vel(3) =  0.0d0;

% convert position and velocity components to au and au/day

for j = 1:1:3

    pos(j) = pos(j) / aukm;

    vel(j) = vel(j) / aukm * 86400.0d0;

end


