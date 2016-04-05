function pos2 = wobble (tjd, x, y, pos1)

% this function corrects a vector in the itrs (a rotating earth-
% fixed system) for polar motion, and also corrects the longitude
% origin (by a tiny amount) to the terrestrial intermediate origin
% (tio).  the itrs vector is thereby transformed to the terrestrial
% intermediate system, based on the true (rotational) equator and
% the terrestrial intermediate origin (tio).  since the true equator
% is the plane orthogonal to the direction of the celestial
% intermediate pole (cip), the components of the output vector are
% referred to z and x axes toward the cip and tio, respectively.

% input

%  tjd  = tt or ut1 julian date 

%  x    = conventionally-defined x coordinate of celestial
%         intermediate pole with respect to itrs pole, in arcseconds 

%  y    = conventionally-defined y coordinate of celestial
%         intermediate pole with respect to itrs pole, in arcseconds 

%  pos1 = position vector, geocentric equatorial rectangular
%         coordinates, referred to itrs axes 

% output

%  pos2 = position vector, geocentric equatorial rectangular
%         coordinates, referred to true equator and tio

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.0d0 * 3600.0d0 / pi;

% t0 = tt julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

xpole = x / seccon;

ypole = y / seccon;

t = (tjd - t0) / 36525.0d0;

% compute approximate longitude of tio, using eq. (10) of
% lambert & bizouard (2002), astronomy and astrophysics 394,
% 317-321

sprime = -47.0d-6 * t;

tiolon = -sprime / seccon;

% note that tiolon, the longitude correction, is negligible for
% most astronomical purposes

% compute elements of rotation matrix
% equivalent to r3(-s')r2(x)r1(y) as per iers conventions (2003)

sinx = sin(xpole);
cosx = cos(xpole);

siny = sin(ypole);
cosy = cos(ypole);

sinl = sin(tiolon);
cosl = cos(tiolon);

xx =  cosx * cosl;
yx =  sinx * siny * cosl + cosy * sinl;
zx = -sinx * cosy * cosl + siny * sinl;

xy = -cosx * sinl;
yy =  sinx * siny * sinl + cosy * cosl;
zy =  sinx * cosy * sinl + siny * cosl;

xz =  sinx;
yz = -cosx * siny;
zz =  cosx * cosy;

% perform rotation

pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3);

pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3);

pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3);


