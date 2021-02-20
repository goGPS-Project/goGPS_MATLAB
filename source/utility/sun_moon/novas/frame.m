function pos2 = frame (pos1, k)

% this function transforms a vector from the dynamical reference
% system to the international celestial reference system (icrs),
% or vice versa.  the dynamical reference system is based on the
% dynamical mean equator and equinox of j2000.0.  the icrs is
% based on the space-fixed icrs axes defined by the radio catalog
% positions of several hundred extragalactic objects.  the rotation
% matrix used here is equivalent to that given by hilton and
% hohenkerk (2004), astronomy and astrophysics 413, 765-770,
% eq. (6) and (8).

% input

%  pos1 = position vector, equatorial rectangular coordinates

%  k    = direction of rotation

%         set k < 0 for dynamical to icrs
%         set k > 0 for icrs to dynamical

% output

%  pos2 = position vector, equatorial rectangular coordinates

% note:  for geocentric coordinates, the same transformation is
% used between the dynamical reference system and the gcrs.

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.0d0 * 3600.0d0 / pi;

% xi0, eta0, and da0 are icrs frame biases in arcseconds taken
% from iers conventions (2003), chapter 5

xi0 = -0.0166170d0;

eta0 = -0.0068192d0;

da0 = -0.01460d0;

% compute elements of rotation matrix (to first order)

xx =  1.0d0;
yx = -da0  / seccon;
zx =  xi0  / seccon;

xy =  da0  / seccon;
yy =  1.0d0;
zy =  eta0 / seccon;

xz = -xi0  / seccon;
yz = -eta0 / seccon;
zz =  1.0d0;

% include second-order corrections to diagonal elements

xx = 1.0d0 - 0.5d0 * (yx^2 + zx^2);

yy = 1.0d0 - 0.5d0 * (yx^2 + zy^2);

zz = 1.0d0 - 0.5d0 * (zy^2 + zx^2);

if (k >= 0)

    % perform rotation from icrs to dynamical system

    pos2(1) = xx * pos1(1) + xy * pos1(2) + xz * pos1(3);

    pos2(2) = yx * pos1(1) + yy * pos1(2) + yz * pos1(3);

    pos2(3) = zx * pos1(1) + zy * pos1(2) + zz * pos1(3);

else

    % perform rotation from dynamical system to icrs

    pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3);

    pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3);

    pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3);

end



