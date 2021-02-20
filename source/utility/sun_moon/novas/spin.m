function pos2 = spin (angl, pos1)

% this function transforms a vector from one coordinate system
% to another with same origin and axes rotated about the z axis.

% input

%  angl = angle of coordinate system rotation, positive
%         counterclockwise when viewed from +z, in degrees

%  pos1 = position vector

% output

%  pos2 = position vector expressed in new coordinate
%         system rotated about z by angle ang

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

ang = angl / 180.0d0 * pi;

cosang = cos(ang);

sinang = sin(ang);

% rotation matrix follows

xx =  cosang;
yx =  sinang;
zx =  0.0d0;

xy = -sinang;
yy =  cosang;
zy =  0.0d0;

xz =  0.0d0;
yz =  0.0d0;
zz =  1.0d0;

% perform rotation

pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3);

pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3);

pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3);


