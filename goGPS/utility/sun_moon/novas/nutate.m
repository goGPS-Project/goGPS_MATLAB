function pos2 = nutate (tjd, pos1)

% this function nutates equatorial rectangular coordinates from
% the mean dynamical equator and equinox of epoch to the true
% equator and equinox of epoch. see explanatory supplement to the
% astronomical almanac, pp. 114-115.

% input

%  tjd  = tdb julian date of epoch

%  pos1 = position vector, geocentric equatorial rectangular
%         coordinates, referred to mean dynamical equator and
%         equinox of epoch

% output

%  pos2 = position vector, geocentric equatorial rectangular
%         coordinates, referred to true equator and equinox of epoch

% note:  if tjd is negative, inverse nutation (true to mean) is applied.

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.d0 * 3600.d0 / pi;

tjd1 = abs(tjd);

[oblm, oblt, eqeq, dpsi, deps] = etilt (tjd1);

oblm = oblm * 3600.d0 / seccon;

oblt = oblt * 3600.d0 / seccon;

dpsi = dpsi / seccon;

deps = deps / seccon;

cobm = cos(oblm);

sobm = sin(oblm);

cobt = cos(oblt);

sobt = sin(oblt);

cpsi = cos(dpsi);

spsi = sin(dpsi);

% compute elements of nutation rotation matrix

xx =  cpsi;
yx = -spsi * cobm;
zx = -spsi * sobm;

xy =  spsi * cobt;
yy =  cpsi * cobm * cobt + sobm * sobt;
zy =  cpsi * sobm * cobt - cobm * sobt;

xz =  spsi * sobt;
yz =  cpsi * cobm * sobt - sobm * cobt;
zz =  cpsi * sobm * sobt + cobm * cobt;

if (tjd < 0.0d0)

    % perform rotation from true to mean

    pos2(1) = xx * pos1(1) + xy * pos1(2) + xz * pos1(3);

    pos2(2) = yx * pos1(1) + yy * pos1(2) + yz * pos1(3);

    pos2(3) = zx * pos1(1) + zy * pos1(2) + zz * pos1(3);

else

    % perform rotation from mean to true

    pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3);

    pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3);

    pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3);
    
end

