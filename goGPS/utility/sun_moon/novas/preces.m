function pos2 = preces (tjd1, pos1, tjd2)

% this function precesses equatorial rectangular coordinates from
% one epoch to another.  the coordinates are referred to the mean
% dynamical equator and equinox of the two respective epochs.  see
% explanatory supplement to the astronomical almanac, pp. 103-104,
% and capitaine et al. (2003), astronomy and astrophysics 412,
% 567-586.

% input

%  tjd1 = tdb julian date of first epoch

%  pos1 = position vector, geocentric equatorial rectangular
%         coordinates, referred to mean dynamical equator and
%         equinox of first epoch

%  tjd2 = tdb julian date of second epoch

% output

%  pos2 = position vector, geocentric equatorial rectangular
%         coordinates, referred to mean dynamical equator and
%         equinox of second epoch

% note: either tjd1 or tjd2 must be 2451545.0 (j2000.0) tdb

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.d0 * 3600.d0 / pi;

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

if (tjd1 ~= t0 && tjd2 ~= t0 )

    fprintf('\npreces error: precession from jd %14.8f to %14.8f not to/from j2000\n', tjd1, tjd2);

    return
    
end

% t is time in tdb centuries between the two epochs

t = (tjd2 - tjd1) / 36525.0d0;

if (tjd2 == t0)
    
    t = -t;
    
end

% numerical coefficients of psi_a, omega_a, and chi_a, along with
% epsilon_0, the obliquity at j2000.0, are 4-angle formulation
% from capitaine et al. (2003), eqs. (4), (37), & (39)

eps0 = 84381.406d0;

psia   = ( ( ( ( -    0.0000000951d0   * t ...
    +    0.000132851d0  ) * t ...
    -    0.00114045d0   ) * t ...
    -    1.0790069d0    ) * t ...
    + 5038.481507d0     ) * t;

omegaa = ( ( ( ( +    0.0000003337d0   * t ...
    -    0.000000467d0  ) * t ...
    -    0.00772503d0   ) * t ...
    +    0.0512623d0    ) * t ...
    -    0.025754d0     ) * t + eps0;

chia   = ( ( ( ( -    0.0000000560d0   * t ...
    +    0.000170663d0  ) * t ...
    -    0.00121197d0   ) * t ...
    -    2.3814292d0    ) * t ...
    +   10.556403d0     ) * t;

eps0 = eps0 / seccon;

psia = psia / seccon;

omegaa = omegaa / seccon;

chia = chia / seccon;

sa = sin(eps0);

ca = cos(eps0);

sb = sin(-psia);

cb = cos(-psia);

sc = sin(-omegaa);

cc = cos(-omegaa);

sd = sin(chia);

cd = cos(chia);

% compute elements of precession rotation matrix
% equivalent to r3(chi_a)r1(-omega_a)r3(-psi_a)r1(epsilon_0)

xx =  cd * cb - sb * sd * cc;

yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc;

zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc;

xy = -sd * cb - sb * cd * cc;

yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc;

zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc;

xz =  sb * sc;

yz = -sc * cb * ca - sa * cc;

zz = -sc * cb * sa + cc * ca;

if (tjd2 == t0)

    % perform rotation from epoch to j2000.0

    pos2(1) = xx * pos1(1) + xy * pos1(2) + xz * pos1(3);

    pos2(2) = yx * pos1(1) + yy * pos1(2) + yz * pos1(3);

    pos2(3) = zx * pos1(1) + zy * pos1(2) + zz * pos1(3);

else

    % perform rotation from j2000.0 to epoch

    pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3);

    pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3);

    pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3);

end



