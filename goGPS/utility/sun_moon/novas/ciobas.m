function [x, y, z] = ciobas (tjd, racio, k)

% this function returns the orthonormal basis vectors, with
% respect to the gcrs (geocentric icrs), of the celestial
% intermediate system defined by the celestial intermediate pole
% (cip) (in the z direction) and the celestial intermediate origin
% (cio) (in the x direction).  a tdb julian date and the right
% ascension of the cio at that date is required as input.  the
% right ascension of the cio can be with respect to either the
% gcrs origin or the true equinox of date -- different algorithms
% are used in the two cases.

% input

%  tjd   = tdb julian date

%  racio = right ascension of the cio, in hours

%  k     = reference system in which right ascension is expressed

%          set k = 1 for gcrs

%          set k = 2 for true equator and equinox of date

% output

%   x = unit vector toward the cio, equatorial rectangular
%       coordinates, referred to the gcrs

%   y = unit vector toward the y-direction, equatorial
%       rectangular coordinates, referred to the gcrs

%   z = unit vector toward north celestial pole (cip),
%       equatorial rectangular coordinates, referred to the gcrs

% important: set the values of tlast and klast in the main script

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

global xx yy zz

radcon = pi / 180.0d0;

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

z0(1) = 0.0d0;

z0(2) = 0.0d0;

z0(3) = 1.0d0;

% compute unit vector z toward celestial pole (cip)

w1 = nutate (-tjd, z0);

w2 = preces (tjd, w1, t0);

zz = frame (w2, -1);

% ra of cio expressed in gcrs

if (k == 1)

    % compute vector x toward cio in gcrs

    sinra = sin(racio * 15.0d0 * radcon);

    cosra = cos(racio * 15.0d0 * radcon);

    xx(1) =  zz(3) * cosra;

    xx(2) =  zz(3) * sinra;

    xx(3) = -zz(1) * cosra - zz(2) * sinra;

    % normalize vector x

    xmag = sqrt (xx(1)^2 + xx(2)^2 + xx(3)^2);

    xx(1) = xx(1) / xmag;

    xx(2) = xx(2) / xmag;

    xx(3) = xx(3) / xmag;

    % compute unit vector y orthogonal to x and z (y = z cross x)

    yy(1) = zz(2) * xx(3) - zz(3) * xx(2);

    yy(2) = zz(3) * xx(1) - zz(1) * xx(3);

    yy(3) = zz(1) * xx(2) - zz(2) * xx(1);

    % ra of cio expressed in equator-and-equinox of date system

elseif (k == 2)

    % construct unit vector toward cio in equator-and-equinox-of-date system

    w0(1) = cos(racio * 15.0d0 * radcon);

    w0(2) = sin(racio * 15.0d0 * radcon);

    w0(3) = 0.0d0;

    % rotate the vector into the gcrs to form unit vector x

    w1 = nutate (-tjd, w0);

    w2 = preces (tjd, w1, t0);

    xx = frame (w2, -1);

    % compute unit vector y orthogonal to x and z (y = z cross x)

    yy(1) = zz(2) * xx(3) - zz(3) * xx(2);

    yy(2) = zz(3) * xx(1) - zz(1) * xx(3);

    yy(3) = zz(1) * xx(2) - zz(2) * xx(1);
    
else

    fprintf('\n ciobas error: invalid value for k for jd %f14.8', tjd);
    
    pause;

end

for j = 1:1:3

    x(j) = xx(j);

    y(j) = yy(j);

    z(j) = zz(j);
    
end

