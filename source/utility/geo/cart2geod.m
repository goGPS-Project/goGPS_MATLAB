function [phi, lam, h, phiC] = cart2geod(X, Y, Z)

% SYNTAX:
%   [phi, lam, h, phiC] = cart2geod(X, Y, Z);
%   [phi, lam, h, phiC] = cart2geod(xyz);
%
% INPUT:
%   X = X axis cartesian coordinate
%   Y = Y axis cartesian coordinate
%   Z = Z axis cartesian coordinate
%
% OUTPUT:
%   phi = latitude
%   lam = longitude
%   h = ellipsoidal height
%   phiC = geocentric latitude
%
% DESCRIPTION:
%   Conversion from cartesian coordinates to geodetic coordinates.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

%global a_GPS e_GPS

if nargin == 1
    Z = X(:, 3);
    Y = X(:, 2);
    X = X(:, 1);
end

a = GPS_SS.ELL_A;
e = GPS_SS.ELL_E;

%radius computation
r = sqrt(X.^2 + Y.^2 + Z.^2);

%longitude
lam = atan2(Y,X);

%geocentric latitude
phiC = atan(Z./sqrt(X.^2 + Y.^2));

%coordinate transformation
psi = atan(tan(phiC)/sqrt(1-e^2));

phi = atan((r.*sin(phiC) + e^2*a/sqrt(1-e^2) * (sin(psi)).^3) ./ ...
    			(r.*cos(phiC) - e^2*a * (cos(psi)).^3));

N = a ./ sqrt(1 - e^2 * sin(phi).^2);

%height
h = r .* cos(phiC)./cos(phi) - N;
