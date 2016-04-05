function [phi, lam, h] = cart2geod(X, Y, Z)

% SYNTAX:
%   [phi, lam, h] = cart2geod(X, Y, Z);
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
%
% DESCRIPTION:
%   Conversion from cartesian coordinates to geodetic coordinates.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

%global a_GPS e_GPS

a = goGNSS.ELL_A_GPS;
e = goGNSS.ELL_E_GPS;

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
