function [EAST, NORTH, h] = cart2plan(X, Y, Z)

% SYNTAX:
%   [EAST, NORTH, h] = cart2plan(X, Y, Z);
%
% INPUT:
%   X = X axis cartesian coordinate
%   Y = Y axis cartesian coordinate
%   Z = Z axis cartesian coordinate
%
% OUTPUT:
%   EAST = EAST coordinate
%   NORTH = NORTH coordinate
%   h = ellipsoidal height
%
% DESCRIPTION:
%   Conversion from cartesian coordinates to planimetric coordinates (UTM WGS84).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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

global a e

%radius computation
r = sqrt(X.^2 + Y.^2 + Z.^2);

%longite
lam = atan2(Y,X);

%latitude
phi = atan(Z./sqrt(X.^2 + Y.^2));

%coordinate transformation
[phi, lam, h] = geoc2geod(phi, lam, r, a, e);

%projection to UTM
[EAST, NORTH] = geod2plan(phi, lam);

%-------------------------------------------------------------------------------