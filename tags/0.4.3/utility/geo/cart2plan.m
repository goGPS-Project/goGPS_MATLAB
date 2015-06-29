function [EAST, NORTH, h, utm_zone] = cart2plan(X, Y, Z)

% SYNTAX:
%   [EAST, NORTH, h, utm_zone] = cart2plan(X, Y, Z);
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
%   utm_zone = UTM zone (example: '32 T')
%
% DESCRIPTION:
%   Conversion from cartesian coordinates to planimetric coordinates (UTM WGS84).

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

%coordinate transformation
[phi, lam, h] = cart2geod(X, Y, Z);

%projection to UTM
[EAST, NORTH, utm_zone] = geod2plan(phi, lam);
