function [X,Y,Z] = geod2cart (phi, lam, h, a, f)

% SYNTAX:
%   [X,Y,Z] = geod2cart (phi, lam, h, a, f);
%
% INPUT:
%   phi     = geodetic latitude                  [rad]
%   lam     = geodetic longitude                 [rad]
%   h       = ellipsoid height                   [m]
%   a       = ellipsoid semi-major axis          [m]
%   f       = ellipsoid flattening
%
% OUTPUT:
%   X       = X cartesian coordinate
%   Y       = Y cartesian coordinate
%   Z       = Z cartesian coordinate
%
% DESCRIPTION:
%   Conversion from geodetic to geocentric cartesian coordinates.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
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

if (nargin == 3)
    cc = Constellation_Collector('G');
    a = cc.gps.ELL_A;
    e = cc.gps.ELL_E;
    e2 = e^2;
else
    e = sqrt(1-(1-f)^2);
    e2 = 1 - (1 - f)^2;
end
N = a ./ sqrt(1 - e2 * sin(phi).^2);

X = (N + h) .* cos(lam) .* cos(phi);
Y = (N + h) .* sin(lam) .* cos(phi);
Z = (N * (1 - e2) + h) .* sin(phi);
