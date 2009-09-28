function [phiD, lambdaD, h] = geoc2geod(phiC, lambdaC, r, a, e)

% SYNTAX:
%   [phiD, lambdaD, h] = geoc2geod(phiC, lambdaC, r, a, e);
%
% INPUT:
%   phiC    = geocentric latitude                [rad]
%   lambdaC = longitude                          [rad]
%   r       = radius                             [m]
%   a       = ellipsoid semi-major axis          [m]
%   e       = ellipsoid eccentricity
%
% OUTPUT:
%   phiD    = geodetic latitude                  [rad]
%   lambdaD = longitude                          [rad]
%   h       = ellipsoidal height                 [m]
%
% DESCRIPTION:
%   Conversion from geocentric to geodetic coordinates.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
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

psi = atan(tan(phiC)/sqrt(1-e^2));

phiD = atan((r.*sin(phiC) + e^2*a/sqrt(1-e^2) * (sin(psi)).^3) ./ ...
    			(r.*cos(phiC) - e^2*a * (cos(psi)).^3));

lambdaD = lambdaC;

N = a ./ sqrt(1 - e^2 * sin(phiD).^2);

h = r .* cos(phiC)./cos(phiD) - N;

%----------------------------------------------------------------------------------------------