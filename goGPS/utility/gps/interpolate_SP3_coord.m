function [pos_S, vel_S] = interpolate_SP3_coord(time, SP3, sat)

% SYNTAX:
%   [pos_S, vel_S] = interpolate_SP3_coord(time, SP3, sat);
%
% INPUT:
%   time = interpolation time (GPS time, continuous since 6-1-1980)
%   SP3  = structure containing precise ephemeris data
%   sat = satellite PRN
%
% OUTPUT:
%   pos_S = interpolated satellite coordinates
%   vel_S = satellite velocity
%
% DESCRIPTION:
%   SP3 (precise ephemeris) coordinates 1-second interpolation by Lagrange
%   polynomials. Satellite velocity computation. Relativistic correction.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 2
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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

SP3_time  = SP3.time;
SP3_coord = SP3.coord(:, sat, :);
antPCO    = SP3.antPCO(:, :, sat)';

%degree of interpolation polynomial (Lagrange)
n = 10;

%number of seconds in a quarter of an hour
interval = SP3.coord_rate;

%find the SP3 epoch closest to the interpolation time
[~, p] = min(abs(SP3_time - time));

b = SP3_time(p) - time;

pos_S = zeros(3,1);

%extract the SP3 coordinates
SP3_X = SP3_coord(1,p+(-n/2:n/2));
SP3_Y = SP3_coord(2,p+(-n/2:n/2));
SP3_Z = SP3_coord(3,p+(-n/2:n/2));

%Lagrange interpolation (coordinates)
s = 1/interval;
d = n/2+1 - b/interval;
u = d - s : s : d + s;

LI_SP3_X = fastLI(SP3_X, u);
LI_SP3_Y = fastLI(SP3_Y, u);
LI_SP3_Z = fastLI(SP3_Z, u);

X_sat = [LI_SP3_X(2); LI_SP3_Y(2); LI_SP3_Z(2)];

%apply satellite antenna phase center correction
[i, j, k] = satellite_fixed_frame(time, X_sat, SP3);
X_sat = X_sat + [i j k]*antPCO;

pos_S(1,1) = X_sat(1);
pos_S(2,1) = X_sat(2);
pos_S(3,1) = X_sat(3);

pos_S_v = zeros(3,2);
pos_S_v(1,1:2) = [LI_SP3_X(1) LI_SP3_X(3)];
pos_S_v(2,1:2) = [LI_SP3_Y(1) LI_SP3_Y(3)];
pos_S_v(3,1:2) = [LI_SP3_Z(1) LI_SP3_Z(3)];

%compute velocity
vel_S = (pos_S_v(:,2) - pos_S_v(:,1)) / 2;
