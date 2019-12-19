function [pos_S, vel_S] = interpolate_SP3_coord(time, SP3, sat, p_rate)
% SYNTAX:
%   [pos_S, vel_S] = interpolate_SP3_coord(time, SP3, sat);
%
% INPUT:
%   time       = interpolation time (GPS time, continuous since 6-1-1980)
%   SP3        = structure containing precise ephemeris data
%   sat        = satellite PRN
%   p_rate     = processing interval [s]
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
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
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

SP3_time  = SP3.time;
SP3_coord = SP3.coord(:, sat, :);
antPCO    = SP3.antPCO(:, :, sat)';

%number of seconds in a quarter of an hour
interval = SP3.coord_rate;

%find the SP3 epoch closest to the interpolation time
%[~, p] = min(abs(SP3_time - time));
% speed improvement of the above line
% supposing SP3_time regularly sampled
p = round((time - SP3_time(1)) / interval) + 1;

b = SP3_time(p) - time;

pos_S = zeros(3,1);

%Lagrange interpolation
%degree of interpolation polynomial (Lagrange)
% n = 10;
%u = (n/2+1) + (- b + (-1:1))/interval;
u = 6 + (- b + (-1:1))/interval;    % using 6 since n = 10;
%X_sat = fastLI(SP3_coord(:,p+(-n/2:n/2)), u);
X_sat = fastLI(SP3_coord(:,p + (-5:5)), u);

%apply satellite antenna phase center correction
[i, j, k] = satellite_fixed_frame(time, X_sat(:,2), SP3, p_rate);
X_sat(:,2) = X_sat(:,2) + [i j k]*antPCO;

pos_S(1,1) = X_sat(1,2);
pos_S(2,1) = X_sat(2,2);
pos_S(3,1) = X_sat(3,2);

%compute velocity
vel_S = (X_sat(:,3) - X_sat(:,1)) / 2;
