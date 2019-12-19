function [i, j, k] = satellite_fixed_frame(time, X_sat, SP3, p_rate)

% SYNTAX:
%   [i, j, k] = satellite_fixed_frame(time, X_sat, SP3, p_rate);
%
% INPUT:
%   time     = GPS time
%   X_sat    = satellite position (X,Y,Z)
%   SP3      = structure containing precise ephemeris data
%   p_rate   = processing interval [s]
%
% OUTPUT:
%   i = unit vector that completes the right-handed system
%   j = resulting unit vector of the cross product of k vector with the unit vector from the satellite to Sun
%   k = unit vector pointing from the Satellite Mass Centre (MC) to the Earth's centre
%
% DESCRIPTION:
%   Computation of the unit vectors defining the satellite-fixed frame.

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
%  Contributors:     Andrea Gatti...
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

t_sun = SP3.t_sun;
X_sun = SP3.X_sun;

%[~, q] = min(abs(t_sun - time));
% speed improvement of the above line
% supposing t_sun regularly sampled
q = round((time - t_sun(1)) / p_rate) + 1;

X_sun = X_sun(:,q);
e = (X_sun - X_sat) / norm(X_sun - X_sat);
k = -X_sat/norm(X_sat);
%j = cross(k,e);
j = [k(2).*e(3)-k(3).*e(2);
     k(3).*e(1)-k(1).*e(3);
     k(1).*e(2)-k(2).*e(1)];
%i = cross(j,k);
i = [j(2).*k(3)-j(3).*k(2);
     j(3).*k(1)-j(1).*k(3);
     j(1).*k(2)-j(2).*k(1)];

j = j / norm(j);
i = i / norm(i);


