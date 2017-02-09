function [i, j, k] = satellite_fixed_frame(time, X_sat, SP3)

% SYNTAX:
%   [i, j, k] = satellite_fixed_frame(time, X_sat, SP3);
%
% INPUT:
%   time  = GPS time
%   X_sat = satellite position (X,Y,Z)
%   SP3   = structure containing precise ephemeris data
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
%    |___/                    v 0.5.0
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

t_sun = SP3.t_sun;
X_sun = SP3.X_sun;

[~, q] = min(abs(t_sun - time));
X_sun = X_sun(:,q);
e = (X_sun-X_sat)/norm(X_sun-X_sat);
k = -X_sat/norm(X_sat);
j = cross(k,e);
i = cross(j,k);
j = j/norm(j);
i = i/norm(i);
