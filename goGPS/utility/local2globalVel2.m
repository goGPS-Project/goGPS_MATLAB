function [y] = local2globalVel2(V, lon,lat)

% SYNTAX:
%   [y] = local2globalVel2(V, X);
%
% INPUT:
%   V = local position vector(s)
%   lon = longitude of orifin vector in radians
%   lat = latirude of orifin vector in radians
%
% OUTPUT:
%   y = global position vector(s)
%
% DESCRIPTION:
%   Rototation from local-level reference frame to Earth-fixed reference frame

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
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

%initialize new position vector
y = zeros(size(V));

for i = 1 : size(V,2)
    %rotation matrix from global to local reference system
    R = [-sin(lon) cos(lon) 0;
         -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
         +cos(lat)*cos(lon) +cos(lat)*sin(lon) sin(lat)];

    %rototraslation
    y(:,i) = R\V(:,i);
end
