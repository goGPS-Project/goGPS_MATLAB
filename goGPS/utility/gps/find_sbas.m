function [sbas_t] = find_sbas(sbas, t)

% SYNTAX:
%   [sbas_t] = find_sbas(sbas, t);
%
% INPUT:
%   sbas = full SBAS data structure
%   t    = current epoch index
%
% OUTPUT:
%   sbas_t = selected SBAS data
%
% DESCRIPTION:
%   Extract the SBAS data referred to the current epoch.

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


if (~isempty(sbas))

    sbas_t = struct('prc',     sbas.prc(t,:), ...
                    'dx',      sbas.dx(t,:), ...
                    'dy',      sbas.dy(t,:), ...
                    'dz',      sbas.dz(t,:), ...
                    'doffset', sbas.doffset(t,:), ...
                    'iode',    sbas.iode(t,:), ...
                    'ivd',     sbas.ivd(t,:), ...
                    'igp',     sbas.igp, ...
                    'lat_igp', sbas.lat_igp, ...
                    'lon_igp', sbas.lon_igp);
else
    sbas_t = [];
end
