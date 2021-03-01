function [ZHD] = saast_dry(P, h, lat)

% SYNTAX:
%   [ZHD] = saast_dry(P, h, lat);
%
% INPUT:
%   P = atmospheric pressure [hPa]
%   h = orthometric height [m]
%   lat = latitude [deg]
%
% OUTPUT:
%   ZHD = Zenith Hydrostatic Delay
%
% DESCRIPTION:
%   Zenith Hydrostatic Delay (ZHD) computation by Saastamoinen model.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
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

%ZHD (Saastamoinen model)
ZHD = 0.0022768 * P(:) .* (1 + 0.00266 * cosd(2*lat(:)) + 0.00000028 * h(:));
