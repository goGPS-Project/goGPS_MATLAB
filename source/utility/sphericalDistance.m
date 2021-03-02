function dist = sphericalDistance(latd_a, lond_a, latd_b, lond_b)
% SYNTAX:
%   dist = sphericalDistance(latd_a, lond_a, latd_b, lond_b);
%
% INPUT:
%   latd_a  = latitude of point a  [degree]
%   lond_a  = longitude of point a [degree]
%   latd_b  = latitude of point b  [degree]
%   lond_b  = longitude of point b [degree]
%
% OUTPUT:
%   dist = spherical distance in degrees
%
% DESCRIPTION:
%   Compute spherical distance among points (a or b can be an array)
%

    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 1.0RC1
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
    %  Written by: Andrea Gatti
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
    
    dist = 2 * asind(sqrt(sind((latd_a - latd_b) / 2) .^2 + cosd(latd_a) .* cosd(latd_b) .* sind((lond_a - lond_b) / 2).^2));
end
