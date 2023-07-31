function [corrTime] = check_t(time)

% SYNTAX:
%   [corrTime] = check_t(time);
%
% INPUT:
%   time = GPS time
%
% OUTPUT:
%   corrTime = corrected GPS time
%
% DESCRIPTION:
%   Function accounting for beginning or end of week crossover. From the
%   Interface Specification document revision E (IS-GPS-200E), page 93.

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


%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) Kai Borre
% Kai Borre 04-01-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

half_week = 302400;     % seconds

corrTime = time;
corrTime(time > half_week) =time(time > half_week) - 2*half_week;
corrTime(time < - half_week) =time(time < - half_week) + 2*half_week;
% if time > half_week
%     corrTime = time - 2*half_week;
% elseif time < -half_week
%     corrTime = time + 2*half_week;
% end
