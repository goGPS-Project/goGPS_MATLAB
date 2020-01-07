function [eph_t] = rt_find_eph(eph_in, time, n_sat)

% SYNTAX:
%   [eph_t] = rt_find_eph(eph_in, time, n_sat);
%
% INPUT:
%   eph_in = ephemerides in input
%   time   = GPS time
%   n_sat  = total number of satellites (depending on enabled constellations)
%
% OUTPUT:
%   eph_t = selected ephemerides
%
% DESCRIPTION:
%   Extract the ephemerides referred to the current epoch.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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

eph_t     = zeros(33, n_sat);

for sv = 1 : n_sat
    icol = find_eph(eph_in, sv, time);
    if (~isempty(icol))
        eph_t(:, sv) = eph_in(:, icol);
    end
end
