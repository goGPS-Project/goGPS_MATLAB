function [Eph_t] = rt_find_eph(Eph_in, time, nsat)

% SYNTAX:
%   [Eph_t] = rt_find_eph(Eph_in, time, nsat);
%
% INPUT:
%   Eph_in = ephemerides in input
%   time = GPS time
%   nsat = total number of satellites (depending on enabled constellations)
%
% OUTPUT:
%   Eph_t = selected ephemerides
%
% DESCRIPTION:
%   Extract the ephemerides referred to the current epoch.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___ 
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta
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

empty_col = zeros(33,1);
Eph_t     = zeros(33,nsat);

for sv = 1 : nsat
    icol = find_eph(Eph_in, sv, time);
    if (~isempty(icol))
        Eph_t(:,sv) = Eph_in(:,icol);
    else
        Eph_t(:,sv) = empty_col;
    end
end
