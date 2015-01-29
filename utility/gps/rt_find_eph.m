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

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

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
