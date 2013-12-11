function [Eph_t] = rt_find_eph (Eph_in, time)

% SYNTAX:
%   [Eph_t] = rt_find_eph (Eph_in, time);
%
% INPUT:
%   Eph_in = ephemerides in input
%   time = GPS time
%
% OUTPUT:
%   Eph_t = selected ephemerides
%
% DESCRIPTION:
%   Extract the ephemerides referred to the current epoch
%   and save them in the global variable 'eph' (size 21x32).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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

for sv = 1 : 32
    icol = find_eph(Eph_in, sv, time);
    %dt = Eph_in(21,icol) - time;
    %if (~isempty(icol) & (dt < 0))
    if (~isempty(icol))
        Eph_t(:,sv) = Eph_in(:,icol);
    end
end
