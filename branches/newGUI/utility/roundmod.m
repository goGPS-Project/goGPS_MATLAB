function [rx] = roundmod(x,y)

% SYNTAX:
%   [rx] = roundmod(x,y);
%
% INPUT:
%   x = values to be rounded
%   y = resolution
%
% OUTPUT:
%   rx = rounded values
%
% DESCRIPTION:
%   Rounds the input values to the nearest float determined by the
%   resolution in input.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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

remainder = mod(x,y);
rx = zeros(size(remainder));
pos1 = find(remainder<=y/2);
pos2 = find(remainder>y/2);
rx(pos1) = x(pos1) - remainder(pos1);
rx(pos2) = x(pos2) - remainder(pos2) + y;
