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
%  Contributors:     Andrea Gatti, ...
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

rx = round(x./y).*y;

% Does this piece of code make sense?
% This should do the same 10 times faster: rx = round(x./y).*y;
%
% remainder = mod(x,y);
% rx = zeros(size(remainder));
% pos1 = find(remainder<=y/2);
% pos2 = find(remainder>y/2);
% rx(pos1) = x(pos1) - remainder(pos1);
% rx(pos2) = x(pos2) - remainder(pos2) + y;
