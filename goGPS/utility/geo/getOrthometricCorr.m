function N = getOrthometricCorr(phi, lam, geoid)

% SYNTAX:
%   N = getOrthometricCorr(phi, lam, geoid);
%
% EXAMPLE:
%   gs = Go_State.getInstance;
%   gs.initGeoid();
%   getOrthometricCorr(45.69 ./ 180*pi, 9.03 ./ 180*pi, gs.getRefGeoid())
%
% INPUT:
%   phi     = geocentric latitude                [rad]
%   lam     = geocentric longitude               [rad]
%   geoid   = regular map in geocentric coordinates <default = EGM08 0.5x0.5 deg>
%
% OUTPUT:
%   N       = geoid ondulation
%
% DESCRIPTION:
%   Get the geoid ondulation (orthometric correction)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
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

if (nargin == 2)
    gs = Go_State.getInstance();
    geoid = gs.getRefGeoid();
end

N = zeros(numel(lam), 1);
for i = 1 : numel(lam)
    N(i) = grid_bilin_interp(lam(i) / pi * 180, phi(i) / pi * 180, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
end


