function N = getOrthometricCorr(phi, lam, geoid, method)

% SYNTAX:
%   N = getOrthometricCorr(phi, lam, geoid);
%
% EXAMPLE:
%   core = Core.getInstance;
%   core.initGeoid();
%   getOrthometricCorr(45.69 ./ 180*pi, 9.03 ./ 180*pi, gs.getRefGeoid())
%   % answer should be 46.1767008
%
% INPUT:
%   phi     = geocentric latitude                [deg (rad for legacy method)]
%   lam     = geocentric longitude               [deg (rad for legacy method)]
%   geoid   = regular map in geocentric coordinates <default = EGM08 0.5x0.5 deg>
%   method  = interpolation approach:
%              - legacy
%              - grid
%              - grid_cubic  ( phi, lam are array of grid coordinates )
%              - grid_akima  ( phi, lam are array of grid coordinates )
%              - linear
%              - natural
%
% OUTPUT:
%   N       = geoid ondulation [m]
%
% DESCRIPTION:
%   Get the geoid ondulation (orthometric correction)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
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

if (nargin < 3) || isempty(geoid)
    geoid = Core.getRefGeoid();
end
if (nargin < 4) || isempty(method)
    method = 'legacy';
end

if  (geoid.ncols == 0 || geoid.nrows == 0) 
    Core.initGeoid();
    geoid = Core.getRefGeoid();
end
N = zeros(numel(lam), 1);

switch method
    case 'legacy'
        for i = 1 : numel(lam)
            N(i) = grid_bilin_interp(lam(i) / pi * 180, phi(i) / pi * 180, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);            
        end
    case 'grid'
        x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
        y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
        
        [xmg, ymg] = meshgrid(x_grid, y_grid);
        N = interp2(xmg, ymg, geoid.grid, lam, phi, 'linear');
    case 'grid_cubic'
        x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
        y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
        
        [xmg, ymg] = meshgrid(x_grid, y_grid);
        N = interp2(xmg, ymg, geoid.grid, lam, phi, 'cubic');
    case 'grid_akima'
        x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
        y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
        
        [xmg, ymg] = meshgrid(x_grid, y_grid);
        N = interp2(xmg, ymg, geoid.grid, lam, phi, 'makima');
    case 'linear'
        x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
        y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
        [xmg, ymg] = meshgrid(x_grid, y_grid);
        finterp = scatteredInterpolant(xmg(:), ymg(:), geoid.grid(:), 'linear');
        
        N = finterp(lam, phi);        
    case 'natural'
        x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
        y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
        
        [xmg, ymg] = meshgrid(x_grid, y_grid);
        finterp = scatteredInterpolant(xmg(:), ymg(:), geoid.grid(:), 'natural');        
        N = finterp(lam, phi);
end

