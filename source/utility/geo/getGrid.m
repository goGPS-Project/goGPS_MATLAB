function [phi_grid, lambda_grid] = getGrid(step, phi_min, phi_max, lambda_min, lambda_max)
%  Get the array of knots of the grid
%
% INPUT:
%   step:   step of the grid
%   ...     limits of the grid
%
% OUTPUT:
%   phiGrid    = array of knots of the grid (phi) as
%                phiMax - step/2 : -step : phiMin + step / 2;
%   lambdaGrid = array of knots of the grid (lambda) as
%                lambdaMin + step / 2 : step : lambdaMax - step / 2;
%
% SYNTAX:
%   [phi_grid, lambda_grid] = getGrid(step, phi_min, phi_max, lambda_min, lambda_max)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2011 Andrea Gatti
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

    if nargin == 1
        phi_min = -90;
        phi_max = 90;
        lambda_min = -180;
        lambda_max = 180;
    end

    if phi_min > phi_max
        tmp = phi_min;
        phi_min = phi_max;
        phi_max = tmp;
    end

    if lambda_min > lambda_max
        tmp = lambda_min;
        lambda_min = lambda_max;
        lambda_max = tmp;
    end
    
    phi_grid = (phi_max - step(1)/2 : -step(1) : phi_min + step(1) / 2)';
    lambda_grid = (lambda_min + step(end) / 2 : step(end) : lambda_max - step(end) / 2)';
end
