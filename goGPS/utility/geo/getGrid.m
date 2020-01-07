function [phiGrid, lambdaGrid] = getGrid(step, phiMin, phiMax, lambdaMin, lambdaMax)

% SYNTAX:
%   [phiGrid, lambdaGrid] = getGrid(step, phiMin, phiMax, lambdaMin, lambdaMax)
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
% DESCRIPTION:
%   Get the array of knots of the grid

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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

    if nargin ==1
        phiMin = -90;
        phiMax = 90;
        lambdaMin = -180;
        lambdaMax = 180;
    end

    if phiMin>phiMax
        tmp = phiMin;
        phiMin = phiMax;
        phiMax = tmp;
    end

    if lambdaMin > lambdaMax
        tmp = lambdaMin;
        lambdaMin = lambdaMax;
        lambdaMax = tmp;
    end

    phiGrid = (phiMax - step/2 : -step : phiMin + step / 2)';
    lambdaGrid = (lambdaMin + step / 2 : step : lambdaMax - step / 2)';
end
