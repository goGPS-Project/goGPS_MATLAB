function f = figureColor(n_col, varargin)
% Open figure, amnd set color order with the requested amount of colors
% To have multiple lines (in plot) of different colors!
% Using function Core_UI.getColor
%
% SYNTAX:
%   f = figureColor(n_col, <varargin>)
%
% EXAMPLE:
%   figureColor(200)
%
% INPUT:
%   n_col = maximum number of colors for the ColorOrder
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti ...
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

if ~isempty(varargin)
    f = figure(varargin);
else
    f = figure;
end
ax = gca;
set(ax,'ColorOrder', Core_UI.getColor(1:n_col,n_col), 'NextPlot', 'replacechildren');

