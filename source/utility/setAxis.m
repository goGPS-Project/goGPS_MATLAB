function ax = setAxis(fh, axis_id)
% Enabled a certain axis in the figure, update MATLAB curfigure status
% useful to reset the current axis in case of new plotting
% 
% INPUT
%   fh       figurehandle
%   axis_id  number of the axis (1 is the first inserted)
%
% SYNTAX
%   setAxis(fh, axis_id)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by: Giulio Tagliaferro
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
    if nargin == 0
        fh = gcf;
        axis_id = 1;
    elseif nargin == 1
        axis_id = 1;
    end
    
    if not(isa(fh, 'matlab.ui.Figure'))
        axis_id = fh;
        fh = gcf;
    end
    set(0, 'CurrentFigure', fh);
    ax_list = fh.Children;
    % filter Axes
    id_ax = false(numel(ax_list),1);
    for a = 1:numel(ax_list)
        id_ax(a) = isa(ax_list(a), 'matlab.graphics.axis.Axes');
    end
    ax_list = ax_list(id_ax);
    try
        ax = ax_list(max(1, end + 1 - axis_id));
        subplot(ax);
    catch
        ax = axes();
    end
end
