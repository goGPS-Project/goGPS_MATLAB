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
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
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

    if not(isa(fh, 'matlab.ui.Figure'))
        axis_id = fh;
        fh = gcf;
    end
    set(0, 'CurrentFigure', fh);
    try
        ax = fh.Children(max(1, end + 1 - axis_id));
        subplot(ax);
    catch
        ax = axes();
    end
end