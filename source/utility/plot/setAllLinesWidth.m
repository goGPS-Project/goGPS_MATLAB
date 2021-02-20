function setAllLinesWidth(h, width)
% setAllLinesWidth get all the lines handler and change them to a <width> value
%
% SINTAX:
%   setAllLinesWidth(<h>,width)
%   setAllLinesWidth(width)
%
% EXAMPLE:
%   setAllLinesWidth(gcf,2);
%   setAllLinesWidth(2);
%
% INPUT:
%   h       = handler to the figure to modify           <optional argument>
%   width   = requested width of all the lines
%
% DEFAULT VALUES:
%   h       = gcf
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
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

if nargin < 2
    width = h;
    h = gcf;
end
hline = findobj(h, 'type', 'line');
set(hline, 'LineWidth', width);
