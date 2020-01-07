function maximizeFig(fig_handle)
% Maximize the figure
% SYNTAX:
%   maximizeFig(<fig_handle>);
%
% EXAMPLE:
%   maximizeFig(gcf);
%
% INPUT:
%   fig_handle = handler to the figure to modify        <optional argument>
%
% DEFAULT VALUES:
%   fig_handle = gcf;
%

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
if nargin < 1
    fig_handle = gcf;
end
    set(fig_handle, 'units','normalized', 'outerposition',[0 0 1 1]);
%     try
%         if isa(fig_handle, 'matlab.ui.Figure')
%             fig_handle = get(fig_handle,'Number');
%         end
%         drawnow % Required to avoid Java errors
%         j_frame = get(fig_handle, 'JavaFrame');
%         j_frame.setMaximized(true);
%     catch
%         % no more supported :-(
%     end
end
