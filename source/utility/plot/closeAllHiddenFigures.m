function closeAllHiddenFigures(flag_dockable)
% Close all the hidden figures 
% goGPS main windows are not dockable and will remain open ^_^
%
% SINTAX:
%   closeAllHiddenFigures(<flag_dockable = false>);
%
% EXAMPLE:
%   dockAllFigures();
%
% INPUT:
%   flag_dockable = close also hidden dockable figures
%
% DEFAULT VALUES:
%   flag_dockable = false
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
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
    if (nargin == 0) || isempty(flag_dockable)
        flag_dockable = false;
    end

    all_fh = findall(0,'Type','figure');
    for fid = 1:length(all_fh)
        if (~flag_dockable || ((islogical(all_fh(fid).DockControls) && (all_fh(fid).DockControls == true)) || ...
                (ischar(all_fh(fid).DockControls) && (all_fh(fid).DockControls(2) == 'n')) || ...
                (isnumeric(all_fh(fid).DockControls) && (numel(all_fh(fid).DockControls) >= 2)))) && ...
                all_fh(fid).Visible == false
            fprintf('Closing "%s"\n', all_fh(fid).Name);
            delete(all_fh(fid));
        end
    end
end
