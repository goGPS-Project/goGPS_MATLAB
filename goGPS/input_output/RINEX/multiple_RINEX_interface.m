function [filename_obs] = multiple_RINEX_interface(filename_R_obs, filename_M_obs, mode)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
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

%if multiple RINEX files in input
if (iscell(filename_R_obs))
    %if mode is not multi-receiver, fall back to the first rover filename
    if (~goGNSS.isMR(mode) && numel(filename_R_obs)>1)
        filename_obs = filename_R_obs(1);
        fprintf('... WARNING: multiple rover RINEX files in input for a single rover mode; using the first file.\n');
    else
        nFiles = length(filename_R_obs);
        filename_obs = cell(nFiles,1);
        for f = 1 : nFiles
            filename_obs{f,1} = filename_R_obs{f};
        end
    end
else
    filename_obs{1} = filename_R_obs;
end

%if relative positioning
if (goGNSS.isDD(mode))
    filename_obs = [filename_obs; filename_M_obs];
end
