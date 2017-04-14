function goGPSgo(ini_settings_file, use_gui) %#ok<INUSD>
% SYNTAX:
%   goGPSgo(<ini_settings_file>, <use_gui =false>);
%
% INPUT:
%   ini_settings_file       path to the settings file
%   use_gui                 flag to activate GUI editing of the settings
%                           default = false
%  
% DESCRIPTION:
%   function launcher for goGPS
%
% EXAMPLE:
%   goGPSgo('../data/project/default_PPP/config/settings.ini');

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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
    if nargin >= 1
        ini_settings_file = checkPath(ini_settings_file);
        if  ~exist(ini_settings_file,'file')
            clear ini_settings_file;
        end
    end
    if nargin < 2
        use_gui = false; %#ok<NASGU>
    end
    goGPS

end
