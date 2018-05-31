function goGPS(ini_settings_file, use_gui)
% SYNTAX:
%   goGPS(<ini_settings_file>, <use_gui =false>);
%
% INPUT:
%   ini_settings_file       path to the settings file
%   use_gui                 (0/1) flag to activate GUI editing of the settings
%                           default = 0 (false
%
% DESCRIPTION:
%   function launcher for goGPS
%   
% OUTPUT:
%   goGPS creates a singleton object, for debug purposes it is possibble to obtain it typing:
%   core = Core.getInstance();
%
% EXAMPLE:
%   goGPS('../data/project/default_PPP/config/settings.ini');
%
% COMPILATION STRING:
%   tic; mcc -v -d ./bin/ -m goGPS -a tai-utc.dat -a cls.csv -a icpl.csv -a nals.csv -a napl.csv; toc;
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
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

% if the plotting gets slower than usual, there might be problems with the
% Java garbage collector. In case, you can try to use the following
% command:
%
% java.lang.System.gc() %clear the Java garbage collector
%
% or:
%
% clear java

% add all the subdirectories to the search path

    %% Preparing execution and settings
    if (~isdeployed)
        addpath(genpath(pwd));
    end
    
    core = Core.getInstance(true); % Init Core

    if nargin >= 1
        if  ~exist(ini_settings_file,'file')
            clear ini_settings_file;
        else            
            core.importIniFile(ini_settings_file);
        end
    end
    if nargin < 2
        if isdeployed
            use_gui = true;
       else
            use_gui = true;
        end
    end

    % Every parameters when the application is deployed are strings
    if isdeployed
        if (ischar(use_gui) && (use_gui == '1'))
            use_gui = true;
        end
    end
    
    if use_gui
        ui = Core_UI.getInstance();
        ui.openGUI();
        
        if ~ui.isGo()
            return
        end
    end
    
    %% GO goGPS - here the computations start
    ok_go = true; % here a check on the validity of the parameters should be done

    core.prepareProcessing(true); % download important files
    
    ok_go = true; % here a check on the validity of the resources should be done

    if ok_go
        core.go(); % execute all
    end

    %% Closing all
    if ~use_gui
        close all;
    end
    
    if ~isdeployed && ok_go
        log = Logger.getInstance();
        log.addMessage('Execute the script "getResults", to load the object created during the processing');
    end
end

