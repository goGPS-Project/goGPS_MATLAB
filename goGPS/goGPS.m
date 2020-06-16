function goGPS(ini_settings, use_gui, flag_online)
% SYNTAX:
%   goGPS(<ini_settings_file>, <use_gui = 1 >);
%
% INPUT:
%   ini_settings_file       path to the settings file
%   use_gui                 flag to activate GUI editing of the settings
%                           0   MATLAB Command Window only
%                           1   GUI only
%                           -1  GUI but do not open settings GUI 
%                           2   GUI + MATLAB Command Window only
%
% DESCRIPTION:
%   function launcher for goGPS
%   
% OUTPUT:
%   goGPS creates a singleton object, for debug purposes it is possible to obtain it typing:
%   core = Core.getInstance();
%
% EXAMPLE:
%   goGPS('../data/project/default_PPP/config/settings.ini');
%
% COMPILATION STRING:
%   tic; mcc -v -d ./bin/ -m goGPS -a tai-utc.dat -a cls.csv -a icpl.csv -a nals.csv -a napl.csv -a remote_resource.ini -a icons/*.png -a utility/thirdParty/guiLayoutToolbox/layout/+uix/Resources/*.png; toc;
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
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

    %% Preparing execution and settings
    
    if (~isdeployed)
        % add all the subdirectories to the search path
        addPathGoGPS;
    end    
    
    if nargin < 2 || isempty(use_gui)
        if isdeployed
            use_gui = 1;
        else
            use_gui = 1;
        end
    end
    
    if ischar(use_gui)
        use_gui = str2num(use_gui); %#ok<ST2NM>
    end    
        
    log = Logger.getInstance();
    log.setColorMode(Go_Settings.getInstance.isLogColorMode);
    log.disableFileOut();
    if use_gui && Go_Settings.getInstance.getLogMode() == 1
        log.enableGUIOut();
        log.disableScreenOut();
        Core.getMsgGUI(true);
    else
        log.enableScreenOut();
        log.disableGUIOut();
    end
    
    % Show coloured header
    cm = log.getColorMode();
    log.setColorMode(true);    
    Core_UI.showTextHeader();
    log.setColorMode(cm);
    
    if nargin >= 1 && ~isempty(ini_settings)
        core = Core.getInstance(true, false, ini_settings);    
    else
        state = Core.getCurrentSettings();
        core = Core.getInstance(true, false, state); % Init Core
    end        
    core.setModeGUI(use_gui);
    if nargin < 3 || isempty(flag_online)
        flag_online = true;
    end
    
    if ischar(flag_online)
        flag_online = iif(lower(flag_online(1)) == 'f' || lower(flag_online(1)) == '0', false, true);
    end
    
    % Every parameters when the application is deployed are strings
    if isdeployed
        if (ischar(use_gui) && (use_gui == '1'))
            use_gui = true;
        end
    end
    
    if use_gui > 0
        ui = Core_UI.getInstance();
        flag_wait = true;
        ui.openGUI(flag_wait);
        
        if ~ui.isGo()
            return
        end
        ok_go = true;
    end
    
    % Enable file logging
    if core.state.isLogOnFile()
        log.newLine();
        log.enableFileOut();
        log.setOutFile(fullfile(core.state.getOutDir, 'log', 'goGPS_run_${NOW}.log')); % <= to enable project logging
        % log.setOutFile(); <= to enable system logging (save into system folder)
        log.disableFileOut();
        fnp = File_Name_Processor();
        log.simpleSeparator()
        log.addMarkedMessage(sprintf('Logging to: "%s"\n', fnp.getRelDirPath(log.getFilePath, core.state.getHomeDir)));
        log.simpleSeparator()
        log.enableFileOut();
        log.addMessageToFile(Core_UI.getTextHeader());
        core.logCurrentSettings();
    end
    
    if use_gui <= 0 % if this check is not performed by GUI
        % GO goGPS - here the computations start
        err_code = core.checkValidity();        
        ok_go = err_code.go == 0; % here a check on the validity of the parameters should be done
    end            
    
    if ~ok_go
        log.addError('Invalid configuration found! Check the log messages above.');
    else
        
        if nargin >= 3 && ~isempty(flag_online)
            core.prepareProcessing(flag_online); % download important files
        else
            core.prepareProcessing();
        end
        
        ok_go = true; % here a check on the validity of the resources should be done
        
        if ok_go
            core.go(); % execute all
        end
        
        %% Closing all
        if ~use_gui && isdeployed
            close all;
        end
    end
    
    % Stop logging
    if core.state.isLogOnFile()
        log.disableFileOut();
        log.closeFile();
    end
    
    if ~isdeployed && ok_go
        % Do not export to workspace
        %log.addMessage('Execute the script "getResults", to load the object created during the processing');
        
        % Export into workspace
        rec = core.rec;
        assignin('base', 'core', core);
        assignin('base', 'rec', rec);
        
        log.addMarkedMessage('Now you should be able to see 2 variables in workspace:');
        log.addMessage(log.indent(' - core      the core processor object containing all the goGPS structures'));
        log.addMessage(log.indent(' - rec       the array of Receivers'));
        
        screen_log = log.isScreenOut;
        if ~screen_log; log.enableScreenOut(); end
        log.simpleSeparator();
        log.addStatusOk(['Execution completed at ' GPS_Time.now.toString('yyyy/mm/dd HH:MM:SS') '      ^_^']);       
        log.simpleSeparator();
        if ~screen_log; log.disableScreenOut(); end
        
    end
    
    if use_gui
        goInspector;
    end
    % Core_Utils.playAlert();
end

