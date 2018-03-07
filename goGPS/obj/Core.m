%   CLASS Core
% =========================================================================
%
% DESCRIPTION
%   Collector of settings to manage a useful parameters of goGPS
%   This singleton class collects multiple objects containing various
%   parameters
%
% EXAMPLE
%   core = Core();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 1 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef Core < handle

    %% PROPERTIES CONSTANTS
    % ==================================================================================================================================================
    properties (Constant)
        GO_GPS_VERSION = '0.6.0 alpha 1 - nightly';
        GUI_MODE = 0; % 0 means text, 1 means GUI, 5 both
    end

    %% PROPERTIES SINGLETON POINTERS
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons
        log
        gc
        state
        w_bar
        sky
        cmd
    end

    %% PROPERTIES RECEIVERS    
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons
        rec      % List of all the receiver used in a session
        
        rec_list % List of all the receiver used in all the session
    end

    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = private)
        % Concrete implementation.  See Singleton superclass.
        function this = Core()
            % Core object creator
            this.log = Logger.getInstance();
            this.init();
            
            this.gc = Global_Configuration.getInstance();
            this.state = Global_Configuration.getCurrentSettings();
            this.w_bar = Go_Wait_Bar.getInstance(100,'Welcome to goGPS', Core.GUI_MODE);  % 0 means text, 1 means GUI, 5 both
            this.sky = Core_Sky.getInstance();
            this.cmd = Command_Interpreter.getInstance;
        end
    end
    
    %% METHODS UI
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_core__
            unique_instance_core__ = [];

            if isempty(unique_instance_core__)
                this = Core();
                unique_instance_core__ = this;
            else
                this = unique_instance_core__;
                this.init();
            end
        end

        function ok_go = openGUI()
            ok_go = gui_goGPS;
        end
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods      
        function init(this)
            this.log.setColorMode(true);
            Core_UI.showTextHeader();
            fclose('all');
        end
        
        function importIniFile(this, ini_settings_file)
            this.state.importIniFile(ini_settings_file);
        end
        
        function prepareProcessing(this)
            this.log.newLine();
            this.log.addMarkedMessage(sprintf('PROJECT: %s', this.state.getPrjName()));

            this.gc.initConfiguration(); % Set up / download observations and navigational files
            this.log.addMarkedMessage('Conjuring all the files!');
            fw = File_Wizard;
            c_mode = this.log.getColorMode();
            this.log.setColorMode(0);
            fw.conjureFiles();
            this.log.setColorMode(c_mode);
        end
    end
    
    %% METHODS RUN
    % ==================================================================================================================================================
    methods      
        function prepareSession(this, session_number)
            %-------------------------------------------------------------------------------------------
            % SESSION START
            %-------------------------------------------------------------------------------------------
            session = session_number;
            
            this.log.newLine;
            this.log.simpleSeparator();
            this.log.addMessage(sprintf('Starting session %d of %d', session, this.state.getSessionCount()));
            this.log.simpleSeparator();
            
            % Init sky
            
            clear rec;  % handle to all the receivers
            for r = 1 : this.state.getRecCount()
                this.log.newLine();
                this.log.addMessage(sprintf('Reading receiver %d of %d', r, this.state.getRecCount()));
                this.log.simpleSeparator();
                
                rec(r) = Receiver(this.state.getConstellationCollector(), this.state.getRecPath(r, session), this.state.getDynMode(r)); %#ok<AGROW>
            end
            this.rec = rec;
            
            % Init sky for this session
            [~, time_lim_large] = rec.getTimeSpan();
            this.sky.initSession(time_lim_large.first, time_lim_large.last);
            this.log.simpleSeparator();
        end
        
        function go(this)
            for s = 1 : this.state.getSessionCount()
                this.prepareSession(s);
                this.cmd.exec(this.rec);
                
                if this.state.isKeepRecList()
                    if isempty(this.rec_list)
                        clear rec_list;
                        rec_list(:,s) = this.rec;
                        this.rec_list = rec_list;
                    else
                        this.rec_list(:,s) = this.rec;
                    end
                end
            end
        end
        
        function exec(this, cmd)
            this.cmd.exec(this.rec, cmd);
        end
    end
        
    %% METHODS UTILITIES
    % ==================================================================================================================================================
    methods
        function [state, log, w_bar] = getUtilities(this)
            state = this.state;
            log = this.log;
            w_bar = this.w_bar;
        end        
    end

    methods % Public Access (Legacy support)
        
    end

end
