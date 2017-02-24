%   CLASS IO_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the file_paths and folders
%
% EXAMPLE
%   settings = IO_Settings();
%
% FOR A LIST OF CONSTANTS and METHODS use doc IO_Settings

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___ 
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.9.1
% 
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
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

classdef IO_Settings < Settings_Interface
    
    properties (Constant, Access = 'protected')
        % id to string of out modes
        OUT_MODE = {'0: old goGPS naming - each file saved with run number in out_dir', ...
                    '1: run naming - each file saved in a folder containing the run number in out_dir' ...
                    '2... other values for future implementation, e.g. each output in a folder with a certain format (doy_hh-hh)'}
                
        DEFAULT_DIR_IN = ['..' filesep 'data' filesep];
        DEFAULT_DIR_OUT = ['..' filesep 'data' filesep];  
        
        % Location of the latest project (the ini contains just a reference to the default project ini file - that is actually a settings file
        last_prj_ini_path = [IO_Settings.DEFAULT_DIR_IN 'last_prj.ini'];
    end
    
    properties (SetAccess = private, GetAccess = public)
        %------------------------------------------------------------------
        % PROJECT
        %------------------------------------------------------------------        

        % Name of the project
        prj_name = 'Default PPP project';
                
        % Location of the project <relative path from goGPS folder>
        prj_home = [IO_Settings.DEFAULT_DIR_IN 'default_ppp_prj' filesep];
        
        cur_ini = [IO_Settings.DEFAULT_DIR_IN 'default_PPP_prj' filesep 'Config' filesep 'config.ini'];        
        
        %------------------------------------------------------------------
        % DEPRECATE
        %------------------------------------------------------------------
        % deprecate INI - it contains some additional setting (yet not imported in the new settings system)
        deprecate_ini_path = '../data/deprecate_settings/testTropo_ZIMM_InputFiles.ini';

        %------------------------------------------------------------------
        % RECEIVERS
        %------------------------------------------------------------------

        % Observation files of the Receivers
        %receiver_rinex_file = Receiver_File;        

        %------------------------------------------------------------------
        % SATELLITES
        %------------------------------------------------------------------
        
        % Path to Navigational files folder
        nav_dir = [IO_Settings.DEFAULT_DIR_IN 'SAT' filesep 'NAV' filesep];
        
        % Path to Clock Offset files folder
        clk_dir = [IO_Settings.DEFAULT_DIR_IN 'SAT' filesep 'CLK' filesep];
        
        % Path to CRX folder containing files of Satellites problems
        crx_dir = [IO_Settings.DEFAULT_DIR_IN 'SAT' filesep 'CRX' filesep];

        % Path to DCB folder containing files of Differential Code Biases
        dcb_dir = [IO_Settings.DEFAULT_DIR_IN 'SAT' filesep 'DCB' filesep];

        %------------------------------------------------------------------
        % REFERENCE
        %------------------------------------------------------------------
        
        % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        geoid_dir = [IO_Settings.DEFAULT_DIR_IN 'geoid' filesep];
        
        %------------------------------------------------------------------
        % DTM (SET PATH AND LOAD PARAMETER FILES)
        %------------------------------------------------------------------

        % Path to DTM folder containing DTM files
        dtm_dir = [IO_Settings.DEFAULT_DIR_IN 'dtm' filesep];
        
        %------------------------------------------------------------------
        % UI IMAGES
        %------------------------------------------------------------------

        % Path to images used by the interface
        img_dir = [IO_Settings.DEFAULT_DIR_IN 'img' filesep];
        
        % Location of the goGPS logo 64x64
        img_logo64;        
        
        %------------------------------------------------------------------
        % OUTPUT
        %------------------------------------------------------------------
        
        % Output style
        out_style = 0;
        %  - out_style = 0: old goGPS naming: each file saved with run number in out_dir 
        %  - out_style = 1: run naming: each file saved in a folder containing the run number in out_dir 
        %  - out_style = 2... other values for future implementation, e.g. each output in a folder with a certain format (doy_hh-hh)
                
        % Directory containing the output of the project
        out_dir = 'out'; % location relative to the project home
        
        % Every time a solution is computed a folder with prefix followed by the run number is created
        prefix = 'run_'
        
        % This parameter store the current run number
        run_counter = 0;
    end
            
    % =========================================================================
    %  INIT
    % =========================================================================
    methods
        function this = IO_Settings()
            % Creator of IO_settings - verbosity level (true/false) can be set or ini file
            this.img_logo64 = [this.img_dir 'goGPS_logo_64.png'];
            if exist(this.last_prj_ini_path, 'file')
                ini = Ini_Manager(this.last_prj_ini_path);
                this.cur_ini = ini.getData('CurrentProject', 'last_prj_ini');
            end
        end                
    end
    
    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods
        function import(this, settings)
            % This function import IO (only) settings from another setting object
            if isa(settings, 'Ini_Manager')
                % PROJECT
                this.prj_name   = settings.getData('prj_name');
                this.prj_home   = settings.getData('prj_home');
                %this.cur_ini    = settings.getData('cur_ini');
                % DEPRECATE
                this.deprecate_ini_path = settings.getData('deprecate_ini_path');
                % RECEIVERS
                % SATELLITES
                this.nav_dir    = settings.getData('nav_dir');
                this.clk_dir    = settings.getData('clk_dir');
                this.crx_dir    = settings.getData('crx_dir');
                this.dcb_dir    = settings.getData('dcb_dir');
                % REFERENCE
                this.geoid_dir  = settings.getData('geoid_dir');
                this.dtm_dir    = settings.getData('dtm_dir');
                % UI
                this.img_dir    = settings.getData('img_dir');
                this.img_logo64 = settings.getData('img_logo64');
                % OUTPUT
                this.out_style  = settings.getData('out_style');
                this.out_dir = settings.getData('out_dir');
                this.prefix = settings.getData('prefix');
                this.run_counter = settings.getData('run_counter');
            else
                % PROJECT
                this.prj_name   = settings.prj_name;
                this.prj_home   = settings.prj_home;
                %this.cur_ini   = settings.cur_ini;
                % DEPRECATE
                this.deprecate_ini_path = settings.deprecate_ini_path;
                % RECEIVERS
                % SATELLITES
                this.nav_dir    = settings.nav_dir;
                this.clk_dir    = settings.clk_dir;
                this.crx_dir    = settings.crx_dir;
                this.dcb_dir    = settings.dcb_dir;
                % REFERENCE
                this.geoid_dir  = settings.geoid_dir;
                this.dtm_dir    = settings.dtm_dir;
                % UI
                this.img_dir    = settings.img_dir;
                this.img_logo64 = settings.img_logo64;
                % OUTPUT
                this.out_style  = settings.out_style;
                this.out_dir = settings.out_dir;
                this.prefix = settings.prefix;
                this.run_counter = settings.run_counter;                
            end
        end
        
        function str = toString(this, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end            
            str = [str '---- PROJECT --------------------------------------------------------------' 10 10];
            str = [str sprintf(' Project name:                                     %s\n', this.prj_name)];
            str = [str sprintf(' Project home:                                     %s\n', this.prj_home)];
            str = [str sprintf(' Path to the current project ini file:             %s\n\n', this.cur_ini)];
            str = [str '---- DEPRECATE ------------------------------------------------------------' 10 10];
            str = [str sprintf(' Deprecate ini (of additional parameters):         %s\n\n', this.deprecate_ini_path)];            
            str = [str '---- INPUT FOLDERS: SATELLITE ---------------------------------------------' 10 10];
            str = [str sprintf(' Directory of Navigational Files:                  %s\n', this.nav_dir)];
            str = [str sprintf(' Directory of Satellite clock offsets:             %s\n', this.clk_dir)];
            str = [str sprintf(' Directory of CRX (satellite problems):            %s\n', this.crx_dir)];
            str = [str sprintf(' Directory of DCB (Differential Code Biases):      %s\n\n', this.dcb_dir)];
            str = [str '---- INPUT FOLDERS: REFERENCE ---------------------------------------------' 10 10];
            str = [str sprintf(' Directory of Geoid models:                        %s\n', this.geoid_dir)];
            str = [str sprintf(' Directory of DTM data:                            %s\n\n', this.dtm_dir)];
            str = [str '---- INPUT FOLDERS: UI ----------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of images for UI:                       %s\n', this.img_dir)];
            str = [str sprintf('  - Image of the goGPS logo:                       %s\n\n', this.img_logo64)];
            str = [str '---- OUTPUT SETTINGS ------------------------------------------------------' 10 10];
            str = [str sprintf(' Output style %s\n\n', this.OUT_MODE{this.out_style + 1})];
            str = [str sprintf(' Directory containing the output of the project:   %s\n', this.out_dir)];
            str = [str sprintf(' Prefix of each run:                               %s\n', this.prefix)];
            str = [str sprintf(' Run counter:                                      %d\n\n', this.run_counter)];
        end
        
        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the this            
            if (nargin == 1)
                str_cell = {};
            end
            % PROJECT
            str_cell = Ini_Manager.toIniStringSection('PROJECT', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of the project', str_cell);
            str_cell = Ini_Manager.toIniString('prj_name', this.prj_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Home of the project', str_cell);
            str_cell = Ini_Manager.toIniString('prj_home', this.prj_home, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % DEPRECATE
            str_cell = Ini_Manager.toIniStringSection('DEPRECATE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Deprecate ini - path - it contains some additional setting (yet not imported in the new settings system)', str_cell);
            str_cell = Ini_Manager.toIniString('deprecate_ini_path', this.deprecate_ini_path, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % RECEIVERS
            % SATELLITES
            str_cell = Ini_Manager.toIniStringSection('INPUT_SAT', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Navigational files', str_cell);
            str_cell = Ini_Manager.toIniString('nav_dir', this.nav_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of clock offset files', str_cell);
            str_cell = Ini_Manager.toIniString('clk_dir', this.clk_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of CRX files (containing satellite problems)', str_cell);
            str_cell = Ini_Manager.toIniString('crx_dir', this.crx_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DCB files (Differential Code Biases)', str_cell);
            str_cell = Ini_Manager.toIniString('dcb_dir', this.dcb_dir, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % REFERENCE
            str_cell = Ini_Manager.toIniStringSection('INPUT_REF', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Geoid files', str_cell);
            str_cell = Ini_Manager.toIniString('geoid_dir', this.geoid_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DTM data', str_cell);
            str_cell = Ini_Manager.toIniString('dtm_dir', this.dtm_dir, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % UI
            str_cell = Ini_Manager.toIniStringSection('INPUT_UI', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of images for UI', str_cell);
            str_cell = Ini_Manager.toIniString('img_dir', this.img_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Path to the image of the logo 64x64 px', str_cell);
            str_cell = Ini_Manager.toIniString('img_logo64', this.img_logo64, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % OUTPUT
            str_cell = Ini_Manager.toIniStringSection('OUTPUT', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Output style', str_cell);
            str_cell = Ini_Manager.toIniString('out_style', this.out_style, str_cell);
            for i = 1 : numel(this.OUT_MODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.OUT_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('out_style = 0 -> Directory containing the basic output folder of the project - path relative to goGPS folder', str_cell);
            str_cell = Ini_Manager.toIniStringComment('out_style = 1 -> Directory containing the basic output folder of the project - path relative to home folder', str_cell);
            str_cell = Ini_Manager.toIniString('out_dir', this.out_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('out_style = 0 -> Every time a solution is computed a file with prefix followed by the run number is created', str_cell);
            str_cell = Ini_Manager.toIniStringComment('out_style = 1 -> Every time a solution is computed a folder with prefix followed by the run number is created and the outputs are stored into it', str_cell);
            str_cell = Ini_Manager.toIniString('prefix', this.prefix, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Current run number', str_cell);
            str_cell = Ini_Manager.toIniString('run_counter', this.run_counter, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);            
        end
    end    
    
    % =========================================================================
    %  LEGACY IMPORT
    % =========================================================================
    methods (Access = 'public')
        function legacyImport(this, state)
            % import from the state variable (saved into the old interface mat file of goGPS)
            % If a group of imports fails display a warning but continue the
            % import of other groups
            % SYNTAX: this.legacyImport(state)
            
            % RECEIVER DEFAULT PARAMETERS ---------------------------------            
            try 
                this.deprecate_ini_path = state.INIsettings;
                this.out_dir = state.gogps_data_output;
                this.prefix = state.gogps_data_output_prefix;
            catch ex
                this.logger.addWarning(['Legacy import "IO file / folders" failed - ', ex.message])
            end
        end
    end
    
    % =========================================================================
    %  GETTERS
    % =========================================================================
    methods
        function dir = getImgDir(this)
            dir = this.img_dir;
        end
        
        function dir = getLogo(this)
            dir = this.img_logo64;
        end
        
    end
    
    % =========================================================================
    %  TEST
    % =========================================================================
    methods (Static, Access = 'public')
        function test()
            % test the class
            s = IO_Settings();
            s.testInterfaceRoutines();
        end
    end
end
