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
%    |___/                    v 0.5.0
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
    
    % Default values for each field - useful to restore corrupted field
    properties (Constant, Access = 'protected')
        
        % PROJECT                
        PRJ_NAME = 1;  % Name of the project
        PRJ_HOME = [IO_Settings.DEFAULT_DIR_IN 'project' filesep 'default_DD' filesep]; % Location of the project <relative path from goGPS folder> 
        CUR_INI = [IO_Settings.DEFAULT_DIR_IN 'project' filesep 'default_DD' filesep 'Config' filesep 'settings.ini']; % Location of the current ini file
        % DEPRECATE                
        INPUT_FILE_INI_PATH = '../data/project/default_DD/config/inputFiles.ini'; % deprecate INI - it contains some additional setting (yet not imported in the new settings system)
        % RECEIVERS
        % SATELLITES                
        EPH_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'EPH' filesep]; % Path to Ephemeris files folder
        CLK_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'CLK' filesep]; % Path to Clock Offset files folder
        CRX_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'CRX' filesep]; % Path to CRX folder containing files of Satellites problems
        DCB_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'DCB' filesep]; % Path to DCB folder containing files of Differential Code Biases
        EMS_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'SBAS' filesep 'EMS' filesep]; % Path to EMS folder containing files of EGNOS Message Server.
        % STATIONS
        CRD_DIR = [IO_Settings.DEFAULT_DIR_IN 'stations' filesep 'CRD' filesep]; % Path to Ephemeris files folder
        MET_DIR = [IO_Settings.DEFAULT_DIR_IN 'stations' filesep 'MET' filesep]; % Path to Clock Offset files folder
        OCEAN_DIR = [IO_Settings.DEFAULT_DIR_IN 'stations' filesep 'ocean' filesep]; % Path to CRX folder containing files of Satellites problems
        % REFERENCE        
        GEOID_DIR = [IO_Settings.DEFAULT_DIR_IN 'reference' filesep 'geoid' filesep]; % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs        
        % DTM (SET PATH AND LOAD PARAMETER FILES)        
        DTM_DIR = [IO_Settings.DEFAULT_DIR_IN 'reference' filesep 'DTM' filesep]; % Path to DTM folder containing DTM files
        % UI IMAGES
        IMG_DIR = [IO_Settings.DEFAULT_DIR_IN 'img' filesep];  % Path to images used by the interface
        IMG_LOGO64 = [IO_Settings.IMG_DIR 'goGPS_logo_64.png'];     % Location of the goGPS logo 64x64
        OUT_STYLE = 0;       % Output style
                             %  - out_style = 0: old goGPS naming: each file saved with run number in out_dir 
                             %  - out_style = 1: run naming: each file saved in a folder containing the run number in out_dir 
                             %  - out_style = 2... other values for future implementation, e.g. each output in a folder with a certain format (doy_hh-hh)
        OUT_DIR = [IO_Settings.DEFAULT_DIR_OUT  'project' filesep 'default_DD' filesep 'out' filesep]; % Directory containing the output of the project
        OUT_PREFIX = 'out';  % Every time a solution is computed a folder with prefix followed by the run number is created
        RUN_COUNTER = 0;     % This parameter store the current run number
    end
    
    properties (Constant, Access = 'protected')
        % id to string of out modes
        OUT_MODE = {'0: old goGPS naming - each file saved with run number in out_dir', ...
                    '1: run naming - each file saved in a folder containing the run number in out_dir' ...
                    '2... other values for future implementation, e.g. each output in a folder with a certain format (doy_hh-hh)'}                
        DEFAULT_DIR_IN = ['..' filesep 'data' filesep];
        DEFAULT_DIR_OUT = ['..' filesep 'data' filesep];  
    end
    
    properties (Constant, Access = 'public')
        % Location of the latest project (the ini contains just a reference to the default project ini file - that is actually a settings file
        LAST_SETTINGS = [IO_Settings.DEFAULT_DIR_IN 'last_settings.ini'];
    end
    
    properties (SetAccess = protected, GetAccess = public)
        %------------------------------------------------------------------
        % PROJECT
        %------------------------------------------------------------------        

        % Name of the project
        prj_name = IO_Settings.PRJ_NAME;
                
        % Location of the project <relative path from goGPS folder>
        prj_home = IO_Settings.PRJ_HOME;
        
        % Location of the current ini file
        cur_ini = IO_Settings.CUR_INI;
        
        %------------------------------------------------------------------
        % DEPRECATE
        %------------------------------------------------------------------
        % deprecate INI - it contains some additional setting (yet not imported in the new settings system)
        input_file_ini_path = IO_Settings.INPUT_FILE_INI_PATH;

        %------------------------------------------------------------------
        % RECEIVERS
        %------------------------------------------------------------------

        % Observation files of the Receivers
        % receiver_rinex_file = Receiver_File;        

        %------------------------------------------------------------------
        % SATELLITES
        %------------------------------------------------------------------
        
        % Path to Ephemeris files folder
        eph_dir = IO_Settings.EPH_DIR;
        % Path to Clock Offset files folder
        clk_dir = IO_Settings.CLK_DIR;
        % Path to CRX folder containing files of Satellites problems
        crx_dir = IO_Settings.CRX_DIR;
        % Path to DCB folder containing files of Differential Code Biases
        dcb_dir = IO_Settings.DCB_DIR;
        % Path to EMS folder containing files of EGNOS Message Server.
        ems_dir = IO_Settings.EMS_DIR;

        %------------------------------------------------------------------
        % STATIONS
        %------------------------------------------------------------------

        % Path to stations coordinates files
        crd_dir = IO_Settings.CRD_DIR;
        % Path to stations metereological files
        met_dir = IO_Settings.MET_DIR;
        % Path to stations ocean loading files
        ocean_dir = IO_Settings.OCEAN_DIR;

        %------------------------------------------------------------------
        % REFERENCE
        %------------------------------------------------------------------
        
        % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        geoid_dir = IO_Settings.GEOID_DIR;
        
        %------------------------------------------------------------------
        % DTM (SET PATH AND LOAD PARAMETER FILES)
        %------------------------------------------------------------------

        % Path to DTM folder containing DTM files
        dtm_dir = IO_Settings.DTM_DIR;
        
        %------------------------------------------------------------------
        % UI IMAGES
        %------------------------------------------------------------------

        % Path to images used by the interface
        img_dir = IO_Settings.IMG_DIR;
        
        % Location of the goGPS logo 64x64
        img_logo64 = IO_Settings.IMG_LOGO64;      
        
        %------------------------------------------------------------------
        % OUTPUT
        %------------------------------------------------------------------
        
        % Output style
        out_style = IO_Settings.OUT_STYLE;
        %  - out_style = 0: old goGPS naming: each file saved with run number in out_dir 
        %  - out_style = 1: run naming: each file saved in a folder containing the run number in out_dir 
        %  - out_style = 2... other values for future implementation, e.g. each output in a folder with a certain format (doy_hh-hh)
                
        % Directory containing the output of the project
        out_dir = IO_Settings.OUT_DIR;
        
        % Every time a solution is computed a folder with prefix followed by the run number is created
        out_prefix = IO_Settings.OUT_PREFIX;
        
        % This parameter store the current run number
        run_counter = IO_Settings.RUN_COUNTER;
    end
            
    % =========================================================================
    %  INIT
    % =========================================================================
    methods
        function this = IO_Settings()
            % Creator of IO_settings - verbosity level (true/false) can be set or ini file
            this.img_logo64 = [this.img_dir 'goGPS_logo_64.png'];
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
                this.prj_name   = checkPath(settings.getData('prj_name'));
                this.prj_home   = checkPath(settings.getData('prj_home'));
                %this.cur_ini    = settings.getData('cur_ini');
                % DEPRECATE
                this.input_file_ini_path = checkPath(settings.getData('input_file_ini_path'));
                % RECEIVERS
                % SATELLITES
                this.eph_dir    = checkPath(settings.getData('eph_dir'));
                this.clk_dir    = checkPath(settings.getData('clk_dir'));
                this.crx_dir    = checkPath(settings.getData('crx_dir'));
                this.dcb_dir    = checkPath(settings.getData('dcb_dir'));
                this.ems_dir    = checkPath(settings.getData('ems_dir'));
                % STATIONS
                this.crd_dir    = checkPath(settings.getData('crd_dir'));
                this.met_dir    = checkPath(settings.getData('met_dir'));
                this.ocean_dir  = checkPath(settings.getData('ocean_dir'));
                % REFERENCE
                this.geoid_dir  = checkPath(settings.getData('geoid_dir'));
                this.dtm_dir    = checkPath(settings.getData('dtm_dir'));
                % UI
                this.img_dir    = checkPath(settings.getData('img_dir'));
                this.img_logo64 = checkPath(settings.getData('img_logo64'));
                % OUTPUT
                this.out_style  = settings.getData('out_style');
                this.out_dir = checkPath(settings.getData('out_dir'));
                this.out_prefix = checkPath(settings.getData('out_prefix'));
                this.run_counter = settings.getData('run_counter');
            else
                % PROJECT
                this.prj_name   = settings.prj_name;
                this.prj_home   = settings.prj_home;
                %this.cur_ini   = settings.cur_ini;
                % DEPRECATE
                this.input_file_ini_path = settings.input_file_ini_path;
                % RECEIVERS
                % SATELLITES
                this.eph_dir    = settings.eph_dir;
                this.clk_dir    = settings.clk_dir;
                this.crx_dir    = settings.crx_dir;
                this.dcb_dir    = settings.dcb_dir;
                this.ems_dir    = settings.ems_dir;
                % STATIONS
                this.crd_dir    = settings.crd_dir;
                this.met_dir    = settings.met_dir;
                this.ocean_dir    = settings.ocean_dir;
                % REFERENCE
                this.geoid_dir  = settings.geoid_dir;
                this.dtm_dir    = settings.dtm_dir;
                % UI
                this.img_dir    = settings.img_dir;
                this.img_logo64 = settings.img_logo64;
                % OUTPUT
                this.out_style  = settings.out_style;
                this.out_dir = settings.out_dir;
                this.out_prefix = settings.out_prefix;
                this.run_counter = settings.run_counter;                
            end
            this.check();
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
            str = [str sprintf(' Deprecate ini (of additional parameters):         %s\n\n', this.input_file_ini_path)];            
            str = [str '---- INPUT FOLDERS: SATELLITE ---------------------------------------------' 10 10];
            str = [str sprintf(' Directory of Ephemeris files:                     %s\n', this.eph_dir)];
            str = [str sprintf(' Directory of Satellite clock offsets:             %s\n', this.clk_dir)];
            str = [str sprintf(' Directory of CRX (satellite problems):            %s\n', this.crx_dir)];
            str = [str sprintf(' Directory of DCB (Differential Code Biases):      %s\n', this.dcb_dir)];
            str = [str sprintf(' Directory of EMS (EGNOS Message Server):          %s\n\n', this.ems_dir)];
            str = [str '---- INPUT FOLDERS: STATIONS ----------------------------------------------' 10 10];
            str = [str sprintf(' Directory of coordinates files:                   %s\n', this.crd_dir)];
            str = [str sprintf(' Directory of metereological data:                 %s\n', this.met_dir)];
            str = [str sprintf(' Directory of ocean loading files:                 %s\n', this.ocean_dir)];
            str = [str '---- INPUT FOLDERS: REFERENCE ---------------------------------------------' 10 10];
            str = [str sprintf(' Directory of Geoid models:                        %s\n', this.geoid_dir)];
            str = [str sprintf(' Directory of DTM data:                            %s\n\n', this.dtm_dir)];
            str = [str '---- INPUT FOLDERS: UI ----------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of images for UI:                       %s\n', this.img_dir)];
            str = [str sprintf('  - Image of the goGPS logo:                       %s\n\n', this.img_logo64)];
            str = [str '---- OUTPUT SETTINGS ------------------------------------------------------' 10 10];
            str = [str sprintf(' Output style %s\n\n', this.OUT_MODE{this.out_style + 1})];
            str = [str sprintf(' Directory containing the output of the project:   %s\n', this.out_dir)];
            str = [str sprintf(' Prefix of each run:                               %s\n', this.out_prefix)];
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
            str_cell = Ini_Manager.toIniString('input_file_ini_path', this.input_file_ini_path, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % RECEIVERS
            % SATELLITES
            str_cell = Ini_Manager.toIniStringSection('INPUT_SAT', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Ephemeris files', str_cell);
            str_cell = Ini_Manager.toIniString('eph_dir', this.eph_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of clock offset files', str_cell);
            str_cell = Ini_Manager.toIniString('clk_dir', this.clk_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of CRX files (containing satellite problems)', str_cell);
            str_cell = Ini_Manager.toIniString('crx_dir', this.crx_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DCB files (Differential Code Biases)', str_cell);
            str_cell = Ini_Manager.toIniString('dcb_dir', this.dcb_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of EMS files (EGNOS Message Server).', str_cell);
            str_cell = Ini_Manager.toIniString('ems_dir', this.ems_dir, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % STATIONS
            str_cell = Ini_Manager.toIniStringSection('INPUT_STATIONS', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of coordinates files', str_cell);
            str_cell = Ini_Manager.toIniString('crd_dir', this.crd_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of metereological data', str_cell);
            str_cell = Ini_Manager.toIniString('met_dir', this.met_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of ocean loading files', str_cell);
            str_cell = Ini_Manager.toIniString('ocean_dir', this.ocean_dir, str_cell);
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
            str_cell = Ini_Manager.toIniString('out_prefix', this.out_prefix, str_cell);
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
                this.input_file_ini_path = state.INIsettings;
                this.out_dir = state.gogps_data_output;
                this.out_prefix = state.gogps_data_output_prefix;
            catch ex
                this.logger.addWarning(['Legacy import "IO file / folders" failed - ', ex.message])
            end
            this.check();
        end
    end
    
    % =========================================================================
    %  GETTERS
    % =========================================================================
    methods
        function dir = getImgDir(this)
            % Get the directory of UI images
            % SYNTAX: dir = this.getImgDir()
            dir = this.img_dir;
        end
        
        function dir = getLogo(this)
            % Get the logo of goGPS to be shown on the UI
            % SYNTAX: dir = this.getLogo()
            dir = this.img_logo64;
        end        
    end
    
    % =========================================================================
    %  SETTERS
    % =========================================================================
    methods
        function setDeprecateIniPath(this, new_path)
            % Get the directory of UI images
            % SYNTAX: dir = this.setDeprecateIniPath(new_path)
            this.input_file_ini_path = new_path;
        end
        
        function setFilePath(this, file_path)
            % Set the file name of the current settings
            % SYNTAX: this.setFilePath(file_path)
            [path_str, name, ~] = fileparts(file_path);
            this.cur_ini = [path_str filesep name '.ini'];
        end
        
        function updatePrj(this, file_path)
            % Set the file project name / home / file_path from file_path
            % SYNTAX: this.autoUpdatePrj(this, file_path)
            [path_str, name, ~] = fileparts(checkPath(file_path));
            this.cur_ini = [path_str filesep name '.ini'];
            path_parts = strsplit(path_str,filesep);
            if numel(path_parts) > 3
                this.prj_home = fullfile(path_parts{1:end-1}, filesep);
                this.prj_name = path_parts{end-1};
                this.logger.addMessage('Trying to guess project name / home / ini');
                this.logger.addMessage(sprintf(' name: %s', this.prj_name));
                this.logger.addMessage(sprintf(' home: %s', this.prj_home));
                this.logger.addMessage(sprintf(' ini:  %s', this.cur_ini));                
            end            
        end
    end

    % =========================================================================
    %  TEST PARAMETERS VALIDITY
    % =========================================================================    
    methods (Access = 'protected')
        function checkLogicalField(this, field_name)
            % Check if a logical field of the object is a valid logical number
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.checkLogicalField(string_field_name);
            this.(field_name) = this.checkLogical(field_name, this.(field_name), this.(upper(field_name)));
        end

        function checkStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.checkStringField(string_field_name, <empty_is_valid == false>, <check_existence == false>);            
            switch nargin
                case 2, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid);
                case 4, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid, check_existence);
                otherwise, error('Settings checkStringField called with the wrong number of parameters');
            end
        end
        
        function checkNumericField(this, field_name, limits, valid_val)
            % Check if a numeric field of the object is valid
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.checkNumericField(string_field_name, <limits>, <valid_values>);
            switch nargin
                case 2, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits);
                case 4, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits, valid_val);
                otherwise, error('Settings checkNumericField called with the wrong number of parameters');
            end
        end
    end
    
    % =========================================================================
    %  TEST PARAMETERS VALIDITY
    % =========================================================================       
    
    methods (Access = 'public')
        function check(this)            
            % Check the validity of the fields
            % SYNTAX: this.check();
            
            this.checkStringField('prj_name', false);
            this.checkStringField('prj_home', false, true);
            this.checkStringField('cur_ini', false);

            this.checkStringField('input_file_ini_path', false);
            
            this.checkStringField('eph_dir', false, true);
            this.checkStringField('clk_dir', false, true);
            this.checkStringField('crx_dir', false);
            this.checkStringField('dcb_dir', false);
            this.checkStringField('ems_dir', true);

            this.checkStringField('crd_dir', false);
            this.checkStringField('met_dir', false);
            this.checkStringField('ocean_dir', false);

            this.checkStringField('geoid_dir', true);
            this.checkStringField('dtm_dir', true);

            this.checkStringField('img_dir', false, true);
            this.checkStringField('img_logo64', false, true);

            this.checkNumericField('out_style',[0 numel(this.OUT_MODE)-2]);
            this.checkStringField('out_prefix', true);
            this.checkNumericField('run_counter',[0 1e6]);
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
