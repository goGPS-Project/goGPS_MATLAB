%   CLASS App_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the Processing parameters
%
% EXAMPLE
%   state = App_Settings();
%
% FOR A LIST OF CONSTANTs and METHODS use doc App_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD) 
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, Giulio Taliaferro, ...
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

classdef App_Settings < Settings_Interface

    properties (Constant, Access = 'protected')
    end
    
    % Real constant
    properties(Constant, Access = 'private')
        CUR_INI = 'app_settings.ini';

        LOG_DEFAULT_MODE = 1; % 0 text mode, 1 graphic mode
        LOG_COLOR_MODE = 1;
        
        GUI_DEFAULT_MODE = 'light';
        GUI_DEFAULT_EXPORT_MODE = 'light';

        FLAG_EXPORT_TRANSPARENT = true;
    end
    
    properties (Constant, Access = 'public')
    end

    properties (SetAccess = protected, GetAccess = protected)
        % Location of the current ini file        
        log_default_mode = App_Settings.LOG_DEFAULT_MODE; % 0 text mode, 1 graphic mode
        log_color_mode = App_Settings.LOG_COLOR_MODE;
        
        gui_default_mode = App_Settings.GUI_DEFAULT_MODE;
        gui_default_export_mode = App_Settings.GUI_DEFAULT_EXPORT_MODE;

        flag_export_transparent = App_Settings.FLAG_EXPORT_TRANSPARENT;
        
        flag_deployed_slaves = false;
        bin_path = '../bin/${ARC}/';        
    end

    properties (SetAccess = public, GetAccess = public)
        cur_ini = Core.getFilePath('app_settings');
        creation_time = GPS_Time(now);
    end

    % =========================================================================
    %%  INIT
    % =========================================================================
    methods(Static,Access = private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        
        function this = App_Settings(ini_settings_file)
            % Creator
            %
            % SYNTAX
            %   s_obj = App_Settings(<ini_settings_file>, <prj_home>);

            log = Core.getLogger();
            log.addMarkedMessage('Preparing goGPS settings...');
            log.newLine();
                            
            if nargin == 1 && not(isempty(ini_settings_file))
                this.cur_ini = ini_settings_file;
            end
            
            if (exist(this.cur_ini, 'file') == 2)
                this.import(this.cur_ini);
            else log.addMarkedMessage('Using default settings');
                log.newLine();
            end
        end
    end

    methods (Static)
        function this = getInstance(ini_settings_file)
            % Concrete implementation.  See Singleton superclass.
            persistent unique_instance_gocfg__
            if isempty(unique_instance_gocfg__)
                if nargin == 1
                    this = App_Settings(ini_settings_file);
                else
                    this = App_Settings();
                end
                unique_instance_gocfg__ = this;
            else
                this = unique_instance_gocfg__;
                if nargin == 1
                    this.import(ini_settings_file)
                end
            end
        end
    end

    % =========================================================================
    %%  SET-UP
    % =========================================================================
    
    methods (Static, Access = 'public')
        function setUpSlaves(slaves_folder)
            % Set the binary folder containing slaves
            %
            % SYNTAX
            %   App_Settings.setUpSlaves(slaves_folder)
            app_settings = App_Settings.getInstance;
            if isunix
                % Compiled slaves at the moment works only for unix
                use_deployed = app_settings.useDeployedSlaves;

                app_settings.import();
                if nargin < 1 || isempty(slaves_folder)
                    if use_deployed
                        slaves_folder = app_settings.getBinDir;
                        if ~exist(fullfile(slaves_folder, 'run_startSlave.sh'), 'file')
                            Core.getLogger.addWarning(sprintf('Slave executable not found at "%s", please select the directory containing the binaries', fullfile(slaves_folder, 'run_startSlave.sh')));
                            msg = Core.getMsgGUI();
                            slaves_folder = uigetdir(slaves_folder, 'Select the directory containing binaries files');
                            app_settings.setBinDir(slaves_folder);
                            close(msg.win);
                        end
                        fnp = File_Name_Processor;
                        slaves_folder = fnp.getFullDirPath(slaves_folder, pwd);
                    else
                        slaves_folder = pwd;
                    end
                else
                    app_settings.setBinDir(slaves_folder);
                end
                if ~use_deployed
                    uiwait(errordlg(sprintf(['Set use_deployed in app_settings to use binaries'])));
                else
                    if exist(fullfile(slaves_folder, 'run_startSlave.sh'), 'file')
                        uiwait(msgbox(sprintf('The slave binary have been found at "%s"\goGPS can now work with parallel workers', slaves_folder), 'goGPS'));
                    else
                        uiwait(errordlg(sprintf('The slave binary have not been found at "%s".\nRerun this setup and set the right folder!', slaves_folder), 'goGPS'));
                    end
                end
            else
                uiwait(errordlg(sprintf('Compiled slaves at the moment works only for unix systems')));
            end
        end
    end
    
    % =========================================================================
    %%  INTERFACE REQUIREMENTS
    % =========================================================================
    methods (Access = 'public')
        
        function reload(this)
            this.import();
        end
        
        function importIniFile(this, file_name)
            this.import(file_name);
        end
        
        function import(this, file_name)
            % This function import processing settings from another setting object or ini file
            %
            % SYNTAX
            %   s_obj.import(state)

            if nargin < 2 || isempty(file_name)
                file_name = this.cur_ini;
            end
            
            ini = Ini_Manager(file_name);
            
            try this.log_default_mode = ini.getData('LOG', 'default_mode'); catch, log = Core.getLogger; log.addWarning(sprintf('LOG default_mode not found in app_settings')); end
            try this.log_color_mode = ini.getData('LOG', 'color_mode'); catch, log = Core.getLogger; log.addWarning(sprintf('LOG color_mode not found in app_settings')); end
            
            try this.gui_default_mode = ini.getData('GUI', 'default_mode'); catch, log = Core.getLogger; log.addWarning(sprintf('GUI default_mode not found in app_settings')); end
            try this.gui_default_export_mode = ini.getData('GUI', 'default_export_mode'); catch, log = Core.getLogger; log.addWarning(sprintf('GUI default_export_mode not found in app_settings')); end
            
            try this.flag_export_transparent = ini.getData('EXPORT', 'flag_export_transparent'); catch, log = Core.getLogger; log.addWarning(sprintf('EXPORT flag_export_transparent not found in app_settings')); end
            
            try this.flag_deployed_slaves = logical(ini.getData('PARALLELISM', 'use_deployed_slaves')); catch, log = Core.getLogger; log.addWarning(sprintf('PARALLELISM use_deployed_slaves not found in app_settings')); end
            if isempty(this.flag_deployed_slaves)
                this.flag_deployed_slaves = false;
            end
            
            try this.bin_path = ini.getData('PARALLELISM', 'bin_path'); catch, log = Core.getLogger; log.addWarning(sprintf('PARALLELISM bin_path not found in app_settings')); end

            this.check(); % check after import            
        end
        
        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the this
            %
            % SYNTAX
            %   s_obj.export(str_cell)

            log.addMarkedMessage('Export of goGPS settings - to do...');
        end
        
        function toString(this)
            ini = Ini_Manager(this.cur_ini);
            ini.showData;
        end
    end

    %%  CHECKING FUNCTIONS - to be kept here to access private parameters
    % =========================================================================
    methods (Access = 'protected')
        function checkLogicalField(this, field_name)
            % Check if a logical field of the object is a valid logical number
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % SYNTAX
            %   this.checkLogicalField(string_field_name);
            this.(field_name) = this.checkLogical(field_name, this.(field_name), this.(upper(field_name)));
        end

        function checkCellStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % SYNTAX
            %   this.Field(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            switch nargin
                case 2, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid);
                case 4, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid, check_existence);
                otherwise, error('Settings checkCellStringField called with the wrong number of parameters');
            end
        end

        function is_existing = checkStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % check_existence 0: do not check
            %                 1: check and do nothing
            %                 2: check and                                try                                                                to correct
            %
            % SYNTAX
            %
            %   this.checkStringField(string_field_name, <empty_is_valid == false>, <check_existence == false>);            
            switch nargin
                case 2, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), this.(upper(field_name)));
                case 3, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid);
                case 4, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid, check_existence);
                otherwise, error('Settings checkStringField called with the wrong number of parameters');
            end
        end

        function is_existing = checkPathField(this, field_name, empty_is_valid, check_existence)
            % Check if a string path field of the object is a valid path
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % SYNTAX
            %   this.checkPathField(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            fnp = File_Name_Processor();
            tmp_path = fnp.checkPath(this.(field_name), this.getHomeDir);
            tmp_default_path = fnp.getFullDirPath(fnp.checkPath(this.(upper(field_name)), this.getHomeDir), this.getHomeDir);
            if ~isempty(tmp_path) || empty_is_valid
                this.(field_name) = fnp.getFullDirPath(tmp_path, this.prj_home, [], fnp.getFullDirPath(tmp_default_path));
                switch nargin
                    case 2, [this.(field_name), is_existing] = this.checkString(field_name, tmp_path, tmp_default_path);
                    case 3, [this.(field_name), is_existing] = this.checkString(field_name, tmp_path, tmp_default_path, empty_is_valid);
                    case 4, [this.(field_name), is_existing] = this.checkString(field_name, tmp_path, tmp_default_path, empty_is_valid, check_existence);
                    otherwise, error('Settings checkStringField called with the wrong number of parameters');
                end
                this.(field_name) = fnp.getFullDirPath(this.(field_name), this.getHomeDir);
            else
                this.(field_name) = tmp_default_path;
            end
        end

        function checkNumericField(this, field_name, limits, valid_val)
            % Check if a numeric field of the object is valid
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % SYNTAX
            %   this.checkNumericField(string_field_name, <limits>, <valid_values>);
            switch nargin
                case 2, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits);
                case 4, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits, valid_val);
                otherwise, error('Settings checkNumericField called with the wrong number of parameters');
            end
        end        
    end
   
    % =========================================================================
    %%  TEST PARAMETERS VALIDITY
    % =========================================================================

    methods (Access = 'public')
        function check(this)
            % Check the validity of the fields
            %
            % SYNTAX
            %   this.check();

            this.checkNumericField('log_default_mode', [0 1]);
            this.checkLogicalField('log_color_mode');
            this.checkLogicalField('flag_export_transparent');   
            
            if isdeployed
                this.log_color_mode = false;
            end
        end        
    end

    % =========================================================================
    %%  GETTERS
    % =========================================================================
    methods  
        function out = getLogMode(this)
            out = this.log_default_mode;
        end
        
        function out = isLogColorMode(this)
            out = logical(this.log_color_mode);
        end
        
        function out = getGUIMode(this)
            out = this.gui_default_mode;
        end
        
        function out = getGUIModeExport(this)
            out = this.gui_default_export_mode;
        end
        
        function out = isExportTransparent(this)
            out = logical(this.flag_export_transparent);
        end
        
        function out = useDeployedSlaves(this)
            out = this.flag_deployed_slaves || isdeployed;
        end
        
        function out = getBinDir(this)
            % Get the dir containing binary files
            %
            % SYNTAX
            %   bin_dir = this.getBinDir()
            out = strrep(this.bin_path, '${ARC}',computer('arch'));
        end
               
        function setBinDir(this, bin_path)
            % Set the dir containing binary files
            % Overwrite cur_ini file to save the setting
            %
            % SYNTAX
            %   this.setBinDir(bin_dir)
            this.bin_path = bin_path;
            % modify the settings file
            % read the file
            txt = Core_Utils.readTextFile(this.cur_ini);
            bin_path_id = strfind(txt,'bin_path = "');
            if not(isempty(bin_path_id))
                bin_path_id = bin_path_id(1) + 11;
                end_id = bin_path_id + regexp(txt((bin_path_id + 1):end),'"','once');
                if not(isempty(end_id))
                    txt = [txt(1:bin_path_id) bin_path txt(end_id:end)];
                    fid = fopen(this.cur_ini, 'wt');
                    if fid > 0
                        fwrite(fid, txt);
                        fclose(fid);
                    end
                end
            end
        end
    end  
    %%  SETTERS
    % =========================================================================
    methods
        function out = setLogMode(this, val)
            this.log_default_mode = val;
        end
        
        function out = setLogColorMode(this, val)
            this.log_color_mode = logical(val);
        end
        
        function out = setModeGUI(this, val)
            this.gui_default_mode = val;
        end
        
        function out = setGUIModeExport(this, val)
            this.gui_default_export_mode = val;
        end
        
        function out = setExportTransparent(this, val)
            this.flag_export_transparent = logical(val);
        end
    end
    
end
