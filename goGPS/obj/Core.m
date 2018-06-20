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
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
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
        GO_GPS_VERSION = '0.6.0 alpha 2 - nightly';
        GUI_MODE = 0; % 0 means text, 1 means GUI, 5 both
    end

    %% PROPERTIES SINGLETON POINTERS
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons
        log         % Logger
        gc          % global configuration
        state       % state (gc.state)
        w_bar       % WaitBar
        sky         % Core_Sky 
        cmd         % Command_Interpreter
    end

    %% PROPERTIES RECEIVERS
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons
        cur_session     % id of the current session
        
        rin_list        % List of observation file (as File_Rinex objects) to store minimal information on the input files
        met_list        % List of meteorological file (as File_Rinex objects) to store minimal information on the input files
        
        rec             % List of all the receiver used
    end

    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = private)
        % Concrete implementation.  See Singleton superclass.
        function this = Core(force_clean)
            if nargin < 1
                force_clean = false;
            end
            % Core object creator
            this.log = Logger.getInstance();
            this.init(force_clean);            
        end
    end
    
    %% METHODS UI
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = getInstance(force_clean)
            if nargin < 1
                force_clean = false;
            end
            % Get the persistent instance of the class
            persistent unique_instance_core__

            if isempty(unique_instance_core__)
                this = Core(force_clean);
                unique_instance_core__ = this;
            else
                this = unique_instance_core__;
                this.init(force_clean);
            end            
        end

        function ok_go = openGUI()
            ok_go = gui_goGPS;
        end
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this, force_clean)
            % Get instances for:
            %   - Global_Configuration
            %   - Settings
            %   - Wait_Bar
            %   - Sky
            %   - Command_Interpreter
            %
            % SYNTAX:
            %   this.init(<force_clean>)
            
            if nargin < 2
                force_clean = false;
            end
            this.log.setColorMode(true);
            Core_UI.showTextHeader();
            fclose('all');
            this.gc = Global_Configuration.getInstance();
            this.state = Global_Configuration.getCurrentSettings();
            this.w_bar = Go_Wait_Bar.getInstance(100,'Welcome to goGPS', Core.GUI_MODE);  % 0 means text, 1 means GUI, 5 both
            this.sky = Core_Sky.getInstance(force_clean);
            this.cmd = Command_Interpreter.getInstance(this);
            if force_clean
                this.cur_session = 0;
                this.rec = [];
            end
        end
        
        function import(this, state)
            % Import Settings from ini files
            % 
            % INPUT:
            %   state       can be a string to the path of the settings file
            %               or a state object
            %
            % SYNTAX:
            %   this.import(state)
            
            if ischar(state)
                this.importIniFile(state);
            else
                this.importState(state);
            end
        end
            
        function importIniFile(this, ini_settings_file)
            % Import Settings from ini files
            %
            % SYNTAX:
            %   this.importIniFile(ini_settings_file)
            if  exist(ini_settings_file, 'file')
                this.importIniFile(ini_settings_file);
            end
        end
        
        function importState(this, state)
            % Import Settings from state
            %
            % SYNTAX:
            %   this.importState(state)
            
            this.state.import(state);
        end
        
        function prepareProcessing(this, flag_rem_check)
            % Init settings, and download necessary files
            %
            % SYNTAX:
            %   this.prepareProcessing(flag_rem_check)
            if nargin == 2
                this.state.setRemCheck(flag_rem_check);
            end
            
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
            % Check the time-limits for the files in the session
            % Init the Sky and Meteo object
            %
            % SYNTAX
            %   this.prepareSession(session_number)
            
            this.cur_session = session_number;
            session = session_number;
            
            this.log.newLine;
            this.log.simpleSeparator();
            this.log.addMessage(sprintf('Starting session %d of %d', session, this.state.getSessionCount()));
            this.log.simpleSeparator();
            
            if isempty(this.rec)
                rec = GNSS_Station;
            else
                rec = this.rec;
            end
            this.log.newLine();
            if ~this.state.isRinexSession()
                rin_list = this.getRinFileList();
                if ~isempty(rin_list)
                    [~, out_limits] = this.state.getSessionLimits(session);
                    time_lim_large = out_limits;
                end
            else
                [~, time_lim_large] = this.getRecTimeSpan(session);
            end

            for r = 1 : this.state.getRecCount()
                this.log.addMarkedMessage(sprintf('Preparing receiver %d of %d', r, this.state.getRecCount()));
                if (numel(rec) < r) || rec(r).isEmpty
                    rec(r) = GNSS_Station(this.state.getConstellationCollector(), this.state.getDynMode(r) == 0); %#ok<AGROW>
                else
                    if this.state.isRinexSession()
                        rec(r).work.resetWorkSpace();
                    end
                end
                if ~this.state.isRinexSession()
                    rec(r).work.setOutLimits(out_limits.first, out_limits.last);
                end
            end
            if numel(rec) > 0
                this.log.newLine();
                this.rec = rec;
                % Init Meteo and Sky objects
                this.initSkySession(time_lim_large);
                this.log.newLine();
                this.initMeteoNetwork(time_lim_large);
                this.log.simpleSeparator();
                
                % inti atmo object
                if this.state.isVMF()
                    atmo = Atmosphere.getInstance();
                    atmo.initVMF(time_lim_large.first,time_lim_large.last);
                end
                if this.state.needIonoMap()
                    atmo = Atmosphere.getInstance();
                    atmo.initIonex(time_lim_large.first,time_lim_large.last);
                end
            end
        end  
        
        function initSkySession(this, time_lim)
            % Init sky for this session
            this.sky = Core_Sky.getInstance();
            this.sky.initSession(time_lim.first, time_lim.last);
        end
        
        function initMeteoNetwork(this, time_lim)
            % Init the meteo network            
            mn = Meteo_Network.getInstance();
            mn.initSession(time_lim.first, time_lim.last);            
        end
        
        function go(this, session_num)
            % Run a session and execut the command list in the settings
            %
            % SYNTAX
            %   this.go(session_num)
            
            t0 = tic;
            if nargin == 1
                session_list = 1 : this.state.getSessionCount();
            else
                session_list = session_num;
            end
            % init refererecne frame object
            rf = Core_Reference_Frame.getInstance();
            rf.init();
            for s = session_list
                this.prepareSession(s);
                this.cmd.exec(this.rec);                   
                if ~this.state.isKeepRecList()
                    for r = 1 : numel(this.rec)
                        this.rec(r).resetOut();
                    end                    
                end
                if ~isunix, fclose('all'); end
            end
            this.log.newLine;
            this.log.addMarkedMessage(sprintf('Computation done in %.2f seconds', toc(t0)));
            this.log.newLine;
        end
        
        function exec(this, cmd)
            this.cmd.exec(this.rec, cmd);
        end
    end
    
    %% CHECK VALIDITY METHODS
    methods
        function err_code = checkValidity(this, level, flag_verbose)
            % Check validity of requiremets
            %
            % INPUT
            %   level        0 - check before conjure, 
            %                1 - single files check (but not downloadable resurces)
            %                5 - check remote resources only
            %                15 - check all necessary files
            %   
            % SYNTAX
            %   core.checkValidity(level, flag_verbosity);
            
            
            if nargin < 3
                flag_verbose = true;
            end         
            
            if nargin < 2 
                level = 1;
            end
            
            this.log.newLine();
            this.log.addMessage('Checking input files and folders...');
            this.log.newLine();
            
            err_code.go = 0; % Global ok check
            
            state = this.state;
            if (level < 5) || (level == 10)
                err_code.home  = state.checkDir('prj_home', 'Home dir', flag_verbose);
                err_code.obs   = state.checkDir('obs_dir', 'Observation dir', flag_verbose);                
            end
            
            if (level == 1) || (level > 10)
                [n_ok, n_ko] = this.checkRinFileList();
                if sum(n_ok) > 0
                    if sum(n_ko) > 0
                        if flag_verbose
                            this.log.addWarning('Some observation files are missing');
                        end
                        err_code.obs_f = sum(n_ko);
                    else
                        if flag_verbose
                            this.log.addStatusOk('Observation rinex are present');
                        end
                        err_code.obs_f = 0;
                    end
                else
                    if flag_verbose
                        this.log.addError('Observation files are missing!!!');
                    end
                    err_code.obs_f = -sum(n_ko);
                end
            end
            
            if (level < 5) || (level > 10)
                err_code.crd   = state.checkDir('crd_dir', 'Coordinate dir', flag_verbose);
            end
            
            if (level < 5) || (level > 10)
                if this.state.isMet()
                    err_code.met   = state.checkDir('met_dir', 'Meteorological dir', flag_verbose);
                    
                    if (level == 1) || (level > 10)
                        %[n_ok, n_ko] = this.checkMetFileList();
                        [n_ok, n_ko] = this.checkFileList(this.state.met_dir, this.state.met_name, [], 0);
                        if sum(n_ok) > 0
                            if sum(n_ko) > 0
                                if flag_verbose
                                    this.log.addWarning('Some met files are missing');
                                end
                                err_code.obs_f = sum(n_ko);
                            else
                                if flag_verbose
                                    this.log.addStatusOk('Met rinex are present');
                                end
                                err_code.obs_f = 0;
                            end
                        else
                            if flag_verbose
                                this.log.addError('Met files are missing!!!');
                            end
                            err_code.obs_f = -sum(n_ko);
                        end
                    end
                else
                    if flag_verbose
                        this.log.addStatusDisabled('Meteorological data not requested');
                    end
                    err_code.met = 0;
                end
            end
            
            if state.isOceanLoading
                err_code.ocean = state.checkDirErr('ocean_dir', 'Ocean loading dir', flag_verbose);
            else
                if flag_verbose
                    this.log.addStatusDisabled('Ocean loading disabled');
                end
                err_code.ocean = 0;
            end
            
            if state.isAtmLoading
                err_code.atm = state.checkDir('atm_load_dir', 'Atmospheric loading dir', flag_verbose);
            else
                if flag_verbose
                    this.log.addStatusDisabled('Athmospheric loading disabled');
                end
                err_code.atm = 0;
            end
            
            err_code.atx   = state.checkDir('atx_dir', 'Antenna dir', flag_verbose);
            err_code.eph   = state.checkDir('eph_dir', 'Ephemerides dir', flag_verbose);
            err_code.clk   = state.checkDir('clk_dir', 'Clock Offset dir', flag_verbose);
            err_code.erp   = state.checkDir('erp_dir', 'Earth Rotation Parameters dir', flag_verbose);
            err_code.crx   = state.checkDir('crx_dir', 'Satellite Manouvers dir', flag_verbose);
            err_code.dcb   = state.checkDir('dcb_dir', 'Differential Code Biases dir', flag_verbose);
            err_code.ems   = state.checkDir('ems_dir', 'EGNOS Message Center dir', flag_verbose);
            
            %err_code.geoid = state.checkDir('geoid_dir', 'Geoid loading dir', flag_verbose);
            err_code.geoid_f = state.checkFile({'geoid_dir', 'geoid_name'}, 'Geoid file', flag_verbose);
            
            geoid = this.gc.getRefGeoid();
            if isempty(geoid) || isempty(geoid.grid)
                err_code.geoid = -1;
                if flag_verbose
                    this.log.addWarning('The Geoid is missing');
                end
            else
                err_code.geoid = 0;
                if flag_verbose
                    this.log.addStatusOk('The Geoid is loaded');
                end
            end
            
            if state.isHOI
                err_code.iono = state.checkDir('iono_dir', 'Ionospheric map dir', flag_verbose);
                err_code.igrf = state.checkDir('igrf_dir', 'International Geomagnetic Reference Frame dir', flag_verbose);
            else
                if flag_verbose
                    this.log.addStatusDisabled('High order ionospheric corrections disabled');
                end
                err_code.iono = 1;
                err_code.igrf = 1;
            end
            
            if state.isVMF
                err_code.vmf = state.checkDir('vmf_dir', 'Vienna Mapping Function dir', flag_verbose);
            else
                if flag_verbose
                    this.log.addStatusDisabled('Not using Vienna Mapping functions');
                end
                err_code.vmf = 1;
            end
            
            this.log.newLine();
            
            % Checking folder that are not created in conjure phase
            err_code.go = err_code.home + ...
                err_code.obs + ...
                (err_code.obs_f < 0) + ...
                (err_code.ocean < 0) + ...
                (err_code.atm < 0) + ...
                err_code.atx;
        end
    end
    
    %% RIN FILE LIST
    % ==================================================================================================================================================
    methods
        function [n_ok, n_ko] = updateRinFileList(this, force_update, verbosity)
            % Update and check rinex list
            %
            % SYNTAX
            %   this.updateRinFileList()
            if nargin < 2 || isempty(force_update)
                force_update = true;
            end
            if nargin < 3 || isempty(verbosity)
                verbosity = false;
            end
            if verbosity
                this.log.addMarkedMessage('Checking RINEX input files');
            end
            this.state.updateObsFileName;
            n_rec = this.state.getRecCount;
            rec_path = this.state.getRecPath;
            
            clear fr
            if ~force_update
                fr = this.getRinFileList();
            end
            n_ok = zeros(n_rec, 1);
            n_ko = zeros(n_rec, 1);
            for r = 1 : n_rec
                name = File_Name_Processor.getFileName(rec_path{r}{1});
                if verbosity
                    this.log.addMessage(this.log.indent(sprintf('- %s ...checking...', upper(name(1:4)))));
                end
                if force_update
                    fr(r) = File_Rinex(rec_path{r}, 100);
                end
                n_ok(r) = sum(fr(r).is_valid_list);
                n_ko(r) = sum(~fr(r).is_valid_list);
                if verbosity
                    this.log.addMessage(sprintf('%s (%d ok - %d ko)', char(8 * ones(1, 2 + 14)), n_ok(r), n_ko(r)));
                end
            end            
            
            this.rin_list = fr;
            if n_rec == 0
                % 'No receivers found';
            end
        end
        
        function [n_ok, n_ko] = updateMetFileList(this, force_update, verbosity)
            % Update and check met RINEX list
            %
            % SYNTAX
            %   this.updateMetFileList()
            if nargin < 2 || isempty(force_update)
                force_update = true;
            end
            if nargin < 3 || isempty(verbosity)
                verbosity = false;
            end
            if verbosity
                this.log.addMarkedMessage('Checking Meteorological RINEX input files');
            end
            
            station_file_list = this.getRequiredMetFile();
            if isempty(station_file_list)
                n_ok = 0;
                n_ko = 0;
                this.met_list = [];
            else
                n_rec = size(station_file_list, 2);
                n_ok = zeros(n_rec, 1);
                n_ko = zeros(n_rec, 1);
                for r = 1 : n_rec
                    name = File_Name_Processor.getFileName(station_file_list{r}{1});
                    if verbosity
                        this.log.addMessage(this.log.indent(sprintf('- %s ...checking...', upper(name(1:4)))));
                    end
                    if force_update
                        fr(r) = File_Rinex(station_file_list{r}, 100);
                    end
                    n_ok(r) = sum(fr(r).is_valid_list);
                    n_ko(r) = sum(~fr(r).is_valid_list);
                    if verbosity
                        this.log.addMessage(sprintf('%s (%d ok - %d ko)', char(8 * ones(1, 2 + 14)), n_ok(r), n_ko(r)));
                    end
                end
                
                this.met_list = fr;
                if n_rec == 0
                    % 'No receivers found';
                end
            end
        end
        
        function [n_ok, n_ko] = checkFileList(this, dir_path, file_name, time_par, margin)
            % Check file list
            %
            % SYNTAX
            %   [n_ok, n_ko] = this.checkRinFileList(dir_path, file_name)
            %   [n_ok, n_ko] = this.checkRinFileList(dir_path, file_name, margin)
            %   [n_ok, n_ko] = this.checkRinFileList(dir_path, file_name, time_limits)                        
            
            verbosity = true;
            
            if nargin < 4 || isempty(time_par) || ~isclass(time_par, 'GPS_Time') 
                [~, time_lim_large, is_empty] = this.getRecTimeSpan();
                time_lim = time_lim_large;
            else
                time_lim = time_par.getCopy();
            end
            
            if ~isa(time_par, 'GPS_Time')
                % is a margin!!!
                time_lim.first.addIntSeconds(margin(1));
                time_lim.first.addIntSeconds(margin(end));
            end
                        
            file_list = {};
            fnp = File_Name_Processor();
            if ~iscell(file_name)
                file_name = {file_name};
            end
            for i = 1 : numel(file_name)
                file_list{i} = fnp.dateKeyRepBatch(fnp.checkPath(strcat(dir_path, filesep, file_name{i})), time_lim.first,  time_lim.last, this.state.sss_id_list, this.state.sss_id_start, this.state.sss_id_stop);
            end
            
            if isempty(file_list)
                n_ok = 0;
                n_ko = 0;
                this.met_list = [];
            else
                n_rec = size(file_list, 2);
                n_ok = zeros(n_rec, 1);
                n_ko = zeros(n_rec, 1);
                for r = 1 : n_rec
                    name = File_Name_Processor.getFileName(file_name{r});
                    name = name(1:4);
                    if verbosity
                        this.log.addMessage(this.log.indent(sprintf('- %s ...checking...', name)));
                    end
                    for s = 1 : numel(file_list{r})
                        if (exist(file_list{r}{s}, 'file') == 2)
                            n_ok(r) = n_ok(r) + 1;
                        else
                            n_ko(r) = n_ko(r) + 1;
                        end
                    end                    
                    if verbosity
                        this.log.addMessage(sprintf('%s (%d ok - %d ko)', char(8 * ones(1, 2 + 14)), n_ok(r), n_ko(r)));
                    end
                end
                
                if n_rec == 0
                    % 'No receivers found';
                end
            end
        end

        
        function [n_ok, n_ko] = checkRinFileList(this, force_update)
            % Update and check rinex list
            %
            % SYNTAX
            %   [n_ok, n_ko] = this.checkRinFileList()
            if nargin < 2 || isempty(force_update)
                force_update = true;
            end            
            [n_ok, n_ko] = this.updateRinFileList(force_update, true);
        end
        
        function [n_ok, n_ko] = checkMetFileList(this, force_update)
            % Update and check meteorological rinex list
            %
            % SYNTAX
            %   [n_ok, n_ko] = this.checkMetFileList()
            if nargin < 2 || isempty(force_update)
                force_update = true;
            end
            [n_ok, n_ko] = this.updateMetFileList(force_update, true);
        end
        
        function rin_list = getRinFileList(this)
            % Update and check rinex list
            %
            % SYNTAX
            %   rin_list = this.getRinFileList()            
            if isempty(this.rin_list)
                this.updateRinFileList()
            end            
            rin_list = this.rin_list;
        end
        
        function [time_lim_small, time_lim_large, is_empty] = getRecTimeSpan(this, session)
            % return a GPS_Time containing the first and last epoch for a session
            %
            % OUTPUT:
            %   time_lim_small     GPS_Time (first and last) epoch of the smaller interval
            %   time_lim_large     GPS_Time (first and last) epoch of the larger interval
            %
            % SYNTAX:
            %   [time_lim_small, time_lim_large] = this.getRecTimeSpan(session)
            %   [time_lim_small, time_lim_large] = this.getRecTimeSpan()
            
            fr = this.getRinFileList();
            is_empty = ~fr.isValid;
            if nargin == 1 % Start and stop limits of all the sessions                
                time_lim_small = fr(1).first_epoch.first;
                tmp_small = fr(1).last_epoch.last;
                time_lim_large = time_lim_small.getCopy;
                tmp_large = tmp_small.getCopy;
                for r = 2 : numel(fr)
                    if fr(r).isValid()
                        if isempty(time_lim_small) || time_lim_small < fr(r).getFirstEpoch.first
                            time_lim_small = fr(r).getFirstEpoch.first;
                        end
                        if isempty(time_lim_large) || time_lim_large > fr(r).getFirstEpoch.first
                            time_lim_large = fr(r).getFirstEpoch.first;
                        end
                        
                        if isempty(tmp_small) || tmp_small > fr(r).getLastEpoch.last
                            tmp_small = fr(r).getLastEpoch.last;
                        end
                        if isempty(tmp_large) || tmp_large < fr(r).getLastEpoch.last
                            tmp_large = fr(r).getLastEpoch.last;
                        end
                    end
                end
                time_lim_small.append(tmp_small);
                time_lim_large.append(tmp_large);
            else % Start and stop of a certain session
                time_lim_small = fr(1).getFirstEpoch(session);
                tmp_small = fr(1).getLastEpoch(session);
                time_lim_large = time_lim_small.getCopy;
                tmp_large = tmp_small.getCopy;
                for r = 2 : numel(fr)
                    if fr(r).isValid(session)
                        if isempty(time_lim_small) || time_lim_small < fr(r).getFirstEpoch(session)
                            time_lim_small = fr(r).getFirstEpoch(session);
                        end
                        if isempty(time_lim_large) || time_lim_large > fr(r).getFirstEpoch(session)
                            time_lim_large = fr(r).getFirstEpoch(session);
                        end
                        
                        if isempty(tmp_small) || tmp_small > fr(r).getLastEpoch(session)
                            tmp_small = fr(r).getLastEpoch(session);
                        end
                        if isempty(tmp_large) || tmp_large < fr(r).getLastEpoch(session)
                            tmp_large = fr(r).getLastEpoch(session);
                        end
                    end
                end
                time_lim_small.append(tmp_small);
                time_lim_large.append(tmp_large);
            end
        end
        
        function cur_session = getCurSession(this)
            % Get the id of the current session
            %
            % SYNTAX
            %   cur_session = this.getCurSession()            
            cur_session = this.cur_session;
        end
    end
    
    methods
        function file_list = getRequiredMetFile(this)
            % Return the list of met_file local paths that will be loaded
            % 
            % SYNTAX 
            %    file_list = getRequiredMetFile(this)
            [~, time_lim_large, is_empty] = this.getRecTimeSpan();
            if ~is_empty
                file_list = this.state.getMetFileName(time_lim_large.first, time_lim_large.last);
            else
                file_list = {};
            end
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
