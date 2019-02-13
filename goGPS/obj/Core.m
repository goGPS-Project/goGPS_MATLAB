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
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 2
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
        GO_GPS_VERSION = '1.0 beta 2';
        GUI_MODE = 0; % 0 means text, 1 means GUI, 5 both
    end

    %% PROPERTIES MISC
    % ==================================================================================================================================================
    properties (GetAccess = private, SetAccess = private) % Public Access
        local_storage = '';
        
        creation_time = GPS_Time(now);
        is_advanced = true;
        
        session_list = [];  % in the go function this variable is update to keep the list of session to compute
        preloaded_session = 0;  % Current preloaded session in Core (for meteo, orbits, and other data)
    end
    
    %% PROPERTIES SINGLETON POINTERS
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons
        log         % Logger handler
        w_bar       % Wait_Bar handler
        
        state       % state
        sky         % Core_Sky handler
        atmo        % Atmosphere handler
        mn          % Meteorological Network handler
        rf          % Reference Frame handler
        cmd         % Command_Interpreter handler
        
        geoid = struct('file', [], 'grid', 0, 'cellsize', 0, 'Xll', 0, 'Yll', 0, 'ncols', 0, 'nrows', 0); % parameters of the reference geoid
        
        gom         % Parallel controller
    end

    %% PROPERTIES RECEIVERS
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons        
        rin_list        % List of observation file (as File_Rinex objects) to store minimal information on the input files
        met_list        % List of meteorological file (as File_Rinex objects) to store minimal information on the input files
        
        rec             % List of all the receiver used
        
        net             % List of all the network used
    end

    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = private)
        % Concrete implementation.  See Singleton superclass.
        function this = Core()
            % Core object creator
            
            % init logger
            this.log = Logger.getInstance();
            
            if ispc()
                home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
                this.local_storage = [home '\AppData\Local\goGPS'];
            else
                home = getenv('HOME');
                if ismac()
                    this.local_storage = [home '/Library/Application Support/goGPS'];
                else
                    this.local_storage = [home '/.goGPS'];
                end
            end
            if ~(exist(this.local_storage, 'dir'))
                mkdir(this.local_storage)
            end
        end
    end    
    
    %% METHODS INIT & STATIC GETTERS & SETTERS
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = getInstance(force_clean, skip_init)                        
            % Get the persistent instance of the class
            %
            % INPUT
            %   force_clean  remove rec and reset state from core
            %   skip_init    take the object core as it is
            %
            % SYNTAX
            %    this = getInstance(force_clean, skip_init)
            %
            % ALTERNATIVE
            % Set the persistent instance of the class
            % substitute the current persistent link to the core
            %
            % INPUT 
            %   core        existent core
            %
            % SYNTAX
            %   core = getInstance(core)
            
            
            persistent unique_instance_core__

            if (nargin == 1) && isa(force_clean, 'Core')
                unique_instance_core__ = force_clean;
                this = force_clean;
            else                
                if nargin < 1 || isempty(force_clean)
                    force_clean = false;
                end
                if nargin < 2 || isempty(skip_init)
                    skip_init = false;
                end
                
                
                if isempty(unique_instance_core__)
                    this = Core();
                    unique_instance_core__ = this;
                    if ~skip_init
                        this.init(force_clean);
                    end
                else
                    this = unique_instance_core__;
                    if ~skip_init
                        this.init(force_clean);
                    end
                end
            end
        end
        
        
        function this = setInstance(core)
            % Set the persistent instance of the class
            % substitute the current persistent link to the core
            %
            % INPUT
            %   core        existent core
            %
            % SYNTAX
            %   core = Core.setInstance(core)
            
            this = Core.getInstance(core);
        end
        
        function initGeoid(geoid)
            % load external geoid
            %
            % SYNTAX
            %   this.initGeoid(); -> load from file
            %   this.initGeoid(geoid); -> import from obj
            core = Core.getCurrentCore;
            if nargin == 1
                core.geoid = geoid;
            else
                try
                    geoid_file = core.state.getGeoidFile();
                    if ~exist(core.state.getGeoidDir, 'file')
                        core.state.geoid_dir = File_Name_Processor.getFullDirPath('../data/reference/geoid', pwd);
                        geoid_file = core.state.getGeoidFile();
                    end
                    if ~strcmp(core.geoid.file, geoid_file) || (core.geoid.ncols == 0 || core.geoid.nrows == 0) 
                        g = load(geoid_file);
                        fn = fieldnames(g);
                        % geoid grid and parameters
                        core.geoid.file = geoid_file;
                        core.geoid.grid = g.(fn{1});
                        core.geoid.cellsize = 360 / size(core.geoid.grid, 2);
                        core.geoid.Xll = -180 + core.geoid.cellsize / 2;
                        core.geoid.Yll = -90 + core.geoid.cellsize / 2;
                        core.geoid.ncols = size(core.geoid.grid, 2);
                        core.geoid.nrows = size(core.geoid.grid, 1);
                        clear g
                    end
                catch
                    core.log.addWarning('Reference geoid not found');
                    % geoid unavailable
                    core.geoid.grid = 0;
                    core.geoid.cellsize = 0;
                    core.geoid.Xll = 0;
                    core.geoid.Yll = 0;
                    core.geoid.ncols = 0;
                    core.geoid.nrows = 0;
                end
            end
        end
        
        function [go_dir] = getLocalStorageDir()
            % Get local storage
            core = Core.getInstance(false, true);
            
            go_dir = core.local_storage;
        end       

        function geoid = getRefGeoid()
            % Get reference geoid
            core = Core.getInstance(false, true);
            
            if isempty(core.geoid)
                core.initGeoid();
            end
            geoid = core.geoid;
        end
        
        function [is_adv] = isAdvanced()
            % Get the status of usage (normal/advanced)
            %
            % SYNTAX:
            %   this.isAdvanced()
            core = Core.getInstance(false, true);
            
            is_adv = core.is_advanced;
        end    
    
        function ok_go = openGUI()
            ok_go = gui_goGPS;
        end
              
        function log = getLogger()
            % Return the pointer to the Logger Object
            %
            % SYNTAX
            %   log = Core.getLogger()
            
            core = Core.getInstance(false, true);
            log = core.log;
            if isempty(log)
                log = Logger.getInstance();
                core.log = log;
            end            
        end
        
        function wb = getWaitBar()
            % Return the pointer to the Go_Wait_Bar
            %
            % SYNTAX
            %   log = Core.getWaitBar()
            
            core = Core.getInstance(false, true);
            wb = core.w_bar;
            if isempty(wb)
                wb = Go_Wait_Bar.getInstance();
                core.wb = wb;
            end            
        end
        
        function atmo = getAtmosphere()
            % Return the pointer to the Atmosphere Object
            %
            % SYNTAX
            %   atmo = Core.getAtmosphere()
            
            core = Core.getInstance(false, true);
            atmo = core.atmo;
            if isempty(atmo)
                atmo = Atmosphere();
                core.mn = atmo;
            end            
        end
        
        function mn = getMeteoNetwork()
            % Return the pointer to the Meteorological Network Object
            %
            % SYNTAX
            %   mn = Core.getMeteoNetwork()
            
            core = Core.getInstance(false, true);
            mn = core.mn;
            if isempty(mn)
                mn = Meteo_Network();
                core.mn = mn;
            end
        end
        
        function rf = getReferenceFrame(flag_reset)
            % Return the pointer to the Reference Frame
            %
            % SYNTAX
            %   rf = Core.getReferenceFrame()
            
            core = Core.getInstance(false, true);
            rf = core.rf;
            if isempty(rf)
                rf = Core_Reference_Frame();
                core.rf = rf;
            end
            if nargin == 1 && flag_reset
            %    rf.reset();
            end
        end
        
        function sky = getCoreSky()
            % Return the pointer to the Core Sky Object
            %
            % SYNTAX
            %   sky = Core.getCoreSky()
            
            core = Core.getInstance(false, true);
            sky = core.sky;
            if isempty(sky)
                sky = Core_Sky();
                core.sky = sky;
            end
            
        end
        
        function rec_list = getRecList()
            % Get receiver list
            %
            % SYNTAX
            %   rec_list = Core.getRecList()          
            core = Core.getInstance(false, true);
            rec_list = core.rec;
        end
                     
        function [cur_session, session_list, cur_pos] = getCurrentSession()
            % Get the id of the current session
            %
            % OUTPUT
            %   cur_session     current session
            %   session_list    list of active sessions to process (in go)
            %   cur_pos         pos of cur_session in session_list
            %
            % SYNTAX
            %   [cur_session, session_list, cur_pos] = Core.getCurrentSession()
            core = Core.getInstance(false, true);
            cur_session = core.state.getCurSession;
            if nargout >= 2
                session_list = core.session_list;
                if nargout >= 3
                    cur_pos = find(session_list ==cur_session);
                end
            end
        end

        function [rin_list, met_list] = getRinLists()
            % Get the rinex lists (GNSS observations and meteorological)
            %
            % SYNTAX
            %   [rin_list, met_list] = Core.getRinLists()
            core = Core.getInstance(false, true);
            rin_list = core.rin_list;
            met_list = core.met_list;
        end              
        
        function state = getState()
            % Return the pointer to the State pointed by Core
            %
            % SYNTAX
            %   state = Core.getState()
            
            state = Core.getCurrentSettings;
        end
        
        function state = getCurrentSettings(ini_settings_file)
            % Get the persistent settings
            %
            % SYNTAX
            %   state = Core.getCurrentSettings(<ini_settings_file>)
            
            core = Core.getInstance(false, true);
            if isempty(core.state)
                core.state = Main_Settings();
            end
            if nargin == 1 && ~isempty(ini_settings_file)
                this.state.importIniFile(ini_settings_file);
            end
            % Return the handler to the object containing the current settings
            state = handle(core.state);
        end
        
        function cmd = getCommandInterpreter()
            % Return the pointer to the Command Interpreter
            %
            % SYNTAX
            %   cmd = Core.getCommandInterpreter()
            
            core = Core.getInstance(false, true);
            if isempty(core.cmd)
                core.cmd = Command_Interpreter(core);
            end
            cmd = core.cmd;
        end

        function core = getCurrentCore()
            % Return the pointer to the actual Core instance
            %
            % SYNTAX
            %   core = Core.getCurrentCore()
            
            core = Core.getInstance(false, true);
        end
        
        function core = setCurrentCore(core)
            % Set the pointer to the actual Core instance
            %
            % SYNTAX
            %   core = Core.setCurrentCore(core)
            
            core = Core.getInstance(core);
        end

        function setAdvanced(mode)
            % Set the status of usage (normal/advanced)
            %
            % SYNTAX:
            %   this.setAdvanced(<mode>)
            core = Core.getInstance(false, true);
            if nargin == 0
                mode = true;
            end
            core.is_advanced = mode;
        end        

        function setAtmosphere(atmo)
            % Set the pointer to the Atmosphere Object
            %
            % SYNTAX
            %   Core.setAtmosphere(atmo)
            
            core = Core.getInstance(false, true);
            core.atmo = atmo;
        end

        function setMeteoNetwork(mn)
            % Set the pointer to the Meteorological Network Object
            %
            % SYNTAX
            %   Core.setMeteoNetwork(mn)
            
            core = Core.getInstance(false, true);
            core.mn = mn;
        end

        function setReferenceFrame(rf)
            % Return the pointer to the Reference Frame
            %
            % SYNTAX
            %   Core.setReferenceFrame(rf)
            
            core = Core.getInstance(false, true);
            core.rf = rf;
        end

        function setCoreSky(sky)
            % Set the pointer to the Core Sky Object
            %
            % SYNTAX
            %   Core.setCoreSky(sky)
            
            core = Core.getInstance(false, true);
            core.sky = sky;
        end
        
        function setState(state)
            % Set the pointer to the State pointed by Core
            %
            % SYNTAX
            %   Core.setState(state)
            
            core = Core.getInstance(false, true);
            core.state = state;
        end
        
        function setCurrentSettings(cs)
            % Set the persistent settings
            %
            % SYNTAX
            %   this.setCurrentSettings(ini_settings_file)
            core = Core.getInstance(false, true);
            core.state = cs;
        end
        
        function setCommandInterpreter(cmd)
            % Return the pointer to the Command Interpreter
            %
            % SYNTAX
            %   Core.setCommandInterpreter(cmd)
            
            core = Core.getInstance(false, true);
            core.cmd = cmd;
        end
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this, force_clean)
            % Get instances for:
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
            this.log.setOutMode([], false);                
            if ispc, fclose('all'); end

            cm = this.log.getColorMode();
            this.log.setColorMode(true);   
            
            Core_UI.showTextHeader();
            this.log.setColorMode(cm);            
            this.state = Main_Settings();
            this.w_bar = Go_Wait_Bar.getInstance(100,'Welcome to goGPS', Core.GUI_MODE);  % 0 means text, 1 means GUI, 5 both
            this.sky = Core_Sky(force_clean);
            this.atmo = Atmosphere();
            this.rf = Core_Reference_Frame();
            this.cmd = Command_Interpreter(this);
            if force_clean
                this.preloaded_session = 0;
                this.state.setCurSession(0);
                this.rec = [];
            end
        end

        function initSimpleHandlers(this)
            % Get instances for:
            %   - Logger
            %   - Wait_Bar
            %   - Command_Interpreter
            %
            % SYNTAX:
            %   this.initSimpleHandlers()
            
            this.log = Logger.getInstance;
            this.log.setOutMode([], false);
            if ispc, fclose('all'); end
            this.w_bar = Go_Wait_Bar.getInstance(100, 'Welcome to goGPS', Core.GUI_MODE);  % 0 means text, 1 means GUI, 5 both
            this.cmd = Command_Interpreter(this);
        end
        
        function initConfiguration(this)
            % Load all the files necessary to the functioning of a goGPS session
            % SYNTAX:   this.initProcessing()
            
            % Load external resources and update
            fnp = File_Name_Processor();
            out_dir = fnp.checkPath(this.state.getOutDir());
            if ~exist(out_dir, 'dir')
                this.log.newLine();
                this.log.addWarning(sprintf('The out folder does not exists\n Creating %s', out_dir));
                mkdir(out_dir);
            end

            this.log.addMessage(this.log.indent(this.state.cc.toString, 5));

            this.initGeoid();
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
                this.state.importIniFile(ini_settings_file);
            end
        end
        
        function importState(this, state)
            % Import Settings from state
            %
            % SYNTAX:
            %   this.importState(state)
            
            this.state.import(state);
        end
        
        function prepareProcessing(this, flag_download)
            % Init settings, and download necessary files
            %
            % SYNTAX:
            %   this.prepareProcessing(flag_download)
            if nargin == 2
                this.state.setAutomaticDownload(flag_download);
            end
            
            this.log.newLine();
            this.log.addMarkedMessage(sprintf('PROJECT: %s', this.state.getPrjName()));

            this.initConfiguration(); % Set up / download observations and navigational files
            this.log.addMarkedMessage('Conjuring all the auxilliary files!');
            rin_list = this.getRinFileList();
            
            fw = File_Wizard;
            c_mode = this.log.getColorMode();
            this.log.setColorMode(false);
            if ~this.state.isRinexSession()

                    [buff_lim, ~] = this.state.getSessionLimits(1);
                    time_lim_large = buff_lim.getEpoch(1);
                    [buff_lim, ~] = this.state.getSessionLimits(this.state.getSessionCount());
                    time_lim_large.append(buff_lim.getEpoch(2));

            else
              [~, time_lim_large, is_empty] = this.getRecTimeSpan();
            end
            fw.conjureFiles(time_lim_large.first, time_lim_large.last);
            this.log.setColorMode(c_mode);
        end        
        
        function activateParallelWorkers(this, id_rec2pass)
            % Call the activation of the receivers
            %
            % INPUT
            %   id_rec2pass     id of the work receiver to broadcast to the receivers
            %
            % SYNTAX
            %   this.activateParallelWorkers(id_rec2pass)
            this.gom = Parallel_Manager.getInstance;
            this.gom.activateWorkers(id_rec2pass);
        end
    end
    
    %% METHODS EXPORT
    methods
        function logCurrentSettings(this)
            log = Core.getLogger();
            log.addMessageToFile('\n============================================================================================');
            log.addMessageToFile('== CONFIG FILE =============================================================================');
            log.addMessageToFile('============================================================================================\n');
            log.addMessageToFile(this.state.export);
            log.addMessageToFile('\n============================================================================================');
            log.addMessageToFile('== END OF CONFIG FILE ======================================================================');
            log.addMessageToFile('============================================================================================\n');
        end
    end
    
    %% METHODS RUN
    % ==================================================================================================================================================
    methods
        function is_empty = prepareSession(this, session_number, flag_preload)
            % Check the time-limits for the files in the session
            % Init the Sky and Meteo object
            %
            % SYNTAX
            %   is_empty = this.prepareSession(session_number)
            
            if nargin < 3 || isempty(flag_preload)
                flag_preload = false;
            end
            
            session = session_number;
            if ~flag_preload
                this.state.setCurSession(session_number);
                
                this.log.newLine;
                this.log.simpleSeparator();
                this.log.addMessage(sprintf('Starting session %d of %d', session, this.state.getSessionCount()));
                
                if isempty(this.rec)
                    rec = GNSS_Station;
                else
                    rec = this.rec;
                end
                
                this.log.newLine();                
            elseif (session ~= this.preloaded_session)
                this.log.addMarkedMessage(sprintf('Preloading session %d', session), 0);
            end
            
            rin_list = this.getRinFileList();
            if ~this.state.isRinexSession()
                if ~isempty(rin_list)
                    [buff_lim, out_limits] = this.state.getSessionLimits(session);
                    time_lim_large = buff_lim;
                end
            else
                [out_limits, time_lim_large] = this.getRecTimeSpan(session);
            end
                        
            if out_limits.length < 2 || (~this.state.isRinexSession() && ~rin_list.isValid) || ((time_lim_large.last - time_lim_large.first) < 0) % && ~rin_list_chk.isValid
                is_empty = true;
                this.log.addMessage(sprintf('No valid receivers are present / session not valid %d', session));
            else

                is_empty = false;
                this.log.addMessage(sprintf('Begin %s', out_limits.first.toString()));
                this.log.addMessage(sprintf('End   %s', out_limits.last.toString()));
                this.log.simpleSeparator();
                
                if ~flag_preload
                    for r = 1 : this.state.getRecCount()
                        this.log.addMessage(sprintf('[ -- ] Preparing receiver %d of %d', r, this.state.getRecCount()));
                        if (numel(rec) < r) || rec(r).isEmpty
                            rec(r) = GNSS_Station(this.state.getConstellationCollector(), this.state.getDynMode() == 0);
                        else
                            buf = min(round(rec(r).work.state.getBuffer / rec(r).work.getRate), rec(r).work.length);
                            if buf == 0
                                rec(r).work.resetWorkSpace();
                                rec(r).old_work = Receiver_Work_Space(rec(r).work.cc, rec(r).work.parent);
                            else
                                % keep the old work
                                rec(r).old_work = rec(r).work;
                                rec(r).old_work.keepEpochs(rec(r).work.length + ((-buf + 1) : 0));
                                % reset the new work
                                rec(r).work = Receiver_Work_Space(rec(r).work.cc, rec(r).work.parent);
                            end
                        end
                    end
                end
                if flag_preload || numel(rec) > 0
                    this.log.newLine();
                    if ~flag_preload
                        this.rec = rec;
                    end
                                        
                    % Init Meteo and Sky objects
                    if ~this.state.isNoResources() && (session ~= this.preloaded_session)
                        this.initSkySession(time_lim_large);
                        this.log.newLine();
                        if this.state.isMet()
                            this.initMeteoNetwork(time_lim_large);
                        end
                        this.log.simpleSeparator();
                        
                        % init atmo object
                        if this.state.isVMF()
                            this.atmo.initVMF(time_lim_large.first,time_lim_large.last);
                        end
                        if this.state.needIonoMap() && ~this.state.isIonoBroadcast()
                            this.atmo.initIonex(time_lim_large.first,time_lim_large.last);
                        end
                        this.preloaded_session = session;
                    end
                end
            end
        end  
        
        function initSkySession(this, time_lim)
            % Init sky for this session
            if isempty(this.sky)
                this.sky = Core_Sky();
            end
            if ~this.state.isNoResources()
                this.sky.initSession(time_lim.first, time_lim.last);
            end
        end
        
        function initMeteoNetwork(this, time_lim)
            % Init the meteo network            
            this.mn = Meteo_Network();
            this.mn.initSession(time_lim.first, time_lim.last);            
        end
        
        function go(this, session_num)
            % Run a session and execute the command list in the settings
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
            this.rf.init();
            
            [cmd_list, ~, execution_block, sss_list, ~, sss_level, flag_push] = this.cmd.fastCheck(this.state.cmd_list);
            
            cmd_line = 1;
            last_sss = 0;
            for eb = unique(execution_block) % for each execution block
                t1 = tic;
                if cmd_line <= numel(sss_list)
                    sessions = sss_list{cmd_line};
                    this.session_list = sessions;
                    this.net = [];
                    for s = sessions
                        if s ~= last_sss
                            is_empty = this.prepareSession(s);
                        end
                        last_sss = s;
                        if ~is_empty
                            this.exec(cmd_list(execution_block == eb), sss_level(execution_block == eb));
                            cmd_line = find(execution_block == eb, 1, 'last') + 1;
                            
                            % to eventually reset the Out
                            % for r = 1 : numel(this.rec)
                            %     this.rec(r).resetOut();
                            % end
                            if flag_push(eb+1)
                                for i = 1 : length(this.rec)
                                    this.rec(i).work.pushResult();
                                end
                            end
                        end
                    end
                    
                    this.log.newLine;
                    this.log.simpleSeparator([], 'Green');
                    if ~isempty(sessions)
                        this.log.addMessage(sprintf('End of session loop from %d to %d in %.3f seconds', sessions(1), sessions(end), toc(t1)));
                    end
                end
            end            
            this.log.newLine;
            this.log.addMarkedMessage(sprintf('Computation done in %.2f seconds', toc(t0)));
            this.log.newLine;          
        end
        
        function exec(this, cmd_list, level)
            % Execute a list of commands on the current session
            %
            % SYNTAX
            %   this.exec(cmd_list, level)
            % 
            % EXAMPLE
            %   core.exec({'LOAD T*', 'PREPRO T*', 'PPP T*'})
            if nargin < 3 || isempty(level)
                level = zeros(size(cmd_list,1), 1);
            end
            this.cmd.exec(this.rec, cmd_list, level);
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
            
            err_code.atx   = state.checkDirErr('atx_dir', 'Antenna dir', flag_verbose);
            if state.isHOI
                err_code.hoi   = state.checkDirErr('igrf_dir', 'International Geomagnetic reference dir', flag_verbose);
            else
                err_code.hoi   = state.checkDir('igrf_dir', 'International Geomagnetic reference dir', flag_verbose);
            end
            err_code.eph   = state.checkDir('eph_dir', 'Ephemerides dir', flag_verbose);
            err_code.clk   = state.checkDir('clk_dir', 'Clock Offset dir', flag_verbose);
            err_code.erp   = state.checkDir('erp_dir', 'Earth Rotation Parameters dir', flag_verbose);
            err_code.crx   = state.checkDir('crx_dir', 'Satellite Manouvers dir', flag_verbose);
            err_code.dcb   = state.checkDir('dcb_dir', 'Differential Code Biases dir', flag_verbose);
            err_code.ems   = state.checkDir('ems_dir', 'EGNOS Message Center dir', flag_verbose);
            
            %err_code.geoid = state.checkDir('geoid_dir', 'Geoid loading dir', flag_verbose);
            err_code.geoid_f = state.checkFile({'geoid_dir', 'geoid_name'}, 'Geoid file', flag_verbose);
            
            geoid = Core.getRefGeoid();
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
                err_code.atx + ...
                err_code.hoi * state.isHOI;
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
            
            if isempty(file_name)
                n_ok = 0;
                n_ko = 0;
            else
                if nargin < 4 || isempty(time_par) || ~isclass(time_par, 'GPS_Time')
                    %[~, time_lim_large, is_empty] = this.getRecTimeSpan();
                    time_lim_large = this.state.getSessionsStartExt;
                    time_lim_large.append(this.state.getSessionsStopExt);
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
                    file_list{i} = fnp.dateKeyRepBatch(fnp.checkPath(strcat(dir_path, filesep, file_name{i})), this.state.getSessionsStartExt,  this.state.getSessionsStopExt, this.state.sss_id_list, this.state.sss_id_start, this.state.sss_id_stop);
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
            end
            
            % Check start limits out of sessions
            if time_lim_small < this.state.getSessionsStartExt()
                time_lim_small = this.state.getSessionsStartExt();
            end
            if time_lim_large < this.state.getSessionsStartExt()
                time_lim_large = this.state.getSessionsStartExt();
            end
            
            % Check stop limits out of sessions
            if tmp_small > this.state.getSessionsStopExt()
                tmp_small = this.state.getSessionsStopExt();
            end
            if tmp_large > this.state.getSessionsStopExt()
                tmp_large = this.state.getSessionsStopExt();
            end
            time_lim_small.append(tmp_small);
            time_lim_large.append(tmp_large);
        end
        
        function net = getNetwork(this, rid, rec_list)
            % get a network based on the receiver id, 'if is not present crate a new one
            %
            % SYNTAX:
            %    net = this.getNetwork(rid)
            if nargin < 2
                rec_list = [];
            end
            for i = 1: length(this.net)
                if Core_Utils.permutedEqual(this.net(i).net_id,rid)
                    % Update rec_list
                    this.net(i).rec_list = rec_list(this.net(i).net_id); 
                    net = this.net(i);
                    return
                end
            end
            % if not found
            if isempty(this.net)
                this.net = Network(rec_list(rid), rid);
            else
                this.net(end+1) = Network(rec_list(rid), rid);
            end
            net = this.net(end);
        end
        
        function cur_session = getCurSession(this)
            % Get the id of the current session
            %
            % SYNTAX
            %   cur_session = this.getCurSession()             
            cur_session = this.state.getCurSession;
        end
        
        function id = getStationId(this, marker_name)
            % Given a marker_name get the sequencial id of a station
            %
            % SYNTAX
            %   id = findStationId(this, marker_name)
            marker4ch_list = '';
            for r = 1 : numel(this.rin_list)
                try
                    marker4ch_list(r, :) = char(this.rin_list(r).file_name_list{1}(1 : 4));
                catch
                    % the name is shorter or missing => ignore
                end
            end
            id = find(Core_Utils.code4Char2Num(upper(marker4ch_list)) == Core_Utils.code4Char2Num(upper(marker_name)));
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
        
        function toString(this)
            % Display on screen information about the core object
            %
            % INPUT
            %   this    core object
            %
            % SYNTAX
            %   this.toString();
            fprintf('----------------------------------------------------------------------------------\n')
            this.log.addMarkedMessage(sprintf('Core object created at %s\n', this.creation_time.toString));
            fprintf('----------------------------------------------------------------------------------\n')
            this.log.newLine();
        end
    end

    methods (Static) % Public Access
        
        function clearSingletons()
            % clear all singletons
            %
            % SYNTAX:
            % Core.clearSingletons
            clear Core ...
                Core_UI ...
                GUI_Main ...
                Remote_Resource_Management ...
                Parallel_Manager ...
                Go_Slave ...
                Go_Wait_Bar ...
                Logger;
        end
    end

end
