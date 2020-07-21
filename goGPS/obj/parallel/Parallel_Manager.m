%   CLASS Parallel_Manager
% =========================================================================
%
% DESCRIPTION
%   Master controller for parallel goGPS computation
%
% EXAMPLE
%   gom = Parallel_Manager
%
% FOR A LIST OF CONSTANTs and METHODS use doc Parallel_Manager
%
% NAMING CONVENTIONs
%   the master  is the parallel manager unique instance that manages slaves
%   a slave     is an idle process waiting to be activated
%   a worker    is an activated slave (it has been verified by the Master that it's alive and running)
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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

classdef Parallel_Manager < Com_Interface
    %% PROPERTIES CONSTANTS
    % ==================================================================================================================================================
    properties (Constant, GetAccess = public)
        ID = 'MASTER'
        
        MSG_KILLALL = 'KILLALL_'
        MSG_RESTART = 'RESTART_'
        MSG_ASKACK = 'ASKACK_'
        MSG_ASKWORK = 'ASKWORK_'
        
        MSG_DO = 'DO_'
        
        BRD_STATE = 'BRD_STATE_'
        BRD_SKY = 'BRD_SKY_'
        BRD_CMD = 'BRD_CMD_'
        BRD_REC = 'BRD_REC_'
    end
    
    properties (GetAccess = private, SetAccess = private)
        worker_id = {}; % list of active workers
        timeout = 25; % additional wait time for fast remote answers
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = private)
        % Concrete implementation.  See Singleton superclass.
        function this = Parallel_Manager(com_dir)
            % Core object creator
            this.id = this.ID;
            this.initComDir(com_dir);
            this.initHandles();
            this.log.addMarkedMessage('Creating goGPS Master Parallel Manager');
            if ~exist(this.getComDir, 'file')
                mkdir(this.getComDir);
            end
        end
    end
    
    %% METHOD DESTRUCTOR
    % ==================================================================================================================================================
    methods (Access = private)
        function delete(this)
            % delete
            Logger.getInstance.addMarkedMessage('Closing goGPS Master');
        end
    end
    %
    %% METHOD INTERFACE STATIC
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = getInstance(com_dir, destroy_this)
            % Get the persistent instance of the class
            % INPUT
            %   com_dir        comunication dir path
            %   destroy_this   flag, if true delete the persistent handler
            %
            % SYNTAX
            %   this.getInstance(com_dir, destroy_this)
            
            persistent unique_instance_gom__
            if nargin >= 2 && ~isempty(destroy_this) && destroy_this
                if ~isempty(unique_instance_gom__)
                    unique_instance_gom__.delete();
                    clear unique_instance_gom__
                end
            else
                if isempty(unique_instance_gom__)
                    if nargin < 1 || isempty(com_dir)
                        com_dir = Core.getState.getComDir;
                    end
                    this = Parallel_Manager(com_dir);
                    unique_instance_gom__ = this;
                else
                    this = unique_instance_gom__;
                    this.initHandles();
                    if nargin > 1 && ~isempty(com_dir)
                        this.initComDir(com_dir);
                    end
                end
            end
        end
        
        function killAll(com_dir)
            % Kill all the slaves
            %
            % SYNTAX
            %   Parallel_Manager.killAll(com_dir);
            if nargin >= 1
                gom = Parallel_Manager.getInstance(com_dir);
            else
                gom = Parallel_Manager.getInstance;
            end
            gom.setComDir(Core.getState.getComDir);
            gom.killThemAll
        end
        
        function restartWorkers(gom, verbose)
            % Force all the slaves to die and resurge
            %
            % SYNTAX
            %   Parallel_Manager.restartWorkers();
            if nargin < 1
                gom = Parallel_Manager.getInstance;
            end
            if nargin < 2 || isempty(verbose)
                verbose = true;
            end
            gom.resurgit([], verbose);
        end
        
        function requestSlaves(n_slaves, com_dir)
            % Start new sessions of matlab executing the slave worker
            %
            % INPUT
            %   n_slaves     number of slaves to launch
            %
            % SYNTAX
            %   this.requestSlaves(n_slaves, com_dir);
            
            t0 = tic();
            if nargin >= 2
                gom = Parallel_Manager.getInstance(com_dir);
            else
                gom = Parallel_Manager.getInstance;
            end
            gom.setComDir(Core.getState.getComDir);
            gom.log.addMessage(gom.log.indent(sprintf('Checking for slaves into "%s"', gom.getComDir)));
            n_living_slaves = gom.checkLivingSlaves();
            if n_living_slaves < n_slaves
                gom.createSlaves(n_slaves - n_living_slaves, gom.getComDir());
            end
            
            % Wait for slaves to appear (for max 25 seconds)
            if n_living_slaves < n_slaves
                n_living_slaves = gom.waitForSlaves(n_slaves, 25);
            end
            gom.log.addMarkedMessage(sprintf('At the moment I see %d living %s ready for processing - after %.2f seconds', n_living_slaves, iif(n_living_slaves == 1, 'slave', 'slaves'), toc(t0)));
        end
        
        function testSlaves(com_dir)
            % Check the number of slaves available
            %
            % INPUT
            %   n_slaves     number of slaves available
            %
            % SYNTAX
            %   this.testSlaves(com_dir);
            
            t0 = tic();
            if nargin >= 1
                gom = Parallel_Manager.getInstance(com_dir);
            else
                gom = Parallel_Manager.getInstance;
            end
            gom.log.addMessage(gom.log.indent(sprintf('Checking for slaves into "%s"', gom.getComDir)));
            n_living_slaves = gom.checkLivingSlaves();
            gom.resurgit([], false);
            gom.log.addMarkedMessage(sprintf('At the moment I see %d living %s ready for processing - after %.2f seconds', n_living_slaves, iif(n_living_slaves == 1, 'slave', 'slaves'), toc(t0)));
        end
    end
    
    methods (Static, Access = private)
        function createSlaves(n_slaves, com_dir, goGPS_folder)
            % Start new sessions of matlab executing the slave worker
            %
            % INPUT
            %   n_slaves      number of slaves to launch
            %   com_dir       directory of comunication
            %   goGPS_folder  goGPS source directory
            %
            % DEFAULT
            %   com_dir       fullfile(pwd, 'com')
            %   goGPS_folder  pwd
            %
            % SYNTAX
            %   this.createSlaves(n_slaves, com_dir, goGPS_folder);
            if nargin < 2 || isempty(com_dir)
                com_dir = fullfile(pwd, 'com');
            end
            if nargin < 3 || isempty(goGPS_folder)
                goGPS_folder = pwd;
            end
            
            slave_cmd = [...
                'cd ' goGPS_folder '; ' ...
                'addPathGoGPS;', ...
                'log = Logger.getInstance;' ...
                'log.setColorMode(0);' ...
                'log.setVerbosityLev(0);' ...
                'gos = Go_Slave.getInstance(''' com_dir ''');' ...
                'gos.live; exit'];
            
            if isunix
                if ismac
                    mat_exe = [matlabroot '/bin/maci64/matlab'];
                else
                    mat_exe = [matlabroot '/bin/matlab'];
                end
                run_cmd = [mat_exe ' -singleCompThread -nodisplay -nosplash -r "' slave_cmd '" &'];
            elseif ispc
                mat_exe = [matlabroot '/bin/matlab.exe'];
                % In windows I need to create a bat to be able to run different matlab in background
                fid = fopen('win_create_worker.bat','w');
                run_cmd = ['"' mat_exe '" -singleCompThread -nodisplay -nosplash -r "not_ready = true; while(not_ready); pause(0.1); jmf = com.mathworks.mde.desk.MLDesktop.getInstance.getMainFrame; if not(isempty(jmf)); not_ready = false; end; end; jmf.hide, ' slave_cmd '" &'];
                fwrite(fid, run_cmd);
                fclose(fid);
            end
            
            log = Core.getLogger();
            log.addMarkedMessage(sprintf('Creating %03d new %s\n - MATLAB executable path: %s\n - goGPS source folder:    %s\n - comunication folder:    %s', n_slaves, iif(n_slaves == 1, 'slave', 'slaves'), mat_exe, goGPS_folder, com_dir));
            log.newLine;
            
            for i = 1 : n_slaves
                log.addMessage(log.indent(sprintf('Creating slave worker %03d / %03d', i, n_slaves)));
                if ispc
                    system('start /b win_create_worker.bat');
                else
                    dos(run_cmd);
                end
                pause(0.5); % pause is necessary to avoid twin slaves (slaves with the same id)
            end
            
            if ispc
                % Under windows I created a file to be removed
                delete('win_create_worker.bat');
            end
        end
    end
    %
    %% METHOD INTERFACE
    % ==================================================================================================================================================
    methods (Access = public)
        function n_workers = getNumWorkers(this)
            % Get the active workers number
            %
            % SYNTAX
            %   n_workers = this.getNumWorkers();
            n_workers = numel(this.worker_id);
        end
        
        function n_slaves = checkLivingSlaves(this)
            % Check how many slaves already existing are alive
            %
            % SYNTAX
            %   n_living_slaves = this.checkLivingSlaves(n_wanted_slaves, max_time)
            
            % Check for already existing slaves stuck in wait
            % and restart them
            this.resurgit([], false);
            % Check for life
            n_searched = this.testWorkers();
            if isunix() && ~ismac()
                % On our Linux box, this update takes longer
                pause(0.5);
            end
            
            % Check for slaves
            slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_BORN Go_Slave.SLAVE_WAIT_PREFIX '*']));
            i = 0;
            while (n_searched > numel(slave_list)) && (i < 2)
                i = i + 1;
                pause(1);
                slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_BORN Go_Slave.SLAVE_WAIT_PREFIX '*']));
            end
            n_slaves = numel(slave_list);
            if n_slaves > 0
                this.log.addMessage(this.log.indent(sprintf('I already have %d living %s', n_slaves, iif(n_slaves == 1, 'slave', 'slaves'))));
            end
        end
        
        function n_workers = activateWorkers(this, flag_send_sky, id_rec2send)
            % Get the available workers,
            % Activate them by sending state and core_sky objects
            %
            % OUTPUT
            %   n_workers       number of available workers
            %   this.worker_id  (implicit) id of the available ready workers
            %
            % SYNTAX
            %   n_workers = this.activateWorkers(flag_send_sky, id_rec2send);
            
            % Checking available slaves
            t0 = tic();
            if nargin < 1 || isempty(flag_send_sky)
                flag_send_sky = true;
            end
            this.deleteMsg('*'); % delete all master massages
            this.deleteMsg(Go_Slave.MSG_DIE, true);
            this.deleteMsg(Go_Slave.MSG_ACK, true);
            this.deleteMsg(Go_Slave.MSG_JOBREADY, true);
            this.deleteMsg([Go_Slave.MSG_BORN, Go_Slave.SLAVE_READY_PREFIX '*'], true);
            
            slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_BORN '*']));
            n_slaves = numel(slave_list);
            this.log.addMessage(this.log.indent(sprintf('%d slaves found', n_slaves)));
            
            % Check the life of slaves (ask ack)
            for w = 1 : n_slaves
                slave_id = regexp(slave_list(w).name, [Go_Slave.SLAVE_WAIT_PREFIX, '[0-9]*'], 'match');
                slave_id = slave_id{1};
                msg = [slave_id, '_' this.MSG_ASKACK];
                this.sendMsg(msg, sprintf('"%s" are you ready to work?', slave_id));
            end
            
            % Wait 2 seconds max for slave answers or till all the slaves have responded
            n_workers = 0;
            elapsed_time = 0;
            while elapsed_time < (2 + this.timeout) && (n_workers < n_slaves)
                pause(0.1);
                elapsed_time = elapsed_time + 0.1;
                slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_ACK '*']));
                n_workers = numel(slave_list);
            end
            this.log.addMessage(this.log.indent(sprintf('%d workers found', n_workers)));
            this.deleteMsg([Go_Slave.MSG_BORN, Go_Slave.SLAVE_WAIT_PREFIX '*'], true);
            this.deleteMsg([Parallel_Manager.MSG_ASKACK, Go_Slave.SLAVE_WAIT_PREFIX '*'], true);
            
            % Activate answering workers
            this.worker_id = {};
            if flag_send_sky
                slave_type = Go_Slave.SLAVE_TARGET_PREFIX;
            else
                slave_type = Go_Slave.SLAVE_SESSION_PREFIX;
            end
            
            for w = 1 : n_workers
                delete(fullfile(this.getComDir, slave_list(w).name));
                slave_id = regexp(slave_list(w).name, [Go_Slave.SLAVE_WAIT_PREFIX, '[0-9]*'], 'match');
                slave_id = slave_id{1};
                this.worker_id{w} = [Go_Slave.SLAVE_READY_PREFIX slave_type sprintf('%03d_', w)];
                msg = [slave_id, '_' this.MSG_ASKWORK this.worker_id{w} '_'];
                this.sendMsg(msg, sprintf('"%s" you will be known as "%s', slave_id, this.worker_id{w}));
            end
            
            if n_workers == 0
                this.log.addError('I need slaves! Please create some slaves!\nsee Parallel_Manager.createWorkers(n);');
            else
                if ~isempty(id_rec2send)
                    this.sendReceiver(id_rec2send);
                end
                this.sendState(flag_send_sky);
                if flag_send_sky
                    this.sendSkyData();
                end
                n_workers = this.waitForWorkerAck(n_workers);
                this.deleteMsg('*');
                this.deleteMsg([Go_Slave.MSG_ACK, Go_Slave.SLAVE_READY_PREFIX, slave_type '*'], true);
                delete(fullfile(this.getComDir, 's*.mat'));
                delete(fullfile(this.getComDir, 'rec_list.mat'));
            end
            this.log.addMarkedMessage(sprintf('Parallel init took %.3f seconds', toc(t0)));
        end
        
        function orderProcessing(this, cmd_list, par_type, list)
            % order to the slaves to execute some work
            %
            % SINTAX
            %   orderProcessing(this, cmd_list, sss_list, trg_list)
            %
            % EXAMLE
            %   gom.orderProcessing(par_cmd_list, sss_list, trg_list{l});
            
            this.sendCommandList(cmd_list);
            
            %             switch par_type
            %                 case 1
            %                     sss_list = list;
            %                 case 2
            %                     trg_list = list;
            %                     % Get info on the curret session
            %                     [cur_sss, sss_list, sss_id] = Core.getCurrentSession();
            %                     log = core.getLogger();
            %             end
            
            worker_stack = this.worker_id;
            worker2job = zeros(size(worker_stack)); % keep track of which job is assigned to which worker
            worker2jobstart = zeros(size(worker_stack),'uint64');
            missing_job = list;
            completed_job = [];
            active_jobs = 0;
            
            while numel(completed_job) < numel(list) && ~isempty(worker_stack)
                % Send work to slaves
                parallel_job = 1 : min(numel(worker_stack), numel(missing_job));
                t = 0;
                for w = parallel_job
                    t = t + 1;
                    
                    % send an order to a worker
                    msg = [worker_stack{w}, this.MSG_DO num2str(missing_job(t), '%04d') '_'];
                    if par_type == 1
                        this.sendMsg(msg, sprintf('"%s" process session %d', worker_stack{w}(1 : end-1), missing_job(t)));
                    else
                        this.sendMsg(msg, sprintf('"%s" process rec %d', worker_stack{w}(1 : end-1), missing_job(t)));
                    end
                    worker2job(str2double(worker_stack{w}(10:12))) = missing_job(t);
                    worker2jobstart(str2double(worker_stack{w}(10:12))) = tic();
                    active_jobs = active_jobs + 1;
                end
                missing_job(1 : t) = []; % remove the currently executing jobs
                worker_stack(parallel_job) = []; %  remove the active worker from the worker stack
                
                % Before checking for finished job I can maybe start computing the orbits for the next session!!!
                % EXPERIMENTAL: internal usage only for target parallelism
                %test_parallel_session_load = false;
                %if (sss_id < numel(sss_list)) && test_parallel_session_load
                %    vl = log.getVerbosityLev(); log.setVerbosityLev(log.WARNING_VERBOSITY_LEV);
                %    core.prepareSession(sss_list(sss_id + 1), true);
                %    log.setVerbosityLev(vl);
                %else
                % wait for complete job
                pause(0.1);
                %end
                
                % check for complete jobs
                [active_jobs, completed_job, worker_stack] = this.waitCompletedJob(active_jobs, completed_job, worker_stack, worker2job, worker2jobstart, par_type);
            end
            if isempty(worker_stack)
                Core.getLogger.addError('All workers has been lost. Stopping parallel execution');
            end
            delete(fullfile(this.getComDir, 'cmd_list.mat'));
            this.sendMsg(this.MSG_RESTART, 'My slaves, your job is done!\n        > There is no time for resting!\n        > Wait for other jobs!');
        end
        
        function importParallelSessions(this)
            core = Core.getCurrentCore();
            log = core.getLogger();
            state = Core.getState();
            n_rec = core.state.getRecCount;
            
            % get available computed sessions
            sss_path = fullfile(this.getComDir, ['job*' Go_Slave.SLAVE_READY_PREFIX Go_Slave.SLAVE_SESSION_PREFIX '*.mat']);
            sss_file = dir(sss_path);
            % Create receiver list
            if isempty(core.rec)
                clear rec;
                for r = 1 : n_rec
                    rec(r) = GNSS_Station(state.getDynMode() == 0); %#ok<AGROW>
                end
                core.rec = rec;
                clear rec;
            end
            
            % Assign all the new work to current rec;
            for s = 1 : numel(sss_file)
                %%
                sss_id = regexp(sss_file(s).name, '(?<=job)[0-9]*', 'match', 'once');
                
                if ~isempty(sss_id)
                    sss_id = str2double(sss_id);
                    w_id = regexp(sss_file(s).name, ['(?<=' Go_Slave.SLAVE_READY_PREFIX Go_Slave.SLAVE_SESSION_PREFIX ')[0-9]*'], 'match', 'once');
                    if ~isempty(w_id)
                        w_id = str2double(w_id);
                    end
                    not_corrupted = true;
                    try
                        tmp = load(fullfile(this.getComDir, sss_file(s).name));
                    catch
                        log.addWarning(sprintf('File %s seems corrupted', sss_file(s).name));
                        not_corrupted = false;
                    end
                    if not_corrupted
                        core.state.setCurSession(sss_id); % load the current session number
                        % Check that all the results are present in the session file
                        if isfield(tmp, 'rec') && (numel(tmp.rec) == n_rec) && isfield(tmp, 'atmo')
                            % Import atmosphere as computed by rec
                            core.setAtmosphere(tmp.atmo);
                            log.addMessage(log.indent(sprintf('%s - Importing session %d computed by worker %d', GPS_Time.now.toString('HH:MM:SS'), sss_id, w_id)));
                            drawnow;
                            for r = 1 : n_rec
                                if core.rec(r).isEmpty
                                    % The first time (computed session) import the entire object
                                    core.rec(r) = tmp.rec(r);
                                else
                                    core.rec(r).work = tmp.rec(r).work;
                                    core.rec(r).work.parent = core.rec(r); % reset handler to the parent object
                                end
                                core.rec(r).initHandles(); % reset other handles
                                core.rec(r).work.initHandles(); % reset other handles
                                core.rec(r).work.pushResult;
                            end
                        else
                            log.addWarning(sprintf('Session %d have been computed by worker %d but seems empty or corrupted', sss_id, w_id));
                        end
                    end
                    
                end
            end
            log.addMarkedMessage('All the parallel session present in the com folder have been imported');
            
            % delete all the imported files!
            delete(sss_path);
        end
        
    end
    
    methods (Access = private)
        function [n_old_slaves, slave_list] = getNumSlaves(this)
            % Check how many slaves already exist
            %
            % SYNTAX
            %   n_living_slaves = this.waitForSlaves(n_wanted_slaves, max_time)
            
            slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_BORN Go_Slave.SLAVE_WAIT_PREFIX '*']));
            n_old_slaves = numel(slave_list);
        end
        
        function n_living_slaves = waitForSlaves(this, n_wanted_slaves, max_time)
            % Check how many slaves exist, and wait for them
            %
            % INPUT
            %   n_wanted_slaves    number of slaves to wait for
            %   max_time           time in seconds of maximum wait
            %
            % SYNTAX
            %   n_living_slaves = this.waitForSlaves(n_wanted_slaves, max_time)
            
            % Wait for slaves to appear (for max 1 second)
            n_living_slaves = 0;
            reset_count = 0;
            twait = tic;
            while (n_living_slaves < n_wanted_slaves) && (toc(twait) < max_time)
                slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_BORN Go_Slave.SLAVE_WAIT_PREFIX '*']));
                n_living_slaves = numel(slave_list);
                reset_count = this.pause(0.05, reset_count);
            end
            this.pause(0, -1);
        end
        
        function die(this)
            % Kill the go Master
            %
            % SYNTAX:
            %   this.die();
            this.log.addMarkedMessage('Bye bye!');
            Parallel_Manager.getInstance([], true)
        end
        
        function n_slaves = testWorkers(this)
            % Get the available workers,
            % Activate them by sending state and core_sky objects
            %
            % OUTPUT
            %   n_workers       number of available workers
            %   this.worker_id  (implicit) id of the available ready workers
            %
            % SYNTAX
            %   n_workers = this.testWorkers();
            
            % Checking available slaves
            t0 = tic();
            this.deleteMsg('*'); % delete all master messages
            this.deleteMsg(Go_Slave.MSG_DIE, true);
            this.deleteMsg(Go_Slave.MSG_ACK, true);
            this.deleteMsg([Go_Slave.MSG_BORN, Go_Slave.SLAVE_READY_PREFIX '*'], true);
            
            [n_slaves, slave_list] = this.getNumSlaves();
            this.log.addMessage(this.log.indent(sprintf('%d %s found', n_slaves, iif(n_slaves == 1, 'slave', 'slaves'))));
            
            if n_slaves > 0
                % Check the life of slaves (ask ack)
                for w = 1 : n_slaves
                    slave_id = regexp(slave_list(w).name, [Go_Slave.SLAVE_WAIT_PREFIX, '[0-9]*'], 'match');
                    slave_id = slave_id{1};
                    msg = [slave_id, '_' this.MSG_ASKACK];
                    this.sendMsg(msg);
                end
                
                % Wait 1 seconds for slave answers or till all the slaves have responded
                n_workers = 0;
                elapsed_time = 0;
                while elapsed_time < (1 + this.timeout) && (n_workers < n_slaves)
                    pause(0.1);
                    elapsed_time = elapsed_time + 0.1;
                    slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_ACK '*']));
                    n_workers = numel(slave_list);
                end
                this.log.addMessage(this.log.indent(sprintf('%d skilled %s found', n_workers, iif(n_slaves == 1, 'slave', 'slaves'))));
                this.deleteMsg([Go_Slave.MSG_BORN, Go_Slave.SLAVE_WAIT_PREFIX '*'], true);
                this.deleteMsg([Parallel_Manager.MSG_ASKACK, Go_Slave.SLAVE_WAIT_PREFIX '*'], true);
                
                % Dispose trash
                this.worker_id = {};
                for w = 1 : n_workers
                    delete(fullfile(this.getComDir, slave_list(w).name));
                end
                % Restart active workers
                if n_workers > 0
                    this.resurgit(n_workers, false);
                end
                n_slaves = n_workers;
            end
        end
        
        function sendCommandList(this, cmd_list)
            % Send command list to the slaves
            %
            % SYNTAX
            %   this.sendCommandList()
            %
            
            % Save state on file
            state = Core.getCurrentSettings();
            save(fullfile(this.getComDir, 'cmd_list.mat'), 'cmd_list');
            this.sendMsg(this.BRD_CMD, 'Broadcast state');
        end
        
        function sendState(this, slave_type)
            % Send state to the slaves
            %
            % INPUT
            %   slave_type  1 for target workers
            %               0 for session workers
            %
            % SYNTAX
            %   this.sendState(slave_type)
            %
            
            % Save state on file
            geoid = Core.getRefGeoid();
            state = Core.getCurrentSettings();
            atx = Core.getAntennaManager();
            cur_session = Core.getCurrentSession();
            [rin_list, met_list] = Core.getRinLists();
            save(fullfile(this.getComDir, 'state.mat'), 'geoid', 'state', 'atx', 'cur_session', 'rin_list', 'met_list', 'slave_type');
            this.sendMsg(this.BRD_STATE, 'Broadcast state');
        end
        
        function sendSkyData(this)
            % Send core sky to the slaves
            %
            % SYNTAX
            %   this.sendSkyData()
            %
            
            % Save core sky on file
            sky = Core.getCoreSky();
            atmo = Core.getAtmosphere();
            mn = Core.getMeteoNetwork();
            atx = Core.getAntennaManager();
            save(fullfile(this.getComDir, 'sky.mat'), 'sky', 'atmo', 'mn');
            this.sendMsg(this.BRD_SKY, 'Broadcast core sky');
        end
        
        function sendReceiver(this, rec_num)
            % Send a receiver to the slave
            %
            % SYNTAX
            %   this.sendReceiver(rec_num)
            %
            
            % Note:
            %   rec_num{1} - target loop
            %   rec_num{2} - pass receivers (they actually pass only the work obj)
            %   rec_num{3} - work receivers
            %   rec_num{4} - out receivers
            
            % Save work receivers to pass to all the slaves
            rec_list = Core.getRecList();
            log = Core.getLogger();
            
            id_work = min(numel(rec_list), max(1, unique(setdiff([rec_num{2} rec_num{3}], 0))));
            id_out = min(numel(rec_list), max(1, unique(setdiff(rec_num{4}, 0))));
            id_target = min(numel(rec_list), max(1, unique(setdiff(rec_num{4}, 0))));
            
            % Pass receivers basic informations
            clear rec_info
            for r = 1 : numel(rec_list)
                rec_info(r) = struct('marker_name', rec_list(r).marker_name, ...
                    'marker_type', rec_list(r).marker_type, ...
                    'number', rec_list(r).number, ...
                    'type', rec_list(r).type, ...
                    'version', rec_list(r).version, ...
                    'observer', rec_list(r).observer, ...
                    'agency', rec_list(r).agency, ...
                    'ant_serial', rec_list(r).ant_serial, ...
                    'ant_type', rec_list(r).ant_type, ...
                    'ant_delta_h', rec_list(r).ant_delta_h, ...
                    'ant_delta_en', rec_list(r).ant_delta_en, ...
                    ... % 'ant_mp', rec_list(r).ant_mp, ...
                    'static', rec_list(r).static);
            end
            if numel(rec_list) == 0
                rec_info = struct();
            end
            
            if isempty(id_work)
                rec_work = [];
            else
                rec_work = [rec_list(id_work).work]; %#ok<NASGU>
            end
            if isempty(id_out)
                rec_out = [];
            else
                rec_out = [rec_list(id_out).out]; %#ok<NASGU>
            end
            log.addMessage(log.indent(' - Saving receivers to pass as broadcast'));
            save(fullfile(this.getComDir, 'rec_list.mat'), 'rec_info', 'rec_work', 'rec_out', 'id_target', 'id_work', 'id_out');
            this.sendMsg(this.BRD_REC, 'Broadcast receiver');
            
            % These other files are needed only after the activation of the workers
            % Save work receivers to pass to the respective slave slaves
            if ismember(0, [rec_num{2} rec_num{3}])
                for r = unique(rec_num{1})
                    log.addMessage(log.indent(sprintf(' - Saving receiver %d workspace for my slaves', r)));
                    rec_work = rec_list(r).work;
                    save(fullfile(this.getComDir, sprintf('rec_work_%04d.mat', r)), 'rec_work');
                end
            end
            
            % Save work receivers to pass to the respective slave slaves
            if ismember(0, rec_num{4})
                for r = unique(rec_num{1})
                    log.addMessage(log.indent(sprintf(' - Saving receiver %d output for my slaves', r)));
                    rec_out = rec_list(r).out;
                    save(fullfile(this.getComDir, sprintf('rec_out_%04d.mat', r)), 'rec_out');
                end
            end
        end
        
        function n_workers = waitForWorkerAck(this, n_slaves)
            % Wait for the workers to load the state and core sky
            %
            % INPUT
            %   n_slaves    number of slave process ready to work
            %
            % OUTPUT
            %   n_workers   number of workers ready to work
            %
            % SYNTAX
            %   n_workers = this.waitForWorkerAck(n_slaves)
            %
            
            % Keep workers that have loaded the files
            % Wait 10 seconds for slave answers or till all the slaves have responded
            n_workers = 0;
            elapsed_time = 0;
            while elapsed_time < (10 + this.timeout) && (n_workers < n_slaves)
                pause(0.1);
                elapsed_time = elapsed_time + 0.1;
                slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_ACK '*']));
                n_workers = numel(slave_list);
            end
            this.log.addMarkedMessage(sprintf('%d workers ready', n_workers));
            this.deleteMsg([Go_Slave.MSG_ACK, Go_Slave.SLAVE_READY_PREFIX '*'], true);
        end
        
        function [active_jobs, completed_job, worker_stack] = waitCompletedJob(this, active_jobs, completed_job, worker_stack, worker2job, worker2jobstart, par_type)
            % Wait for the workers for sending new jobs
            %
            % INPUT
            %   n_slaves    number of slave process ready to work
            %
            % OUTPUT
            %   n_workers   number of workers ready to work
            %
            % SYNTAX
            %   n_workers = this.waitCompletedJob(n_slaves)
            %
            
            if active_jobs > 0
                core = Core.getCurrentCore();
                
                elapsed_time = 0;
                n_job_done = 0;
                while (n_job_done < 1) || (active_jobs == 0)
                    pause(0.1);
                    % To send jobs if anything goes wrong use:
                    % this.sendMsg('WORKER_S_001_DO_0020_', 'helooo');
                    % that line ask to worker 001 for the job 0020
                    elapsed_time = elapsed_time + 0.1;
                    % Search for finished jobs
                    slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_JOBREADY '*']));
                    n_job_done = numel(slave_list);
                    for j = 1 : numel(slave_list)
                        worker_id = regexp(slave_list(j).name, [ Go_Slave.SLAVE_READY_PREFIX '[T|S]_[0-9]*'], 'match', 'once');
                        % get the result stored into jobXXXX_WORKER_T|S_YYY.mat
                        job_file = dir(fullfile(this.getComDir, ['job*' worker_id '.mat']));
                        switch par_type
                            case 1 % if it's a session parallelism:
                                worker_id_num = str2double(regexp(worker_id, [ '(?<=' Go_Slave.SLAVE_READY_PREFIX '[T|S]_)[0-9]*'], 'match', 'once'));
                                job_id = worker2job(worker_id_num);
                                if isempty(job_file)
                                    this.log.addError(sprintf('Worker %d completed its job with no results', worker_id_num));
                                else
                                    this.log.addMessage(sprintf('Worker %d completed its job %d', worker_id_num, job_id));
                                end
                            case 2 % if it's a target parallelism:
                                % Start to import the processed target while data is still computed
                                
                                if isempty(job_file)
                                    tmp = GNSS_Station(Core.getState.getDynMode() == 0);
                                    worker_id_num = str2double(regexp(worker_id, [ '(?<=' Go_Slave.SLAVE_READY_PREFIX '[T|S]_)[0-9]*'], 'match', 'once'));
                                    this.log.addError(sprintf('Worker %d completed its job with no results', worker_id_num));
                                    job_id = worker2job(worker_id_num);
                                    tmp.work.flag_currupted = true;
                                else
                                    job_id = str2double(regexp(job_file(1).name, '(?<=job)[0-9]*', 'match', 'once'));
                                    try
                                        tmp = load(fullfile(this.getComDir(), job_file(1).name));
                                    catch
                                        pause(1); % due to sync problems the file could not be immediately ready,
                                        % try to read it again
                                        tmp = load(fullfile(this.getComDir(), job_file(1).name));
                                    end
                                    s0 = tmp.rec.work.sat.res.getStd();
                                    if s0 * 1e2 > 2
                                        this.log.addWarning(sprintf('s0 = %.3f of the residuals for parallel job %d (session %d)', s0, job_id, tmp.rec.state.getCurSession()));
                                    else
                                        this.log.addMessage(this.log.indent(sprintf('job %d residuals s0 = %.3f from parallel execution (session %d)', job_id, s0, tmp.rec.state.getCurSession), 9));
                                    end
                                    if core.rec(job_id).out.isEmpty
                                        % import all
                                        tmp.rec.out = core.rec(job_id).out;
                                        % Detect if the parallel GNSS_Station is empty
                                        par_sta_is_empty = (strcmp(tmp.rec.marker_name, 'unknown') && strcmp(tmp.rec.number, '000') && strcmp(tmp.rec.version, '000'));
                                        local_sta_is_empty = (strcmp(core.rec(job_id).marker_name, 'unknown') && strcmp( core.rec(job_id).number, '000') && strcmp(core.rec(job_id).version, '000'));
                                        if ~par_sta_is_empty || local_sta_is_empty
                                            % import the entire receiver
                                            core.rec(job_id) = tmp.rec;
                                        else
                                            % Import only work and out
                                            core.rec(job_id).work = tmp.rec.work;
                                            core.rec(job_id).out = tmp.rec.out;
                                        end
                                        % Remember to out (Luke) who's the father!
                                        % Luke, I'm your father!!!
                                        core.rec(job_id).out.parent = core.rec(job_id);
                                        core.rec(job_id).out.initHandles();
                                        % relink singletons
                                        core.rec(job_id).log = Core.getLogger;
                                        core.rec(job_id).state = Core.getState;
                                        % import results in out
                                    else
                                        % import only work
                                        tmp.rec.work.parent = core.rec(job_id);
                                        tmp.rec.work.initHandles();
                                        core.rec(job_id).work = tmp.rec.work;
                                    end
                                    % Remember to work (Leia) who's the father!
                                    core.rec(job_id).work.parent = core.rec(job_id);
                                    if core.rec(job_id).work.flag_currupted
                                        this.log.addError(sprintf('Something bad appened in parallel job %d - processing corrupted', job_id));
                                    end
                                    this.log.addMessage(this.log.indent(sprintf('Deleting %s',fullfile(this.getComDir(), job_file(1).name))),200);
                                    delete(fullfile(this.getComDir(), job_file(1).name));
                                end
                        end
                        % core.rec(job_id).work.pushResult();
                        this.deleteMsg([Go_Slave.MSG_JOBREADY, worker_id], true);
                        
                        
                        completed_job = [completed_job; job_id]; %#ok<AGROW>
                        worker_stack = [worker_stack {[worker_id '_']}]; %#ok<AGROW>
                        
                        
                    end
                    % check if executiong time has not passed max time
                    % (sometimes a parallel job die without notice)
                    for w = 1:length(worker2jobstart)
                        if ~ismember(worker2job(w),completed_job)
                            if toc(worker2jobstart(w)) > Core.getCurrentSettings.getMaxExecutionTimePar()
                                completed_job = [completed_job; worker2job(w)]; %#ok<AGROW>
                                this.log.addWarning(sprintf('Worker %d is taking too long to complete, signaling job %d as complete, results will be missing',w, numel(completed_job)));
                                this.log.addWarning(sprintf('Removing worker %d from workers'' pool',w));
                                n_job_done = n_job_done +1;
                            end
                        end
                    end
                end
                
                active_jobs = active_jobs - n_job_done;
                this.log.addMarkedMessage(sprintf('%d jobs completed @ %s', numel(completed_job), GPS_Time.now.toString('yyyy-mm-dd HH:MM:SS')));
                this.deleteMsg([Go_Slave.MSG_ACK, Go_Slave.SLAVE_READY_PREFIX '*'], true);
            end
        end
        
        function cleanDummy(obj, rootname)
            dos(['del ' obj.COM_DIR 'dummy' rootname '* > nil']);
        end
        
        function rmComDir(this)
            % Delete the comunication dir
            %
            % SYNTAX:
            %   this.rmComDir();
            if exist(this.getComDir, 'file')
                rmdir(this.getComDir, 's');
            end
        end
        
        function killThemAll(this)
            % Send a kill message to close all the parallel workers/slaves
            %
            % SYNTAX
            %   this.killThemAll
            this.sendMsg(this.MSG_KILLALL, 'As the mad king said: Kill them all!!!');
            pause(1);
            this.deleteMsg(Go_Slave.MSG_DIE, true);
            this.deleteMsg(Go_Slave.MSG_BORN, true);
            this.deleteMsg('*');
            delete(fullfile(this.COM_DIR,'*.mat'));
        end
        
        function resurgit(this, n_workers, verbosity)
            % Try to revive crashed and idle workers -> free them
            % and put them back into the stack of slaves
            %
            % SYNTAX
            %   this.resurgit(n_workers, verbosity)
            if nargin < 3 || isempty(verbosity)
                verbosity = true;
            end
            n_old_slaves = 0;
            if nargin < 2 || isempty(n_workers)
                % Check for already existing slaves
                n_old_slaves = this.getNumSlaves();
                
                % Check for workers stuck in the work loop
                % (they are not doing anything, probably they are in this state due to a crash of the old Master)
                worker_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_BORN Go_Slave.SLAVE_READY_PREFIX  '*']));
                n_workers = numel(worker_list);
                this.log.addMessage(this.log.indent(sprintf('%d idler worker found', n_workers)));
                if n_workers > 0
                    % If I found some stuck workers I can try to revive it!
                    this.log.addMessage(this.log.indent('Unuseful scum restart your life!'));
                end
            end
            this.sendMsg(this.MSG_RESTART, iif(verbosity, 'There is no time for resting!\n          Wait for other jobs!', ''));
            if n_workers == 0
                % Wait for slaves whos id have been deleted, but they're still alive
                pause(0.2);
            end
            n_workers = n_workers + n_old_slaves;
            this.waitForSlaves(n_workers, 1);
            this.deleteMsg(Go_Slave.MSG_DIE, true);
            this.deleteMsg('*');
        end
    end
    %
    %% METHOD INIT
    % ==================================================================================================================================================
    methods
        function initHandles(this)
            this.log = Core.getLogger();
        end
    end
end
