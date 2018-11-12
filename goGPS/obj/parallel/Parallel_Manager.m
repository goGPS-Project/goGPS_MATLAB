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
%    |___/                    v 0.999.0 - nightly
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
            this.log.addMarkedMessage('Closing goGPS Master');
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
                        com_dir = fullfile(pwd, 'com');
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
        
        function killAll()
            % Kill all the slaves
            %
            % SYNTAX
            %   Parallel_Manager.killAll();
            gom = Parallel_Manager.getInstance;
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
        
        function requestSlaves(n_slaves)
            % Start new sessions of matlab executing the slave worker
            %
            % INPUT
            %   n_slaves     number of slaves to launch
            %
            % SYNTAX
            %   this.requestSlaves(n_slaves);
            
            t0 = tic();
            gom = Parallel_Manager.getInstance;
            n_living_slaves = gom.checkLivingSlaves();
            if n_living_slaves < n_slaves
                gom.createSlaves(n_slaves - n_living_slaves, gom.getComDir());
            end
            
            % Wait for slaves to appear (for max 20 seconds)
            if n_living_slaves < n_slaves
                n_living_slaves = gom.waitForSlaves(n_slaves, 20);
            end
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
            %   this.gimmeWorkers(n_workers, com_dir, goGPS_folder);
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
                run_cmd = [mat_exe ' -singleCompThread -nodisplay -nosplashmaxNumCompThreadsmaxNumCompThreads -r "' slave_cmd '" &'];
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
            this.testWorkers();
            
            % Check for slaves
            slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_BORN Go_Slave.SLAVE_WAIT_PREFIX '*']));
            n_slaves = numel(slave_list);
            if n_slaves > 0
                this.log.addMessage(this.log.indent(sprintf('I alredy have %d living %s', n_slaves, iif(n_slaves == 1, 'slave', 'slaves'))));
            end
        end
        
        function n_workers = activateWorkers(this, id_rec2pass)
            % Get the available workers,
            % Activate them by sending state and core_sky objects
            %
            % OUTPUT
            %   n_workers       number of available workers
            %   this.worker_id  (implicit) id of the available ready workers
            %
            % SYNTAX
            %   n_workers = this.activateWorkers();
            
            % Checking available slaves
            t0 = tic();
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
            
            % Wait 1 seconds for slave answers or till all the slaves have responded
            n_workers = 0;
            elapsed_time = 0;
            while elapsed_time < 1 && (n_workers < n_slaves)
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
            for w = 1 : n_workers
                delete(fullfile(this.getComDir, slave_list(w).name));
                slave_id = regexp(slave_list(w).name, [Go_Slave.SLAVE_WAIT_PREFIX, '[0-9]*'], 'match');
                slave_id = slave_id{1};
                this.worker_id{w} = [Go_Slave.SLAVE_READY_PREFIX sprintf('%03d_', w)];
                msg = [slave_id, '_' this.MSG_ASKWORK this.worker_id{w} '_'];
                this.sendMsg(msg, sprintf('"%s" you will be known as "%s', slave_id, this.worker_id{w}));
            end
            
            if n_workers == 0
                this.log.addError('I need slaves! Please create some slaves!\nsee Parallel_Manager.createWorkers(n);');
            else
                if ~isempty(id_rec2pass)
                    this.sendReceiver(id_rec2pass);
                end
                this.sendState();
                this.sendSkyData();
                n_workers = this.waitForWorkerAck(n_workers);
                this.deleteMsg('*');
                this.deleteMsg([Go_Slave.MSG_ACK, Go_Slave.SLAVE_READY_PREFIX '*'], true);
                delete(fullfile(this.getComDir, '*.mat'));
            end
            this.log.addMarkedMessage(sprintf('Parallel init took %.3f seconds', toc(t0)));
        end
        
        function orderProcessing(this, cmd_list, trg_list)
            % order to the slaves to execute some work
            %
            % SINTAX
            %   orderProcessing(this, cmd_list, trg_list)
            %
            % EXAMLE
            %   gom.orderProcessing(par_cmd_list, trg_list{l});
            this.sendCommandList(cmd_list);
            
            worker_stack = this.worker_id;
            
            missing_job = trg_list;
            completed_job = [];
            active_jobs = 0;
            while numel(completed_job) < numel(trg_list)
                % Send work to slaves
                parallel_job = 1 : min(numel(worker_stack), numel(missing_job));
                t = 0;
                for w = parallel_job
                    t = t + 1;
                    
                    % send an order to a worker
                    msg = [worker_stack{w}, this.MSG_DO num2str(missing_job(t), '%04d') '_'];
                    this.sendMsg(msg, sprintf('"%s" process rec %d', worker_stack{w}(1 : end-1), missing_job(t)));
                    active_jobs = active_jobs + 1;
                end
                missing_job(1 : t) = []; % remove the currently executing jobs
                worker_stack(parallel_job) = []; %  remove the active worker from the worker stack
                
                % wait for complete job
                pause(0.1);
                
                % check for complete jobs
                [active_jobs, completed_job, worker_stack] = this.waitCompletedJob(active_jobs, completed_job, worker_stack);
            end
            delete(fullfile(this.getComDir, 'cmd_list.mat'));
            this.sendMsg(this.MSG_RESTART, 'My slaves, your job is done!\n        > There is no time for resting!\n        > Wait for other jobs!');
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
            %   n_workers = this.activateWorkers();
            
            % Checking available slaves
            t0 = tic();
            this.deleteMsg('*'); % delete all master massages
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
                while elapsed_time < 1 && (n_workers < n_slaves)
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
            gc = Core.getGlobalConfig();
            state = gc.getCurrentSettings();
            save(fullfile(this.getComDir, 'cmd_list.mat'), 'cmd_list');
            this.sendMsg(this.BRD_CMD, 'Broadcast state');
        end
        
        function sendState(this)
            % Send state to the slaves
            %
            % SYNTAX
            %   this.sendState()
            %
            
            % Save state on file
            gc = Core.getGlobalConfig();
            geoid = gc.getRefGeoid();
            state = gc.getCurrentSettings();
            cur_session = Core.getCurrentSession();
            [rin_list, met_list] = Core.getRinLists();
            save(fullfile(this.getComDir, 'state.mat'), 'geoid', 'state', 'cur_session', 'rin_list', 'met_list');
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
            save(fullfile(this.getComDir, 'sky.mat'), 'sky', 'atmo', 'mn');
            this.sendMsg(this.BRD_SKY, 'Broadcast core sky');
        end
        
        function sendReceiver(this, rec_num)
            % Send a receiver to the slave
            %
            % SYNTAX
            %   this.sendReceiver(rec_num)
            %
            
            % Save state on file
            rec_list = Core.getRecList();
            rec_work = [rec_list(rec_num).work]; %#ok<NASGU>
            save(fullfile(this.getComDir, 'rec_list.mat'), 'rec_work', 'rec_num');
            this.sendMsg(this.BRD_REC, 'Broadcast receiver');
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
            while elapsed_time < 10 && (n_workers < n_slaves)
                pause(0.1);
                elapsed_time = elapsed_time + 0.1;
                slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_ACK '*']));
                n_workers = numel(slave_list);
            end
            this.log.addMarkedMessage(sprintf('%d workers ready', n_workers));
            this.deleteMsg([Go_Slave.MSG_ACK, Go_Slave.SLAVE_READY_PREFIX '*'], true);
        end
        
        function [active_jobs, completed_job, worker_stack] = waitCompletedJob(this, active_jobs, completed_job, worker_stack)
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
                    elapsed_time = elapsed_time + 0.1;
                    % Search for finished jobs
                    slave_list = dir(fullfile(this.getComDir, [Go_Slave.MSG_JOBREADY '*']));
                    n_job_done = numel(slave_list);
                    for j = 1 : numel(slave_list)
                        worker_id = regexp(slave_list(j).name, [ Go_Slave.SLAVE_READY_PREFIX '[0-9]*'], 'match', 'once');
                        % get the result stored into jobXXXX_WORKER_YYYY.mat
                        job_file = dir(fullfile(this.getComDir, ['job*' worker_id '.mat']));
                        job_id = str2double(regexp(job_file(1).name, '(?<=job)[0-9]*', 'match', 'once'));
                        tmp = load(fullfile(this.getComDir(), job_file(1).name));
                        if std(zero2nan(tmp.rec.work.sat.res(:)), 'omitnan') * 1e2 > 2
                            this.log.addWarning(sprintf('s0 = %.3f of the residuals for job %d', std(zero2nan(tmp.rec.work.sat.res(:)), 'omitnan'), job_id));
                        end
                        if core.rec(job_id).out.isEmpty
                            % import all
                            tmp.rec.out = core.rec(job_id).out;
                            tmp.rec.out.parent = tmp.rec;
                            tmp.rec.out.initHandles();
                            core.rec(job_id) = tmp.rec;
                            % relink singletons
                            core.rec(job_id).log = Core.getLogger;
                            core.rec(job_id).state = Core.getState;
                            % import results in out
                        else
                            % import only work
                            tmp.rec.work.initHandles();
                            core.rec(job_id).work = tmp.rec.work;
                            core.rec(job_id).work.parent = core.rec(job_id);
                        end
                        %core.rec(job_id).work.pushResult();
                        this.deleteMsg([Go_Slave.MSG_JOBREADY, worker_id], true);
                        delete(fullfile(this.getComDir(), job_file(1).name));
                        completed_job = [completed_job; job_id]; %#ok<AGROW>
                        worker_stack = [worker_stack {[worker_id '_']}]; %#ok<AGROW>
                    end
                end
                active_jobs = active_jobs - n_job_done;
                this.log.addMarkedMessage(sprintf('%d jobs completed', numel(completed_job)));
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
            this.sendMsg(this.MSG_KILLALL, 'As the mad king said: Kill them all!!!');
            pause(1);
            this.deleteMsg(Go_Slave.MSG_DIE, true);
            this.deleteMsg('*');
            delete(fullfile(this.COM_DIR,'*.mat'));
        end
        
        function resurgit(this, n_workers, verbosity)
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
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function initHandles(this)
            this.log = Core.getLogger();
        end
    end
end
