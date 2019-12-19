%   CLASS Go_Slave
% =========================================================================
%
% DESCRIPTION
%   Slave controller for parallel goGPS computation
%
% EXAMPLE
%   gos = Go_Slave
%
% FOR A LIST OF CONSTANTs and METHODS use doc Go_Slave


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
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

classdef Go_Slave < Com_Interface
    %% PROPERTIES CONSTANTS
    % ==================================================================================================================================================
    properties (Constant, GetAccess = public)
        SLAVE_WAIT_PREFIX = 'SLAVE_'
        SLAVE_SESSION_PREFIX = 'S_'
        SLAVE_TARGET_PREFIX = 'T_'
        SLAVE_READY_PREFIX = 'WORKER_'
                
        MSG_BORN = 'HELO_'
        MSG_DIE = 'ADIOS_'
        MSG_ACK = 'ACK_'
        MSG_JOBREADY = 'DONE_'
        
        EXIT = -1;
        RESTART = 0;
    end
    
    properties (SetAccess = private, GetAccess = private)
        rnd_id = 0;     % personal serial number of the slave
        slave_type = 1; % 0 for session worker
                        % 1 for target worker
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = private)
        % Concrete implementation.  See Singleton superclass.
        function this = Go_Slave(com_dir)
            % Core object creator
            this.id = '';
            this.initComDir(com_dir);
            this.init();
            this.log.addMarkedMessage('Creating goGPS Slave');
        end
    end
    
    %% METHOD DESTRUCTOR
    % ==================================================================================================================================================
    methods (Access = private)        
        function delete(this)
            % delete
        end
    end
    
    %% METHOD INTERFACE
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = getInstance(com_dir, destroy_this)
            % Get the persistent instance of the class
            %
            % SYNTAX
            %   gos = Go_Slave.getInstance(com_dir, destroy_this)
            persistent unique_instance_gos__
            if nargin >= 2 && ~isempty(destroy_this) && destroy_this
                if ~isempty(unique_instance_gos__)
                    unique_instance_gos__.delete();
                    clear unique_instance_gos__
                end
            else
                if isempty(unique_instance_gos__)
                    if nargin < 1 || isempty(com_dir)
                        com_dir = fullfile(pwd, 'com');
                    end
                    this = Go_Slave(com_dir);
                    unique_instance_gos__ = this;
                else
                    this = unique_instance_gos__;
                    this.init();
                    if nargin > 0 && ~isempty(com_dir)
                        this.initComDir(com_dir);
                    end
                end
            end
        end        
    end
    
    methods (Access = public)
        function is_trg_worker = isTargetWorker(this)
            is_trg_worker = this.slave_type == 1;
        end
        
        function is_sss_worker = isSessionWorker(this)
            is_sss_worker = this.slave_type == 0;
        end
        
        function revive(this)
            % Kill the go Slave but reborn with new id
            %
            % SYNTAX:
            %   this.revive();
            this.die(true);
        end
        
        function die(this, reborn)
            % Kill the go Slave
            %
            % SYNTAX:
            %   this.die();
            
            if ~isempty(this.id)
                this.log.addMarkedMessage('Closing goGPS Slave');
                this.deleteMsg();
                this.sendMsg(this.MSG_DIE, sprintf('I was "%s", goodbye cruel world!', this.id));
            end
            if nargin > 1 && ~isempty(reborn) && reborn
                if this.rnd_id == 0
                    rng('shuffle', 'simdTwister');
                    this.rnd_id = randi(1e15);
                end
                this.id = [this.SLAVE_WAIT_PREFIX sprintf('%015d', this.rnd_id)];
            else
                Go_Slave.getInstance([], true)
            end
        end
        
        function msg = checkMsg(this, msg, remove_msg, check_revive, wait_until_received)
            % Wait for a specific order
            %
            % SYNTAX:
            %   this.live();
            if nargin < 3 || isempty(remove_msg)
                remove_msg = false;
            end
            if nargin < 4 || isempty(check_revive)
                check_revive = false;
            end
            if nargin < 5 || isempty(wait_until_received)
                wait_until_received = true;
            end
            is_active = false;
            reset_count = 1;
            while ~is_active
                if msg(1) == 'S'
                    % Check that there is my helo message!
                    slave_list = dir(fullfile(this.getComDir, [this.MSG_BORN this.id]));
                    if isempty(slave_list)
                        % Resend helo message
                        this.sendMsg(this.MSG_BORN);
                    end
                end
                slave_list = dir(fullfile(this.getComDir, [Parallel_Manager.MSG_KILLALL Parallel_Manager.ID]));
                if ~isempty(slave_list)
                    is_active = true;
                    msg = this.EXIT; % This means exit and die
                else
                    slave_list = dir(fullfile(this.getComDir, [Parallel_Manager.MSG_RESTART Parallel_Manager.ID]));
                    if check_revive && ~isempty(slave_list)
                        is_active = true;
                        msg = this.RESTART; % This means exit and die
                    else
                        slave_list = dir(fullfile(this.getComDir, msg));
                        if ~isempty(slave_list)
                            if remove_msg
                                delete(fullfile(this.getComDir, slave_list(1).name))
                            end
                            is_active = true;
                            msg = slave_list(1).name;
                        else
                            if ~wait_until_received
                                % if I'm doing just a check on the existence of the message do not wait for it
                                is_active = true;
                                msg = [];
                            else
                                reset_count = this.pause(0.1, reset_count);
                            end
                        end
                    end
                end
            end
            this.pause(0, -1);
        end
        
        function live(this)
            % Start the life of the slave
            %
            % SYNTAX:
            %   this.live();        
            
            % No color mode for slaves
            this.log.setColorMode(0);
            
            % Handshake
            stay_alive = true;
            while stay_alive
                this.revive();
                this.sendMsg(this.MSG_BORN, sprintf('Helo! My name is "%s"', this.id));
                this.checkMsg([this.id, '_' Parallel_Manager.MSG_ASKACK Parallel_Manager.ID], true, false); % WAIT ACK MESSAGE
                this.sendMsg(this.MSG_ACK, sprintf('I''m ready to work!'));
                msg = this.checkMsg([this.id, '_' Parallel_Manager.MSG_ASKWORK '*' Parallel_Manager.ID], true, true); % WAIT WORK MESSAGE
                this.deleteMsg();
                if ~(isnumeric(msg))
                    
                    % Creating worker
                    core = Core.getCurrentCore(); % Init Core
                    core.initPreloadSession();
                    core.clearSingletons();
                    core.initSimpleHandlers();
                    core.initLocalPath();
                    this.checkMsg([Parallel_Manager.BRD_STATE Parallel_Manager.ID], false, false); % WAIT WORK MESSAGE
                    tmp = load(fullfile(this.getComDir, 'state.mat'), 'geoid', 'state', 'atx', 'cur_session', 'rin_list', 'met_list', 'slave_type');
                    this.slave_type = tmp.slave_type;
                    core.state = tmp.state; % load the state
                    core.atx = tmp.atx;   % load the antenna manager
                    core.setCurrentSettings(tmp.state); % load the state
                    core.initGeoid(tmp.geoid); % load the geoid
                    core.state.setCurSession(tmp.cur_session); % load the current session number
                    core.rin_list = tmp.rin_list; % load the rinex list of files
                    core.met_list = tmp.met_list; % load the meteorological list of files
                    
                    % No colormode for slaves (faster execution)
                    core.log.setColorMode(0);
                    Logger.getInstance.setColorMode(0);
                    
                    this.id = regexp(msg, [Go_Slave.SLAVE_READY_PREFIX iif(this.slave_type == 1, this.SLAVE_TARGET_PREFIX, this.SLAVE_SESSION_PREFIX) '[0-9]*'], 'match', 'once');
                    
                    this.log.addMarkedMessage('State updated');
                    clear tmp;
                    if this.isTargetWorker()
                        this.checkMsg([Parallel_Manager.BRD_SKY Parallel_Manager.ID], false, true); % WAIT WORK MESSAGE
                        tmp = load(fullfile(this.getComDir, 'sky.mat'), 'sky', 'atmo', 'mn');
                        core.sky  = tmp.sky;  % load the state
                        core.atmo = tmp.atmo; % load the atmosphere
                        core.mn = tmp.mn; % load the meteorological network
                        clear tmp;                        
                        this.log.addMarkedMessage('Sky updated');
                    else
                        core.atmo = Atmosphere();
                    end
                    
                    % Check for receiver to load
                    msg = this.checkMsg([Parallel_Manager.BRD_REC Parallel_Manager.ID], false, true, false); % CHECK REC PASSING MESSAGE
                    rec_pass = [];
                    if ~isempty(msg) && ~isnumeric(msg)
                        % I received a rec_list to load
                        rec_pass = load(fullfile(this.getComDir, 'rec_list.mat'), 'rec_work', 'rec_num');
                        this.log.addMarkedMessage('Passed receiver have been read');
                    end
                    this.sendMsg(this.MSG_ACK, sprintf('Everything loaded'));
                    this.sendMsg(this.MSG_BORN, sprintf('Helo! My new name is "%s", gimme work', this.id));
                    
                    % Waiting work
                    this.checkMsg([Parallel_Manager.BRD_CMD Parallel_Manager.ID], false, true); % WAIT ACK MESSAGE
                    
                    active_ps = true;
                    while active_ps
                        try
                            msg = this.checkMsg([this.id '_' Parallel_Manager.MSG_DO '*' Parallel_Manager.ID], true, true); % WAIT ACK MESSAGE
                            if isnumeric(msg)
                                active_ps = false;
                            else
                                cmd_file = load(fullfile(this.getComDir, 'cmd_list.mat'));
                                % get request id
                                req_id = str2double(regexp(msg, '[0-9]*(?=_MASTER)', 'match', 'once'));
                                
                                % prepare receivers
                                state = core.getState();
                                clear rec
                                n_rec = core.state.getRecCount;
                                for r = 1 : n_rec
                                    rec(r) = GNSS_Station(state.getDynMode() == 0); %#ok<AGROW>
                                end
                                core.rec = rec;
                                if ~isempty(rec_pass)
                                    for r = 1 : numel(rec_pass.rec_num)
                                        core.rec(rec_pass.rec_num(r)).work = rec_pass.rec_work(r);
                                        core.rec(rec_pass.rec_num(r)).work.parent = core.rec(rec_pass.rec_num(r));
                                        core.rec(rec_pass.rec_num(r)).work.initHandles();
                                    end
                                end
                                
                                % Substitute key "$" into command list with the one from PAR target
                                if this.isTargetWorker()
                                    for c = 1 : numel(cmd_file.cmd_list)
                                        cmd_file.cmd_list{c} = strrep(cmd_file.cmd_list{c},'$', num2str(req_id));
                                    end
                                end
                                
                                if this.isSessionWorker()
                                    core.state.setCurSession(req_id);
                                    core.prepareSession(req_id);
                                    %cmd_file.cmd_list = [{sprintf('FOR S%d', req_id);} cmd_file.cmd_list(:) {'ENDPAR'}];
                                end
                                
                                core.exec(cmd_file.cmd_list);
                                
                                if this.isTargetWorker()
                                    % Save the output ar rec
                                    % store command list in the rec as it was executed
                                    rec(req_id).work.state.cmd_list = cmd_file.cmd_list;
                                    
                                    % Export work
                                    rec = core.rec(req_id);
                                    rec.out = []; % do not want to save out
                                    save(fullfile(this.getComDir, sprintf('job%04d_%s.mat', req_id, this.id)), 'rec');
                                elseif this.isSessionWorker()
                                    % Export all the rec work spaces of the session
                                    rec = core.rec;
                                    % Make the receiver lighter to save
                                    for r = 1 : numel(rec)
                                        rec(r).out = []; % do not want to save out
                                        rec(r).clearHandles(); % do not want to save handles
                                        rec(r).work.clearHandles(); % do not want to save handles
                                    end
                                    atmo = Core.getAtmosphere;
                                    save(fullfile(this.getComDir, sprintf('job%04d_%s.mat', req_id, this.id)), 'rec', 'atmo');
                                end
                                pause(0.1); % be sure that the file is saved correctly
                                core.rec = []; % empty space
                                clear rec;
                                this.sendMsg(this.MSG_JOBREADY, sprintf('Work done!'));
                            end
                        catch ex
                            % Export work
                            try
                                rec = core.rec(req_id);
                            catch
                                % I'm going to create an empty rec if something
                                % goes wrong
                            end
                            try
                                if isempty(rec)
                                    rec = GNSS_Station(state.getDynMode() == 0);
                                end
                                rec.out = []; % do not want to save out
                                rec.work.flag_currupted = true;
                                save(fullfile(this.getComDir, sprintf('job%04d_%s.mat', req_id, this.id)), 'rec');
                                pause(0.1); % be sure that the file is saved correctly
                            catch
                                % try to send the receiver, if something goes bad,
                                % the master with deal with it
                            end
                            core.rec = []; % empty space
                            clear rec;
                            
                            % If something bad happen during work restart
                            this.sendMsg(this.MSG_JOBREADY, sprintf('Work done!'));
                            this.log.addError(sprintf('Something bad happened: %s\n', ex.message));
                        end
                    end
                    clear cmd_file rec_pass;
                end
                if isnumeric(msg)
                    switch msg
                        case this.EXIT
                            % It means that my services are no more needed
                            stay_alive = false;
                            this.log.addMarkedMessage('Thank you Master!\n Have been a pleasure living for you');
                        otherwise % RESTART
                            this.log.addMarkedMessage('Thank you Master!');
                    end
                    clear core state rec_id r msg
                end
            end
            
            this.die();
            exit
        end
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this)
            % Init the slave object and delete old messages of its previous life
            %
            % SYNTAX
            %   this.init();
            this.initLogger();
        end
    end
    
end
