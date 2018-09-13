%   CLASS Go_Slave
% =========================================================================
%
% DESCRIPTION
%   Slave controller for parallel goGPS computation
%
% EXAMPLE
%   gom = Go_Slave
%
% FOR A LIST OF CONSTANTs and METHODS use doc Command_Interpreter


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 4 - nightly
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

classdef Go_Slave < Com_Interface        
    %% PROPERTIES CONSTANTS
    % ==================================================================================================================================================
    properties (Constant, GetAccess = public)        
        SLAVE_WAIT_PREFIX = 'SLAVE_'
        SLAVE_READY_PREFIX = 'WORKER_'
        
        MSG_BORN = 'HELO_'
        MSG_DIE = 'ADIOS_'
        MSG_ACK = 'ACK_'
        MSG_JOBREADY = 'DONE_'
        
        EXIT = -1;
        RESTART = 0;
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
                    if nargin > 1 && ~isempty(com_dir)
                        this.initComDir(com_dir);
                    end
                end
            end
        end
    end
    
    methods (Access = public)
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
                rng shuffle;
                this.id = [this.SLAVE_WAIT_PREFIX sprintf('%06d', randi(1000000))];
            else
                Go_Slave.getInstance([], true)
            end
        end
        
        function msg = waitMsg(this, msg, remove_msg, check_revive)
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
            is_active = false;
            reset_count = 1;
            while ~is_active
                slave_list = dir(fullfile(this.getComDir, [Go_Master.MSG_KILLALL Go_Master.ID]));
                if ~isempty(slave_list)
                   is_active = true;
                   msg = this.EXIT; % This means exit and die
                else
                    slave_list = dir(fullfile(this.getComDir, [Go_Master.MSG_RESTART Go_Master.ID]));
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
                            reset_count = this.pause(0.1, reset_count);
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
            
            % Handshake
            stay_alive = true;
            while stay_alive
                this.revive();
                this.sendMsg(this.MSG_BORN, sprintf('Helo! My name is "%s"', this.id));
                this.waitMsg([this.id, '_' Go_Master.MSG_ASKACK Go_Master.ID], true); % WAIT ACK MESSAGE
                this.deleteMsg();
                this.sendMsg(this.MSG_ACK, sprintf('I''m ready to work!'));
                msg = this.waitMsg([this.id, '_' Go_Master.MSG_ASKWORK '*' Go_Master.ID], true); % WAIT WORK MESSAGE
                if ~(isnumeric(msg))
                    this.id = regexp(msg, [Go_Slave.SLAVE_READY_PREFIX '[0-9]*'], 'match', 'once');
                    
                    % Creating worker
                    core = Core.getInstance(); % Init Core
                    this.waitMsg([Go_Master.BRD_STATE Go_Master.ID], false, true); % WAIT WORK MESSAGE
                    tmp = load(fullfile(this.getComDir, 'state.mat'), 'geoid', 'state', 'cur_session', 'rin_list', 'met_list');
                    core.state = tmp.state; % load the state
                    core.gc.cur_settings = tmp.state; % load the state
                    core.gc.initGeoid(tmp.geoid); % load the geoid
                    core.cur_session = tmp.cur_session; % load the current session number
                    core.rin_list = tmp.rin_list; % load the rinex list of files
                    core.met_list = tmp.met_list; % load the meteorological list of files
                    this.log.addMarkedMessage('State updated');
                    clear tmp;
                    this.waitMsg([Go_Master.BRD_SKY Go_Master.ID]); % WAIT WORK MESSAGE
                    tmp = load(fullfile(this.getComDir, 'sky.mat'), 'sky', 'atmo', 'mn');
                    core.sky  = tmp.sky;  % load the state
                    core.atmo = tmp.atmo; % load the atmosphere
                    core.mn = tmp.mn; % load the meteorological network
                    clear tmp;
                    this.sendMsg(this.MSG_ACK, sprintf('Sky loaded'));
                    this.sendMsg(this.MSG_BORN, sprintf('Helo! My new name is "%s", gimme work', this.id));
                    
                    % Waiting work
                    this.waitMsg([Go_Master.BRD_CMD Go_Master.ID], false); % WAIT ACK MESSAGE
                    
                    active_ps = true;
                    while active_ps
                        msg = this.waitMsg([this.id '_' Go_Master.MSG_DO '*' Go_Master.ID], true, true); % WAIT ACK MESSAGE
                        if isnumeric(msg)
                            active_ps = false;
                        else
                            cmd_file = load(fullfile(this.getComDir, 'cmd_list.mat'));
                            rec_id = str2double(regexp(msg, '[0-9]*(?=_MASTER)', 'match', 'once'));
                            
                            % prepare receiver
                            state = core.getState();
                            clear rec
                            for r = 1 : rec_id
                                rec(r) = GNSS_Station(state.getConstellationCollector(), state.getDynMode() == 0); %#ok<AGROW>
                            end
                            core.rec = rec;
                            
                            for c = 1 : numel(cmd_file.cmd_list)
                                cmd_file.cmd_list{c} = strrep(cmd_file.cmd_list{c},'$', num2str(rec_id));
                            end
                            core.exec(cmd_file.cmd_list);
                            
                            % Export work
                            rec = core.rec(rec_id);
                            rec.out = []; % do not want to save out
                            save(fullfile(this.getComDir, sprintf('job%04d_%s.mat', rec_id, this.id)), 'rec');
                            clear rec;
                            core.rec = []; % empty space
                            this.sendMsg(this.MSG_JOBREADY, sprintf('Work done!'));
                        end
                    end
                    clear cmd_file;
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
