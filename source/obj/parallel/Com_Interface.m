%   CLASS Com_Interface
% =========================================================================
%
% DESCRIPTION
%   Abstract class with basic properties and methods that a setting class
%   must have
%
% COMMENTS
% Settings have been build with multiple inheritance
% A standard abstract interface have been created: Setting_Interface
% it add to each subclass the object log
% force the subclasses to implement three basic methods:
%  - import:       when a setting object of identical class, or that inherits from the same class is passed to this function, the relative parameters are copied in the calling object
%                  when an ini file is the input of the function, the object is updated with the settings contained into the file
%  - toString:     display the content of the object, in a human readable way, a goGPS user can "ask" for the value of the settings on screen
%  - export:       create a cell array of strings containing the settings in plain text ini format. The variable it's the raw data format of Ini_Manager
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 (GReD srl) Andrea Gatti, Giulio Tagliaferro, Eugenio Realini
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti
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

classdef Com_Interface < handle
    properties (SetAccess = protected, GetAccess = protected)
        log;        % Handler to the log object
        id;         % Id of the sender        

        nch_wait = 0;
    end

    properties
        COM_DIR = './com';
        working_dir = '';
    end    

    methods  (Abstract)
    end

    methods (Access = 'public')
        function com_dir = getComDir(this)
            % Get the comunication dir
            % 
            % OUTPUT
            %   com_dir = fullfile(pwd, 'com') if empty
            %
            % SYNTAX
            %   com_dir = this.getComDir()
            %
            if isempty(this.COM_DIR)
                com_dir = Core.getState.getComDir();
            else
                com_dir = this.COM_DIR;
            end
        end

        function setComDir(this, com_dir)
            % Set the comunication dir
            % (must be used on object creation)
            % 
            % SYNTAX
            %   this.initComDir(com_dir)
            %
            % EXAMPLE default)
            %   this.initComDir(fullfile(pwd, 'com'));
            
            this.COM_DIR = com_dir;
        end
        
        function initComDir(this, com_dir)
            % alias of setComDir
            % Set the comunication dir
            % (must be used on object creation)
            % 
            % SYNTAX
            %   this.initComDir(com_dir)
            %
            % EXAMPLE default)
            %   this.initComDir(fullfile(pwd, 'com'));
            
            this.setComDir(com_dir);
        end
        
        function initLogger(this)
            % Init the log object
            %
            % SYNTAX
            %   this.initLogger();
            this.log = Core.getLogger();
            this.log.setOutMode(1,[],0);
            this.log.setColorMode(0);
        end
        
        function sendMsg(this, msg, msg_feedback)
            % Send a message to the comunication interface
            %
            % SINTAX
            %   this.sendMsg(msg, msg_feedback)
            %
            % EXAMPLE
            %   this.sendMsg(this.idRep(this.MSG_BORN), sprintf('Helo! I''m "%s"', this.id))
            if ~exist(this.getComDir, 'file')
                mkdir(this.getComDir);
            end
            fid = fopen(fullfile(this.getComDir, this.working_dir, [msg, this.id]), 'w');
            if fid > 0
                fclose(fid);
                if nargin > 2 && ~isempty(msg_feedback)
                    this.log.addMessage(this.log.indent([' > ' msg_feedback]));
                end
            end
        end
        
        function is_master = isMaster(this)
            % tell if the caller is the master
            %
            % SINTAX
            %   this.isMaster()
            is_master = strcmp(this.id, Parallel_Manager.ID);
        end
        
        function deleteAllMsg(this)
            % Delete all the messages of the caller
            % (of everyone in case of Master calling)
            %
            % SINTAX
            %   this.deleteAllMsg()

            if this.isMaster
                this.deleteMsg('*', true);
            else
                this.deleteMsg('*');
            end
        end
        
        function deleteMsg(this, msg_type, all_ids)
            % Delete all the messages with a certain msg_type
            %
            % SINTAX
            %   this.deleteMsg(<msg_type>, <all_ids>)
            %
            % EXAMPLE
            %   this.deleteMsg()
            
            if exist(this.getComDir, 'file')
                if nargin < 2 || isempty(msg_type)
                    msg_type = '*';
                end
                if nargin >= 3 && ~isempty(all_ids) && all_ids
                    msg_type = [msg_type '*'];
                else
                    msg_type = [msg_type this.id];
                end
                msg_type = strrep(msg_type, '**', '*');
                warning off; % I don't care if the file does not exist
                delete(fullfile(this.getComDir, this.working_dir, msg_type));
                warning on;
            end            
        end

        function reset_count = pause(this, seconds, reset_count)
            % Display wait message
            %
            % SINTAX
            %   this.deleteMsg(<msg_type>, <all_ids>)
            
            if nargin > 2 && ~isempty(reset_count) && reset_count > 0
               this.nch_wait = 0;
            end
            
            if this.nch_wait == 20 || nargin > 2 && ~isempty(reset_count) && reset_count ~= 0
                fprintf(char(ones(1, this.nch_wait) * 8)); % delete last chars
                this.nch_wait = 0;
            end
            
            if nargin < 3 || (~isempty(reset_count) && reset_count >= 0)
                switch (this.nch_wait)
                    case 0, fprintf('Waiting'); this.nch_wait = 7;
                    otherwise, fprintf('.'); this.nch_wait = this.nch_wait + 1;
                end
                pause(seconds);
            end 
            reset_count = 0;            
        end
    end

    % =========================================================================
    %%  SETTERS IO
    % =========================================================================
    methods
        function value = setProperty(this, prop_name, value)
            this.(prop_name) = value;
        end
    end
    
    % =========================================================================
    %%  GETTERS IO
    % =========================================================================
    methods
        function value = getProperty(this, prop_name)
            props = fieldnames(this);
            if any(strcmp(props,prop_name))
                value = this.(prop_name);
            else
                value = [];
            end
        end
    end
    
end
