%   CLASS Command_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the goGPS command
%
% EXAMPLE
%   settings = Command_Settings();
%
% FOR A LIST OF CONSTANTS and METHODS use doc Command_Settings

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti
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

classdef Command_Settings < Settings_Interface

    
    % Default values for each field - useful to restore corrupted fields
    properties (Constant, GetAccess = public)
        CMD_LIST = {};
    end

    properties (Constant, GetAccess = protected)
        CMD_SECTION = 'COMMANDS';
    end

    properties (SetAccess = public, GetAccess = public)
        cmd_list = Command_Settings.CMD_LIST;
    end

    % =========================================================================
    %  INIT
    % =========================================================================
    methods
        function this = Command_Settings()
            % Creator of Command_Settings
        end
    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods
        function importPlainCommands(this, txt)
            if isempty(txt)
                this.cmd_list = [];
            else
                log = Core.getLogger();
                ci = Core.getCommandInterpreter();
                n_cmd = size(txt, 1);
                [cmd, err_list] = ci.fastCheck(txt);
                this.cmd_list = cmd;
            end
        end
        
        function import(this, settings)
            % This function import Mode (only) settings from another setting object
            %
            % SYNTAX:
            %   this.import(settings);
            
            if isa(settings, 'Ini_Manager')
                cmd_keys = settings.getKeys(this.CMD_SECTION);
                this.cmd_list = {};
                % cleaning list
                l = 0;
                while l < numel(cmd_keys)
                    l = l + 1;
                    if (numel(cmd_keys{l}) < 3) || ~strcmpi(cmd_keys{l}(1:3), 'cmd')
                        Core.getLogger.addWarning(sprintf('%s command unrecognized\nit should start with "cmd" (e.g. cmd_001)', cmd_keys{l}));
                        cmd_keys(l) = [];
                        l = l - 1;
                    else
                        this.cmd_list{l} = settings.getData(this.CMD_SECTION, cmd_keys{l});
                        if iscell(this.cmd_list{l}) && (numel(this.cmd_list{l}) > 1)
                            this.cmd_list{l} = [this.cmd_list{l}{1} '"' this.cmd_list{l}{2:end} '"'];
                        end
                    end
                end
            else
                this.cmd_list = settings.cmd_list;
            end
            this.check();
        end

        function str = toString(this, str)
            if (nargin == 1)
                str = '';
            end
            str = this.cmdToString(str);
        end
        
        function str = cmdToString(this, str)
            % Display the command list
            %
            % SYNTAX:
            %   str = this.cmdToString(str);
            if (nargin == 1)
                str = '';
            end
            if ischar(str)
                str = [str '---- CMD LIST ------------------------------------------------------------' 10 10];
                str = [str sprintf('Execution list:\n')];
                cmd_list = this.cmd_list;
            else
                cmd_list = str;
                str = '';
            end
            for l = 1 : numel(cmd_list)
                str = [str sprintf(' %03d %s\n', l, strrep(cmd_list{l}, Command_Interpreter.SUB_KEY, ' '))];
            end
            %str = [str 10];
        end

        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the object
            %
            % SYNTAX:
            %   str_cell = this.export(str_cell);
            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniStringSection(this.CMD_SECTION, str_cell);
            str_cell = Ini_Manager.toIniStringComment('goGPS command list', str_cell);
            str_cell = Ini_Manager.toIniStringComment('NOTE: All the commands will be executed for each session', str_cell);
            cmd = Core.getCommandInterpreter();
            % To be moved in the manual in the future
            str_cell = Ini_Manager.toIniStringComment(cmd.getHelp, str_cell);
            for l = 1 : numel(this.cmd_list)
                str_cell = Ini_Manager.toIniString(sprintf('cmd_%03d', l), strrep(this.cmd_list{l}, Command_Interpreter.SUB_KEY, ' '), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end
        
        function str_cell = exportCmdListHelp(this, str_cell)
            % Conversion to string ini format of the help
            %
            % SYNTAX:
            %   str_cell = this.exportCmdListHelp(str_cell);
            
            if (nargin == 1)
                str_cell = {};
            end
            
            cmd = Core.getCommandInterpreter();
            % To be moved in the manual in the future
            str_cell = Ini_Manager.toIniStringComment(cmd.getHelp, str_cell);
        end
        
        function str_cell = exportCmdListExamples(this, str_cell)
            % Conversion to string ini format of the example list
            %
            % SYNTAX:
            %   str_cell = this.exportCmdListExamples(str_cell);
            
            if (nargin == 1)
                str_cell = {};
            end
            
            cmd = Core.getCommandInterpreter();
            % To be moved in the manual in the future
            str_cell = strrep(cmd.getExamples,'#','%');
        end
        
        function str_cell = exportCmdList(this, str_cell)
            % Conversion to string ini format of the command list
            %
            % SYNTAX:
            %   str_cell = this.export(str_cell);
            
            if (nargin == 1)
                str_cell = {};
            end
            
            cmd = Core.getCommandInterpreter();
            [cmd_list, ~, loop_lev] = cmd.fastCheck(this.cmd_list);
            % find how to indent commands:
            loop_lev = loop_lev - (diff([0 loop_lev]) > 0);
            this.cmd_list = cmd_list;
            % To be moved in the manual in the future            
            for l = 1 : numel(this.cmd_list)
                str_cell = [str_cell; {sprintf('%s%s', char(32 * ones(1,3 * loop_lev(l))), strtrim(strrep(strrep(this.cmd_list{l}, Command_Interpreter.SUB_KEY, ' '), '''', '"')))}];
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end
    end
   
    % =========================================================================
    %  CHECKER CALL
    % =========================================================================
    methods (Access = 'private')
        function check(this)
            if isempty(this.cmd_list)
            % Minimal check on cmd_list entries validity
            % Modify and (try to) correct cmd_list
            %
            % SYNTAX:
            %   this.check();
            %
                this.cmd_list = this.CMD_LIST;
                Core.getLogger.addWarning(sprintf('Command list seems to be empty\n'));
            else
                this.cmd_list = Core.getCommandInterpreter.fastCheck(this.cmd_list);
            end
        end
    end
    
    methods (Access = 'public')
        function cmd_list = getCommandList(this)
            cmd_list = this.cmd_list;
        end
    end

    methods (Access = 'public')
        function setCommandList(this, cmd_list)
            this.cmd_list = cmd_list;
        end
    end

    % =========================================================================
    %  TEST
    % =========================================================================
    methods (Static, Access = 'public')
        function test()
            % test the class
            %
            % SYNTAX:
            %   this.test();
            s = Command_Settings();
            s.testInterfaceRoutines();
        end
    end
end
