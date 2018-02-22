%   CLASS Remote_Resource_Manager
% =========================================================================
%
% DESCRIPTION
%
% EXAMPLE
%
% FOR A LIST OF CONSTANTs and METHODS use doc Remote_Resource_Manager

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 1 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
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

classdef Remote_Resource_Manager < Ini_Manager
    properties (Constant, Access = private)
        DEFAULT_RESOURCE_FILE = '../data/goGPSconfig/remote_resource.ini';
        DEFAULT_RESOURCE_TXT = ['to be filled'];
    end
    
    properties (Access = private)
    end
    
    properties (Access = private)
        log
    end
    
    methods
        function this = Remote_Resource_Manager(file_name)
            if (nargin == 0)
                file_name = Remote_Resource_Manager.DEFAULT_RESOURCE_FILE;
                if exist(file_name, 'file') ~= 2
                    %Remote_Resource_Manager.writeDefault(); %deafult file not
                    %stable enough
                end
            end
            
            this = this@Ini_Manager(file_name);
            this.readFile();
            this.log = Logger.getInstance();
        end
        
        function [ip, port] = getServerIp(this, name)
            % return the ip of a server given the server name
            ip = [];
            port = [];
            ip_port = this.getData('SERVER', name);
            ip = ip_port{1};
            port = ip_port{2};
        end
        
        function f_struct = getFileLoc(this, file_name)
            % return the remote path of the file
            f_struct.filename = this.getData(['f_' file_name],'filename');
            f_struct.const = this.getData(['f_' file_name],'sys');
            locations = this.getData(['f_' file_name],'location');
            if ~iscell(locations)
                locations = {locations};
            end
            f_struct.loc_number = length(locations);
            for i = 1 : f_struct.loc_number
                f_struct.(['loc' sprintf('%03d',i)]) = this.getData('LOCATION',locations{i});
            end
            
            
        end
        
        
        function [file_structure, latency] = getFileStr(this, center_name, resource_name)
            str = this.getData(['c_' center_name], resource_name);
            if isempty(str)
                this.log.addWarning(sprintf('No resource %s for center %s',resource_name, center_name))
                file_structure = [];
                latency = [];
            else
                file_structure = this.parseLogicTree(str);
                latency = this.getData(['c_' center_name], [resource_name '_latency']);
            end
        end
    end
    
    methods ( Access = private)
        function file_structure = parseLogicTree(this, str)
            [status, list] = this.findElements(str);
            file_structure = [];
            if status == 0
                file_structure  = {strtrim(str), false};
                return
            else
                if status == 1
                    cond = 'or';
                else
                    cond = 'and';
                end
                for i = 1 : length(list)
                    file_structure.(cond).(['f' num2str(i)]) = this.parseLogicTree(list{i});
                end
            end
            
        end
        
        function [status, list] = findElements(this, str)
            % OUTPUT:
            % status: -1 and 0 nothing 1 or
            % list: list of string parts
            [matches] = regexp(str, '\&|\|', 'match');
            if isempty(matches)
                status = 0;
                list = {};
                return
            else
                open = 0; %number if open pharentesis
                status = '';
                index = [];
                for i = 1:length(str)
                    if open == 0
                        if str(i) == '&' | str(i) == '|'
                            if isempty(status)
                                status = str(i);
                                index = [index; i];
                            else
                                if status ~= str(i)
                                    this.log.addWarning('| and & can not exist at the same level, check parenthesis')
                                    status = 0;
                                    return
                                else
                                    index = [index; i];
                                end
                            end
                        end
                    end
                    if str(i) == '('
                        open = open + 1;
                    end
                    if str(i) == ')'
                        open = open - 1;
                    end
                end
                if status == '|'
                    status = 1;
                elseif status == '&'
                    status = -1;
                end
                list = {};
                for i = 1 : length(index)
                    if i == 1
                        list{end + 1} = this.removeTrailingPar(str(1:index(i)-1));
                    else
                        list{end + 1} = this.removeTrailingPar(str(index(i-1)+1: index(i)-1));
                    end
                end
                list{end + 1} = this.removeTrailingPar(str(index(end)+1 : end));
                
            end
        end
        
        function str = removeTrailingPar(this, str)
            for i =1 :length(str)
                if str(i)~=' '
                    if str(i)=='('
                        str(1:i) = [];
                    end
                    break
                end
            end
            for i =length(str) : -1:1
                if str(i)~=' '
                    if str(i)==')'
                        str(i:end) = [];
                    end
                    break
                end
            end
        end
    end
    methods (Static)
        function writeDefault(this)
            fid = fopen(Remote_Resource_Manager.DEFAULT_RESOURCE_FILE,'w+');
            fprintf(fid, Remote_Resource_Manager.DEFAULT_RESOURCE_TXT);
            fclose(fid);
        end
    end
end
