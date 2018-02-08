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
    end
    
    properties (Access = private)
        servers
        file_resources
        remote_sources
        computational_centers
    end
    
    properties (Access = private)
        log        
    end
    
    methods
        function this = Remote_Resource_Manager(file_name)
            if (nargin == 0)
                file_name = Remote_Resource_Manager.DEFAULT_RESOURCE_FILE;
            end
            
            this = this@Ini_Manager(file_name);
            this.readFile();
            for i = 1 : length(this.section)
                if strcmp(this.section{i}.name, 'SERVER')
                    this.servers = this.section{i}.key;
                end
                if strcmp(this.section{i}.name(1:2), 'f_')
                    this.file_resources{end+1} = this.section{i};
                    this.file_resources{end}.name(1:2) = [];
                end
                if strcmp(this.section{i}.name(1:2), 'r_')
                    this.remote_sources{end+1} = this.section{i};
                    this.remote_sources{end}.name(1:2) = [];
                end
                if strcmp(this.section{i}.name(1:2), 'c_')
                    this.computational_centers{end+1} = this.section{i};
                    this.computational_centers{end}.name(1:2) = [];
                end
            end
            
            this.log = Logger.getInstance();
        end
        
        function [ip, port] = getServerIp(this, name)
            % return the ip of a server given the server name
            ip = [];
            port = [];
            for i = 1 : length(this.servers)
                if strcmp(this.servers{i}.name, name)
                    ip = this.servers{i}.data{1};
                    port = this.servers{i}.data{2};
                end
            end
        end
        
        function f_struct = getFileLoc(this, file_name, sys_c)
            % return the ip of a server given the server name
            for i = 1 : length(this.file_resources)
                name_part = strsplit(this.file_resources{i}.name,'@');
                name = name_part{1};
                if length(name_part) > 1
                    const = name_part{2};
                else
                    const = 'GRECJIS';
                end
                cond_const = true;
                if nargin > 2
                    cond_const = ~isempty(strfind(const, sys_c));
                end
                if strcmp(name, file_name) && cond_const
                    f_struct = struct();
                    f_struct.name = file_name;
                    f_struct.const = const;
                    for j = 1 : length(this.file_resources{i}.key)
                        name_k = this.file_resources{i}.key{j}.name;
                        f_struct.(name_k) = this.file_resources{i}.key{j}.data;
                    end
                    
                end
            end
        end
        
        function center_code = getCenterCode(this, center_name, resource_name, sys_c)
            % return the center code given a resource name and desired
            % constelltion
            for i = 1 :length(this.computational_centers)
                if strcmp(this.computational_centers{i}.name, center_name)
                    for j = 1 : length(this.computational_centers{i}.key)
                        resource = this.computational_centers{i}.key{j};
                        if strcmp(resource.name, resource_name)
                            idx = [];
                            if nargin > 3
                                valid = [];
                                if ~iscell(resource.data)
                                    resource.data = {resource.data};
                                end
                                for k = 1:length(resource.data)
                                    center_code_part = strsplit(resource.data{k},'@');
                                    if length(center_code_part) > 1
                                        consts = center_code_part{1};
                                        found = true;
                                        for l = 1 : length(sys_c)
                                            found = found && ~isempty(strfind(consts, sys_c(l)));
                                        end
                                        if found
                                            valid = [valid ; [k, (length(consts) -length(sys_c))]];
                                        end
                                    else
                                        valid = [valid; [k,Inf]];
                                    end
                                    if ~isempty(valid)
                                        idx = valid(valid(:,2) == min(valid(:, 2)), 1); %select the center_code that has all contellations ad give priority to the ones that has the minum number of other constellations
                                    end
                                end
                            else
                                idx = 1;
                            end
                            if isempty(idx)
                                this.log.addWarning('No vaild center code found for the desidered combination')
                            end
                            center_code_part = strsplit(resource.data{idx},'@');
                            if length(center_code_part) > 1
                                center_code = center_code_part{2};
                            else
                                center_code = center_code_part{1};
                            end
                            
                        end
                    end
                end
            end
        end
        
        function file_structure = getFileStr(this, resource_name)
            for i = 1 : length(this.remote_sources)
                if strcmp(this.remote_sources{i}.name, resource_name)
                    string = this.remote_sources{i}.key{1}.data;
                    file_structure = this.parseLogicTree(string);
                end
            end
        end        
    end
    
    methods ( Access = private)
        function file_structure = parseLogicTree(this, str)
            [status, list] = this.findElements(str);
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
end
