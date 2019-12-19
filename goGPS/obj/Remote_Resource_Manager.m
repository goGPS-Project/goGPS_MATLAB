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
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
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

classdef Remote_Resource_Manager < Ini_Manager
    properties (Constant, Access = private)
        DEFAULT_RESOURCE_FILE = 'remote_resource.ini';
        DEFAULT_RESOURCE_TXT = 'to be filled';
    end
    
    properties (Access = private)
        local_storage = ''
    end
    
    properties (Access = private)
        log
    end
    
    % =========================================================================
    %  INIT
    % =========================================================================
    methods (Access = public)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function this = Remote_Resource_Manager(file_name)
            if (nargin == 0)
                file_name = Remote_Resource_Manager.DEFAULT_RESOURCE_FILE;
                if exist(file_name, 'file') ~= 2
                    %Remote_Resource_Manager.writeDefault(); %deafult file not
                    %stable enough
                end
            end
            this = this@Ini_Manager(file_name);
            
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
            
            this.readFile();
            this.log = Core.getLogger();
        end
                
    end

    % =========================================================================
    %  SINGLETON GETTERS
    % =========================================================================
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function this = getInstance(resources_file, force_clean)
            % Get the persistent instance of the class
            persistent unique_instance_resource_manager__
            ini_is_present = false;
            switch nargin
                case 0
                    force_clean = false;
                case 1
                    if ischar(resources_file)
                        ini_is_present = true;
                        force_clean = false;
                    else
                        force_clean = resources_file;
                    end
                case 2
                    ini_is_present = true;
            end

            if force_clean
                log = Core.getLogger();
                clear unique_instance_resource_manager__;
                unique_instance_resource_manager__ = [];
            end

            if isempty(unique_instance_resource_manager__)
                if ini_is_present
                    this = Remote_Resource_Manager(resources_file);
                else
                    this = Remote_Resource_Manager();
                end
                unique_instance_resource_manager__ = this;
            else
                if ini_is_present
                    this = Remote_Resource_Manager(resources_file);
                    unique_instance_resource_manager__ = this;
                else
                    this = unique_instance_resource_manager__;
                end
            end
        end
    end
    
    methods        
        function [ip, port] = getServerIp(this, name)
            % Return the ip of a server given the server name
            %
            % SYNTAX:
            %   [ip, port] = this.getServerIp(name)
            ip = [];
            port = [];
            ip_port = this.getData('SERVER', name);
            ip = ip_port{1};
            port = ip_port{2};
        end
        
        function f_struct = getFileLoc(this, file_name)
            % Return the remote path of the file
            %
            % SYNTAX:
            %   f_struct = this.getFileLoc(file_name)
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
            % Get the logical file structure for the desidered center and
            % resource the latecncy of the resource
            %
            % SYNTAX
            %   [file_structure, latency] = this.getFileStr(center_name, resource_name)
            % 
            % OUTPUT:
            % file_strcuture
            %           the structure is a tree and can cointains fields
            %           named in 3 way 'fn' (where n is progressive number at
            %           the same level of the structure) , 'and' and 'or'.
            %           Or means that all sub field of the structure has to
            %           be found or means at least one.
            %           Leaves of the tree are cell containing the file code to 
            %           be found in the remote resource ini file and a boolean to
            %           tell if the file has been found or not
            %
            %           example: fs.and.f1
            %                          .f2.or.f1
            %                                .f2
            %                          .f3
            %                     
            %                    f1 = {'cnes_erp' , 0}
            % latency
            %           [h1 h2] h1 -> hours before which we now the
            %                         resource is not there
            %                   h2 -> hours after which we are sure we will
            %                         found the resource
            %           
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
    
    % Specific file query methods
    methods 
        function [center, center_ss] = getCenterList(this)
            % Get the list of available centers and the supported constellations
            %
            % SYNTAX
            %   [center, center_ss] = this.getCenterList()
            tmp = this.getData('CENTER', 'available');
            n_center = numel(tmp);
            center = cell(n_center, 1);
            center_ss = cell(n_center, 1);
            for c = 1 : n_center
                center{c} = regexp(tmp{c}, '(?<=@).*', 'match', 'once');
                if isempty(center{c})
                    center{c} = tmp{c};
                end
                center_ss{c} = regexp(tmp{c}, '.*(?=@)', 'match', 'once');;
            end            
        end
        
        function [center, center_ss] = getCenterListExtended(this)
            % Get the list of available centers and the supported constellations
            %
            % SYNTAX
            %   [center, center_ss] = this.getCenterList()
            tmp = this.getData('CENTER', 'available');
            n_center = numel(tmp);
            center = cell(n_center, 1);
            center_ss = cell(n_center, 1);
            for c = 1 : n_center
                center{c} = regexp(tmp{c}, '(?<=@).*', 'match', 'once');
                if isempty(center{c})
                    center{c} = tmp{c};
                else
                    center{c} = [upper(center{c}) ' - ' this.getData(['c_' center{c}], 'description')];
                end
                center_ss{c} = regexp(tmp{c}, '.*(?=@)', 'match', 'once');;
            end            
        end        
        
        function [descr] = centerToString(this, center)
            % Get the center description
            %
            % SYNTAX
            %   descr = this.centerToString(center)
            descr = sprintf('Center: "%s" - %s\n\nAvailable resources:\n', center, this.getData(['c_' center], 'description'));
            
            % get resource_list
            key = this.getKeys(['c_' center]);
            for k = 1 : numel(key)
                if isempty(regexp(key{k}, '(description)|(.*_latency)', 'once'))
                    descr = sprintf('%s%s', descr, this.resourceTreeToString(center, key{k}));
                end                
            end
        end
        
        function [flag_frub] = getOrbitType(this, center)
            % Get the orbit type availability
            %
            % SYNTAX
            %   [flag_frub] = this.getOrbitType(center)            
            flag_frub(1) = ~isempty(this.getData(['c_' center], 'final'));
            flag_frub(2) = ~isempty(this.getData(['c_' center], 'rapid'));
            flag_frub(3) = ~isempty(this.getData(['c_' center], 'ultra'));
            flag_frub(4) = ~isempty(this.getData(['c_' center], 'broadcast'));
        end
        
        function [flag_fp1p2b] = getIonoType(this, center)
            % Get the iono type availability
            %
            % SYNTAX
            %   [flag_frub] = this.getIonoType(center)            
            flag_fp1p2b(1) = ~isempty(this.getData(['c_' center], 'iono_final'));
            flag_fp1p2b(2) = ~isempty(this.getData(['c_' center], 'iono_predicted1'));
            flag_fp1p2b(3) = ~isempty(this.getData(['c_' center], 'iono_predicted2'));
            flag_fp1p2b(4) = ~isempty(this.getData(['c_' center], 'iono_broadcast'));
            
            if ~any(flag_fp1p2b)
                % Switch to default center
                flag_fp1p2b(1) = ~isempty(this.getData(['c_default'], 'iono_final'));
                flag_fp1p2b(2) = ~isempty(this.getData(['c_default'], 'iono_predicted1'));
                flag_fp1p2b(3) = ~isempty(this.getData(['c_default'], 'iono_predicted2'));
                flag_fp1p2b(4) = ~isempty(this.getData(['c_default'], 'iono_broadcast'));
            end
        end
        
        function [tree_str] = resourceTreeToString(this, center, resource_type)
            % Parse a resource tree
            %
            % SYNTAX
            %   [tree_str] = this.resourceTreeToString(center, resource_type)
            tmp = this.getData(['c_' center], resource_type);
            tree = this.parseLogicTree(tmp);
            if iscell(tree)
                clear tmp
                tmp.f1 = tree;
                tree = tmp;
            end
            tree_str = sprintf(' - %s\n', resource_type);
            if ~isempty(tree)
                tree_str = sprintf('%s%s', tree_str, this.treeToString(tree, 4));
            end
        end
        
        function remote_path = resourceToString(this, resource)
            if strcmp(resource, 'null')
                remote_path = 'none';
            else
                try
                    file_name = this.getData(['f_', resource], 'filename');
                    location = this.getData(['f_', resource], 'location');
                    if iscell(location)
                        location = location{1};
                    end
                    if isempty(location)
                        location = '';
                        server = cell(2,1);
                    else
                        location = this.getData('LOCATION', location);
                        % remove the server
                        server = regexp(location, '(?<=\?\{)[a-z\_A-Z]*(?=\})', 'match');
                        if isempty(server)
                            server = cell(2,1);
                        else
                            location = regexp(location, '(?<=(\?\{[a-z_A-Z]*\})).*', 'match');
                        end
                        server = this.getData('SERVER', server{1});
                    end
                    % The protocol could be improved woth more values
                    if iscell(server) && numel(server) == 2
                        switch server{2}
                            case '21'
                                protocol = 'ftp://';
                            otherwise
                                protocol = 'http://';
                        end
                    else
                        protocol = '';
                    end
                    remote_path = [protocol server{1} ':' server{2} location{1} file_name];
                catch
                    remote_path = sprintf('corrupted resource "%s" -> check resource file', resource);
                end
            end
        end
        
        function tree_str = treeToString(this, tree, indentation_lev)
            % Parse a resource tree
            %
            % SYNTAX
            %   [tree_str] = this.treeToString(tree, indentation_lev)
            if nargin == 1
                indentation_lev = 0;
            end
            str_padding = char(32*ones(1, indentation_lev + 1));
            tree_str = '';
            tree_names = fieldnames(tree);
            for b = 1 : numel(tree_names)
                if this.isNode(tree_names(b))
                    % is a node
                    tree_str = sprintf('%s%s|- %s\n%s', tree_str, str_padding, upper(tree_names{b}), this.treeToString(tree.(tree_names{b}), indentation_lev + 4));
                else
                    % is a leaf
                    if iscell(tree.(tree_names{b}))
                        % tree_str = sprintf('%s%s|- %s\n', tree_str, str_padding, tree.(tree_names{b}){1});
                        tree_str = sprintf('%s%s|- %s\n', tree_str, str_padding, this.resourceToString(tree.(tree_names{b}){1}));
                    elseif isstruct(tree.(tree_names{b}))
                        tree_str = sprintf('%s%s', tree_str, this.treeToString(tree.(tree_names{b}), indentation_lev));                        
                    else
                        % NON VALID
                    end
                end
            end
        end
    end
    
    methods ( Access = public)
        function file_structure = parseLogicTree(this, str)
            % Description parse the logic structure found in
            % remote_resource.ini to get the file structure descripbed in
            % this.getFileStr
            %
            % SYNTAX:
            %   file_structure = this.parseLogicTree(this)
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
            % Return the element of the string separated by | and &
            %
            % SYNTAX:
            %       [status, list] = this.findElements(str)
            %
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
            % Remove trailing parenthesis
            %
            % SYNTAX:
            %   str = this.removeTrailingPar(str)
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
            % Write the deafut remote resource ini file if it is not found 
            %
            % SYNTAX:
            %       Remote_Resource_Manager.writeDefault()
            fid = fopen(Remote_Resource_Manager.DEFAULT_RESOURCE_FILE,'w+');
            fprintf(fid, Remote_Resource_Manager.DEFAULT_RESOURCE_TXT);
            fclose(fid);
        end
        
        function is_node = isNode(tree_el)
            % Local usage return if a node of node_tree is an "and" or an "or"
            if numel(tree_el) == 1 && iscell(tree_el) && ~isempty(regexp(tree_el{1}, '(and)*|(or)*'))
                is_node = true;
            else
                is_node = false;
            end
        end
        
    end
end
