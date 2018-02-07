classdef Remote_Resource_Manager < Ini_Manager
    properties
        servers
        file_resources
        remote_sources
        computational_centers
    end
    properties (Access = private)
        log
        
end
    methods
        function this = Remote_Resource_Manager(filename)
            this = this@Ini_Manager(filename);
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
        function f_struct = getFileLoc(this, fname, sys_c)
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
                if strcmp(name, fname) && cond_const
                    f_struct = struct();
                    f_struct.name = fname;
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
        function file_structure = parseLogicTree(this, string)
            [status, list] = this.findElements(string);
            if status == 0
                file_structure  = {strtrim(string), false};
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
        function [status, list] = findElements(this, string)
            % OUTPUT:
            % status: -1 and 0 nothing 1 or
            % list: list of string parts
            [matches] = regexp(string, '\&|\|', 'match');
            if isempty(matches)
                status = 0;
                list = {};
                return 
            else
                open = 0; %number if open pharentesis
                status = '';
                index = [];
                for i = 1:length(string)
                    if open == 0
                        if string(i) == '&' | string(i) == '|'
                            if isempty(status)
                                status = string(i);
                                index = [index; i];
                            else
                                if status ~= string(i)
                                    this.log.addWarning('| and & can not exist at the same level, check parenthesis')
                                    status = 0;
                                    return
                                else
                                    index = [index; i];
                                end
                            end
                        end
                    end
                    if string(i) == '('
                        open = open + 1;
                    end
                    if string(i) == ')'
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
                    list{end + 1} = this.removeTrailingPar(string(1:index(i)-1));
                    else
                        list{end + 1} = this.removeTrailingPar(string(index(i-1)+1: index(i)-1));
                    end
                end
                 list{end + 1} = this.removeTrailingPar(string(index(end)+1 : end));
                 
            end
        end
        function string = removeTrailingPar(this,string)
            for i =1 :length(string)
                if string(i)~=' '
                    if string(i)=='('
                        string(1:i) = [];
                    end
                    break
                end
            end
            for i =length(string) : -1:1
                if string(i)~=' '
                    if string(i)==')'
                        string(i:end) = [];
                    end
                    break
                end
            end
        end
    end
end