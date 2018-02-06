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
        function resources_name = getResourcesPath(this, resource_name)
        end       
        function ip = getServerIp(this, name)
            % return the ip of a server given the server name
            ip = [];
            for i = 1 : length(this.servers)
                if strcmp(this.servers{i}.name, name)
                    ip = this.servers{i}.data;
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
                    const = '';
                end
                cond_const = true;
                if nargin > 2
                    cond_const = strfind(const, sys_c);
                end
                if strcmp(name, fname) && cond_const
                    f_struct = struct();
                    f_struct.name = fname;
                    f_struct.const = const;
                    for j = 1 : length(this.file_resources{i}.key)
                        name_k = this.file_resources{i}.key{j}.name;
                        f_struct.(name_k) = this.file_resources{i}.key{j}.data;
                        if strcmp(name_k(1:3),'loc') && ~(name_k(4) == '_')
                            % substitute server name with server address
                            [matches] = regexp(f_struct.(name_k),'(?<=\?{)\w*(?=})','match'); % saerch for ?{server_name} in paths
                            for k = 1 : length(matches)
                                s_ip = this.getServerIp(matches{k});
                                if ~isempty(s_ip)
                                    f_struct.(name_k) = strrep(f_struct.(name_k),['?{' matches{k} '}'], s_ip);
                                else
                                    this.log.addWarning(['No server called ' matches{k} ' found.'])
                                end
                            end
                        end
                    end
                    
                end
            end
        end
        function center_code = getCenterCode(this, center_name, resource_name, sys_c)
            % return the cneter code given a resource name and desidered
            % constalltion
            for i = 1 :length(this.computational_centers)
                if strcmp(this.computational_centers{i}.name, center_name)
                    for j = 1 : length(this.computational_centers{i}.key)
                        resource = this.computational_centers{i}.key{j};
                        if strcmp(resource.name, resource_name)
                            idx = [];
                            if nargin > 3
                                valid = [];
                                for k = 1:length(resource.data)
                                    center_code_part = strsplit(resource.data{k},'@');
                                    if length(center_code_part) > 1
                                        consts = center_code_part{1};
                                        found = true;
                                        for l = 1 : length(sys_c)
                                            found = found && ~isempty(strfind(consts, sys_c(l)));
                                        end
                                        if found
                                            valid = [valid ; [k, length(consts) -length(sys_c)]];
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
    end
end