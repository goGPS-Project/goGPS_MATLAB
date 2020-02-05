classdef Core_Reference_Frame < handle
    % This class contains properties and methods to manage reference frames
    % and station coordinates
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
    %  Written by: Giulio Tagliaferro
    %  Contributors:     ...
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
    
    
    properties (Constant)
        FLAG_STRING = {'0) rough position', '1) a-priori','2) Fixed', '3) Fixed for PREPRO'};
    end
    
    properties        
        station_code
        xyz
        vxvyvz
        
        end_validity_epoch
        start_validity_epoch
        
        flag
        is_valid
    end
    
    methods
        % Creator
        function this = Core_Reference_Frame()
            % Core object creator
        end
    end
    
    methods
        function init(this, crd_file)
            % initilize the reference frame object loading the crd file specified in state
            %
            % SYNTAX:
            % this.init()
            
            if nargin == 2
                Core.getState.setCrdFile(crd_file);
            end
            
            this.clear();
            this.is_valid = true;
            try
                fname = Core.getState.getCrdFile();
                this.load(fname);
            catch ex
                this.is_valid = false;
                log = Logger.getInstance();
                log.addError('CRD file seems empty or corrupted');
            end
        end
        
        function load(this, crd_file)
            % initilize the reference frame object loading coordinates from
            % crd file
            %
            % SYNTAX:
            % this.init()
            
            this.clear();
            this.is_valid = true;
            if ~isempty(crd_file)
                fid = fopen([crd_file],'r');
                if fid == -1
                    Core.getLogger.addWarning(sprintf('Core RF: File %s not found', crd_file));
                    return
                end
                Core.getLogger.addMessage(Core.getLogger.indent(sprintf('Opening file %s for reading', crd_file)));
                txt = fread(fid,'*char')';
                fclose(fid);
                
                % get new line separators
                nl = regexp(txt, '\n')';
                if nl(end) <  numel(txt)
                    nl = [nl; numel(txt)];
                end
                lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
                lim = [lim lim(:,2) - lim(:,1)];
                if lim(end,3) < 3
                    lim(end,:) = [];
                end
                header_line = find(txt(lim(:,1)) == '#');
                lim(header_line,:) = [];
                %initilaize array
                n_sta = size(lim,1);
                
                this.station_code = {};
                this.xyz = zeros(n_sta,3);
                this.vxvyvz = zeros(n_sta,3);
                this.flag = zeros(n_sta,1);
                re_date = zeros(n_sta,1);
                st_date = zeros(n_sta,1);
                en_date = zeros(n_sta,1);
                for i = 1:n_sta
                    line = txt(lim(i,1):lim(i,2));
                    parts = strsplit(line);
                    l = length(parts);
                    this.station_code(i) = parts(1);
                    this.xyz(i,:) = [str2double(parts{2}) str2double(parts{3}) str2double(parts{4})];
                    if l > 4
                        this.flag(i) = str2double(parts{5});
                        if l > 7
                            this.vxvyvz(i,:) = [str2double(parts{6}) str2double(parts{7}) str2double(parts{8})];
                            if l > 9
                                st_date(i) = datenum([parts{9} ' ' parts{10}]);
                                if l > 11
                                    en_date(i) = datenum([parts{11} ' ' parts{12}]);
                                end
                            end
                        end
                    end
                end
                this.start_validity_epoch = GPS_Time(st_date);
                en_date(en_date == 0)     = datenum(2099, 1,1);
                this.end_validity_epoch   = GPS_Time(en_date);
            end
        end
        
        function clear(this)
            % clear the object
            %
            % SYNTAX:
            % this.clear()
            this.is_valid = false;
            this.xyz = [];
            this.vxvyvz = [];
            this.station_code = {};
            this.flag = [];
            this.start_validity_epoch = [];
            this.end_validity_epoch = [];
        end
        
        function [xyz, is_valid] = getCoo(this, sta_name, epoch)
            % get the coordinates of the station defined by marker name at the desidered epoch
            %
            % SYNTAX:
            %  [xyz, is_valid] = this.getCoo(sta_name, epoch)
            xyz = [];
            
            % load RF if not loaded
            is_valid = false;
            if ~this.isValid()
                if exist(Core.getState.getCrdFile, 'file') == 2
                    this.init();
                end
            end
            
            if this.isValid()
                
                if ~isempty(this.station_code) && length(sta_name) == 4
                    idx_sta = find(strcmpi(this.station_code, sta_name));
                    if sum(idx_sta) > 0
                        st_validity_time = this.start_validity_epoch.getSubSet(idx_sta).getGpsTime();
                        end_validity_time = this.end_validity_epoch.getSubSet(idx_sta).getGpsTime();
                        epoch_gps = epoch.getGpsTime();
                        idx_sta2 = st_validity_time < epoch_gps & end_validity_time > epoch_gps;
                        idx_sta = idx_sta(idx_sta2);
                        idx_sta = idx_sta(1);
                        dt = epoch - this.start_validity_epoch.getEpoch(idx_sta);
                        xyz = this.xyz(idx_sta,:) + (this.vxvyvz(idx_sta,:)' * (dt./(365.25 * 86400))')';
                        is_valid = true;
                    end
                end
            end
        end
        
        function setCoo(this, sta_name, xyz, flag, vxvyvz, start_validity_epoch, end_validity_epoch, flag_overwrite)
            % set the coordiates at the reference epoch
            %
            % SYNTAX:
            %  this.setCoo(sta_name, xyz)
            % load RF if not loaded            
            
            if ~this.isValid()
                if exist(Core.getState.getCrdFile, 'file') == 2
                    this.init();
                end
            end
            if length(sta_name) == 4
                idx_sta = find(strcmpi(this.station_code, sta_name),1 , 'first');
                if sum(idx_sta) > 0
                    this.xyz(idx_sta,:) = xyz;
                    if nargin > 3 && ~isempty(flag)
                        this.flag(idx_sta) = flag;
                    end
                    if nargin > 4 && ~isempty(vxvyvz)
                        this.vxvyvz(idx_sta,:) = vxvyvz;
                    end
                    
                    if nargin > 5 && ~isempty(start_validity_epoch)
                        this.start_validity_epoch.setEpoch(idx_sta,start_validity_epoch);
                    end
                    if nargin > 6 && ~isempty(end_validity_epoch)
                        this.end_validity_epoch.setEpoch(idx_sta,end_validity_epoch);
                    end
                else
                    this.xyz = [this.xyz; xyz];
                    this.station_code{end+1} = sta_name;
                    this.flag = [this.flag; flag];
                    this.vxvyvz = [this.vxvyvz; vxvyvz];
                    if isempty(this.start_validity_epoch)
                        this.start_validity_epoch = start_validity_epoch;
                    else
                        this.start_validity_epoch.append(start_validity_epoch);
                        
                    end
                    if isempty(this.end_validity_epoch)
                        this.end_validity_epoch = end_validity_epoch;
                    else
                        this.end_validity_epoch.append(end_validity_epoch);
                        
                    end
                end
            end
        end
        
        function crx_list = getEntryCell(this)
            % Get the list of CORD entries in the format:
            % 'Marker Name'; 'X'; 'Y'; 'Z'; 'type'; 'start', 'stop'
            %
            % SYNTAX:
            %  [xyz, is_valid] = this.getCoo(sta_name, epoch)
            % load RF if not loaded
            is_valid = false;
            if ~this.isValid()
                if exist(Core.getState.getCrdFile, 'file') == 2
                    this.init();
                end
            end
            
            crx_list = {};
            if this.isValid()
                crx_list = cell(numel(this.station_code), 7);
                for i = 1 : numel(this.station_code)
                    crx_list{i, 1} = this.station_code{i};
                    crx_list{i, 2} = this.xyz(i, 1);
                    crx_list{i, 3} = this.xyz(i, 2);
                    crx_list{i, 4} = this.xyz(i, 3);
                    crx_list{i, 5} = this.FLAG_STRING{nan2zero(this.flag(i)) + 1};
                    crx_list{i, 6} = this.start_validity_epoch.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS');
                    crx_list{i, 7} = this.end_validity_epoch.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS');
                    crx_list{i, 8} = this.vxvyvz(i,1);
                    crx_list{i, 9} = this.vxvyvz(i,2);
                    crx_list{i, 10} = this.vxvyvz(i,3);
                end
            end
        end
        
        function [status] = isValid(this)
            % Tell if station coordiantes are loaded
            %
            % SYNTAX:
            %  [status] = this.isValid()
            status = ~isempty(this.is_valid) && this.is_valid;
        end
        
        function [status] = isFixed(this, sta_code)
            % tell if station coordiantes are meant to be fixed
            % in case sation not sound return false
            %
            % SYNTAX:
            %  [status] = this.isFixed(sta_code)
            status = false;
            if size(this.station_code) > 0
                sta_idx = find(strcmpi(this.station_code, sta_code), 1, 'first');
                if sum(sta_idx) > 0
                    status  = this.flag(sta_idx(1)) == 2;
                end
            end
        end
        
        function [status] = hasAPriori(this, sta_code)
            % tell if station coordiantes are meant to be fixed
            % in case sation not sound return false
            %
            % SYNTAX:
            %  [status] = this.hasAPriori(sta_code)
            status = false;
            if numel(this.station_code) > 0
                sta_idx = find(strcmpi(this.station_code, sta_code), 1, 'first');
                if sum(sta_idx) > 0
                    status  = this.flag(sta_idx(1)) == 2 || this.flag(sta_idx(1)) == 1  || this.flag(sta_idx(1)) == 3;
                end
            end
        end
        
        function [status] = hasGoodAPriori(this, sta_code)
            % tell if station coordiantes are good enough to skip initla code postioining
            % in case sation not sound return false
            %
            % SYNTAX:
            %  [status] = this.hasAPriori(sta_code)
            status = false;
            if numel(this.station_code) > 0
                sta_idx = find(strcmpi(this.station_code, sta_code), 1, 'first');
                if sum(sta_idx) > 0
                    status  = this.flag(sta_idx(1)) == 3 || this.flag(sta_idx(1)) == 2;
                end
            end
        end
        
        function setFlag(this, sta_code, flag)
            % Set the flag [fixed approximate] for the station sta code
            %
            % SYNTAX:
            %  this.setFlag(sta_code,flag)
            sta_idx = find(strcmpi(this.station_code, sta_code), 1, 'first');
            if sum(sta_idx) > 0
                this.flag(sta_idx) = flag;
            end
        end
        
        function flag = getFlag(this, sta_code)
            % Get the flag [fixed / approximate] for the station sta code
            %
            % SYNTAX:
            %  flag = this.getFlag(sta_code)
            sta_idx = find(strcmpi(this.station_code, sta_code), 1, 'first');
            flag = this.flag(sta_idx);
        end
        
        function importTableData(this, data)
            % Import from table (GUI) format to Core_Reference_Frame
            %
            % SYNTAX:
            %  this.importTableData(data)
            
            % get marker names:
            if ~isempty(data)
                name = {};
                for i = 1 : size(data, 1)
                    if ischar(data{i,1})
                        name_start = find(data{i,1} == '>', 1, 'last');
                        name_start = iif(isempty(name_start), 1, name_start + 1);
                        name{i} = data{i,1}(name_start : min(name_start+4, numel(data{i,1}))) ;
                    else
                        name{i} = 'NAME';
                    end
                end
                
                this.station_code = name;
                
                % get location
                this.xyz = [[data{:,2}]' [data{:,3}]' [data{:,4}]'];
                
                % get speed
                this.vxvyvz = [[data{:,8}]' [data{:,9}]' [data{:,10}]'];
                
                % epochs
                date = [];
                for i = 1 : size(data,1)
                    try
                        date(i) = datenum(iif(isempty(data{i,6}), '1980/01/01 00:00:00', data{i,6}), 'yyyy/mm/dd HH:MM:SS');
                    catch ex
                        % not valid epoch
                        date(i) = datenum('1980/01/01 00:00:00', 'yyyy/mm/dd HH:MM:SS');
                    end
                end
                this.start_validity_epoch = GPS_Time(date');
                
                date = [];
                for i = 1 : size(data,1)
                    try
                        date(i) = datenum(iif(isempty(data{i,7}), '2099/01/01 00:00:00', data{i,7}), 'yyyy/mm/dd HH:MM:SS');
                    catch ex
                        % not valid epoch
                        date(i) = datenum('2099/01/01 00:00:00', 'yyyy/mm/dd HH:MM:SS');
                    end
                end
                this.end_validity_epoch = GPS_Time(date');
                
                flag = []; for i = 1 : size(data, 1); flag(i) = iif(isempty(data{i,5}), 0, str2double(data{i,5}(1))); end
                this.flag = flag;
                this.is_valid = 1;
            else
                this.is_valid = 0;
            end
        end
        
        function str = toCrdString(this)
            % Create the string to Export the object
            %
            % SYNTAX:
            %  str = this.toCrdString()
            str = sprintf('#goGPS Coordinate file\n');
            str = sprintf('%s#This file contains position and velocity for multiple stations\n', str);
            str = sprintf('%s#F = FLAG: %s\n', str, sprintf('%s    ', Core_Reference_Frame.FLAG_STRING{:}));
            str = sprintf('%s#-------------------------------------------------------------------------------------------------------------------------------\n', str);
            str = sprintf('%s#STA       X [m]          Y [m]          Z [m]    F   dx [m/y]   dy [m/y]   dz [m/y]   date validity start    date validity stop\n', str);
            str = sprintf('%s#-------------------------------------------------------------------------------------------------------------------------------\n', str);
            for i = 1 : size(this.xyz, 1)
                str = sprintf('%s%4s %+14.5f %+14.5f %+14.5f %1d %+10.5f %+10.5f %+10.5f %s %s\n', str, this.station_code{i}, this.xyz(i, 1), this.xyz(i, 2), this.xyz(i, 3), this.flag(i), this.vxvyvz(i, 1), this.vxvyvz(i, 2), this.vxvyvz(i, 3), this.start_validity_epoch.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS.s'), this.end_validity_epoch.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS.s'));
            end
        end
        
        function export(this, file_path)
            % Export the object in a CRD file
            %
            % SYNTAX:
            %  this.export(file_name)
            
            [path, ~] = fileparts(file_path);
            if ~(exist(path, 'file') == 7)
                mkdir(path);
            end
            fid = fopen(file_path, 'Wb');
            if fid > 0
                fwrite(fid, this.toCrdString, 'char');
                fclose(fid);
            else
                Core.getLogger.addError(sprintf('"%s" cannot be saved', file_path));
            end
        end
    end
    
end
