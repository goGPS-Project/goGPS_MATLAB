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
    %  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
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
        state
        log
       
        station_code
        xyz
        vxvyvz
        
        ref_epoch
        end_validity_epoch
        start_validity_epoch
        
        flag
        is_valid
    end
    
    methods
        % Creator
        function this = Core_Reference_Frame()
            % Core object creator
            this.state = Core.getCurrentSettings();
            this.log = Logger.getInstance();
        end
    end
    
    methods
        function init(this, crd_file)
            % initilize the reference frame object loading the crd file specified in state
            %
            % SYNTAX:
            % this.init()
            
            this.state = Core.getCurrentSettings();
            if nargin == 2
                this.state.setCrdFile(crd_file);
            end
            
            this.clear();
            this.is_valid = true;
            try
                fname = this.state.getCrdFile();
                if ~isempty(fname)
                    fid = fopen([fname],'r');
                    if fid == -1
                        this.log.addWarning(sprintf('Core RF: File %s not found', fname));
                        return
                    end
                    this.log.addMessage(this.log.indent(sprintf('Opening file %s for reading', fname)));
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
                    for i = header_line
                        line = txt(lim(i,1):lim(i,2));
                        idx = strfind(line,'Reference epoch');
                        if ~isempty(idx)
                            this.ref_epoch = GPS_Time( sscanf(line((idx+15):end),'%f %f %f %f %f %f')');
                        end
                    end
                    lim(header_line,:) = [];
                    %initilaize array
                    n_sta = size(lim,1);
                    
                    this.station_code=char(' '*ones(n_sta,4));
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
                        this.station_code(i,:) = parts{1};
                        this.xyz(i,:) = [str2double(parts{2}) str2double(parts{3}) str2double(parts{4})];
                        if l > 4
                            this.flag(i) = str2double(parts{5});
                            if l > 7
                                this.vxvyvz(i,:) = [str2double(parts{6}) str2double(parts{7}) str2double(parts{8})];
                                if l > 9
                                    re_date(i) = datenum([parts{9} ' ' parts{10}]);
                                    if l > 11
                                        st_date(i) = datenum([parts{11} ' ' parts{12}]);
                                        if l > 13
                                            en_date(i) = datenum([parts{13} ' ' parts{14}]);
                                        end
                                    end
                                end
                            end
                        end
                    end
                    this.ref_epoch            = GPS_Time(re_date);
                    this.start_validity_epoch = GPS_Time(st_date);
                    en_date(en_date == 0)     = datenum(2099, 1,1);
                    this.end_validity_epoch   = GPS_Time(en_date);
                end
            catch
                this.is_valid = false;
                log = Logger.getInstance();
                log.addError('CRD file seems to be corrupted');
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
            this.station_code = [];
            this.flag = [];
            this.start_validity_epoch = [];
            this.end_validity_epoch = [];
            this.ref_epoch = [];
        end
        
        function [xyz, is_valid] = getCoo(this, sta_name, epoch)
            % get the coordinates of the station defined by marker name at the  desidered epoch
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
                    idx_sta = strLineMatch(lower(this.station_code), lower(sta_name));
                    if sum(idx_sta) > 0
                        st_validity_time = this.start_validity_epoch.getSubSet(idx_sta).getGpsTime();
                        end_validity_time = this.end_validity_epoch.getSubSet(idx_sta).getGpsTime();
                        epoch_gps = epoch.getGpsTime();
                        idx_sta2 = st_validity_time < epoch_gps & end_validity_time > epoch_gps;
                        idx_sta = find(idx_sta);
                        idx_sta = idx_sta(idx_sta2);
                        dt = epoch - this.ref_epoch.getEpoch(idx_sta);
                        xyz = this.xyz(idx_sta,:) + (this.vxvyvz(idx_sta,:)' * (dt./(365.25 * 86400))')';
                        is_valid = true;
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
                crx_list = cell(size(this.station_code, 1), 7);
                for i = 1 : size(this.station_code, 1)
                    crx_list{i, 1} = this.station_code(i, :);
                    crx_list{i, 2} = this.xyz(i, 1);
                    crx_list{i, 3} = this.xyz(i, 2);
                    crx_list{i, 4} = this.xyz(i, 3);
                    crx_list{i, 5} = this.FLAG_STRING{this.flag(i) + 1};
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
                sta_idx  = strLineMatch(lower(this.station_code), lower(sta_code));
                if sum(sta_idx) > 0
                    status  = this.flag(sta_idx) == 2;
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
            if size(this.station_code) > 0
                sta_idx  = strLineMatch(lower(this.station_code), lower(sta_code));
                if sum(sta_idx) > 0
                    status  = this.flag(sta_idx) == 2 || this.flag(sta_idx) == 1  || this.flag(sta_idx) == 3;
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
            if size(this.station_code) > 0
                sta_idx  = strLineMatch(lower(this.station_code), lower(sta_code));
                if sum(sta_idx) > 0
                    status  = this.flag(sta_idx) == 3 || this.flag(sta_idx) == 2;
                end
            end
        end
        
        function setFlag(this, sta_code, flag)
            % Set the falg [fixed approximate] for the station sta code
            %
            % SYNTAX:
            %  this.setFlag(sta_code,flag)
            sta_idx  = strLineMatch(this.station_code, sta_code);
            if sum(sta_idx) > 0
                this.flag(sta_idx) = flag;
            end
        end
    end
    
end
