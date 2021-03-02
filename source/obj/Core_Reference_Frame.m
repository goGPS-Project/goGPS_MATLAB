% This class contains properties and methods to manage reference frames
% and station coordinates

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:        Giulio Tagliaferro
%  Contributors:      Giulio Tagliaferro, Andrea Gatti ...
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

classdef Core_Reference_Frame < handle
    properties (Constant)
        FLAG_STRING = {'0) rough position', '1) a-priori','2) Fixed', '3) Fixed for PREPRO'};
        FLAG_ROUGH = 0;
        FLAG_APRIORI = 1
        FLAG_FIXED = 2
        FLAG_FIXED4PREPRO = 3
    end
    
    properties
        station_code
        xyz     % XYZ coordinates
        vxvyvz  % XYZ velocities
        std_pup % std planar up
        
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
                log.addWarning('CRD file seems empty or corrupted');
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
                fid = fopen([crd_file],'rt');
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
                this.std_pup = zeros(n_sta,2);
                this.flag = zeros(n_sta,1);
                st_date = zeros(n_sta,1);
                en_date = zeros(n_sta,1);
                for i = 1:n_sta
                    line = txt(lim(i,1):lim(i,2));
                    parts = strsplit(line);
                    l = length(parts);
                    this.station_code(i) = parts(1);
                    this.xyz(i,:) = [str2double(parts{2}) str2double(parts{3}) str2double(parts{4})];
                    flag_std = l >= 10 && ~isnan(str2double(parts{9}));
                    if l > 4
                        this.flag(i) = str2double(parts{5});
                        if l > 7
                            this.vxvyvz(i,:) = [str2double(parts{6}) str2double(parts{7}) str2double(parts{8})];
                            if flag_std
                                if l > 9
                                    this.std_pup(i,:) = [str2double(parts{9}) str2double(parts{10})];
                                    if l > 11
                                        st_date(i) = datenum([parts{11} ' ' parts{12}]);
                                        if l > 13
                                            en_date(i) = datenum([parts{13} ' ' parts{14}]);
                                        end
                                    end
                                end
                            end
                        elseif l > 9
                            st_date(i) = datenum([parts{9} ' ' parts{10}]);
                            if l > 11
                                en_date(i) = datenum([parts{11} ' ' parts{12}]);
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
        
        function [xyz, is_valid, std_pup] = getCoo(this, sta_name, epoch)
            % get the coordinates of the station defined by marker name at the desidered epoch
            %
            % SYNTAX:
            %  [xyz, is_valid] = this.getCoo(sta_name, epoch)
            xyz = [];
            std_pup = [];
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
                        if nargin > 2 && not(isempty(epoch))
                            epoch_gps = epoch.getGpsTime();
                            idx_sta2 = st_validity_time < epoch_gps & end_validity_time > epoch_gps;
                            if not(any(idx_sta2))
                                log = Core.getLogger;
                                log.addMessage(log.indent('No valid a-priori coordinate found.'));
                            end
                            idx_sta = idx_sta(idx_sta2);
                        end
                        if not(isempty(idx_sta))
                            idx_sta = idx_sta(end); % take only the last coordinate matching
                            
                            % remove velocity
                            if nargin > 2 && not(isempty(epoch))
                                dt = epoch - this.start_validity_epoch.getEpoch(idx_sta);
                                xyz = this.xyz(idx_sta,:) + (this.vxvyvz(idx_sta,:)' * (dt./(365.25 * 86400))')';
                            else
                                xyz = this.xyz(idx_sta,:);
                            end
                            
                            is_valid = true;
                            std_pup = this.std_pup(idx_sta,:);
                        end
                    end
                end
            end
        end
        
        function setCoo(this, sta_name, xyz, flag, vxvyvz, std_pup, start_validity_epoch, end_validity_epoch, flag_overwrite)
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
                    
                    if nargin > 5 && ~isempty(std_pup)
                        this.std_pup(idx_sta,:) = std_pup;
                    end
                    if nargin > 6 && ~isempty(start_validity_epoch)
                        this.start_validity_epoch.setEpoch(idx_sta,start_validity_epoch);
                    end
                    if nargin > 7 && ~isempty(end_validity_epoch)
                        this.end_validity_epoch.setEpoch(idx_sta,end_validity_epoch);
                    end
                else
                    this.xyz = [this.xyz; xyz];
                    this.station_code{end+1} = sta_name;
                    this.flag = [this.flag; flag];
                    this.vxvyvz = [this.vxvyvz; vxvyvz];
                    this.std_pup = [this.std_pup; std_pup];
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
            % 'Marker Name'; 'X'; 'Y'; 'Z'; 'type'; 'std Planar'; 'std Up'; 'start'; 'stop'
            %
            % SYNTAX:
            %  [xyz, is_valid] = this.getCoo(sta_name, epoch)
            % load RF if not loaded
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
                    crx_list{i, 6} = this.std_pup(i, 1);
                    crx_list{i, 7} = this.std_pup(i, 2);
                    crx_list{i, 8} = this.start_validity_epoch.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS');
                    crx_list{i, 9} = this.end_validity_epoch.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS');
                    crx_list{i, 10} = this.vxvyvz(i,1);
                    crx_list{i, 11} = this.vxvyvz(i,2);
                    crx_list{i, 12} = this.vxvyvz(i,3);
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
        
        function [status] = isFixedPrepro(this, sta_code)
            % tell if station coordiantes are meant to be fixed
            % in case sation not sound return false
            %
            % SYNTAX:
            %  [status] = this.isFixed(sta_code)
            status = false;
            if size(this.station_code) > 0
                sta_idx = find(strcmpi(this.station_code, sta_code), 1, 'first');
                if sum(sta_idx) > 0
                    status  = this.flag(sta_idx(1)) == 3;
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
                this.vxvyvz = [[data{:,10}]' [data{:,11}]' [data{:,12}]'];
                
                % get speed
                this.std_pup = [[data{:,6}]' [data{:,7}]'];
                
                % epochs
                date = [];
                for i = 1 : size(data,1)
                    try
                        date(i) = datenum(strrep(iif(isempty(data{i,8}), '2099/01/01 00:00:00', data{i,8}), '-', '/'), 'yyyy/mm/dd HH:MM:SS');
                    catch ex
                        % not valid epoch
                        Core_Utils.printEx(ex);
                        date(i) = datenum('1980/01/01 00:00:00', 'yyyy/mm/dd HH:MM:SS');
                    end
                end
                this.start_validity_epoch = GPS_Time(date');
                
                date = [];
                for i = 1 : size(data,1)
                    try
                        date(i) = datenum(strrep(iif(isempty(data{i,9}), '2099/01/01 00:00:00', data{i,9}), '-', '/'), 'yyyy/mm/dd HH:MM:SS');
                    catch ex
                        % not valid epoch
                        Core_Utils.printEx(ex);
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
            str = sprintf('%s#--------------------------------------------------------------------------------------------------------------------------------------------------\n', str);
            str = sprintf('%s#STA             X             Y              Z  F       dX       dY        dZ   std Planar    std Up   date validity start   date validity stop\n', str);
            str = sprintf('%s#               [m]           [m]            [m]       [m/y]    [m/y]     [m/y]         [m]       [m]  yyyy-mm-dd HH:MM:SS.s yyyy-mm-dd HH:MM:SS.s\n', str);
            str = sprintf('%s#--------------------------------------------------------------------------------------------------------------------------------------------------\n', str);
            for i = 1 : size(this.xyz, 1)
                str = sprintf('%s%4s %+14.5f %+14.5f %+14.5f %1d %+9.5f %+9.5f %+9.5f  %+9.5f %+9.5f  %s %s\n', str, this.station_code{i}, this.xyz(i, 1), this.xyz(i, 2), this.xyz(i, 3), this.flag(i), this.vxvyvz(i, 1), this.vxvyvz(i, 2), this.vxvyvz(i, 3),  this.std_pup(i, 1),  this.std_pup(i, 2),  this.start_validity_epoch.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS.s'), this.end_validity_epoch.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS.s'));
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
            Core.getLogger.addMessage(Core.getLogger.indent(sprintf('Opening file %s for writing', file_path)));
            fid = fopen(file_path, 'Wb');
            if fid > 0
                fwrite(fid, this.toCrdString, 'char');
                fclose(fid);
            else
                Core.getLogger.addError(sprintf('"%s" cannot be saved', file_path));
            end
        end
        
        function fh_list = showMapGoogle(this)
            % Show Google Map of the stations
            %
            % CITATION:
            %   Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software],
            %   available online at www.eoas.ubc.ca/~rich/map.html.
            %
            % INPUT
            %   new_fig     open a new figure
            %
            % SYNTAX
            %   sta_list.showMapGoogle(new_fig);
           
            N_MAX_STA = 50;
            
            flag_labels = true;
            flag_large_points = true;
            point_size = 15;
            point_color = [0, 255, 10]/256;
            
            f = figure('Visible', 'off');
            
            fh_list = f;
            fig_name = sprintf('RFMapGoogle');
            f.UserData = struct('fig_name', fig_name);

            Core.getLogger.addMarkedMessage('Preparing map, please wait...');
            
            f.Color = [1 1 1];
            [lat, lon] = Coordinates.fromXYZ(this.xyz).getGeodetic;
            lat = lat / pi * 180;
            lon = lon / pi * 180;
            
            lat_tmp = lat;
            lon_tmp = lon;
            
            % set map limits
            if numel(lon_tmp) == 1
                lon_lim = minMax(lon_tmp) + [-0.05 0.05];
                lat_lim = minMax(lat_tmp) + [-0.05 0.05];
            else
                lon_lim = minMax(lon_tmp); lon_lim = lon_lim + [-1 1] * diff(lon_lim) / 6;
                lat_lim = minMax(lat_tmp); lat_lim = lat_lim + [-1 1] * diff(lat_lim) / 6;
            end
            lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
            lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
            nwse = [lat_lim(2), lon_lim(1), lat_lim(1), lon_lim(2)];
            clon = nwse([2 4]) + [-0.001 0.001];
            clat = max(-90, min(90, nwse([3 1]) + [-0.001 0.001]));

            axes
            xlim(clon);
            ylim(clat);
            [lon_ggl,lat_ggl, img_ggl] = Core_Utils.addGoogleMaps('alpha', 0.95, 'maptype','satellite','refresh',0,'autoaxis',0);
            xlim(lon_lim);
            ylim(lat_lim);
            
            m_proj('equidistant','lon',clon,'lat',clat);   % Projection
            %m_proj('utm', 'lon',lon_lim,'lat',lat_lim);   % Projection
            drawnow
            m_image(lon_ggl, lat_ggl, img_ggl);
            
            % read shapefile
            shape = 'none';
            if (~strcmp(shape,'none'))
                if (~strcmp(shape,'coast')) && (~strcmp(shape,'fill'))
                    if (strcmp(shape,'10m'))
                        M = m_shaperead('countries_10m');
                    elseif (strcmp(shape,'30m'))
                        M = m_shaperead('countries_30m');
                    else
                        M = m_shaperead('countries_50m');
                    end
                    [x_min, y_min] = m_ll2xy(min(lon_lim), min(lat_lim));
                    [x_max, y_max] = m_ll2xy(max(lon_lim), max(lat_lim));
                    for k = 1 : length(M.ncst)
                        lam_c = M.ncst{k}(:,1);
                        ids = lam_c <  min(lon);
                        lam_c(ids) = lam_c(ids) + 360;
                        phi_c = M.ncst{k}(:,2);
                        [x, y] = m_ll2xy(lam_c, phi_c);
                        if sum(~isnan(x))>1
                            x(find(abs(diff(x)) >= abs(x_max - x_min) * 0.90) + 1) = nan; % Remove lines that occupy more than th 90% of the plot
                            line(x,y,'color', [0.3 0.3 0.3]);
                        end
                    end
                else
                    if (strcmp(shape,'coast'))
                        m_coast('line','color', lineCol);
                    else
                        m_coast('patch',lineCol);
                    end
                end
            end
            hold on;
            
            m_grid('box','fancy','tickdir','in', 'fontsize', 16);
            % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
            drawnow
            %m_ruler([.7 1], -0.05, 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
            [x, y] = m_ll2xy(lon, lat);
            
            
            %point_color = Cmap.get('viridis', numel(x));
            %point_size = 25;
            if size(point_color, 1) > 1
                scatter(x(:), y(:), point_size, 1:numel(x), 'filled'); hold on;
                colormap(point_color);
            else
                if numel(this.station_code) > N_MAX_STA
                    plot(x(:), y(:),'.', 'MarkerSize', point_size, 'Color', Core_UI.BLACK); hold on;
                    plot(x(:), y(:),'.', 'MarkerSize', point_size-3, 'Color', point_color); hold on;
                else
                    plot(x(:), y(:),'.', 'MarkerSize', 8, 'Color', Core_UI.BLACK); hold on;
                end
            end
            if flag_labels
                % Label BG (in background w.r.t. the point)
                if numel(this.station_code) < N_MAX_STA
                    for r = 1 : numel(this.station_code)
                        name = upper(this.station_code{r});
                        txt = text(x(r), y(r), ['   ' name], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                            'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                    end
                end
            end
            
            if flag_large_points && numel(this) < N_MAX_STA
                for r = 1 : numel(this.station_code)
                    plot(x(r), y(r), '.', 'MarkerSize', 45, 'Color', Core_UI.getColor(r, numel(this)), 'UserData', 'GNSS_point');
                end
                plot(x(:), y(:), '.k', 'MarkerSize', 5);
                plot(x(:), y(:), 'ko', 'MarkerSize', 15, 'LineWidth', 2);
            end
            if flag_labels
                if numel(this.station_code) < N_MAX_STA
                    for r = 1 : numel(this.station_code)
                        name = upper(this.station_code{r});
                        t = text(x(r), y(r), ['   ' name], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                            ...%'FontWeight', 'bold', 'FontSize', 10, 'Color', [0 0 0], ...
                            ...%'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                        %t.Units = 'pixels';
                        %t.Position(1) = t.Position(1) + 20 + 10 * double(numel(sta_list) == 1);
                        %t.Units = 'data';
                    end
                end
            end
            
            Core_UI.addExportMenu(f); Core_UI.addBeautifyMenu(f); Core_UI.beautifyFig(f);
            f.Visible = 'on'; drawnow;
            title(sprintf('Map of GNSS stations\\fontsize{5} \n'), 'FontSize', 16);
            %xlabel('Longitude [deg]');
            %ylabel('Latitude [deg]');
            ax = gca; ax.FontSize = 16;
            Core.getLogger.addStatusOk('The map is ready ^_^');
        end

        
    end
    
    methods (Static)
        function test()
            rf = Core_Reference_Frame;      
            rf.load('../data/project/default_PPP/station/CRD/stations.crd');
            rf.export('../data/project/default_PPP/station/CRD/stations_test.crd')
        end
    end
    
end
