%   CLASS Core_Antenna
% =========================================================================
%
% DESCRIPTION
%   Collector of antenna models
%   Container of the ATX file()
%
% EXAMPLE
%   atx = Core_Antenna();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_Antenna

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 (GReD srl) Andrea Gatti, Giulio Tagliaferro, Eugenio Realini
%  Written by:        Andrea Gatti, Giulio Tagliaferro ...
%  Contributors:      Andrea Gatti, Giulio Tagliaferro ...
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

classdef Core_Antenna < handle
    
    %% PROPERTIES CONSTANTS
    % ==================================================================================================================================================
    properties (Constant)
    end
    
    %% PROPERTIES MISC
    % ==================================================================================================================================================
    properties (GetAccess = private, SetAccess = private) % Public Access
        creation_time = GPS_Time(now);
    end
    
    %% PROPERTIES ANTENNA
    % ==================================================================================================================================================
    properties
        list        % Array of Antenna objects
        type        % List of antenna names present in ant_list (mirrored copy) (for look-up)
        serial      % List of Antennas serial number (for look-up)
        start       % Validity interval (MATLAB TIME)
        stop        % Validity interval (MATLAB TIME)
        
        mnt_tbl     % monument table (Used only for GeoNet)
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = public)
        % Concrete implementation.  See Singleton superclass.
        function this = Core_Antenna()
            % Core object creator
        end
    end
    
    %% METHODS INIT & STATIC GETTERS & SETTERS
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = fromAntex(file_name)
            this = Core_Antenna();
            this.init();
            this.importAntex(file_name);
        end
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this, atx_file_name)
            this.type = '';
            this.serial = '';
            this.list = [];
            if nargin == 2
                % Load custom files:
                file_list = dir(fullfile(Core.getState.getHomeDir, 'station', 'ATX', 'custom', '*.ATX'));
                file_list = [file_list; dir(fullfile(Core.getState.getHomeDir, 'station', 'ATX', 'custom', '*.atx'))];
                file_name = {};
                log = Logger.getInstance();
                for f = 1 : numel(file_list)
                    file_name{f} = fullfile(Core.getState.getHomeDir, 'station', 'ATX', 'custom', file_list(f).name); %#ok<AGROW>
                    log.addMessage('Custom ANTEX file found in "%s"\n Loading...', file_name{f});
                    this.importAntex(file_name{f});
                end
                
                % Load UI defined file
                this.importAntex(atx_file_name);
            end
            this.importMonumentTable()
            this.importGeoNetAntex()
        end
        
        function importAntex(this, file_name)
            if nargin == 0
                file_name = Core.getState.getAtxFile();
            end
            log = Core.getLogger();
            
            % open RINEX observation file
            fid = fopen(file_name,'rt');
            if fid > 0
                txt = fread(fid,'*char')';
                % try to see if carriage return is present in the file (Windows stupid standard)
                % On Windows file lines ends with char(13) char(10)
                % instead of just using char(10)
                if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                    has_cr = true;  % The file has carriage return - I hate you Bill!
                else
                    has_cr = false;  % The file is UNIX standard
                end
                % txt = txt(txt ~= 13);  % remove carriage return - I hate you Bill!
                fclose(fid);
                
                % get new line separators
                nl = regexp(txt, '\n')';
                if nl(end) <  (numel(txt) - double(has_cr))
                    nl = [nl; numel(txt)];
                end
                lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                lim = [lim (lim(:,2) - lim(:,1) + 1)];
                lim(lim(:,3) < 64, :) = []; % Do not consider any line shorter than 65 characters
                lim((txt((lim(:,1) + 60)) == 'C'), :) = []; % detect COMMENT and ignore them
                
                id_start = find((txt((lim(:,1) + 60)) == 'S') & (txt((lim(:,1) + 69)) == 'A'));  % detect START OF ANTENNA
                id_stop = find((txt((lim(:,1) + 60)) == 'E') & (txt((lim(:,1) + 69)) == 'T'));   % detect END OF ANTENNA
                id_ant = find((txt((lim(:,1) + 60)) == 'T'));
                n_ant = numel(id_start);
                if n_ant ~= numel(id_stop) || n_ant ~= numel(id_ant)
                    log.addError('Corrupted antex file!\n antennas are not properly opened or closed by (START / END)');
                else
                    ant_type = txt(repmat(lim(id_ant, 1), 1, 20) + repmat(0:19, n_ant, 1));
                    ant_serial = txt(repmat(lim(id_ant, 1), 1, 20) + repmat(20:39, n_ant, 1)); % Used for Satellites antenna
                    ant_list = Antenna();
                    for a = 1 : n_ant
                        ant_list(a) = Antenna.fromText(txt(lim(id_start(a) + 1, 1) : lim(id_stop(a) - 1, 2)), lim((id_start(a) + 1) : (id_stop(a) - 1), :), has_cr);
                    end
                end
                
                % Append imported antennas
                if isempty(this.type)
                    this.type = ant_type;
                    this.serial = ant_serial;
                    this.list = ant_list;
                    tmp = [ant_list.start];
                    this.start = tmp.getMatlabTime;
                    tmp = [ant_list.stop];
                    this.stop = tmp.getMatlabTime;
                else
                    this.type = [this.type; ant_type];
                    this.serial = [this.serial; ant_serial];
                    this.list = [this.list ant_list];
                    tmp = [ant_list.start];
                    this.start = [this.start; tmp.getMatlabTime];
                    tmp = [ant_list.stop];
                    this.stop = [this.stop; tmp.getMatlabTime];
                end
            end
        end
        
        function importGeoNetAntex(this, file_name)
            if nargin == 1
                file_list = dir(fullfile(Core.getState.getHomeDir, 'station', 'ATX', 'custom', 'GEONET*'));
                file_name = {};
                for f = 1 : numel(file_list)
                    file_name{f} = fullfile(Core.getState.getHomeDir, 'station', 'ATX', 'custom', file_list(f).name); %#ok<AGROW>
                end
            else
                if ~iscell(file_name)
                    file_name = {file_name};
                end
            end
            log = Core.getLogger();
            
            for f = 1 : numel(file_name)
                % open RINEX observation file
                fid = fopen(file_name{f},'rt');
                if fid > 0
                    txt = fread(fid, '*char')';
                    % try to see if carriage return is present in the file (Windows stupid standard)
                    % On Windows file lines ends with char(13) char(10)
                    % instead of just using char(10)
                    if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                        has_cr = true;  % The file has carriage return - I hate you Bill!
                    else
                        has_cr = false;  % The file is UNIX standard
                    end
                    % txt = txt(txt ~= 13);  % remove carriage return - I hate you Bill!
                    fclose(fid);
                    
                    % get new line separators
                    nl = regexp(txt, '\n')';
                    if nl(end) <  (numel(txt) - double(has_cr))
                        nl = [nl; numel(txt)];
                    end
                    lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                    lim = [lim (lim(:,2) - lim(:,1) + 1)];
                    lim(lim(:,3) < 21, :) = []; % Do not consider any line shorter than 65 characters
                    
                    % in character 39 of the line PCO lines contains the frequuency number
                    pco_id = find(txt(lim(:,1) + 38)' > '0' & txt(lim(:,1) + 38)' < '9'); % lines containing PCO
                    type_id = find(txt(lim(:,1) + 38)' == '1'); % lines antenna names
                    ant_type = txt(repmat(lim(type_id,1), 1, 20) + repmat(0:19, size(type_id, 1), 1));
                    pco_size_per_type = diff([type_id; (pco_id(end) +1)]);
                    % Before an antenna there is always a line of "*"                        
                    % get antenna PCV start line
                    ant_start = (find(txt(lim(:,1)) == '*') + 1); ant_start = ant_start(2:end)';
                    
                    if size(ant_type, 1) ~= size(ant_start, 1) && ... % number of PCV == PCO
                        sum(pco_size_per_type == 2) == size(ant_type, 1) % I should have 2 frequencies per antenna
                        % the file seems to be corrupted 
                        log.addWarning('The GeoNet antenna file seems to be corrupted, skipping it');
                    else
                        n_ant = size(ant_start, 1);
                        % import PCO
                        serial_pco = reshape(sscanf(serialize(txt(repmat(lim(pco_id, 1), 1, 24) + repmat(41:64, size(pco_id, 1), 1))')', '%f'), 3, 2 * n_ant)'  * 1e3;
                        for a = 1 : n_ant
                            % PCO in the geonet config file is in meters => convert it in mm (*1e3)
                            pco{a} = [{serial_pco(2*a -1,:)} {serial_pco(2*a,:)}]; %#ok<AGROW> 
                        end      
                        
                        % import antenna calibration date (when available)
                        cal_date = txt(repmat(lim(ant_start, 1), 1, 9) + repmat(111:119, n_ant, 1)); cal_date((lim(ant_start, 3) < 120), :) = ' ';
                        cal_flag = find(sum(cal_date, 2) ~= ' ' * size(cal_date, 2))';
                        for a = find(sum(cal_date, 2) == ' ' * size(cal_date, 2))'
                            cal_date(a, :) = '06-JAN-80';
                        end
                        cal_date = datenum(cal_date, 'dd-mmm-yy');
                        
                        % read dazi
                        dazi = str2num(txt(repmat(lim(ant_start, 1), 1, 3) + repmat(69:71, n_ant, 1)));
                        dazi(dazi == 360) = 0;
                        az_grid = cell(n_ant, 1);
                        for a = find(dazi > 0)'
                            az_grid{a} = (0 : dazi(a) : 360)';
                        end
                        
                        % read zen2 dzen
                        zen1 = zeros(n_ant, 1);
                        zen2 = str2num(txt(repmat(lim(ant_start, 1), 1, 2) + repmat(75:76, n_ant, 1)));
                        dzen = str2num(txt(repmat(lim(ant_start, 1), 1, 2) + repmat(65:66, n_ant, 1)));
                        
                        % n_freq
                        n_freq = pco_size_per_type;
                        
                        % f_code
                        f_code = ['G01';'G02']; % at the moment I know and I checked that all the antenna have just two frequencies of calibration                        end
                        
                        % Start-Stop
                        start = cal_date;
                        stop = cal_date + 365*100; % One hundred years of validity is more than enough
                        
                        % Sinex_code
                        sinex_code = cell(n_ant, 1);
                        for a = cal_flag
                            sinex_code{a} = txt(repmat(lim(ant_start(a), 1), 1, 10) + repmat(79:88, 1, 1));
                        end                        
                        
                        % Read PCV
                        pcv_noaz = cell(n_ant, 1);
                        pcv = cell(n_ant, 1);
                        for a = 1 : n_ant
                            if dazi(a) > 0
                                rows = lim(ant_start(a) + 1 + (1 : n_freq(a) * numel(az_grid{a})), 1);
                                data = sscanf(serialize((txt(repmat(rows, 1, 132) + repmat(8:139, size(rows, 1), 1)))')', '%f');
                                data = reshape(data, size(zen1(a) : dzen(a) : 90, 2), size(rows, 1))';
                                data = data(:, 1 : (zen2(a)/dzen(a) + 1));
                                pcv{a} = [{data(1 : 2 : end, :)} {data(2 : 2 : end, :)}];
                                pcv_noaz{a} = [{mean(pcv{a}{1})} {mean(pcv{a}{2})}];
                            else
                                rows = lim(ant_start(a) + 1 + (1 : n_freq(a)), 1);
                                data = sscanf(serialize((txt(repmat(rows, 1, 132) + repmat(8:139, size(rows, 1), 1)))')', '%f');
                                data = reshape(data, size(zen1(a) : dzen(a) : 90, 2), size(rows, 1))';
                                data = data(:, 1 : (zen2(a)/dzen(a) + 1));
                                pcv_noaz{a} = [{data(1,:)} {data(2,:)}];
                            end
                        end
                        
                        ant_serial = 'CUSTOM              ';
                        ant_list(n_ant) = Antenna; %#ok<AGROW>
                        for a = 1 : n_ant
                            ant_list(a).import(ant_type(a, :), ...
                                ant_serial, ...
                                GPS_Time(cal_date(a)), ...
                                dazi(a), ...
                                az_grid{a}, ...
                                zen1(a), ...
                                zen2(a), ...
                                dzen(a), ...
                                n_freq(a), ...
                                f_code, ...
                                GPS_Time(start(a)), ...
                                GPS_Time(stop(a)), ...
                                sinex_code{a}, ...
                                pco{a}, ...
                                pcv_noaz{a}, ... 
                                pcv{a});
                        end
                        
                        % Append imported antennas
                        if isempty(this.type)
                            this.type = ant_type;
                            this.serial = repmat(ant_serial, n_ant, 1);
                            this.list = ant_list;
                            this.start = start;
                            this.stop = tmp.getMatlabTime;
                        else
                            this.type = [this.type; ant_type];
                            this.serial = [this.serial; repmat(ant_serial, n_ant, 1)];
                            this.list = [this.list ant_list];
                            tmp = [ant_list.start];
                            this.start = [this.start; start];
                            this.stop = [this.stop; stop];
                        end
                    end
                end
            end
        end
        
        function importMonumentTable(this, file_name)
            % Import monument table for GEONET Japan receivers
            %
            % SYNTAX
            %   this.importMonumentTable(file_name)            
            
            if nargin < 2
                file_name = fullfile(Core.getState.getHomeDir, 'station', 'ATX', 'custom', 'MONUMENT.TBL');
            end
            %log = Core.getLogger();
            
            % open RINEX observation file
            fid = fopen(file_name,'rt');
            if fid > 0
                txt = fread(fid,'*char')';
                % try to see if carriage return is present in the file (Windows stupid standard)
                % On Windows file lines ends with char(13) char(10)
                % instead of just using char(10)
                if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                    has_cr = true;  % The file has carriage return - I hate you Bill!
                else
                    has_cr = false;  % The file is UNIX standard
                end
                % txt = txt(txt ~= 13);  % remove carriage return - I hate you Bill!
                fclose(fid);
                
                % get new line separators
                nl = regexp(txt, '\n')';
                if nl(end) <  (numel(txt) - double(has_cr))
                    nl = [nl; numel(txt)];
                end
                lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                lim = [lim (lim(:,2) - lim(:,1) + 1)];
                
                % Start to parse monument.tbl
                sep = regexp(txt, ';')';
                txt(lim(2:end,1) - 1) = ';';
                tmp = textscan(txt(lim(2,1) : end), '%s %s %s %s %d %d %f %f %f %s', 'Delimiter', ';');
                                
                n_entry = size(tmp{1}, 1);
                % Supposing station marker starting from character 2 to 5
                %sta_name = txt(repmat(lim(2 : end, 1), 1, 4) + repmat((2:5), size(lim,1) - 1, 1));
                % Convert station ID to marker name (keep last 4 char
                sta_name = char(zeros(n_entry, 4, 'uint8'));
                for e = 1 : n_entry
                    sta_name(e, :) = tmp{1}{e}(end-3:end);
                end
                monument = txt(repmat(sep(12:9:end), 1, 4) + repmat((1:4), size(lim,1) - 1, 1));
                monument(monument(:,1) == '?', :) = ' ';
                time = GPS_Time.fromDoySod(double(tmp{5}), double(tmp{6}), double(0*tmp{5}));                
                
                % compose antenna name
                antenna_type = char(32 * ones(n_entry, 20, 'uint8'));
                for a = 1 : n_entry
                    antenna_type(a, 1 : length(tmp{3}{a})) = tmp{3}{a};
                    antenna_type(a, 17:20) = monument(a, :);
                end
                this.mnt_tbl = struct('marker', sta_name, 'type', antenna_type, 'time', time);
            end
        end
    end
    
    %% METHODS UTILITIES
    % ==================================================================================================================================================
    methods
        function sat_type = getSatAntennaType(this, serial, time)
            % get the satellite type given a 3ch antenna (e.g. G01, E12, ...)
            %
            % INPUT
            %   serial      3ch serial value
            %   time        GPS_Time reference time for antenna validity
            %
            % SYNTAX
            %   sat_type = this.getSatAntennaType(serial, time)
            if ~iscell(serial)
                serial = {serial};
            end
            time = time.getMatlabTime;
            if isempty(this.serial)
                Core.getLogger.addError('goGPS requires a satellite antenna file to work properly!'); 
                error('No antennas no fun :-(')
            else
                ss = Core_Utils.code3Char2Num(this.serial(:, 1:3)); % stored serials expressed as a number
                sat_type = cell(size(serial(:), 1), 1);
                for a = 1 : size(serial(:), 1)
                    id_sat = find(ss == Core_Utils.code3Char2Num(serial{a}(1:3)) & ...
                        this.start <= time & ...
                        this.stop >= time, 1, 'first');
                    sat_type{a} = this.type(id_sat, :);
                end
            end
        end
        
        function ant = getAntenna(this, type, serial, time, log_lev)
            % get the satellite type given a 3ch antenna (e.g. G01, E12, ...)
            %
            % INPUT
            %   serial      3ch serial value
            %   time        GPS_Time reference time for antenna validity
            %
            % SYNTAX
            %   sat_type = this.getAntenna(type, serial, time)
            if ~isempty(type) && ~iscell(type)
                type = {type};
            end
            if ~isempty(serial) && ~iscell(serial)
                serial = {serial};
            end
            
            log = Core.getLogger;
            
            time = time.getMatlabTime;
            % time check
            t_check = find(this.start <= time & this.stop >= time);
            
            n_ant = max(size(serial(:), 1), size(type(:), 1));
            none_code = Core_Utils.code4Char2Num('NONE');
            ant(n_ant + 1) = Antenna(); ant(end) = []; % init output, define empty array of Antennas
            for a = 1 : n_ant
                if isempty(type) || isempty(type{a})
                    cmp_str = this.serial;
                    ant_name = serial{a};
                elseif isempty(serial) || isempty(serial{a})
                    cmp_str = this.type;
                    ant_name = type{a};
                else
                    cmp_str = [this.type this.serial];
                    ant_name = [type{a} serial{a}];
                end
                tmp = ant_name;
                ant_name = char(32*ones(1,20));
                ant_name(1:length(tmp)) = tmp;
                if (numel(ant_name) >= 4) && Core_Utils.code4Char2Num(ant_name(1:4)) == none_code
                    % ignoring antenna type == 'NONE'
                    id_ant = [];
                else
                    % Get the antenna id containing the searched antenna file
                    id_ant = ((strfind(serialize(cmp_str')', ant_name) - 1) / size(cmp_str, 2) + 1);
                    id_ant(rem(id_ant,1) > 0) = []; % remove wrong matches
                    if isempty(id_ant) && ~(isempty(type) || isempty(type{a})) && (Core_Utils.code4Char2Num(ant_name(17:20)) ~= none_code) % search for same antenna with no radome
                        % Check with no radome
                        if isempty(serial) || isempty(serial{a})
                            cmp_str = this.type;
                            ant_name = type{a};
                            ant_name(17:20) = 'NONE';
                        else
                            cmp_str = [this.type this.serial];
                            ant_name = [type{a} serial{a}];
                            ant_name(17:20) = 'NONE';
                        end
                        id_ant = ((strfind(serialize(cmp_str')', ant_name) - 1) / size(cmp_str, 2) + 1);
                        id_ant(rem(id_ant,1) > 0) = []; % remove wrong matches
                        
                        id_ant = intersect(id_ant, t_check);
                        if ~isempty(id_ant)
                            if (nargin > 4)
                                log.addWarning(sprintf('"%s" antenna PCO/PCV found but with no radome', ant_name), log_lev);
                            else
                                log.addWarning(sprintf('"%s" antenna PCO/PCV found but with no radome', ant_name));
                            end
                        end
                    else
                        if isempty(id_ant)
                            if (nargin > 4)
                                log.addWarning(sprintf('"%s" antenna PCO/PCV not found', ant_name), log_lev);
                            else
                                log.addWarning(sprintf('"%s" antenna PCO/PCV not found', ant_name));
                            end
                        else
                            id_ant = intersect(id_ant, t_check);
                            if isempty(id_ant)
                                if (nargin > 4)
                                    log.addWarning(sprintf('"%s" antenna PCO/PCV found but it is out of validity range,\n      not using corrections', ant_name), log_lev);
                                else
                                    log.addWarning(sprintf('"%s" antenna PCO/PCV found but it is out of validity range,\n      not using corrections', ant_name));
                                end
                            else
                                % Keep the first occurence of the antenna
                                id_ant = id_ant(1);
                            end
                        end
                    end
                end
                if ~isempty(id_ant)
                    % Keep the antenna with the latest calibration
                    if numel(id_ant) > 1
                        id_keep = 1;
                        cal_date = 0;
                        for i = 1 : numel(id_ant)
                            tmp = this.list(id_ant(i)).cal_date.getMatlabTime;
                            if cal_date < tmp
                                id_keep = i;
                                cal_date = tmp;
                            end
                        end
                        id_ant = id_ant(id_keep);
                    end
                    ant(a) = this.list(id_ant);
                end
            end
        end
        
        function [type, found] = getTypeFromMarker(this, marker_name, time, log_lev)
            % Get antenna type from Monument File given the marker name of the station and the epoch
            %
            % SYNTAX
            %   [type, found] = this.getTypeFromMarker(marker_name, time, log_lev)
            type = '';
            found = false;
            if ~isempty(this.mnt_tbl)
                % Get the last monument for this station
                id_marker = find(Core_Utils.code4Char2Num(this.mnt_tbl.marker) == Core_Utils.code4Char2Num(marker_name));
                if ~isempty(id_marker)
                    id_marker = id_marker(this.mnt_tbl.time.getEpoch(id_marker) < time); % keep calibration older than the observations
                    if ~isempty(id_marker)
                        found = true;
                        [~, id_max] = max(this.mnt_tbl.time.getEpoch(id_marker).getMatlabTime);
                        id_marker = id_marker(id_max);
                        type = this.mnt_tbl.type(id_marker,:);
                    end
                end
                if ~found
                    log = Core.getLogger();
                    if (nargin > 4)
                        log.addWarning(sprintf('Monument not found for receiver "%s"', marker_name), log_lev);
                    else
                        log.addWarning(sprintf('Monument not found for receiver "%s"', marker_name));
                    end
                end
            end            
        end
    end
end
