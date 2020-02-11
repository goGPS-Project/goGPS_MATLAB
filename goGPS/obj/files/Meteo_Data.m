%   CLASS Meteo_Data
% =========================================================================
%
% DESCRIPTION
%   Class to store receiver data (observations, and characteristics
%
% EXAMPLE
%   settings = Meteo_Data();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Meteo_Data
%
% REFERENCE

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
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
%----------------------------------------------------------------------------------------------
classdef Meteo_Data < handle

    properties (Constant, GetAccess = private)
        PR_ID = sum(uint16('PR').*uint16([1 256]));     % Internal id of Pressure [mbar]
        TD_ID = sum(uint16('TD').*uint16([1 256]));     % Internal id of Dry temperature (deg Celsius)
        HR_ID = sum(uint16('HR').*uint16([1 256]));     % Internal id of Relative humidity (percent)
        ZW_ID = sum(uint16('ZW').*uint16([1 256]));     % Internal id of Wet zenith path delay (mm), (for WVR data)
        ZD_ID = sum(uint16('ZD').*uint16([1 256]));     % Internal id of Dry component of zen.path delay (mm)
        ZT_ID = sum(uint16('ZT').*uint16([1 256]));     % Internal id of Total zenith path delay (mm)
        WD_ID = sum(uint16('WD').*uint16([1 256]));     % Internal id of Wind azimuth (deg) from where the wind blows
        WS_ID = sum(uint16('WS').*uint16([1 256]));     % Internal id of Wind speed (m/s)
        RI_ID = sum(uint16('RI').*uint16([1 256]));     % Internal id of Rain increment (1/10 mm): Rain accumulation since last measurement
        RR_ID = sum(uint16('RR').*uint16([1 256]));     % Internal id of Rain Rate (mm/h): Rain rate
        HI_ID = sum(uint16('HI').*uint16([1 256]));     % Internal id of Hail indicator non-zero: Hail detected since last measurement

        % Array of all the meteorological data types
        DATA_TYPE_ID = [Meteo_Data.PR_ID, ...
                        Meteo_Data.TD_ID, ...
                        Meteo_Data.HR_ID, ...
                        Meteo_Data.ZW_ID, ...
                        Meteo_Data.ZD_ID, ...
                        Meteo_Data.ZT_ID, ...
                        Meteo_Data.WD_ID, ...
                        Meteo_Data.WS_ID, ...
                        Meteo_Data.RI_ID, ...
                        Meteo_Data.RR_ID, ...
                        Meteo_Data.HI_ID];
    end

    properties (Constant)
        DATA_TYPE = ['PR'; 'TD'; 'HR'; 'ZW'; 'ZD'; 'ZT'; 'WD'; 'WS'; 'RI'; 'RR'; 'HI'];
        DATA_TYPE_EXT = { 'Pressure [mbar]', ...
                          'Dry temperature [deg Celsius]', ...
                          'Relative humidity [percent]', ...
                          'Wet zenith path delay [mm], (for WVR data)', ...
                          'Dry component of zen.path delay [mm]', ...
                          'Total zenith path delay [mm]', ...
                          'Wind azimuth [deg] from where the wind blows', ...
                          'Wind speed [m/s]', ...
                          'Rain increment [1/10 mm]: Rain accumulation since last measurement', ...
                          'Rain Rate (mm/h): Rain rate', ...
                          'Hail indicator non-zero: Hail detected since last measurement'};

        PR = 1;     % Numeric id of Pressure [mbar]
        TD = 2;     % Numeric id of Dry temperature (deg Celsius)
        HR = 3;     % Numeric id of Relative humidity (percent)
        ZW = 4;     % Numeric id of Wet zenith path delay (mm), (for WVR data)
        ZD = 5;     % Numeric id of Dry component of zen.path delay (mm)
        ZT = 6;     % Numeric id of Total zenith path delay (mm)
        WD = 7;     % Numeric id of Wind azimuth (deg) from where the wind blows
        WS = 8;     % Numeric id of Wind speed (m/s)
        RI = 9;     % Numeric id of Rain increment (1/10 mm): Rain accumulation since last measurement
        RR = 10;     % Numeric id of Rain Rate (mm/h): Rain rate'
        HI = 11;    % Numeric id of Hail indicator non-zero: Hail detected since last measurement
    end

    properties (SetAccess = protected, GetAccess = protected)
        log = Logger.getInstance(); % Handler to the log object
    end

    properties (SetAccess = private, GetAccess = public)
        % contains an object to read the RINEX file
        file;   % init this with File_Rinex('filename')
        rin_type       % rinex version format

        marker_name = '';   % name of the station
        n_type = 0;         % number of observation types
        type = 1:11;        % supposing to have all the fields

        time = GPS_Time();  % array of observation epochs
        rate                % observations rate;

        data = [];          % Meteorological file

        xyz = [0 0 0];      % geocentric coordinate of the sensor
        amsl = 0;           % hortometric height of the sensor
        is_valid = false;   % Status of valitity of the file;
        
        max_bound = 90;     % Max bound to extrapolate
        
        smoothing = [900 300 300 0 0 0 0 0 0 0 0]; % Spline base for smoothing (in seconds)
    end

    methods (Access = private)
        function parseRinHead(this, txt, lim, eoh)
            % Parse header and update the object, having as input the meteo_file in txt
            %
            % SYNTAX
            %    this.parseRinHead(txt, nl)
            %
            % INPUT
            %    txt    raw txt of the RINEX
            %    lim    indexes to determine start-stop of a line in "txt"  [n_line x 2/<3>]
            %    eoh    end of header line
            
            h_std{1} = 'RINEX VERSION / TYPE';                  %  1
            h_std{2} = 'MARKER NAME';                           %  2
            h_std{3} = 'SENSOR MOD/TYPE/ACC';                   %  4
            h_std{4} = 'SENSOR POS XYZ/H';                      %  5
            h_std{5} = '# / TYPES OF OBSERV';                   %  6
            
            h_opt{1} = 'MARKER NUMBER';                         %  7
            head_field = {h_std{:} h_opt{:}}';
            
            try
                % read RINEX type 3 or 2 ---------------------------------------------------------------------------------------------------------------------------
                l = 0;
                type_found = false;
                while ~type_found && l < eoh
                    l = l + 1;
                    if strcmp(strtrim(txt((lim(l,1) + 60) : lim(l,2))), h_std{1})
                        type_found = true;
                        dataset = textscan(txt(lim(1,1):lim(1,2)), '%f%c%18c%c');
                    end
                end
                this.rin_type = dataset{1};
                if dataset{2} == 'M'
                    % OK
                else
                    throw(MException('VerifyINPUTInvalidObservationFile', 'This RINEX seems corrupted - RINEX header parsing error'));
                end
                
                % parsing ------------------------------------------------------------------------------------------------------------------------------------------
                
                % retriving the kind of header information is contained on each line
                line2head = zeros(eoh, 1);
                l = 0;
                while l < eoh
                    l = l + 1;
                    %DEBUG: txt((lim(l,1) + 60) : lim(l,2))
                    tmp = find(strcmp(strtrim(txt((lim(l,1) + 60) : lim(l,2))), head_field));
                    if ~isempty(tmp)
                        % if the field have been recognized (it's not a comment)
                        line2head(l) = tmp;
                    end
                end
                
                % reading parameters -------------------------------------------------------------------------------------------------------------------------------
                
                % 1) 'RINEX VERSION / TYPE'
                % already parsed
                % 2) 'MARKER NAME'
                fln = find(line2head == 2, 1, 'first'); % get field line
                if isempty(fln)
                    this.marker_name = 'NO_NAME';
                else
                    this.marker_name = strtrim(txt(lim(fln, 1) + (0:59)));
                end
                % 3) 'SENSOR MOD/TYPE/ACC'
                % ignoring
                % 4) 'SENSOR POS XYZ/H'
                fln = find(line2head == 4, 1, 'first'); % get field line
                if isempty(fln)
                    this.xyz = [0 0 0];
                else
                    tmp = sscanf(txt(lim(fln, 1) + (0:41)),'%f')';                                               % read value
                    this.xyz = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 3), [0 0 0], tmp);          % check value integrity
                    if sum(abs(this.xyz)) > 0
                        [~, lam, h, phiC] = cart2geod(this.xyz(1), this.xyz(2), this.xyz(3));
                        this.amsl = h - getOrthometricCorr(phiC, lam);
                    else
                        this.amsl = 0;
                    end
                end
                % 5) '# / TYPES OF OBSERV'
                fln = find(line2head == 5, 1, 'first'); % get field line
                % Read types of obs
                if isempty(fln)
                    throw(MException('MeteorologicalFile:InvalidHeader', 'unrecognized data type'));
                else
                    this.n_type = sscanf(txt(lim(fln, 1) + (0 : 5)),'%d');
                    if (this.n_type >= 10)
                        this.log.addWarning('Reading more than 10 fields is not yet supported');
                    end
                    this.type = [];
                    for t = 0 : (this.n_type - 1)
                        str_type = txt(lim(fln, 1) + 7 + (3 : 4) + t * 6);
                        this.type = [this.type find(this.DATA_TYPE_ID == (sum(uint16(str_type) .* uint16([1 256]))))];
                    end
                    if this.n_type > numel(this.type)
                        throw(MException('MeteorologicalFile:InvalidHeader', 'unrecognized data type'));
                    end
                    if isempty(this.type)
                        throw(MException('MeteorologicalFile:InvalidHeader', 'file does not contain the header field "# / TYPES OF OBSERV"'));
                    end
                end
            catch ex
                this.log.addWarning(sprintf('Problem detected in the header: %s', ex.message))
                this.log.addWarning(sprintf('try to read it as it contains date (6 fields) + %s', sprintf('%c%c ', this.getType()')));
            end
        end
        
        function parseRin2Data(this, txt, lim, eoh)
            % Parse the data part of a RINEX 2 file -  the header must already be parsed
            %
            % SYNTAX
            %   this.parseRin2Data(txt, lim, eoh)
            
            % Read the data
            try
                % search for lines containing comments (to be ignored)
                comment_line = sum(txt(repmat(lim(1:end-2,1),1,7) + repmat(60:66, size(lim,1)-2, 1)) == repmat('COMMENT', size(lim,1)-2, 1),2) == 7;
                comment_line(1 : eoh) = false;
                lim(comment_line, :) = [];
                
                % find all the observation lines
                t_line = find([false(eoh, 1); (txt(lim(eoh+1:end,1) + 2) ~= ' ')' & (txt(lim(eoh+1:end,1) + 3) == ' ')' & lim(eoh+1:end,3) > 25]);
                t_line(lim(t_line,3) ~= lim(t_line(1),3)) = [];
                
                n_epo = numel(t_line);
                % extract all the epoch lines
                string_time = txt(repmat(lim(t_line,1),1,18) + repmat(1:18, n_epo, 1))';
                % convert the times into a 6 col time
                date = cell2mat(textscan(string_time,'%f %f %f %f %f %f'));
                after_70 = (date(:,1) < 70); date(:, 1) = date(:, 1) + 1900 + after_70 * 100; % convert to 4 digits
                % import it as a GPS_Time obj
                this.time = GPS_Time(date, [], this.file.first_epoch.is_gps);
                this.rate = this.time.getRate();
                
                this.data = nan(this.n_type, n_epo);
                line = txt(repmat(lim(t_line,1),1, lim(t_line(1),3) - 17) + 17 + repmat(1 : (lim(t_line,3) - 17), n_epo, 1))';
                data = sscanf(line, '%f');
                this.data(1:end) = data;
            catch ex
                this.log.addWarning(sprintf('Problem detected while reading meteorological data: %s', ex.message));
                this.is_valid = false;
            end
            this.data = this.data'; % keep one column per data type
        end
        
        function parseRin3Data(this, txt, lim, eoh)
            % Parse the data part of a RINEX 2 file -  the header must already be parsed
            %
            % SYNTAX 
            %   this.parseRin2Data(txt, lim, eoh)
            
            % Read the data
            try
                % search for lines containing comments (to be ignored)
                comment_line = sum(txt(repmat(lim(1:end-2,1),1,7) + repmat(60:66, size(lim,1)-2, 1)) == repmat('COMMENT', size(lim,1)-2, 1),2) == 7;
                comment_line(1 : eoh) = false;
                lim(comment_line, :) = [];
                
                % find all the observation lines
                t_line = find([false(eoh, 1); (txt(lim(eoh+1:end,1) + 3) ~= ' ')' & (txt(lim(eoh+1:end,1) + 5) == ' ')' & lim(eoh+1:end,3) > 25]);
                t_line(lim(t_line,3) ~= lim(t_line(1),3)) = [];
                
                n_epo = numel(t_line);
                % extract all the epoch lines
                string_time = txt(repmat(lim(t_line,1),1,19) + repmat(1:19, n_epo, 1))';
                % convert the times into a 6 col time
                date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %2f'));
                % import it as a GPS_Time obj
                this.time = GPS_Time(date, [], this.file.first_epoch.is_gps);
                this.rate = this.time.getRate();
                
                this.data = nan(this.n_type, n_epo);
                line = txt(repmat(lim(t_line,1),1, lim(t_line(1),3) - 19) + 19 + repmat(1 : (lim(t_line,3) - 19), n_epo, 1))';
                data = sscanf(line, '%f');
                this.data(1:end) = data;
            catch ex
                this.log.addWarning(sprintf('Problem detected while reading meteorological data: %s', ex.message));
                this.is_valid = false;
            end
            this.data = this.data'; % keep one column per data type
        end
        
        function parseData(this, meteo_file)
            % Parse data and update the object, having as input the meteo_file
            % as read with textscan (cell array of string lines)
            %
            % SYNTAX
            %  this.parseData(meteo_file)

            % Read the data
            try
                eoh = this.file.getEOH();
                if (numel(meteo_file) < eoh)
                    eoh = 0; % file with no header
                end
                n_epoch = numel(meteo_file) - eoh;

                % Guess the rate of the data (in seconds)
                % rate = round((this.file.last_epoch.getMatlabTime-this.file.first_epoch.getMatlabTime)*86400) / n_epochs;

                this.data = nan(this.n_type, n_epoch);
                % try to guess the time format
                [id_start, id_stop] = regexp(meteo_file{eoh+1}, '[.0-9]*');
                id_date = id_start(1) : id_stop(6); % save first and last char limits of the date in the line -> suppose it composed by 6 fields

                for l = (eoh + 1) : numel(meteo_file)
                    line = meteo_file{l};
                    this.time.addEpoch(line(id_date), [], true);
                    [value, ~, id_stop] = regexp(line((id_date(end) + 1):end), '[.0-9]*', 'match');
                    col_id = round(id_stop / 7);
                    this.data(col_id, l - eoh) = str2double(value);
                end
            catch ex
                this.log.addWarning(sprintf('Problem detected while reading meteorological data: %s', ex.message));
                this.is_valid = false;
            end
            this.data = this.data'; % keep one column per data type
        end

        function init(this, file_name, type, verbosity_lev)
            % Try to read the file
            %
            % SYNTAX
            %   this.init(file_name, type, verbosity_lev)
            if (nargin < 4)
                verbosity_lev = Logger.DEFAULT_VERBOSITY_LEV;
            end

            this.log.addMessage(sprintf('Loading Meteorological data from "%s"', file_name), iif(verbosity_lev < 50, Logger.DEFAULT_VERBOSITY_LEV, verbosity_lev));
            try
                this.time = GPS_Time(); % empty the time
                this.file = File_Rinex(file_name, verbosity_lev);
                
                fid = fopen(file_name,'r');
                txt = fread(fid,'*char')';
                txt(txt == 13) = []; % remove carriage return - I hate you Bill!
                fclose(fid);
                
                % get new line separators
                nl = regexp(txt, '\n')';
                if nl(end) <  numel(txt)
                    nl = [nl; numel(txt)];
                end
                lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
                lim = [lim lim(:,2) - lim(:,1)];
                while lim(end,3) < 3
                    lim(end,:) = [];
                end
                
                % removing empty lines at end of file
                while (lim(end,1) - lim(end-1,1))  < 2
                    lim(end,:) = [];
                end
                
                % importing header informations
                eoh = this.file.eoh;
            catch ex
                this.log.addError(sprintf('Error reading meteorological file "%s" (%s)', file_name, ex.message));
                Core_Utils.printEx(ex);
                return
            end

            % Parse the header and detect the types of data contained in the meteorological file
            this.parseRinHead(txt, lim, eoh);

            if (nargin == 3) && ~isempty(type)
                this.log.addWarning('Meteorological file - overriding the file data types with custom types');
                this.type = type;
                this.n_type = numel(type);
                assert(sum(type < 1) + sum(type > 10) == 0, 'Invalid custom types');
            end

            this.log.addMessage(sprintf('\nThe following meteorological data are present in the header:'), verbosity_lev);
            for t = 1 : this.n_type
                this.log.addMessage([' - ' this.getTypeExt{t}], verbosity_lev);
            end
            this.is_valid = this.file.isValid();
            % Parse the data
            if this.rin_type < 3
                this.parseRin2Data(txt, lim, eoh);
            else
                this.parseRin3Data(txt, lim, eoh);
            end
            if ~any(this.xyz)
                this.log.addWarning(sprintf('No position found in meteorological file "%s"\n this meteorological station cannot be used correctly', File_Name_Processor.getFileName(file_name)), verbosity_lev);                
            end
        end
    end
    
    methods (Access = public)
        
        function inject(this, md)
            % inject a  mete file into the other meteo file
            %
            % SYNTAX : this.inject(md
            
            % NOTE: xyz are kept the ones of the first meteo file, it does not manage overlapping data
            % neither object containign different datat types
            
            if isequal(this.type, md.type) 
            [this.time, idx1, idx2] = this.time.injectBatch(md.time);
            
            this.data = Core_Utils.injectData(this.data, md.data, idx1, idx2);
            else
                this.log.addWarning('Meteo files contains different data types, inject skipped')
            end
        end
    end
    
    % =========================================================================
    %  SETTER
    % =========================================================================
    methods
        function setMaxBound(this, max_bound)
            % Set the maximum extrapolation span
            %
            % SYNTAX
            %   this.setMaxBound(max_bound)
            this.max_bound = max_bound;
        end      
    end

    % =========================================================================
    %  INIT / READER
    % =========================================================================
    methods
        function this = Meteo_Data(file_name, type, verbosity_lev)
            % Creator Meteo_Data(file_name, <type = empty>, <verbosity_lev>)

            % The function calls all its creation methods within try and
            % catch statements, reading the Meteo file should not be
            % blocking for the processing, even if the data are not good
            % SYNTAX
            %   this = Meteo_Data(file_name, type, verbosity_lev)
            switch nargin
                case 1
                    this.init(file_name);
                case 2
                    this.init(file_name, type);
                case 3
                    this.init(file_name, type, verbosity_lev);
            end
        end
    end

    % =========================================================================
    %  IMPORT / EXPORT / TOSTRING
    % =========================================================================
    methods
        function importRaw(this, obs_time, data, type, marker_name, pos_xyz)
            % Import a meteorological file
            % EXAMPLE: this.importRaw(GPS_Time(time - 1/12), [pres temp hum rain], [Meteo_Data.PR Meteo_Data.TD Meteo_Data.HR Meteo_Data.RT], 'GReD', xyz);
            
            % SYNTAX
            %   this.importRaw(obs_time, data, type, marker_name, pos_xyz)
            narginchk(6, 6);

            % Skip NaN epochs
            this.marker_name = marker_name;
            this.xyz = pos_xyz;
            [~, lam, h, phiC] = cart2geod(this.xyz(1), this.xyz(2), this.xyz(3));
            this.amsl = h - getOrthometricCorr(phiC, lam);
            ok = ~obs_time.isnan();
            this.data = data(ok, :);
            this.type = type;

            % Cut empty epochs
            invalid_epoch = sum(isnan(this.data),2) == numel(this.type);
            ok(ok) = ok(ok) & ~invalid_epoch;
            this.time = obs_time.getEpoch(ok);
            this.data(invalid_epoch, :) = [];

            % Cut empty data types
            if ~isempty(this.data)
                invalid_data = sum(isnan(this.data)) == size(this.data, 1);
                this.data(:, invalid_data) = [];
                this.type(:, invalid_data) = [];
                this.n_type = numel(this.type);
                this.is_valid = true;
            end
        end

        function import(this, file_name, type)
            % import a meteorological file
            
            % SYNTAX
            %   this.import(file_name, type)
            narginchk(2,3);

            switch nargin
                case 1
                    this.init(file_name);
                case 2
                    this.init(file_name, type);
            end
        end

        function export(this, file_name)
            % Export data to a meteorological RINEX
            
            % SYNTAX
            %   this.export(file_name)
            
            narginchk(1,2);
            if this.time.isempty()
                this.log.addError(sprintf('Export failed - missing data - %s', this.marker_name));
            else
                fnp = File_Name_Processor;
                
                % Find the time span of the observations
                [yyyy, doy] =  this.time.getDOY;
                [~, day_start, day_id] = unique(yyyy*1e4+doy);
                
                if (nargin == 1)
                    state = Core.getCurrentSettings();
                    file_name =  this.marker_name;
                    % generate short 4 letters name
                    if numel(file_name) < 4
                        file_name = sprintf(['%0' num2str(4-numel(file_name)) 'd%s'], 0, file_name);
                    else
                        file_name = upper(file_name(1:4));
                    end
                    file_name = [state.getMetDir() filesep '${YYYY}_${DOY}' filesep file_name '_${DOY}0.${YY}m'];
                end
                
                for d = 1 : numel(day_start)
                    id = day_id == d;
                    cur_file_name = fnp.dateKeyRep(fnp.checkPath(file_name), this.time.getEpoch(day_start(d)));
                    
                    dir_container = fileparts(cur_file_name);
                    if ~isempty(dir_container) && ~exist(dir_container, 'dir')
                        mkdir(dir_container);
                    end
                    
                    this.log.addMessage(sprintf('Exporting met data to %s', cur_file_name));
                    try
                        fid = fopen(cur_file_name, 'Wb');
                        str = ['     3.03           METEOROLOGICAL DATA                     RINEX VERSION / TYPE', 10 ...
                            'EXPORTED MET FILE FROM METEO_DATA MATLAB CLASS              COMMENT', 10];
                        if ~isempty(this.marker_name)
                            str = sprintf(['%s%s%' num2str(59-numel(this.marker_name)) 's MARKER_NAME\n'], str, this.marker_name, '');
                        end
                        
                        line = sprintf('%6d', this.n_type);
                        for t = 1 : this.n_type
                            line = sprintf('%s%6s', line, this.DATA_TYPE(this.type(t), :));
                        end
                        str = sprintf(['%s%s%' num2str(60 - (this.n_type + 1) * 6) 's# / TYPES OF OBSERV\n'], str, line, '');
                        str = sprintf('%s%14.4f%14.4f%14.4f%14.4f PR SENSOR POS XYZ/H\n', str, this.xyz, this.amsl);
                        str = [str '                                                            END OF HEADER' 10]; %#ok<*AGROW>
                        fwrite(fid, str);
                        epochs = this.time.getEpoch(id).toString(' yyyy mm dd HH MM SS ')';
                        str = [epochs; reshape(sprintf('%7.1f', this.data(id,:)'), 7 * size(this.data(id,:),2), size(this.data(id,:),1)); 10 * ones(1, size(this.data(id,:),1))];
                        fwrite(fid, str);
                        fclose(fid);
                    catch ex
                        this.log.addError(sprintf('Export failed - %s', ex.message));
                    end
                end
            end
        end

        function str = toString(this, str)
            % Display the loaded meteorological data
            
            % SYNTAX
            %   this.toString(str)

            if (nargin == 1)
                str = '';
            end

            str = [str '---- METEOROLOGICAL DATA -------------------------------------------------' 10 10];
            if ~isempty(this.marker_name)
                str = [str sprintf(' Station %s\n\n', this.marker_name)];
            end

            str = [str sprintf(' Location (XYZ)  %f %f %f\n', this.xyz)];
            str = [str sprintf('          (amsl) %f\n\n', this.amsl)];

            str = [str sprintf(' Data available from %s\n                  to %s\n\n', this.time.first.toString('dd mmm yyyy HH:MM:SS'), this.time.last.toString('dd mmm yyyy HH:MM:SS'))];
            str = [str sprintf(' The following meteorological data are present:\n')];
            type_ext = this.getTypeExt();
            for t = 1 : this.n_type
                str = [str sprintf('  - %s\n', type_ext{t})]; %#ok<AGROW>
            end
            str = [str 10];
        end
    end

    % =========================================================================
    %  GETTERS
    % =========================================================================
    
    methods
        function is_empty = isEmpty(md_list)
            % Return a logical array contain which station is 
            is_empty = false(numel(md_list), 1);
            for i = 1 : numel(md_list)
                is_empty(i) = isempty(md_list(i)) | isempty(md_list(i).data);
            end
        end
        
        function validity = isValid(this)
            % Get the validity of a RINEX file or the object
            
            % SYNTAX
            %   validity = isValid()
            validity = this.file.isValid();
        end
                
        function name = getMarkerName(this)
            % Get the name of the station
            %
            % SYNTAX
            %   time = this.getMarkerName()
            name = this.marker_name;
        end
        
        function time = getTime(this)
            % Get the epochs of the data
            %
            % SYNTAX
            %   time = this.getTime()
            time = this.time;
        end
        
        function type = getType(this)
            % Get the types of data stored in the RINEX
            %
            % SYNTAX
            %   
            id = this.getTypeId();
            type = this.DATA_TYPE(id,:);
        end
        
        function max_bound = getMaxBound(this)
            % Get the maximum extrapolation
            %
            % SYNTAX
            %   max_bound = getMaxBound(this)
            max_bound = this.max_bound;
        end
        
        function type = getTypeExt(this)
            % Get the description of the types of data stored in the RINEX
            %
            % SYNTAX
            %   type = getTypeExt(this)
            id = this.getTypeId();
            type = this.DATA_TYPE_EXT(id);
        end
        
        function id = getTypeId(this)
            % Get the id of the types of data stored in the RINEX
            %
            % SYNTAX
            %   id = getTypeId(this)
            id = this.type;
        end
        
        function data = getComponent(this, data_id, time)
            % Get the data with id of the type wanted
            % Passing a time array as GPS_Time the object interpolate the
            % data contained in the meteorological file
            % SYNTAX
            %   data = this.getComponent(id, <time>)
            id = find(this.type == data_id);
            if isempty(id)
                if nargin == 3
                    data = nan(time.length(), 1);
                else
                    data = [];
                end
            else
                if size(this.data, 2) >= id
                    data_in = this.data(:,id);
                    if nargin == 3
                        data = nan(time.length(), 1);
                        time_data = this.time.getMatlabTime();
                        time_pred = time.getMatlabTime();
                        if (sum(~isnan(data_in)) > 0)
                            if numel(data_in(~isnan(data_in))) == 1
                                [~, id] = min(time_data(~isnan(data_in)) - time_pred);
                                data(id) = data_in(~isnan(data_in));
                            else
                                time_data = time_data(~isnan(data_in));
                                data_in = data_in(~isnan(data_in));
                                data_in = [data_in(1); data_in; data_in(end)];
                                time_data = [ (min(time_pred(1), time_data(1)) - 1/86400); time_data; (max(time_pred(end), time_data(end)) + 1/86400)];
                                data = interp1(time_data, data_in, time_pred, 'pchip','extrap');
                                if this.smoothing(data_id) > 0
                                    data = splinerMat(time_pred * 86400, data - mean(data), (this.smoothing(data_id)), 0) + mean(data);
                                end
                                % do not extrapolate further than 20 minutes in time
                                data((time_pred < time_data(1) - this.getMaxBound / 1440) | (time_pred > time_data(end) + this.getMaxBound / 1440)) = NaN;
                                % extrapoleted value
                            end
                        end
                    else
                        data = data_in;
                    end
                else
                    data = [];
                end
            end
        end
        
        function t_dist = getTimeInterpDistance(this, time, time2)
            % Get the time distance in seconds from the observed data
            %
            % SYNTAX
            %   t_dist = getTimeInterpDistance(this, time)
            
            t_dist = nan(time.length, 1);
            i = 1;
            j = 1;
            % moving on prediction time
            t0_mat = time.getMatlabTime;
            if nargin == 3
                t1_mat = time2.getMatlabTime();
            else
                t1_mat = this.time.getMatlabTime();
            end
            while (i <= time.length())
                t0 = t0_mat(i);
                % moving on data time
                while j < length(t1_mat) && (abs(t0 - t1_mat(j)) > (abs(t0 - t1_mat(j + 1)))); j = j + 1; end
                t_dist(i) = abs(t0 - t1_mat(j)) * 86400;
                i = i + 1;
            end
        end
        
        function data = getPressure(this, time, amsl)
            % Get the pressure data
            % SYNTAX
            %   data = this.getPressure()
            if (nargin == 1)
                data = this.getComponent(1);
            else
                data = this.getComponent(1, time);
            end
            
            if (nargin == 3)
                data = Meteo_Data.pressureAdjustment(data, this.amsl, amsl);
            end
        end
        
        function data = getTemperature(this, time, amsl)
            % Get the temperature data
            % SYNTAX
            %   data = this.getTemperature()
            if (nargin == 1)
                data = this.getComponent(2);
            else
                data = this.getComponent(2, time);
            end
            
            if (nargin == 3)
                data = Meteo_Data.temperatureAdjustment(data, this.amsl, amsl);
            end
        end
        
        function data = getHumidity(this, time, amsl)
            % Get the humidity data
            % SYNTAX
            %   data = this.getHumidity()
            if (nargin == 1)
                data = this.getComponent(3);
            else
                data = this.getComponent(3, time);
            end
            
            if (nargin == 3)
                data = Meteo_Data.humidityAdjustment(data, this.amsl, amsl);
            end
        end
        
        function [x, y, z, amsl] = getLocation(md_list)
            % Get meteo station location
            % SINTAX
            %   [x, y, z, amsl] = this.getLocation();
            if nargout == 1
                x = nan(numel(md_list), 3);
                for i = 1 : numel(md_list)
                    x(i, :) = md_list(i).xyz;
                end
            else                
                x = nan(numel(md_list), 1);
                y = nan(numel(md_list), 1);
                z = nan(numel(md_list), 1);
                amsl = nan(numel(md_list), 1);
                for i = 1 : numel(md_list)
                    x(i) = md_list(i).xyz(1);
                    y(i) = md_list(i).xyz(2);
                    z(i) = md_list(i).xyz(3);
                    amsl(i) = md_list(i).amsl;
                end
            end                
        end
        
        function time = getObsTime(this)
            % Get meteo station observatipon time
            % SINTAX
            %   time = this.getObsTime();
            time = this.time;
        end
    end

    % =========================================================================
    %  STATIC
    % =========================================================================
    methods (Static)
        function data = getMeteoData(station, st_type, fun, time, amsl, d_oo, d_op, type)
            % Get interpolated data of meteo parameter
            % pressure / temperaure / humidity
            %
            % SYNTAX
            %   data = Meteo_Data.getMeteoData(station, st_type, fun, time, amsl, d_oo, d_op, type)
            
            t_thr = 3000;
            
            q_fun_obs = fun(d_oo) .* repmat(fun(d_op)', size(d_oo,1), 1);
            q_fun_obs = triu(q_fun_obs) + triu(q_fun_obs,1)';
            
            % getting data
            id_data = find(st_type(:, type));
            % id_pr = id_pr(real_dist(id_pr) < 20e3);
            
            data_obs = zeros(numel(id_data), time.length());
            switch type
                case Meteo_Data.PR
                    for s = 1 : numel(id_data)
                        data_obs(s, :) = station(id_data(s)).getPressure(time, amsl);
                    end
                case Meteo_Data.TD
                    for s = 1 : numel(id_data)
                        data_obs(s, :) = station(id_data(s)).getTemperature(time, amsl);
                    end
                case Meteo_Data.HR
                    for s = 1 : numel(id_data)
                        data_obs(s, :) = station(id_data(s)).getHumidity(time, amsl);
                    end
            end
            
            id_data(sum(isnan(data_obs),2) > 0) = [];
            data_obs(sum(isnan(data_obs),2) > 0, :) = [];
            
            % get the time distance from true observations
            t_dist = zeros(numel(id_data), time.length());
            for s = 1 : numel(id_data)
                t_dist(s, :) = station(id_data(s)).getTimeInterpDistance(time, station(id_data(s)).time.getEpoch(find(~isnan(station(id_data(s)).getPressure()))));
            end
            
            if isempty(id_data)
                log =  Core.getLogger();
                switch type
                    case Meteo_Data.PR
                        log.addWarning('There are no station to get pressure information', 100);
                    case Meteo_Data.TD
                        log.addWarning('There are no station to get temperature information', 100);
                    case Meteo_Data.HR
                        log.addWarning('There are no station to get humidity information', 100);
                end
                data = nan(time.length,1);
            else
                % A = ones(size(id_pr));
                % Q = d2(id_pr, id_pr);
                % AinvQ =  A'/Q;
                % w = (AinvQ*A)\AinvQ;
                
                trans = sum(q_fun_obs(id_data, id_data));
                w = sum(trans)\trans;
                
                % Rescale weigths epoch by epoch
                w_all = repmat(w', 1, size(t_dist,2));
                w_all(t_dist > t_thr) = 0; % eliminate interpolations too much distant in time
                lid_best = (sum(w_all > 0.8)) >= 1;
                if sum(lid_best) < 2
                    [~, id_min] = min(d_op + 1e300 * double(~(sum(w_all, 2) > 0)));
                    lid_best = w_all(id_min, :) > 0;
                end
                w_all = bsxfun(@rdivide, w_all, sum(w_all,1));
                data0 = sum(w_all .* data_obs, 1);
                
                if ~any(data0(:)) % there are no valid data
                    data = [];
                else
                    lim = getOutliers(lid_best);
                    
                    % adjust pres0 and avoid discontinuities
                    % temporary approach
                    ddata = Core_Utils.diffAndPred(data0'); sensor = abs(ddata - medfilt_mat(ddata, 3));
                    id_jmp = sensor > 1e-3;
                    ddata_fill = simpleFill1D(ddata, id_jmp, 'linear');
                    data_fill = cumsum(ddata_fill);
                    data_fill(lid_best) = data0(lid_best);
                    
                    % first block bias
                    if ~isempty(lim) && (lim(1) > 1)
                        data_fill(1 : lim(1) - 1)  = data_fill(1 : lim(1) - 1) + data0(lim(1)) - data_fill(lim(1) - 1) - ddata_fill(lim(1) - 1);
                    end
                    
                    % middle blocks linear interpolations
                    for l = 2 : size(lim, 1)
                        m = (data_fill(lim(l,1)) - ddata_fill(lim(l,1) + 1) - data_fill(lim(l,1) - 1) ...
                            - ( data_fill(lim(l-1,2)) + ddata_fill(lim(l-1,2) + 1) - data_fill(lim(l-1,2) + 1))) / ...
                            (lim(l,1) - lim(l-1,2) + 1);
                        data_fill((lim(l-1, 2) + 1) : (lim(l, 1) - 1)) = data_fill((lim(l-1, 2) + 1) :  (lim(l, 1) - 1)) + ...
                            m * (0 : (lim(l,1) - 2 - lim(l-1, 2)))' + ...
                            ( data_fill(lim(l-1,2)) + ddata_fill(lim(l-1,2) + 1) - data_fill(lim(l-1,2) + 1));
                    end
                    
                    % last block bias
                    
                    if ~isempty(lim) && (lim(end) < size(data_fill,1))
                        data_fill((lim(end) + 1) : end) = data_fill((lim(end) + 1) : end) - data_fill(lim(end) + 1) + data_fill(lim(end)) + ddata_fill(lim(end) + 1);
                    end
                    
                    data = data_fill;
                    
                    if type == Meteo_Data.HR
                        inan = isnan(data);
                        data = min(100, max(0, data));
                        data(inan) = nan;
                    end
                end
            end
        end
                
        function md = getVMS(name, xyz, time, station)
            % Get Virtual Meteo Station
            %
            % INPUT
            %   name    name of the new Meteo_Data virtual station
            %   xyz     coordinates of the new meteo station
            %   time    time of interpolation
            %
            %   station list of Meteo_Data station to use for the interpolation
            %   
            % OUTPUT
            %   md      virtual Meteo_Data station generated at xyz coordinates
            %   
            % SYNTAX
            %   md = this.getVMS(name, xyz, time, station)
            %
            % EXAMPLE
            %   [x, y, z, amsl] = station(1).getLocation();
            %   md1 = Meteo_Data.getVMS('test', [x y z], station(1).getObsTime, md)

            md = Meteo_Data();
            [~, lam, h, phiC] = cart2geod(xyz(1), xyz(2), xyz(3));
            [e, n] = cart2plan(xyz(1), xyz(2), xyz(3));

            amsl = h - getOrthometricCorr(phiC, lam);

            n_station = numel(station);

            % In a VMS I keep only PR TD HR
            st_type = false(3, numel(station));
            e_obs = zeros(numel(station), 1);
            n_obs = zeros(numel(station), 1);
            for s = 1 : n_station
                % Get the supported type by the input stations
                type = station(s).getTypeId;
                st_type(type(type <= 3), s) = true;

                % Get the stations location
                [x, y, z] = station(s).getLocation();
                [e_obs(s), n_obs(s)] = cart2plan(x, y, z);
            end
            st_type = st_type';

            [e_mesh, n_mesh] = meshgrid(e_obs, n_obs);
            d_oo = sqrt(abs(e_mesh - e_mesh').^2 + abs(n_mesh - n_mesh').^2); % distance obs obs
            d_op = sqrt(abs(e_obs - e).^2 + abs(n_obs - n).^2);               % distance obs o prediction point

            % fun for pressure
            fun = @(dist) 0.2 * exp(-(dist/0.8e4)) + exp(-(dist/6e3).^2);
            %fun = @(dist) min(1, 1 ./ dist.);
            pres = Meteo_Data.getMeteoData(station, st_type, fun, time, amsl, d_oo, d_op, Meteo_Data.PR);
                        
            % fun for temperature
            fun = @(dist) 0.2 * exp(-(dist/1e4)) + exp(-(dist/6e3).^2);
            %fun = @(dist) min(1, 1 ./ dist);
            temp = Meteo_Data.getMeteoData(station, st_type, fun, time, amsl, d_oo, d_op, Meteo_Data.TD);
            
            % fun for humidity
            fun = @(dist) exp(-(dist/1e4)) + exp(-(dist/8e3).^2);
            % fun = @(dist) min(1, 1 ./ dist);
            hum = Meteo_Data.getMeteoData(station, st_type, fun, time, amsl, d_oo, d_op, Meteo_Data.HR);
            
            data = [pres temp hum];
            if sum(any(data))
                id_ok = (sum(isnan(data)) < time.length());
                data = data(:, id_ok);
                type = [Meteo_Data.PR Meteo_Data.TD Meteo_Data.HR];
                type = type(:, id_ok);
                md.importRaw(time, data, type, name, xyz);
            end
        end

        function [ temperature_adj ] = temperatureAdjustment( temperature , obs_h, pred_h)
            % Barometric formula taken from Bai and Feng, 2003.
            % The parameter value is taken from Realini et al., 2014
            % Parameter definition
            %
            % SYNTAX
            %   [ temperature_adj ] = temperatureAdjustment( temperature , obs_h, pred_h)
            grad = 0.0065 ; % * C / m gravitational acceleration constant
            temperature_adj = temperature + grad * (obs_h - pred_h);
        end

        function [ pressure_adj ] = pressureAdjustment( pressure , obs_h, pred_h)
            % Barometric formula taken from Berberan-Santos et al., 1997
            % The parameter values are taken from Realini et al., 2014
            % Parameters definition
            %
            % SYNTAX
            % [ pressure_adj ] = pressureAdjustment( pressure , obs_h, pred_h)
            g = 9.80665 ;    % m / s^2 gravitational acceleration constant
            Md = 0.0289644 ; % kg / mol molar mass of dry air
            R = 8.31432 ;    % J(mol * K) gas constant for air
            Tisa = 288.15 ;  % K international standard temperature of the atmosphere at the sea level
            pressure_adj = pressure .* exp(-(g * Md * (pred_h - obs_h))/(R * Tisa ));
        end

        function [ humidity_adj ] = humidityAdjustment( humidity , obs_h, pred_h)
            % Empirically derived from here, to be substituted with something better...
            % % http://www.engineeringtoolbox.com/relative-humidity-air-d_687.html
            %
            % tmp = [0 1.000; 108 0.987; 200 0.976; 400 0.953; 600 0.931; 800 0.909; 1000 0.887; 1500 0.835; 2000 0.785];
            %
            % % quadratics fitting
            % A = [tmp(:, 1).^2 tmp(:, 1) ones(size(tmp,1),1)];
            % y0 = tmp(:,2);
            % x = (A' * A) \ A' * y0
            % figure; plot(tmp(:,1),tmp(:,2),'o'); hold on; amsl = -10 : 3000; plot(amsl, x(1) .* amsl.^2 + x(2) .* amsl + x(3) .* ones(size(amsl))); setAllLinesWidth(2)
            %x = [ 5.25984524194874e-09; -0.000117788989855394; 0.999649177981675 ];
            %humidity_adj = humidity * ([pred_h^2 pred_h 1] * x) / ([obs_h^2 obs_h 1] * x);
            %
            % SYNTAX
            %   [ humidity_adj ] = humidityAdjustment( humidity , obs_h, pred_h)
            humidity_adj = humidity / exp(-6.396e-4 * (obs_h - pred_h));
        end
        
        function importSTEAM(data_path)
            % Import data from STEAM and create MET files
            %
            % SYNTAX
            %   Meteo_Data.importSteam(data_path)            
            if nargin == 0
                data_path = '/Volumes/Data/goGPS_data/project/STEAM_Pilot/Stations/MET';
            end
            
            do_export = true;
            show_fig = true;
            
            state = Core.getCurrentSettings();
            state.setMetDir(data_path);
            log = Logger.getInstance();
            v_lev = log.getVerbosityLev();
            
            pres = load([data_path filesep 'pressione.mat']);
            pres.a2dOsservazioni(pres.a2dOsservazioni < -9000) = NaN;
            
            if show_fig
                figure; plot(pres.a2dOsservazioni'); title('Pressure before outDetection');
            end
            pr_orig = pres.a2dOsservazioni';
            
            % removing Outliers
            pres.a2dOsservazioni(pres.a2dOsservazioni < 200) = NaN;
            pres.a2dOsservazioni(pres.a2dOsservazioni > 1100) = NaN;
            pres.a2dOsservazioni(std(pres.a2dOsservazioni', 'omitnan') > 10, :) = NaN;
            if show_fig
                figure; plot(pres.a2dOsservazioni'); title('Pressure after outDetection');
            end
            
            id_out = (pres.a2dOsservazioni' - pr_orig) ~= 0;
            if show_fig
                figure; imagesc(id_out); title('Changed pressure');
                dockAllFigures;
            end
            
            sta_name_list = char(pres.a1sNomiStazioni);

            n_ok = sum(~isnan(pr_orig));
            n_ko = sum(isnan(pr_orig));
            n_out = sum(id_out);
            n_tot = size(id_out, 1);
              
            fprintf('%d Stations do not contain data\n', sum(n_ko == n_tot));
            for s = find(n_ok == 0)
                fprintf('No data found for %4d %s\n', s, sta_name_list(s,:));
            end
            for s = find(n_ko == n_out & n_ok > 0 & n_ko > 0)
                fprintf('%2d / %2d Invalid data found (%8.2f, %8.2f) for %4d %s\n', n_out(s), n_tot, mean(pr_orig(:,s), 'omitnan'),  std(pr_orig(:,s), 'omitnan'), s, sta_name_list(s,:));
            end
            for s = find(((n_out - n_ko) > 0))
                fprintf('%2d / %2d outliers found (%8.2f, %8.2f) for %4d %s\n', (n_out(s) - n_ko(s)), n_ok(s), mean(pr_orig(:,s), 'omitnan'),  std(pr_orig(:,s), 'omitnan'), s, sta_name_list(s,:));
            end
            
            time = GPS_Time(datenum(char(pres.a1sTempi),'dd/mm/yyyy HH:MM'));            
            for s = 1 : size(sta_name_list, 1) * double(do_export)
                if any(pres.a2dOsservazioni(s,:))
                    sta_name = strrep(strtrim(sta_name_list(s, :)),' ', '_');
                    log.addMessage(sprintf('Exporting pressure of %s', strtrim(sta_name_list(s, :))));
                    log.setVerbosityLev(1);
                    file_name = [data_path filesep '${YYYY}_${DOY}' filesep sta_name '_PR_${DOY}.${YY}m'];
                    
                    ondu = getOrthometricCorr(pres.a1dLatitudini(s)./180*pi, pres.a1dLongitudini(s)./180*pi, Core.getRefGeoid());
                    h_ellips = pres.a1dAltitudini(s) + ondu;
                    [x, y, z] = geod2cart(pres.a1dLatitudini(s)./180*pi, pres.a1dLongitudini(s)./180*pi, h_ellips);
                    
                    md = Meteo_Data();
                    md.importRaw(time, [pres.a2dOsservazioni(s,:)' pres.a2dOsservazioni(s,:)'], [Meteo_Data.PR Meteo_Data.TD], sta_name, [x y z]);
                    md.export(file_name);
                    log.setVerbosityLev(v_lev);
                end
            end
            
            temp = load([data_path filesep 'temperatura.mat']);
            % Find Outliers
            
            temp.a2dOsservazioni(temp.a2dOsservazioni < -55) = NaN;
            tm_orig = temp.a2dOsservazioni';
            if show_fig
                figure; plot(temp.a2dOsservazioni'); title('Temperatures before outDetection');
            end
            tmp = temp.a2dOsservazioni;
            tmp = simpleFill1D(tmp', flagExpand(~isnan(tmp'),3) & isnan(tmp'), 'linear');
            dtmp = Core_Utils.diffAndPred(tmp)';            
            sensor = abs(bsxfun(@minus, dtmp', median(dtmp, 'omitnan')'));
            thr = mean(movmax(sensor,1)'*10,'omitnan');
            id_ko = bsxfun(@minus, sensor, thr') > 0;
            temp.a2dOsservazioni(id_ko' | isnan(sensor')) = NaN;
            if show_fig
                figure; plot(temp.a2dOsservazioni'); title('Temperatures after outDetection');
            end
            
            id_out = (temp.a2dOsservazioni' - tm_orig) ~= 0;
            if show_fig
                figure; imagesc(id_out); title('Changed pressure');
                dockAllFigures;
            end
            
            sta_name_list = char(temp.a1sNomiStazioni);
            
            n_ok = sum(~isnan(tm_orig));
            n_ko = sum(isnan(tm_orig));
            n_out = sum(id_out);
            n_tot = size(id_out, 1);
              
            fprintf('%d Stations do not contain data\n', sum(n_ko == n_tot));
            for s = find(n_ok == 0)
                fprintf('No data found for %4d %s\n', s, sta_name_list(s,:));
            end
            for s = find(n_ko == n_out & n_ok > 0 & n_ko > 0)
                fprintf('%2d / %2d Invalid data found (%8.2f, %8.2f) for %4d %s\n', n_out(s), n_tot, mean(tm_orig(:,s), 'omitnan'),  std(tm_orig(:,s), 'omitnan'), s, sta_name_list(s,:));
            end
            for s = find(((n_out - n_ko) > 0))
                fprintf('%2d / %2d outliers found (%8.2f, %8.2f) for %4d %s\n', (n_out(s) - n_ko(s)), n_ok(s), mean(tm_orig(:,s), 'omitnan'),  std(tm_orig(:,s), 'omitnan'), s, sta_name_list(s,:));
            end            
            
            sta_name_list = char(temp.a1sNomiStazioni);
            time = GPS_Time(datenum(char(temp.a1sTempi),'dd/mm/yyyy HH:MM'));            
            for s = 1 : size(sta_name_list, 1) * double(do_export)
                if any(temp.a2dOsservazioni(s,:))
                    sta_name = strrep(strtrim(sta_name_list(s, :)),' ', '_');
                    log.addMessage(sprintf('Exporting temperature of %4d/%4d %s', s, size(sta_name_list, 1), strtrim(sta_name_list(s, :))));
                    log.setVerbosityLev(1);
                    file_name = [data_path filesep '${YYYY}_${DOY}' filesep sta_name '_TM_${DOY}.${YY}m'];
                    
                    ondu = getOrthometricCorr(temp.a1dLatitudini(s)./180*pi, temp.a1dLongitudini(s)./180*pi, Core.getRefGeoid());
                    h_ellips = temp.a1dAltitudini(s) + ondu;
                    [x, y, z] = geod2cart(temp.a1dLatitudini(s)./180*pi, temp.a1dLongitudini(s)./180*pi, h_ellips);
                    
                    md = Meteo_Data();
                    md.importRaw(time, [temp.a2dOsservazioni(s,:)' temp.a2dOsservazioni(s,:)'], [Meteo_Data.PR Meteo_Data.TD], sta_name, [x y z]);
                    md.export(file_name);
                    log.setVerbosityLev(v_lev);
                end
            end
        end
        
        function met2Csv(file_list, file_out)
            if ~iscell(file_list)
                file_list = {file_list};
            end
            fid = fopen(file_out, 'Wb');
            if (fid < 0)
                Logger.getInstance.addError(sprintf('Not possible to open "%s" for writing', file_out));
            else
                fprintf(fid, 'time;temperature;pressure\n');
                for i = 1 : numel(file_list)
                    md = Meteo_Data(file_list{i});
                    time = md.getTime.toString('yyyy-mm-dd HH:MM:SS');
                    temp = num2str(md.getTemperature,'%+6.2f');
                    pres = num2str(md.getPressure,'%7.2f');
                    %uint8(';') = 59
                    fprintf(fid, '%s', [time char(ones(size(time,1),1)*59) temp char(ones(size(time,1),1)*59) pres char(ones(size(time,1),1)*10)]');
                end
                fclose(fid);
            end
        end
    end
end
