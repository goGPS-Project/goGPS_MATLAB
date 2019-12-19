%   CLASS Meteo_Data_Old
% =========================================================================
%
% DESCRIPTION
%   Class to store receiver data (observations, and characteristics
%
% EXAMPLE
%   settings = Meteo_Data_Old();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Meteo_Data_Old
%
% REFERENCE

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
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
classdef Meteo_Data_Old < handle

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
        DATA_TYPE_ID = [Meteo_Data_Old.PR_ID, ...
                        Meteo_Data_Old.TD_ID, ...
                        Meteo_Data_Old.HR_ID, ...
                        Meteo_Data_Old.ZW_ID, ...
                        Meteo_Data_Old.ZD_ID, ...
                        Meteo_Data_Old.ZT_ID, ...
                        Meteo_Data_Old.WD_ID, ...
                        Meteo_Data_Old.WS_ID, ...
                        Meteo_Data_Old.RI_ID, ...
                        Meteo_Data_Old.RR_ID, ...
                        Meteo_Data_Old.HI_ID];
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

    properties (SetAccess = private, GetAccess = protected)
        % contains an object to read the RINEX file
        file;   % init this with File_Rinex('filename')

        marker_name = '';   % name of the station
        n_type = 0;         % number of observation types
        type = 1:11;        % supposing to have all the fields

        time = GPS_Time();  % array of observation epochs

        data = [];          % Meteorological file

        xyz = [0 0 0];      % geocentric coordinate of the sensor
        amsl = 0;           % hortometric height of the sensor
        is_valid = false;   % Status of valitity of the file;
        
        max_bound = 90;     % Max bound to extrapolate
        
        smoothing = [900 300 300 0 0 0 0 0 0 0 0]; % Spline base for smoothing (in seconds)
    end

    methods (Access = private)
        function parseHeader(this, meteo_file)
            % Parse header and update the object, having as input the meteo_file
            % as read with textscan (cell array of string lines)
            % SYNTAX: this.parseHeader(meteo_file)

            if ~this.file.isValid(1)
                this.log.addWarning('Meteorological file with no header or corrupted');
                this.log.addWarning(sprintf('Try to read it as it contains date (6 fields) + %s', sprintf('%c%c ', this.getType()')));
            else
                % Try to parse header
                try
                    % Check for file description to verify the header
                    if isempty(strfind(meteo_file{1}(61:end),'VERSION / TYPE')) || isempty(strfind(meteo_file{1},'METEOROLOGICAL DATA'))
                        throw(MException('MeteorologicalFile:InvalidHeader', 'file type is not described as "METEOROLOGICAL DATA"'));
                    end
                catch ex
                    this.log.addWarning(sprintf('Problem detected in header of %s: %s', file.getFileName(), ex.message))
                end

                % Scan the file header
                try
                    l = 1;
                    this.type = [];
                    type_is_present = false;
                    pos_is_present = false;
                    marker_name_is_present = false;
                    while (l < this.file.getEOH)
                        l = l + 1;
                        line = meteo_file{l};
                        if (~type_is_present)
                            % Check for type of observations
                            type_is_present = ~isempty(strfind(line(61:end),'# / TYPES OF OBSERV'));

                            if type_is_present
                                % Read types of obs
                                this.n_type = sscanf(line(1:6),'%d');
                                if (this.n_type >= 10)
                                    this.log.addWarning('Reading more than 10 fields is not yet supported');
                                end
                                for t = 0 : (this.n_type - 1)
                                    str_type = line(8 + (3:4) + t * 6);
                                    this.type = [this.type find(this.DATA_TYPE_ID == (sum(uint16(str_type) .* uint16([1 256]))))];
                                end
                                if this.n_type > numel(this.type)
                                    throw(MException('MeteorologicalFile:InvalidHeader', 'unrecognized data type'));
                                end
                            end
                        end
                        if (~pos_is_present)
                            % Check for sensor position
                            pos_is_present = ~isempty(strfind(line(61:end),'SENSOR POS XYZ/H'));

                            if pos_is_present
                                this.xyz = sscanf(line(1:42), '%14f%14f%14f');
                                if sum(abs(this.xyz)) > 0
                                    [~, lam, h, phiC] = cart2geod(this.xyz(1), this.xyz(2), this.xyz(3));
                                    this.amsl = h - getOrthometricCorr(phiC, lam);
                                else
                                    this.amsl = 0;
                                end
                            end
                        end
                        if (~marker_name_is_present)
                            % Check for sensor position
                            marker_name_is_present = ~isempty(strfind(line(61:end),'MARKER NAME'));
                            if marker_name_is_present
                                this.marker_name = strtrim(line(1:60));
                            end
                        end
                    end

                    % set default (empty) values if the following paremeters are not found in header
                    if ~pos_is_present
                        this.xyz = [0 0 0];
                        this.amsl = 0;
                    end

                    if ~marker_name_is_present
                        this.marker_name = '';
                    end

                    if isempty(this.type)
                        this.type = 1 : 10;
                        this.n_type = numel(this.type);
                    end
                    if ~type_is_present
                        throw(MException('MeteorologicalFile:InvalidHeader', 'file does not contain the header field "# / TYPES OF OBSERV"'));
                    end
                catch ex
                    this.log.addWarning(sprintf('Problem detected in the header: %s', ex.message))
                    this.log.addWarning(sprintf('try to read it as it contains date (6 fields) + %s', sprintf('%c%c ', this.getType()')));
                end
            end
        end

        function parseData(this, meteo_file)
            % Parse data and update the object, having as input the meteo_file
            % as read with textscan (cell array of string lines)
            % SYNTAX: this.parseData(meteo_file)

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
            if (nargin < 4)
                verbosity_lev = Logger.DEFAULT_VERBOSITY_LEV;
            end

            this.log.addMessage(sprintf('Loading Meteorological data from "%s"', file_name), iif(verbosity_lev < 50, Logger.DEFAULT_VERBOSITY_LEV, verbosity_lev));
            try
                this.time = GPS_Time(); % empty the time
                this.file = File_Rinex(file_name, verbosity_lev);
                fid = fopen(this.file.getFileName(), 'r');
                meteo_file = textscan(fid,'%s','Delimiter', '\n', 'whitespace', '');
                fclose(fid);
                meteo_file = meteo_file{1};
            catch ex
                this.log.addError(sprintf('Error reading meteorological file "%s" (%s)', file_name, ex.message));
                return
            end

            % Parse the header and detect the types of data contained in the meteorological file
            this.parseHeader(meteo_file);

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
            this.parseData(meteo_file);
            if ~any(this.xyz)
                this.log.addWarning(sprintf('No position found in meteorological file "%s"\n this meteorological station cannot be used correctly', File_Name_Processor.getFileName(file_name)), verbosity_lev);                
            end
        end
    end
    
    % =========================================================================
    %  SETTER
    % =========================================================================
    methods
        function setMaxBound(this, max_bound)
            % Set the maximum extrapolation span
            this.max_bound = max_bound;
        end
        
    end

    % =========================================================================
    %  INIT / READER
    % =========================================================================
    methods
        function this = Meteo_Data_Old(file_name, type, verbosity_lev)
            % Creator Meteo_Data_Old(file_name, <type = empty>, <verbosity_lev>)

            % The function calls all its creation methods within try and
            % catch statements, reading the Meteo file should not be
            % blocking for the processing, even if the data are not good

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
            % EXAMPLE: this.importRaw(GPS_Time(time - 1/12), [pres temp hum rain], [Meteo_Data_Old.PR Meteo_Data_Old.TD Meteo_Data_Old.HR Meteo_Data_Old.RT], 'GReD', xyz);
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
            ok = ok & ~invalid_epoch;
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
            % SYNTAX: this.toString(str)
            
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
                        file_name = file_name(1:4);
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
                        fid = fopen(cur_file_name, 'w');
                        str = ['     3.03           METEOROLOGICAL DATA                     RINEX VERSION / TYPE', 10 ...
                            'EXPORTED MET FILE FROM Meteo_Data_Old MATLAB CLASS              COMMENT', 10];
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
            % SYNTAX: this.toString(str)

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
        function validity = isValid(this)
            % Get the validity of a RINEX file or the object
            % SYNTAX: validity = isValid()
            validity = this.file.isValid();
        end
                
        function name = getMarkerName(this)
            % Get the name of the station
            % SYNTAX: time = this.getMarkerName()
            name = this.marker_name;
        end
        
        function time = getTime(this)
            % Get the epochs of the data
            % SYNTAX: time = this.getTime()
            time = this.time;
        end
        
        function type = getType(this)
            % Get the types of data stored in the RINEX
            id = this.getTypeId();
            type = this.DATA_TYPE(id,:);
        end
        
        function max_bound = getMaxBound(this)
            % Get the maximum extrapolation
            max_bound = this.max_bound;
        end
        
        function type = getTypeExt(this)
            % Get the description of the types of data stored in the RINEX
            id = this.getTypeId();
            type = this.DATA_TYPE_EXT(id);
        end
        
        function id = getTypeId(this)
            % Get the id of the types of data stored in the RINEX
            id = this.type;
        end
        
        function data = getComponent(this, data_id, time)
            % Get the data with id of the type wanted
            % Passing a time array as GPS_Time the object interpolate the
            % data contained in the meteorological file
            % SYNTAX: data = this.getComponent(id, <time>)
            id = find(this.type == data_id);
            if isempty(id)
                if nargin == 3
                    data = nan(time.length(), 1);
                else
                    data = [];
                end
            else
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
                            data_in = data_in(~isnan(data_in));
                            time_data = time_data(~isnan(data_in));
                            data_in = [data_in(1); data_in; data_in(end)];
                            time_data = [ (min(time_pred(1), time_data(1)) - 1/86400); time_data; (max(time_pred(end), time_data(end)) + 1/86400)];
                            data = interp1(time_data, data_in, time_pred, 'pchip','extrap');
                            if this.smoothing(data_id) > 0
                                data = splinerMat(time_pred, data - mean(data), this.smoothing(data_id) / 86400, 0) + mean(data);
                            end
                            % do not extrapolate further than 20 minutes in time
                            data((time_pred < time_data(1) - this.getMaxBound / 1440) | (time_pred > time_data(end) + this.getMaxBound / 1440)) = NaN;
                            % extrapoleted value 
                        end
                    end
                else
                    data = data_in;
                end
            end
        end
        
        function data = getPressure(this, time, amsl)
            % Get the pressure data
            % SYNTAX: data = this.getPressure()
            if (nargin == 1)
                data = this.getComponent(1);
            else
                data = this.getComponent(1, time);
            end
            
            if (nargin == 3)
                data = Meteo_Data_Old.pressure_adjustment(data, this.amsl, amsl);
            end
        end
        
        function data = getTemperature(this, time, amsl)
            % Get the temperature data
            % SYNTAX: data = this.getTemperature()
            if (nargin == 1)
                data = this.getComponent(2);
            else
                data = this.getComponent(2, time);
            end
            
            if (nargin == 3)
                data = Meteo_Data_Old.temperature_adjustment(data, this.amsl, amsl);
            end
        end
        
        function data = getHumidity(this, time, amsl)
            % Get the humidity data
            % SYNTAX: data = this.getHumidity()
            if (nargin == 1)
                data = this.getComponent(3);
            else
                data = this.getComponent(3, time);
            end
            
            if (nargin == 3)
                data = Meteo_Data_Old.humidity_adjustment(data, this.amsl, amsl);
            end
        end
        
        function [x, y, z, amsl] = getLocation(this)
            % Get meteo station location
            % SINTAX: [x, y, z, amsl] = this.getLocation();
            x = this.xyz(1);
            y = this.xyz(2);
            z = this.xyz(3);
            amsl = this.amsl;
        end
        
        function time = getObsTime(this)
            % Get meteo station observatipon time
            % SINTAX: time = this.getObsTime();
            time = this.time;
        end
    end

    % =========================================================================
    %  STATIC
    % =========================================================================
    methods (Static)
        function md = getVMS(name, xyz, time, station)
            % Get Virtual Meteo Station
            %
            % INPUT
            %   name    name of the new Meteo_Data_Old virtual station
            %   xyz     coordinates of the new meteo station
            %   time    time of interpolation
            %
            %   station list of Meteo_Data_Old station to use for the interpolation
            %   
            % OUTPUT
            %   md      virtual Meteo_Data_Old station generated at xyz coordinates
            %   
            % SYNTAX
            %   md = this.getVMS(name, xyz, time, station)
            %
            % EXAMPLE
            %   [x, y, z, amsl] = station(1).getLocation();
            %   md1 = Meteo_Data_Old.getVMS('test', [x y z], station(1).getObsTime, md)

            log = Logger.getInstance();

            md = Meteo_Data_Old();
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
            d2 = sqrt(abs(e_mesh - e_mesh').^2 + abs(n_mesh - n_mesh').^2);
            d = sqrt(abs(e_obs - e).^2 + abs(n_obs - n).^2);

            % fun for pressure
            fun = @(dist) 0.2 * exp(-(dist/0.8e4)) + exp(-(dist/6e3).^2);
            q_fun_obs = fun(d2) .* repmat(fun(d)', size(d2,1), 1);
            q_fun_obs = triu(q_fun_obs) + triu(q_fun_obs,1)';

            % getting pressure
            id_pr = find(st_type(:,1) == 1);
            % id_pr = id_pr(real_dist(id_pr) < 20e3);

            pr_obs = zeros(numel(id_pr), time.length());
            for s = 1 : numel(id_pr)
                pr_obs(s, :) = station(id_pr(s)).getPressure(time, amsl);
            end

            id_pr(sum(isnan(pr_obs),2) > 0) = [];
            pr_obs(sum(isnan(pr_obs),2) > 0, :) = [];

            if isempty(id_pr)
                log.addWarning('There are no station to get pressure information', 100);
                pres = nan(time.length,1);
            else
                %A = ones(size(id_pr));
                %Q = d2(id_pr, id_pr);
                %AinvQ =  A'/Q;
                %w = (AinvQ*A)\AinvQ;
                trans = sum(q_fun_obs(id_pr, id_pr));
                w = sum(trans)\trans;
                pres = (w * pr_obs)';
            end

            % fun for temperature
            fun = @(dist) 0.2 * exp(-(dist/1e4)) + exp(-(dist/6e3).^2);
            q_fun_obs = fun(d2) .* repmat(fun(d)', size(d2,1), 1);
            q_fun_obs = triu(q_fun_obs) + triu(q_fun_obs,1)';

            % getting temperature
            id_td = find(st_type(:,2) == 1);

            td_obs = zeros(numel(id_td), time.length);
            for s = 1 : numel(id_td)
                td_obs(s, :) = station(id_td(s)).getTemperature(time, amsl);
            end
            id_td(sum(isnan(td_obs),2) > 1) = [];
            td_obs(sum(isnan(td_obs),2) > 1, :) = [];

            if isempty(id_td)
                log.addWarning('There are no station to get temperature information', 100);
                temp = nan(time.length,1);
            else
                trans = sum(q_fun_obs(id_td, id_td));
                w = sum(trans)\trans;
                temp = (w * td_obs)';
            end

            % fun for humidity
            fun = @(dist) exp(-(dist/1e4)) + exp(-(dist/8e3).^2);
            q_fun_obs = fun(d2) .* repmat(fun(d)', size(d2,1), 1);
            q_fun_obs = triu(q_fun_obs) + triu(q_fun_obs,1)';

            % getting humidity
            id_hr = find(st_type(:,2) == 1);

            hr_obs = zeros(numel(id_hr), time.length);
            for s = 1 : numel(id_hr)
                hr_obs(s, :) = station(id_hr(s)).getHumidity(time, amsl);
            end

            id_hr(sum(isnan(hr_obs),2) > 1) = [];
            hr_obs(sum(isnan(hr_obs),2) > 1, :) = [];

            if isempty(id_hr)
                log.addWarning('There are no station to get relative humidity information', 100);
                hum = nan(time.length,1);
            else
                trans = sum(q_fun_obs(id_hr, id_hr));
                w = sum(trans)\trans;
                hum = (w * hr_obs)';
            end

            data = [pres temp hum];
            id_ok = (sum(isnan(data)) < time.length());
            data = data(:, id_ok);
            type = [Meteo_Data_Old.PR Meteo_Data_Old.TD Meteo_Data_Old.HR];
            type = type(:, id_ok);
            md.importRaw(time, data, type, name, xyz);
        end

        function [ temperature_adj ] = temperature_adjustment( temperature , obs_h, pred_h)
            % Barometric formula taken from Bai and Feng, 2003.
            % The parameter value is taken from Realini et al., 2014
            % Parameter definition
            grad = 0.0065 ; % * C / m gravitational acceleration constant
            temperature_adj = temperature + grad * (obs_h - pred_h);
        end

        function [ pressure_adj ] = pressure_adjustment( pressure , obs_h, pred_h)
            % Barometric formula taken from Berberan-Santos et al., 1997
            % The parameter values are taken from Realini et al., 2014
            % Parameters definition
            g = 9.80665 ;    % m / s^2 gravitational acceleration constant
            Md = 0.0289644 ; % kg / mol molar mass of dry air
            R = 8.31432 ;    % J(mol * K) gas constant for air
            Tisa = 288.15 ;  % K international standard temperature of the atmosphere at the sea level
            pressure_adj = pressure .* exp(-(g * Md * (pred_h - obs_h))/(R * Tisa ));
        end

        function [ humidity_adj ] = humidity_adjustment( humidity , obs_h, pred_h)
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
            humidity_adj = humidity / exp(-6.396e-4 * (obs_h - pred_h));
        end
    end
end
