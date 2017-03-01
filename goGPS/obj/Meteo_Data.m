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
% ftp://igs.org/pub/data/format/rinex303.pdf

%--------------------------------------------------------------------------
%               ___ ___ ___ 
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.0
% 
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
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
%----------------------------------------------------------------------------------------------
classdef Meteo_Data < handle
    
    properties (Constant, GetAccess = private);
        PR_ID = sum(uint16('PR').*uint16([1 256]));     % Internal id of Pressure [mbar]
        TD_ID = sum(uint16('TD').*uint16([1 256]));     % Internal id of Dry temperature (deg Celsius)        
        HR_ID = sum(uint16('HR').*uint16([1 256]));     % Internal id of Relative humidity (percent)
        ZW_ID = sum(uint16('ZW').*uint16([1 256]));     % Internal id of Wet zenith path delay (mm), (for WVR data)
        ZD_ID = sum(uint16('ZD').*uint16([1 256]));     % Internal id of Dry component of zen.path delay (mm)
        ZT_ID = sum(uint16('ZT').*uint16([1 256]));     % Internal id of Total zenith path delay (mm)
        WD_ID = sum(uint16('WD').*uint16([1 256]));     % Internal id of Wind azimuth (deg) from where the wind blows
        WS_ID = sum(uint16('WS').*uint16([1 256]));     % Internal id of Wind speed (m/s)
        RI_ID = sum(uint16('RI').*uint16([1 256]));     % Internal id of Rain increment (1/10 mm): Rain accumulation since last measurement
        HI_ID = sum(uint16('HI').*uint16([1 256]));     % Internal id of Hail indicator non-zero: Hail detected since last measurement                
    end
        properties (Constant, GetAccess = public);
        DATA_TYPE = ['PR'; 'TD'; 'HR'; 'ZW'; 'ZD'; 'ZT'; 'WD'; 'WS'; 'RI'; 'HI'];
        DATA_TYPE_EXT = { 'Pressure [mbar]', ...
                          'Dry temperature (deg Celsius)', ...
                          'Relative humidity (percent)', ...
                          'Wet zenith path delay (mm), (for WVR data)', ...
                          'Dry component of zen.path delay (mm)', ...
                          'Total zenith path delay (mm)', ...
                          'Wind azimuth (deg) from where the wind blows', ...
                          'Wind speed (m/s)', ...
                          'Rain increment (1/10 mm): Rain accumulation since last measurement', ...
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
        HI = 10;    % Numeric id of Hail indicator non-zero: Hail detected since last measurement
        
        DATA_TYPE_ID = [Meteo_Data.PR_ID, ... 
                        Meteo_Data.TD_ID, ... 
                        Meteo_Data.HR_ID, ... 
                        Meteo_Data.ZW_ID, ... 
                        Meteo_Data.ZD_ID, ... 
                        Meteo_Data.ZT_ID, ... 
                        Meteo_Data.WD_ID, ... 
                        Meteo_Data.WS_ID, ... 
                        Meteo_Data.RI_ID, ... 
                        Meteo_Data.HI_ID];
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        logger = Logger.getInstance(); % Handler to the logger object
    end
    
    properties (SetAccess = private, GetAccess = public)
        % contains an object to read the RINEX file 
        file;   % init this with File_Rinex('filename')
        
        n_type = 0;         % number of observation types
        type = 1:10;        % supposing to have all the fields
        
        obs_epoc = GPS_Time(); % array of observation epochs
        
        pr = 0;     % Column of Pressure [mbar]
        td = 0;     % Column of Dry temperature (deg Celsius)        
        hr = 0;     % Column of Relative humidity (percent)
        zw = 0;     % Column of Wet zenith path delay (mm), (for WVR data)
        zd = 0;     % Column of Dry component of zen.path delay (mm)
        zt = 0;     % Column of Total zenith path delay (mm)
        wd = 0;     % Column of Wind azimuth (deg) from where the wind blows
        ws = 0;     % Column of Wind speed (m/s)
        ri = 0;     % Column of Rain increment (1/10 mm): Rain accumulation since last measurement
        hi = 0;     % Column of Hail indicator non-zero: Hail detected since last measurement
        
        data = [];        
    end
    
    methods (Access = private)
        function parseHeader(this, meteo_file)
            % Parse header and update the object, having as input the meteo_file 
            % as read with textscan (cell array of string lines)
            % SYNTAX: this.parseHeader(meteo_file)
            
            if ~this.file.isValid(1)
                this.logger.addWarning('Meteorological file with no header or corrupted');
                this.logger.addWarning(sprintf('Try to read it as it contains date (6 fields) + %s', sprintf('%c%c ', this.getType()')));
            else           
                % Try to parse header
                try
                    % Check for file description
                    if isempty(strfind(meteo_file{1}(61:end),'VERSION / TYPE')) || isempty(strfind(meteo_file{1},'METEOROLOGICAL DATA'))
                        throw(MException('MeteorologicalFile:InvalidHeader', 'file type is not described as "METEOROLOGICAL DATA"'));
                    end
                catch ex
                    this.logger.addWarning(sprintf('Problem detected in header of %s: %s', file.getFileName(), ex.message))
                end
            
                % Scan the file header
                try
                    l = 1;
                    this.type = [];
                    type_is_present = false;
                    while (l < this.file.getEOH)
                        l = l + 1;
                        line = meteo_file{l};
                        if (~type_is_present)
                            % Check for type of observations
                            type_is_present = ~isempty(strfind(line(61:end),'# / TYPES OF OBSERV'));
                            
                            if type_is_present
                                % Read types of obs
                                this.n_type = sscanf(line(1:6),'%d');
                                for t = 0 : (this.n_type - 1)
                                    str_type = line(8 + (3:4) + t * 6);
                                    this.type = [this.type find(this.DATA_TYPE_ID == (sum(uint16(str_type) .* uint16([1 256]))))];
                                end
                                if this.n_type > numel(this.type)
                                    throw(MException('MeteorologicalFile:InvalidHeader', 'unrecognized data type'));
                                end
                            end
                        end
                    end
                    if isempty(this.type)
                        this.type = [1:10];
                        this.n_type = numel(this.type);
                    end
                    if ~type_is_present
                        throw(MException('MeteorologicalFile:InvalidHeader', 'file does not contain the header field "# / TYPES OF OBSERV"'));
                    end  
                catch ex
                    this.logger.addWarning(sprintf('Problem detected in the header: %s', ex.message))
                    this.logger.addWarning(sprintf('try to read it as it contains date (6 fields) + %s', sprintf('%c%c ', this.getType()')));
                end  
            end
        end
    end
    
    methods
        function this = Meteo_Data(file_name, type) 
            % Try to read the file
            this.logger.addMessage(sprintf('Loading Metereological data from %s\n', file_name));
            try
                this.file = File_Rinex(file_name);
                fid = fopen(this.file.getFileName(), 'r');
                meteo_file = textscan(fid,'%s','Delimiter', '\n', 'whitespace', '');
                fclose(fid);
                meteo_file = meteo_file{1};
            catch ex
                this.logger.addError(sprintf('Error reading metereological file %s (%s)', file_name, ex.message));
                return
            end
            
            % Parse the header and detect the types of data contained in the metereological file
            this.parseHeader(meteo_file);
            
            if (nargin == 2)
                this.logger.addWarning('Overriding the file type with custom types');
                this.type = type;
                this.n_type = numel(type);
                assert(sum(type < 1) + sum(type > 10) == 0, 'Invalid custom types');
            end
            
            this.logger.addMessage(sprintf('\nThe following meteorological data are present in the header:'));
            for t = 1:this.n_type
                this.logger.addMessage([' - ' this.getTypeExt{t}]);
            end
            
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
                    this.obs_epoc.addEpoch(line(id_date), [], true);
                    [value, ~, id_stop] = regexp(line((id_date(end) + 1):end), '[.0-9]*', 'match');
                    col_id = round(id_stop / 7);
                    this.data(col_id, l - eoh) = str2double(value);
                end
            catch ex 
                this.logger.addWarning(sprintf('Problem detected while reading metereological data: %s', ex.message));
            end
            this.data = this.data'; % keep one column per data type
        end
    end
    
    methods        
        function type = getType(this)
            % Get the types of data stored in the RINEX
            id = this.getTypeId();
            type = this.DATA_TYPE(id,:);
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
        
    end        
end
