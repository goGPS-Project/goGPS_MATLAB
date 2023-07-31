%   CLASS Radiosonde
% =========================================================================
%
% DESCRIPTION
%   Class to manage the download and comparison of Radiosondes around the World
%   This class currently uses data collected from the University of Wyoming at 
%   - http://weather.uwyo.edu/upperair/sounding.html
%   
%   It could be useful in a future to use/add the NOAA raob archives, bigger but every file contains the entire history of the soundings
%   - https://www.ncdc.noaa.gov/data-access/weather-balloon/integrated-global-radiosonde-archive
%
%
% FOR A LIST OF CONSTANTs and METHODS use doc Radiosonde

classdef Radiosonde < handle

    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 1.0
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
    %  Written by:        Andrea Gatti, Alice Bonfiglio, Stefano Barindelli
    %  Contributors:      Andrea Gatti, Alice Bonfiglio, Stefano Barindelli, Alessandra Mascitelli
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
    
    properties
        name            % station name
        sta_num         % station number
        
        lat             % latitude  [rad]
        lon             % logitude  [rad]
        elevation       % orthometric height [m]
        
        data_time       % time as datetime                  datetime [n_records x 1]
        pressure        % pressure [hPa]                    double   [n_records x 1]
        height          % height [m]                        double   [n_records x 1]
        temperature     % temperature [ÿC]                  double   [n_records x 1]
        rel_humidity    % relatibe humidity [%]             double   [n_records x 1]
        
        ref_time        % launch epoch;
        pwv             % precipitable water vapor [mm]     double   [n_launches x 1]
        ztd             % zenith tropospheric delay
    end
    
    properties (Constant)
        MAX_DIST = 100;  % Maximum distance in Km from a station to consider it valid
        
        JAPAN_STATION = {'47401', '47418', '47412', '47580', '47582', '47600', '47646', '47681', '47678', '47741', '47778', '47807', '47827', '47909', '47945', '47918'};
        ITALY_STATION = {'16045', '16064', '16113', '16245', '16320', '16429', '16546'};
        ITALY_NOVARA = '16064';
    end
    
    methods (Static)
        % Creator
        function this = Radiosonde(region, sta_num, date_start, date_stop) %valid for one station for one month
            % Core object creator
            
            log = Logger.getInstance();
            this.reset();
            if nargin == 0
                % init only;
            else
                if nargin < 4
                    log.addMessage(sprintf('Error downloading the files, not enough input arguments'));
                else
                    log.addMessage(log.indent(sprintf(' - downloading files for %s-%s', region, sta_num)));
                end
                
                this.download(region, sta_num, date_start, date_stop);
            end
        end
        
        function rds_list = fromList(sta_num, date_start, date_stop)
            % Get a list of radiosondes objects downloading data from
            %  - http://weather.uwyo.edu/upperair/sounding.html
            %
            % SINTAX
            %   rds_list = Radiosonde.fromList(region, sta_num, date_start, date_stop)
            %
            % EXAMPLE
            %   rds_list = Radiosonde.fromList(Radiosonde.JAPAN_STATION, date_start, date_stop)
            
            % http://weather.uwyo.edu/cgi-bin/sounding?region=seasia&TYPE=TEXT%3ALIST&YEAR=2019&MONTH=05&FROM=0600&TO=0600&STNM=32150
            if ~iscell(sta_num)
                sta_num = {sta_num};
            end
            region = 'europe'; % it is not really necessary
            
            % Init out
            rds_list(numel(sta_num)) = Radiosonde();
            for s = 1 : numel(sta_num)
                try
                    rds_list(s) = Radiosonde(region, sta_num{s}, date_start, date_stop);
                catch
                    % catch timeout
                    pause(1 + rand(1) * 3);
                    rds_list(s) = Radiosonde(region, sta_num{s}, date_start, date_stop);
                end
            end
        end
        
        function rds_list = fromJapan(date_start, date_stop)
            % http://weather.uwyo.edu/cgi-bin/sounding?region=seasia&TYPE=TEXT%3ALIST&YEAR=2019&MONTH=05&FROM=0600&TO=0600&STNM=32150
            region = 'seasia';
            % there is also '32150' very north
            sta_num = {'47401', ...
                '47418', ...
                '47412', ...
                '47580', ...
                '47582', ...
                '47600', ...
                '47646', ...
                '47681', ...
                '47678', ...
                '47741', ...
                '47778', ...
                '47807', ...
                '47827', ...
                '47909', ...
                '47945', ...
                '47918', ...
                };
            % Init out
            rds_list(numel(sta_num)) = Radiosonde();
            for s = 1 : numel(sta_num)
                try
                    rds_list(s) = Radiosonde(region, sta_num{s}, date_start, date_stop);
                catch
                    % catch timeout
                    pause(1 + rand(1) * 3);
                    rds_list(s) = Radiosonde(region, sta_num{s}, date_start, date_stop);
                end
            end
        end
    end
    
    % =========================================================================
    %  METHODS
    % =========================================================================
    
    methods % Public Access
        function reset(this)
            this.name = '';
            this.sta_num = 0;
            
            this.lat = 0;
            this.lon = 0;
            this.elevation = 0;
            
            this.data_time = [];
            this.pressure  = [];
            this.height  = [];
            this.temperature  = [];
            this.rel_humidity  = [];
            
            this.ref_time = [];
            this.pwv  = [];
            this.ztd = [];
        end
                
        function rds_pwv = download(this, region, sta_num, date_start, date_stop)
            date_list = datevec(floor(date_start.getMatlabTime() * 2)/2 : 0.5 : floor(date_stop.getMatlabTime() * 2)/2);
            [~, id_unique] = unique(date_list(:,1)*1e5 + date_list(:,2));
            year = date_list(id_unique, 1);
            month = date_list(id_unique, 2);
            
            % sta_num_cell=cellstr(sta_num);
            j = 1;
            for e = 1 : size(month, 1)
                from = '0100';
                if e == 1
                    from = sprintf('%02d%02d', date_list(1,3), date_list(1,4));
                end
                if e == size(month, 1)
                    to = sprintf('%02d%02d', date_list(end,3), date_list(end,4));
                else
                    switch sprintf('%02d',month(e,:))
                        case {'01', '03', '05', '07', '08', '10', '12'}
                            to = '3112';
                        case {'04', '06', '09' '11'}
                            to = '3012';
                        case '02'
                            to = '2812';
                    end
                end
                
                [sta_path, sta_dir] = getLocalFilePath(sta_num, year(e), month(e,:), from, to);
                flag_download = false;
                if exist(sta_path, 'file')                    
                    try
                        fid = fopen(sta_path, 'rt');
                        char_array = fread(fid, '*char')';
                        fclose(fid);
                        flag_retry = 0;
                    catch
                        flag_download = true;
                    end
                else
                    try
                        if ~exist(sta_dir, 'dir')
                            mkdir(sta_dir);
                        end
                    catch
                        Core.getLogger.addWarning(sprintf('RAOB dir cannot be created at "%s"', sta_dir));
                    end
                    flag_download = true;
                end
                   
                if flag_download
                    options = weboptions;
                    options.Timeout = 10;
                    plot_type = 'TEXT';
                    % http://weather.uwyo.edu/cgi-bin/sounding?region=seasia&TYPE=TEXT%3ALIST&YEAR=2019&MONTH=05&FROM=0600&TO=0600&STNM=32150
                    address = ['http://weather.uwyo.edu/cgi-bin/sounding?region=' region '&TYPE=' plot_type '%3ALIST&YEAR=' sprintf('%04d', year(e)) '&MONTH=' sprintf('%04d', month(e,:)) '&FROM=' from '&TO=' to '&STNM=' sta_num(:,:)];
                    fprintf('           get %s\n', address);
                    flag_retry = 4;
                    while flag_retry > 1
                        try
                            char_array = webread(address);
                            flag_retry = 0;
                        catch
                            flag_retry = flag_retry - 1;
                            pause(0.5 + randi(1,1) * 5);
                        end
                    end
                    % Last try
                    if flag_retry == 1
                        % try again after some seconds
                        try
                            fprintf('           retry access to %s\n', address);
                            pause(6 + randi(1,1) * 5);
                            char_array = webread(address);
                            flag_retry = 0;
                        catch ex
                            fprintf('           download failed\n');
                            %Core_Utils.printEx(ex);
                        end
                    end
                    
                    % Try to save the file
                    
                    try
                        if ~isempty(char_array)
                            fid = fopen(sta_path, 'wb');
                            fwrite(fid, char_array, '*char');
                            fclose(fid);
                            Core.getLogger.addStatusOk(sprintf('RAOB file cached at "%s"', sta_path));
                        end
                    catch
                        Core.getLogger.addWarning(sprintf('RAOB file cannot be written at "%s"', sta_path));
                    end
                end
                if (flag_retry == 0)
                    char_array = regexprep(char_array,'<script.*?/script>','');
                    char_array = regexprep(char_array,'<style.*?/style>','');
                    char_array = regexprep(char_array,'<.*?>','');
                    if ~instr(char_array, 'Observations') || instr(char_array, 'Can''t get')
                        Core.getLogger.addWarning(sprintf('no data available in "" :-(', sta_path));
                        rds_pwv(1,j).pwv = [];
                        rds_pwv(1,j).datetime = [];
                        rds_pwv(1,j).sta_numer = sta_num;
                    else
                        this.lat = str2double(regexp(char_array,'(?<=(Station latitude\: ))([\-0-9]|\.)*','match', 'once'));
                        this.lon = str2double(regexp(char_array,'(?<=(Station longitude\: ))([\-0-9]|\.)*','match', 'once'));
                        this.elevation = str2double(regexp(char_array,'(?<=(Station elevation\: ))([0-9]|\.)*','match', 'once'));
                        this.name = regexp(char_array,['(?<=' sta_num ' ).+?(?=( Observations))'],'match', 'once');
                        this.sta_num = sta_num;

                        %pwv
                        str_pw = 'Precipitable water [mm] for entire sounding:';
                        str_idx_pw = strfind(char_array,str_pw);
                        pw_vec_d = [];
                        for i=1:length(str_idx_pw)
                            pw_vec(i,:) = char_array(1, str_idx_pw(i)+length(str_pw)+1 : str_idx_pw(i)+length(str_pw)+5);
                            pw_vec_d(i,:) = str2double(pw_vec(i,:));
                        end
                        
                        %time
                        str_obs_time = 'Observation time:';
                        str_idx_time = strfind(char_array, str_obs_time);                        
                        if isempty(str_idx_time)
                            Core.getLogger.addWarning(sprintf('no data available in "%s" :-(', sta_path));
                        else                            
                            for i=1:length(str_idx_time)
                                time_vec(i,:)=[char_array(1,str_idx_time(i)+length(str_obs_time)+1:str_idx_time(i)+length(str_obs_time)+7) char_array(1,str_idx_time(i)+length(str_obs_time)+8:str_idx_time(i)+length(str_obs_time)+11) '00'];
                            end
                            datetime_vec = datetime(time_vec,'InputFormat','yyMMdd/HHmmSS');
                            
                            str_idx_header = strfind(char_array,'-----------------------------------------------------------------------------');
                            str_idx_header = str_idx_header(2:2:end);
                            str_idx_header_fin = strfind(char_array,'Station information and sounding indices');
                            for i = 1 : length(str_idx_header)
                                id_1 = length(this.data_time) + 1;
                                try
                                    temp = char_array(str_idx_header(i)+79:str_idx_header_fin(i)-2);
                                    err_code = false;
                                catch
                                    err_code = true;
                                    log = Logger.getInstance();
                                    log.addError(sprintf('The radiosonde "%s" launched at %s have a corrupted entry\nDeleting "%s", retry later', sta_num, datestr(date_time, 'yyyy-mmm-dd HH:MM'), sta_path));
                                    try
                                        delete(sta_path);
                                    catch
                                    end
                                end
                                if not(err_code)
                                    data(:,i) = textscan(temp,'%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%[^\n\r]','delimiter','','multipleDelimsAsOne',false,'TreatAsEmpty',{'[]'},'EmptyValue',NaN,'whitespace','','CollectOutput',true);
                                    %pressure (hpa)
                                    pres = str2double(data{1,i}(:,1));
                                    %height (m)
                                    height = str2double(data{1,i}(:,2));
                                    %temperature (Celsius degree)
                                    temp = str2double(data{1,i}(:,3));
                                    %rel humidity (%)
                                    relh = str2double(data{1,i}(:,5));
                                    %radiosonde(i,1).name=sta_num(s,:);
                                    date_time = datetime_vec(i,:);
                                    
                                    if max(pres) < 800 % this is an outlier!!!!
                                        log = Logger.getInstance();
                                        log.addError(sprintf('The radiosonde "%s" launched at %s have a maximum invalid pressure (%.f mbar)', sta_num, datestr(date_time, 'yyyy-mmm-dd HH:MM'), max(pres)));
                                    else
                                        this.data_time = [this.data_time; repmat(date_time, length(pres), 1)];
                                        this.pressure = [this.pressure; pres];
                                        this.height = [this.height; height];
                                        this.temperature = [this.temperature; temp];
                                        this.rel_humidity = [this.rel_humidity; relh];
                                        [ ztd, err_code] = Radiosonde.rad2ztd(this.temperature(id_1:end), this.rel_humidity(id_1:end), this.pressure(id_1:end), this.height(id_1:end), this.lat);
                                        
                                        if err_code
                                            log = Logger.getInstance();
                                            log.addError(sprintf('The radiosonde "%s" launched at %s did not collected enough data\nthe minimum valid altitude logging all the data needs to be 8900m (99/100 of water vapour)', sta_num, datestr(date_time, 'yyyy-mmm-dd HH:MM')));
                                        end
                                        this.ztd = [this.ztd ztd];
                                        this.ref_time = [this.ref_time this.data_time(id_1)];
                                    end
                                end
                            end
                            
                            rds_pwv(1,j).pwv = pw_vec_d;
                            rds_pwv(1,j).datetime = datetime_vec;
                            rds_pwv(1,j).sta_numer = sta_num(:,:);
                            this.pwv = [this.pwv, rds_pwv(e).pwv'];
                            clear datetime_vec pw_vec_d time_vec
                        end
                    end
                end
                j = j + 1;
            end
            
            function [sta_path, sta_dir] = getLocalFilePath(sta_num, year, month, from, to)
                % Get the RAOB local path
                %
                % SYNTAX
                %   [sta_path, sta_dir] = getLocalFilePath(sta_num, year, month, from, to)
                data_dir = fullfile(Core.getInstallDir, '..' , 'data');
                fnp = File_Name_Processor;
                sta_dir = fullfile(sprintf('%s/station/RAOB/%s/', data_dir, sta_num));
                sta_path = fullfile(sta_dir, sprintf('%s_%04d_%02d_%s_%s.raob', sta_num, year, month, from, to));
            end            
        end
    end
    
    methods % Public Access
        function coo = getPos(this)
            % Get coordinate of the radiosonde launch site
            %
            % INPUT
            %   this    single radiosonde object
            %
            % OUTPUT 
            %   coo     Coordinate object
            %
            % SYNTAX
            %   coo = this.getPos()
            
            coo = Coordinates.fromGeodetic(this.lat, this.lon, [], this.elevation, GPS_Time(this.ref_time));
        end
        
        function height = getHeight(rds_list)
            % Get launch site orthometric height
            %
            % INPUT
            %   rds_list    radiosonde list
            %
            % OUTPUT 
            %   height      orthometric height [m]
            %
            % SYNTAX
            %   height = rds_list.getLat()
            
            height = [rds_list.elevation]';
        end

        function lat = getLat(rds_list)
            % Get latitude
            %
            % INPUT
            %   rds_list    radiosonde list
            %
            % OUTPUT 
            %   lat         latitude [deg]
            %
            % SYNTAX
            %   lat = rds_list.getLat()
            
            lat = [rds_list.lat]';
        end

        function lon = getLon(rds_list)
            % Get latitude
            %
            % INPUT
            %   rds_list    radiosonde list
            %
            % OUTPUT 
            %   lon         longitude [deg]
            %
            % SYNTAX
            %   lon = rds_list.getLon()
                        
            lon = [rds_list.lon]';
        end
        
        function height = getLaunchHeight(rds_list)
            % Get height (orthometric height)
            %
            % INPUT
            %   rds_list    radiosonde list
            %
            % OUTPUT 
            %   height      orthometric height [deg]
            %
            % SYNTAX
            %   heigth = rds_list.getLaunchHeight()
            
            height = [rds_list.elevation]';
        end
        
        function [ztd, time] = getZtd(this)
            % Get the zenit total delay computed from sounding
            %
            % SYNTAX
            %   [ztd, time] = getZtd(this)
            ztd = this.ztd * 1e2;
            time = GPS_Time(datenum(this.ref_time)');
        end
        
        function name = getName(sta_list)
            name = {};
            for s = 1 : numel(sta_list)
                name{s} = sta_list(s).name;
                if isempty(name{s})
                    name{s} = Radiosonde.getRaobName(sta_list(s).sta_num);
                end
                name{s} = strtrim(name{s});
            end
            if numel(sta_list) == 1
                name = name{1};
            end
        end
    end
    
    methods (Static)
        function [ztd, err_status] = rad2ztd(temperature, relative_humidity, pressure, height, lat)
            % Compute ZTD from RAOB (no saastamoinen)
            %
            % For question on how this script works ask to: stefano.barindelli (at) polimi.it
            %
            % OUTPUT
            %   ztd     zenith total delay [m]
            %   err     boolean, represent if the radiosonde data are valid
            %
            % SYNTAX
            %  [ztd, err] = Radiosonde.rad2ztd(temperature, relative_humidity, pressure, height, lat)
            
            err_status = 0;
            
            celsius_to_kelvin = 273.16; %tempearture conversion
            R = 8.31432;  %universal gas constant for air
            dry_molar_mass = 0.02896; % molar mass of dry air [kg mole-1]
            % wet_molar_mass = 0.0180153; % molar mass of H2o [kg mole-1]
            g = 9.80665; % Gravitational acceleration [m]
            
            if isnan(temperature(end)) && isnan(height(end)) && isnan(relative_humidity(end))
                temperature = temperature(1 : end - 1);
                relative_humidity = relative_humidity(1 : end - 1);
                pressure = pressure(1 : end - 1);
                height = height(1 : end - 1);
            end
            
            nan_T = (isnan(temperature));
            nan_RH = (isnan(relative_humidity));
            nan_P = (isnan(pressure));
            nan_height = (isnan(height));
            
            % find the last valid epoch
            lid_nan = nan_T | nan_RH | nan_P | nan_height;
            first_ok = find(~lid_nan, 1, 'first');
            last_ok = find(~lid_nan, 1, 'last');
            
            %"analyzing the radiosonde accumulated water vapor as a function of altitude a threshold value was determined; 99% of the total accumulated water vapor
            %was reached at an altitude between 8 and 9 km. A conservative threshold value of 10 km was thus chosen for tests in order to make a precise comparison
            %between radiosonde- and GPS-derived PWV."
            %Sato, K., Realini, E., Tsuda, T., Oigawa, M., Iwaki, Y., Shoji, Y., & Seko, H. (2013).
            % A high-resolution, precipitable water vapor monitoring system using a dense network of GNSS receivers.
            % cit. https://www.fujipress.jp/jdr/dr/dsstr000800010037/
            %In this case a threshold value of 8.9 km was set in order to consider empirical evidences
            
            if isempty(last_ok) || height(last_ok) < 8900
                err_status = 1;
                [ztd] = deal(nan);
            else
                if last_ok < length(lid_nan) || (first_ok > 1)
                    % cut last nan values
                    temperature = temperature(first_ok : last_ok);
                    relative_humidity = relative_humidity(first_ok : last_ok);
                    pressure = pressure(first_ok : last_ok);
                    height = height(first_ok : last_ok);
                    
                    nan_T = nan_T(first_ok : last_ok);
                    nan_RH = nan_RH(first_ok : last_ok);
                    nan_P = nan_P(first_ok : last_ok);
                    nan_height = nan_height(first_ok : last_ok);
                end
                
                if sum(nan_T)
                    flag_interval_ref = getFlagsLimits(nan_T);
                    temperature = simpleFill1D(temperature, uint32(flag_interval_ref), 0 );
                end
                
                if sum(nan_RH)
                    flag_interval_ref = getFlagsLimits(nan_RH);
                    relative_humidity = simpleFill1D(relative_humidity, uint32(flag_interval_ref), 0 );
                end
                
                if sum(nan_P)
                    flag_interval_ref = getFlagsLimits(nan_P);
                    pressure = simpleFill1D(pressure, uint32(flag_interval_ref), 0 );
                end
                
                if sum(nan_height)
                    %barometric formula for finding h
                    p0 = 101325; %Pa -> Sea level standard atmospheric pressure
                    %L = 0.0065; %K/m -> Temperature lapse rate, = g/c_p for dry air
                    %c_p = 1007; %J/(kgÿK) -> Constant-pressure specific heat
                    t0 = 288.15; %K -> Sea level standard temperature
                    flag_interval_ref = getFlagsLimits(nan_height);
                    for h = 1:size(flag_interval_ref,1)
                        for hh = flag_interval_ref(h,1):flag_interval_ref(h,2)
                            height(hh) = -log(pressure(hh)*100/p0)*t0*R/g/dry_molar_mass;
                        end
                    end
                    %height = fill1D(height, uint32(flag_interval_ref), 0 );
                end
                % height=[height];
                T = (temperature(1:end-1)+temperature(2:end))/2;
                RH = (relative_humidity(1:end-1)+relative_humidity(2:end))/2;
                P = (pressure(1:end-1)+pressure(2:end))/2;
                % delta_h = double(height(2:end)-height(1:end-1));
                
                
                %https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html(bolton,1980)
                %From Bevis et al. 1994 k1,k2,k3
                k1 = 7.76 * 10^-1; %[K/Pa] % from Vedel,2001
                k2 = 7.04 * 10^-1; %[K/Pa] % from Vedel,2001
                k3 = 3.739 * 10^3; %[K^2/Pa] % from Vedel,2001
                %http://glossary.ametsoc.org/wiki/Gas_constant
                % R_w = 461.5; %[J Kg^-1 K^-1]% Realini,2014
                R_d = 287; %[J Kg^-1 K^-1]% Realini,2014
                
                %eps is the ratio of the molecular weight of water vapour to that of dry air
                eps=18/28.94;
                
                
                % Saturation vapor pressure (e_s)
                e_s = 6.1078*exp((17.27*T)./(237.3+T));
                % The partial pressure of water vapor (e)
                e = double(RH)/100.*e_s;
                % Mixing ratio (r)
                r = 0.622*(e)./P;
                % dry air density (rho_d)
                rho_d = (dry_molar_mass/R)./((T+celsius_to_kelvin)).*P;
                % wet air density
                rho_w = rho_d.*r;
                %specific humidity (q)
                q = ((rho_w)./(rho_w+rho_d));
                
                %see Vedel et al. 2000 Conversion of WGS84 geometric heights to NWP model HIRLAM geopotential heights, Danish Meteorological Institute.
                
                delta_P = double(pressure(1:end-1)-pressure(2:end));
                g_e = 9.780356; %[m/s^2] %vedel et al. 2000
                a1 = 5.2885 * 10^-3; %vedel et al. 2000
                a2 = -5.9 * 10^-6; %vedel et al. 2000
                % lat = 45.4605 * pi/180; %latitude of the radiosonde station of Spino d'Adda (http://www.radiosonde.eu.bonplans.info/RS00-I/RS02C-I.html)
                lat = lat / 180 * pi;
                R_e = 6378.1; %[km] % average equatorial radius
                R_p = R_e-21.5;%[km] %average pole radius
                g_s = g_e*(1+a1*sin(lat)^2+a2*sin(2*lat)^2); %An approximate expression for the value of g at the geoid surface as function of latitude
                R_s = R_e / sqrt((R_e/R_p)^2*sin(lat)^2+cos(lat)^2); %distance to the center of the earth from that point of the geoid surface
                %height_g = (height(1:end-1)+height(2:end))/2;
                g = g_s * (R_s./(R_s+(height)*(10^-3))).^2; %variation of g with height
                
                g_fin = (g(1:end-1)+g(2:end))/2;
                
                % cit: Calculation of zenith delays from meteorological
                % data, comparison of NWP model, radiosonde and GPS delays
                % Vedel et al. 2001
                % http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.542.6157&rep=rep1&type=pdf
                % ZTD in meters
                ztd = ((10^-6)*sum((k1*(R_d./g_fin).*delta_P))+(10^-6)*sum(((R_d./(g_fin*eps)).*q.*((k2-(k1*eps))+(k3./(T+celsius_to_kelvin)))).*delta_P)) * 1e2;
                ztd = (ztd + (10^-4)*(k1*R_d*P(end)/g_s(end))*(1+2*(R_d*(T(end)+celsius_to_kelvin))/(R_s(end)*10^3*g_s(end))+2*((R_d*(T(end)+celsius_to_kelvin))/(R_s(end)*10^3*g_s(end)))^2));
                
            end
        end
        
        function fh_list = showMap(sta_list, new_fig, flag_labels, flag_polar)
            % Show Map of the stations
            % downloading the DTM and showing it
            %
            % CITATION:
            %   Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software],
            %   available online at www.eoas.ubc.ca/~rich/map.html.
            %
            % INPUT
            %   new_fig     open a new figure
            %   flag_labels show/not show labels
            %   flag_polar  could be 'N' / 'S' / false
            %
            % SYNTAX
            %   Radiosonde.showMap(sta_list, new_fig, flag_labels, flag_polar);
            %
            % EXAMPLE
            %   Radiosonde.showMap(Radiosonde.getRaobList())
            
            if nargin < 3 || isempty(flag_labels)
                flag_labels = false;
            end
            if nargin < 4 || isempty(flag_polar)
                flag_polar = false;
            end
            flag_large_points = true;
            point_size = 20;
            large_point_size = iif(flag_polar, 50, 25);
            point_color = [14, 25, 41]/256;
            resolution = 'high';
            
            if nargin < 2 || isempty(new_fig)
                new_fig = true;
            end
            if new_fig
                f = figure('Visible', 'off');
            else
                f = gcf;
                hold on;
            end
            
            fh_list = f;
            fig_name = sprintf('RaobMapDtm');
            f.UserData = struct('fig_name', fig_name);
            
            Logger.getInstance.addMarkedMessage('Preparing map, please wait...');
            
            maximizeFig(f);
            f.Color = [1 1 1];
            
            if isa(sta_list, 'Radiosonde')
                lat = sta_list.getLat();
                lon = sta_list.getLon();
                id_ko = isnan(lat);
                lat(id_ko) = [];
                lon(id_ko) = [];
                
                name = sta_list(~id_ko).getName;
                if ~iscell(name)
                    name = {name};
                end
                
                % set map limits
                if numel(sta_list) == 1
                    lon_lim = minMax(lon) + [-0.05 0.05];
                    lat_lim = minMax(lat) + [-0.05 0.05];
                else
                    lon_lim = minMax(lon); lon_lim = lon_lim + [-1 1] * diff(lon_lim)/15;
                    lat_lim = minMax(lat); lat_lim = lat_lim + [-1 1] * diff(lat_lim)/15;
                end
            else
                tmp = struct2array(sta_list);
                lat = [tmp.lat]'; lon = [tmp.lon]';
                lon_lim = [-180 180];
                lat_lim = [-90 90];
                name = {tmp.name};
                clear tmp;
            end
            
            lon_lim = max(-179.999, min(179.999, lon_lim));
            lat_lim = max(-89.999, min(89.999, lat_lim));
            
            nwse = [lat_lim(2), lon_lim(1), lat_lim(1), lon_lim(2)];
            clon = nwse([2 4]) + [-0.02 0.02];
            clat = nwse([3 1]) + [-0.02 0.02];
            clon = max(-180, min(180, clon));
            clat = max(-90, min(90, clat));
            
            if (flag_polar)
                if flag_polar == 'S'                    
                    m_proj('stereographic','lat',-90,'long',0,'radius', 25);
                    id_ko = lat > -65;
                else
                    m_proj('stereographic','lat',90,'long',0,'radius', 25);
                    id_ko = lat < 65;
                end
                lat(id_ko) = [];
                lon(id_ko) = [];
                name = name(~id_ko);
            else
                if isa(sta_list, 'Radiosonde')
                    m_proj('equidistant','lon',clon,'lat',clat);   % Projection
                else
                    m_proj('miller', 'lon', clon, 'lat', clat);   % Projection
                end
            end
            
            axes
            cmap = flipud(gray(1000)); colormap(cmap(150: end, :));
            
            % retrieve external DTM
            if flag_polar
                % use ETOPO instead
                colormap(gca, cmap(100 : end - 100, :))
                m_etopo2('shadedrelief', 'gradient', 3);
            else
                try
                    cmap = flipud(gray(1000)); colormap(gca, cmap(150: end, :));
                    [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', resolution);

                    sensor_pole = sum(~isnan(dtm) / size(dtm, 2),2);
                    id_pole = (find(sensor_pole > 0.5, 1, 'first') - 1);
                    if any(id_pole)
                        dtm(1 : id_pole, :) = repmat(dtm(id_pole + 1, : ), id_pole, 1);
                    end
                    dtm = flipud(dtm);
                    id_pole = (find(flipud(sensor_pole) > 0.5, 1, 'first') - 1);
                    if any(id_pole)
                        dtm(1 : id_pole, :) = repmat(dtm(id_pole + 1, : ), id_pole, 1);
                    end
                    
                    % comment the following line to have bathimetry
                    dtm(dtm < 0) = nan; % - 1/3 * max(dtm(:));
                    
                    % uncomment the following line to have colors
                    %colormap(gca, Cmap.adaptiveTerrain(minMax(dtm(:))));
                    drawnow;
                    
                    [shaded_dtm, x, y] = m_shadedrelief(lon_dtm, lat_dtm, dtm, 'nan', [0.98, 0.98, 1]);
                    %h_dtm = m_pcolor(lon_dtm, lat_dtm, dtm);
                    %h_dtm.CData = shaded_dtm;
                    
                    m_image(lon_dtm, lat_dtm, shaded_dtm);
                catch ex
                    % use ETOPO1 instead
                    colormap(gca, cmap(100 : end, :))
                    m_etopo2('shadedrelief','gradient', 3);
                end
            end
            
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
            
            if (flag_polar)
                if flag_polar == 'S'
                    m_grid('XaxisLocation', 'top', 'tickdir','in', 'fontsize', 16);
                else
                    m_grid('tickdir','in', 'fontsize', 16);
                end
                % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                ylim([-0.47 0.47]);
                drawnow
            else
                m_grid('box','fancy','tickdir','in', 'fontsize', 16);
                % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                drawnow
                if isa(sta_list, 'Radiosonde')
                    m_ruler([.7 1], -0.05, 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                end
            end
            [x, y] = m_ll2xy(lon, lat);
            
            %point_color = Cmap.get('viridis', numel(x));
            %point_size = 25;
            if size(point_color, 1) > 1
                scatter(x(:), y(:), point_size, 1:numel(x), 'filled'); hold on;
                colormap(point_color);
            else
                plot(x(:), y(:),'.', 'MarkerSize', point_size, 'Color', point_color); hold on;
            end
            if flag_labels
                % Label BG (in background w.r.t. the point)
                for r = 1 : numel(x)
                    tmp = name{r};
                    text(x(r), y(r), char(32 * ones(1, 4 + 2 * length(tmp), 'uint8')), ...
                        'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                        'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                        'Margin', 2, 'LineWidth', 2, ...
                        'HorizontalAlignment','left');
                end
            end
            
            if flag_large_points
                plot(x(:), y(:), 'ko', 'MarkerSize', large_point_size/3, 'LineWidth', 1 + 1 * logical(flag_polar));
                for r = 1 : numel(x)
                    plot(x(r), y(r), '.', 'MarkerSize', large_point_size, 'Color', Core_UI.getColor(r, numel(x)));
                end
                plot(x(:), y(:), '.k', 'MarkerSize', large_point_size/7);
            end
            
            if flag_labels
                for r = 1 : numel(x)
                    tmp = name{r};
                    t = text(x(r), y(r), ['   ' tmp], ...
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
            Core_UI.addExportMenu(f); Core_UI.addBeautifyMenu(f); Core_UI.beautifyFig(f, 'light');
            f.Visible = 'on';
            title(sprintf('Radiosondes position\\fontsize{5} \n'), 'FontSize', 16);
            %xlabel('Longitude [deg]');
            %ylabel('Latitude [deg]');
            ax = gca; ax.FontSize = 16;
            Logger.getInstance.addStatusOk('The map is ready ^_^');
        end
        
        function [name, raob_code] = getRaobName(sta_num)
            all_raob = {'16429  Trapani/Birgi (LICT)';
                '16546  Decimomannu (LIED)';
                '16754  Heraklion (Airport) (LGIR)';
                '17196  Kayseri (LTAU)';
                '17220  Izmir/Guzelyali';
                '17240  Isparta (LTBM)';
                '17281  Diyarbakir';
                '17351  Adana/Bolge';
                '40179  Bet Dagan';
                '40265  Mafraq (OJMF)';
                '40373  Al-Qaisumah (OEPA)';
                '40375  Tabuk (OETB)';
                '40394  Hail (OEHL)';
                '40417  K.F.I.A.-Dammam (OEDF)';
                '40430  Al-Madinah (OEMA)';
                '40437  King Khaled Intl Arpt (OERK)';
                '40706  Tabriz (OITT)';
                '40745  Mashhad (OIMM)';
                '40754  Tehran-Mehrabad (OIII)';
                '40766  Kermanshah (OICC)';
                '40800  Esfahan (OIFM)';
                '40848  Shiraz (OISS)';
                '40856  Zahedan (OIZH)';
                '40948  Kabul Airport (OAKB)';
                '41024  Jeddah (King Abdul Aziz) (OEJN)';
                '41112  Abha (OEAB)';
                '41217  Abu Dhabi Inter Arpt (OMAA)';
                '41256  Seeb; Intl Airport (OOMS)';
                '41316  Salalah (OOSA)';
                '41624  Dera Ismail Khan (OPDI)';
                '41715  Jacobabad (OPJA)';
                '41718  Khanpur';
                '41749  Nawabshah (OPNH)';
                '60018  Guimar-Tenerife';
                '60155  Casablanca (GMMC)';
                '60390  Dar-El-Beida (DAAG)';
                '60571  Bechar (DAOR)';
                '60680  Tamanrasset';
                '60715  Tunis-Carthage (DTTA)';
                '60760  Tozeur (DTTZ)';
                '63741  Nairobi/Dagoretti (HKNC)';
                '64400  Pointe-Noire (FCPP)';
                '64450  Brazzaville/Maya-Maya (FCBB)';
                '65578  Abidjan (DIAP)';
                '67083  Antananarivo/Ivato (FMMI)';
                '68263  Pretoria (Irene) (FAIR)';
                '68424  Upington (FAUP)';
                '68842  Port Elizabeth (FAPE)';
                '16429  Trapani/Birgi (LICT)';
                '16546  Decimomannu (LIED)';
                '16754  Heraklion (Airport) (LGIR)';
                '17196  Kayseri (LTAU)';
                '17220  Izmir/Guzelyali';
                '17240  Isparta (LTBM)';
                '17281  Diyarbakir';
                '17351  Adana/Bolge';
                '40179  Bet Dagan';
                '40265  Mafraq (OJMF)';
                '40373  Al-Qaisumah (OEPA)';
                '40375  Tabuk (OETB)';
                '40394  Hail (OEHL)';
                '40417  K.F.I.A.-Dammam (OEDF)';
                '40430  Al-Madinah (OEMA)';
                '40437  King Khaled Intl Arpt (OERK)';
                '40706  Tabriz (OITT)';
                '40745  Mashhad (OIMM)';
                '40754  Tehran-Mehrabad (OIII)';
                '40766  Kermanshah (OICC)';
                '40800  Esfahan (OIFM)';
                '40848  Shiraz (OISS)';
                '40856  Zahedan (OIZH)';
                '40948  Kabul Airport (OAKB)';
                '41024  Jeddah (King Abdul Aziz) (OEJN)';
                '41112  Abha (OEAB)';
                '41217  Abu Dhabi Inter Arpt (OMAA)';
                '41256  Seeb; Intl Airport (OOMS)';
                '41316  Salalah (OOSA)';
                '41624  Dera Ismail Khan (OPDI)';
                '41715  Jacobabad (OPJA)';
                '41718  Khanpur';
                '41749  Nawabshah (OPNH)';
                '60018  Guimar-Tenerife';
                '60155  Casablanca (GMMC)';
                '60390  Dar-El-Beida (DAAG)';
                '60571  Bechar (DAOR)';
                '60680  Tamanrasset';
                '60715  Tunis-Carthage (DTTA)';
                '60760  Tozeur (DTTZ)';
                '63741  Nairobi/Dagoretti (HKNC)';
                '64400  Pointe-Noire (FCPP)';
                '64450  Brazzaville/Maya-Maya (FCBB)';
                '65578  Abidjan (DIAP)';
                '67083  Antananarivo/Ivato (FMMI)';
                '68263  Pretoria (Irene) (FAIR)';
                '68424  Upington (FAUP)';
                '68842  Port Elizabeth (FAPE)';
                '12374  Legionowo';
                '13388  Nis (LYNI)';
                '15420  Bucuresti Inmh-Banesa (LRBS)';
                '16622  Thessaloniki (Airport) (LGTS)';
                '16754  Heraklion (Airport) (LGIR)';
                '17030  Samsun';
                '17064  Istanbul/Kartal';
                '17095  Erzurum (ERZM)';
                '17130  Ankara/Central';
                '17196  Kayseri (LTAU)';
                '17220  Izmir/Guzelyali';
                '17240  Isparta (LTBM)';
                '17281  Diyarbakir';
                '17351  Adana/Bolge';
                '27962  Penza; MS (UWPP)';
                '27995  Samara (Bezencuk); MS';
                '28951  Kostanai; AK';
                '33041  Gomel; MI';
                '33345  Kyiv; KI (UKKK)';
                '33393  Lviv; KI (UKLL)';
                '34009  Kursk; MS';
                '34122  Voronez; MS (UUOO)';
                '34172  Saratov; MS';
                '34247  Kalac; MS';
                '34300  Kharkiv; KI (UKHH)';
                '34467  Volgograd; TB (URWW)';
                '34731  Rostov-Na-Donu; TB (URRR)';
                '34858  Divnoe; TB';
                '34882  Astrakhan; TB';
                '35121  Orenburg; AL';
                '35229  Aktjubinsk; AL (UATT)';
                '35700  Atyran; AL';
                '37011  Tuapse; TB';
                '37055  Mineralnye Vody; TB (URMM)';
                '37259  Mahachkala';
                '37789  Yerevan/Yerevan-Arabkir; TB (UGEE)';
                '40179  Bet Dagan';
                '40265  Mafraq (OJMF)';
                '40373  Al-Qaisumah (OEPA)';
                '40375  Tabuk (OETB)';
                '40394  Hail (OEHL)';
                '40417  K.F.I.A.-Dammam (OEDF)';
                '40430  Al-Madinah (OEMA)';
                '40437  King Khaled Intl Arpt (OERK)';
                '40706  Tabriz (OITT)';
                '40745  Mashhad (OIMM)';
                '40754  Tehran-Mehrabad (OIII)';
                '40766  Kermanshah (OICC)';
                '40800  Esfahan (OIFM)';
                '40848  Shiraz (OISS)';
                '40856  Zahedan (OIZH)';
                '41024  Jeddah (King Abdul Aziz) (OEJN)';
                '41112  Abha (OEAB)';
                '41217  Abu Dhabi Inter Arpt (OMAA)';
                '41256  Seeb; Intl Airport (OOMS)';
                '41316  Salalah (OOSA)';
                '29839  Barnaul; IR';
                '29862  Hakasskaja; NO';
                '30635  Ust-Barguzin; IR';
                '30715  Angarsk; IR';
                '30758  Chita; IR (UIAA)';
                '30935  Krasnyj Chikoj; IR';
                '30965  Borzja; IR';
                '31510  Blagovescensk; HA';
                '31538  Sutur';
                '31736  Habarovsk; HA';
                '31873  Dalnerechensk; HA';
                '31977  Vladivostok (Sad Gorod); HA';
                '32150  Juzhno-Sahalinsk; HA (UHSS)';
                '35394  Karaganda; AL';
                '35671  Zhezkazgan; AL';
                '36003  Pavlodar; AL';
                '36096  Kyzyl; NO';
                '36872  Almaty; AL';
                '38064  Kyzylorda';
                '38341  Zhambyl; AL';
                '40745  Mashhad (OIMM)';
                '40856  Zahedan (OIZH)';
                '40948  Kabul Airport (OAKB)';
                '41624  Dera Ismail Khan (OPDI)';
                '41715  Jacobabad (OPJA)';
                '41718  Khanpur';
                '41749  Nawabshah (OPNH)';
                '41859  Rangpur';
                '41883  Bogra';
                '41891  Sylhet (VGSY)';
                '41907  Ishurdi (VGIS)';
                '41923  Dhaka (VGTJ)';
                '41936  Jessore (VGJR)';
                '41943  Feni';
                '41950  Barisal';
                '41977  Chittagong (Ambagan)';
                '41992  Coxs Bazar (VGCB)';
                '42027  Srinagar';
                '42101  Patiala';
                '42182  New Delhi/Safdarjung (VIDD)';
                '42299  Gangtok';
                '42314  Dibrugarh/Mohanbari (VEMN)';
                '42339  Jodhpur (VIJO)';
                '42361  Gwalior (VIGR)';
                '42369  Lucknow/Amausi (VILK)';
                '42379  Gorakhpur (VEGK)';
                '42410  Gauhati (VEGT)';
                '42492  Patna (VEPT)';
                '42623  Imphal (VEIM)';
                '42647  Ahmadabad (VAAH)';
                '42667  Bhopal/Bairagarh (VABP)';
                '42701  M.O. Ranchi (VERC)';
                '42724  Agartala (VEAT)';
                '42809  Calcutta/Dum Dum (VECC)';
                '42867  Nagpur Sonegaon (VANP)';
                '42874  Pbo Raipur';
                '42971  Bhubaneswar (VEBS)';
                '43003  Bombay/Santacruz (VABB)';
                '43014  Aurangabad Chikalthan (VAAU)';
                '43041  Jagdalpur';
                '43128  Hyderabad Airport (VOHY)';
                '43150  Vishakhapatnam/Waltair';
                '43192  Goa/Panjim';
                '43279  Madras/Minambakkam (VOMM)';
                '43285  Mangalore/Panambur';
                '43295  Bangalore';
                '43311  Amini Divi';
                '43333  Port Blair (VEPB)';
                '43346  Karaikal';
                '43353  Cochin/Willingdon (VOCC)';
                '43369  Minicoy';
                '43371  Thiruvananthapuram';
                '43466  Colombo';
                '43497  Hambantota';
                '44231  Muren';
                '44292  Ulaan-Baator';
                '44373  Dalanzadgad';
                '45004  Kings Park';
                '47102  Baengnyeongdo';
                '47104  Bukgangneung';
                '47122  Osan Ab (RKSO)';
                '47138  Pohang';
                '47158  Kwangju Ab (RKJJ)';
                '47169  Heuksando';
                '47186  National Typhoon Centre';
                '47401  Wakkanai';
                '47412  Sapporo';
                '47418  Kushiro';
                '47580  Misawa Ab (RJSM)';
                '47582  Akita';
                '47600  Wajima';
                '47646  Tateno';
                '47678  Hachijyojima/Omure';
                '47741  Matsue';
                '47778  Shionomisaki';
                '47807  Fukuoka';
                '47827  Kagoshima';
                '47909  Naze/Funchatoge';
                '47918  Ishigakijima (ROIG)';
                '47945  Minamidaitojima (ROMD)';
                '48453  Bangna';
                '48477  Sattahip';
                '48500  Prachuap Khirikhan (VTBP)';
                '48601  Penang/Bayan Lepas (WMKP)';
                '48615  Kota Bharu (WMKC)';
                '48650  Sepang';
                '48657  Kuantan (WMKD)';
                '48698  Singapore/Changi Arpt (WSSS)';
                '48811  Dien Bien Phu';
                '48820  Ha Noi (VVNB)';
                '48839  Bach Longvi';
                '48845  Vinh (VVVH)';
                '48855  Da Nang (VVDN)';
                '48870  Qui Nhon';
                '48877  Nha Trang';
                '48887  Phan Thiet';
                '48900  Ho Chi Minh (VVTS)';
                '48914  Ca Mau';
                '50527  Hailar; SY';
                '50557  Nenjiang; SY';
                '50774  Yichun; SY';
                '50953  Harbin; SY';
                '51076  Altay; UQ';
                '51431  Yining; UQ (ZWYN)';
                '51463  Urumqi; UQ';
                '51644  Kuqa; UQ';
                '51709  Kashi; UQ (ZWSH)';
                '51777  Ruoqiang; UQ';
                '51828  Hotan; UQ (ZWTN)';
                '51839  Minfeng; UQ';
                '52203  Hami; UQ (ZWHM)';
                '52267  Ejin Qi; LZ';
                '52323  Mazong Shan; LZ';
                '52418  Dunhuang; LZ';
                '52533  Jiuquan; LZ (ZLJQ)';
                '52681  Minqin; LZ';
                '52818  Golmud; LZ';
                '52836  Dulan; LZ';
                '52866  Xining; LZ (ZLXN)';
                '52983  Yu Zhong; LZ';
                '53068  Erenhot; BJ';
                '53463  Hohhot; BJ (ZBHH)';
                '53513  Linhe; BJ';
                '53614  Yinchuan; LZ (ZLIC)';
                '53772  Taiyuan; BJ (ZBYN)';
                '53845  Yan An; LZ (ZLYA)';
                '53915  Pingliang; LZ';
                '54102  Xilin Hot; BJ';
                '54135  Tongliao; SY';
                '54161  Changchun; SY (ZYCC)';
                '54218  Chifeng; SY';
                '54292  Yanji; SY';
                '54374  Linjiang; SY';
                '54511  Beijing; BJ (ZBAA)';
                '54662  Dalian; BJ (ZYTL)';
                '54727  Zhangqiu';
                '54857  Qingdao; BJ (ZSQD)';
                '55299  Nagqu; CD';
                '55591  Lhasa; CD (ZULS)';
                '56029  Yushu; LZ';
                '56080  Hezuo; LZ';
                '56137  Qamdo; CD';
                '56146  Garze; CD';
                '56187  Wenjiang';
                '56571  Xichang; CD';
                '56691  Weining; CD';
                '56739  Tengchong; CD';
                '56778  Kunming; CD (ZPPP)';
                '56964  Simao; CD';
                '56985  Mengzi; CD';
                '57083  Zhengzhou; BJ (ZHCC)';
                '57127  Hanzhong; LZ';
                '57131  Jinghe';
                '57178  Nanyang; BJ';
                '57447  Enshi; HK';
                '57461  Yichang; HK';
                '57494  Wuhan; HK (ZHHH)';
                '57516  Chongqing; CD (ZUCK)';
                '57749  Huaihua; HK';
                '57816  Guiyang; CD (ZUGY)';
                '57957  Guilin; GZ (ZGKL)';
                '57972  Chenzhou; HK';
                '57993  Ganzhou; HK (ZSGZ)';
                '58027  Xuzhou; SH';
                '58150  Sheyang; SH';
                '58203  Fuyang; HK';
                '58238  Nanjing; SH (ZSNJ)';
                '58362  Shanghai; SH';
                '58424  Anqing; HK';
                '58457  Hangzhou; SH (ZSHC)';
                '58606  Nanchang; HK (ZSCN)';
                '58633  Qu Xian; SH';
                '58665  Hongjia; SH';
                '58725  Shaowu; SH';
                '58847  Fuzhou; SH (ZSFZ)';
                '58968  Taibei; SH';
                '59134  Xiamen; SH (ZSAM)';
                '59211  Baise; GZ';
                '59265  Wuzhou; GZ';
                '59280  Qing Yuan; GZ';
                '59316  Shantou; GZ (ZGOW)';
                '59431  Nanning; GZ (ZGNN)';
                '59758  Haikou; GZ (ZGHK)';
                '59981  Xisha Dao; GZ';
                '96009  Lhokseumawe/Malikussaleh (WITM)';
                '96015  Meulaboh/Cut Nyak Dhien (WITC)';
                '96035  Medan/Polonia (WIMM)';
                '96091  Tanjung Pinang/Kijang (WIKN)';
                '96109  Pekan Baru/Simpangtiga (WIBB)';
                '96147  Ranai (WION)';
                '96163  Padang/Tabing (WIMG)';
                '96171  Rengat/Japura (WIPR)';
                '96179  Singkep/Dabo (WIKS)';
                '96195  Jambi/Sultan Taha (WIPA)';
                '96221  Palembang/St. Badarudin (WIPP)';
                '96237  Pangkal Pinang (WIKK)';
                '96249  Tanjung Pandan/Buluh (WIKD)';
                '96253  Bengkulu/Padang Kemiling (WIPL)';
                '96295  Tanjung Karang/Radin (WIIT)';
                '96315  Brunei Airport (WBSB)';
                '96413  Kuching (WBGG)';
                '96441  Bintulu (WBGB)';
                '96471  Kota Kinabalu (WBKK)';
                '96481  Tawau (WBKW)';
                '96509  Tarakan/Juwata (WRLR)';
                '96535  Paloh';
                '96581  Pontianak/Supadio (WIOO)';
                '96633  Balikpapan/Sepinggan (WRLL)';
                '96645  Pangkalan Bun/Iskandar (WRBI)';
                '96685  Banjarmasin/Syamsudin (WRBB)';
                '96737  Serang';
                '96739  Curug/Budiarto (WIIA)';
                '96749  Jakarta/Soekarno-Hatta (WIII)';
                '96791  Cirebon/Jatiwangi';
                '96797  Tegal';
                '96805  Cilacap (WIIL)';
                '96839  Semarang/Ahmad Yani (WIIS)';
                '96935  Surabaya/Juanda (WRSJ)';
                '96987  Banyuwangi';
                '97014  Menado/ Sam Ratulangi (WAMM)';
                '97028  Toli-Toli/Lalos (WAMI)';
                '97048  Gorontalo/Jalaluddin (WAMG)';
                '97072  Palu/Mutiara (WAML)';
                '97086  Luwuk/Bubung (WAMW)';
                '97180  Ujung Pandang/Hasanuddin (WAAA)';
                '97230  Denpasar/Ngurah Rai (WRRR)';
                '97240  Mataram/Selaparang (WRRA)';
                '97260  Sumbawa Besar/Brangbiji (WRRS)';
                '97270  Bima/M.Salahuddin (WRRB)';
                '97300  Maumere/Wai Oti (WRKC)';
                '97340  Waingapu/Mau Hau (WRRW)';
                '97372  Kupang/Eltari (WRKK)';
                '97430  Ternate/Babullah (WAMT)';
                '97600  Sanana (WAPN)';
                '98223  Laoag (RPLI)';
                '98328  Baguio (RPUB)';
                '98433  Tanay';
                '98444  Legaspi (RPMP)';
                '98618  Puerto Princesa (RPVP)';
                '98646  Mactan (RPMT)';
                '98747  Cagayan De Oro';
                '98753  Davao Airport (RPMD)';
                '01028  Bjornoya (ENBJ)';
                '01415  Stavanger/Sola (ENZV)';
                '03005  Lerwick';
                '03238  Albermarle';
                '03354  Nottingham/Watnall';
                '03882  Herstmonceux';
                '04018  Keflavikurflugvollur (BIKF)';
                '04220  Aasiaat (Egedesminde) (BGEM)';
                '04270  Narsarsuaq (BGBW)';
                '04320  Danmarkshavn (BGDH)';
                '04339  Ittoqqortoormiit (BGSC)';
                '04360  Tasiilaq (Ammassalik) (BGAM)';
                '04417  Geosummit';
                '06011  Torshavn';
                '10035  Schleswig';
                '10113  Norderney';
                '10184  Greifswald';
                '10393  Lindenberg';
                '12120  Leba';
                '12374  Legionowo';
                '20046  Polargmo Im. Krenkelja; DK';
                '20292  Gmo Im.E.K. Fedorova; DK';
                '20674  Ostrov Dikson; DK';
                '20744  Malye Karmakuly; DK';
                '21432  Ostrov Kotelnyj; DK';
                '21824  Tiksi; TK';
                '21946  Chokurdah; TK';
                '22008  Murmansk';
                '22217  Kandalaksa; AR';
                '22271  Sojna; AR';
                '22522  Kem; LE';
                '22543  Arhangelsk; AR';
                '22820  Petrozavodsk; LE';
                '22845  Kargopol; AR';
                '23078  Norilsk';
                '23205  Narjan-Mar; AR';
                '23330  Salehard; NO';
                '23415  Pechora; AR';
                '23472  Turuhansk; NO';
                '23802  Syktyvkar; AR';
                '23884  Bor; NO';
                '23921  Ivdel; SV';
                '23933  Hanty-Mansijsk; NO (USHH)';
                '23955  Aleksandrovskoe; NO';
                '24122  Olenek; HA';
                '24266  Verhojansk; HA';
                '24343  Zhigansk; HA';
                '24507  Tura; NO';
                '24641  Viljujsk; HA';
                '24688  Ojmjakon; HA';
                '24726  Mirnvy; HB';
                '24908  Vanavara; NO';
                '24947  Olekminsk; HA';
                '24959  Jakutsk; HA (UEEE)';
                '25123  Cherskij; HA';
                '25403  Zyrjanka; HA';
                '25428  Omolon; HB';
                '25703  Sejmchan; HA';
                '25913  Magadan; HA (UHMM)';
                '26038  Tallinn; LE (ULTT)';
                '26075  St.Petersburg(Voejkovo); St.Petersburg(Voejkovo); LE (ULLI)';
                '26298  Bologoe; LE';
                '26702  Kaliningrad; MI';
                '26781  Smolensk; MI';
                '27038  Vologda; AR (ULWW)';
                '27199  Kirov; MS';
                '27459  Niznij Novgorod; MS';
                '27594  Kazan (Vyzaovye); MS';
                '27713  Moskva (Dolgoprudnyj); MS';
                '27707  Suhinici; MS';
                '27730  Rjazan; MS';
                '27962  Penza; MS (UWPP)';
                '27995  Samara (Bezencuk); MS';
                '28225  Perm; SV';
                '28275  Tobolsk; NO';
                '28445  Verhnee Dubrovo; SV';
                '28661  Kurgan; SV';
                '28722  Ufa; SV';
                '28951  Kostanai; AK';
                '29231  Kolpasevo; NO';
                '29263  Enisejsk; NO (UNII)';
                '29282  Bogucany; NO';
                '29572  Emeljanovo; NO';
                '29612  Barabinsk; NO';
                '29634  Novosibirsk; NO (UNNN)';
                '29698  Nizhneudinsk; NO (UINN)';
                '29839  Barnaul; IR';
                '29862  Hakasskaja; NO';
                '30054  Vitim; HA';
                '30309  Bratsk; IR';
                '30372  Chara; IR';
                '30635  Ust-Barguzin; IR';
                '30715  Angarsk; IR';
                '30758  Chita; IR (UIAA)';
                '30965  Borzja; IR';
                '31004  Aldan; HA';
                '31088  Ohotsk; HA';
                '31168  Ajan; HA';
                '31369  Nikolaevsk-Na-Amure; HA';
                '31510  Blagovescensk; HA';
                '31538  Sutur';
                '31736  Habarovsk; HA';
                '31770  Sovetskay Gavan';
                '31873  Dalnerechensk; HA';
                '31977  Vladivostok (Sad Gorod); HA';
                '32061  Aleksandrovsk-Sahalnskij; HA';
                '32098  Poronajsk; HA';
                '32150  Juzhno-Sahalinsk; HA (UHSS)';
                '32215  Severo-Kurilsk; HA';
                '32477  Sobolevo';
                '32540  Kamchatskij; HA (UHPP)';
                '32618  Nikolskoe; HA';
                '33041  Gomel; MI';
                '33345  Kyiv; KI (UKKK)';
                '33393  Lviv; KI (UKLL)';
                '34009  Kursk; MS';
                '34122  Voronez; MS (UUOO)';
                '34172  Saratov; MS';
                '34247  Kalac; MS';
                '34300  Kharkiv; KI (UKHH)';
                '34467  Volgograd; TB (URWW)';
                '34731  Rostov-Na-Donu; TB (URRR)';
                '34858  Divnoe; TB';
                '34882  Astrakhan; TB';
                '35121  Orenburg; AL';
                '35229  Aktjubinsk; AL (UATT)';
                '35700  Atyran; AL';
                '36003  Pavlodar; AL';
                '37011  Tuapse; TB';
                '37055  Mineralnye Vody; TB (URMM)';
                '37259  Mahachkala';
                '37789  Yerevan/Yerevan-Arabkir; TB (UGEE)';
                '47401  Wakkanai';
                '50527  Hailar; SY';
                '50557  Nenjiang; SY';
                '50774  Yichun; SY';
                '50953  Harbin; SY';
                '70026  Barrow/W. Post W.Rogers; AK (PABR)';
                '70133  Kotzebue; Ralph Wien; AK (PAOT)';
                '70200  Nome; AK (PAOM)';
                '70219  Bethel/Bethel Airport; AK (PABE)';
                '70231  Mcgrath; AK (PAMC)';
                '70261  Fairbanks/Int; AK (PAFA)';
                '70273  Anchorage/Int; AK (PANC)';
                '70308  St. Paul; AK (PASN)';
                '70326  King Salmon; AK (PAKN)';
                '70350  Kodiak; AK (PADQ)';
                '70361  Yakutat; AK (PAYA)';
                '70398  Annette Island; AK (PANT)';
                '71043  Norman Wells Ua; NT (YVQ)';
                '71081  Hall Beach; NT (YUX)';
                '71082  Alert; NT (WLT)';
                '71109  Port Hardy; BC (YZT)';
                '71119  Edmonton Stony Plain; AB (WSE)';
                '71126  Edmonton; AB (ZED)';
                '71802  Mt Pearl; NF (AYT)';
                '71811  Sept-Iles; QB (YZV)';
                '71815  Stephenville; NF (YJT)';
                '71816  Goose Bay; NF (YYR)';
                '71823  La Grande Iv; QB (YAH)';
                '71867  The Pas; MB (YQD)';
                '71906  Kuujjuaq; QB (YVP)';
                '71907  Inukjuak; QB (WPH)';
                '71908  Prince George; BC (ZXS)';
                '71909  Iqaluit; NT (YFB)';
                '71913  Churchill; MB (YYQ)';
                '71917  Eureka; NT (WEU)';
                '71924  Resolute; NT (YRB)';
                '71925  Cambridge Bay; NT (YCB)';
                '71926  Baker Lake; NT (YBK)';
                '71934  Fort Smith; NT (YSM)';
                '71945  Fort Nelson; BC (YYE)';
                '71957  Inuvik; NT (YEV)';
                '72797  Quillayute; WA (UIL)';
                '73033  Vernon; BC (CWVK)';
                '89009  Amundsen-Scott';
                '89532  Syowa';
                '89571  Davis';
                '89592  Mirnyj';
                '89611  Casey';
                '89664  Mcmurdo';
                '93112  Whenuapai (NZWP)';
                '93417  Paraparaumu Aerodrome (NZPP)';
                '93844  Invercargill Aerodrome (NZNV)';
                '45004  Kings Park';
                '47909  Naze/Funchatoge';
                '47918  Ishigakijima (ROIG)';
                '47945  Minamidaitojima (ROMD)';
                '47971  Chichijima (RJAO)';
                '47991  Minamitorishima (RJAM)';
                '48839  Bach Longvi';
                '48855  Da Nang (VVDN)';
                '48870  Qui Nhon';
                '48877  Nha Trang';
                '48887  Phan Thiet';
                '48900  Ho Chi Minh (VVTS)';
                '57516  Chongqing; CD (ZUCK)';
                '57749  Huaihua; HK';
                '57816  Guiyang; CD (ZUGY)';
                '57957  Guilin; GZ (ZGKL)';
                '57972  Chenzhou; HK';
                '57993  Ganzhou; HK (ZSGZ)';
                '58606  Nanchang; HK (ZSCN)';
                '58633  Qu Xian; SH';
                '58665  Hongjia; SH';
                '58725  Shaowu; SH';
                '58847  Fuzhou; SH (ZSFZ)';
                '58968  Taibei; SH';
                '59134  Xiamen; SH (ZSAM)';
                '59211  Baise; GZ';
                '59265  Wuzhou; GZ';
                '59280  Qing Yuan; GZ';
                '59316  Shantou; GZ (ZGOW)';
                '59431  Nanning; GZ (ZGNN)';
                '59758  Haikou; GZ (ZGHK)';
                '59981  Xisha Dao; GZ';
                '91165  Lihue; HI (PHLI)';
                '91212  Guam Intl Arpt (PGAC)';
                '91285  Hilo/Gen; HI (PHTO)';
                '91334  Truk (PTKK)';
                '91348  Ponape (PTPN)';
                '91366  Kwajalein/Bucholz (PKWA)';
                '91376  Majuro/Marshall Is (PKMJ)';
                '91408  Koror; Palau Is (PTRO)';
                '91413  Yap (PTYA)';
                '91610  Tarawa (NGTA)';
                '91643  Funafuti (NGFU)';
                '91765  Pago Pago/Int.Airp. (NSTU)';
                '94120  Darwin Airport; NT (YPDN)';
                '94203  Broome Amo; WE (YBRM)';
                '94294  Townsville Aero; QU (YBTL)';
                '94299  Willis Island; QU';
                '94312  Port Hedland Amo; WE (YPPD)';
                '94326  Alice Springs Aero; NT (YBAS)';
                '94332  Mount Isa Amo; QU (YBMA)';
                '94374  Rockhampton Aero; QU (YBRK)';
                '94403  Geraldton Amo; WE (YPGN)';
                '94430  Meekatharra Amo; WE (YPMR)';
                '94461  Giles; WE';
                '94510  Charleville Amo; QU (YBCV)';
                '94578  Brisbane Airport Aero; QU (YBBN)';
                '94610  Perth Airport; WE (YPPH)';
                '94638  Esperance Mo; WE';
                '94653  Ceduna Amo; SA (YPCD)';
                '94659  Woomera Aerodrome Mo; SA (YPWR)';
                '94672  Adelaide Airport; SA (YPAD)';
                '94711  Cobar Mo; NW';
                '94776  Williamtown Amo Raaf; NW (YSWM)';
                '94802  Albany Airport; WE (YPAL)';
                '94866  Melbourne Airport; VC (YMML)';
                '94910  Wagga Wagga Amo/Aws; NW (YSWG)';
                '94975  Hobart Airport; TA (YMHB)';
                '94995  Lord Howe Island';
                '94996  Norfolk Island Aero (YSNF)';
                '95527  Moree Mo; NW';
                '96147  Ranai (WION)';
                '96249  Tanjung Pandan/Buluh (WIKD)';
                '96315  Brunei Airport (WBSB)';
                '96413  Kuching (WBGG)';
                '96441  Bintulu (WBGB)';
                '96471  Kota Kinabalu (WBKK)';
                '96481  Tawau (WBKW)';
                '96509  Tarakan/Juwata (WRLR)';
                '96535  Paloh';
                '96581  Pontianak/Supadio (WIOO)';
                '96633  Balikpapan/Sepinggan (WRLL)';
                '96645  Pangkalan Bun/Iskandar (WRBI)';
                '96685  Banjarmasin/Syamsudin (WRBB)';
                '96739  Curug/Budiarto (WIIA)';
                '96749  Jakarta/Soekarno-Hatta (WIII)';
                '96791  Cirebon/Jatiwangi';
                '96797  Tegal';
                '96805  Cilacap (WIIL)';
                '96839  Semarang/Ahmad Yani (WIIS)';
                '96935  Surabaya/Juanda (WRSJ)';
                '96987  Banyuwangi';
                '97014  Menado/ Sam Ratulangi (WAMM)';
                '97028  Toli-Toli/Lalos (WAMI)';
                '97048  Gorontalo/Jalaluddin (WAMG)';
                '97072  Palu/Mutiara (WAML)';
                '97086  Luwuk/Bubung (WAMW)';
                '97180  Ujung Pandang/Hasanuddin (WAAA)';
                '97230  Denpasar/Ngurah Rai (WRRR)';
                '97240  Mataram/Selaparang (WRRA)';
                '97260  Sumbawa Besar/Brangbiji (WRRS)';
                '97270  Bima/M.Salahuddin (WRRB)';
                '97300  Maumere/Wai Oti (WRKC)';
                '97340  Waingapu/Mau Hau (WRRW)';
                '97372  Kupang/Eltari (WRKK)';
                '97430  Ternate/Babullah (WAMT)';
                '97502  Sorong/Jefman (WASS)';
                '97560  Biak/Frans Kaisiepo (WABB)';
                '97600  Sanana (WAPN)';
                '97724  Ambon/Pattimura (WAPP)';
                '97748  Geser';
                '97810  Tual/Dumatubun';
                '97900  Saumlaki/Olilit (WAPI)';
                '97980  Merauke/Mopah (WAKK)';
                '98223  Laoag (RPLI)';
                '98328  Baguio (RPUB)';
                '98433  Tanay';
                '98444  Legaspi (RPMP)';
                '98618  Puerto Princesa (RPVP)';
                '98646  Mactan (RPMT)';
                '98747  Cagayan De Oro';
                '98753  Davao Airport (RPMD)';
                '78954  Grantley Adams (TBPB)';
                '78970  Piarco Int. Airport (TTPP)';
                '78988  Hato Airport; Curacao (TNCC)';
                '80001  San Andres Isl (SKSP)';
                '82022  Boa Vista (SBBV)';
                '82099  Macapa (SBMQ)';
                '82026  Tirios';
                '82107  Sao Gabriel Da Cachoeira;';
                '82193  Belem (Aeroporto) (SBBE)';
                '82244  Santarem;';
                '82281  Sao Luiz';
                '82332  Manaus (Aeroporto) (SBMN)';
                '82400  Fernando De Noronha (SBFN)';
                '82411  Tabatinga;';
                '82532  Manicore';
                '82599  Natal Aeroporto (SBNT)';
                '82824  Porto Velho (Aeroporto) (SBPV)';
                '82917  Rio Branco';
                '82965  Alta Floresta (Aero) (SBAT)';
                '83208  Vilhena (Aeroporto) (SBVH)';
                '83362  Cuiaba (Aeroporto) (SBCY)';
                '83378  Brasilia (Aeroporto) (SBBR)';
                '83525  Uberlandia (SBUL)';
                '83554  Corumba';
                '83566  Confis Intnl Arpt';
                '83612  Campo Grande (Aero) (SBCG)';
                '83649  Vitoria';
                '83746  Galeao (SBGL)';
                '83768  Londrina (SBLO)';
                '83779  Marte Civ/Mil (SBMT)';
                '83827  Foz Do Iguacu (Aero) (SBFI)';
                '83840  Curitiba (Aeroporto) (SBCT)';
                '83899  Florianopolis (SBFL)';
                '83928  Uruguaniana (SBUG)';
                '83937  Santa Maria (SBSM)';
                '83971  Porto Alegre (Aero) (SBPA)';
                '85586  Santo Domingo (SCSN)';
                '88889  Mount Pleasant Airport (EGYP)';
                '03953  Valentia Observatory';
                '04220  Aasiaat (Egedesminde) (BGEM)';
                '04270  Narsarsuaq (BGBW)';
                '04360  Tasiilaq (Ammassalik) (BGAM)';
                '32618  Nikolskoe; HA';
                '70026  Barrow/W. Post W.Rogers; AK (PABR)';
                '70133  Kotzebue; Ralph Wien; AK (PAOT)';
                '70200  Nome; AK (PAOM)';
                '70219  Bethel/Bethel Airport; AK (PABE)';
                '70231  Mcgrath; AK (PAMC)';
                '70261  Fairbanks/Int; AK (PAFA)';
                '70273  Anchorage/Int; AK (PANC)';
                '70308  St. Paul; AK (PASN)';
                '70326  King Salmon; AK (PAKN)';
                '70350  Kodiak; AK (PADQ)';
                '70361  Yakutat; AK (PAYA)';
                '70398  Annette Island; AK (PANT)';
                '71043  Norman Wells Ua; NT (YVQ)';
                '71081  Hall Beach; NT (YUX)';
                '71109  Port Hardy; BC (YZT)';
                '71119  Edmonton Stony Plain; AB (WSE)';
                '71126  Edmonton; AB (ZED)';
                '71603  Yarmouth; NS (YQI)';
                '71722  Maniwaki; QB (WMW)';
                '71802  Mt Pearl; NF (AYT)';
                '71811  Sept-Iles; QB (YZV)';
                '71815  Stephenville; NF (YJT)';
                '71816  Goose Bay; NF (YYR)';
                '71823  La Grande Iv; QB (YAH)';
                '71836  Moosonee; ON (YMO)';
                '71845  Pickle Lake; ON (WPL)';
                '71867  The Pas; MB (YQD)';
                '71906  Kuujjuaq; QB (YVP)';
                '71907  Inukjuak; QB (WPH)';
                '71908  Prince George; BC (ZXS)';
                '71909  Iqaluit; NT (YFB)';
                '71913  Churchill; MB (YYQ)';
                '71924  Resolute; NT (YRB)';
                '71925  Cambridge Bay; NT (YCB)';
                '71926  Baker Lake; NT (YBK)';
                '71934  Fort Smith; NT (YSM)';
                '71945  Fort Nelson; BC (YYE)';
                '71957  Inuvik; NT (YEV)';
                '72201  Key West/Int; FL (KEY)';
                '72202  Miami; FL (MFL)';
                '72206  Jacksonville Intl; FL (JAX)';
                '72208  Charleston/Muni; SC (CHS)';
                '72210  Tampa Bay Area; FL (TBW)';
                '72214  Tallahassee Fsu; FL (TLH)';
                '72215  Peachtree City; GA (FFC)';
                '72230  Shelby County Airport; AL (BMX)';
                '72233  Slidell Muni; LA (LIX)';
                '72235  Jackson Thompson Fld; MS (JAN)';
                '72240  Lake Charles/Muni; LA (LCH)';
                '72248  Shreveport Reg; LA (SHV)';
                '72249  Ft Worth; TX (FWD)';
                '72250  Brownsville Intl; TX (BRO)';
                '72251  Corpus Christi Intl; TX (CRP)';
                '72261  Del Rio/Int; TX (DRT)';
                '72265  Midland/Midland Reg; TX (MAF)';
                '72274  Tucson; AZ (TUS)';
                '72293  San Diego/Miramar; CA (NKX)';
                '72305  Newport; NC (MHX)';
                '72317  Greensboro/High Pt; NC (GSO)';
                '72318  Blacksburg; VA (RNK)';
                '72327  Nashville/Old Hickory; TN (BNA)';
                '72340  Little Rock/Adams; AR (LZK)';
                '72357  Norman/Westheimer; OK (OUN)';
                '72363  Amarillo Arpt(Awos); TX (AMA)';
                '72364  Santa Teresa; NM (EPZ)';
                '72365  Albuquerque/Int; NM (ABQ)';
                '72376  Flagstaff; AZ (FGZ)';
                '72388  Las Vegas; NV (VEF)';
                '72393  Vandenberg Afb; Vandenberg Afb; CA (VBG)';
                '72402  Wallops Island; VA (WAL)';
                '72403  Sterling; VA (IAD)';
                '72426  Wilmington; OH (ILN)';
                '72440  Springfield/Muno.; MO (SGF)';
                '72451  Dodge City(Awos); KS (DDC)';
                '72456  Topeka/Billard Muni; KS (TOP)';
                '72469  Denver/Stapleton; CO (DNR)';
                '72476  Grand Junction/Walker; CO (GJT)';
                '72489  Reno; NV (REV)';
                '72493  Oakland Int; CA (OAK)';
                '72501  Upton; NY (OKX)';
                '72518  Albany; NY (ALB)';
                '72520  Pittsburgh/Moon; PA (PIT)';
                '72528  Buffalo Int; NY (BUF)';
                '72558  Omaha/Valley; NE (OAX)';
                '72562  North Platte/Lee Bird; NE (LBF)';
                '72572  Salt Lake City/Intnl; UT (SLC)';
                '72582  Elko; NV (LKN)';
                '72597  Medford/Jackson; OR (MFR)';
                '72632  White Lake; MI (DTX)';
                '72634  Gaylord; MI (APX)';
                '72645  Green Bay/Straubel; WI (GRB)';
                '72649  Chanhassen; MN (MPX)';
                '72659  Aberdeen/Reg; SD (ABR)';
                '72662  Rapid City; SD (RAP)';
                '72672  Riverton; WY (RIW)';
                '72681  Boise/Mun; ID (BOI)';
                '72694  Salem/Mcnary; OR (SLE)';
                '72712  Caribou/Mun; ME (CAR)';
                '72747  Int.Falls/Falls Int; MN (INL)';
                '72764  Bismarck/Mun; ND (BIS)';
                '72768  Glasgow/Int; MT (GGW)';
                '72776  Great Falls; MT (TFX)';
                '72786  Spokane; WA (OTX)';
                '72797  Quillayute; WA (UIL)';
                '73033  Vernon; BC (CWVK)';
                '74389  Gray; ME (GYX)';
                '74455  Davenport; IA (DVN)';
                '74560  Lincoln; IL (ILX)';
                '74646  Lamont Oklahoma; OK (LMN)';
                '74794  Cape Kennedy; FL (XMR)';
                '76225  Chihuahua; Chih.';
                '76458  Colonia Juancarrasco';
                '76526  Guadalupe';
                '76612  Guadalajara; Jal.';
                '76654  Manzanillo; Col.';
                '76679  Aerop. Intl Mexico; D.F.';
                '78073  Nassau Airport (MYNN)';
                '78384  Owen Roberts Arpt (MWCR)';
                '78486  Santo Domingo (MDSD)';
                '78526  San Juan/Int (TJSJ)';
                '78583  Phillip Goldston Intl. (MZBZ)';
                '78954  Grantley Adams (TBPB)';
                '78970  Piarco Int. Airport (TTPP)';
                '78988  Hato Airport; Curacao (TNCC)';
                '80001  San Andres Isl (SKSP)';
                '82022  Boa Vista (SBBV)';
                '91165  Lihue; HI (PHLI)';
                '91285  Hilo/Gen; HI (PHTO)';
                '01415  Stavanger/Sola (ENZV)';
                '02365  Sundsvall-Harnosand Fpl';
                '02591  Visby Aerologiska Stn (ESQV)';
                '03005  Lerwick';
                '03808  Camborne';
                '03953  Valentia Observatory';
                '06011  Torshavn';
                '06458  Beauvecchain (EBBE)';
                '08579  Lisboa/Gago Coutinho';
                '10035  Schleswig';
                '10113  Norderney';
                '10184  Greifswald';
                '10393  Lindenberg';
                '10410  Essen (EDZE)';
                '10548  Meiningen';
                '10739  Stuttgart/Schnarrenberg';
                '10868  Muenchen-Oberschlssheim';
                '11520  Praha-Libus';
                '11747  Prostejov';
                '11952  Poprad-Ganovce';
                '12120  Leba';
                '12374  Legionowo';
                '12425  Wroclaw I';
                '12843  Budapest/Lorinc';
                '12982  Szeged (LHUD)';
                '13275  Beograd/Kosutnjak';
                '13388  Nis (LYNI)';
                '14240  Zagreb/Maksimir (LDDD)';
                '14430  Zadar';
                '15420  Bucuresti Inmh-Banesa (LRBS)';
                '15614  Sofia (Observ) (LBSF)';
                '16045  Rivolto (LIPI)';
                '16064  Novara/Cameri (LIMN)';
                '16113  Cuneo-Levaldigi';
                '16245  Pratica Di Mare (LIRE)';
                '16320  Brindisi (LIBR)';
                '16429  Trapani/Birgi (LICT)';
                '16546  Decimomannu (LIED)';
                '16716  Athinai (Airport) (LGAT)';
                '17030  Samsun';
                '17064  Istanbul/Kartal';
                '17130  Ankara/Central';
                '17196  Kayseri (LTAU)';
                '17220  Izmir/Guzelyali';
                '17240  Isparta (LTBM)';
                '17351  Adana/Bolge';
                '17516  Nicosia';
                '22820  Petrozavodsk, LE';
                '22845  Kargopol, AR';
                '26075  St.Petersburg(Voejkovo), St.Petersburg(Voejkovo), LE (ULLI)';
                '26298  Bologoe, LE';
                '26477  Velikie Luki, LE (ULOL)';
                '26702  Kaliningrad, MI';
                '26781  Smolensk, MI';
                '27038  Vologda, AR (ULWW)';
                '27199  Kirov, MS';
                '27459  Niznij Novgorod, MS';
                '27594  Kazan (Vyzaovye), MS';
                '27713  Moskva (Dolgoprudnyj), MS';
                '27707  Suhinici, MS';
                '27730  Rjazan, MS';
                '27962  Penza, MS (UWPP)';
                '27995  Samara (Bezencuk), MS';
                '33317  Shepetivka, KI';
                '33345  Kyiv, KI (UKKK)';
                '34009  Kursk, MS';
                '34122  Voronez, MS (UUOO)';
                '34172  Saratov, MS';
                '34247  Kalac, MS';
                '34467  Volgograd, TB (URWW)';
                '34731  Rostov-Na-Donu, TB (URRR)';
                '37011  Tuapse, TB';
                '01004  ENAS';
                '94998  Macquarie Island';
                '78016  Bermuda Nvl Stn Kindley (TXKF)';
                '91938  Tahiti-Faaa (NTAA)'};
            
            all_raob = sort(unique(all_raob));
            
            sta_id = find(contains(all_raob, sta_num), 1);
            name = sta_num;
            raob_code = sta_num;
            if ~isempty(sta_id)
                name = all_raob{sta_id}(8 : end);
                raob_code = all_raob{sta_id}(1:5);
            end
        end
        
        function raob_list = getRaobList()
            % Get a static manually inserted struct of all the Radiosondes sites available
            %
            % SYNTAX
            %   raob_list = Radiosondes.getRaobList();
            raob_list = struct();
            raob_list.s01004 = struct('lat',   78.91, 'lon',   11.93, 'name', 'ENAS Ny-Alesund Ii');
            raob_list.s01028 = struct('lat',   74.50, 'lon',   19.00, 'name', 'ENBJ Bjornoya');
            raob_list.s01415 = struct('lat',   58.87, 'lon',    5.67, 'name', 'ENZV Stavanger');
            raob_list.s03005 = struct('lat',   60.13, 'lon',   -1.18, 'name', 'Lerwick');
            raob_list.s03238 = struct('lat',   55.01, 'lon',   -1.52, 'name', 'Albermarle');
            raob_list.s03354 = struct('lat',   53.00, 'lon',   -1.25, 'name', 'Nottingham');            
            raob_list.s03882 = struct('lat',   50.90, 'lon',    0.32, 'name', 'Herstmonceux');
            raob_list.s03808 = struct('lat',   50.22, 'lon',   -5.32, 'name', 'Camborne');
            raob_list.s03953 = struct('lat',   51.93, 'lon',  -10.25, 'name', 'Valentia Observatory');
            raob_list.s04018 = struct('lat',   63.96, 'lon',  -22.60, 'name', 'BIKF Keflavikurflugvollur');
            raob_list.s04220 = struct('lat',   68.70, 'lon',  -52.85, 'name', 'BGEM Aasiaat (Egedesminde)');
            raob_list.s04270 = struct('lat',   61.15, 'lon',  -45.43, 'name', 'BGBW Narsarsuaq');
            raob_list.s04320 = struct('lat',   76.76, 'lon',  -18.66, 'name', 'BGDH Danmarkshavn');
            raob_list.s04339 = struct('lat',   70.48, 'lon',  -21.95, 'name', 'BGSC Ittoqqortoormiit');
            raob_list.s04360 = struct('lat',   65.60, 'lon',  -37.63, 'name', 'BGAM Tasiilaq (Ammassalik)');
            raob_list.s04417 = struct('lat',   72.57, 'lon',  -38.45, 'name', 'Geosummit');
            raob_list.s06011 = struct('lat',   62.01, 'lon',   -6.76, 'name', 'Torshavn');
            raob_list.s06458 = struct('lat',   50.75, 'lon',    4.77, 'name', 'EBBE Beauvecchain');
            raob_list.s08579 = struct('lat',   38.76, 'lon',   -9.13, 'name', 'Lisboa');
            raob_list.s10035 = struct('lat',   54.53, 'lon',    9.55, 'name', 'Schleswig');
            raob_list.s10113 = struct('lat',   53.71, 'lon',    7.15, 'name', 'Norderney');
            raob_list.s10184 = struct('lat',   54.10, 'lon',   13.40, 'name', 'Greifswald');
            raob_list.s10393 = struct('lat',   52.21, 'lon',   14.12, 'name', 'Lindenberg');
            raob_list.s10410 = struct('lat',   51.40, 'lon',    6.97, 'name', 'EDZE Essen');
            raob_list.s10548 = struct('lat',   50.56, 'lon',   10.38, 'name', 'Meiningen');
            raob_list.s10739 = struct('lat',   48.83, 'lon',    9.20, 'name', 'Stuttgart');
            raob_list.s10868 = struct('lat',   48.25, 'lon',   11.55, 'name', 'Muenchen-Oberschlssheim');
            raob_list.s11520 = struct('lat',   50.01, 'lon',   14.45, 'name', 'Praha-Libus');
            raob_list.s11747 = struct('lat',   49.45, 'lon',   17.13, 'name', 'Prostejov');
            raob_list.s11952 = struct('lat',   49.03, 'lon',   20.31, 'name', 'Poprad-Ganovce');
            raob_list.s12120 = struct('lat',   54.75, 'lon',   17.53, 'name', 'Leba');
            raob_list.s12374 = struct('lat',   52.40, 'lon',   20.96, 'name', 'Legionowo');
            raob_list.s12425 = struct('lat',   51.13, 'lon',   16.98, 'name', 'Wroclaw I');
            raob_list.s12843 = struct('lat',   47.43, 'lon',   19.18, 'name', 'Budapest');
            raob_list.s12982 = struct('lat',   46.25, 'lon',   20.10, 'name', 'LHUD Szeged');
            raob_list.s13275 = struct('lat',   44.76, 'lon',   20.42, 'name', 'Beograd');
            raob_list.s13388 = struct('lat',   43.33, 'lon',   21.90, 'name', 'LYNI Nis');
            raob_list.s14240 = struct('lat',   45.82, 'lon',   16.03, 'name', 'LDDD Zagreb');
            raob_list.s14430 = struct('lat',   44.10, 'lon',   15.34, 'name', 'Zadar');
            raob_list.s15420 = struct('lat',   44.50, 'lon',   26.13, 'name', 'LRBS Bucuresti Inmh-Banesa');
            raob_list.s15614 = struct('lat',   42.65, 'lon',   23.38, 'name', 'LBSF Sofia (Observ)');
            raob_list.s16045 = struct('lat',   45.97, 'lon',   13.05, 'name', 'LIPI Rivolto');
            raob_list.s16064 = struct('lat',   45.43, 'lon',    9.28, 'name', 'LIML Milano');
            raob_list.s16113 = struct('lat',   44.53, 'lon',    7.61, 'name', 'Cuneo-Levaldigi');
            raob_list.s16245 = struct('lat',   41.65, 'lon',   12.43, 'name', 'LIRE Pratica Di Mare');
            raob_list.s16320 = struct('lat',   40.65, 'lon',   17.95, 'name', 'LIBR Brindisi');
            raob_list.s16429 = struct('lat',   37.91, 'lon',   12.50, 'name', 'LICT Trapani');
            raob_list.s16546 = struct('lat',   39.35, 'lon',    8.85, 'name', 'LIED Decimomannu');
            raob_list.s16622 = struct('lat',   40.51, 'lon',   22.96, 'name', 'LGTS Thessaloniki (Airport)');
            raob_list.s16716 = struct('lat',   37.90, 'lon',   23.73, 'name', 'LGAT Athinai (Airport)');
            raob_list.s16754 = struct('lat',   35.33, 'lon',   25.18, 'name', 'LGIR Heraklion (Airport)');
            raob_list.s17030 = struct('lat',   41.28, 'lon',   36.30, 'name', 'Samsun');
            raob_list.s17064 = struct('lat',   40.90, 'lon',   29.15, 'name', 'Istanbul');
            raob_list.s17095 = struct('lat',   41.28, 'lon',   39.90, 'name', 'ERZM Erzurum');
            raob_list.s17196 = struct('lat',   38.69, 'lon',   35.50, 'name', 'LTAU Kayseri');
            raob_list.s17220 = struct('lat',   38.43, 'lon',   27.16, 'name', 'Izmir');
            raob_list.s17240 = struct('lat',   37.75, 'lon',   30.55, 'name', 'LTBM Isparta');
            raob_list.s17281 = struct('lat',   37.54, 'lon',   40.12, 'name', 'Diyarbakir');
            raob_list.s17351 = struct('lat',   36.98, 'lon',   35.35, 'name', 'Adana');
            raob_list.s17516 = struct('lat',   35.10, 'lon',   33.30, 'name', 'Nicosia');
            raob_list.s20046 = struct('lat',   80.61, 'lon',   58.05, 'name', 'Polargmo Im. Krenkelja');
            raob_list.s20292 = struct('lat',   77.71, 'lon',  104.30, 'name', 'Gmo Im.E.K. Fedorova');
            raob_list.s20674 = struct('lat',   73.50, 'lon',   80.40, 'name', 'Ostrov Dikson');
            raob_list.s20744 = struct('lat',   72.36, 'lon',   52.70, 'name', 'Malye Karmakuly');
            raob_list.s21432 = struct('lat',   76.00, 'lon',  137.86, 'name', 'Ostrov Kotelnyj');
            raob_list.s21824 = struct('lat',   71.58, 'lon',  128.91, 'name', 'Tiksi');
            raob_list.s21946 = struct('lat',   70.61, 'lon',  147.88, 'name', 'Chokurdah');
            raob_list.s22008 = struct('lat',   68.10, 'lon',   33.11, 'name', 'Murmansk');
            raob_list.s22217 = struct('lat',   67.15, 'lon',   32.35, 'name', 'Kandalaksa');
            raob_list.s22271 = struct('lat',   67.88, 'lon',   44.13, 'name', 'Sojna');
            raob_list.s22522 = struct('lat',   64.95, 'lon',   34.65, 'name', 'Kem');
            raob_list.s22543 = struct('lat',   64.62, 'lon',   40.51, 'name', 'Arhangelsk');
            raob_list.s22820 = struct('lat',   61.81, 'lon',   34.26, 'name', 'Petrozavodsk');
            raob_list.s22845 = struct('lat',   61.50, 'lon',   38.93, 'name', 'Kargopol');
            raob_list.s23078 = struct('lat',   69.32, 'lon',   88.22, 'name', 'Norilsk');
            raob_list.s23205 = struct('lat',   67.63, 'lon',   53.03, 'name', 'Narjan-Mar');
            raob_list.s23330 = struct('lat',   66.53, 'lon',   66.66, 'name', 'Salehard');
            raob_list.s23415 = struct('lat',   65.12, 'lon',   57.10, 'name', 'Pechora');
            raob_list.s23472 = struct('lat',   65.78, 'lon',   87.93, 'name', 'Turuhansk');
            raob_list.s23802 = struct('lat',   61.68, 'lon',   50.78, 'name', 'Syktyvkar');
            raob_list.s23884 = struct('lat',   61.60, 'lon',   90.01, 'name', 'Bor');
            raob_list.s23921 = struct('lat',   60.68, 'lon',   60.45, 'name', 'Ivdel');
            raob_list.s23933 = struct('lat',   61.01, 'lon',   69.03, 'name', 'USHH Hanty-Mansijsk');
            raob_list.s23955 = struct('lat',   60.43, 'lon',   77.86, 'name', 'Aleksandrovskoe');
            raob_list.s24122 = struct('lat',   68.50, 'lon',  112.43, 'name', 'Olenek');
            raob_list.s24266 = struct('lat',   67.56, 'lon',  133.40, 'name', 'Verhojansk');
            raob_list.s24343 = struct('lat',   66.76, 'lon',  123.40, 'name', 'Zhigansk');
            raob_list.s24507 = struct('lat',   64.26, 'lon',  100.23, 'name', 'Tura');
            raob_list.s24641 = struct('lat',   63.76, 'lon',  121.61, 'name', 'Viljujsk');
            raob_list.s24688 = struct('lat',   63.25, 'lon',  143.15, 'name', 'Ojmjakon');
            raob_list.s24726 = struct('lat',   62.53, 'lon',  113.86, 'name', 'Mirnvy');
            raob_list.s24908 = struct('lat',   60.33, 'lon',  102.26, 'name', 'Vanavara');
            raob_list.s24947 = struct('lat',   60.37, 'lon',  120.42, 'name', 'Olekminsk');
            raob_list.s24959 = struct('lat',   62.01, 'lon',  129.71, 'name', 'UEEE Jakutsk');
            raob_list.s25123 = struct('lat',   68.75, 'lon',  161.28, 'name', 'Cherskij');
            raob_list.s25403 = struct('lat',   65.72, 'lon',  150.90, 'name', 'Zyrjanka');
            raob_list.s25428 = struct('lat',   65.23, 'lon',  160.53, 'name', 'Omolon');
            raob_list.s25703 = struct('lat',   62.91, 'lon',  152.41, 'name', 'Sejmchan');
            raob_list.s25913 = struct('lat',   59.55, 'lon',  150.78, 'name', 'UHMM Magadan');
            raob_list.s26038 = struct('lat',   59.38, 'lon',   24.58, 'name', 'ULTT Tallinn');
            raob_list.s26075 = struct('lat',   59.95, 'lon',   30.70, 'name', 'ULLI St.Petersburg(Voejkovo)');
            raob_list.s26298 = struct('lat',   57.90, 'lon',   34.05, 'name', 'Bologoe');
            raob_list.s26477 = struct('lat',   56.35, 'lon',   30.61, 'name', 'Velikie Luki');
            raob_list.s26702 = struct('lat',   54.71, 'lon',   20.55, 'name', 'Kaliningrad');
            raob_list.s26781 = struct('lat',   54.75, 'lon',   32.06, 'name', 'Smolensk');
            raob_list.s27038 = struct('lat',   59.32, 'lon',   39.93, 'name', 'ULWW Vologda');
            raob_list.s27199 = struct('lat',   58.60, 'lon',   49.63, 'name', 'Kirov');
            raob_list.s27459 = struct('lat',   56.26, 'lon',   44.00, 'name', 'Niznij Novgorod');
            raob_list.s27594 = struct('lat',   55.82, 'lon',   48.52, 'name', 'Kazan (Vyzaovye)');
            raob_list.s27707 = struct('lat',   54.10, 'lon',   35.35, 'name', 'Suhinici');
            raob_list.s27713 = struct('lat',   55.93, 'lon',   37.52, 'name', 'Moskva (Dolgoprudnyj)');
            raob_list.s27730 = struct('lat',   54.63, 'lon',   39.70, 'name', 'Rjazan');
            raob_list.s27962 = struct('lat',   53.11, 'lon',   45.01, 'name', 'UWPP Penza');
            raob_list.s27995 = struct('lat',   52.98, 'lon',   49.43, 'name', 'Samara (Bezencuk)');
            raob_list.s28225 = struct('lat',   57.95, 'lon',   56.20, 'name', 'Perm');
            raob_list.s28275 = struct('lat',   58.15, 'lon',   68.25, 'name', 'Tobolsk');
            raob_list.s28445 = struct('lat',   56.73, 'lon',   61.06, 'name', 'Verhnee Dubrovo');
            raob_list.s28661 = struct('lat',   55.46, 'lon',   65.40, 'name', 'Kurgan');
            raob_list.s28722 = struct('lat',   54.71, 'lon',   55.83, 'name', 'Ufa');
            raob_list.s28951 = struct('lat',   53.23, 'lon',   63.62, 'name', 'Kostanai');
            raob_list.s29231 = struct('lat',   58.31, 'lon',   82.95, 'name', 'Kolpasevo');
            raob_list.s29263 = struct('lat',   58.45, 'lon',   92.15, 'name', 'UNII Enisejsk');
            raob_list.s29282 = struct('lat',   58.38, 'lon',   97.45, 'name', 'Bogucany');
            raob_list.s29572 = struct('lat',   56.18, 'lon',   92.61, 'name', 'Emeljanovo');
            raob_list.s29612 = struct('lat',   55.33, 'lon',   78.36, 'name', 'Barabinsk');
            raob_list.s29634 = struct('lat',   54.96, 'lon',   82.95, 'name', 'UNNN Novosibirsk');
            raob_list.s29698 = struct('lat',   54.88, 'lon',   99.03, 'name', 'UINN Nizhneudinsk');
            raob_list.s29839 = struct('lat',   53.35, 'lon',   83.81, 'name', 'Barnaul');
            raob_list.s29862 = struct('lat',   53.76, 'lon',   91.31, 'name', 'Hakasskaja');
            raob_list.s30054 = struct('lat',   59.45, 'lon',  112.58, 'name', 'Vitim');
            raob_list.s30309 = struct('lat',   56.28, 'lon',  101.75, 'name', 'Bratsk');
            raob_list.s30372 = struct('lat',   56.90, 'lon',  118.26, 'name', 'Chara');
            raob_list.s30635 = struct('lat',   53.41, 'lon',  109.01, 'name', 'Ust-Barguzin');
            raob_list.s30715 = struct('lat',   52.48, 'lon',  103.85, 'name', 'Angarsk');
            raob_list.s30758 = struct('lat',   52.08, 'lon',  113.48, 'name', 'UIAA Chita');
            raob_list.s30935 = struct('lat',   50.36, 'lon',  108.75, 'name', 'Krasnyj Chikoj');
            raob_list.s30965 = struct('lat',   50.40, 'lon',  116.51, 'name', 'Borzja');
            raob_list.s31004 = struct('lat',   58.61, 'lon',  125.36, 'name', 'Aldan');
            raob_list.s31088 = struct('lat',   59.36, 'lon',  143.20, 'name', 'Ohotsk');
            raob_list.s31168 = struct('lat',   56.45, 'lon',  138.15, 'name', 'Ajan');
            raob_list.s31369 = struct('lat',   53.15, 'lon',  140.70, 'name', 'Nikolaevsk-Na-Amure');
            raob_list.s31510 = struct('lat',   50.53, 'lon',  127.50, 'name', 'Blagovescensk');
            raob_list.s31538 = struct('lat',   50.07, 'lon',  132.13, 'name', 'Sutur');
            raob_list.s31736 = struct('lat',   48.53, 'lon',  135.23, 'name', 'Habarovsk');
            raob_list.s31770 = struct('lat',   49.00, 'lon',  140.27, 'name', 'Sovetskay Gavan');
            raob_list.s31873 = struct('lat',   45.86, 'lon',  133.73, 'name', 'Dalnerechensk');
            raob_list.s31977 = struct('lat',   43.26, 'lon',  132.05, 'name', 'Vladivostok (Sad Gorod)');
            raob_list.s32061 = struct('lat',   50.90, 'lon',  142.16, 'name', 'Aleksandrovsk-Sahalnskij');
            raob_list.s32098 = struct('lat',   49.22, 'lon',  143.10, 'name', 'Poronajsk');
            raob_list.s32150 = struct('lat',   46.95, 'lon',  142.71, 'name', 'UHSS Juzhno-Sahalinsk');
            raob_list.s32215 = struct('lat',   50.68, 'lon',  156.13, 'name', 'Severo-Kurilsk');
            raob_list.s32477 = struct('lat',   54.30, 'lon',  155.97, 'name', 'Sobolevo');
            raob_list.s32540 = struct('lat',   53.08, 'lon',  158.58, 'name', 'UHPP Kamchatskij');
            raob_list.s32618 = struct('lat',   55.20, 'lon',  165.98, 'name', 'Nikolskoe');
            raob_list.s33041 = struct('lat',   52.40, 'lon',   30.95, 'name', 'Gomel');
            raob_list.s33317 = struct('lat',   50.16, 'lon',   27.03, 'name', 'Shepetivka');
            raob_list.s33345 = struct('lat',   50.40, 'lon',   30.56, 'name', 'UKKK Kyiv');
            raob_list.s33393 = struct('lat',   49.81, 'lon',   23.95, 'name', 'UKLL Lviv');
            raob_list.s34009 = struct('lat',   51.76, 'lon',   36.16, 'name', 'Kursk');
            raob_list.s34122 = struct('lat',   51.65, 'lon',   39.25, 'name', 'UUOO Voronez');
            raob_list.s34172 = struct('lat',   51.56, 'lon',   46.03, 'name', 'Saratov');
            raob_list.s34247 = struct('lat',   50.41, 'lon',   41.05, 'name', 'Kalac');
            raob_list.s34300 = struct('lat',   49.96, 'lon',   36.13, 'name', 'UKHH Kharkiv');
            raob_list.s34467 = struct('lat',   48.68, 'lon',   44.35, 'name', 'URWW Volgograd');
            raob_list.s34731 = struct('lat',   47.25, 'lon',   39.81, 'name', 'URRR Rostov-Na-Donu');
            raob_list.s34858 = struct('lat',   45.91, 'lon',   43.35, 'name', 'Divnoe');
            raob_list.s34882 = struct('lat',   46.27, 'lon',   48.03, 'name', 'Astrakhan');
            raob_list.s35121 = struct('lat',   51.68, 'lon',   55.10, 'name', 'Orenburg');
            raob_list.s35229 = struct('lat',   50.28, 'lon',   57.15, 'name', 'UATT Aktjubinsk');
            raob_list.s35394 = struct('lat',   49.80, 'lon',   73.15, 'name', 'Karaganda');
            raob_list.s35671 = struct('lat',   47.80, 'lon',   67.71, 'name', 'Zhezkazgan');
            raob_list.s35700 = struct('lat',   47.11, 'lon',   51.91, 'name', 'Atyran');
            raob_list.s36003 = struct('lat',   52.30, 'lon',   76.93, 'name', 'Pavlodar');
            raob_list.s36096 = struct('lat',   51.71, 'lon',   94.50, 'name', 'Kyzyl');
            raob_list.s36872 = struct('lat',   43.36, 'lon',   77.00, 'name', 'Almaty');
            raob_list.s37011 = struct('lat',   44.10, 'lon',   39.07, 'name', 'Tuapse');
            raob_list.s37055 = struct('lat',   44.22, 'lon',   43.10, 'name', 'URMM Mineralnye Vody');
            raob_list.s37259 = struct('lat',   43.01, 'lon',   47.48, 'name', 'Mahachkala');
            raob_list.s37789 = struct('lat',   40.13, 'lon',   44.46, 'name', 'UGEE Yerevan');
            raob_list.s38064 = struct('lat',   44.77, 'lon',   65.52, 'name', 'Kyzylorda');
            raob_list.s38341 = struct('lat',   42.85, 'lon',   71.38, 'name', 'Zhambyl');
            raob_list.s40179 = struct('lat',   32.00, 'lon',   34.81, 'name', 'Bet Dagan');
            raob_list.s40265 = struct('lat',   32.36, 'lon',   36.25, 'name', 'OJMF Mafraq');
            raob_list.s40373 = struct('lat',   28.31, 'lon',   46.13, 'name', 'OEPA Al-Qaisumah');
            raob_list.s40375 = struct('lat',   28.38, 'lon',   36.60, 'name', 'OETB Tabuk');
            raob_list.s40394 = struct('lat',   27.43, 'lon',   41.68, 'name', 'OEHL Hail');
            raob_list.s40417 = struct('lat',   26.45, 'lon',   49.81, 'name', 'OEDF K.F.I.A.-Dammam');
            raob_list.s40430 = struct('lat',   24.55, 'lon',   39.70, 'name', 'OEMA Al-Madinah');
            raob_list.s40437 = struct('lat',   24.93, 'lon',   46.71, 'name', 'OERK King Khaled Intl Arpt');
            raob_list.s40706 = struct('lat',   38.08, 'lon',   46.28, 'name', 'OITT Tabriz');
            raob_list.s40745 = struct('lat',   36.26, 'lon',   59.63, 'name', 'OIMM Mashhad');
            raob_list.s40754 = struct('lat',   35.68, 'lon',   51.35, 'name', 'OIII Tehran-Mehrabad');
            raob_list.s40766 = struct('lat',   34.26, 'lon',   47.11, 'name', 'OICC Kermanshah');
            raob_list.s40800 = struct('lat',   32.46, 'lon',   51.71, 'name', 'OIFM Esfahan');
            raob_list.s40848 = struct('lat',   29.53, 'lon',   52.58, 'name', 'OISS Shiraz');
            raob_list.s40856 = struct('lat',   29.46, 'lon',   60.88, 'name', 'OIZH Zahedan');
            raob_list.s40948 = struct('lat',   34.55, 'lon',   69.21, 'name', 'OAKB Kabul Airport');
            raob_list.s41024 = struct('lat',   21.70, 'lon',   39.18, 'name', 'OEJN Jeddah (King Abdul Aziz)');
            raob_list.s41112 = struct('lat',   18.23, 'lon',   42.65, 'name', 'OEAB Abha');
            raob_list.s41217 = struct('lat',   24.43, 'lon',   54.65, 'name', 'OMAA Abu Dhabi Inter Arpt');
            raob_list.s41256 = struct('lat',   23.58, 'lon',   58.28, 'name', 'OOMS Seeb, Intl Airport');
            raob_list.s41316 = struct('lat',   17.03, 'lon',   54.08, 'name', 'OOSA Salalah');
            raob_list.s41883 = struct('lat',   24.85, 'lon',   89.36, 'name', 'Bogra');
            raob_list.s41891 = struct('lat',   24.90, 'lon',   91.88, 'name', 'VGSY Sylhet');
            raob_list.s41923 = struct('lat',   23.76, 'lon',   90.38, 'name', 'VGTJ Dhaka');
            raob_list.s42027 = struct('lat',   34.08, 'lon',   74.83, 'name', 'Srinagar');
            raob_list.s42101 = struct('lat',   30.33, 'lon',   76.46, 'name', 'Patiala');
            raob_list.s42182 = struct('lat',   28.58, 'lon',   77.20, 'name', 'VIDD New Delhi');
            raob_list.s42299 = struct('lat',   27.33, 'lon',   88.61, 'name', 'Gangtok');
            raob_list.s42314 = struct('lat',   27.48, 'lon',   95.01, 'name', 'VEMN Dibrugarh');
            raob_list.s42339 = struct('lat',   26.30, 'lon',   73.01, 'name', 'VIJO Jodhpur');
            raob_list.s42361 = struct('lat',   26.23, 'lon',   78.25, 'name', 'VIGR Gwalior');
            raob_list.s42369 = struct('lat',   26.75, 'lon',   80.88, 'name', 'VILK Lucknow');
            raob_list.s42379 = struct('lat',   26.75, 'lon',   83.36, 'name', 'VEGK Gorakhpur');
            raob_list.s42410 = struct('lat',   26.10, 'lon',   91.58, 'name', 'VEGT Gauhati');
            raob_list.s42492 = struct('lat',   25.60, 'lon',   85.10, 'name', 'VEPT Patna');
            raob_list.s42623 = struct('lat',   24.66, 'lon',   93.90, 'name', 'VEIM Imphal');
            raob_list.s42647 = struct('lat',   23.06, 'lon',   72.63, 'name', 'VAAH Ahmadabad');
            raob_list.s42667 = struct('lat',   23.28, 'lon',   77.35, 'name', 'VABP Bhopal');
            raob_list.s42701 = struct('lat',   23.31, 'lon',   85.31, 'name', 'VERC M.O. Ranchi');
            raob_list.s42724 = struct('lat',   23.88, 'lon',   91.25, 'name', 'VEAT Agartala');
            raob_list.s42809 = struct('lat',   22.65, 'lon',   88.45, 'name', 'VECC Calcutta');
            raob_list.s42867 = struct('lat',   21.10, 'lon',   79.05, 'name', 'VANP Nagpur Sonegaon');
            raob_list.s42874 = struct('lat',   21.23, 'lon',   81.65, 'name', 'Pbo Raipur');
            raob_list.s42971 = struct('lat',   20.25, 'lon',   85.83, 'name', 'VEBS Bhubaneswar');
            raob_list.s43003 = struct('lat',   19.11, 'lon',   72.85, 'name', 'VABB Bombay');
            raob_list.s43014 = struct('lat',   19.85, 'lon',   75.40, 'name', 'VAAU Aurangabad Chikalthan');
            raob_list.s43041 = struct('lat',   19.08, 'lon',   82.03, 'name', 'Jagdalpur');
            raob_list.s43128 = struct('lat',   17.45, 'lon',   78.46, 'name', 'VOHY Hyderabad Airport');
            raob_list.s43150 = struct('lat',   17.70, 'lon',   83.30, 'name', 'Vishakhapatnam');
            raob_list.s43192 = struct('lat',   15.48, 'lon',   73.81, 'name', 'Goa');
            raob_list.s43279 = struct('lat',   13.00, 'lon',   80.18, 'name', 'VOMM Madras');
            raob_list.s43285 = struct('lat',   12.95, 'lon',   74.83, 'name', 'Mangalore');
            raob_list.s43295 = struct('lat',   12.96, 'lon',   77.58, 'name', 'Bangalore');
            raob_list.s43311 = struct('lat',   11.12, 'lon',   72.73, 'name', 'Amini Divi');
            raob_list.s43333 = struct('lat',   11.66, 'lon',   92.71, 'name', 'VEPB Port Blair');
            raob_list.s43346 = struct('lat',   10.91, 'lon',   79.83, 'name', 'Karaikal');
            raob_list.s43353 = struct('lat',    9.95, 'lon',   76.26, 'name', 'VOCC Cochin');
            raob_list.s43369 = struct('lat',    8.30, 'lon',   73.15, 'name', 'Minicoy');
            raob_list.s43371 = struct('lat',    8.48, 'lon',   76.95, 'name', 'Thiruvananthapuram');
            raob_list.s43466 = struct('lat',    6.90, 'lon',   79.86, 'name', 'Colombo');
            raob_list.s44231 = struct('lat',   49.56, 'lon',  100.16, 'name', 'Muren');
            raob_list.s44292 = struct('lat',   47.55, 'lon',  106.52, 'name', 'Ulaan-Baator');
            raob_list.s44373 = struct('lat',   43.35, 'lon',  104.25, 'name', 'Dalanzadgad');
            raob_list.s45004 = struct('lat',   22.31, 'lon',  114.17, 'name', 'Kings Park');
            raob_list.s47102 = struct('lat',   37.97, 'lon',  124.63, 'name', 'Baengnyeongdo');
            raob_list.s47104 = struct('lat',   37.81, 'lon',  128.85, 'name', 'Bukgangneung');
            raob_list.s47122 = struct('lat',   37.10, 'lon',  127.03, 'name', 'RKSO Osan Ab');
            raob_list.s47138 = struct('lat',   36.03, 'lon',  129.38, 'name', 'Pohang');
            raob_list.s47158 = struct('lat',   35.11, 'lon',  126.81, 'name', 'RKJJ Kwangju Ab');
            raob_list.s47169 = struct('lat',   34.68, 'lon',  125.45, 'name', 'Heuksando');
            raob_list.s47186 = struct('lat',   33.33, 'lon',  126.68, 'name', 'National Typhoon Centre');
            raob_list.s47401 = struct('lat',   45.41, 'lon',  141.68, 'name', 'Wakkanai');
            raob_list.s47412 = struct('lat',   43.05, 'lon',  141.33, 'name', 'Sapporo');
            raob_list.s47418 = struct('lat',   42.98, 'lon',  144.40, 'name', 'Kushiro');
            raob_list.s47580 = struct('lat',   40.70, 'lon',  141.38, 'name', 'RJSM Misawa Ab');
            raob_list.s47582 = struct('lat',   39.71, 'lon',  140.10, 'name', 'Akita');
            raob_list.s47600 = struct('lat',   37.38, 'lon',  136.90, 'name', 'Wajima');
            raob_list.s47646 = struct('lat',   36.05, 'lon',  140.13, 'name', 'Tateno');
            raob_list.s47678 = struct('lat',   33.11, 'lon',  139.78, 'name', 'Hachijyojima');
            raob_list.s47741 = struct('lat',   35.45, 'lon',  133.07, 'name', 'Matsue');
            raob_list.s47778 = struct('lat',   33.45, 'lon',  135.76, 'name', 'Shionomisaki');
            raob_list.s47807 = struct('lat',   33.58, 'lon',  130.38, 'name', 'Fukuoka');
            raob_list.s47827 = struct('lat',   31.55, 'lon',  130.55, 'name', 'Kagoshima');
            raob_list.s47909 = struct('lat',   28.38, 'lon',  129.55, 'name', 'Naze');
            raob_list.s47918 = struct('lat',   24.33, 'lon',  124.16, 'name', 'ROIG Ishigakijima');
            raob_list.s47945 = struct('lat',   25.83, 'lon',  131.23, 'name', 'ROMD Minamidaitojima');
            raob_list.s47971 = struct('lat',   27.08, 'lon',  142.18, 'name', 'RJAO Chichijima');
            raob_list.s47991 = struct('lat',   24.30, 'lon',  153.96, 'name', 'RJAM Minamitorishima');
            raob_list.s48601 = struct('lat',    5.30, 'lon',  100.26, 'name', 'WMKP Penang');
            raob_list.s48615 = struct('lat',    6.16, 'lon',  102.28, 'name', 'WMKC Kota Bharu');
            raob_list.s48650 = struct('lat',    2.71, 'lon',  101.70, 'name', 'Sepang');
            raob_list.s48657 = struct('lat',    3.78, 'lon',  103.21, 'name', 'WMKD Kuantan');
            raob_list.s48698 = struct('lat',    1.36, 'lon',  103.98, 'name', 'WSSS Singapore');
            raob_list.s48811 = struct('lat',   21.40, 'lon',  103.02, 'name', 'Dien Bien Phu');
            raob_list.s48820 = struct('lat',   21.01, 'lon',  105.80, 'name', 'VVNB Ha Noi');
            raob_list.s48839 = struct('lat',   20.13, 'lon',  107.72, 'name', 'Bach Longvi');
            raob_list.s48845 = struct('lat',   18.68, 'lon',  105.67, 'name', 'VVVH Vinh');
            raob_list.s48855 = struct('lat',   16.03, 'lon',  108.20, 'name', 'VVDN Da Nang');
            raob_list.s48900 = struct('lat',   10.81, 'lon',  106.66, 'name', 'VVTS Ho Chi Minh');
            raob_list.s50527 = struct('lat',   49.21, 'lon',  119.75, 'name', 'Hailar');
            raob_list.s50557 = struct('lat',   49.16, 'lon',  125.23, 'name', 'Nenjiang');
            raob_list.s50774 = struct('lat',   47.71, 'lon',  128.90, 'name', 'Yichun');
            raob_list.s50953 = struct('lat',   45.75, 'lon',  126.76, 'name', 'Harbin');
            raob_list.s51076 = struct('lat',   47.73, 'lon',   88.08, 'name', 'Altay');
            raob_list.s51431 = struct('lat',   43.95, 'lon',   81.33, 'name', 'ZWYN Yining');
            raob_list.s51463 = struct('lat',   43.78, 'lon',   87.62, 'name', 'Urumqi');
            raob_list.s51644 = struct('lat',   41.71, 'lon',   82.95, 'name', 'Kuqa');
            raob_list.s51709 = struct('lat',   39.46, 'lon',   75.98, 'name', 'ZWSH Kashi');
            raob_list.s51777 = struct('lat',   39.03, 'lon',   88.16, 'name', 'Ruoqiang');
            raob_list.s51828 = struct('lat',   37.13, 'lon',   79.93, 'name', 'ZWTN Hotan');
            raob_list.s51839 = struct('lat',   37.06, 'lon',   82.71, 'name', 'Minfeng');
            raob_list.s52203 = struct('lat',   42.81, 'lon',   93.51, 'name', 'ZWHM Hami');
            raob_list.s52267 = struct('lat',   41.95, 'lon',  101.06, 'name', 'Ejin Qi');
            raob_list.s52323 = struct('lat',   41.80, 'lon',   97.03, 'name', 'Mazong Shan');
            raob_list.s52418 = struct('lat',   40.15, 'lon',   94.68, 'name', 'Dunhuang');
            raob_list.s52533 = struct('lat',   39.76, 'lon',   98.48, 'name', 'ZLJQ Jiuquan');
            raob_list.s52681 = struct('lat',   38.63, 'lon',  103.08, 'name', 'Minqin');
            raob_list.s52818 = struct('lat',   36.41, 'lon',   94.90, 'name', 'Golmud');
            raob_list.s52836 = struct('lat',   36.30, 'lon',   98.10, 'name', 'Dulan');
            raob_list.s52866 = struct('lat',   36.71, 'lon',  101.75, 'name', 'ZLXN Xining');
            raob_list.s52983 = struct('lat',   35.87, 'lon',  104.15, 'name', 'Yu Zhong');
            raob_list.s53068 = struct('lat',   43.65, 'lon',  112.00, 'name', 'Erenhot');
            raob_list.s53463 = struct('lat',   40.81, 'lon',  111.68, 'name', 'ZBHH Hohhot');
            raob_list.s53513 = struct('lat',   40.76, 'lon',  107.40, 'name', 'Linhe');
            raob_list.s53614 = struct('lat',   38.48, 'lon',  106.21, 'name', 'ZLIC Yinchuan');
            raob_list.s53772 = struct('lat',   37.78, 'lon',  112.55, 'name', 'ZBYN Taiyuan');
            raob_list.s53845 = struct('lat',   36.60, 'lon',  109.50, 'name', 'ZLYA Yan An');
            raob_list.s53915 = struct('lat',   35.55, 'lon',  106.66, 'name', 'Pingliang');
            raob_list.s54102 = struct('lat',   43.95, 'lon',  116.06, 'name', 'Xilin Hot');
            raob_list.s54135 = struct('lat',   43.60, 'lon',  122.26, 'name', 'Tongliao');
            raob_list.s54161 = struct('lat',   43.90, 'lon',  125.21, 'name', 'ZYCC Changchun');
            raob_list.s54218 = struct('lat',   42.26, 'lon',  118.96, 'name', 'Chifeng');
            raob_list.s54292 = struct('lat',   42.88, 'lon',  129.46, 'name', 'Yanji');
            raob_list.s54374 = struct('lat',   41.71, 'lon',  126.91, 'name', 'Linjiang');
            raob_list.s54511 = struct('lat',   39.93, 'lon',  116.28, 'name', 'ZBAA Beijing');
            raob_list.s54662 = struct('lat',   38.90, 'lon',  121.63, 'name', 'ZYTL Dalian');
            raob_list.s54727 = struct('lat',   36.70, 'lon',  117.55, 'name', 'Zhangqiu');
            raob_list.s54857 = struct('lat',   36.06, 'lon',  120.33, 'name', 'ZSQD Qingdao');
            raob_list.s55299 = struct('lat',   31.48, 'lon',   92.06, 'name', 'Nagqu');
            raob_list.s55591 = struct('lat',   29.66, 'lon',   91.13, 'name', 'ZULS Lhasa');
            raob_list.s56029 = struct('lat',   33.01, 'lon',   97.01, 'name', 'Yushu');
            raob_list.s56080 = struct('lat',   35.00, 'lon',  102.90, 'name', 'Hezuo');
            raob_list.s56137 = struct('lat',   31.15, 'lon',   97.16, 'name', 'Qamdo');
            raob_list.s56146 = struct('lat',   31.61, 'lon',  100.00, 'name', 'Garze');
            raob_list.s56187 = struct('lat',   30.70, 'lon',  103.83, 'name', 'Wenjiang');
            raob_list.s56571 = struct('lat',   27.90, 'lon',  102.26, 'name', 'Xichang');
            raob_list.s56691 = struct('lat',   26.86, 'lon',  104.28, 'name', 'Weining');
            raob_list.s56739 = struct('lat',   25.11, 'lon',   98.48, 'name', 'Tengchong');
            raob_list.s56778 = struct('lat',   25.01, 'lon',  102.68, 'name', 'ZPPP Kunming');
            raob_list.s56964 = struct('lat',   22.76, 'lon',  100.98, 'name', 'Simao');
            raob_list.s56985 = struct('lat',   23.38, 'lon',  103.38, 'name', 'Mengzi');
            raob_list.s57083 = struct('lat',   34.71, 'lon',  113.65, 'name', 'ZHCC Zhengzhou');
            raob_list.s57127 = struct('lat',   33.06, 'lon',  107.03, 'name', 'Hanzhong');
            raob_list.s57131 = struct('lat',   34.43, 'lon',  108.97, 'name', 'Jinghe');
            raob_list.s57178 = struct('lat',   33.03, 'lon',  112.58, 'name', 'Nanyang');
            raob_list.s57447 = struct('lat',   30.28, 'lon',  109.46, 'name', 'Enshi');
            raob_list.s57461 = struct('lat',   30.70, 'lon',  111.30, 'name', 'Yichang');
            raob_list.s57494 = struct('lat',   30.61, 'lon',  114.13, 'name', 'ZHHH Wuhan');
            raob_list.s57516 = struct('lat',   29.51, 'lon',  106.48, 'name', 'ZUCK Chongqing');
            raob_list.s57749 = struct('lat',   27.56, 'lon',  110.00, 'name', 'Huaihua');
            raob_list.s57816 = struct('lat',   26.48, 'lon',  106.65, 'name', 'ZUGY Guiyang');
            raob_list.s57957 = struct('lat',   25.33, 'lon',  110.30, 'name', 'ZGKL Guilin');
            raob_list.s57972 = struct('lat',   25.80, 'lon',  113.03, 'name', 'Chenzhou');
            raob_list.s57993 = struct('lat',   25.85, 'lon',  114.95, 'name', 'ZSGZ Ganzhou');
            raob_list.s58027 = struct('lat',   34.28, 'lon',  117.15, 'name', 'Xuzhou');
            raob_list.s58150 = struct('lat',   33.76, 'lon',  120.25, 'name', 'Sheyang');
            raob_list.s58203 = struct('lat',   32.86, 'lon',  115.73, 'name', 'Fuyang');
            raob_list.s58238 = struct('lat',   32.00, 'lon',  118.80, 'name', 'ZSNJ Nanjing');
            raob_list.s58362 = struct('lat',   31.40, 'lon',  121.46, 'name', 'Shanghai');
            raob_list.s58424 = struct('lat',   30.53, 'lon',  117.05, 'name', 'Anqing');
            raob_list.s58457 = struct('lat',   30.23, 'lon',  120.16, 'name', 'ZSHC Hangzhou');
            raob_list.s58606 = struct('lat',   28.60, 'lon',  115.91, 'name', 'ZSCN Nanchang');
            raob_list.s58633 = struct('lat',   28.96, 'lon',  118.86, 'name', 'Qu Xian');
            raob_list.s58665 = struct('lat',   28.61, 'lon',  121.41, 'name', 'Hongjia');
            raob_list.s58725 = struct('lat',   27.33, 'lon',  117.46, 'name', 'Shaowu');
            raob_list.s58847 = struct('lat',   26.08, 'lon',  119.28, 'name', 'ZSFZ Fuzhou');
            raob_list.s58968 = struct('lat',   25.03, 'lon',  121.51, 'name', 'Taibei');
            raob_list.s59134 = struct('lat',   24.48, 'lon',  118.08, 'name', 'ZSAM Xiamen');
            raob_list.s59211 = struct('lat',   23.90, 'lon',  106.60, 'name', 'Baise');
            raob_list.s59265 = struct('lat',   23.48, 'lon',  111.30, 'name', 'Wuzhou');
            raob_list.s59280 = struct('lat',   23.66, 'lon',  113.05, 'name', 'Qing Yuan');
            raob_list.s59316 = struct('lat',   23.35, 'lon',  116.66, 'name', 'ZGOW Shantou');
            raob_list.s59431 = struct('lat',   22.63, 'lon',  108.21, 'name', 'ZGNN Nanning');
            raob_list.s59758 = struct('lat',   20.03, 'lon',  110.35, 'name', 'ZGHK Haikou');
            raob_list.s59981 = struct('lat',   16.83, 'lon',  112.33, 'name', 'Xisha Dao');
            raob_list.s60018 = struct('lat',   28.47, 'lon',  -16.38, 'name', 'Guimar-Tenerife');
            raob_list.s60155 = struct('lat',   33.56, 'lon',   -7.66, 'name', 'GMMC Casablanca');
            raob_list.s60390 = struct('lat',   36.68, 'lon',    3.21, 'name', 'DAAG Dar-El-Beida');
            raob_list.s60571 = struct('lat',   31.50, 'lon',   -2.25, 'name', 'DAOR Bechar');
            raob_list.s60680 = struct('lat',   22.78, 'lon',    5.52, 'name', 'Tamanrasset');
            raob_list.s60715 = struct('lat',   36.83, 'lon',   10.23, 'name', 'DTTA Tunis-Carthage');
            raob_list.s60760 = struct('lat',   33.91, 'lon',    8.10, 'name', 'DTTZ Tozeur');
            raob_list.s63741 = struct('lat',   -1.30, 'lon',   36.75, 'name', 'HKNC Nairobi');
            raob_list.s64400 = struct('lat',   -4.81, 'lon',   11.90, 'name', 'FCPP Pointe-Noire');
            raob_list.s67083 = struct('lat',  -18.80, 'lon',   47.48, 'name', 'FMMI Antananarivo');
            raob_list.s68263 = struct('lat',  -25.91, 'lon',   28.21, 'name', 'FAIR Pretoria (Irene)');
            raob_list.s68424 = struct('lat',  -28.41, 'lon',   21.26, 'name', 'FAUP Upington');
            raob_list.s68842 = struct('lat',  -33.98, 'lon',   25.61, 'name', 'FAPE Port Elizabeth');
            raob_list.s70026 = struct('lat',   71.28, 'lon', -156.79, 'name', 'PABR Barrow');
            raob_list.s70133 = struct('lat',   66.86, 'lon', -162.63, 'name', 'PAOT Kotzebue, Ralph Wien');
            raob_list.s70200 = struct('lat',   64.50, 'lon', -165.43, 'name', 'PAOM Nome');
            raob_list.s70219 = struct('lat',   60.78, 'lon', -161.84, 'name', 'PABE Bethel');
            raob_list.s70231 = struct('lat',   62.96, 'lon', -155.61, 'name', 'PAMC Mcgrath');
            raob_list.s70261 = struct('lat',   64.81, 'lon', -147.88, 'name', 'PAFA Fairbanks');
            raob_list.s70273 = struct('lat',   61.16, 'lon', -150.01, 'name', 'PANC Anchorage');
            raob_list.s70308 = struct('lat',   57.16, 'lon', -170.22, 'name', 'PASN St. Paul');
            raob_list.s70326 = struct('lat',   58.68, 'lon', -156.67, 'name', 'PAKN King Salmon');
            raob_list.s70350 = struct('lat',   57.75, 'lon', -152.50, 'name', 'PADQ Kodiak');
            raob_list.s70361 = struct('lat',   59.51, 'lon', -139.67, 'name', 'PAYA Yakutat');
            raob_list.s70398 = struct('lat',   55.05, 'lon', -131.59, 'name', 'PANT Annette Island');
            raob_list.s71043 = struct('lat',   65.28, 'lon', -126.75, 'name', 'YVQ Norman Wells Ua');
            raob_list.s71081 = struct('lat',   68.76, 'lon',  -81.21, 'name', 'YUX Hall Beach');
            raob_list.s71082 = struct('lat',   82.50, 'lon',  -62.35, 'name', 'WLT Alert');
            raob_list.s71109 = struct('lat',   50.68, 'lon', -127.36, 'name', 'YZT Port Hardy');
            raob_list.s71119 = struct('lat',   53.53, 'lon', -114.10, 'name', 'WSE Edmonton Stony Plain');
            raob_list.s71126 = struct('lat',   53.55, 'lon', -114.10, 'name', 'ZED Edmonton');
            raob_list.s71603 = struct('lat',   43.87, 'lon',  -66.11, 'name', 'YQI Yarmouth');
            raob_list.s71722 = struct('lat',   46.30, 'lon',  -76.01, 'name', 'WMW Maniwaki');
            raob_list.s71802 = struct('lat',   47.51, 'lon',  -52.78, 'name', 'AYT Mt Pearl');
            raob_list.s71811 = struct('lat',   50.21, 'lon',  -66.24, 'name', 'YZV Sept-Iles');
            raob_list.s71815 = struct('lat',   48.56, 'lon',  -58.56, 'name', 'YJT Stephenville');
            raob_list.s71816 = struct('lat',   53.30, 'lon',  -60.36, 'name', 'YYR Goose Bay');
            raob_list.s71823 = struct('lat',   53.75, 'lon',  -73.66, 'name', 'YAH La Grande Iv');
            raob_list.s71836 = struct('lat',   51.29, 'lon',  -80.62, 'name', 'YMO Moosonee');
            raob_list.s71845 = struct('lat',   51.45, 'lon',  -90.22, 'name', 'WPL Pickle Lake');
            raob_list.s71867 = struct('lat',   53.97, 'lon', -101.09, 'name', 'YQD The Pas');
            raob_list.s71906 = struct('lat',   58.11, 'lon',  -68.41, 'name', 'YVP Kuujjuaq');
            raob_list.s71907 = struct('lat',   58.45, 'lon',  -78.11, 'name', 'WPH Inukjuak');
            raob_list.s71908 = struct('lat',   53.90, 'lon', -122.79, 'name', 'ZXS Prince George');
            raob_list.s71909 = struct('lat',   63.75, 'lon',  -68.55, 'name', 'YFB Iqaluit');
            raob_list.s71913 = struct('lat',   58.73, 'lon',  -94.08, 'name', 'YYQ Churchill');
            raob_list.s71917 = struct('lat',   79.98, 'lon',  -85.93, 'name', 'WEU Eureka');
            raob_list.s71924 = struct('lat',   74.70, 'lon',  -94.97, 'name', 'YRB Resolute');
            raob_list.s71925 = struct('lat',   69.13, 'lon', -105.06, 'name', 'YCB Cambridge Bay');
            raob_list.s71926 = struct('lat',   64.31, 'lon',  -96.00, 'name', 'YBK Baker Lake');
            raob_list.s71934 = struct('lat',   60.03, 'lon', -111.93, 'name', 'YSM Fort Smith');
            raob_list.s71945 = struct('lat',   58.83, 'lon', -122.60, 'name', 'YYE Fort Nelson');
            raob_list.s71957 = struct('lat',   68.31, 'lon', -133.53, 'name', 'YEV Inuvik');
            raob_list.s72201 = struct('lat',   24.55, 'lon',  -81.79, 'name', 'KEY Key West');
            raob_list.s72202 = struct('lat',   25.75, 'lon',  -80.38, 'name', 'MFL Miami');
            raob_list.s72206 = struct('lat',   30.50, 'lon',  -81.70, 'name', 'JAX Jacksonville Intl');
            raob_list.s72208 = struct('lat',   32.90, 'lon',  -80.03, 'name', 'CHS Charleston');
            raob_list.s72210 = struct('lat',   27.70, 'lon',  -82.40, 'name', 'TBW Tampa Bay Area');
            raob_list.s72214 = struct('lat',   30.45, 'lon',  -84.30, 'name', 'TLH Tallahassee Fsu');
            raob_list.s72215 = struct('lat',   33.36, 'lon',  -84.57, 'name', 'FFC Peachtree City');
            raob_list.s72230 = struct('lat',   33.16, 'lon',  -86.76, 'name', 'BMX Shelby County Airport');
            raob_list.s72233 = struct('lat',   30.34, 'lon',  -89.83, 'name', 'LIX Slidell Muni');
            raob_list.s72235 = struct('lat',   32.32, 'lon',  -90.08, 'name', 'JAN Jackson Thompson Fld');
            raob_list.s72240 = struct('lat',   30.11, 'lon',  -93.21, 'name', 'LCH Lake Charles');
            raob_list.s72248 = struct('lat',   32.46, 'lon',  -93.78, 'name', 'SHV Shreveport Reg');
            raob_list.s72249 = struct('lat',   32.83, 'lon',  -97.30, 'name', 'FWD Ft Worth');
            raob_list.s72250 = struct('lat',   25.91, 'lon',  -97.41, 'name', 'BRO Brownsville Intl');
            raob_list.s72251 = struct('lat',   27.76, 'lon',  -97.50, 'name', 'CRP Corpus Christi Intl');
            raob_list.s72261 = struct('lat',   29.36, 'lon', -100.91, 'name', 'DRT Del Rio');
            raob_list.s72265 = struct('lat',   31.95, 'lon', -102.18, 'name', 'MAF Midland');
            raob_list.s72274 = struct('lat',   32.23, 'lon', -110.96, 'name', 'TUS Tucson');
            raob_list.s72293 = struct('lat',   32.85, 'lon', -117.12, 'name', 'NKX San Diego');
            raob_list.s72305 = struct('lat',   34.78, 'lon',  -76.88, 'name', 'MHX Newport');
            raob_list.s72317 = struct('lat',   36.08, 'lon',  -79.95, 'name', 'GSO Greensboro');
            raob_list.s72318 = struct('lat',   37.20, 'lon',  -80.41, 'name', 'RNK Blacksburg');
            raob_list.s72327 = struct('lat',   36.25, 'lon',  -86.57, 'name', 'BNA Nashville');
            raob_list.s72340 = struct('lat',   34.84, 'lon',  -92.26, 'name', 'LZK Little Rock');
            raob_list.s72357 = struct('lat',   35.18, 'lon',  -97.44, 'name', 'OUN Norman');
            raob_list.s72363 = struct('lat',   35.23, 'lon', -101.70, 'name', 'AMA Amarillo Arpt(Awos)');
            raob_list.s72364 = struct('lat',   31.86, 'lon', -106.70, 'name', 'EPZ Santa Teresa');
            raob_list.s72365 = struct('lat',   35.04, 'lon', -106.62, 'name', 'ABQ Albuquerque');
            raob_list.s72376 = struct('lat',   35.23, 'lon', -111.82, 'name', 'FGZ Flagstaff');
            raob_list.s72388 = struct('lat',   36.05, 'lon', -115.18, 'name', 'VEF Las Vegas');
            raob_list.s72393 = struct('lat',   34.75, 'lon', -120.56, 'name', 'VBG Vandenberg Afb');
            raob_list.s72402 = struct('lat',   37.93, 'lon',  -75.47, 'name', 'WAL Wallops Island');
            raob_list.s72403 = struct('lat',   38.98, 'lon',  -77.46, 'name', 'IAD Sterling');
            raob_list.s72426 = struct('lat',   39.41, 'lon',  -83.81, 'name', 'ILN Wilmington');
            raob_list.s72440 = struct('lat',   37.23, 'lon',  -93.38, 'name', 'SGF Springfield');
            raob_list.s72451 = struct('lat',   37.76, 'lon',  -99.97, 'name', 'DDC Dodge City(Awos)');
            raob_list.s72456 = struct('lat',   39.07, 'lon',  -95.62, 'name', 'TOP Topeka');
            raob_list.s72469 = struct('lat',   39.77, 'lon', -104.87, 'name', 'DNR Denver');
            raob_list.s72476 = struct('lat',   39.11, 'lon', -108.53, 'name', 'GJT Grand Junction');
            raob_list.s72489 = struct('lat',   39.56, 'lon', -119.80, 'name', 'REV Reno');
            raob_list.s72493 = struct('lat',   37.73, 'lon', -122.21, 'name', 'OAK Oakland Int');
            raob_list.s72501 = struct('lat',   40.87, 'lon',  -72.86, 'name', 'OKX Upton');
            raob_list.s72518 = struct('lat',   42.69, 'lon',  -73.83, 'name', 'ALB Albany');
            raob_list.s72520 = struct('lat',   40.53, 'lon',  -80.22, 'name', 'PIT Pittsburgh');
            raob_list.s72528 = struct('lat',   42.93, 'lon',  -78.73, 'name', 'BUF Buffalo Int');
            raob_list.s72558 = struct('lat',   41.31, 'lon',  -96.36, 'name', 'OAX Omaha');
            raob_list.s72562 = struct('lat',   41.13, 'lon', -100.68, 'name', 'LBF North Platte');
            raob_list.s72572 = struct('lat',   40.77, 'lon', -111.95, 'name', 'SLC Salt Lake City');
            raob_list.s72582 = struct('lat',   40.86, 'lon', -115.73, 'name', 'LKN Elko');
            raob_list.s72597 = struct('lat',   42.36, 'lon', -122.86, 'name', 'MFR Medford');
            raob_list.s72632 = struct('lat',   42.70, 'lon',  -83.46, 'name', 'DTX White Lake');
            raob_list.s72634 = struct('lat',   44.91, 'lon',  -84.71, 'name', 'APX Gaylord');
            raob_list.s72645 = struct('lat',   44.50, 'lon',  -88.11, 'name', 'GRB Green Bay');
            raob_list.s72649 = struct('lat',   44.85, 'lon',  -93.56, 'name', 'MPX Chanhassen');
            raob_list.s72659 = struct('lat',   45.45, 'lon',  -98.43, 'name', 'ABR Aberdeen');
            raob_list.s72662 = struct('lat',   44.08, 'lon', -103.21, 'name', 'RAP Rapid City');
            raob_list.s72672 = struct('lat',   43.06, 'lon', -108.48, 'name', 'RIW Riverton');
            raob_list.s72681 = struct('lat',   43.56, 'lon', -116.21, 'name', 'BOI Boise');
            raob_list.s72694 = struct('lat',   44.91, 'lon', -123.00, 'name', 'SLE Salem');
            raob_list.s72712 = struct('lat',   46.86, 'lon',  -68.01, 'name', 'CAR Caribou');
            raob_list.s72747 = struct('lat',   48.56, 'lon',  -93.40, 'name', 'INL Int.Falls');
            raob_list.s72764 = struct('lat',   46.76, 'lon', -100.75, 'name', 'BIS Bismarck');
            raob_list.s72768 = struct('lat',   48.21, 'lon', -106.61, 'name', 'GGW Glasgow');
            raob_list.s72776 = struct('lat',   47.46, 'lon', -111.39, 'name', 'TFX Great Falls');
            raob_list.s72786 = struct('lat',   47.68, 'lon', -117.63, 'name', 'OTX Spokane');
            raob_list.s72797 = struct('lat',   47.95, 'lon', -124.55, 'name', 'UIL Quillayute');
            raob_list.s73033 = struct('lat',   50.28, 'lon', -119.27, 'name', 'CWVK Vernon');
            raob_list.s74389 = struct('lat',   43.90, 'lon',  -70.25, 'name', 'GYX Gray');
            raob_list.s74455 = struct('lat',   41.61, 'lon',  -90.58, 'name', 'DVN Davenport');
            raob_list.s74560 = struct('lat',   40.15, 'lon',  -89.33, 'name', 'ILX Lincoln');
            raob_list.s74646 = struct('lat',   36.62, 'lon',  -97.48, 'name', 'LMN Lamont Oklahoma');
            raob_list.s74794 = struct('lat',   28.46, 'lon',  -80.55, 'name', 'XMR Cape Kennedy');
            raob_list.s76225 = struct('lat',   28.67, 'lon', -106.02, 'name', 'Chihuahua, Chih.');
            raob_list.s76458 = struct('lat',   23.20, 'lon', -106.42, 'name', 'Colonia Juancarrasco');
            raob_list.s76526 = struct('lat',   22.75, 'lon', -102.50, 'name', 'Guadalupe');
            raob_list.s76612 = struct('lat',   20.67, 'lon', -103.38, 'name', 'Guadalajara, Jal.');
            raob_list.s76654 = struct('lat',   19.05, 'lon', -104.32, 'name', 'Manzanillo, Col.');
            raob_list.s76679 = struct('lat',   19.40, 'lon',  -99.20, 'name', 'Aerop. Intl Mexico, D.F.');
            raob_list.s78016 = struct('lat',   32.37, 'lon',  -64.68, 'name', 'TXKF Bermuda Nvl Stn Kindley');
            raob_list.s78073 = struct('lat',   25.05, 'lon',  -77.46, 'name', 'MYNN Nassau Airport');
            raob_list.s78384 = struct('lat',   19.30, 'lon',  -81.35, 'name', 'MWCR Owen Roberts Arpt');
            raob_list.s78486 = struct('lat',   18.43, 'lon',  -69.88, 'name', 'MDSD Santo Domingo');
            raob_list.s78526 = struct('lat',   18.43, 'lon',  -65.99, 'name', 'TJSJ San Juan');
            raob_list.s78583 = struct('lat',   17.50, 'lon',  -88.33, 'name', 'MZBZ Phillip Goldston Intl.');
            raob_list.s78954 = struct('lat',   13.07, 'lon',  -59.50, 'name', 'TBPB Grantley Adams');
            raob_list.s78970 = struct('lat',   10.58, 'lon',  -61.35, 'name', 'TTPP Piarco Int. Airport');
            raob_list.s78988 = struct('lat',   12.20, 'lon',  -68.96, 'name', 'TNCC Hato Airport, Curacao');
            raob_list.s80001 = struct('lat',   12.58, 'lon',  -81.72, 'name', 'SKSP San Andres Isl');
            raob_list.s82022 = struct('lat',    2.83, 'lon',  -60.70, 'name', 'SBBV Boa Vista');
            raob_list.s82026 = struct('lat',    2.48, 'lon',  -55.98, 'name', 'Tirios');
            raob_list.s82099 = struct('lat',    0.05, 'lon',  -51.07, 'name', 'SBMQ Macapa');
            raob_list.s82107 = struct('lat',   -1.30, 'lon',  -67.05, 'name', 'Sao Gabriel Da Cachoeira');
            raob_list.s82193 = struct('lat',   -1.38, 'lon',  -48.48, 'name', 'SBBE Belem (Aeroporto)');
            raob_list.s82244 = struct('lat',   -2.43, 'lon',  -54.72, 'name', 'Santarem');
            raob_list.s82281 = struct('lat',   -2.60, 'lon',  -44.23, 'name', 'Sao Luiz');
            raob_list.s82332 = struct('lat',   -3.15, 'lon',  -59.98, 'name', 'SBMN Manaus (Aeroporto)');
            raob_list.s82411 = struct('lat',   -3.67, 'lon',  -69.67, 'name', 'Tabatinga');
            raob_list.s82532 = struct('lat',   -5.82, 'lon',  -61.28, 'name', 'Manicore');
            raob_list.s82599 = struct('lat',   -5.91, 'lon',  -35.25, 'name', 'SBNT Natal Aeroporto');
            raob_list.s82824 = struct('lat',   -8.76, 'lon',  -63.91, 'name', 'SBPV Porto Velho (Aeroporto)');
            raob_list.s82917 = struct('lat',  -10.00, 'lon',  -67.80, 'name', 'Rio Branco');
            raob_list.s82965 = struct('lat',   -9.86, 'lon',  -56.10, 'name', 'SBAT Alta Floresta (Aero)');
            raob_list.s83208 = struct('lat',  -12.70, 'lon',  -60.10, 'name', 'SBVH Vilhena (Aeroporto)');
            raob_list.s83362 = struct('lat',  -15.65, 'lon',  -56.10, 'name', 'SBCY Cuiaba (Aeroporto)');
            raob_list.s83378 = struct('lat',  -15.86, 'lon',  -47.93, 'name', 'SBBR Brasilia (Aeroporto)');
            raob_list.s83525 = struct('lat',  -18.87, 'lon',  -48.22, 'name', 'SBUL Uberlandia');
            raob_list.s83566 = struct('lat',  -19.62, 'lon',  -43.57, 'name', 'Confis Intnl Arpt');
            raob_list.s83612 = struct('lat',  -20.46, 'lon',  -54.66, 'name', 'SBCG Campo Grande (Aero)');
            raob_list.s83649 = struct('lat',  -20.27, 'lon',  -40.28, 'name', 'Vitoria');
            raob_list.s83746 = struct('lat',  -22.81, 'lon',  -43.25, 'name', 'SBGL Galeao');
            raob_list.s83768 = struct('lat',  -23.33, 'lon',  -51.13, 'name', 'SBLO Londrina');
            raob_list.s83779 = struct('lat',  -23.52, 'lon',  -46.63, 'name', 'SBMT Marte Civ');
            raob_list.s83827 = struct('lat',  -25.51, 'lon',  -54.58, 'name', 'SBFI Foz Do Iguacu (Aero)');
            raob_list.s83840 = struct('lat',  -25.51, 'lon',  -49.16, 'name', 'SBCT Curitiba (Aeroporto)');
            raob_list.s83899 = struct('lat',  -27.67, 'lon',  -48.55, 'name', 'SBFL Florianopolis');
            raob_list.s83928 = struct('lat',  -29.78, 'lon',  -57.03, 'name', 'SBUG Uruguaniana');
            raob_list.s83937 = struct('lat',  -29.72, 'lon',  -53.70, 'name', 'SBSM Santa Maria');
            raob_list.s83971 = struct('lat',  -30.00, 'lon',  -51.18, 'name', 'SBPA Porto Alegre (Aero)');
            raob_list.s85586 = struct('lat',  -33.65, 'lon',  -71.61, 'name', 'SCSN Santo Domingo');
            raob_list.s88889 = struct('lat',  -51.81, 'lon',  -58.45, 'name', 'EGYP Mount Pleasant Airport');
            raob_list.s89009 = struct('lat',  -90.00, 'lon',    0.00, 'name', 'Amundsen-Scott');
            raob_list.s89532 = struct('lat',  -69.00, 'lon',   39.58, 'name', 'Syowa');
            raob_list.s89571 = struct('lat',  -68.57, 'lon',   77.97, 'name', 'Davis');
            raob_list.s89592 = struct('lat',  -66.55, 'lon',   93.01, 'name', 'Mirnyj');
            raob_list.s89611 = struct('lat',  -66.28, 'lon',  110.52, 'name', 'Casey');
            raob_list.s89664 = struct('lat',  -77.85, 'lon',  166.66, 'name', 'McMurdo');
            raob_list.s91165 = struct('lat',   21.99, 'lon', -159.34, 'name', 'PHLI Lihue');
            raob_list.s91212 = struct('lat',   13.48, 'lon',  144.79, 'name', 'PGAC Guam Intl Arpt');
            raob_list.s91285 = struct('lat',   19.72, 'lon', -155.05, 'name', 'PHTO Hilo');
            raob_list.s91334 = struct('lat',    7.46, 'lon',  151.84, 'name', 'PTKK Truk');
            raob_list.s91348 = struct('lat',    6.96, 'lon',  158.21, 'name', 'PTPN Ponape');
            raob_list.s91366 = struct('lat',    8.73, 'lon',  167.73, 'name', 'PKWA Kwajalein');
            raob_list.s91376 = struct('lat',    7.07, 'lon',  171.29, 'name', 'PKMJ Majuro');
            raob_list.s91408 = struct('lat',    7.34, 'lon',  134.48, 'name', 'PTRO Koror, Palau Is');
            raob_list.s91413 = struct('lat',    9.50, 'lon',  138.08, 'name', 'PTYA Yap');
            raob_list.s91610 = struct('lat',    1.35, 'lon',  172.91, 'name', 'NGTA Tarawa');
            raob_list.s91643 = struct('lat',   -8.51, 'lon',  179.21, 'name', 'NGFU Funafuti');
            raob_list.s91765 = struct('lat',  -14.33, 'lon', -170.71, 'name', 'NSTU Pago Pago');
            raob_list.s91938 = struct('lat',  -17.55, 'lon', -149.61, 'name', 'NTAA Tahiti-Faaa');
            raob_list.s93112 = struct('lat',  -36.78, 'lon',  174.63, 'name', 'NZWP Whenuapai');
            raob_list.s93417 = struct('lat',  -40.90, 'lon',  174.98, 'name', 'NZPP Paraparaumu Aerodrome');
            raob_list.s93844 = struct('lat',  -46.41, 'lon',  168.31, 'name', 'NZNV Invercargill Aerodrome');
            raob_list.s94120 = struct('lat',  -12.42, 'lon',  130.89, 'name', 'YPDN Darwin Airport');
            raob_list.s94203 = struct('lat',  -17.95, 'lon',  122.23, 'name', 'YBRM Broome Amo');
            raob_list.s94294 = struct('lat',  -19.25, 'lon',  146.76, 'name', 'YBTL Townsville Aero');
            raob_list.s94299 = struct('lat',  -16.30, 'lon',  149.98, 'name', 'Willis Island');
            raob_list.s94312 = struct('lat',  -20.36, 'lon',  118.63, 'name', 'YPPD Port Hedland Amo');
            raob_list.s94326 = struct('lat',  -23.80, 'lon',  133.88, 'name', 'YBAS Alice Springs Aero');
            raob_list.s94332 = struct('lat',  -20.68, 'lon',  139.48, 'name', 'YBMA Mount Isa Amo');
            raob_list.s94374 = struct('lat',  -23.38, 'lon',  150.48, 'name', 'YBRK Rockhampton Aero');
            raob_list.s94403 = struct('lat',  -28.80, 'lon',  114.70, 'name', 'YPGN Geraldton Amo');
            raob_list.s94430 = struct('lat',  -26.61, 'lon',  118.55, 'name', 'YPMR Meekatharra Amo');
            raob_list.s94461 = struct('lat',  -25.03, 'lon',  128.28, 'name', 'Giles');
            raob_list.s94510 = struct('lat',  -26.41, 'lon',  146.26, 'name', 'YBCV Charleville Amo');
            raob_list.s94578 = struct('lat',  -27.38, 'lon',  153.13, 'name', 'YBBN Brisbane Airport Aero');
            raob_list.s94610 = struct('lat',  -31.93, 'lon',  115.96, 'name', 'YPPH Perth Airport');
            raob_list.s94638 = struct('lat',  -33.83, 'lon',  121.88, 'name', 'Esperance Mo');
            raob_list.s94659 = struct('lat',  -31.15, 'lon',  136.81, 'name', 'YPWR Woomera Aerodrome Mo');
            raob_list.s94672 = struct('lat',  -34.95, 'lon',  138.53, 'name', 'YPAD Adelaide Airport');
            raob_list.s94711 = struct('lat',  -31.48, 'lon',  145.83, 'name', 'Cobar Mo');
            raob_list.s94776 = struct('lat',  -32.80, 'lon',  151.83, 'name', 'YSWM Williamtown Amo Raaf');
            raob_list.s94802 = struct('lat',  -34.93, 'lon',  117.80, 'name', 'YPAL Albany Airport');
            raob_list.s94866 = struct('lat',  -37.66, 'lon',  144.85, 'name', 'YMML Melbourne Airport');
            raob_list.s94910 = struct('lat',  -35.16, 'lon',  147.45, 'name', 'YSWG Wagga Wagga Amo');
            raob_list.s94975 = struct('lat',  -42.83, 'lon',  147.50, 'name', 'YMHB Hobart Airport');
            raob_list.s94995 = struct('lat',  -31.53, 'lon',  159.06, 'name', 'Lord Howe Island');
            raob_list.s94996 = struct('lat',  -29.03, 'lon',  167.93, 'name', 'YSNF Norfolk Island Aero');
            raob_list.s94998 = struct('lat',  -54.50, 'lon',  158.95, 'name', 'YMMQ Macquarie Island');
            raob_list.s95527 = struct('lat',  -29.48, 'lon',  149.83, 'name', 'Moree Mo');
            raob_list.s96035 = struct('lat',    3.56, 'lon',   98.68, 'name', 'WIMM Medan');
            raob_list.s96163 = struct('lat',   -0.88, 'lon',  100.35, 'name', 'WIMG Padang');
            raob_list.s96237 = struct('lat',   -2.16, 'lon',  106.13, 'name', 'WIKK Pangkal Pinang');
            raob_list.s96253 = struct('lat',   -3.88, 'lon',  102.33, 'name', 'WIPL Bengkulu');
            raob_list.s96413 = struct('lat',    1.48, 'lon',  110.33, 'name', 'WBGG Kuching');
            raob_list.s96441 = struct('lat',    3.20, 'lon',  113.03, 'name', 'WBGB Bintulu');
            raob_list.s96471 = struct('lat',    5.93, 'lon',  116.05, 'name', 'WBKK Kota Kinabalu');
            raob_list.s96481 = struct('lat',    4.26, 'lon',  117.88, 'name', 'WBKW Tawau');
            raob_list.s96509 = struct('lat',    3.33, 'lon',  117.56, 'name', 'WRLR Tarakan');
            raob_list.s96581 = struct('lat',   -0.15, 'lon',  109.40, 'name', 'WIOO Pontianak');
            raob_list.s96645 = struct('lat',   -2.70, 'lon',  112.70, 'name', 'WRBI Pangkalan Bun');
            raob_list.s96685 = struct('lat',   -3.43, 'lon',  114.75, 'name', 'WRBB Banjarmasin');
            raob_list.s96749 = struct('lat',   -6.11, 'lon',  106.65, 'name', 'WIII Jakarta');
            raob_list.s96805 = struct('lat',   -7.73, 'lon',  109.01, 'name', 'WIIL Cilacap');
            raob_list.s96935 = struct('lat',   -7.36, 'lon',  112.76, 'name', 'WRSJ Surabaya');
            raob_list.s97014 = struct('lat',    1.53, 'lon',  124.91, 'name', 'WAMM Menado');
            raob_list.s97072 = struct('lat',   -0.68, 'lon',  119.73, 'name', 'WAML Palu');
            raob_list.s97180 = struct('lat',   -5.06, 'lon',  119.55, 'name', 'WAAA Ujung Pandang');
            raob_list.s97372 = struct('lat',  -10.16, 'lon',  123.66, 'name', 'WRKK Kupang');
            raob_list.s97502 = struct('lat',   -0.93, 'lon',  131.11, 'name', 'WASS Sorong');
            raob_list.s97560 = struct('lat',   -1.18, 'lon',  136.11, 'name', 'WABB Biak');
            raob_list.s97724 = struct('lat',   -3.70, 'lon',  128.08, 'name', 'WAPP Ambon');
            raob_list.s97900 = struct('lat',   -7.98, 'lon',  131.30, 'name', 'WAPI Saumlaki');
            raob_list.s97980 = struct('lat',   -8.46, 'lon',  140.38, 'name', 'WAKK Merauke');
            raob_list.s98223 = struct('lat',   18.18, 'lon',  120.53, 'name', 'RPLI Laoag');
            raob_list.s98328 = struct('lat',   16.41, 'lon',  120.60, 'name', 'RPUB Baguio');
            raob_list.s98433 = struct('lat',   14.56, 'lon',  121.36, 'name', 'Tanay');
            raob_list.s98444 = struct('lat',   13.13, 'lon',  123.73, 'name', 'RPMP Legaspi');
            raob_list.s98618 = struct('lat',    9.75, 'lon',  118.73, 'name', 'RPVP Puerto Princesa');
            raob_list.s98646 = struct('lat',   10.30, 'lon',  123.96, 'name', 'RPMT Mactan');
            raob_list.s98747 = struct('lat',    8.41, 'lon',  124.61, 'name', 'Cagayan De Oro');
            raob_list.s98753 = struct('lat',    7.11, 'lon',  125.65, 'name', 'RPMD Davao Airport');
        end
        
        function [sta_list, lat_sta, lon_sta] = getCloseStations(lat, lon, n_stations)
            % Get the sta id of the radiosondes in a certain area,
            % or close to a point
            %
            % SYNTAX
            %   [sta_list, lat_sta, lon_sta] = getCloseStations(lat_lim, lon_lim);
            %   [sta_list, lat_sta, lon_sta] = getCloseStations(lat, lon, n_stations)
            %
            % EXAMPLE
            %   [sta_list, lat_sta, lon_sta] = Radiosonde.getCloseStations([43 46], [5 10]);
            %   [sta_list, lat_sta, lon_sta] = Radiosonde.getCloseStations(44, 7, n_stations)
            
            raob_list = Radiosonde.getRaobList;
            sta_list = fieldnames(raob_list);
            raob_list = struct2array(raob_list);
            if nargin < 3 || isempty(n_stations)
                n_stations = numel(raob_list);
            end
            lat_sta = [raob_list.lat]';
            lon_sta = [raob_list.lon]';
            
            if numel(lat) > 1
                % Use a box % it does not work at longitute 180
                % ToDo: manage it
                
                % lat is a limit
                lat = sort(lat);
                lon = sort(lon);
                
                id_ok = (lat_sta >= lat(1)) & (lat_sta <= lat(2)) & ...
                    (lon_sta >= lon(1)) & (lon_sta <= lon(2));
                lat_sta = lat_sta(id_ok);
                lon_sta = lon_sta(id_ok);
                sta_list = sta_list(id_ok);
            else
                % distance from a coordinate
                d = sphericalDistance(lat_sta, lon_sta, lat, lon);
                [d, ids] = sort(d);
                lat_sta = lat_sta(ids);
                lon_sta = lon_sta(ids);
                sta_list = sta_list(ids);
            end
            
            if nargin == 3 && ~isempty(n_stations)
                lat_sta = lat_sta(1 : n_stations);
                lon_sta = lon_sta(1 : n_stations);
                sta_list = sta_list(1 : n_stations);
            end
            
            for s = 1 : numel(sta_list)
                sta_list{s} = sta_list{s}(2 : end);
            end            
        end
        
        % RESERVED
        function rds = getAllStations()
            % RESERVED FUNCTION
            % From the list of station hand inserted here get station position
            % and generate a struct to be used within this Class
            %
            % SYNTAX
            %   Radiosonde.getAllStations()
            
            %% page = 'http://weather.uwyo.edu/upperair/sounding.html';
            raob_code = {'16429'; '16546'; '16754'; '17196'; '17220'; '17240'; '17281'; '17351'; '40179'; '40265'; '40373'; '40375'; '40394'; '40417'; '40430';
                '40437'; '40706'; '40745'; '40754'; '40766'; '40800'; '40848'; '40856'; '40948'; '41024'; '41112'; '41217'; '41256'; '41316'; '41624'; '41715';
                '41718'; '41749'; '60018'; '60155'; '60390'; '60571'; '60680'; '60715'; '60760'; '63741'; '64400'; '64450'; '65578'; '67083'; '68263'; '68424';
                '68842'; '16429'; '16546'; '16754'; '17196'; '17220'; '17240'; '17281'; '17351'; '40179'; '40265'; '40373'; '40375'; '40394'; '40417'; '40430';
                '40437'; '40706'; '40745'; '40754'; '40766'; '40800'; '40848'; '40856'; '40948'; '41024'; '41112'; '41217'; '41256'; '41316'; '41624'; '41715';
                '41718'; '41749'; '60018'; '60155'; '60390'; '60571'; '60680'; '60715'; '60760'; '63741'; '64400'; '64450'; '65578'; '67083'; '68263'; '68424';
                '68842'; '12374'; '13388'; '15420'; '16622'; '16754'; '17030'; '17064'; '17095'; '17130'; '17196'; '17220'; '17240'; '17281'; '17351'; '27962';
                '27995'; '28951'; '33041'; '33345'; '33393'; '34009'; '34122'; '34172'; '34247'; '34300'; '34467'; '34731'; '34858'; '34882'; '35121'; '35229';
                '35700'; '37011'; '37055'; '37259'; '37789'; '40179'; '40265'; '40373'; '40375'; '40394'; '40417'; '40430'; '40437'; '40706'; '40745'; '40754';
                '40766'; '40800'; '40848'; '40856'; '41024'; '41112'; '41217'; '41256'; '41316'; '29839'; '29862'; '30635'; '30715'; '30758'; '30935'; '30965';
                '31510'; '31538'; '31736'; '31873'; '31977'; '32150'; '35394'; '35671'; '36003'; '36096'; '36872'; '38064'; '38341'; '40745'; '40856'; '40948';
                '41624'; '41715'; '41718'; '41749'; '41859'; '41883'; '41891'; '41907'; '41923'; '41936'; '41943'; '41950'; '41977'; '41992'; '42027'; '42101';
                '42182'; '42299'; '42314'; '42339'; '42361'; '42369'; '42379'; '42410'; '42492'; '42623'; '42647'; '42667'; '42701'; '42724'; '42809'; '42867';
                '42874'; '42971'; '43003'; '43014'; '43041'; '43128'; '43150'; '43192'; '43279'; '43285'; '43295'; '43311'; '43333'; '43346'; '43353'; '43369';
                '43371'; '43466'; '43497'; '44231'; '44292'; '44373'; '45004'; '47102'; '47104'; '47122'; '47138'; '47158'; '47169'; '47186'; '47401'; '47412';
                '47418'; '47580'; '47582'; '47600'; '47646'; '47678'; '47741'; '47778'; '47807'; '47827'; '47909'; '47918'; '47945'; '48453'; '48477'; '48500';
                '48601'; '48615'; '48650'; '48657'; '48698'; '48811'; '48820'; '48839'; '48845'; '48855'; '48870'; '48877'; '48887'; '48900'; '48914'; '50527';
                '50557'; '50774'; '50953'; '51076'; '51431'; '51463'; '51644'; '51709'; '51777'; '51828'; '51839'; '52203'; '52267'; '52323'; '52418'; '52533';
                '52681'; '52818'; '52836'; '52866'; '52983'; '53068'; '53463'; '53513'; '53614'; '53772'; '53845'; '53915'; '54102'; '54135'; '54161'; '54218';
                '54292'; '54374'; '54511'; '54662'; '54727'; '54857'; '55299'; '55591'; '56029'; '56080'; '56137'; '56146'; '56187'; '56571'; '56691'; '56739';
                '56778'; '56964'; '56985'; '57083'; '57127'; '57131'; '57178'; '57447'; '57461'; '57494'; '57516'; '57749'; '57816'; '57957'; '57972'; '57993';
                '58027'; '58150'; '58203'; '58238'; '58362'; '58424'; '58457'; '58606'; '58633'; '58665'; '58725'; '58847'; '58968'; '59134'; '59211'; '59265';
                '59280'; '59316'; '59431'; '59758'; '59981'; '96009'; '96015'; '96035'; '96091'; '96109'; '96147'; '96163'; '96171'; '96179'; '96195'; '96221';
                '96237'; '96249'; '96253'; '96295'; '96315'; '96413'; '96441'; '96471'; '96481'; '96509'; '96535'; '96581'; '96633'; '96645'; '96685'; '96737';
                '96739'; '96749'; '96791'; '96797'; '96805'; '96839'; '96935'; '96987'; '97014'; '97028'; '97048'; '97072'; '97086'; '97180'; '97230'; '97240';
                '97260'; '97270'; '97300'; '97340'; '97372'; '97430'; '97600'; '98223'; '98328'; '98433'; '98444'; '98618'; '98646'; '98747'; '98753'; '01028';
                '01415'; '03005'; '03238'; '03354'; '04018'; '04220'; '04270'; '04320'; '04339'; '04360'; '04417'; '06011'; '10035'; '10113'; '10184'; '10393';
                '12120'; '12374'; '20046'; '20292'; '20674'; '20744'; '21432'; '21824'; '21946'; '22008'; '22217'; '22271'; '22522'; '22543'; '22820'; '22845';
                '23078'; '23205'; '23330'; '23415'; '23472'; '23802'; '23884'; '23921'; '23933'; '23955'; '24122'; '24266'; '24343'; '24507'; '24641'; '24688';
                '24726'; '24908'; '24947'; '24959'; '25123'; '25403'; '25428'; '25703'; '25913'; '26038'; '26075'; '26298'; '26702'; '26781'; '27038'; '27199';
                '27459'; '27594'; '27713'; '27707'; '27730'; '27962'; '27995'; '28225'; '28275'; '28445'; '28661'; '28722'; '28951'; '29231'; '29263'; '29282';
                '29572'; '29612'; '29634'; '29698'; '29839'; '29862'; '30054'; '30309'; '30372'; '30635'; '30715'; '30758'; '30965'; '31004'; '31088'; '31168';
                '31369'; '31510'; '31538'; '31736'; '31770'; '31873'; '31977'; '32061'; '32098'; '32150'; '32215'; '32477'; '32540'; '32618'; '33041'; '33345';
                '33393'; '34009'; '34122'; '34172'; '34247'; '34300'; '34467'; '34731'; '34858'; '34882'; '35121'; '35229'; '35700'; '36003'; '37011'; '37055';
                '37259'; '37789'; '47401'; '50527'; '50557'; '50774'; '50953'; '70026'; '70133'; '70200'; '70219'; '70231'; '70261'; '70273'; '70308'; '70326';
                '70350'; '70361'; '70398'; '71043'; '71081'; '71082'; '71109'; '71119'; '71126'; '71802'; '71811'; '71815'; '71816'; '71823'; '71867'; '71906';
                '71907'; '71908'; '71909'; '71913'; '71917'; '71924'; '71925'; '71926'; '71934'; '71945'; '71957'; '72797'; '73033'; '89009'; '89532'; '89571';
                '89592'; '89611'; '89664'; '93112'; '93417'; '93844'; '45004'; '47909'; '47918'; '47945'; '47971'; '47991'; '48839'; '48855'; '48870'; '48877';
                '48887'; '48900'; '57516'; '57749'; '57816'; '57957'; '57972'; '57993'; '58606'; '58633'; '58665'; '58725'; '58847'; '58968'; '59134'; '59211';
                '59265'; '59280'; '59316'; '59431'; '59758'; '59981'; '91165'; '91212'; '91285'; '91334'; '91348'; '91366'; '91376'; '91408'; '91413'; '91610';
                '91643'; '91765'; '94120'; '94203'; '94294'; '94299'; '94312'; '94326'; '94332'; '94374'; '94403'; '94430'; '94461'; '94510'; '94578'; '94610';
                '94638'; '94653'; '94659'; '94672'; '94711'; '94776'; '94802'; '94866'; '94910'; '94975'; '94995'; '94996'; '95527'; '96147'; '96249'; '96315';
                '96413'; '96441'; '96471'; '96481'; '96509'; '96535'; '96581'; '96633'; '96645'; '96685'; '96739'; '96749'; '96791'; '96797'; '96805'; '96839';
                '96935'; '96987'; '97014'; '97028'; '97048'; '97072'; '97086'; '97180'; '97230'; '97240'; '97260'; '97270'; '97300'; '97340'; '97372'; '97430';
                '97502'; '97560'; '97600'; '97724'; '97748'; '97810'; '97900'; '97980'; '98223'; '98328'; '98433'; '98444'; '98618'; '98646'; '98747'; '98753';
                '78954'; '78970'; '78988'; '80001'; '82022'; '82099'; '82026'; '82107'; '82193'; '82244'; '82281'; '82332'; '82400'; '82411'; '82532'; '82599';
                '82824'; '82917'; '82965'; '83208'; '83362'; '83378'; '83525'; '83554'; '83566'; '83612'; '83649'; '83746'; '83768'; '83779'; '83827'; '83840';
                '83899'; '83928'; '83937'; '83971'; '85586'; '88889'; '03953'; '04220'; '04270'; '04360'; '32618'; '70026'; '70133'; '70200'; '70219'; '70231';
                '70261'; '70273'; '70308'; '70326'; '70350'; '70361'; '70398'; '71043'; '71081'; '71109'; '71119'; '71126'; '71603'; '71722'; '71802'; '71811';
                '71815'; '71816'; '71823'; '71836'; '71845'; '71867'; '71906'; '71907'; '71908'; '71909'; '71913'; '71924'; '71925'; '71926'; '71934'; '71945';
                '71957'; '72201'; '72202'; '72206'; '72208'; '72210'; '72214'; '72215'; '72230'; '72233'; '72235'; '72240'; '72248'; '72249'; '72250'; '72251';
                '72261'; '72265'; '72274'; '72293'; '72305'; '72317'; '72318'; '72327'; '72340'; '72357'; '72363'; '72364'; '72365'; '72376'; '72388'; '72393';
                '72402'; '72403'; '72426'; '72440'; '72451'; '72456'; '72469'; '72476'; '72489'; '72493'; '72501'; '72518'; '72520'; '72528'; '72558'; '72562';
                '72572'; '72582'; '72597'; '72632'; '72634'; '72645'; '72649'; '72659'; '72662'; '72672'; '72681'; '72694'; '72712'; '72747'; '72764'; '72768';
                '72776'; '72786'; '72797'; '73033'; '74389'; '74455'; '74560'; '74646'; '74794'; '76225'; '76458'; '76526'; '76612'; '76654'; '76679'; '78073';
                '78384'; '78486'; '78526'; '78583'; '78954'; '78970'; '78988'; '80001'; '82022'; '91165'; '91285'; '01415'; '02365'; '02591'; '03005'; '03808';
                '03953'; '06011'; '06458'; '08579'; '10035'; '10113'; '10184'; '10393'; '10410'; '10548'; '10739'; '10868'; '11520'; '11747'; '11952'; '12120';
                '12374'; '12425'; '12843'; '12982'; '13275'; '13388'; '14240'; '14430'; '15420'; '15614'; '16045'; '16064'; '16113'; '16245'; '16320'; '16429';
                '16546'; '16716'; '17030'; '17064'; '17130'; '17196'; '17220'; '17240'; '17351'; '17516'; '22820'; '22845'; '26075'; '26298'; '26702'; '26781';
                '27038'; '27199'; '27459'; '27594'; '27713'; '27707'; '27730'; '27962'; '27995'; '33317'; '33345'; '34009'; '34122'; '34172'; '34247'; '34467';
                '34731'; '37011'; '01004'; '94998'; '78016'; '91938'};
            
            raob_code = sort(unique(raob_code));
            clear rds;
            for i = 1 : numel(raob_code)
                fprintf(' Checking %03d/%03d:\n', i, numel(raob_code));
                rds(i) = Radiosonde.fromList(raob_code{i}, GPS_Time.now.addIntSeconds(-86400), GPS_Time.now.addIntSeconds(-43200));
                %rds(i) = Radiosonde.fromList(raob_code{i}, GPS_Time('2019-01-01'), GPS_Time('2019-01-02'));
            end
            
            clc
            fprintf('raob_list = struct();\n');
            for i = 1 : numel(raob_code)
                % If the station have a position
                if ~isnan(rds(i).getLat())
                    fprintf('raob_list.s%s = struct(''lat'', %7.2f, ''lon'', %7.2f, ''name'', ''%s'');\n', rds(i).sta_num, rds(i).getLat(), rds(i).getLon(), rds(i).getName());
                    %fprintf('%5s  %+7.2f  %+7.2f  %s\n', rds(i).sta_num, rds(i).getLat(), rds(i).getLon(), rds(i).getName());
                end
            end
        end
    end
end
