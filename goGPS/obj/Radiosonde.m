classdef Radiosonde < handle
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
    %  Written by:       Alice Bonfiglio, Stefano Barindelli
    %  Contributors:     Andrea Gatti, Alessandra Mascitelli
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
        sta_num        % station number
        
        lat             % latitude
        lon             % logitude
        elevation       % elevation
        
        data_time       % time as datetime                  datetime [n_records x 1]
        pressure        % pressure [hPa]                    double   [n_records x 1]
        height          % height [m]                        double   [n_records x 1]
        temperature     % temperature [°C]                  double   [n_records x 1]
        rel_humidity    % relatibe humidity [%]             double   [n_records x 1]
        
        ref_time        % launch epoch;
        pwv             % precipitable water vapor [mm]     double   [n_launches x 1]
        ztd             % zenith tropospheric delay
    end
    
    properties (Constant)
        JAPAN_STATION = {'47401', '47418', '47412', '47580', '47582', '47600', '47646', '47681', '47678', '47741', '47778', '47807', '47827', '47909', '47945', '47918'};
        ITALY_STATION = {'16045', '16080', '16113', '16245', '16320', '16429', '16546'};
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
            date_list = datevec(floor(date_start.getMatlabTime() * 2)/2 : 0.5 : ceil(date_stop.getMatlabTime() * 2)/2);
            [~, id_unique] = unique(date_list(:,1)*1e5 + date_list(:,2));
            year = date_list(id_unique, 1);
            month = date_list(id_unique, 2);            
            
            % sta_num_cell=cellstr(sta_num);
            sta_num_struct = horzcat(repmat('s_', size(sta_num, 1), 1), sta_num);
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
                
                options = weboptions;
                options.Timeout = 10;
                plot_type = 'TEXT';
                % http://weather.uwyo.edu/cgi-bin/sounding?region=seasia&TYPE=TEXT%3ALIST&YEAR=2019&MONTH=05&FROM=0600&TO=0600&STNM=32150
                address = ['http://weather.uwyo.edu/cgi-bin/sounding?region=' region '&TYPE=' plot_type '%3ALIST&YEAR=' sprintf('%04d', year(e)) '&MONTH=' sprintf('%04d', month(e,:)) '&FROM=' from '&TO=' to '&STNM=' sta_num(:,:)];
                fprintf('           get %s\n', address);
                try
                    char_array = webread(address);
                catch
                    % try again after some seconds
                    try
                        fprintf('           hretry access to %s\n', address);
                        pause(3 + randi(1,1) * 5);
                        char_array = webread(address);
                    catch
                    end
                end
                char_array = regexprep(char_array,'<script.*?/script>','');
                char_array = regexprep(char_array,'<style.*?/style>','');
                char_array = regexprep(char_array,'<.*?>','');
                this.lat = str2double(regexp(char_array,'(?<=(Station latitude\: ))([0-9]|\.)*','match', 'once'));
                this.lon = str2double(regexp(char_array,'(?<=(Station longitude\: ))([0-9]|\.)*','match', 'once'));
                this.elevation = str2double(regexp(char_array,'(?<=(Station elevation\: ))([0-9]|\.)*','match', 'once'));
                this.name = regexp(char_array,['(?<=(' sta_num '  )).+?(?=( Observations))'],'match', 'once');
                this.sta_num = sta_num;
                if ~instr(char_array, 'Observations')
                    warning('no data available :-(')
                    rds_pwv(1,j).pwv = [];
                    rds_pwv(1,j).datetime = [];
                    rds_pwv(1,j).sta_numer = sta_num;
                else
                    
                    %pwv
                    str_pw = 'Precipitable water [mm] for entire sounding:';
                    str_idx_pw = strfind(char_array,str_pw);
                    for i=1:length(str_idx_pw)
                        pw_vec(i,:) = char_array(1, str_idx_pw(i)+length(str_pw)+1 : str_idx_pw(i)+length(str_pw)+5);
                        pw_vec_d(i,:) = str2double(pw_vec(i,:));
                    end
                    
                    %time
                    str_obs_time = 'Observation time:';
                    str_idx_time = strfind(char_array,str_obs_time);
                    
                    for i=1:length(str_idx_time)
                        time_vec(i,:)=[char_array(1,str_idx_time(i)+length(str_obs_time)+1:str_idx_time(i)+length(str_obs_time)+7) char_array(1,str_idx_time(i)+length(str_obs_time)+8:str_idx_time(i)+length(str_obs_time)+11) '00'];
                    end
                    datetime_vec = datetime(time_vec,'InputFormat','yyMMdd/HHmmSS');
                    
                    str_idx_header = strfind(char_array,'-----------------------------------------------------------------------------');
                    str_idx_header = str_idx_header(2:2:end);
                    str_idx_header_fin = strfind(char_array,'Station information and sounding indices');
                    for i = 1 : length(str_idx_header)
                        id_1 = length(this.data_time) + 1;
                        temp = char_array(str_idx_header(i)+79:str_idx_header_fin(i)-2);
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
                    
                    rds_pwv(1,j).pwv = pw_vec_d;
                    rds_pwv(1,j).datetime = datetime_vec;
                    rds_pwv(1,j).sta_numer = sta_num(:,:);
                    this.pwv = [this.pwv, rds_pwv(e).pwv'];
                    clear datetime_vec pw_vec_d time_vec                    
                end
                j = j + 1;
            end
        end
    end
    
    methods % Public Access
        function lat = getLat(rds_list)
            lat = [rds_list.lat]';
        end
        function lon = getLon(rds_list)
            lon = [rds_list.lon]';
        end
        function el = getElevation(rds_list)
            el = [rds_list.elevation]';
        end
        
        function [ztd, time] = getZtd(this)
            % Get the zenit total delay computed from sounding
            %
            % SYNTAX
            %   [ztd, time] = getZtd(this)
            ztd = this.ztd * 1e2;
            time = GPS_Time(datenum(this.ref_time)');
        end
        
        function name = getName(this)
            name = this.name;
            if isempty(name)
                name = this.sta_num;
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
                    flag_interval_ref = getOutliers(nan_T);
                    temperature = simpleFill1D(temperature, uint32(flag_interval_ref), 0 );
                end
                
                if sum(nan_RH)
                    flag_interval_ref = getOutliers(nan_RH);
                    relative_humidity = simpleFill1D(relative_humidity, uint32(flag_interval_ref), 0 );
                end
                
                if sum(nan_P)
                    flag_interval_ref = getOutliers(nan_P);
                    pressure = simpleFill1D(pressure, uint32(flag_interval_ref), 0 );
                end
                
                if sum(nan_height)
                    %barometric formula for finding h
                    p0 = 101325; %Pa -> Sea level standard atmospheric pressure
                    %L = 0.0065; %K/m -> Temperature lapse rate, = g/c_p for dry air
                    %c_p = 1007; %J/(kg�K) -> Constant-pressure specific heat
                    t0 = 288.15; %K -> Sea level standard temperature
                    flag_interval_ref = getOutliers(nan_height);
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
    end
end
