%   CLASS Atmosphere
% =========================================================================
%
% DESCRIPTION
%   Class to store static method to compute atmospheric corrections using
%   different models
%
% EXAMPLE
%   ls = LS();
%


%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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
classdef Atmosphere < handle
    properties  (Constant)
        STD_TEMP =  291.15;
        STD_PRES = 1013.25;
        STD_HUMI = 50;
    end
    
    properties  (SetAccess = private, GetAccess = public)
        state
        geoid
        ionex = struct( ...
            'data',       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
            'first_lat',    [], ...    % first latitude
            'first_lon',    [], ...    % first longitude
            'd_lat',      [], ...    %lat spacing
            'd_lon',      [], ...    % lon_spacing
            'n_lat',      [], ...    % num lat
            'n_lon',      [], ...    % num lon
            'first_time', [], ...    % times [time] of the maps
            'first_time_double', [], ...    % times [time] of the maps [seconds from GPS zero]
            'dt',         [], ...    % time spacing
            'n_t',        [], ...    % num of epocvhs
            'heigth',     []  ...    % heigh of the layer
            )
    end
    
    properties  (SetAccess = private, GetAccess = private)
        lat = 1e3;    % current values for lat
        lon = 1e3;    % current values for lon
        undu;         % current values for undu
        
        P       % value of legendre polynomial at the latitude
        V       % V   saved value for GMF computation
        W       % W   saved value for GMF computation
        ahm     % ahm saved value for GMF computation
        aha     % aha saved value for GMF computation
        awm     % awm saved value for GMF computation
        awa     % awa saved value for GMF computation
        
        apm     % awa saved value for GPT computation
        apa     % apa saved value for GPT computation
        atm     % atm saved value for GPT computation
        ata     % ata saved value for GPT computation
        
        emf     % current Earth geomagnetic field object
        log
    end
    
    methods
        function this = Atmosphere()
            % Initialisation of the variables
            gs = Global_Configuration.getInstance;
            this.geoid = gs.getRefGeoid;
            this.state = gs.getCurrentSettings();
            this.log = Logger.getInstance();
        end
        function importIonex(this, filename)
            fid = fopen([filename],'r');
            if fid == -1
                this.log.addWarning(sprintf('      File %s not found', filename));
                return
            end
            this.log.addMessage(sprintf('      Opening file %s for reading', filename));
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
            % read header
            for l = 1:size(lim,1)
                line = txt(lim(l,1):lim(l,2));
                if strfind(line,'END OF HEADER')
                    break
                elseif strfind(line,'EPOCH OF FIRST MAP')
                    first_epoch = GPS_Time(sscanf(line(1:60),'%f %f %f %f %f %f')');
                elseif strfind(line,'EPOCH OF LAST MAP')
                    last_epoch = GPS_Time(sscanf(line(1:60),'%f %f %f %f %f %f')');
                elseif strfind(line,'INTERVAL')
                    interval = sscanf(line(1:60),'%f')';
                elseif strfind(line,'HGT1 / HGT2 / DHGT')
                    height = sscanf(line(1:60),'%f %f %f')';
                elseif strfind(line,'LAT1 / LAT2 / DLAT')
                    lats = sscanf(line(1:60),'%f %f %f')';
                elseif strfind(line,'LON1 / LON2 / DLON')
                    lons = sscanf(line(1:60),'%f %f %f')';
                end
            end
            lim(1:l,:) = [];
            txt(1:(lim(1,1)-1)) = [];
            lim(:,1:2) = lim(:,1:2) - lim(1,1) +1;
            if isempty(this.ionex.data)
                this.ionex.first_time = first_epoch;
                this.ionex.first_time_double = first_epoch.getGpsTime();
                this.ionex.d_t = interval;
                this.ionex.n_t =  round((last_epoch - first_epoch) / interval);
                this.ionex.first_lat= lats(1);
                this.ionex.d_lat = lats(3);
                this.ionex.n_lat = round((lats(2)-lats(1))/lats(3))+1;
                this.ionex.first_lon= lons(1);
                this.ionex.d_lon = lons(3);
                this.ionex.n_lon = round((lons(2)-lons(1))/lons(3))+1;
                this.ionex.height = height;
                this.ionex.data = zeros(round((lats(2)-lats(1))/lats(3))+1, round((lons(2)-lons(1))/lons(3))+1, this.ionex.n_t);
            end
            n_line_1_lat = ceil(size(this.ionex.data,2)*5 / 80);
            n_lat = size(this.ionex.data,1);
            n_lon = size(this.ionex.data,2);
            nt = this.ionex.n_t;
            lines = repmat([false; false; repmat([false; true; false(n_line_1_lat-1,1) ],n_lat,1); false],nt,1);
            st_l  = lim(lines, 1);
            cols = [0:(n_lon*5+n_line_1_lat-2)];
            idx = repmat(cols,length(st_l),1) + repmat(st_l,1,length(cols));
            idx(:,81:81:length(cols))   = [];
            idx(:,366:end) = []; %% trial and error fix bug fix
            vals = txt(serialize(idx'));
            %vals = serialize(vals');
            vals = reshape(vals,5,length(vals)/5);
            nums = sscanf(vals,'%f');
            this.ionex.data = reshape(nums,n_lat,n_lon,nt);
        end
        
        function tec = interpolateTEC(this, gps_time, lat, lon)
            % find indexes and interpolating length
            %time
            dt = this.ionex.d_t;
            nt = this.ionex.n_t;
            ion_gps_time = this.ionex.first_time_double;
            it = max(min(floor((gps_time - ion_gps_time)/ dt)+1,nt-1),1);
            st = max(min(gps_time - ion_gps_time - (it-1)*dt, dt), 0) / dt;
            
            %lat
            dlat = this.ionex.d_lat;
            nlat = this.ionex.n_lat;
            ilat = max(min(floor((lat - this.ionex.first_lat)/ dlat)+1,nlat-1),1);
            slat = min(max(lat - this.ionex.first_lat - (ilat-1)*dlat, dlat), 0) / dlat;
            %lon
            dlon = this.ionex.d_lon;
            nlon = this.ionex.n_lon;
            ilon = max(min(floor((lon - this.ionex.first_lon)/ dlon)+1,nlon-1),1);
            slon = max(min(lon - this.ionex.first_lon - (ilon-1)*dlon, dlon), 0) / dlon;
            
            % interpolate along time
            % [ 1 2  <- index of the cell at the smae time
            %   3 4]
            data = this.ionex.data;
            tec1 = data(ilat   , ilon   , it)*(1-st) + this.ionex.data(ilat   , ilon   , it+1)*st;
            tec2 = data(ilat   , ilon+1 , it)*(1-st) + this.ionex.data(ilat   , ilon+1 , it+1)*st;
            tec3 = data(ilat+1 , ilon   , it)*(1-st) + this.ionex.data(ilat+1 , ilon   , it+1)*st;
            tec4 = data(ilat+1 , ilon+1 , it)*(1-st) + this.ionex.data(ilat+1 , ilon+1 , it+1)*st;
            
            %interpolate along long
            tecn = tec1*(1-slon) + tec2*slon;
            tecs = tec3*(1-slon) + tec4*slon;
            
            %interpolate along lat
            tec = tecn*(1-slat) + tecs*slat;
            
            
        end
        function [stec, pp,mfpp, k] = getSTEC(this,lat,lon, az,el,h, time)
            % get slant total electron component
            
            
            thin_shell_height = this.ionex.height(1)*1000;       %ionopshere thin shell height [km]
            % get piercing point and mapping function
            [latpp, lonpp, mfpp, k] = this.getPiercePoint( lat/180*pi, lon/180*pi, h, az/180*pi, el/180*pi, thin_shell_height,6371000);
            %inetrpolate TEC at piercing point
            tec = this.interpolateTEC( time, latpp*180/pi, lonpp*180/pi);
            
            %apply mapping function
            stec = tec.* mfpp;
            if nargout > 1
                pp = [latpp , lonpp];
            end
            
            
        end
        
        function foi_delay = getFOIdelay(this,lat,lon, az,el,h,time,f)
            %get first horder ionosphere delay
            [stec] = getSTEC(this,lat,lon, az,el,h, time);
            foi_delay = 40.3 * 1e16/ f^2 .* stec;
        end
        
        function [hoi_delay2, hoi_delay3, bending, ppo] = getHOIdelay(this,lat,lon, az,el,h,time,lambda)
            % [1] Fritsche, M., R. Dietrich, C. Knÿfel, A. Rÿlke, S. Vey, M. Rothacher, and P. Steigenberger. Impact
            % of higher-order ionospheric terms on GPS estimates. Geophysical Research Letters, 32(23),
            % 2005. doi: 10.1029/2005GL024342.
            % [2] Odijk, Dennis. "Fast precise GPS positioning in the presence of ionospheric delays." (2002).
            % get High order ionophere delays -- Return phase group
            hoi_delay2 = zeros(size(el));
            hoi_delay3 = zeros(size(el));
            bending = zeros(size(el));
            ppo = zeros([size(el)]);
            gps_time = time.getGpsTime();
            for t = 1: size(el,1)
                idx_ep = find(el(t,:) ~= 0);
                t_time= gps_time(t);
                for s = idx_ep
                    A = 80.6;
                    [stec, pp,mfpp, k] = this.getSTEC(lat,lon, az(t,s),el(t,s),h, t_time);
                    if isempty(this.emf) % do not reaload the model each time
                        this.emf = Earth_Magnetic_Field();
                    end
                    b = this.emf.getB(t_time, GPS_SS.ELL_A/1000 + this.ionex.height(1), pp(2), pp(1));
                    bok = b'*k; %to Tesla
                    c = Global_Configuration.V_LIGHT ;
                    Nemax = (3* 1e12);
                    vtec =  1e18;
                    ni = 0.66;
                    zi = acos(1/mfpp);
                    hoi_delay2(t,s) = - 1 / 2 * 7527 / c^2 * lambda(s)^3 * bok * stec * 1e16;% Eq (10) (11) in [1]
                    hoi_delay3(t,s) = - 1 / 3 * 2437 / c^4 * lambda(s)^4 * Nemax/vtec * ni * (stec  * 1e16)^2;% Eq (1g) (15) (14) in [1]
                    bending(t,s) = A^2 / (8 * c^4) *lambda(s)^4 * tan(zi)^2 * ni * Nemax * stec * 1e16;% Eq(4.34) in [2]
                    ppo(t,s) = stec;
                end
            end
        end
    end
    
    methods
        %-----------------------------------------------------------
        % TROPO
        %-----------------------------------------------------------
        function [delay] = saastamoinenModel (this, h, undu, el)
            % SYNTAX:
            %   [delay] = Atmosphere.tropo_error_correction(lat, lon, h, el);
            %
            % INPUT:
            %   time_rx = receiver reception time
            %   lat = receiver latitude          [degrees]
            %   lon = receiver longitude         [degrees]
            %   h  = receiver ellipsoidal height [meters]  [nx1]
            %   el = satellite elevation         [degrees] [nx1]
            %
            % OUTPUT:
            %   corr = tropospheric error correction
            %
            % DESCRIPTION:
            %   Computation of the pseudorange correction due to tropospheric refraction.
            %   Saastamoinen algorithm using standard atmosphere accounting for
            %   height gradient for temperature pressure and humidity.
            %   --> multi epoch for static receiver
            
            % Saastamoinen model requires (positive) orthometric height
            % ÿÿ undulation is never less than 300 m (on Earth)
            %h(undu > -300) = h(undu > -300) - undu(undu > -300);
            h = h - undu;
            h(h < 0) = 0;
            
            if (h < 5000)
                
                % conversion to radians
                el = abs(el) * pi/180;
                
                % Standard atmosphere - Berg, 1948 (Bernese)
                % pressure [mbar]
                Pr = this.STD_PRES;
                % temperature [K]
                Tr = this.STD_TEMP;
                % humidity [%]
                Hr = this.STD_HUMI;
                
                P = Pr * (1-0.0000226*h).^5.225;
                T = Tr - 0.0065*h;
                H = Hr * exp(-0.0006396*h);
                
                %----------------------------------------------------------------------
                
                %linear interpolation
                h_a = [0; 500; 1000; 1500; 2000; 2500; 3000; 4000; 5000];
                B_a = [1.156; 1.079; 1.006; 0.938; 0.874; 0.813; 0.757; 0.654; 0.563];
                
                t = zeros(length(T),1);
                B = zeros(length(T),1);
                
                for i = 1 : length(T)
                    
                    d = h_a - h(i);
                    [~, j] = min(abs(d));
                    if (d(j) > 0)
                        index = [j-1; j];
                    else
                        index = [j; j+1];
                    end
                    
                    t(i) = (h(i) - h_a(index(1))) ./ (h_a(index(2)) - h_a(index(1)));
                    B(i) = (1-t(i))*B_a(index(1)) + t(i)*B_a(index(2));
                end
                
                %----------------------------------------------------------------------
                
                e = 0.01 * H .* exp(-37.2465 + 0.213166*T - 0.000256908*T.^2);
                
                % tropospheric delay
                w_fun = (1-tan(el).^2./(tan(el).^2+tan(el+2).^2));
                %delay = ((0.002277 ./ sin(el)) .* (P - (B ./ tan(el).^2)) + (0.002277 ./ sin(el)) .* (1255./T + 0.05) .* e);
                delay = ((0.002277 ./ sin(el)) .* (P - (B ./ (tan(el).^2+ 0.01.*w_fun))) + (0.002277 ./ sin(el)) .* (1255./T + 0.05) .* e); % max to eliminate numeric instability near 0
            else
                delay = zeros(size(el));
            end
        end
        
        function [delay] = saastamoinenModelGPT (this, gps_time, lat, lon, h, undu, el)
            % SYNTAX:
            %   [delay] = Atmosphere.saastamoinen_model_GPT(time_rx, lat, lon, h, undu, el)
            %
            % INPUT:
            %   time_rx = receiver reception time
            %   lat = receiver latitude          [degrees]
            %   lon = receiver longitude         [degrees]
            %   h  = receiver ellipsoidal height [meters]
            %   el = satellite elevation         [degrees] [nx1]
            %
            % OUTPUT:
            %   corr = tropospheric error correction
            %
            % DESCRIPTION:
            %   Computation of the pseudorange correction due to tropospheric refraction.
            %   Saastamoinen algorithm using P T from Global Pressure and Temperature
            %   (GPT), and H from standard atmosphere accounting for humidity height gradient.
            %   --> single epoch
            [pres, temp] = this.gpt(gps_time, lat*pi/180, lon*pi/180, h, undu);
            
            t_h = h;
            t_h(undu > -300) = t_h(undu > -300) - undu(undu > -300);
            ZHD_R = saast_dry(pres, t_h, lat);
            ZWD_R = saast_wet(temp, this.STD_HUMI, t_h);
            
            [gmfh_R, gmfw_R] = this.gmf(gps_time, lat*pi/180, lon*pi/180, h - undu, (90-el)*pi/180);
            delay = gmfh_R .* ZHD_R + gmfw_R .* ZWD_R;
        end
        
        function [pres,temp,undu] = gpt(this, gps_time, dlat, dlon, dhgt, undu)
            % This subroutine determines Global Pressure and Temperature
            % based on Spherical Harmonics up to degree and order 9
            %
            % input data
            % ----------
            % dmjd: modified julian date
            % dlat: latitude in radians
            % dlon: longitude in radians
            % dhgt: ellipsoidal height in m
            %
            % output data
            % -----------
            % pres: pressure in hPa
            % temp: temperature in Celsius
            % undu: Geoid undulation in m (from a 9x9 EGM based model)
            
            
            %
            % Johannes Boehm, 2006 June 12
            % rev 2006 June 16: geoid undulation is accounted for
            %
            % Reference:
            % J. Boehm, R. Heinkelmann, and H. Schuh,
            % Global Pressure and Temperature (GPT): A spherical harmonics expansion
            % of annual pressure and temperature variations for geodetic applications,
            % to be submitted to Journal of Geodesy, 2006
            
            % reference day is 28 January
            % this is taken from Niell (1996) to be consistent
            doy = (gps_time/86400 - 22) / 365.25d0;
            
            cached = (dlon == this.lon) && (dlat == this.lat) && ~isempty(this.undu);
            if cached % no cache (debugging purpouse)
                apm = this.apm;
                apa = this.apa;
                atm = this.atm;
                ata = this.ata;
                P   = this.P;
                undu= this.undu;
            else
                a_geoid = [ ...
                    -5.6195d-001,-6.0794d-002,-2.0125d-001,-6.4180d-002,-3.6997d-002, ...
                    +1.0098d+001,+1.6436d+001,+1.4065d+001,+1.9881d+000,+6.4414d-001, ...
                    -4.7482d+000,-3.2290d+000,+5.0652d-001,+3.8279d-001,-2.6646d-002, ...
                    +1.7224d+000,-2.7970d-001,+6.8177d-001,-9.6658d-002,-1.5113d-002, ...
                    +2.9206d-003,-3.4621d+000,-3.8198d-001,+3.2306d-002,+6.9915d-003, ...
                    -2.3068d-003,-1.3548d-003,+4.7324d-006,+2.3527d+000,+1.2985d+000, ...
                    +2.1232d-001,+2.2571d-002,-3.7855d-003,+2.9449d-005,-1.6265d-004, ...
                    +1.1711d-007,+1.6732d+000,+1.9858d-001,+2.3975d-002,-9.0013d-004, ...
                    -2.2475d-003,-3.3095d-005,-1.2040d-005,+2.2010d-006,-1.0083d-006, ...
                    +8.6297d-001,+5.8231d-001,+2.0545d-002,-7.8110d-003,-1.4085d-004, ...
                    -8.8459d-006,+5.7256d-006,-1.5068d-006,+4.0095d-007,-2.4185d-008];
                
                b_geoid = [ ...
                    +0.0000d+000,+0.0000d+000,-6.5993d-002,+0.0000d+000,+6.5364d-002, ...
                    -5.8320d+000,+0.0000d+000,+1.6961d+000,-1.3557d+000,+1.2694d+000, ...
                    +0.0000d+000,-2.9310d+000,+9.4805d-001,-7.6243d-002,+4.1076d-002, ...
                    +0.0000d+000,-5.1808d-001,-3.4583d-001,-4.3632d-002,+2.2101d-003, ...
                    -1.0663d-002,+0.0000d+000,+1.0927d-001,-2.9463d-001,+1.4371d-003, ...
                    -1.1452d-002,-2.8156d-003,-3.5330d-004,+0.0000d+000,+4.4049d-001, ...
                    +5.5653d-002,-2.0396d-002,-1.7312d-003,+3.5805d-005,+7.2682d-005, ...
                    +2.2535d-006,+0.0000d+000,+1.9502d-002,+2.7919d-002,-8.1812d-003, ...
                    +4.4540d-004,+8.8663d-005,+5.5596d-005,+2.4826d-006,+1.0279d-006, ...
                    +0.0000d+000,+6.0529d-002,-3.5824d-002,-5.1367d-003,+3.0119d-005, ...
                    -2.9911d-005,+1.9844d-005,-1.2349d-006,-7.6756d-009,+5.0100d-008];
                
                ap_mean = [ ...
                    +1.0108d+003,+8.4886d+000,+1.4799d+000,-1.3897d+001,+3.7516d-003, ...
                    -1.4936d-001,+1.2232d+001,-7.6615d-001,-6.7699d-002,+8.1002d-003, ...
                    -1.5874d+001,+3.6614d-001,-6.7807d-002,-3.6309d-003,+5.9966d-004, ...
                    +4.8163d+000,-3.7363d-001,-7.2071d-002,+1.9998d-003,-6.2385d-004, ...
                    -3.7916d-004,+4.7609d+000,-3.9534d-001,+8.6667d-003,+1.1569d-002, ...
                    +1.1441d-003,-1.4193d-004,-8.5723d-005,+6.5008d-001,-5.0889d-001, ...
                    -1.5754d-002,-2.8305d-003,+5.7458d-004,+3.2577d-005,-9.6052d-006, ...
                    -2.7974d-006,+1.3530d+000,-2.7271d-001,-3.0276d-004,+3.6286d-003, ...
                    -2.0398d-004,+1.5846d-005,-7.7787d-006,+1.1210d-006,+9.9020d-008, ...
                    +5.5046d-001,-2.7312d-001,+3.2532d-003,-2.4277d-003,+1.1596d-004, ...
                    +2.6421d-007,-1.3263d-006,+2.7322d-007,+1.4058d-007,+4.9414d-009];
                
                bp_mean = [ ...
                    +0.0000d+000,+0.0000d+000,-1.2878d+000,+0.0000d+000,+7.0444d-001, ...
                    +3.3222d-001,+0.0000d+000,-2.9636d-001,+7.2248d-003,+7.9655d-003, ...
                    +0.0000d+000,+1.0854d+000,+1.1145d-002,-3.6513d-002,+3.1527d-003, ...
                    +0.0000d+000,-4.8434d-001,+5.2023d-002,-1.3091d-002,+1.8515d-003, ...
                    +1.5422d-004,+0.0000d+000,+6.8298d-001,+2.5261d-003,-9.9703d-004, ...
                    -1.0829d-003,+1.7688d-004,-3.1418d-005,+0.0000d+000,-3.7018d-001, ...
                    +4.3234d-002,+7.2559d-003,+3.1516d-004,+2.0024d-005,-8.0581d-006, ...
                    -2.3653d-006,+0.0000d+000,+1.0298d-001,-1.5086d-002,+5.6186d-003, ...
                    +3.2613d-005,+4.0567d-005,-1.3925d-006,-3.6219d-007,-2.0176d-008, ...
                    +0.0000d+000,-1.8364d-001,+1.8508d-002,+7.5016d-004,-9.6139d-005, ...
                    -3.1995d-006,+1.3868d-007,-1.9486d-007,+3.0165d-010,-6.4376d-010];
                
                ap_amp = [ ...
                    -1.0444d-001,+1.6618d-001,-6.3974d-002,+1.0922d+000,+5.7472d-001, ...
                    -3.0277d-001,-3.5087d+000,+7.1264d-003,-1.4030d-001,+3.7050d-002, ...
                    +4.0208d-001,-3.0431d-001,-1.3292d-001,+4.6746d-003,-1.5902d-004, ...
                    +2.8624d+000,-3.9315d-001,-6.4371d-002,+1.6444d-002,-2.3403d-003, ...
                    +4.2127d-005,+1.9945d+000,-6.0907d-001,-3.5386d-002,-1.0910d-003, ...
                    -1.2799d-004,+4.0970d-005,+2.2131d-005,-5.3292d-001,-2.9765d-001, ...
                    -3.2877d-002,+1.7691d-003,+5.9692d-005,+3.1725d-005,+2.0741d-005, ...
                    -3.7622d-007,+2.6372d+000,-3.1165d-001,+1.6439d-002,+2.1633d-004, ...
                    +1.7485d-004,+2.1587d-005,+6.1064d-006,-1.3755d-008,-7.8748d-008, ...
                    -5.9152d-001,-1.7676d-001,+8.1807d-003,+1.0445d-003,+2.3432d-004, ...
                    +9.3421d-006,+2.8104d-006,-1.5788d-007,-3.0648d-008,+2.6421d-010];
                
                bp_amp = [ ...
                    +0.0000d+000,+0.0000d+000,+9.3340d-001,+0.0000d+000,+8.2346d-001, ...
                    +2.2082d-001,+0.0000d+000,+9.6177d-001,-1.5650d-002,+1.2708d-003, ...
                    +0.0000d+000,-3.9913d-001,+2.8020d-002,+2.8334d-002,+8.5980d-004, ...
                    +0.0000d+000,+3.0545d-001,-2.1691d-002,+6.4067d-004,-3.6528d-005, ...
                    -1.1166d-004,+0.0000d+000,-7.6974d-002,-1.8986d-002,+5.6896d-003, ...
                    -2.4159d-004,-2.3033d-004,-9.6783d-006,+0.0000d+000,-1.0218d-001, ...
                    -1.3916d-002,-4.1025d-003,-5.1340d-005,-7.0114d-005,-3.3152d-007, ...
                    +1.6901d-006,+0.0000d+000,-1.2422d-002,+2.5072d-003,+1.1205d-003, ...
                    -1.3034d-004,-2.3971d-005,-2.6622d-006,+5.7852d-007,+4.5847d-008, ...
                    +0.0000d+000,+4.4777d-002,-3.0421d-003,+2.6062d-005,-7.2421d-005, ...
                    +1.9119d-006,+3.9236d-007,+2.2390d-007,+2.9765d-009,-4.6452d-009];
                
                at_mean = [ ...
                    +1.6257e+001,+2.1224e+000,+9.2569e-001,-2.5974e+001,+1.4510e+000, ...
                    +9.2468e-002,-5.3192e-001,+2.1094e-001,-6.9210e-002,-3.4060e-002, ...
                    -4.6569e+000,+2.6385e-001,-3.6093e-002,+1.0198e-002,-1.8783e-003, ...
                    +7.4983e-001,+1.1741e-001,+3.9940e-002,+5.1348e-003,+5.9111e-003, ...
                    +8.6133e-006,+6.3057e-001,+1.5203e-001,+3.9702e-002,+4.6334e-003, ...
                    +2.4406e-004,+1.5189e-004,+1.9581e-007,+5.4414e-001,+3.5722e-001, ...
                    +5.2763e-002,+4.1147e-003,-2.7239e-004,-5.9957e-005,+1.6394e-006, ...
                    -7.3045e-007,-2.9394e+000,+5.5579e-002,+1.8852e-002,+3.4272e-003, ...
                    -2.3193e-005,-2.9349e-005,+3.6397e-007,+2.0490e-006,-6.4719e-008, ...
                    -5.2225e-001,+2.0799e-001,+1.3477e-003,+3.1613e-004,-2.2285e-004, ...
                    -1.8137e-005,-1.5177e-007,+6.1343e-007,+7.8566e-008,+1.0749e-009];
                
                bt_mean = [ ...
                    +0.0000e+000,+0.0000e+000,+1.0210e+000,+0.0000e+000,+6.0194e-001, ...
                    +1.2292e-001,+0.0000e+000,-4.2184e-001,+1.8230e-001,+4.2329e-002, ...
                    +0.0000e+000,+9.3312e-002,+9.5346e-002,-1.9724e-003,+5.8776e-003, ...
                    +0.0000e+000,-2.0940e-001,+3.4199e-002,-5.7672e-003,-2.1590e-003, ...
                    +5.6815e-004,+0.0000e+000,+2.2858e-001,+1.2283e-002,-9.3679e-003, ...
                    -1.4233e-003,-1.5962e-004,+4.0160e-005,+0.0000e+000,+3.6353e-002, ...
                    -9.4263e-004,-3.6762e-003,+5.8608e-005,-2.6391e-005,+3.2095e-006, ...
                    -1.1605e-006,+0.0000e+000,+1.6306e-001,+1.3293e-002,-1.1395e-003, ...
                    +5.1097e-005,+3.3977e-005,+7.6449e-006,-1.7602e-007,-7.6558e-008, ...
                    +0.0000e+000,-4.5415e-002,-1.8027e-002,+3.6561e-004,-1.1274e-004, ...
                    +1.3047e-005,+2.0001e-006,-1.5152e-007,-2.7807e-008,+7.7491e-009];
                
                at_amp = [ ...
                    -1.8654e+000,-9.0041e+000,-1.2974e-001,-3.6053e+000,+2.0284e-002, ...
                    +2.1872e-001,-1.3015e+000,+4.0355e-001,+2.2216e-001,-4.0605e-003, ...
                    +1.9623e+000,+4.2887e-001,+2.1437e-001,-1.0061e-002,-1.1368e-003, ...
                    -6.9235e-002,+5.6758e-001,+1.1917e-001,-7.0765e-003,+3.0017e-004, ...
                    +3.0601e-004,+1.6559e+000,+2.0722e-001,+6.0013e-002,+1.7023e-004, ...
                    -9.2424e-004,+1.1269e-005,-6.9911e-006,-2.0886e+000,-6.7879e-002, ...
                    -8.5922e-004,-1.6087e-003,-4.5549e-005,+3.3178e-005,-6.1715e-006, ...
                    -1.4446e-006,-3.7210e-001,+1.5775e-001,-1.7827e-003,-4.4396e-004, ...
                    +2.2844e-004,-1.1215e-005,-2.1120e-006,-9.6421e-007,-1.4170e-008, ...
                    +7.8720e-001,-4.4238e-002,-1.5120e-003,-9.4119e-004,+4.0645e-006, ...
                    -4.9253e-006,-1.8656e-006,-4.0736e-007,-4.9594e-008,+1.6134e-009];
                
                bt_amp = [ ...
                    +0.0000e+000,+0.0000e+000,-8.9895e-001,+0.0000e+000,-1.0790e+000, ...
                    -1.2699e-001,+0.0000e+000,-5.9033e-001,+3.4865e-002,-3.2614e-002, ...
                    +0.0000e+000,-2.4310e-002,+1.5607e-002,-2.9833e-002,-5.9048e-003, ...
                    +0.0000e+000,+2.8383e-001,+4.0509e-002,-1.8834e-002,-1.2654e-003, ...
                    -1.3794e-004,+0.0000e+000,+1.3306e-001,+3.4960e-002,-3.6799e-003, ...
                    -3.5626e-004,+1.4814e-004,+3.7932e-006,+0.0000e+000,+2.0801e-001, ...
                    +6.5640e-003,-3.4893e-003,-2.7395e-004,+7.4296e-005,-7.9927e-006, ...
                    -1.0277e-006,+0.0000e+000,+3.6515e-002,-7.4319e-003,-6.2873e-004, ...
                    -8.2461e-005,+3.1095e-005,-5.3860e-007,-1.2055e-007,-1.1517e-007, ...
                    +0.0000e+000,+3.1404e-002,+1.5580e-002,-1.1428e-003,+3.3529e-005, ...
                    +1.0387e-005,-1.9378e-006,-2.7327e-007,+7.5833e-009,-9.2323e-009];
                
                % Computing Legendre Polynomial
                %parameter t
                t = sin(dlat);
                
                % degree n and order m
                n = 9;
                m = 9;
                
                
                % determine n!  (faktorielle)  moved by 1
                dfac(1) = 1;
                for i = 1:(2*n + 1)
                    dfac(i+1) = dfac(i)*i;
                end
                
                % determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
                for i = 0:n
                    for j = 0:min(i,m)
                        ir = floor((i - j)/2);
                        sum_t = 0;
                        for k = 0:ir
                            sum_t = sum_t + (-1)^k*dfac(2*i - 2*k + 1)/dfac(k + 1)/dfac(i - k + 1)/dfac(i - j - 2*k + 1)*t^(i - j - 2*k);
                        end
                        % Legendre functions moved by 1
                        P(i + 1,j + 1) = 1.d0/2^i*sqrt((1 - t^2)^(j))*sum_t;
                    end
                end
                
                % spherical harmonics
                i = 0;
                for n = 0:9
                    for m = 0:n
                        i = i + 1;
                        aP(i) = P(n+1,m+1)*cos(m*dlon);
                        bP(i) = P(n+1,m+1)*sin(m*dlon);
                    end
                end
                % vectorial computation
                %                 md = (0 : 9)' * dlon;
                %                 aP = bsxfun(@times, P, cos(md)); aP = aP(triu(true(10)))';
                %                 bP = bsxfun(@times, P, sin(md)); bP = bP(triu(true(10)))';
                
                % Geoidal height
                % undu = 0.d0;
                % for i = 1:55
                %    undu = undu + (a_geoid(i)*aP(i) + b_geoid(i)*bP(i));
                % end
                % vectorial computation
                
                if nargin < 6 || isempty(undu)
                    undu = sum(a_geoid .* aP + b_geoid .* bP);
                end
                % now this function get directly the orthometric height from input
                
                % Surface pressure on the geoid
                % apm = 0.d0;
                % apa = 0.d0;
                % for i = 1:55
                %     apm = apm + (ap_mean(i)*aP(i) + bp_mean(i)*bP(i));
                %     apa = apa + (ap_amp(i) *aP(i) + bp_amp(i) *bP(i));
                % end
                % vectorial computation
                apm = sum(ap_mean .* aP + bp_mean .* bP);
                apa = sum(ap_amp  .* aP + bp_amp  .* bP);
                
                % Surface temperature on the geoid
                % atm = 0.d0;
                % ata = 0.d0;
                % for i = 1:55
                %     atm = atm + (at_mean(i)*aP(i) + bt_mean(i)*bP(i));
                %     ata = ata + (at_amp(i) *aP(i) + bt_amp(i) *bP(i));
                % end
                % vectorial computation
                atm = sum(at_mean .* aP + bt_mean .* bP);
                ata = sum(at_amp  .* aP + bt_amp  .* bP);
                
                this.apm = apm;
                this.apa = apa;
                this.atm = atm;
                this.ata = ata;
                this.P = P;
                this.lon = dlon;
                this.lat = dlat;
                this.undu = undu;
            end
            
            % orthometric height
            h_ort = dhgt - undu;
            
            % height correction for pressure
            pres0  = apm + apa*cos(doy*2.d0*pi);
            pres = pres0*(1.d0-0.0000226d0*h_ort)^5.225d0;
            
            % height correction for temperature
            temp0 =  atm + ata*cos(doy*2*pi);
            temp = temp0 - 0.0065d0 * h_ort;
        end
        
        function [gmfh, gmfw] = gmf (this, gps_time, dlat, dlon, dhgt, zd)
            % This subroutine determines the Global Mapping Functions GMF
            % Reference: Boehm, J., A.E. Niell, P. Tregoning, H. Schuh (2006),
            % Global Mapping Functions (GMF): A new empirical mapping function based on numerical weather model data,
            % Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.
            %
            % input data
            % ----------
            % time: GPS_Time
            % dlat: ellipsoidal latitude in radians
            % dlon: longitude in radians
            % dhgt: orthometric height in m
            % zd:   zenith distance in radians
            %
            % output data
            % -----------
            % gmfh: hydrostatic mapping function
            % gmfw: wet mapping function
            %
            % Johannes Boehm, 2005 August 30
            %
            % ref 2006 Aug. 14: recursions for Legendre polynomials (O. Montenbruck)
            % ref 2011 Jul. 21: latitude -> ellipsoidal latitude (J. Boehm)
            
            % reference day is 28 January
            % this is taken from Niell (1996) to be consistent
            doy = (gps_time/86400 - 22) / 365.25d0; % years from 28 jan 1980
            
            pi = 3.14159265359d0;
            
            cached = (dlon == this.lon) && (dlat == this.lat) && ~isempty(this.ahm) && ~isempty(this.aha);
            if cached
                V = this.V;
                W = this.W;
            else
                
                ah_mean = ...
                    [+1.2517d+02, +8.503d-01, +6.936d-02, -6.760d+00, +1.771d-01, ...
                    +1.130d-02, +5.963d-01, +1.808d-02, +2.801d-03, -1.414d-03, ...
                    -1.212d+00, +9.300d-02, +3.683d-03, +1.095d-03, +4.671d-05, ...
                    +3.959d-01, -3.867d-02, +5.413d-03, -5.289d-04, +3.229d-04, ...
                    +2.067d-05, +3.000d-01, +2.031d-02, +5.900d-03, +4.573d-04, ...
                    -7.619d-05, +2.327d-06, +3.845d-06, +1.182d-01, +1.158d-02, ...
                    +5.445d-03, +6.219d-05, +4.204d-06, -2.093d-06, +1.540d-07, ...
                    -4.280d-08, -4.751d-01, -3.490d-02, +1.758d-03, +4.019d-04, ...
                    -2.799d-06, -1.287d-06, +5.468d-07, +7.580d-08, -6.300d-09, ...
                    -1.160d-01, +8.301d-03, +8.771d-04, +9.955d-05, -1.718d-06, ...
                    -2.012d-06, +1.170d-08, +1.790d-08, -1.300d-09, +1.000d-10];
                
                bh_mean = ...
                    [+0.000d+00, +0.000d+00, +3.249d-02, +0.000d+00, +3.324d-02, ...
                    +1.850d-02, +0.000d+00, -1.115d-01, +2.519d-02, +4.923d-03, ...
                    +0.000d+00, +2.737d-02, +1.595d-02, -7.332d-04, +1.933d-04, ...
                    +0.000d+00, -4.796d-02, +6.381d-03, -1.599d-04, -3.685d-04, ...
                    +1.815d-05, +0.000d+00, +7.033d-02, +2.426d-03, -1.111d-03, ...
                    -1.357d-04, -7.828d-06, +2.547d-06, +0.000d+00, +5.779d-03, ...
                    +3.133d-03, -5.312d-04, -2.028d-05, +2.323d-07, -9.100d-08, ...
                    -1.650d-08, +0.000d+00, +3.688d-02, -8.638d-04, -8.514d-05, ...
                    -2.828d-05, +5.403d-07, +4.390d-07, +1.350d-08, +1.800d-09, ...
                    +0.000d+00, -2.736d-02, -2.977d-04, +8.113d-05, +2.329d-07, ...
                    +8.451d-07, +4.490d-08, -8.100d-09, -1.500d-09, +2.000d-10];
                
                ah_amp = ...
                    [-2.738d-01, -2.837d+00, +1.298d-02, -3.588d-01, +2.413d-02, ...
                    +3.427d-02, -7.624d-01, +7.272d-02, +2.160d-02, -3.385d-03, ...
                    +4.424d-01, +3.722d-02, +2.195d-02, -1.503d-03, +2.426d-04, ...
                    +3.013d-01, +5.762d-02, +1.019d-02, -4.476d-04, +6.790d-05, ...
                    +3.227d-05, +3.123d-01, -3.535d-02, +4.840d-03, +3.025d-06, ...
                    -4.363d-05, +2.854d-07, -1.286d-06, -6.725d-01, -3.730d-02, ...
                    +8.964d-04, +1.399d-04, -3.990d-06, +7.431d-06, -2.796d-07, ...
                    -1.601d-07, +4.068d-02, -1.352d-02, +7.282d-04, +9.594d-05, ...
                    +2.070d-06, -9.620d-08, -2.742d-07, -6.370d-08, -6.300d-09, ...
                    +8.625d-02, -5.971d-03, +4.705d-04, +2.335d-05, +4.226d-06, ...
                    +2.475d-07, -8.850d-08, -3.600d-08, -2.900d-09, +0.000d+00];
                
                bh_amp = ...
                    [+0.000d+00, +0.000d+00, -1.136d-01, +0.000d+00, -1.868d-01, ...
                    -1.399d-02, +0.000d+00, -1.043d-01, +1.175d-02, -2.240d-03, ...
                    +0.000d+00, -3.222d-02, +1.333d-02, -2.647d-03, -2.316d-05, ...
                    +0.000d+00, +5.339d-02, +1.107d-02, -3.116d-03, -1.079d-04, ...
                    -1.299d-05, +0.000d+00, +4.861d-03, +8.891d-03, -6.448d-04, ...
                    -1.279d-05, +6.358d-06, -1.417d-07, +0.000d+00, +3.041d-02, ...
                    +1.150d-03, -8.743d-04, -2.781d-05, +6.367d-07, -1.140d-08, ...
                    -4.200d-08, +0.000d+00, -2.982d-02, -3.000d-03, +1.394d-05, ...
                    -3.290d-05, -1.705d-07, +7.440d-08, +2.720d-08, -6.600d-09, ...
                    +0.000d+00, +1.236d-02, -9.981d-04, -3.792d-05, -1.355d-05, ...
                    +1.162d-06, -1.789d-07, +1.470d-08, -2.400d-09, -4.000d-10];
                
                aw_mean = ...
                    [+5.640d+01, +1.555d+00, -1.011d+00, -3.975d+00, +3.171d-02, ...
                    +1.065d-01, +6.175d-01, +1.376d-01, +4.229d-02, +3.028d-03, ...
                    +1.688d+00, -1.692d-01, +5.478d-02, +2.473d-02, +6.059d-04, ...
                    +2.278d+00, +6.614d-03, -3.505d-04, -6.697d-03, +8.402d-04, ...
                    +7.033d-04, -3.236d+00, +2.184d-01, -4.611d-02, -1.613d-02, ...
                    -1.604d-03, +5.420d-05, +7.922d-05, -2.711d-01, -4.406d-01, ...
                    -3.376d-02, -2.801d-03, -4.090d-04, -2.056d-05, +6.894d-06, ...
                    +2.317d-06, +1.941d+00, -2.562d-01, +1.598d-02, +5.449d-03, ...
                    +3.544d-04, +1.148d-05, +7.503d-06, -5.667d-07, -3.660d-08, ...
                    +8.683d-01, -5.931d-02, -1.864d-03, -1.277d-04, +2.029d-04, ...
                    +1.269d-05, +1.629d-06, +9.660d-08, -1.015d-07, -5.000d-10];
                
                bw_mean = ...
                    [+0.000d+00, +0.000d+00, +2.592d-01, +0.000d+00, +2.974d-02, ...
                    -5.471d-01, +0.000d+00, -5.926d-01, -1.030d-01, -1.567d-02, ...
                    +0.000d+00, +1.710d-01, +9.025d-02, +2.689d-02, +2.243d-03, ...
                    +0.000d+00, +3.439d-01, +2.402d-02, +5.410d-03, +1.601d-03, ...
                    +9.669d-05, +0.000d+00, +9.502d-02, -3.063d-02, -1.055d-03, ...
                    -1.067d-04, -1.130d-04, +2.124d-05, +0.000d+00, -3.129d-01, ...
                    +8.463d-03, +2.253d-04, +7.413d-05, -9.376d-05, -1.606d-06, ...
                    +2.060d-06, +0.000d+00, +2.739d-01, +1.167d-03, -2.246d-05, ...
                    -1.287d-04, -2.438d-05, -7.561d-07, +1.158d-06, +4.950d-08, ...
                    +0.000d+00, -1.344d-01, +5.342d-03, +3.775d-04, -6.756d-05, ...
                    -1.686d-06, -1.184d-06, +2.768d-07, +2.730d-08, +5.700d-09];
                
                aw_amp = ...
                    [+1.023d-01, -2.695d+00, +3.417d-01, -1.405d-01, +3.175d-01, ...
                    +2.116d-01, +3.536d+00, -1.505d-01, -1.660d-02, +2.967d-02, ...
                    +3.819d-01, -1.695d-01, -7.444d-02, +7.409d-03, -6.262d-03, ...
                    -1.836d+00, -1.759d-02, -6.256d-02, -2.371d-03, +7.947d-04, ...
                    +1.501d-04, -8.603d-01, -1.360d-01, -3.629d-02, -3.706d-03, ...
                    -2.976d-04, +1.857d-05, +3.021d-05, +2.248d+00, -1.178d-01, ...
                    +1.255d-02, +1.134d-03, -2.161d-04, -5.817d-06, +8.836d-07, ...
                    -1.769d-07, +7.313d-01, -1.188d-01, +1.145d-02, +1.011d-03, ...
                    +1.083d-04, +2.570d-06, -2.140d-06, -5.710d-08, +2.000d-08, ...
                    -1.632d+00, -6.948d-03, -3.893d-03, +8.592d-04, +7.577d-05, ...
                    +4.539d-06, -3.852d-07, -2.213d-07, -1.370d-08, +5.800d-09];
                
                bw_amp = ...
                    [+0.000d+00, +0.000d+00, -8.865d-02, +0.000d+00, -4.309d-01, ...
                    +6.340d-02, +0.000d+00, +1.162d-01, +6.176d-02, -4.234d-03, ...
                    +0.000d+00, +2.530d-01, +4.017d-02, -6.204d-03, +4.977d-03, ...
                    +0.000d+00, -1.737d-01, -5.638d-03, +1.488d-04, +4.857d-04, ...
                    -1.809d-04, +0.000d+00, -1.514d-01, -1.685d-02, +5.333d-03, ...
                    -7.611d-05, +2.394d-05, +8.195d-06, +0.000d+00, +9.326d-02, ...
                    -1.275d-02, -3.071d-04, +5.374d-05, -3.391d-05, -7.436d-06, ...
                    +6.747d-07, +0.000d+00, -8.637d-02, -3.807d-03, -6.833d-04, ...
                    -3.861d-05, -2.268d-05, +1.454d-06, +3.860d-07, -1.068d-07, ...
                    +0.000d+00, -2.658d-02, -1.947d-03, +7.131d-04, -3.506d-05, ...
                    +1.885d-07, +5.792d-07, +3.990d-08, +2.000d-08, -5.700d-09];
                
                % degree n and order m
                nmax = 9;
                % mmax = 9;
                
                % unit vector
                x = cos(dlat)*cos(dlon);
                y = cos(dlat)*sin(dlon);
                z = sin(dlat);
                
                V = zeros(nmax+1,nmax+1);
                W = zeros(nmax+1,nmax+1);
                
                % Legendre polynomials
                V(1,1) = 1.0D0;
                W(1,1) = 0.0D0;
                V(2,1) = z * V(1,1);
                W(2,1) = 0.0;
                
                for n = 2:nmax
                    V(n+1,1) = ((2*n-1) * z * V(n,1) - (n-1) * V(n-1,1)) / n;
                    W(n+1,1) = 0.0D0;
                end
                for m = 1:nmax
                    V(m+1,m+1) = (2*m-1) * (x*V(m,m) - y*W(m,m));
                    W(m+1,m+1) = (2*m-1) * (x*W(m,m) + y*V(m,m));
                    if (m < nmax)
                        V(m+2,m+1) = (2*m+1) * z * V(m+1,m+1);
                        W(m+2,m+1) = (2*m+1) * z * W(m+1,m+1);
                    end
                    for n = m+2:nmax
                        V(n+1,m+1) = ((2*n-1)*z*V(n,m+1) - (n+m-1)*V(n-1,m+1)) / (n-m);
                        W(n+1,m+1) = ((2*n-1)*z*W(n,m+1) - (n+m-1)*W(n-1,m+1)) / (n-m);
                    end
                end
                this.V = V;
                this.W = W;
                this.lat = dlat;
                this.lon = dlon;
            end
            
            % hydrostatic
            bh = 0.0029;
            c0h = 0.062;
            if (dlat < 0) % southern hemisphere
                phh  = pi;
                c11h = 0.007;
                c10h = 0.002;
            else                % northern hemisphere
                phh  = 0;
                c11h = 0.005;
                c10h = 0.001;
            end
            ch = c0h + ((cos(doy*2d0*pi + phh)+1)*c11h/2 + c10h)*(1-cos(dlat));
            
            if cached
                ahm = this.ahm;
                aha = this.aha;
            else
                ahm = 0.d0;
                aha = 0.d0;
                i = 0;
                for n = 0:nmax
                    for m = 0:n
                        i = i+1;
                        ahm = ahm + (ah_mean(i)*V(n+1,m+1) + bh_mean(i)*W(n+1,m+1));
                        aha = aha + (ah_amp(i) *V(n+1,m+1) + bh_amp(i) *W(n+1,m+1));
                    end
                end
                this.ahm = ahm;
                this.aha = aha;
            end
            ah  = (ahm + aha*cos(doy*2.d0*pi))*1d-5;
            
            sine   = sin(pi/2 - zd);
            % cose   = cos(pi/2 - zd);
            beta   = bh./( sine + ch  );
            gamma  = ah./( sine + beta);
            topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)));
            gmfh   = topcon ./ (sine+gamma);
            
            % height correction for hydrostatic mapping function from Niell (1996)
            a_ht = 2.53d-5;
            b_ht = 5.49d-3;
            c_ht = 1.14d-3;
            hs_km  = dhgt/1000.d0;
            
            beta   = b_ht ./ ( sine + c_ht );
            gamma  = a_ht ./ ( sine + beta);
            topcon = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)));
            ht_corr_coef = 1 ./ sine - topcon ./ (sine + gamma);
            ht_corr      = ht_corr_coef * hs_km;
            gmfh         = gmfh + ht_corr;
            
            % wet
            bw = 0.00146;
            cw = 0.04391;
            
            if cached
                awm = this.awm;
                awa = this.awa;
            else
                awm = 0.d0;
                awa = 0.d0;
                i = 0;
                for n = 0 : nmax
                    for m = 0 : n
                        i = i+1;
                        awm = awm + (aw_mean(i)*V(n+1,m+1) + bw_mean(i)*W(n+1,m+1));
                        awa = awa + (aw_amp(i) *V(n+1,m+1) + bw_amp(i) *W(n+1,m+1));
                    end
                end
                this.awm = awm;
                this.awa = awa;
            end
            aw =  (awm + awa*cos(doy*2*pi))*1d-5;
            
            beta   = bw ./ ( sine + cw );
            gamma  = aw ./ ( sine + beta);
            topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)));
            gmfw   = topcon ./ (sine+gamma);
        end
    end
    
    methods (Static)
        
        %-----------------------------------------------------------
        % IONO
        %-----------------------------------------------------------
        function [delay] = klobucharModel(lat, lon, az, el, sow, ionoparams)
            % SYNTAX:
            %   [delay] = Atmosphere. klobuchar_model(lat, lon, az, el, sow, ionoparams)
            %
            % INPUT:
            %   lat = receiver latitude          [degrees] [nx1]
            %   lon = receiver longitude         [degrees] [nx1]
            %   az  = satellite azimuth          [degrees] [nx1]
            %   el = satellite elevation         [degrees] [nx1]
            %   sow = second of week                       [nx1]
            % OUTPUT:
            %   corr = tropospheric error correction
            %
            % DESCRIPTION:
            %   Computation of the pseudorange correction due to ionosphere.
            %   --> multiple epoch for both static and moving target
            %-------------------------------------------------------------------------------
            % KLOBUCHAR MODEL
            %
            % Algorithm from Leick, A. (2004) "GPS Satellite Surveying - 2nd Edition"
            % John Wiley & Sons, Inc., New York, pp. 301-303)
            %-------------------------------------------------------------------------------
            %initialization
            delay = zeros(size(el));
            
            
            
            %ionospheric parameters
            a0 = ionoparams(1);
            a1 = ionoparams(2);
            a2 = ionoparams(3);
            a3 = ionoparams(4);
            b0 = ionoparams(5);
            b1 = ionoparams(6);
            b2 = ionoparams(7);
            b3 = ionoparams(8);
            
            %elevation from 0 to 90 degrees
            el = abs(el);
            
            %conversion to semicircles
            lat = lat / 180;
            lon = lon / 180;
            az = az / 180;
            el = el / 180;
            
            f = 1 + 16*(0.53-el).^3;
            
            psi = (0.0137 ./ (el+0.11)) - 0.022;
            
            phi = lat + psi .* cos(az*pi);
            phi(phi > 0.416)  =  0.416;
            phi(phi < -0.416) = -0.416;
            
            lambda = lon + ((psi.*sin(az*pi)) ./ cos(phi*pi));
            
            ro = phi + 0.064*cos((lambda-1.617)*pi);
            
            t = lambda*43200 + sow;
            t = mod(t,86400);
            
            
            a = a0 + a1*ro + a2*ro.^2 + a3*ro.^3;
            a(a < 0) = 0;
            
            p = b0 + b1*ro + b2*ro.^2 + b3*ro.^3;
            p(p < 72000) = 72000;
            
            x = (2*pi*(t-50400)) ./ p;
            
            %ionospheric delay
            index = find(abs(x) < 1.57);
            delay(index,1) = goGNSS.V_LIGHT * f(index) .* (5e-9 + a(index) .* (1 - (x(index).^2)/2 + (x(index).^4)/24));
            
            index = find(abs(x) >= 1.57);
            delay(index,1) = goGNSS.V_LIGHT * f(index) .* 5e-9;
        end
        
        function [lat_pp, lon_pp, iono_mf, k] = getPiercePoint(lat_rad, lon_rad, h_ortho, az_rad, el_rad, thin_shell_height, rcm)
            % Get the pierce point
            % INPUT:
            %   lat_rad             latitude of the receiver           [rad]
            %   lon_rad             longitude of the receiver          [rad]
            %   h_ortho             orthometric height of the receiver [m]
            %   az_rad              azimuth of the satellites          [rad]
            %   el_rad              elevation of the satellites        [rad]
            %   thin_shell_height   height of the pierce point         [m]
            %   rcm                 meridian radius curvature <optional>
            %
            % OUTPUT
            %   latpp               latitude pierce point [rad]
            %   lonpp               longitude pierce point [rad]
            %   iono_mf             iono mapping function
            %   k                   iono k factor ????
            %
            % SYNTAX: 
            %   [latpp, lonpp, mfpp, k] = getPiercePoint(lat_rad, lon_rad, h_ortho, az_rad, el_rad, thin_shell_height, <rcm>)
            
            % Get radius of curvature at lat
            if nargin < 7
                rcm = getMeridianRadiusCurvature(lat_rad);
            end
            
            input_size = (size(az_rad));
            az_rad = az_rad(:);
            el_rad = el_rad(:);
            
            k = ((rcm + h_ortho)/((rcm + h_ortho) + thin_shell_height)) * cos(el_rad);
            psi_pp = (pi/2) - el_rad - asin(k);
            
            %set azimuth from -180 to 180
            az_rad = mod((az_rad+pi),2*pi)-pi;
            
            %latitude of the ionosphere piercing point
            lat_pp = asin(sin(lat_rad) * cos(psi_pp) + cos(lat_rad) * sin(psi_pp) .* cos(az_rad));
            
            %longitude of the ionosphere piercing point
            id_hl = ((lat_pp >  70*pi/180) & (tan(psi_pp).*cos(az_rad)      > tan((pi/2) - lat_rad))) | ...
                ((lat_pp < -70*pi/180) & (tan(psi_pp).*cos(az_rad + pi) > tan((pi/2) + lat_rad)));
            
            lon_pp = zeros(size(az_rad));
            lon_pp(id_hl) = lon_rad + pi - asin(sin(psi_pp(id_hl)) .* sin(az_rad(id_hl)) ./ cos(lat_pp(id_hl)));

            lon_pp(~id_hl) = lon_rad + asin(sin(psi_pp(~id_hl)) .* sin(az_rad(~id_hl)) ./ cos(lat_pp((~id_hl))));
            
            % using thin shell layer mapping function (Handbook of Global
            % Navigation System pp 185)
            if nargout > 2
                iono_mf = (1-(k)^2)^(-1/2);
            end
            
            if nargout > 3
                k = [cos(az_rad).*cos(el_rad);
                    -sin(az_rad).*cos(el_rad);
                    sin(el_rad)];
                % go to global system
                R = [-sin(lat_rad) cos(lon_rad) 0;
                    -sin(lat_rad)*cos(lon_rad) -sin(lat_rad)*sin(lon_rad) cos(lat_rad);
                    +cos(lat_rad)*cos(lon_rad) +cos(lat_rad)*sin(lon_rad) sin(lat_rad)];
                [k] = R'*k;
            end
            
            lat_pp = reshape(lat_pp, input_size(1), input_size(2));
            lon_pp = reshape(lon_pp, input_size(1), input_size(2));
        end
    end
end
