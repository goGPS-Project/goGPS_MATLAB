classdef Core_Sky < handle
    
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
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
    
    properties
        time_ref_coord  % Gps Times of tabulated aphemerids
        time_ref_clock
        coord  % cvoordinates of tabulated aphemerids [times x num_sat x 3]
        coord_type % 0: Center of Mass 1: Antenna Phase Center
        clock  % cloks of tabulated aphemerids [times x num_sat]
        prn
        sys
        %time_hr
        %clock_hr
        coord_rate = 900;
        clock_rate = 900;
        iono
        
        X_sun  % coord of tabulated sun positions ECEF at the same time of coord
        X_moon % coord of tabulated moon positions ECEF at the same time of coord
        sun_pol_coeff % coeff for polynoimial interpolation of tabulated sun positions
        moon_pol_coeff % coeff for polynoimial interpolation of tabulated moon positions
        
        ERP  % EARH rotation parameters
        group_delays_flags = [ 'GC1C' ; 'GC1S' ; 'GC1L' ; 'GC1X' ; 'GC1P' ; 'GC1W' ; 'GC1Y' ; 'GC1M' ; 'GC2C' ; 'GC2D' ; 'GC2S' ; 'GC2L' ; 'GC2X' ; 'GC2P' ; 'GC2W' ; 'GC2Y' ; 'GC2M' ; 'GC5I' ; 'GC5Q' ; 'GC5X' ; ... %% GPS codes
                               'RC1C' ; 'RC1P' ; 'RC2C' ; 'RC2P' ; 'RC3I' ; 'RC3Q' ; 'RC3X' ; ... %GLONASS code
                               'EC1A' ; 'EC1B' ; 'EC1C' ; 'EC1X' ; 'EC1Z' ; 'EC5I' ; 'EC5Q' ; 'EC5X' ; 'EC7I' ; 'EC7Q' ; 'EC7X' ; 'EC8I' ; 'EC8Q' ; 'EC8X' ; ... %GALIELEO codes
                               'BC2I' ; 'BC2Q' ; 'BC2X' ; 'BC7I' ; 'BC7Q' ; 'BC7X' ; 'BC6I' ; 'BC6Q' ; 'BC6X' ; ... %BeiDou codes
                               'QC1C' ; 'QC1S' ; 'QC1L' ; 'QC1X' ; 'QC1Z' ; 'QC2S' ; 'QC2L' ; 'QC2X' ; 'QC2M' ; 'QC5I' ; 'QC5Q' ; 'QC5X' ; 'QC6S' ; 'QC6L' ; 'QC6X' ; ... %% QZSS codes
                               'IC5A' ; 'IC5B' ; 'IC5C' ; 'IC5X' ; 'IC9A' ; 'IC9B' ; 'IC9C' ; 'IC9X' ; ... %% IRNSS codes
                               'SC1C' ; 'SC5I' ; 'SC5Q' ; 'SC5X' % SBAS   
                               ]; % ALL Rinex 3 code observations flags
        group_delays = zeros(32,77); % group delay of code measurements (meters) referenced to their constellation reference:
                                     %    GPS -> Iono free linear combination C1P C2P
                                     %    GLONASS -> Iono free linear combination C1P C2P
                                     %    Galileo -> Iono free linear combination
                                     %    BedDou -> B3 signal
                                     %    QZS -> Iono free linear combination
                                     %    IRNSS -> Iono free linear combination
                                     %    SABS -> Iono free linear combination
        group_delays_times % 77x1 GPS_Time
        antenna_PCO   %% satellites antenna phase center offset
        antenna_PCV  %% satellites antenna phase center variations
        %satType
        avail
        coord_pol_coeff %% coefficient of the polynomial interpolation for coordinates [11,3,num_sat,num_coeff_sets]
        %start_time_idx  %%% index  defininig start of the day (time == 1)
        cc % constellation collector
    end
    properties (Access = private)
        
        log
        state
    end
    methods (Access = 'private')
        % Creator
        function this = Core_Sky()
            % Core object creator
            this.state = Go_State.getCurrentSettings();
            this.log = Logger.getInstance();
            this.cc = Go_State.getCurrentSettings().getConstellationCollector();
            this.antenna_PCO =zeros(1,this.cc.getNumSat(),3);
            
        end
    end
    
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_core_sky__
            
            if isempty(unique_instance_core_sky__)
                this = Core_Sky();
                unique_instance_core_sky__ = this;
            else
                this = unique_instance_core_sky__;
            end
        end
    end
    % =========================================================================
    %  METHODS
    % =========================================================================
    
    methods % Public Access
        function init(this)
            state = Go_State.getCurrentSettings();
            %%% might serve some scope
        end
        function initSession(this, start_date, stop_time)
            state = Go_State.getCurrentSettings();
            
            %%% load Epehemerids
            eph_f_name = state.getEphFileName(start_date, stop_time);
            clock_f_name = state.getClkFileName(start_date, stop_time);
            clock_in_eph = isempty(setdiff(eph_f_name,clock_f_name)); %%% condition to be tested in differnet cases
            this.clearOrbit();
            if strfind(eph_f_name{1},'.sp3') % assuming all files have the same endings
                for i = 1:length(eph_f_name)
                    this.addSp3(eph_f_name{i},clock_in_eph);
                    this.coord_type = 0; %center of mass
                end
            else %% if not sp3 assume is a rinex navigational file
                this.importBrdcs(eph_f_name,start_date, stop_time, clock_in_eph);
            end
            this.coord_type = 1; %Antenna phase center
            if not(clock_in_eph)
                for i = 1:length(clock_f_name)
                    [~,~,ext] = fileparts(clock_f_name{i});
                    str = strsplit(ext,'_');
                    clk_rate = str2num(str{2}(1:end-1));
                    this.addClk(clock_f_name{i});
                end
            end
            %%% compute pol for ephemerids
            this.computeSatPolyCoeff();
            this.computeSMPolyCoeff();
            %%% load PCV
            this.load_antenna_PCV(state.getAtxFile);
            % pass to antenna phase center if necessary
            if this.coord_type == 0
                this.toAPC();
            end
            %%% load ERP
            this.importERP(state.getErpFileName(start_date, stop_time),start_date);
            %%% load DCB
            this.importCODEDCB();
        end
        function clearOrbit(this, gps_date)
            if nargin > 1
                this.clearCoord(gps_date);
                this.clearClock(gps_date);
            else
                this.clearCoord();
                this.clearClock();
                this.clearSunMoon();
            end
        end
        function clearCoord(this, gps_date)
            % DESCRIPTION: clear coord data , if date is provided clear
            % only data before that date
            if nargin > 1
                if this.time_ref_coord < gps_date
                    n_ep = floor((gps_date - this.time_ref_coord)/this.coord_rate);
                    this.coord(1:n_ep,:,:)=[];
                    this.time_ref_coord.addSeconds(n_ep*this.coord_rate);
                    this.coord_pol_coeff = []; %!!! the coefficient have to been recomputed
                    
                    % deleate also sun e moon data
                    if not(isempty(this.X_sun))
                        this.X_sun(1:n_ep,:)=[];
                    end
                    if not(isempty(this.X_moon))
                        this.X_moon(1:n_ep,:)=[];
                    end
                    this.sun_pol_coeff = []; %!!! the coefficient have to been recomputed
                    this.moon_pol_coeff = []; %!!! the coefficient have to been recomputed
                    
                end
            else
                this.coord=[];
                this.time_ref_coord = [];
                this.coord_pol_coeff = [];
            end
        end
        function clearClock(this, gps_date)
            % DESCRIPTION: clear clock data , if date is provided clear
            % only data before that date
            if nargin > 1
                if this.time_ref_clock < gps_date
                    n_ep = floor((gps_date - this.time_ref_clock)/this.clock_rate);
                    this.clock(1:n_ep,:)=[];
                    this.time_ref_clock.addSeconds(n_ep*this.clock_rate);
                    
                    
                end
            else
                this.clock=[];
                this.time_ref_clock = [];
            end
        end
        function clearSunMoon(this, gps_date)
            % DESCRIPTION: clear sun and moon data , if date is provided clear
            % only data before that date
            if nargin > 1
                if this.time_ref_coord > gps_date
                    n_ep = floor((gps_date - this.time_ref_coord)/this.coord_rate);
                    this.X_sun(1:n_ep,:)=[];
                    this.X_moon(1:n_ep,:)=[];
                    this.sun_pol_coeff = []; %!!! the coefficient have to been recomputed
                    this.moon_pol_coeff = []; %!!! the coefficient have to been recomputed
                end
            end
            this.X_sun = [];
            this.X_moon = [];
            this.sun_pol_coeff = [];
            this.moon_pol_coeff = [];
        end
        function orb_time = getCoordTime(this)
            % DESCRIPTION:
            % return the time of coordinates in GPS_Time (unix time)
            orb_time = this.time_ref_coord.getCopy();
            orb_time.toUnixTime();
            [r_u_t , r_u_t_f ] = orb_time.getUnixTime();
            
            dt = [this.coord_rate : this.coord_rate : (size(this.coord,1)-1)*this.coord_rate]';
            
            
            u_t = r_u_t + uint32(fix(dt));
            u_t_f =  r_u_t_f  + rem(dt,1);
            
            idx = u_t_f >= 1;
            
            u_t(idx) = u_t(idx) + 1;
            u_t_f(idx) = u_t_f(idx) - 1;
            
            idx = u_t_f < 0;
            
            u_t(idx) = u_t(idx) - 1;
            u_t_f(idx) = 1 + u_t_f(idx);
            
            orb_time.appendUnixTime(u_t , u_t_f);
        end
        function orb_time = getClockTime(this)
            % DESCRIPTION:
            % return the time of clock corrections in GPS_Time (unix time)
            orb_time = this.time_ref_clock.getCopy();
            orb_time.toUnixTime();
            
            [r_u_t , r_u_t_f ] = orb_time.getUnixTime();
            
            
            dt = [this.clock_rate : this.clock_rate : (size(this.clock,1)-1)*this.clock_rate]';
            
            
            u_t = r_u_t + uint32(fix(dt));
            u_t_f =  r_u_t_f  + rem(dt,1);
            
            idx = u_t_f >= 1;
            
            u_t(idx) = u_t(idx) + 1;
            u_t_f(idx) = u_t_f(idx) - 1;
            
            idx = u_t_f < 0;
            
            u_t(idx) = u_t(idx) - 1;
            u_t_f(idx) = 1 + u_t_f(idx);
            
            orb_time.appendUnixTime(u_t , u_t_f);
            
        end
        function importEph(this, eph, t_st, t_end, step, clock)
            % SYNTAX:
            %   eph_tab.importEph(eph, t_st, t_end, sat, step)
            %
            % INPUT:
            %   eph         = ephemerids matrix
            %   t_st        = start_time
            %   t_end       = start_time
            %   sat         = available satellite indexes
            %
            % OUTPUT:
            %   XS      = satellite position at time in ECEF(time_rx) (X,Y,Z)
            %   VS      = satellite velocity at time in ECEF(time_tx) (X,Y,Z)
            %   dtS     = satellite clock error (vector)
            %
            % DESCRIPTION:
            
            if nargin < 5
                step = 900;
            end
            this.coord_rate = step;
            if clock
                this.clock_rate = step;
            end
            if nargin < 4 || t_end.isempty()
                t_end = t_st;
            end
            times = t_st.getGpsTime -5*step : step : t_end.getGpsTime+5*step; %%% compute 5 step before and after the day to use for polynomila interpolation
            this.time_ref_coord = t_st.getCopy();
            this.time_ref_coord.toUnixTime();
            this.time_ref_coord.addSeconds(-5*step);
            if clock
                this.time_ref_clock = this.time_ref_coord.getCopy();
            end
            
            sat = unique(eph(30,:)); %% keep only satellite also present in eph
            i = 0;
            
            this.coord = zeros(length(times), this.cc.getNumSat,3 );
            this.clock = zeros (this.cc.getNumSat, length(times));
            t_dist_exced=false;
            for t = times
                i=i+1;
                [this.coord(i,:,:), ~, clock_temp, t_d_e]=this.satellitePositions(t, sat, eph); %%%% loss of precision problem should be less tha 1 mm
                if clock
                    this.clock(:,i) = clock_temp;
                end
                t_dist_exced = t_dist_exced || t_d_e;
            end
            if t_dist_exced
                this.log.addWarning(sprintf('One of the time bonds (%s , %s)\ntoo far from valid ephemerids \nPositions might be inaccurate\n ',t_st.toString(0),t_end.toString(0)))
            end
        end
        function importBrdcs(this,f_names, t_st, t_end, clock, step)
            if nargin < 6
                step = 900;
            end
            if nargin < 5
                clock = true
            end
            sat= this.cc.index;
            if nargin < 4 || t_end.isempty()
                t_end = t_st;
            end
            if not(iscell(f_names))
                f_names = {f_names};
            end
            eph = [];
            for i=1:length(f_names)
                [eph_temp, this.iono, flag_return] = load_RINEX_nav(f_names{i},this.cc,0,0);
                eph = [eph eph_temp];
            end
            
            if not(isempty(eph))
                this.importEph(eph, t_st, t_end, step, clock);
            end
            %%% add TGD delay parameter
            for const = unique(eph(31,:))
                eph_const = eph(:,eph(31,:)==const);
                for s = unique(eph_const(1,:))
                    eph_sat = eph_const(:, eph_const(1,:) == s);
                    GD = eph_sat(28,1); % TGD change only every 3 months
                end
                switch char(const) 
                    case 'G'
                            idx_c1w = this.getGroupDelayIdx('GC1W');
                            idx_c2w = this.getGroupDelayIdx('GC2W');
                            this.group_delays(s,idx_c1w) = -GD * goGNSS.V_LIGHT;
                            f = this.cc.getGPS().F_VEC; % frequencies 
                            this.group_delays(s,idx_c2w) = - f(1)^2 / f(2)^2 * GD * goGNSS.V_LIGHT;
                    case 'R'
                            idx_c1p = this.getGroupDelayIdx('RC1P');
                            idx_c2p = this.getGroupDelayIdx('RC2P');
                            this.group_delays(s,idx_c1p) = -GD * goGNSS.V_LIGHT;
                            f = this.cc.getGLONASS().F_VEC; % frequencies 
                            this.group_delays(s,idx_c2p) = - f(1)^2 / f(2)^2 * GD * goGNSS.V_LIGHT;
                    case 'E'
                            idx_c1p = this.getGroupDelayIdx('EC1B');
                            idx_c2p = this.getGroupDelayIdx('EC5I');
                            this.group_delays(s,idx_c1p) = -GD * goGNSS.V_LIGHT;
                            f = this.cc.getGalileo().F_VEC; % frequencies 
                            this.group_delays(s,idx_c2p) = - f(1)^2 / f(2)^2 * GD * goGNSS.V_LIGHT;
                            
                end
            end
                
            
        end
        function [XS,VS,dt_s, t_dist_exced] =  satellitePositions(this,time, sat, eph)
            
            % SYNTAX:
            %   [XS, VS] = satellite_positions(time_rx, sat, eph);
            %
            % INPUT:
            %   time_rx     = reception time
            %   sat         = available satellite indexes
            %   eph         = ephemeris
            %
            % OUTPUT:
            %   XS      = satellite position at time in ECEF(time_rx) (X,Y,Z)
            %   VS      = satellite velocity at time in ECEF(time_tx) (X,Y,Z)
            %   dtS     = satellite clock error (vector)
            %
            % DESCRIPTION:
            nsat = length(sat);
            
            XS = zeros(this.cc.getNumSat(), 3);
            VS = zeros(this.cc.getNumSat(), 3);
            
            
            dt_s = zeros(this.cc.getNumSat(), 1);
            t_dist_exced = false;
            for i = 1 : nsat
                
                k = find_eph(eph, sat(i), time);
                if not(isempty(k))
                    %compute satellite position and velocity
                    [XS(sat(i),:), VS(sat(i),:)] = satellite_orbits(time, eph(:,k), sat(i), []);
                    dt_s(sat(i)) = sat_clock_error_correction(time, eph(:,k));
                    dt_s(sat(i)) = sat_clock_error_correction(time - dt_s(sat(i)), eph(:,k));
                else
                    t_dist_exced = true;
                end
                
            end
            
            
            %XS=XS';
        end
        function addSp3(this, filename_SP3, clock_flag)
            % SYNTAX:
            %   this.addSp3(filename_SP3, clock_flag)
            %
            % INPUT:
            %   filename_SP3 = name of sp3 file
            %   clock_flag   = load also clock? (optional, dafault = true)
            %
            % DESCRIPTION:
            % add satellite and clock postiion contained in the sp3 file to
            % the object if values are contiguos with the ones already in
            % the object add them, otherwise clear the object and add them
            % data that are alrady present are going to be overwritten
            
            if isempty(this.coord)
                empty_file = true;
            else
                empty_file = false;
            end
            if nargin <3
                clock_flag = true;
            end
            k = 0; % current epoch
            
            
            %SP3 file
            f_sp3 = fopen(filename_SP3,'r');
            
            if (f_sp3 ~= -1)
                
                txt = fread(f_sp3,'*char')';
                fclose(f_sp3);
                
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
                % get end pf header
                % coord  rate
                coord_rate = cell2mat(textscan(txt(repmat(lim(2,1),1,11) + [26:36]),'%f'));
                % n epochs
                nEpochs = cell2mat(textscan(txt(repmat(lim(1,1),1,7) + [32:38]),'%f'));
                % find first epoch
                string_time = txt(repmat(lim(1,1),1,28) + [3:30]);
                % convert the times into a 6 col time
                date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %10.8f'));
                % import it as a GPS_Time obj
                sp3_first_ep = GPS_Time(date, [], true);
                if this.coord_rate ~= coord_rate;
                    if empty_file
                        this.coord_rate = coord_rate;
                        if clock_flag
                            this.clock_rate = coord_rate;
                        end
                    else
                        this.log.addWarning(['Coord rate not match: ' num2str(coord_rate)]);
                        return
                    end
                    
                end
                % checking overlapping and same correct syncro
                sp3_last_ep = sp3_first_ep.getCopy();
                sp3_last_ep.addSeconds(coord_rate*nEpochs);
                if ~empty_file
                    idx_first = (sp3_first_ep - this.time_ref_coord)/this.coord_rate;
                    idx_last = (sp3_last_ep - this.time_ref_coord)/this.coord_rate;
                    memb_idx = ismembertol([idx_first idx_last], -1 : (size(this.coord,1)+1) ); %check whether the extend of sp3 file intersect with the current data
                    if sum(memb_idx)==0
                        empty_file = true;
                        this.clearCoord(); %<---- if new sp3 does not match the already present data clear the data and put the new ones
                        %                         elseif sum(memb_idx)==2 %<--- case new data are already in the class, (this leave out the case wether only one epoch more would be added to the current data, extremely unlikely)
                        %                             return
                    end
                end
                %initlaize array size
                if empty_file
                    this.time_ref_coord = sp3_first_ep.getCopy();
                    if clock_flag
                        this.time_ref_clock = sp3_first_ep.getCopy();
                    end
                    this.coord = zeros(nEpochs, this.cc.getNumSat(),3);
                    if clock_flag
                        this.clock = zeros(nEpochs, this.cc.getNumSat());
                    end
                else
                    c_n_sat = size(this.coord,2);
                    if memb_idx(1) == true & memb_idx(2) == false
                        n_new_epochs = idx_last - size(this.coord, 1);
                        this.coord = cat(1,this.coord,zeros(n_new_epochs,c_n_sat,3));
                        if clock_flag
                            this.clock = cat(1,this.clock,zeros(n_new_epochs,c_n_sat));
                        end
                    elseif memb_idx(1) == false & memb_idx(2) == true
                        this.time_ref_coord = sp3_first_ep.getCopy();
                        if clock_flag
                            this.time_ref_clock = sp3_first_ep.getCopy();
                        end
                        n_new_epochs = -idx_first;
                        this.coord = cat(1,zeros(n_new_epochs,c_n_sat,3),this.coord);
                        if clock_flag
                            this.clock = cat(1,zeros(n_new_epochs,c_n_sat),this.clock);
                        end
                    end
                end
                %%%% read data
                %%% raed epochs
                t_line = find(txt(lim(:,1)) == '*');
                string_time = txt(repmat(lim(t_line,1),1,28) + repmat(3:30, length(t_line), 1))';
                % convert the times into a 6 col time
                date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %10.8f'));
                % import it as a GPS_Time obj
                sp3_times = GPS_Time(date, [], true);
                ant_ids = this.cc.getAntennaId;
                for i = 1 : length(ant_ids)
                    ant_id = ant_ids{i};
                    sat_line = find(txt(lim(:,1)) == 'P' & txt(lim(:,1)+1) == ant_id(1) & txt(lim(:,1)+2) == ant_id(2)& txt(lim(:,1)+3) == ant_id(3));
                    if ~isempty((sat_line))
                        c_ep_idx = round((sp3_times - this.time_ref_coord) / this.coord_rate) +1; %current epoch index
                        this.coord(c_ep_idx,i,:) = cell2mat(textscan(txt(repmat(lim(sat_line,1),1,41) + repmat(5:45, length(sat_line), 1))','%f %f %f'))*1e3;
                        if clock_flag
                            this.clock(c_ep_idx,i) = cell2mat(textscan(txt(repmat(lim(sat_line,1),1,13) + repmat(47:59, length(sat_line), 1))','%f'))/10e6;
                        end
                    else
                    end
                end
            else
                this.log.addWarning([ filename_SP3 ' not found.']);
            end
            clear sp3_file;
            this.coord = zero2nan(this.coord);
        end
        
        function addClk(this,filename_clk)
            % SYNTAX:
            %   eph_tab.addClk(filename_clk)
            %
            % INPUT:
            %   filename_clk = name of clk rinex file file (IMPORTANT:the method
            %   assume 1 day clock filen at 5s)
            %
            % DESCRIPTION:
            % add satellites  clock contained in the clk file to
            % the object if values are contiguos with the ones already in
            % the object add them, otherwise clear the object and add them
            % data that are alrady present are going to be overwritten
            f_clk = fopen(filename_clk,'r');
            
            if (f_clk == -1)
                this.log.addWarning(sprintf('No clk files have been found at %s', filename_clk));
            else
                this.log.addMessage(sprintf('Opening file %s for reading', filename_clk));
                t0 = tic;
                if isempty(this.clock)
                    empty_clk = true;
                else
                    empty_clk = false;
                    [ref_week, ref_sow] =this.time_ref_clock.getGpsWeek();
                    
                end
                
                % open RINEX observation file
                fid = fopen(filename_clk,'r');
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
                % get end pf header
                eoh = strfind(txt,'END OF HEADER');
                eoh = find(lim(:,1)>eoh);
                eoh = eoh(1)-1;
                sats_line = find(txt(lim(eoh+1:end,1)) == 'A' & txt(lim(eoh+1:end,1)+1) == 'S') + eoh;
                % clk rate
                clk_rate = [];
                % find first epoch
                string_time = txt(repmat(lim(sats_line(1),1),1,27) + [8:34]);
                % convert the times into a 6 col time
                date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %10.7f'));
                % import it as a GPS_Time obj
                file_first_ep = GPS_Time(date, [], true);
                % find sampling rate
                ant_ids = this.cc.getAntennaId;
                for i = 1 : length(ant_ids)
                    ant_id = ant_ids{i};
                    sat_line = sats_line(txt(lim(sats_line,1)+3) == ant_id(1) & txt(lim(sats_line,1)+4) == ant_id(2)& txt(lim(sats_line,1)+5) == ant_id(3));
                    if not(isempty(sat_line))
                        n_ep_sat = length(sat_line);
                        string_time = txt(repmat(lim(sat_line,1),1,27) + repmat(8:34, n_ep_sat, 1))';
                        % convert the times into a 6 col time
                        date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %10.7f'));
                        % import it as a GPS_Time obj
                        sat_time = GPS_Time(date, [], true);
                        
                        % initilize matrix
                        if isempty(clk_rate)
                            clk_rate = mean(diff(sat_time.getGpsTime()));
                            if not(empty_clk) & clk_rate ~= this.clock_rate
                                this.log.addWarning('Clock rate in file different from one in Core_Sky\n Discarding old data\n');
                                thic.clearClock();
                                empty_clk = true;
                                this.clock_rate = clk_rate;
                            end
                            if empty_clk
                                this.clock_rate = clk_rate;
                                this.time_ref_clock = file_first_ep;
                                [ref_week, ref_sow] =this.time_ref_clock.getGpsWeek();
                                
                                this.clock = zeros(86400 / this.clock_rate,this.cc.getNumSat());
                            else
                                
                                c_ep_idx = round((file_first_ep - this.time_ref_clock) / this.clock_rate) +1; % epoch index
                                if c_ep_idx < 1
                                    this.clock = [zeros(abs(c_ep_idx)+1,size(this.clock,2)); this.clock];
                                    this.time_ref_clock = file_first_ep;
                                    [ref_week, ref_sow] =this.time_ref_clock.getGpsWeek();
                                end
                                c_ep_idx = round((file_first_ep - this.time_ref_clock) / this.clock_rate) +1; % epoch index
                                if c_ep_idx + 86400/this.clock_rate -1 > size(this.clock,1)
                                    this.clock = [this.clock; zeros( c_ep_idx + 86400/this.clock_rate -1 - size(this.clock,1) ,size(this.clock,2)); ];
                                end
                            end
                        end
                        c_ep_idx = round((sat_time - this.time_ref_clock) / this.clock_rate) +1; % epoch index
                        this.clock(c_ep_idx,i) =  cell2mat(textscan(txt(repmat(lim(sat_line,1),1,20) + repmat(41:60, n_ep_sat, 1))','%f'));
                        
                        
                    end
                end
                this.log.addMessage(sprintf('Parsing completed in %.2f seconds', toc(t0)));
                this.log.newLine();
            end
        end
        function importERP(this, f_name, time)
            this.ERP = load_ERP(f_name, time.getGpsTime());
        end
        function importCODEDCB(this)
            [DCB] = load_dcb(this.state.DCB_DIR, double(this.time_ref_coord.getGpsWeek), this.time_ref_coord.getGpsTime, true, goGNSS.initConstellation(true , true, true,true,true,true));
            %%% assume that CODE DCB contains only GPS and GLONASS
            %GPS C1W - C2W
            idx_w1 =  this.getGroupDelayIdx('GC1W');
            idx_w2 =  this.getGroupDelayIdx('GC2W');
            p1p2 = DCB.P1P2.value(DCB.P1P2.sys == 'G');
            iono_free = this.cc.getGPS.getIonoFree();
            this.group_delays(DCB.P1P2.prn(DCB.P1P2.sys == 'G') , idx_w1) = iono_free.alpha2 *p1p2*goGNSS.V_LIGHT;
            this.group_delays(DCB.P1P2.prn(DCB.P1P2.sys == 'G') , idx_w2) = iono_free.alpha1 *p1p2*goGNSS.V_LIGHT;
            % GPS C1W - C1C
            idx_w1 =  this.getGroupDelayIdx('GC1C');
            idx_w2 =  this.getGroupDelayIdx('GC2D');
            p1c1 = DCB.P1P2.value(DCB.P1C1.sys == 'G');
            this.group_delays(DCB.P1P2.prn(DCB.P1P2.sys == 'G') , idx_w1) = (iono_free.alpha2 *p1p2 + p1c1)*goGNSS.V_LIGHT;
            this.group_delays(DCB.P1P2.prn(DCB.P1P2.sys == 'G') , idx_w2) = (iono_free.alpha1 *p1p2 + p1c1)*goGNSS.V_LIGHT; %semi codeless tracking
            %GLONASS C1P - C2P
            idx_w1 =  this.getGroupDelayIdx('RC1P');
            idx_w2 =  this.getGroupDelayIdx('RC2P');
            p1p2 = DCB.P1P2.value(DCB.P1P2.sys == 'R');
            iono_free = this.cc.getGLONASS.getIonoFree();
            this.group_delays(DCB.P1P2.prn(DCB.P1P2.sys == 'R') , idx_w1) = (iono_free.alpha2 *p1p2)*goGNSS.V_LIGHT;
            this.group_delays(DCB.P1P2.prn(DCB.P1P2.sys == 'R') , idx_w2) = (iono_free.alpha1 *p1p2)*goGNSS.V_LIGHT;
            
            
            
        end
        function importSinexDCB(this, filename)
            %DESCRIPTION: import DCB in sinex format 
            % TBD
        end
        function idx = getGroupDelayIdx(this,flag)
            %DESCRIPTION: get the index of the gorup delay for the given
            %flag
            idx = find(sum(this.group_delays_flags == repmat(flag,size(this.group_delays_flags,1),1),2)==4);
        end
        function importIono(this,f_name)
            [~, this.iono, flag_return ] = load_RINEX_nav(f_name,this.cc,0,0);
            if (flag_return)
                return
            end
        end
        function [sx ,sy, sz] = getSatFixFrame(this,time)
            
            % SYNTAX:
            %   [i, j, k] = satellite_fixed_frame(time,X_sat);
            %
            % INPUT:
            %   time     = GPS_Time [nx1]
            %   X_sat    = postition of satellite [n_epoch x n-sat x 3]
            % OUTPUT:
            %   sx = unit vector that completes the right-handed system [n_epoch x n_sat x 3]
            %   sy = resulting unit vector of the cross product of k vector with the unit vector from the satellite to Sun [n_epoch x n_sat x 3]
            %   sz = unit vector pointing from the Satellite Mass Centre (MC) to the Earth's centre [n_epoch x n_sat x 3]
            %
            % DESCRIPTION:
            %   Computation of the unit vectors defining the satellite-fixed frame.
            
            
            t_sun = time;
            X_sun = this.sunMoonInterpolate(t_sun, true);
            
            X_sat = this.coordInterpolate(time);
            n_sat = size(X_sat,2);
            sx = zeros(size(X_sat)); sy = sx; sz = sx;
            for idx = 1 : t_sun.length()
                x_sun = X_sun(idx,:);
                x_sat = X_sat(idx,:,:);
                e = permute(repmat(x_sun,1,1,n_sat),[1 3 2]) - x_sat() ;
                e = e./repmat(normAlngDir(e,3),1,1,3);
                k = -x_sat./repmat(normAlngDir(x_sat,3),1,1,3);
                j=cross(k,e);
                %                 j = [k(2).*e(3)-k(3).*e(2);
                %                     k(3).*e(1)-k(1).*e(3);
                %                     k(1).*e(2)-k(2).*e(1)];
                i=cross(j,k);
                %                 i = [j(2).*k(3)-j(3).*k(2);
                %                     j(3).*k(1)-j(1).*k(3);
                %                     j(1).*k(2)-j(2).*k(1)];
                sx(idx,:,:) = k;
                sy(idx,:,:) = j ./ repmat(normAlngDir(j,3),1,1,3);
                sz(idx,:,:) = i ./ repmat(normAlngDir(i,3),1,1,3);
            end
            function nrm=normAlngDir(A,d)
                nrm=sqrt(sum(A.^2,d));
            end
        end
        function toCOM(this)
            %DESCRIPTION : convert coord to center of mass
            if this.coord_type == 0
                return %already ceneter of amss
            end
            this.coord = this.getCOM();
            this.coord_type = 0;
        end
        function coord = getCOM(this)
            if this.coord_type == 0
                coord = this.coord; %already ceneter of amss
            else
                [sx, sy, sz] = this.getSatFixFrame(this.getCoordTime());
                coord = this.coord - cat(3, sum(repmat(this.antenna_PCO,size(this.coord,1),1,1) .* sx , 3) ...
                    , sum(repmat(this.antenna_PCO,size(this.coord,1),1,1) .* sy , 3) ...
                    , sum(repmat(this.antenna_PCO,size(this.coord,1),1,1) .* sz , 3));
            end
        end
        function toAPC(this)
            %DESCRIPTION : convert coord to center of antenna phase center
            if this.coord_type == 1
                return %already antennna phase center
            end
            this.coord = this.getAPC();
            this.coord_type = 1;
        end
        function coord = getAPC(this)
            if this.coord_type == 1
                coord = this.coord; %already antennna phase center
            end
            [sx, sy, sz] = this.getSatFixFrame(this.getCoordTime());
            coord = this.coord + cat(3, sum(repmat(this.antenna_PCO,size(this.coord,1),1,1) .* sx , 3) ...
                , sum(repmat(this.antenna_PCO,size(this.coord,1),1,1) .* sy , 3) ...
                , sum(repmat(this.antenna_PCO,size(this.coord,1),1,1) .* sz , 3));
        end
        function [dt_S] = clockInterpolate(this,time, sat)
            
            % SYNTAX:
            %   [dt_S_SP3] = interpolate_SP3_clock(time, sat);
            %
            % INPUT:
            %   time  = interpolation timespan GPS_Time
            %   SP3   = structure containing precise ephemeris data
            %   sat   = satellite PRN
            %
            % OUTPUT:
            %   dt_S_SP3  = interpolated clock correction
            %
            % DESCRIPTION:
            %   SP3 (precise ephemeris) clock correction linear interpolation.
            if nargin < 3
                sat= this.cc.index;
            end
            
            
            
            interval = this.clock_rate;
            
            %find the SP3 epoch closest to the interpolation time
            %[~, p] = min(abs(SP3_time - time));
            % speed improvement of the above line
            % supposing SP3_time regularly sampled
            p = (round((time - this.time_ref_clock) / interval) + 1)';
            times = this.getClockTime();
            b =  (times.getSubSet(p)- time)';
            
            %extract the SP3 clocks
            if (b>0)
                SP3_c = cat(3, this.clock(p-1,sat), this.clock(p,sat));
                u = 1 - b/interval;
            else
                SP3_c = cat(3, this.clock(p,sat), this.clock(p+1,sat));
                u = -b/interval;
            end
            
            dt_S  = NaN*ones(size(SP3_c,1),size(SP3_c,2));
            idx=(sum(SP3_c~=0,3) == 2 .* ~any(SP3_c >= 0.999,3))>0;
            dt_S=repmat((1-u)',1,size(SP3_c,2)).*SP3_c(:,:,1) + repmat((u)',1,size(SP3_c,2)).*SP3_c(:,:,2);
            dt_S(not(idx))=NaN;
            
            %             dt_S_SP3=NaN;
            %             if (sum(SP3_c~=0) == 2 && ~any(SP3_c >= 0.999))
            %
            %                 %linear interpolation (clock)
            %                 dt_S_SP3 = (1-u)*SP3_c(1) + u*SP3_c(2);
            %
            %                 %plot([0 1],SP3_c,'o',u,dt_S_SP3,'.')
            %                 %pause
            %             end
        end
        function computeSatPolyCoeff(this)
            % SYNTAX:
            %   this.computeSatPolyCoeff();
            %
            % INPUT:
            %
            % OUTPUT:
            %
            % DESCRIPTION: Precompute the coefficient of the 10th poynomial for all the possible support sets
            n_pol=10;
            n_coeff=n_pol+1;
            A=zeros(n_coeff,n_coeff);
            A(:,1)=ones(n_coeff,1);
            x=[-5:5]; % *this.coord_rat
            for i=1:10
                A(:,i+1)=(x.^i)';
            end
            n_coeff_set= size(this.coord,1)-10;%86400/this.coord_rate+1;
            %this.coord_pol_coeff=zeros(this.cc.getNumSat,n_coeff_set,n_coeff,3)
            this.coord_pol_coeff=zeros(n_coeff,3,this.cc.getNumSat,n_coeff_set);
            for s=1:this.cc.getNumSat
                for i=1:n_coeff_set
                    for j=1:3
                        this.coord_pol_coeff(:,j,s,i)=A\squeeze(this.coord(i:i+10,s,j));
                    end
                end
            end
        end
        function computeSMPolyCoeff(this)
            % SYNTAX:
            %   this.computeSatPolyCoeff();
            %
            % INPUT:
            %
            % OUTPUT:
            %
            % DESCRIPTION: Precompute the coefficient of the 10th poynomial for all the possible support sets
            n_pol=10;
            n_coeff=n_pol+1;
            A=zeros(n_coeff,n_coeff);
            A(:,1)=ones(n_coeff,1);
            x=[-5:5]; % *this.coord_rat
            for i=1:10
                A(:,i+1)=(x.^i)';
            end
            n_coeff_set= size(this.X_sun,1)-10;%86400/this.coord_rate+1;
            %this.coord_pol_coeff=zeros(this.cc.getNumSat,n_coeff_set,n_coeff,3)
            this.sun_pol_coeff=zeros(n_coeff,3,n_coeff_set);
            this.moon_pol_coeff=zeros(n_coeff,3,n_coeff_set);
            for i=1:n_coeff_set
                for j=1:3
                    this.sun_pol_coeff(:,j,i)=A\squeeze(this.X_sun(i:i+10,j));
                    this.moon_pol_coeff(:,j,i)=A\squeeze(this.X_moon(i:i+10,j));
                end
            end
        end
        function [X_sat, V_sat]=coordInterpolate(this,t,sat)
            % SYNTAX:
            %   [X_sat]=Eph_Tab.polInterpolate(t,sat)
            %
            % INPUT:
            %    t = vector of times where to interpolate
            %    sat = satellite to be interpolated (optional)
            % OUTPUT:
            %
            % DESCRIPTION: interpolate coordinates of staellites
            n_sat=this.cc.getNumSat;
            if nargin <3
                sat_idx=ones(n_sat,1)>0;
            else
                sat_idx=sat;
            end
            
            if isempty(this.coord_pol_coeff)
                this.computeSatPolyCoeff();
            end
            n_sat=length(sat_idx);
            nt=t.length();
            %c_idx=round(t_fd/this.coord_rate)+this.start_time_idx;%coefficient set  index
            
            c_idx = round((t - this.time_ref_coord) / this.coord_rate) - 4;
            
            c_idx(c_idx<1) = 1;
            c_idx(c_idx > size(this.coord,1)-10) = size(this.coord,1)-10;
            
            c_times = this.getCoordTime();
            %l_idx=idx-5;
            %u_id=idx+10;
            
            X_sat=zeros(nt,n_sat,3);
            V_sat=zeros(nt,n_sat,3);
            un_idx=unique(c_idx)';
            for idx=un_idx
                t_idx=c_idx==idx;
                times=t.getSubSet(t_idx);
                %t_fct=((times-this.time(5+idx)))';%time from coefficient time
                t_fct =  (times -  c_times.getSubSet(idx+5))/this.coord_rate; %
                %%%% compute position
                eval_vec = [ones(size(t_fct)) ...
                    t_fct ...
                    t_fct.^2 ...
                    t_fct.^3 ...
                    t_fct.^4 ...
                    t_fct.^5 ...
                    t_fct.^6 ...
                    t_fct.^7 ...
                    t_fct.^8 ...
                    t_fct.^9 ...
                    t_fct.^10];
                X_sat(t_idx,:,:) = reshape(eval_vec*reshape(this.coord_pol_coeff(:,:,sat_idx,idx),11,3*n_sat),sum(t_idx),n_sat,3);
                %%% compute velocity
                eval_vec = [ ...
                    ones(size(t_fct))  ...
                    2*t_fct  ...
                    3*t_fct.^2 ...
                    4*t_fct.^3 ...
                    5*t_fct.^4 ...
                    6*t_fct.^5 ...
                    7*t_fct.^6 ...
                    8*t_fct.^7 ...
                    9*t_fct.^8 ...
                    10*t_fct.^9];
                V_sat(t_idx,:,:) = reshape(eval_vec*reshape(this.coord_pol_coeff(2:end,:,sat_idx,idx),10,3*n_sat),sum(t_idx),n_sat,3)/this.coord_rate;
                
            end
            if size(X_sat,2)==1
                X_sat = squeeze(X_sat);
                V_sat = squeeze(V_sat);
                if size(X_sat,2) ==1
                    X_sat = X_sat';
                    V_sat = V_sat';
                end
            end
        end
        function [sun_ECEF,moon_ECEF ] = sunMoonInterpolate(this,t,no_moon);
            % SYNTAX:
            %   [X_sat]=Eph_Tab.sunInterpolate(t,sat)
            %
            % INPUT:
            %    time = vector of times where to interpolate
            %    no_mmon = do not compute moon postion (default false)
            % OUTPUT:
            %
            % DESCRIPTION: interpolate sun and moon positions
            if isempty(this.X_moon) | isempty(this.X_sun)
                this.tabulateSunMoonPos();
            end
            
            if isempty(this.sun_pol_coeff)
                this.computeSMPolyCoeff();
            end
            
            if nargin < 3
                moon = true;
            else
                moon = not(no_moon);
            end
            %c_idx=round(t_fd/this.coord_rate)+this.start_time_idx;%coefficient set  index
            
            c_idx = round((t - this.time_ref_coord) / this.coord_rate) - 4;
            
            c_idx(c_idx<1) = 1;
            c_idx(c_idx > size(this.coord,1)-10) = size(this.coord,1)-10;
            
            c_times = this.getCoordTime();
            %l_idx=idx-5;
            %u_id=idx+10;
            nt = t.length();
            sun_ECEF=zeros(nt,3);
            if moon
                moon_ECEF=zeros(nt,3);
            end
            un_idx=unique(c_idx)';
            for idx=un_idx
                t_idx=c_idx==idx;
                times=t.getSubSet(t_idx);
                %t_fct=((times-this.time(5+idx)))';%time from coefficient time
                t_fct =  (times -  c_times.getSubSet(idx+5))/this.coord_rate; %
                %%%% compute position
                eval_vec = [ones(size(t_fct)) ...
                    t_fct ...
                    t_fct.^2 ...
                    t_fct.^3 ...
                    t_fct.^4 ...
                    t_fct.^5 ...
                    t_fct.^6 ...
                    t_fct.^7 ...
                    t_fct.^8 ...
                    t_fct.^9 ...
                    t_fct.^10];
                sun_ECEF(t_idx,:) = eval_vec*reshape(this.sun_pol_coeff(:,:,idx),11,3);
                if moon
                    moon_ECEF(t_idx,:) = eval_vec*reshape(this.moon_pol_coeff(:,:,idx),11,3);
                end
            end
        end
        function [sun_ECEF , moon_ECEF] = computeSunMoonPos(this,time,no_moon)
            % SYNTAX:
            %   this.computeSunMoonPos(p_time)
            %
            % INPUT:
            %    time = Gps_Time [n_epoch x 1]
            %    no_moon = do not compute moon (Boolena deafult false)
            % OUTPUT:
            % sun_ECEF  : sun  coordinate Earth Centered Earth Fixed [n_epoch x 3]
            % moon_ECEF : moon coordinate Earth Centered Earth Fixed [n_epoch x 3]
            % DESCRIPTION: Compute sun and moon psitions at the time
            % desidered time
            
            global iephem km ephname inutate psicor epscor ob2000
            %time = GPS_Time((p_time(1))/86400+GPS_Time.GPS_ZERO);
            if nargin < 3
                moon = true;
            else
                moon = not(no_moon);
            end
            
            sun_id = 11; moon_id = 10; earth_id = 3;
            
            readleap; iephem = 1; ephname = 'de436.bin'; km = 1; inutate = 1; ob2000 = 0.0d0;
            
            tmatrix = j2000_icrs(1);
            
            setmod(2);
            % setdt(3020092e-7);
            setdt(5.877122033683494);
            xp = 171209e-6; yp = 414328e-6;
            
            gs = Go_State.getInstance();
            go_dir = gs.getLocalStorageDir();
            
            %if the binary JPL ephemeris file is not available, generate it
            if (exist(fullfile(go_dir, 'de436.bin'),'file') ~= 2)
                fprintf('Warning: file "de436.bin" not found in at %s\n         ... generating a new "de436.bin" file\n',fullfile(go_dir, 'de436.bin'));
                fprintf('         (this procedure may take a while, but it will be done only once on each installation):\n')
                fprintf('-------------------------------------------------------------------\n\n')
                asc2eph(436, {'ascp01950.436', 'ascp02050.436'}, fullfile(go_dir, 'de436.bin'));
                fprintf('-------------------------------------------------------------------\n\n')
            end
            
            sun_ECEF = zeros(time.length(), 3);
            moon_ECEF = zeros(time.length(), 3);
            for e = 1 : time.length()
                [year , month ,day,hour,min,sec]= time.getCalEpoch(e);
                %UTC to TDB
                jdutc = julian(month, day+hour/24+min/1440+sec/86400, year);
                jdtdb = utc2tdb(jdutc);
                %this.t_sun(e) = time.getGpsTime();%(datenum( year,month, day) - GPS_Time.GPS_ZERO)*86400;
                %precise celestial pole (disabled)
                [psicor, epscor] = celpol(jdtdb, 1, 0.0d0, 0.0d0);
                
                %compute the Sun position (ICRS coordinates)
                rrd = jplephem(jdtdb, sun_id, earth_id);
                sun_ECI = rrd(1:3);
                sun_ECI = tmatrix*sun_ECI;
                
                %Sun ICRS coordinates to ITRS coordinates
                deltat = getdt;
                jdut1 = jdutc - deltat;
                tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
                sun_ECEF(e,:) = celter(tjdh, tjdl, xp, yp, sun_ECI)*1e3;
                
                if moon
                    %compute the Moon position (ICRS coordinates)
                    rrd = jplephem(jdtdb, moon_id, earth_id);
                    moon_ECI = rrd(1:3);
                    moon_ECI = tmatrix*moon_ECI;
                    
                    %Moon ICRS coordinates to ITRS coordinates
                    deltat = getdt;
                    jdut1 = jdutc - deltat;
                    tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
                    moon_ECEF(e,:) = celter(tjdh, tjdl, xp, yp, moon_ECI)*1e3;
                end
            end
            
            %this.t_sun_rate =
        end
        function tabulateSunMoonPos(this)
            % SYNTAX:
            %   this.computeSunMoonPos(p_time)
            %
            % INPUT:
            %    p_time = Gps_Time [n_epoch x 1]
            % OUTPUT:
            % DESCRIPTION: Compute sun and moon positions at coordinates time and
            % store them in the object (Overwrite previous data)
            
            %his.t_sun = p_time;
            [this.X_sun , this.X_moon] = this.computeSunMoonPos(this.getCoordTime());
            this.computeSMPolyCoeff();
        end
        function [eclipsed] = check_eclipse_condition(this, time, XS, sat)  %%% TO BE CORRECT
            
            % SYNTAX:
            %   [eclipsed] = check_eclipse_condition(time, XS, SP3, sat, p_rate);
            %
            % INPUT:
            %   time     = GPS time
            %   XS       = satellite position (X,Y,Z)
            %   sat      = staellite goIds
            % OUTPUT:
            %   eclipsed = boolean value to define satellite eclipse condition (0: OK, 1: eclipsed)
            %
            % DESCRIPTION:
            %   Check if the input satellite is under eclipse condition.
            
            eclipsed = zeros(time.length(),1) > 1;
            
            X_sun = this.sunMoonInterpolate(time,true);
            
            %satellite geocentric position
            XS_n = sqrt(sum(XS.^2,2));
            XS_u = XS ./ repmat(XS_n,1,3);
            
            %sun geocentric position
            X_sun_n = sqrt(sum(X_sun.^2,2));
            X_sun_u = X_sun ./ repmat(X_sun_n,1,3);
            
            %satellite-sun angle
            %cosPhi = dot(XS_u, X_sun_u);
            % speed improvement of the above line
            cosPhi = sum(XS_u.*X_sun_u,2);
            
            
            %threshold to detect noon/midnight maneuvers
            if sat > 32
                t = 0; %ignore noon/midnight maneuvers for other constellations (TBD)
            elseif (~isempty(strfind(this.antenna_PCV(sat).sat_type,'BLOCK IIA')))
                t = 4.9*pi/180; % maximum yaw rate of 0.098 deg/sec (Kouba, 2009)
                
            elseif (~isempty(strfind(this.antenna_PCV(sat).sat_type,'BLOCK IIR')))
                t = 2.6*pi/180; % maximum yaw rate of 0.2 deg/sec (Kouba, 2009)
            elseif (~isempty(strfind(this.antenna_PCV(sat).sat_type,'BLOCK IIF')))
                t = 4.35*pi/180; % maximum yaw rate of 0.11 deg/sec (Dilssner, 2010)
            end
            
            if (sat <= 32 & ~isempty(strfind(this.antenna_PCV(sat).sat_type,'BLOCK IIA')))
                %shadow crossing affects only BLOCK IIA satellites
                shadowCrossing = cosPhi < 0 & XS_n.*sqrt(1 - cosPhi.^2) < goGNSS.ELL_A_GPS;
                eclipsed(shadowCrossing) = 1;
            end
            %noon/midnight maneuvers affect all satellites
            noonMidnightTurn = acos(abs(cosPhi)) < t;
            eclipsed(noonMidnightTurn) = 3;
            
            
            
        end
        function importSP3Struct(this, sp3)
            this.time = sp3.time;
            this.coord =permute(sp3.coord,[3 2 1]);
            this.clock = sp3.clock',
            this.prn = sp3.prn;
            this.sys = sp3.sys;
            this.time_hr = sp3.time_hr;
            this.clock_hr = sp3.clock_hr;
            this.coord_rate = sp3.coord_rate;
            this.clock_rate = sp3.clock_rate;
            this.t_sun = sp3.t_sun;
            this.X_sun = sp3.X_sun';
            this.X_moon = sp3.X_moon';
            this.ERP = sp3.ERP;
            this.DCB = sp3.DCB;
            this.antenna_PCO = sp3.antPCO;
            this.start_time_idx = find(this.time == 1);
            
            
        end
        function load_antenna_PCV(this, filename_pco)
            antmod_S = this.cc.getAntennaId();
            this.antenna_PCV = read_antenna_PCV(filename_pco, antmod_S, this.time_ref_coord.getMatlabTime());
            this.antenna_PCO = zeros(1,this.cc.getNumSat(),3);
            %this.satType = cell(1,size(this.antenna_PCV,2));
            if isempty(this.avail)
                this.avail = zeros(size(this.antenna_PCV,2),1)
            end
            for sat = 1 : size(this.antenna_PCV,2)
                if (this.antenna_PCV(sat).n_frequency ~= 0)
                    this.antenna_PCO(:,sat,:) = this.antenna_PCV(sat).offset(:,:,1);
                    %this.satType{1,sat} = this.antenna_PCV(sat).sat_type;
                else
                    this.avail(sat) = 0;
                end
            end
        end
        function [stidecorr] = solid_earth_tide_correction(this, time, XR, XS, p_rate, phiC, lam)
            
            % SYNTAX:
            %   [stidecorr] = solid_earth_tide_correction(time, XR, XS,p_rate, phiC, lam);
            %
            % INPUT:
            %   time = GPS time
            %   XR   = receiver position  (X,Y,Z)
            %   XS   = satellite position (X,Y,Z)
            %   p_rate   = processing interval [s]
            %   phiC = receiver geocentric latitude (rad)
            %   lam  = receiver longitude (rad)
            %
            % OUTPUT:
            %   stidecorr = solid Earth tide correction terms (along the satellite-receiver line-of-sight)
            %
            % DESCRIPTION:
            %   Computation of the solid Earth tide displacement terms.
            
            
            if (nargin < 6)
                [~, lam, ~, phiC] = cart2geod(XR(1,1), XR(2,1), XR(3,1));
            end
            %north (b) and radial (c) local unit vectors
            b = [-sin(phiC)*cos(lam) -sin(phiC)*sin(lam) cos(phiC)];
            c = [+cos(phiC)*cos(lam) +cos(phiC)*sin(lam) sin(phiC)];
            
            %Sun and Moon position
            t_sun  = this.t_sun;
            X_sun  = this.X_sun;
            X_moon = this.X_moon;
            %[~, q] = min(abs(t_sun - time));
            % speed improvement of the above line
            % supposing t_sun regularly sampled
            q = round((time - t_sun(1)) / p_rate) + 1;
            X_sun  = X_sun(q,:);
            X_moon = X_moon(q,:);
            
            %receiver geocentric position
            XR_n = norm(XR);
            XR_u = XR / XR_n;
            
            %sun geocentric position
            X_sun_n = norm(X_sun);
            X_sun_u = X_sun / X_sun_n;
            
            %moon geocentric position
            X_moon_n = norm(X_moon);
            X_moon_u = X_moon / X_moon_n;
            
            %latitude dependence
            p = (3*sin(phiC)^2-1)/2;
            
            %gravitational parameters
            GE = goGNSS.GM_GAL; %Earth
            GS = GE*332946.0; %Sun
            GM = GE*0.01230002; %Moon
            
            %Earth equatorial radius
            R = 6378136.6;
            
            %nominal degree 2 Love number
            H2 = 0.6078 - 0.0006*p;
            %nominal degree 2 Shida number
            L2 = 0.0847 + 0.0002*p;
            
            %solid Earth tide displacement (degree 2)
            Vsun  = sum(conj(X_sun_u) .* XR_u);
            Vmoon = sum(conj(X_moon_u) .* XR_u);
            r_sun2  = (GS*R^4)/(GE*X_sun_n^3) *(H2*XR_u*(1.5*Vsun^2  - 0.5) + 3*L2*Vsun *(X_sun_u  - Vsun *XR_u));
            r_moon2 = (GM*R^4)/(GE*X_moon_n^3)*(H2*XR_u*(1.5*Vmoon^2 - 0.5) + 3*L2*Vmoon*(X_moon_u - Vmoon*XR_u));
            r = r_sun2 + r_moon2;
            
            %nominal degree 3 Love number
            H3 = 0.292;
            %nominal degree 3 Shida number
            L3 = 0.015;
            
            %solid Earth tide displacement (degree 3)
            r_sun3  = (GS*R^5)/(GE*X_sun_n^4) *(H3*XR_u*(2.5*Vsun^3  - 1.5*Vsun)  +   L3*(7.5*Vsun^2  - 1.5)*(X_sun_u  - Vsun *XR_u));
            r_moon3 = (GM*R^5)/(GE*X_moon_n^4)*(H3*XR_u*(2.5*Vmoon^3 - 1.5*Vmoon) +   L3*(7.5*Vmoon^2 - 1.5)*(X_moon_u - Vmoon*XR_u));
            r = r + r_sun3 + r_moon3;
            
            %from "conventional tide free" to "mean tide"
            radial = (-0.1206 + 0.0001*p)*p;
            north  = (-0.0252 + 0.0001*p)*sin(2*phiC);
            r = r + radial*c + north*b;
            
            %displacement along the receiver-satellite line-of-sight
            stidecorr = zeros(size(XS,1),1);
            for s = 1 : size(XS,1)
                LOS  = XR - XS(s,:);
                LOSu = LOS / norm(LOS);
                %stidecorr(s,1) = dot(r,LOSu);
                stidecorr(s,1) = sum(conj(r).*LOSu);
            end
        end
        function [poletidecorr] = pole_tide_correction(this, time, XR, XS, phiC, lam)
            
            % SYNTAX:
            %   [poletidecorr] = pole_tide_correction(time, XR, XS, SP3, phiC, lam);
            %
            % INPUT:
            %   time = GPS time
            %   XR   = receiver position  (X,Y,Z)
            %   XS   = satellite position (X,Y,Z)
            %   phiC = receiver geocentric latitude (rad)
            %   lam  = receiver longitude (rad)
            %
            % OUTPUT:
            %   poletidecorr = pole tide correction terms (along the satellite-receiver line-of-sight)
            %
            % DESCRIPTION:
            %   Computation of the pole tide displacement terms.
            if (nargin < 5)
                [~, lam, ~, phiC] = cart2geod(XR(1,1), XR(2,1), XR(3,1));
            end
            
            poletidecorr = zeros(size(XS,1),1);
            
            %interpolate the pole displacements
            if (~isempty(this.ERP))
                if (length(this.ERP.t) > 1)
                    m1 = interp1(this.ERP.t, this.ERP.m1, time, 'linear', 'extrap');
                    m2 = interp1(this.ERP.t, this.ERP.m2, time, 'linear', 'extrap');
                else
                    m1 = this.ERP.m1;
                    m2 = this.ERP.m2;
                end
                
                deltaR   = -33*sin(2*phiC)*(m1*cos(lam) + m2*sin(lam))*1e-3;
                deltaLam =  9* cos(  phiC)*(m1*sin(lam) - m2*cos(lam))*1e-3;
                deltaPhi = -9* cos(2*phiC)*(m1*cos(lam) + m2*sin(lam))*1e-3;
                
                corrENU(1,1) = deltaLam; %east
                corrENU(2,1) = deltaPhi; %north
                corrENU(3,1) = deltaR;   %up
                
                %displacement along the receiver-satellite line-of-sight
                XRcorr = local2globalPos(corrENU, XR);
                corrXYZ = XRcorr - XR;
                for s = 1 : size(XS,1)
                    LOS  = XR - XS(s,:)';
                    LOSu = LOS / norm(LOS);
                    poletidecorr(s,1) = -dot(corrXYZ,LOSu);
                end
            end
            
        end
        
        
        function writeSP3(this, f_name, prec)
            % SYNTAX:
            %   eph_tab.writeSP3(f_name, prec)
            %
            % INPUT:
            %   f_name       = file name of the sp3 file to be written
            %   prec        = precision (cm) of satellite orbit for all
            %   satellites (default 100)
            %
            %
            % DESCRIPTION:
            %   Write the current satellite postions and clocks bias into a sp3
            %   file
            if nargin<3
                prec=99;
            end
            %%% check if clock rate and coord rate are compatible
            rate_ratio=this.coord_rate/this.clock_rate;
            if abs(rate_ratio-round(rate_ratio)) > 0.00000001
                this.log.addWarning(sprintf('Incompatible coord rate (%s) and clock rate (%s) , sp3 not produced',this.coord_rate,this.clock_rate))
                return
            end
            %%% check if sun and moon positions ahve been computed
            if isempty(this.X_sun) || this.X_sun(1,1)==0
                this.sun_moon_pos();
            end
            %%% compute center of mass position (X_sat - PCO)
            switch_back = false;
            if this.coord_type == 1
                this.toCOM();
                switch_back = true;
            end
            %%% write to file
            rate_ratio = round(rate_ratio);
            fid=fopen(f_name,'w');
            this.writeHeader(fid, prec);
            
            for i=1:length(this.coord)
                this.writeEpoch(fid,[squeeze(this.coord(i,:,:)/1000) this.clock((i-1)/rate_ratio+1,:)'*1000000],i); %% convert coord in km and clock in microsecodns
            end
            fprintf(fid,'EOF\n');
            fclose(fid);
            if switch_back
                this.toAPC();
            end
            
            
            
        end
        function writeHeader(this, fid, prec)
            
            if nargin<3
                %%% unknown precision
                prec=99;
            end
            %prec = num2str(prec);
            time=this.time_ref_coord.getCopy();
            str_time = time.toString();
            year = str2num(str_time(1:4));
            month = str2num(str_time(6:7));
            day = str2num(str_time(9:10));
            hour = str2num(str_time(12:13));
            minute = str2num(str_time(15:16));
            second = str2num(str_time(18:27));
            week = time.getGpsWeek();
            sow = time.getGpsTime()-week*7*86400;
            mjd = jd2mjd(cal2jd(year,month,day));
            d_frac = hour/24+minute/24*60+second/86400;
            step = this.coord_rate;
            num_epoch = length(this.time_ref_coord);
            cc = this.cc;
            fprintf(fid,'#cP%4i %2i %2i %2i %2i %11.8f %7i d+D   IGS14 CNV GReD\n',year,month,day,hour,minute,second,num_epoch);
            fprintf(fid,'## %4i %15.8f %14.8f %5i %15.13f\n',week,sow,step,mjd,d_frac);
            
            sats = [];
            pre = [];
            ids = cc.prn;
            for i = 1:length(ids)
                sats=[sats, strrep(sprintf('%s%2i', cc.system(i), ids(i)), ' ', '0')];
                pre=[pre, sprintf('%3i', prec)];
            end
            n_row=ceil(length(sats)/51);
            rows=cell(5,1);
            rows(:)={repmat('  0',1,17)};
            pres=cell(5,1);
            pres(:)={repmat('  0',1,17)};
            for i =1:n_row
                rows{i}=sats((i-1)*51+1:min(length(sats),i*51));
                pres{i}=pre((i-1)*51+1:min(length(pre),i*51));
            end
            last_row_length=length((i-1)*51+1:length(sats));
            rows{n_row}=[rows{n_row} repmat('  0',1,(51-last_row_length)/3)];
            pres{n_row}=[pres{n_row} repmat('  0',1,(51-last_row_length)/3)];
            
            fprintf(fid,'+   %2i   %s\n',sum(cc.n_sat),rows{1});
            for i=2:length(rows)
                fprintf(fid,'+        %s\n',rows{i});
            end
            for i=1:length(rows)
                fprintf(fid,'++       %s\n',pres{i});
            end
            fprintf(fid,'%%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n');
            fprintf(fid,'%%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n');
            fprintf(fid,'%%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n');
            fprintf(fid,'%%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n');
            fprintf(fid,'%%i    0    0    0    0      0      0      0      0         0\n');
            fprintf(fid,'%%i    0    0    0    0      0      0      0      0         0\n');
            fprintf(fid,'/* Produced using goGPS                                     \n');
            fprintf(fid,'/*                 Non                                      \n');
            fprintf(fid,'/*                     Optional                             \n');
            fprintf(fid,'/*                              Lines                       \n');
        end
        function writeEpoch(this,fid,XYZT,epoch)
            t=this.time_ref_coord.getCopy();
            t.addIntSeconds((epoch)*900);
            cc=this.cc;
            str_time=t.toString();
            year=str2num(str_time(1:4));
            month=str2num(str_time(6:7));
            day=str2num(str_time(9:10));
            hour=str2num(str_time(12:13));
            minute=str2num(str_time(15:16));
            second=str2num(str_time(18:27));
            fprintf(fid,'*  %4i %2i %2i %2i %2i %11.8f\n',year,month,day,hour,minute,second);
            for i=1:size(XYZT,1)
                fprintf(fid,'P%s%14.6f%14.6f%14.6f%14.6f\n',strrep(sprintf('%s%2i', cc.system(i), cc.prn(i)), ' ', '0'),XYZT(i,1),XYZT(i,2),XYZT(i,3),XYZT(i,4));
            end
            
        end
        
        
        
    end
    
    
    methods (Static)
    end
end
