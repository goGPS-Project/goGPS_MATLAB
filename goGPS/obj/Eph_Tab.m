classdef Eph_Tab < handle
    
    
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
        time
        coord
        clock
        prn
        sys
        time_hr
        clock_hr
        coord_rate = 900;
        clock_rate = 900;
        iono
        t_sun
        X_sun
        t_sun_rate
        X_moon
        ERP
        antPCO
        satType
        avail
        
    end
    properties (Access = private)
        cc
        logger
    end
    methods
        function this = Eph_Tab(cc)
            if nargin == 0
                this.cc = Go_State.getCurrentSettings().getConstellationCollector();
            else
                this.cc = cc;
            end
            this.logger = Logger.getInstance();
            this.sys = this.cc.system;
            this.prn = this.cc.prn;
            this.antPCO =zeros(1,3,this.cc.getNumSat());
            
        end
        function importEph(this, eph, t_st, t_end, sat, step)
            % SYNTAX:
            %   eph_tab.importEph(eph, t_st, t_end, sat, step)
            %
            % INPUT:
            %   eph         = ephemerids matrix
            %   eph         = ephemerids matrix
            %   sat         = available satellite indexes
            %
            % OUTPUT:
            %   XS      = satellite position at time in ECEF(time_rx) (X,Y,Z)
            %   VS      = satellite velocity at time in ECEF(time_tx) (X,Y,Z)
            %   dtS     = satellite clock error (vector)
            %
            % DESCRIPTION:
            if nargin < 6
                step = 900;
            end
            this.coord_rate = step;
            this.clock_rate = step;
            if nargin < 4 || t_end.isempty()
                t_end = t_st;
            end
            times = t_st.getGpsTime : step : t_end.getGpsTime;
            this.time=times;
            sat = intersect(sat,unique(eph(30,:))); %% keep only satellite also present in eph
            i = 0;
            
            this.coord = zeros(3, this.cc.getNumSat, length(times));
            this.clock = zeros (this.cc.getNumSat, length(times));
            t_dist_exced=false;
            for t = times
                i=i+1;
                [this.coord(:,:,i), ~, this.clock(:,i), t_d_e]=this.satellitePositions(t, sat, eph); %%%% consider loss of precision problem
                t_dist_exced = t_dist_exced || t_d_e;
            end
            if t_dist_exced
                this.logger.addError(sprintf('One of the time bonds (%s , %s)\ntoo far from valid ephemerids ',t_st.toString(0),t_end.toString(0)))
            end
        end
        function importBrdc(this,f_name, t_st, t_end, sat, step)
            if nargin < 6
                step = 900;
            end
            if nargin < 5
                sat= this.cc.index;
            end
            if nargin < 4 || t_end.isempty()
                t_end = t_st;
            end
            [eph, this.iono, flag_return] = load_RINEX_nav(f_name,this.cc,0,0);
            if isempty(eph)
                
            end
            this.importEph(eph,t_st,t_end,sat,step);
            
        end
        function importIono(this,f_name)
            %%% to be implemeted
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
            
            
            XS=XS';
        end
        function [pos_S, vel_S] = interpolate_SP3_coord(this,time, sat)
            % SYNTAX:
            %   [pos_S, vel_S] = interpolate_SP3_coord(time, SP3, sat);
            %
            % INPUT:
            %   time       = interpolation time (GPS time, continuous since 6-1-1980)
            %   sat        = satellite PRN
            %   p_rate     = processing interval [s]
            %
            % OUTPUT:
            %   pos_S = interpolated satellite coordinates
            %   vel_S = satellite velocity
            %
            % DESCRIPTION:
            %   SP3 (precise ephemeris) coordinates 1-second interpolation by Lagrange
            %   polynomials. Satellite velocity computation. Relativistic correction.
            SP3_time  = this.time;
            SP3_coord = squeeze(this.coord(:, sat, :));
            antPCO    = this.antPCO(:, :, sat);
            
            %number of seconds in a quarter of an hour
            interval = this.coord_rate;
            
            %find the SP3 epoch closest to the interpolation time
            %[~, p] = min(abs(SP3_time - time));
            % speed improvement of the above line
            % supposing SP3_time regularly sampled
            p = round((time - SP3_time(1)) / interval) + 1;
            
            b = SP3_time(p) - time;
            
            pos_S = zeros(3,1);
            
            %Lagrange interpolation
            %degree of interpolation polynomial (Lagrange)
            % n = 10;
            %u = (n/2+1) + (- b + (-1:1))/interval;
            u = 6 + (- b + (-1:1))/interval;    % using 6 since n = 10;
            %X_sat = fastLI(SP3_coord(:,p+(-n/2:n/2)), u);
            X_sat = fastLI(SP3_coord(:,p + (-5:5)), u);
            
            %apply satellite antenna phase center correction
            [i, j, k] = this.satellite_fixed_frame(time,X_sat(:,2));
            X_sat(:,2) = X_sat(:,2) + [i j k]*antPCO';
            
            pos_S(1,1) = X_sat(1,2);
            pos_S(2,1) = X_sat(2,2);
            pos_S(3,1) = X_sat(3,2);
            
            %compute velocity
            
            vel_S = (X_sat(:,3) - X_sat(:,1)) / 2;
        end
        function [i, j, k] = satellite_fixed_frame(this,time,X_sat)
            
            % SYNTAX:
            %   [i, j, k] = satellite_fixed_frame(time);
            %
            % INPUT:
            %   time     = GPS time
            %
            % OUTPUT:
            %   i = unit vector that completes the right-handed system
            %   j = resulting unit vector of the cross product of k vector with the unit vector from the satellite to Sun
            %   k = unit vector pointing from the Satellite Mass Centre (MC) to the Earth's centre
            %
            % DESCRIPTION:
            %   Computation of the unit vectors defining the satellite-fixed frame.
            
            
            t_sun = this.t_sun;
            X_sun = this.X_sun;
            
            %[~, q] = min(abs(t_sun - time));
            % speed improvement of the above line
            % supposing t_sun regularly sampled
            q = round((time - t_sun(1)) / this.coord_rate) + 1;
            X_sun = X_sun(:,q);
            e = (X_sun - X_sat) / norm(X_sun - X_sat);
            k = -X_sat/norm(X_sat);
            %j = cross(k,e);
            j = [k(2).*e(3)-k(3).*e(2);
                k(3).*e(1)-k(1).*e(3);
                k(1).*e(2)-k(2).*e(1)];
            %i = cross(j,k);
            i = [j(2).*k(3)-j(3).*k(2);
                j(3).*k(1)-j(1).*k(3);
                j(1).*k(2)-j(2).*k(1)];
            
            j = j / norm(j);
            i = i / norm(i);
            
            
        end
        function [dt_S_SP3] = interpolate_SP3_clock(this,time, sat)
            
            % SYNTAX:
            %   [dt_S_SP3] = interpolate_SP3_clock(time, SP3, sat);
            %
            % INPUT:
            %   time  = interpolation timespan (GPS time, continuous since 6-1-1980)
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
            if (isempty(this.clock_hr))
                SP3_time = this.time;
                SP3_clock = this.clock;
            else
                SP3_time = this.time_hr;
                SP3_clock = this.clock_hr;
            end
            
            interval = this.clock_rate;
            
            %find the SP3 epoch closest to the interpolation time
            %[~, p] = min(abs(SP3_time - time));
            % speed improvement of the above line
            % supposing SP3_time regularly sampled
            p = round((time - SP3_time(1)) / interval) + 1;
            
            b = SP3_time(p) - time;
            
            %extract the SP3 clocks
            if (b>0)
                SP3_c = [SP3_clock(sat,p-1) SP3_clock(sat,p)];
                u = 1 - b/interval;
            else
                SP3_c = [SP3_clock(sat,p) SP3_clock(sat,p+1)];
                u = -b/interval;
            end
            
            dt_S_SP3  = NaN*ones(size(SP3_c,1),length(time));
            idx=(sum(SP3_c~=0,2) == 2 .* ~any(SP3_c >= 0.999,2))>0;
            dt_S_SP3(idx)=(1-u)*SP3_c(idx,1) + u*SP3_c(idx,2);
            
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
        function sun_moon_pos(this)
            
            global iephem km ephname inutate psicor epscor ob2000
            time = GPS_Time((this.time(1))/86400+GPS_Time.GPS_ZERO);
            
            
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
            
            sun_ECEF = zeros(3,length(this.time));
            moon_ECEF = zeros(3,length(this.time));
            this.t_sun = zeros(1,length(this.time));
            for e = 1 : length(this.time)
                [year , month ,day,~,~,~]= time.getCalEpoch();
                %UTC to TDB
                jdutc = julian(month, day, year);
                jdtdb = utc2tdb(jdutc);
                this.t_sun(e) = (datenum( year,month, day) - GPS_Time.GPS_ZERO)*86400;
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
                sun_ECEF(:,e) = celter(tjdh, tjdl, xp, yp, sun_ECI);
                
                if (nargout > 1)
                    %compute the Moon position (ICRS coordinates)
                    rrd = jplephem(jdtdb, moon_id, earth_id);
                    moon_ECI = rrd(1:3);
                    moon_ECI = tmatrix*moon_ECI;
                    
                    %Moon ICRS coordinates to ITRS coordinates
                    deltat = getdt;
                    jdut1 = jdutc - deltat;
                    tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
                    moon_ECEF(:,e) = celter(tjdh, tjdl, xp, yp, moon_ECI);
                end
                time.addSeconds(this.coord_rate);
            end
            
            this.X_sun  = sun_ECEF*1e3;
            this.X_moon = moon_ECEF*1e3;
            
            %this.t_sun_rate =
        end
        function importSP3(this, sp3)
            this.time = sp3.time;
            this.coord =sp3.coord;
            this.clock = sp3.clock,
            this.prn = sp3.prn;
            this.sys = sp3.sys;
            this.time_hr = sp3.time_hr;
            this.clock_hr = sp3.clock_hr;
            this.coord_rate = sp3.coord_rate;
            this.clock_rate = sp3.clock_rate;
            this.t_sun = sp3.t_sun;
            this.X_sun = sp3.X_sun;
            this.X_moon = sp3.X_moon;
            this.ERP = sp3.ERP;
            this.antPCO = sp3.antPCO;
            
            
        end
        function load_antenna_PCO(this, filename_pco)
            antmod_S = this.cc.getAntennaId();
            antenna_PCV_S = read_antenna_PCV(filename_pco, antmod_S, this.time(1));
            this.antPCO = zeros(1,3,size(antenna_PCV_S,2));
            this.satType = cell(1,size(antenna_PCV_S,2));
            if isempty(this.avail)
                this.avail=zeros(size(antenna_PCV_S,2),1)
            end
            for sat = 1 : size(antenna_PCV_S,2)
                if (antenna_PCV_S(sat).n_frequency ~= 0)
                    this.antPCO(:,:,sat) = antenna_PCV_S(sat).offset(:,:,1);
                    this.satType{1,sat} = antenna_PCV_S(sat).type;
                else
                    this.avail(sat) = 0;
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
                prec=100,
            end
            %%% check if clock rate and coord rate are compatible
            rate_ratio=this.coord_rate/this.clock_rate
            if abs(rate_ratio-round(rate_ratio)) > 0.00000001
                this.logger.addWarning(sprintf('Incompatible coord rate (%s) and clock rate (%s) , sp3 not produced',this.coord_rate,this.clock_rate))
                return
            end
            rate_ratio = round(rate_ratio);
            fid=fopen(f_name,'w');
            this.writeHeader(fid, prec);
            
            for i=1:length(this.time)
                this.writeEpoch(fid,[this.coord(:,:,i)'/1000 this.clock(:,(i-1)/rate_ratio+1)],i)
            end
            fclose(fid);
            
            
            
            
        end
        function writeHeader(this, fid, prec)
            
            if nargin<3
                prec=100,
            end
            prec = num2str(prec);
            time = GPS_Time((this.time(1))/86400+GPS_Time.GPS_ZERO);
            str_time = time.toString()
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
            num_epoch = length(this.time);
            cc = this.cc;
            fprintf(fid,'#cP%4i %2i %2i %2i %2i %11.8f %7i d+D   IGS14 CNV GReD\n',year,month,day,hour,minute,second,num_epoch);
            fprintf(fid,'## %4i %14.8f %14.8f %5i %15.13f\n',week,sow,step,mjd,d_frac);
            
            sats = [];
            pre = [];
            ids = cc.prn;
            for i = 1:length(ids)
                sats=[sats, strrep(sprintf('%s%2i', cc.system(i), ids(i)), ' ', '0')];
                pre=[pre, prec];
            end
            n_row=ceil(length(sats)/51);
            rows=cell(n_row,1);
            pres=cell(n_row,1);
            for i =1:n_row
                rows{i}=sats((i-1)*51+1:min(length(sats),i*51));
                pres{i}=pre((i-1)*51+1:min(length(pre),i*51));
            end
            rows{end}=[rows{end} repmat('  0',1,(51-length(rows{end}))/3)];
            pres{end}=[pres{end} repmat('  0',1,(51-length(pres{end}))/3)];
            fprintf(fid,'+   %2i   %s\n',cc.n_sat,rows{1});
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
            fprintf(fid,'/* GoGPS orbits and clocks, converted from broadcast        \n');
        end
        function writeEpoch(this,fid,XYZT,epoch)
            t=this.time(epoch);
            t=GPS_Time(t/86400+GPS_Time.GPS_ZERO);
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