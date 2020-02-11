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
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
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
            'first_lat',  [], ...    % first latitude
            'first_lon',  [], ...    % first longitude
            'd_lat',      [], ...    % lat spacing
            'd_lon',      [], ...    % lon_spacing
            'n_lat',      [], ...    % num lat
            'n_lon',      [], ...    % num lon
            'first_time', [], ...    % times [time] of the maps
            'first_time_double', [], ...    % times [time] of the maps [seconds from GPS zero]
            'dt',         [], ...    % time spacing
            'n_t',        [], ...    % num of epocvhs
            'height',     []  ...    % heigh of the layer
            )
        atm_load_nt = struct( ...
            'data_u',       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
            'data_e',       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
            'data_n',       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
            'first_lat',  [], ...    % first latitude
            'first_lon',  [], ...    % first longitude
            'd_lat',      [], ...    % lat spacing
            'd_lon',      [], ...    % lon_spacing
            'n_lat',      [], ...    % num lat
            'n_lon',      [], ...    % num lon
            'first_time', [], ...    % times [time] of the maps
            'first_time_double', [], ...    % times [time] of the maps [seconds from GPS zero]
            'dt',         [], ...    % time spacing
            'n_t',        [] ...   % num of epocvhs
            )
        atm_load_t = struct( ...
            'harmonics',       [], ...    % ionosphere single layer map [n_lat x _nlon n n_harmonics]
            'first_lat',  89.5, ...    % first latitude
            'first_lon',  0.5, ...    % first longitude
            'd_lat',      1, ...    % lat spacing
            'd_lon',      1, ...    % lon_spacing
            'n_lat',      180, ...    % num lat
            'n_lon',      360, ...    % num lon
            'first_time', 0, ...    % times [time] of the maps
            'first_time_double', 0, ...    % times [time] of the maps [seconds from GPS zero]
            'dt',         1, ...    % time spacing
            'n_t',        0 ...   % num of epocvhs
            )
        vmf_coeff = struct( ...
            'ell_height', [], ... %ellipsoidal height for the vmf values
            'ah',       [], ...    % alpha coefficient dry
            'aw',       [], ...    % alpha coefficent wet
            'zhd',       [], ...    % zhd
            'zwd',       [], ...    % zwd
            'first_lat',  [], ...    % first latitude
            'first_lon',  [], ...    % first longitude
            'd_lat',      [], ...    % lat spacing
            'd_lon',      [], ...    % lon_spacing
            'n_lat',      [], ...    % num lat
            'n_lon',      [], ...    % num lon
            'first_time', [], ...    % times [time] of the maps
            'first_time_double', [], ...    % times [time] of the maps [seconds from GPS zero]
            'dt',         [], ...    % time spacing
            'n_t',        [] ...   % num of epocvhs
            )
        
        emf     % current Earth geomagnetic field object
    end
    
    properties  (SetAccess = private, GetAccess = private)
        lat = 1e3;    % current values for lat
        lon = 1e3;    % current values for lon
        undu         % current values for undu
        
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
        
        
        log
        
        V_LIGHT = Core_Utils.V_LIGHT;
    end
    
    methods
        function this = Atmosphere()
            % Initialisation of the variables
            %
            % SYNTAX
            %   this = Atmosphere()
            this.state = Core.getCurrentSettings();
            this.log = Core.getLogger();
        end
    end
    
    methods
        function importIonex(this, file_name)
            % import IONEX file
            % if the files data overlaps or touch the present data in the object the new data are appended to the old one, otherwise the old data re-discarded
            %
            % SYNTAX
            %   importIonex(this, file_name)
            fid = fopen(file_name,'r');
            fnp = File_Name_Processor();
            if fid == -1
                this.log.addWarning(sprintf('File IONEX %s not found', file_name));
                return
            end
            this.log.addMessage(this.log.indent(sprintf('Opening file %s for reading', fnp.getFileName(file_name))));
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
                if instr(line,'END OF HEADER')
                    break
                elseif instr(line,'EPOCH OF FIRST MAP')
                    first_epoch = GPS_Time(sscanf(line(1:60),'%f %f %f %f %f %f')');
                elseif instr(line,'EPOCH OF LAST MAP')
                    last_epoch = GPS_Time(sscanf(line(1:60),'%f %f %f %f %f %f')');
                elseif instr(line,'INTERVAL')
                    interval = sscanf(line(1:60),'%f')';
                elseif instr(line,'HGT1 / HGT2 / DHGT')
                    height = sscanf(line(1:60),'%f %f %f')';
                elseif instr(line,'LAT1 / LAT2 / DLAT')
                    lats = sscanf(line(1:60),'%f %f %f')';
                elseif instr(line,'LON1 / LON2 / DLON')
                    lons = sscanf(line(1:60),'%f %f %f')';
                elseif instr(line,'EXPONENT')
                    exponent = sscanf(line(1:60),'%f')';
                end
            end
            %-------------------
            isempty_obj = isempty(this.ionex.data);
            if not(isempty_obj) && (first_epoch >= this.ionex.first_time) && (first_epoch - this.ionex.first_time) < this.ionex.d_t *  this.ionex.n_t
                %% file is contained in the data already
                this.log.addMessage(this.log.indent('File already present, skipping'));
            else
                if not(isempty_obj) && ((first_epoch - this.ionex.first_time) > (this.ionex.d_t *  this.ionex.n_t +3600*24) || (first_epoch - this.ionex.first_time) < (-3600*24))
                    %% file too far away emptying the object
                    this.log.addMessage(this.log.indent('File too far away, emptying old atmospheric pressure loading'));
                    this.clearAtmLoad();
                    isempty_obj = true;
                end
                lim(1:l,:) = [];
                txt(1:(lim(1,1)-1)) = [];
                lim(:,1:2) = lim(:,1:2) - lim(1,1) +1;
                n_lat = round((lats(2)-lats(1))/lats(3))+1;
                n_lon = round((lons(2)-lons(1))/lons(3))+1;
                n_line_1_lat = ceil(n_lon*5 / 80);
                
                nt = round((last_epoch - first_epoch) / interval);
                lines = repmat([false; false; repmat([false; true; false(n_line_1_lat-1,1) ],n_lat,1); false],nt,1);
                st_l  = lim(lines, 1);
                cols = [0:(n_lon*5+n_line_1_lat-2)];
                idx = repmat(cols,length(st_l),1) + repmat(st_l,1,length(cols));
                idx(:,81:81:length(cols))   = [];
                idx(:,366:end) = []; %% trial and error fix bug fix
                vals = txt(serialize(idx'));
                vals = reshape(vals,5,length(vals)/5);
                nums = sscanf(vals,'%f');
                data = permute(reshape(nums,n_lon,n_lat,nt),[2,1,3])*10^(exponent);
                if mod(lons(1), 360) == mod(lons(2),360)
                    data(:,end,:) = [];
                    n_lon = n_lon -1;
                end
                if isempty_obj
                    this.ionex.first_time = first_epoch;
                    this.ionex.first_time_double = first_epoch.getGpsTime();
                    this.ionex.d_t = interval;
                    this.ionex.n_t =  nt;
                    this.ionex.first_lat= lats(1);
                    this.ionex.d_lat = lats(3);
                    this.ionex.n_lat = n_lat ;
                    this.ionex.first_lon= lons(1);
                    this.ionex.d_lon = lons(3);
                    this.ionex.n_lon = n_lon;
                    this.ionex.height = height;
                    this.ionex.data = data;
                else
                    if first_epoch < this.ionex.first_time
                        this.ionex.first_time = first_epoch;
                        this.ionex.first_time_double = first_epoch.getGpsTime();
                        this.ionex.data = cat(3,data,this.ionex.data);
                    else
                        this.ionex.data = cat(3,this.ionex.data,data);
                    end
                    this.ionex.n_t = this.ionex.n_t + 24;
                end
            end
        end
        
        function initIonex(this, dsa, dso)
            % initialise Ionosphere map repository importing all the necessary IONEX files
            %
            % SYNTAX
            %   initIonex(this, dsa, dso)
            dso = dso.getCopy();
            dsa = dsa.getCopy();
            dso.addSeconds(6*3600);
            fname = this.state.getIonoFileName( dsa, dso);
            for i = 1 : length(fname)
                this.importIonex(fname{i});
            end
        end
        
        function initVMF(this, dsa, dso)
            % Initialize Vienna Mapping Function
            %
            % SYNTAX
            %   initVMF(this, dsa, dso)
            dso = dso.getCopy();
            dsa = dsa.getCopy();
            %dso.addSeconds(6*3600);
            fname = this.state.getVMFFileName( dsa, dso);
            % import the coefficeints files
            for i = 1 : length(fname)
                this.importVMFCoeffFile(fname{i});
            end
            h_fname = this.state.getVMFHeightFileName();
            fid = fopen(h_fname);
            fgetl(fid); %skip one line header
            formatSpec = [repmat([repmat(' %f',1,10) '\n'],1,14) repmat(' %f',1,5)];
            h_data = cell2mat(textscan(fid,formatSpec,91));
            h_data(:,end) = [];
            this.vmf_coeff.ell_height = h_data;
            fclose(fid);
        end
        
        function importTidalAtmLoadHarmonics(this)
            % importing Tidal Atm and loading Harmonics
            %
            % SYNTAX
            %   importTidalAtmLoadHarmonics(this)
            fname = this.state.getTAtmLoadFileName();
            data = importdata(fname);
            this.atm_load_t.harmonics = permute(reshape(data(:,3:end),360 ,180, 18),[2 1 3])/1e3;
            
        end
        
        function importAtmLoadCoeffFile(this, filename)
            % import data of atmospehric loading file
            %
            % SYNTAX
            %   importAtmLoadCoeffFile(this, filename)
            fid = fopen([filename],'r');
            if fid == -1
                this.log.addWarning(sprintf('Athmosphere: File %s not found', filename));
                return
            else
                this.log.addMessage(this.log.indent(sprintf('Loading  %s', File_Name_Processor.getFileName(filename))));
            end
            [~, file_name, ~] = fileparts(filename);
            year = str2num(file_name(1:4));
            month = str2num(file_name(5:6));
            day = str2num(file_name(7:8));
            hour = str2num(file_name(9:10));
            file_ref_ep = GPS_Time([year month day hour 0 0]);
            isempty_obj = isempty(this.atm_load_nt.data_u);
            if not(isempty_obj) && (file_ref_ep >= this.atm_load_nt.first_time) && (file_ref_ep - this.atm_load_nt.first_time) < this.atm_load_nt.dt *  this.atm_load_nt.n_t
                % file is contained in the data already
                this.log.addMessage(this.log.indent('File already present, skipping'));
                fclose(fid);
            else
                if not(isempty_obj) && ((file_ref_ep - this.atm_load_nt.first_time) > (this.atm_load_nt.dt *  this.atm_load_nt.n_t +3600*6) || (file_ref_ep - this.atm_load_nt.first_time) < (-3600*6))
                    % file too far away emptying the object
                    this.log.addMessage(this.log.indent('File too far away, emptying old atmospheric pressure loading'));
                    this.clearAtmLoad();
                    isempty_obj = true;
                end
                
                txt = fread(fid, '*char')';
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
                
                % searching for header
                eoh = 0;
                while txt(lim(eoh + 1)) == '!'
                    eoh = eoh + 1;
                end
                
                % parsing data
                data = sscanf(txt(lim(eoh + 1, 1) : lim(end, 2)), '%f ');
                data_u = reshape(data(3 : 5 : end), 360, 180)';
                data_e = reshape(data(4 : 5 : end), 360, 180)';
                data_n = reshape(data(5 : 5 : end), 360, 180)';
                if isempty_obj
                    this.atm_load_nt.data_u = data_u;
                    this.atm_load_nt.data_e = data_e;
                    this.atm_load_nt.data_n = data_n;
                    this.atm_load_nt.first_lat = 89.5;
                    this.atm_load_nt.first_lon = 0.5;
                    this.atm_load_nt.d_lat = -1;
                    this.atm_load_nt.d_lon = 1;
                    this.atm_load_nt.n_lat = 180;
                    this.atm_load_nt.n_lon = 360;
                    this.atm_load_nt.first_time = file_ref_ep;
                    this.atm_load_nt.first_time_double = file_ref_ep.getGpsTime();
                    this.atm_load_nt.dt = 3600 * 6;
                    this.atm_load_nt.n_t = 1;
                else
                    if file_ref_ep < this.atm_load_nt.first_time
                        this.atm_load_nt.first_time = file_ref_ep;
                        this.atm_load_nt.first_time_double = file_ref_ep.getGpsTime();
                        this.atm_load_nt.data_u = cat(3, data_u, this.atm_load_nt.data_u);
                        this.atm_load_nt.data_e = cat(3, data_e, this.atm_load_nt.data_e);
                        this.atm_load_nt.data_n = cat(3, data_n, this.atm_load_nt.data_n);
                    else
                        this.atm_load_nt.data_u = cat(3, this.atm_load_nt.data_u, data_u);
                        this.atm_load_nt.data_e = cat(3, this.atm_load_nt.data_e, data_e);
                        this.atm_load_nt.data_n = cat(3, this.atm_load_nt.data_n, data_n);
                    end
                    this.atm_load_nt.n_t = this.atm_load_nt.n_t + 1;
                end
            end
        end
        
        function importVMFCoeffFile(this, ffile_name)
            % import data of atmospehric loading file
            %
            % SYNTAX
            %   importVMFCoeffFile(this, filename)
            
            fnp = File_Name_Processor;
            fid = fopen(ffile_name,'r');
            if fid == -1
                this.log.addWarning(sprintf('Atmosphere: File %s not found', ffile_name));
            else
                % Parsing title
                [~, file_name, ext ] = fileparts(ffile_name);
                year = str2double(file_name(6 : 9));
                month = str2double(file_name(10 : 11));
                day = str2double(file_name(12 : 13));
                hour = str2double(ext(3 : 4));
                
                file_ref_ep = GPS_Time([year month day hour 0 0]);
                isempty_obj = isempty(this.vmf_coeff.ah);
                
                if not(isempty_obj) && (file_ref_ep >= this.vmf_coeff.first_time) && (file_ref_ep - this.vmf_coeff.first_time) < this.vmf_coeff.dt *  this.vmf_coeff.n_t
                    % file is contained in the data already
                    this.log.addMessage(this.log.indent(sprintf('%s already present, skipping', fnp.getFileName(ffile_name))));
                    fclose(fid);
                else
                    if not(isempty_obj) && ((file_ref_ep - this.vmf_coeff.first_time) > (this.vmf_coeff.dt *  this.vmf_coeff.n_t +3600*6) || (file_ref_ep - this.vmf_coeff.first_time) < (-3600*6))
                        % file too far away emptying the object
                        this.log.addMessage(this.log.indent('File too far away, emptying old vmf coefficients'));
                        this.clearVMF();
                        isempty_obj = true;
                    end
                    this.log.addMessage(this.log.indent(sprintf('Opening file %s for reading', fnp.getFileName(ffile_name))));
                    
                    txt = fread(fid, '*char')';
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
                    
                    % searching for header
                    eoh = 0;
                    while txt(lim(eoh + 1)) == '!'
                        eoh = eoh + 1;
                    end
                    
                    % parsing data
                    data = sscanf(txt(lim(eoh + 1, 1) : lim(end,2)), '%f ');
                    ah =  reshape(data(3:6:end), 144, 91)';
                    aw =  reshape(data(4:6:end), 144, 91)';
                    zhd = reshape(data(5:6:end), 144, 91)';
                    zwd = reshape(data(6:6:end), 144, 91)';
                    
                    if isempty_obj
                        this.vmf_coeff.ah         = ah;
                        this.vmf_coeff.aw         = aw;
                        this.vmf_coeff.zhd        = zhd;
                        this.vmf_coeff.zwd        = zwd;
                        this.vmf_coeff.first_lat  = 90;
                        this.vmf_coeff.first_lon  = 0.0;
                        this.vmf_coeff.d_lat      = -2;
                        this.vmf_coeff.d_lon      = 2.5;
                        this.vmf_coeff.n_lat      = 91;
                        this.vmf_coeff.n_lon      = 144;
                        this.vmf_coeff.first_time = file_ref_ep;
                        this.vmf_coeff.first_time_double = file_ref_ep.getGpsTime();
                        this.vmf_coeff.dt         = 3600*6;
                        this.vmf_coeff.n_t        = 1;
                    else
                        if file_ref_ep < this.vmf_coeff.first_time
                            this.vmf_coeff.first_time = file_ref_ep;
                            this.vmf_coeff.first_time_double = file_ref_ep.getGpsTime();
                            this.vmf_coeff.ah  = cat(3, ah,  this.vmf_coeff.ah);
                            this.vmf_coeff.aw  = cat(3, aw,  this.vmf_coeff.aw);
                            this.vmf_coeff.zhd = cat(3, zhd, this.vmf_coeff.zhd);
                            this.vmf_coeff.zwd = cat(3, zwd, this.vmf_coeff.zwd);
                        else
                            this.vmf_coeff.ah  = cat(3, this.vmf_coeff.ah,  ah);
                            this.vmf_coeff.aw  = cat(3, this.vmf_coeff.aw,  aw);
                            this.vmf_coeff.zhd = cat(3, this.vmf_coeff.zhd, zhd);
                            this.vmf_coeff.zwd = cat(3, this.vmf_coeff.zwd, zwd);
                        end
                        this.vmf_coeff.n_t = this.vmf_coeff.n_t + 1;
                    end
                end
            end
        end
        
        function [ah, aw, zhd, zwd] = readVMF(file_name)
            fid = fopen(filename);
            if fid > 0
                txt = fread(fid, '*char')';
                txt(txt == 13) = []; % remove carriage return - I hate you Bill!
                
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
                
                % searching for header
                eoh = 0;
                while txt(lim(eoh + 1)) == '!'
                    eoh = eoh + 1;
                end
                
                data = sscanf(txt(lim(eoh + 1, 1) : lim(end,2)), '%f ');
                ah =  reshape(data(3:6:end), 144, 91)';
                aw =  reshape(data(4:6:end), 144, 91)';
                zhd = reshape(data(5:6:end), 144, 91)';
                zwd = reshape(data(6:6:end), 144, 91)';
                
                fclose(fid);
            else
                this.log.addWarning(sprintf('VMF file not found or invalid "%s"', file_name));
                ah = [];
                aw = [];
                zhd = [];
                zwd = [];
            end
        end
        
        function clearAtmLoad(this)
            % cleaning Atmosphere loading data
            %
            % SYNTAX
            %   clearAtmLoad(this)
            this.atm_load_nt = struct( ...
                'data_u',       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
                'data_e',       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
                'data_n',       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
                'first_lat',  [], ...    % first latitude
                'first_lon',  [], ...    % first longitude
                'd_lat',      [], ...    % lat spacing
                'd_lon',      [], ...    % lon_spacing
                'n_lat',      [], ...    % num lat
                'n_lon',      [], ...    % num lon
                'first_time', [], ...    % times [time] of the maps
                'first_time_double', [], ...    % times [time] of the maps [seconds from GPS zero]
                'dt',         [], ...    % time spacing
                'n_t',        [] ...    % num of epocvhs
                );
        end
        
        function clearIonex(this)
            % cleaning IONEX
            %
            % SYNTAX
            %   clearIonex(this)
            this.ionex = struct( ...
                'data',       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
                'first_lat',  [], ...    % first latitude
                'first_lon',  [], ...    % first longitude
                'd_lat',      [], ...    % lat spacing
                'd_lon',      [], ...    % lon_spacing
                'n_lat',      [], ...    % num lat
                'n_lon',      [], ...    % num lon
                'first_time', [], ...    % times [time] of the maps
                'first_time_double', [], ...    % times [time] of the maps [seconds from GPS zero]
                'dt',         [], ...    % time spacing
                'n_t',        [], ...    % num of epocvhs
                'height',     []  ...    % heigh of the layer
                );
        end
        
        function clearVMF(this)
            % cleaning Vienna Mapping Function data
            %
            % SYNTAX
            %   clearVMF(this)
            this.vmf_coeff = struct( ...
                'ah',       [], ...    % alpha coefficient dry
                'aw',       [], ...    % alpha coefficent wet
                'zhd',       [], ...    % zhd
                'zwd',       [], ...    % zwd
                'first_lat',  [], ...    % first latitude
                'first_lon',  [], ...    % first longitude
                'd_lat',      [], ...    % lat spacing
                'd_lon',      [], ...    % lon_spacing
                'n_lat',      [], ...    % num lat
                'n_lon',      [], ...    % num lon
                'first_time', [], ...    % times [time] of the maps
                'first_time_double', [], ...    % times [time] of the maps [seconds from GPS zero]
                'dt',         [], ...    % time spacing
                'n_t',        [] ...   % num of epocvhs
                );
        end
        
        function is_loaded = isAtmoLoadLoaded(this)
            % check if atmospheric  loading coefficient has been loaded
            %
            % SYNTAX
            %     is_loaded = this.isAtmoLoadLoaded()
            is_loaded = any(this.atm_load_nt.data_u);
        end
        
        function [corrxyz] = getNTAtmLoadCorr(this, lat, lon, h_ellips, time)
            % get non tidal atmospheric loading corrections
            %
            % SYNTAX
            %   [corrxyz] = getNTAtmLoadCorr(this, lat, lon, h_ellips, time)
            gps_time = time.getGpsTime();
            % interpolate the values from the downloaded grids
            up = Core_Utils.linInterpLatLonTime(this.atm_load_nt.data_u, this.atm_load_nt.first_lat, this.atm_load_nt.d_lat, this.atm_load_nt.first_lon, this.atm_load_nt.d_lon, this.atm_load_nt.first_time_double, this.atm_load_nt.dt, lat, lon,gps_time);
            east = Core_Utils.linInterpLatLonTime(this.atm_load_nt.data_e, this.atm_load_nt.first_lat, this.atm_load_nt.d_lat, this.atm_load_nt.first_lon, this.atm_load_nt.d_lon, this.atm_load_nt.first_time_double, this.atm_load_nt.dt, lat, lon,gps_time);
            north = Core_Utils.linInterpLatLonTime(this.atm_load_nt.data_n, this.atm_load_nt.first_lat, this.atm_load_nt.d_lat, this.atm_load_nt.first_lon, this.atm_load_nt.d_lon, this.atm_load_nt.first_time_double, this.atm_load_nt.dt, lat, lon,gps_time);
            % pass from north east up to cartesian coordinates
            [x,y,z] = geod2cart(lat, lon, h_ellips);
            [corrxyz] = local2globalVel([east north up]', repmat([x,y,z]',1,length(east)));
            corrxyz = corrxyz';
        end
        
        function [corrxyz] = getTAtmLoadCorr(this, lat, lon, h_ellips, time)
            % get non tidal atmospheric loading corrections
            %
            % SYNTAX
            %   [corrxyz] = getTAtmLoadCorr(this, lat, lon, h_ellips, time)
            
            %get the harmonics
            [harm_r, harm_e, harm_n] = getAtmLoadHarm(this, lat,lon);
            T = mod(time.getMJD, 1)*2*pi;
            % convert the displcaement from local to cartesian coordinates
            [x,y,z] = geod2cart(lat, lon, h_ellips);
            [harm_xyz] = local2globalVel([harm_n harm_e harm_r]', repmat([x,y,z]',1,length(harm_e)));
            harm_xyz = harm_xyz';
            % compute the displacents using diurnal semidiurnal and terdiurnal componenent
            corrx = harm_xyz(1,1)*sin(T) + harm_xyz(2,1)*cos(T) + harm_xyz(3,1)*sin(2*T) + harm_xyz(4,1)*cos(2*T) + harm_xyz(5,1)*sin(3*T) + harm_xyz(6,1)*cos(3*T);
            corry = harm_xyz(1,2)*sin(T) + harm_xyz(2,2)*cos(T) + harm_xyz(3,2)*sin(2*T) + harm_xyz(4,2)*cos(2*T) + harm_xyz(5,2)*sin(3*T) + harm_xyz(6,2)*cos(3*T);
            corrz = harm_xyz(1,3)*sin(T) + harm_xyz(2,3)*cos(T) + harm_xyz(3,3)*sin(2*T) + harm_xyz(4,3)*cos(2*T) + harm_xyz(5,3)*sin(3*T) + harm_xyz(6,3)*cos(3*T);
            corrxyz = [corrx corry corrz];
        end
        
        function [corrxyz] = getAtmLoadCorr(this, lat, lon, h_ellips, time)
            % get atmospheric loading corrections
            %
            % SYNTAX
            %   [corrxyz] = getAtmLoadCorr(this, lat, lon, h_ellips, time)
            [corrxyz_nt] = getNTAtmLoadCorr(this, lat, lon, h_ellips, time);
            [corrxyz_t] = getTAtmLoadCorr(this, lat, lon, h_ellips, time);
            corrxyz = corrxyz_nt + corrxyz_t;
        end
        
        function [harm_r, harm_e, harm_n] = getAtmLoadHarm(this, lat,lon)
            if isempty(this.atm_load_t.harmonics)
                this.importTidalAtmLoadHarmonics();
            end
            % Interpolate the harmonics at the receiver position
            % sine diurnal radial componenet
            S1SR = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,1), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine diurnal radial componenet
            S1CR = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,2), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % sine semi-diurnal radial componenet
            S2SR = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,3), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine semi-diurnal radial componenet
            S2CR = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,4), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % sine ter-diurnal radial componenet
            S3SR = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,5), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine ter-diurnal radial componenet
            S3CR = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,6), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % sine diurnal east componenet
            S1SE = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,7), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine diurnal east componenet
            S1CE = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,8), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % sine semi-diurnal east componenet
            S2SE = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,9), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine semi-diurnal east componenet
            S2CE = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,10), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % sine ter-diurnal east componenet
            S3SE = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,11), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine ter-diurnal east componenet
            S3CE = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,12), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % sine diurnal north componenet
            S1SN = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,13), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine diurnal north componenet
            S1CN = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,14), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % sine semi-diurnal north componenet
            S2SN = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,15), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine semi-diurnal north componenet
            S2CN = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,16), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % sine ter-diurnal north componenet
            S3SN = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,17), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            % cosine ter-diurnal north componenet
            S3CN = Core_Utils.linInterpLatLonTime(this.atm_load_t.harmonics(:,:,18), this.atm_load_t.first_lat, this.atm_load_t.d_lat, this.atm_load_t.first_lon, this.atm_load_t.d_lon, this.atm_load_t.first_time_double, this.atm_load_t.dt, lat, lon,0);
            
            harm_r = [S1SR S1CR S2SR S2CR S3SR S3CR]';
            harm_e = [S1SE S1CE S2SE S2CE S3SE S3CE]';
            harm_n = [S1SN S1CN S2SN S2CN S3SN S3CN]';
        end
        
        function [ah, aw] = interpolateAlpha(this, gps_time, lat, lon)
            % interpolate a (continued fraction form) coefficent for both wet idrostatit part
            %
            % SYNTAX
            %   [ah, aw] = interpolateAlpha(this, gps_time, lat, lon)
            ah = Core_Utils.linInterpLatLonTime(this.vmf_coeff.ah, this.vmf_coeff.first_lat, this.vmf_coeff.d_lat, this.vmf_coeff.first_lon, this.vmf_coeff.d_lon, this.vmf_coeff.first_time_double, this.vmf_coeff.dt, lat, lon,gps_time);
            aw = Core_Utils.linInterpLatLonTime(this.vmf_coeff.aw, this.vmf_coeff.first_lat, this.vmf_coeff.d_lat, this.vmf_coeff.first_lon, this.vmf_coeff.d_lon, this.vmf_coeff.first_time_double, this.vmf_coeff.dt, lat, lon,gps_time);
        end
        
        function [ it, st, ilons, ilone, slon, ilat, slat] = getVMFIndex(this, gps_time, lat, lon)
            % get the index of the elements of the vmf to be used in the trilinear interpolation
            %
            % SYNTAX
            %    [ it, st, ilons, ilone, slon, ilat, slat] = this.tgetVMFIndex(gps_time, lat, lon)
            [ it, st, ilons, ilone, slon, ilat, slat] = Core_Utils.getIntIdx(this.vmf_coeff.ah, this.vmf_coeff.first_lat, this.vmf_coeff.d_lat, this.vmf_coeff.first_lon, this.vmf_coeff.d_lon, this.vmf_coeff.first_time_double, this.vmf_coeff.dt, lat, lon,gps_time);
        end
        
        function [zhd] = interpolateZhd(this, gps_time, lat, lon)
            % interpolate Zenit hidrostatic delay
            %
            % SYNTAX
            %   [zhd] = interpolateZhd(this, gps_time, lat, lon)
            % REFERENCE:
            % [1] Kouba, Jan. "Implementation and testing of the gridded Vienna Mapping Function 1 (VMF1)." Journal of Geodesy 82.4-5 (2008): 193-205.
            zhd = Core_Utils.linInterpLatLonTime(this.vmf_coeff.zhd, this.vmf_coeff.first_lat, this.vmf_coeff.d_lat, this.vmf_coeff.first_lon, this.vmf_coeff.d_lon, this.vmf_coeff.first_time_double, this.vmf_coeff.dt, lat, lon,gps_time);
        end
        
        function [zhd] = getVmfZhd(this, gps_time, lat, lon, h_ell_sta, interp_first)
            % get vmf Zenit hidrostatic delay
            %
            % SYNTAX
            %   [zhd] = getVmfZhd(this, gps_time, lat, lon, h_ell_sta)
            % REFERENCE:
            % [1] Kouba, Jan. "Implementation and testing of the gridded Vienna Mapping Function 1 (VMF1)." Journal of Geodesy 82.4-5 (2008): 193-205.
            if nargin < 6
                interp_first = false;
            end
            if interp_first
                [zhd] = this.interpolateZhd( gps_time, lat, lon);
                h_ell_vmf = this.interpolateVMFElHeight(lat,lon); % get height of the station
                %             p_h_vmf = (zhd / 0.0022768) .* (1 - 0.00266 * cos(2*lat/180*pi) - 0.28 * 10^-6 * h_ell_vmf);   % formula 3 in [1]
                %             p_h_sta = p_h_vmf .* (1 - 0.0000226 .* (h_ell_sta - h_ell_vmf)).^5.225;   % formula 4 in [1]
                %             zhd     = 0.0022768 * p_h_sta / (1 - 0.00266 * cos(2 * lat/180*pi) - 0.28 * 10^-6 * h_ell_sta);   % formula 3  in [1] again
                
                par1 = 1 - 0.00266 * cos(2*lat/180*pi);
                par2 = 0.28 * 10^-6 ;
                zhd = zhd * (par1 -par2*h_ell_vmf)/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_ell_vmf)).^5.225; % simplification of formula 3 then formual 4 then formula 3 again
            else
                % get the index of the interpolating points
                [ it, st, ilons, ilone, slon, ilat, slat] = this.getVMFIndex(gps_time, lat, lon);
                % time t
                zhd_calc_1 = this.vmf_coeff.zhd([ilat ilat+1], [ilons ilone],it);
                % time t +1
                zhd_calc_2 = this.vmf_coeff.zhd([ilat ilat+1], [ilons ilone],it+1);
                h_calc =  this.vmf_coeff.ell_height([ilat ilat+1], [ilons ilone]);
                % compute the mapping function and the heoght correction for all the inteprolating points
                par1 = 1 - 0.00266 * cos(2*lat/180*pi);
                par2 = 0.28 * 10^-6 ;
                % time t
                zhd_calc_1(1,1,:) = zhd_calc_1(1,1,:) * (par1 -par2*h_calc(1,1))/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_calc(1,1))).^5.225;
                zhd_calc_1(1,2,:) = zhd_calc_1(1,2,:) * (par1 -par2*h_calc(1,2))/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_calc(1,2))).^5.225;
                zhd_calc_1(2,1,:) = zhd_calc_1(2,1,:) * (par1 -par2*h_calc(2,1))/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_calc(2,1))).^5.225;
                zhd_calc_1(2,2,:) = zhd_calc_1(2,2,:) * (par1 -par2*h_calc(2,2))/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_calc(2,2))).^5.225;
                % time t+1
                zhd_calc_2(1,1,:) = zhd_calc_2(1,1,:) * (par1 -par2*h_calc(1,1))/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_calc(1,1))).^5.225;
                zhd_calc_2(1,2,:) = zhd_calc_2(1,2,:) * (par1 -par2*h_calc(1,2))/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_calc(1,2))).^5.225;
                zhd_calc_2(2,1,:) = zhd_calc_2(2,1,:) * (par1 -par2*h_calc(2,1))/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_calc(2,1))).^5.225;
                zhd_calc_2(2,2,:) = zhd_calc_2(2,2,:) * (par1 -par2*h_calc(2,2))/(par1 -par2*h_ell_sta) * (1 - 0.0000226 .* (h_ell_sta - h_calc(2,2))).^5.225;
                
                
                %interpolate along lon
                valbu = zhd_calc_1(1,1,:).*(1-slon) + zhd_calc_1(1,2,:).*slon;
                valau = zhd_calc_2(1,1,:).*(1-slon) + zhd_calc_2(1,2,:).*slon;
                valbd = zhd_calc_1(2,1,:)*(1-slon) + zhd_calc_1(2,2,:).*slon;
                valad = zhd_calc_2(2,1,:)*(1-slon) + zhd_calc_2(2,2,:).*slon;
                
                %interpolate along lat
                valb = valbd.*(1-slat) + valbu.*slat;
                vala = valad.*(1-slat) + valau.*slat;
                
                %interpolate along time
                zhd = squeeze(valb).*(1-st) + squeeze(vala).*st;
                
            end
        end
        
        function [zwd] = interpolateZwd(this, gps_time, lat, lon)
            % interpolate Zenit wet delay
            %
            % SYNTAX
            %   [zwd] = interpolateZwd(this, gps_time, lat, lon)
            
            zwd = Core_Utils.linInterpLatLonTime(this.vmf_coeff.zwd, this.vmf_coeff.first_lat, this.vmf_coeff.d_lat, this.vmf_coeff.first_lon, this.vmf_coeff.d_lon, this.vmf_coeff.first_time_double, this.vmf_coeff.dt, lat, lon,gps_time);
            
        end
        
        function [zwd] = getVmfZwd(this, gps_time, lat, lon, h_ell_sta, interp_first)
            % get vmf Zenit wet delay
            %
            % SYNTAX
            %   [zhd] = getVmfZwd(this, gps_time, lat, lon, h_ell_sta)
            % REFERENCE:
            % [1] Kouba, Jan. "Implementation and testing of the gridded Vienna Mapping Function 1 (VMF1)." Journal of Geodesy 82.4-5 (2008): 193-205.
            if nargin < 6
                interp_first = false;
            end
            if interp_first
                [zwd] = this.interpolateZwd( gps_time, lat, lon);
                h_ell_vmf = this.interpolateVMFElHeight(lat,lon); % get height of the station
                zwd     = zwd * exp(-(h_ell_sta - h_ell_vmf)/2000);   % formula 5  in [1]
            else
                % get the index of the interpolating points
                [ it, st, ilons, ilone, slon, ilat, slat] = this.getVMFIndex(gps_time, lat, lon);
                zwd_calc_1 = this.vmf_coeff.zwd([ilat ilat+1], [ilons ilone],it);
                zwd_calc_2 = this.vmf_coeff.zwd([ilat ilat+1], [ilons ilone],it+1);
                h_calc =  this.vmf_coeff.ell_height([ilat ilat+1], [ilons ilone]);
                % compute the mapping function and the heoght correction for all the inteprolating points
                zwd_calc_1(1,1,:) = zwd_calc_1(1,1,:) * exp(-(h_ell_sta - h_calc(1,1))/2000);
                zwd_calc_1(1,2,:) = zwd_calc_1(1,2,:) * exp(-(h_ell_sta - h_calc(1,2))/2000);
                zwd_calc_1(2,1,:) = zwd_calc_1(2,1,:) * exp(-(h_ell_sta - h_calc(2,1))/2000);
                zwd_calc_1(2,2,:) = zwd_calc_1(2,2,:) * exp(-(h_ell_sta - h_calc(2,2))/2000);
                
                zwd_calc_2(1,1,:) = zwd_calc_2(1,1,:) * exp(-(h_ell_sta - h_calc(1,1))/2000);
                zwd_calc_2(1,2,:) = zwd_calc_2(1,2,:) * exp(-(h_ell_sta - h_calc(1,2))/2000);
                zwd_calc_2(2,1,:) = zwd_calc_2(2,1,:) * exp(-(h_ell_sta - h_calc(2,1))/2000);
                zwd_calc_2(2,2,:) = zwd_calc_2(2,2,:) * exp(-(h_ell_sta - h_calc(2,2))/2000);
                
                
                %interpolate along lon
                valbu = zwd_calc_1(1,1,:).*(1-slon) + zwd_calc_1(1,2,:).*slon;
                valau = zwd_calc_2(1,1,:).*(1-slon) + zwd_calc_2(1,2,:).*slon;
                valbd = zwd_calc_1(2,1,:)*(1-slon) + zwd_calc_1(2,2,:).*slon;
                valad = zwd_calc_2(2,1,:)*(1-slon) + zwd_calc_2(2,2,:).*slon;
                
                %interpolate along lat
                valb = valbd.*(1-slat) + valbu.*slat;
                vala = valad.*(1-slat) + valau.*slat;
                
                %interpolate along time
                zwd =  squeeze(valb).*(1-st) + squeeze(vala).*st;
                
                
                
            end
        end
        
        function [height] = interpolateVMFElHeight(this, lat, lon)
            % interpolate Zenit wet delay
            %
            % SYNTAX
            %   [zwd] = interpolateZwd(this, gps_time, lat, lon)
            height = Core_Utils.linInterpLatLonTime(this.vmf_coeff.ell_height, this.vmf_coeff.first_lat, this.vmf_coeff.d_lat, this.vmf_coeff.first_lon, this.vmf_coeff.d_lon, 1, 1, lat, lon,1);
        end
        
        function tec = interpolateTEC(this, gps_time, lat, lon)
            %interpolate total elecron content
            %
            % SYNTAX
            %   tec = interpolateTEC(this, gps_time, lat, lon)
            tec = Core_Utils.linInterpLatLonTime(this.ionex.data, this.ionex.first_lat, this.ionex.d_lat, this.ionex.first_lon, this.ionex.d_lon, this.ionex.first_time_double, this.ionex.d_t, lat, lon,gps_time);
        end
        
        function thin_shell_height = getThinShellHeight(this)
            % Get Thin Shell Height from Core_Athmosphere
            % 
            % SYNTAX
            %   thin_shell_height = getThinShellHeight(this)

            % ionopshere thin shell height [km]->[m]
            if isempty(this.ionex.height)
                thin_shell_height = 350 * 1e3; % if the ionex is not loaded use 350km
            else
                thin_shell_height = this.ionex.height(1) * 1e3;       % ionopshere thin shell height [km]
            end
        end
        
        function [stec, pp, mfpp, k] = getSTEC(this, lat, lon, az, el, h, time)
            % get slant total electron component
            %
            % SYNTAX
            %   [stec, pp, mfpp, k] = getSTEC(this, lat, lon, az, el, h, time)
            
            thin_shell_height = this.getThinShellHeight();
            
            % get piercing point and mapping function            
            [latpp, lonpp, mfpp, k] = this.getPiercePoint( lat/180*pi, lon/180*pi, h, az(:)/180*pi, el(:)/180*pi, thin_shell_height, 6371000);
            % interpolate TEC at piercing point
            tec = this.interpolateTEC( time, latpp * 180/pi, lonpp * 180/pi);
            
            % apply mapping function
            stec = tec.* mfpp(:);
            if nargout > 1
                pp = [latpp , lonpp];
            end
        end
        
        function foi_delay = getFOIdelayCoeff(this,lat,lon, az,el,h,time)
            % get first horder ionosphere delay
            gps_time = time.getGpsTime();
            foi_delay = zeros(size(el));
            for t = 1: size(el,1)
                idx_sat = find(el(t,:) > 0);
                if length(idx_sat) > 0
                    t_time = gps_time(t);
                    [stec, ~, ~, ~] = this.getSTEC(lat, lon, az(t,idx_sat), el(t,idx_sat), h, t_time);
                    foi_delay(t,idx_sat) = 40.3 * 1e16 .* stec ./ this.V_LIGHT^2; % to be multipleid by wavelength^2
                end
            end
        end
        
        function [hoi_delay2_coeff, hoi_delay3_coeff, bending_coeff, ppo] = getHOIdelayCoeff(this,lat,lon, az,el,h,time)
            % get the coefficient to be multiplied foe a power of the frequency to get the hoi order ionospheric delay and bending
            % hoi_delay2 -> hoi_delay2_coeff * wavelength^3
            % hoi_delay3 -> hoi_delay3_coeff * wavelength^4
            % bending    -> bending_coeff    * wavelength^4
            
            % [1] Fritsche, M., R. Dietrich, C. Kn�fel, A. R�lke, S. Vey, M. Rothacher, and P. Steigenberger. Impact
            % of higher-order ionospheric terms on GPS estimates. Geophysical Research Letters, 32(23),
            % 2005. doi: 10.1029/2005GL024342.
            % [2] Odijk, Dennis. "Fast precise GPS positioning in the presence of ionospheric delays." (2002).
            
            
            hoi_delay2_coeff = zeros(size(el));
            hoi_delay3_coeff = zeros(size(el));
            bending_coeff = zeros(size(el));
            ppo = zeros([size(el)]);
            gps_time = time.getGpsTime();
            t_stec = zeros(size(el));
            t_mfpp = zeros(size(el));
            for t = 1: size(el,1)
                idx_sat = find(el(t,:) > 0);
                if sum(idx_sat > 0)
                    t_time= gps_time(t);
                    [stec, pp, mfpp, k] = this.getSTEC(lat,lon, az(t,idx_sat),el(t,idx_sat),h, t_time);
                    A = 80.6;
                    
                    t_stec(t,idx_sat) = stec;
                    t_mfpp(t,idx_sat) = mfpp;
                    if isempty(this.emf) % do not reload the model each time
                        this.emf = Earth_Magnetic_Field();
                    end
                    b = zeros(length(idx_sat),3);
                    for s = 1:numel(idx_sat)
                        b(s,:) = this.emf.getB(t_time, GPS_SS.ELL_A/1000 + this.ionex.height(1), pp(s,2), pp(s,1));
                    end
                    bok = sum(b.*k,2); %to Tesla
                    c = this.V_LIGHT ;
                    Nemax1 = (20 -6)/(4.55 - 1.38) * 1e12 / 1e18; %in[1]
                    Nemax2 = (3*1e12); % in [2]
                    ni = 0.66;
                    zi = acos(1./mfpp);
                    stec = stec * 1e16;
                    hoi_delay2_coeff(t,idx_sat) =  7527 / c^2  .* bok .* stec;% Eq (10) (11) in [1]
                    hoi_delay3_coeff(t,idx_sat) =  2437 / c^4  .* Nemax1 .* ni .* stec.^2;% Eq (1g) (15) (14) in [1]
                    bending_coeff(t,idx_sat)    =  A^2 ./ (8 .* c^4)  .* tan(zi).^2 .* ni .* Nemax2 .* stec;% Eq(4.34) in [2]
                    ppo(t,idx_sat) = stec;
                end
            end
        end
    end
    
    methods
        %-----------------------------------------------------------
        % TROPO
        %-----------------------------------------------------------
        function [delay] = saastamoinenModel(this, h, undu, el)
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
            % �� undulation is never less than 300 m (on Earth)
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
        
        function [delay] = saastamoinenModelGPT(this, gps_time, lat, lon, h, undu, el)
            % SYNTAX
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
        
        function [delay] = saastamoinenModelPTH(this, gps_time,lat, lon, h, undu, el, P, T, H)
            % SYNTAX
            %   [delay] = Atmosphere.saastamoinen_modelPTH(time_rx, lat, lon, h, undu, el)
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
            
            if isnan(P) || isnan(T) || isnan(H)
                [pres, temp] = this.gpt(gps_time, lat*pi/180, lon*pi/180, h, undu);
                if isnan(P)
                    P = pres;
                end
                if isnan(T)
                    T = temp;
                end
                if isnan(H)
                    H = this.STD_HUMI;
                end
                this.log.addWarning(sprintf('No valid meteo data are present @%s\nUsing standard GPT values \n - %.1f �C\n - %.1f hpa\n - humidity %.1f', datestr(gps_time / 86400 + GPS_Time.GPS_ZERO, 'HH:MM'), T, P, H), 100);
            end
            
            t_h = h;
            ZHD_R = saast_dry(P, t_h, lat);
            ZWD_R = saast_wet(T, H, t_h);
            
            [gmfh_R, gmfw_R] = this.gmf(gps_time, lat*pi/180, lon*pi/180, h - undu, (90-el)*pi/180);
            delay = gmfh_R .* ZHD_R + gmfw_R .* ZWD_R;
        end
        
        function [pres, temp, undu] = gpt(this, gps_time, dlat, dlon, dhgt, undu)
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
        
        function [gmfh, gmfw] = gmf(this, time, lat, lon, hgt, el)
            % This subroutine determines the Global Mapping Functions GMF
            % Reference: Boehm, J., A.E. Niell, P. Tregoning, H. Schuh (2006),
            % Global Mapping Functions (GMF): A new empirical mapping function based on numerical weather model data,
            % Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.
            %
            % input data
            % ----------
            % time: GPS_Time
            % lat: ellipsoidal latitude in radians
            % lon: longitude in radians
            % hgt: orthometric height in m
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
            %
            % SYNTAX
            %   [gmfh, gmfw] = gmf(this, gps_time, dlat, dlon, dhgt, zd)
            doy = time.getMJD()  - 44239 + 1;
            
            pi = 3.14159265359d0;
            
            % -------------------------------------------------------
            % getting ah and aw  from the sperical harmonics expansion
            % --------------------------------------------------------
            cached = (lon == this.lon) && (lat == this.lat) && ~isempty(this.ahm) && ~isempty(this.aha);
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
                x = cos(lat)*cos(lon);
                y = cos(lat)*sin(lon);
                z = sin(lat);
                
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
                this.lat = lat;
                this.lon = lon;
            end
            if cached
                ahm = this.ahm;
                aha = this.aha;
                
                awm = this.awm;
                awa = this.awa;
            else
                % hysdrostatic
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
                
                % wet
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
            ah  = (ahm + aha*cos(doy*2.d0*pi))*1d-5;
            % --------------------------------------------------------
            [bh, bw, ch, cw] = this.getGMFVMFBC(time, lat);
            
            % compute the mapping functions
            [gmfh] = this.mfContinuedFractionForm(repmat(ah,1,size(el,2)),bh,repmat(ch,1,size(el,2)),el);
            [gmfw] = this.mfContinuedFractionForm(repmat(aw,1,size(el,2)),bw,cw,el);
            % correct hydrostatic for height
            [ht_corr] = this.hydrostaticMFHeigthCorrection(hgt,el);
            gmfh       = nan2zero(gmfh + ht_corr);
            gmfh       = nan2zero(gmfh + ht_corr);
        end
        
        function [gmfh, gmfw] = vmf_grd(this, time, lat, lon, el, h_ell, interp_first)
            %angles in radians!!
            %code based on:
            %    [1]  Boehm, J., B. Werl, H. Schuh (2006),  Troposphere mapping functions for GPS and very long baseline interferometry  from European Centre for Medium-Range Weather Forecasts operational analysis data,J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
            %    [2] Kouba, Jan. "Implementation and testing of the gridded Vienna Mapping Function 1 (VMF1)." Journal of Geodesy 82.4-5 (2008): 193-205.
            %    [3] Niell, A. E. "Global mapping functions for the atmosphere delay at radio wavelengths." Journal of Geophysical Research: Solid Earth 101.B2 (1996): 3227-3246.
            %    [4] http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/vmf1_grd.m
            %
            % SYNTAX
            %   [gmfh, gmfw] = vmf(this, gps_time, lat, lon, zd)
            if nargin < 7
                interp_first = false;
            end
            if interp_first % interp first
                [ah, aw] = this.interpolateAlpha(time.getGpsTime(), lat/pi*180, lon/pi*180);
                
                
                %             h_ell_vmf = this.interpolateVMFElHeight(lat*180/pi,lon*180/pi); % get height of the station
                %             % eq (6) in [2]
                %             aw = aw - 4 * 1e-8 * (h_ell - h_ell_vmf);
                
                [bh, bw, ch, cw] = this.getGMFVMFBC(time,lat);
                [gmfh] = this.mfContinuedFractionForm(repmat(ah,1,size(el,2)),bh,repmat(ch,1,size(el,2)),el);
                [gmfw] = this.mfContinuedFractionForm(repmat(aw,1,size(el,2)),bw,cw,el);
                
                [ht_corr] = this.hydrostaticMFHeigthCorrection(h_ell,el);
                gmfh       = nan2zero(gmfh + ht_corr);
            else % compute mapping function first
                % get the index of the interpolating points
                if isempty(this.vmf_coeff.first_lat)
                    this.initVMF(time.first, time.last);
                end
                if isempty(this.vmf_coeff.first_lat)
                    gmfh = [];
                    gmfw = [];
                else
                    [ it, st, ilons, ilone, slon, ilat, slat] = this.getVMFIndex(time.getGpsTime(), lat/pi*180, lon/pi*180);
                    ah_calc_1 = this.vmf_coeff.ah([ilat ilat+1], [ilons ilone],it);
                    aw_calc_1 = this.vmf_coeff.aw([ilat ilat+1], [ilons ilone],it);
                    ah_calc_2 = this.vmf_coeff.ah([ilat ilat+1], [ilons ilone],it+1);
                    aw_calc_2 = this.vmf_coeff.aw([ilat ilat+1], [ilons ilone],it+1);
                    h_calc =  this.vmf_coeff.ell_height([ilat ilat+1], [ilons ilone]);
                    n_sat = size(el,2);
                    % compute the mapping function and the heoght correction for all the inteprolating points
                    [bh, bw, ch, cw] = this.getGMFVMFBC(time, lat);
                    [gmfh_11_1] = this.mfContinuedFractionForm(repmat(squeeze(ah_calc_1(1,1,:)),1,n_sat),bh,repmat(ch,1,n_sat),el);
                    [gmfw_11_1] = this.mfContinuedFractionForm(repmat(squeeze(aw_calc_1(1,1,:)),1,n_sat),bw,cw,el);
                    [gmfh_12_1] = this.mfContinuedFractionForm(repmat(squeeze(ah_calc_1(1,2,:)),1,n_sat),bh,repmat(ch,1,n_sat),el);
                    [gmfw_12_1] = this.mfContinuedFractionForm(repmat(squeeze(aw_calc_1(1,2,:)),1,n_sat),bw,cw,el);
                    [gmfh_21_1] = this.mfContinuedFractionForm(repmat(squeeze(ah_calc_1(2,1,:)),1,n_sat),bh,repmat(ch,1,n_sat),el);
                    [gmfw_21_1] = this.mfContinuedFractionForm(repmat(squeeze(aw_calc_1(2,1,:)),1,n_sat),bw,cw,el);
                    [gmfh_22_1] = this.mfContinuedFractionForm(repmat(squeeze(ah_calc_1(2,2,:)),1,n_sat),bh,repmat(ch,1,n_sat),el);
                    [gmfw_22_1] = this.mfContinuedFractionForm(repmat(squeeze(aw_calc_1(2,2,:)),1,n_sat),bw,cw,el);
                    
                    [gmfh_11_2] = this.mfContinuedFractionForm(repmat(squeeze(ah_calc_2(1,1,:)),1,n_sat),bh,repmat(ch,1,n_sat),el);
                    [gmfw_11_2] = this.mfContinuedFractionForm(repmat(squeeze(aw_calc_2(1,1,:)),1,n_sat),bw,cw,el);
                    [gmfh_12_2] = this.mfContinuedFractionForm(repmat(squeeze(ah_calc_2(1,2,:)),1,n_sat),bh,repmat(ch,1,n_sat),el);
                    [gmfw_12_2] = this.mfContinuedFractionForm(repmat(squeeze(aw_calc_2(1,2,:)),1,n_sat),bw,cw,el);
                    [gmfh_21_2] = this.mfContinuedFractionForm(repmat(squeeze(ah_calc_2(2,1,:)),1,n_sat),bh,repmat(ch,1,n_sat),el);
                    [gmfw_21_2] = this.mfContinuedFractionForm(repmat(squeeze(aw_calc_2(2,1,:)),1,n_sat),bw,cw,el);
                    [gmfh_22_2] = this.mfContinuedFractionForm(repmat(squeeze(ah_calc_2(2,2,:)),1,n_sat),bh,repmat(ch,1,n_sat),el);
                    [gmfw_22_2] = this.mfContinuedFractionForm(repmat(squeeze(aw_calc_2(2,2,:)),1,n_sat),bw,cw,el);
                    
                    
                    % copute the height correction for the mapping function with the height of the four points see [3]
                    [ht_corr_11] = this.hydrostaticMFHeigthCorrection(h_ell,el);
                    [ht_corr_12] = this.hydrostaticMFHeigthCorrection(h_ell,el);
                    [ht_corr_21] = this.hydrostaticMFHeigthCorrection(h_ell,el);
                    [ht_corr_22] = this.hydrostaticMFHeigthCorrection(h_ell,el);
                    
                    gmfh_11_1 = gmfh_11_1 + ht_corr_11;
                    gmfh_12_1 = gmfh_12_1 + ht_corr_12;
                    gmfh_21_1 = gmfh_21_1 + ht_corr_21;
                    gmfh_22_1 = gmfh_22_1 + ht_corr_22;
                    
                    gmfh_11_2 = gmfh_11_2 + ht_corr_11;
                    gmfh_12_2 = gmfh_12_2 + ht_corr_12;
                    gmfh_21_2 = gmfh_21_2 + ht_corr_21;
                    gmfh_22_2 = gmfh_22_2 + ht_corr_22;
                    %interpolate along lon
                    valbu = gmfh_11_1.*(1-slon) + gmfh_12_1.*slon;
                    valau = gmfh_11_2.*(1-slon) + gmfh_12_2.*slon;
                    valbd = gmfh_21_1.*(1-slon) + gmfh_22_1.*slon;
                    valad = gmfh_21_2.*(1-slon) + gmfh_22_2.*slon;
                    
                    %interpolate along lat
                    valb = valbd.*(1-slat) + valbu.*slat;
                    vala = valad.*(1-slat) + valau.*slat;
                    
                    %interpolate along time
                    gmfh = valb.*repmat((1-st),1,n_sat) + vala.*repmat(st,1,n_sat);
                    
                    %interpolate along lon
                    valbu = gmfw_11_1.*(1-slon) + gmfw_12_1.*slon;
                    valau = gmfw_11_2.*(1-slon) + gmfw_12_2.*slon;
                    valbd = gmfw_21_1.*(1-slon) + gmfw_22_1.*slon;
                    valad = gmfw_21_2.*(1-slon) + gmfw_22_2.*slon;
                    
                    %interpolate along lat
                    valb = valbd.*(1-slat) + valbu.*slat;
                    vala = valad.*(1-slat) + valau.*slat;
                    
                    %interpolate along time
                    gmfw = valb.*repmat((1-st),1,n_sat) + vala.*repmat(st,1,n_sat);
                end
            end
        end               
        
        function [gmfh, gmfw] = niell(this, time, lat, el, h_ell)
            %angles in radians!!
            %code based on:
            %    [3] Niell, A. E. "Global mapping functions for the atmosphere delay at radio wavelengths." Journal of Geophysical Research: Solid Earth 101.B2 (1996): 3227-3246.
            %
            % SYNTAX
            %   [gmfh, gmfw] = vmf(this, gps_time, lat, lon, zd)
            T0 = 28  - (lat<0)*365.25/2;
            lat = abs(lat);
            lat_p = [ 15 30 45 60 75]/180*pi;
            ah_coef_lat = [1.2769934e-3 1.2683230e-3 1.2465397e-3 1.2196049e-3 1.2045996e-3];
            bh_coef_lat = [2.9153695e-3 2.9152299e-3 2.9288445e-3 2.9022565e-3 2.9024912e-3];
            ch_coef_lat = [62.610505e-3 62.837393e-3 63.721774e-3 63.824265e-3 64.258455e-3];
            
            ah_coef_amp_lat = [0.0 1.2709626e-5 2.6523662e-5 3.4000452e-5 4.1202191e-5];
            bh_coef_amp_lat = [0.0 2.1414979e-5 3.0160779e-5 7.2562722e-5 11.723375e-5];
            ch_coef_amp_lat = [0.0 9.0128400e-5 4.3497037e-5 84.795348e-5 170.37206e-5];
            
            aw_coef_lat = [5.8021897e-4 5.6794847e-4 5.8118019e-4 5.9727542e-4 6.1641693e-4];
            bw_coef_lat = [1.4275268e-3 1.5138625e-3 1.4572752e-3 1.5007428e-3 1.7599082e-3];
            cw_coef_lat = [4.3472961e-2 4.6729510e-2 4.3908931e-2 4.4626982e-2 5.4736038e-2];
            
            [~,~,idx_lat1] = histcounts(lat,[0 lat_p pi/2]);
            if idx_lat1 == 1
                ll = 1;
                idx_lat2 = idx_lat1;
            elseif idx_lat1 == 6
                ll = 1;
                idx_lat1 = idx_lat1-1;
                idx_lat2 = idx_lat1;
            else
                idx_lat1 = idx_lat1-1;
                idx_lat2 = idx_lat1 + 1;
                ll = (lat - lat_p(idx_lat1))/( lat_p(idx_lat2) -  lat_p(idx_lat1));
            end
            ah_coef = ah_coef_lat(idx_lat1)*ll + (1-ll)*ah_coef_lat(idx_lat2);
            bh_coef = bh_coef_lat(idx_lat1)*ll + (1-ll)*bh_coef_lat(idx_lat2);
            ch_coef = ch_coef_lat(idx_lat1)*ll + (1-ll)*ch_coef_lat(idx_lat2);
            
            ah_coef_amp = ah_coef_amp_lat(idx_lat1)*ll + (1-ll)*ah_coef_amp_lat(idx_lat2);
            bh_coef_amp = bh_coef_amp_lat(idx_lat1)*ll + (1-ll)*bh_coef_amp_lat(idx_lat2);
            ch_coef_amp = ch_coef_amp_lat(idx_lat1)*ll + (1-ll)*ch_coef_amp_lat(idx_lat2);
            
            aw = aw_coef_lat(idx_lat1)*ll + (1-ll)*aw_coef_lat(idx_lat2);
            bw = bw_coef_lat(idx_lat1)*ll + (1-ll)*bw_coef_lat(idx_lat2);
            cw = cw_coef_lat(idx_lat1)*ll + (1-ll)*cw_coef_lat(idx_lat2);
            
            [~,doy] = time.getDOY;
            period = (doy - T0 ) / 365.25 * 2 *pi;
            
            ah = ah_coef + ah_coef_amp * cos(period);
            bh = bh_coef + bh_coef_amp * cos(period);
            ch = ch_coef + ch_coef_amp * cos(period);
            
            aw = aw*ones(size(ah));
            bw = bw*ones(size(bh));
            cw = cw*ones(size(ch));
            n_sat = size(el,2);
            [gmfh] = this.mfContinuedFractionForm(repmat(ah,1,n_sat),repmat(bh,1,n_sat),repmat(ch,1,n_sat),el);
            [gmfw] = this.mfContinuedFractionForm(repmat(aw,1,n_sat),repmat(bw,1,n_sat),repmat(cw,1,n_sat),el);
            
            [h_h_coorection] = this.hydrostaticMFHeigthCorrection(h_ell,el);
            gmfh = gmfh + h_h_coorection;
        end
                        
        %-----------------------------------------------------------
        % IONO
        %-----------------------------------------------------------
        function [iono_mf] = getIonoMF(this, lat_rad, h_ortho, el_rad, rcm)
            % Get the pierce point
            % INPUT:
            %   lat_rad             latitude of the receiver           [rad]
            %   h_ortho             orthometric height of the receiver [m]
            %   el_rad              elevation of the satellites        [rad]
            %   rcm                 meridian radius curvature <optional>
            %
            % OUTPUT
            %   iono_mf             iono mapping function
            % SOURCES
            %   [1] Handbook of Global Navigation Satellite System (2017)
            %
            % SYNTAX
            %   [iono_mf] = getIonoMF(lat_rad, h_ortho, el_rad, rcm)
            
            % Get radius of curvature at lat
            if nargin < 5
                rcm = getMeridianRadiusCurvature(lat_rad);
            end
            
            thin_shell_height = this.getThinShellHeight();
            
            id_ok = el_rad >= 0 & ~isnan(el_rad);
            k = (rcm + h_ortho)/(rcm + h_ortho + thin_shell_height) * cos(el_rad(id_ok));
            iono_mf = nan(size(el_rad));
            iono_mf(id_ok) = (1-(k).^2).^(-1/2); % formula 6.99 in [1]
        end
    end
    
    methods (Static)
        %-----------------------------------------------------------
        % TROPO
        %-----------------------------------------------------------
        function [delay] = mfContinuedFractionForm(a,b,c,el)
            sine = sin(el);
            delay = (1 + (a ./ (1 + (b ./ (1 + c) )))) ...
                ./ (sine + (a ./ (sine + (b ./ (sine + c) ))));
        end
        
        function [ht_corr] = hydrostaticMFHeigthCorrection(h_ell,el)
            % coorect the hysdrostatic mapping functions for the height
            % formulas and paramaters taken from :
            %   Niell, A. E. "Global mapping functions for the atmosphere delay at radio wavelengths." Journal of Geophysical Research: Solid Earth 101.B2 (1996): 3227-3246.
            
            % height correction for the hydrostatic part
            
            % coefficent from tab 3
            a_ht = 2.53d-5;
            b_ht = 5.49d-3;
            c_ht = 1.14d-3;
            h_ell_km     = h_ell/1000;   % convert height to km
            % eq (6)
            ht_corr_coef = 1 ./ sin(zero2nan(el))   -    Atmosphere.mfContinuedFractionForm(a_ht,b_ht,c_ht,el);
            % eq (7)
            ht_corr      = ht_corr_coef * h_ell_km;
        end
        
        function [bh, bw, ch, cw] = getGMFVMFBC(time, lat)
            % get the coefficients b and c of the continued fraction form, both for GMF and VMF
            % INPUT:
            %    time - time GPS_Time object
            %    lat  - latitude in radians
            % coefficients of tab 1 in [1]
            if (lat<0)      % southern hemisphere
                phi_h  = pi;
                c11_h = 0.007;
                c10_h = 0.002;
            else             % northern hemisphere
                phi_h  = 0.d0;
                c11_h = 0.005;
                c10_h = 0.001;
            end
            c0_h = 0.062;
            % hidrostatic b form Isobaric mapping function
            bh = 0.002905;
            
            doy = time.getMJD()  - 44239 + 1;
            % c hydrostatic is taken from equation (7) in [1]
            ch = c0_h + ((cos((doy - 28) / 365.25 * 2 * pi + phi_h) + 1) * c11_h / 2 + c10_h)*(1 - cos(lat));
            % wet b and c form Niell mapping function at 45� lat tab 4 in [3]
            bw = 0.00146;
            cw = 0.04391;
        end
        
        function zhd_corr = getZenithDelayCorrection(coo1, pr1, coo2, pr2)
            %
            % Compute ztd_corr and iwv_corr at one receiver height from another
            %
            % Corresponding about the script: alessandra.mascitelli (at) polimi.it
            %
            % OUTPUT
            %   ztd_corr     zenith total delay correction [mm]
            %   pwv_corr     integrated water vapor correction [mm]
            %
            tc2tk = 273.15; %temperature conversion from Celsius to Kelvin
            Tisa = 288.15; %temperature International Standard Atmosphere [K]
            R = 8.31432;  %universal gas constant for air [J K^-1 mol^-1]
            rv = 461.5; %specific gas constant - Water vapor [N m/kg K]
            md = 0.0289644; %molar mass of dry air [kg mol^-1]
            g = 9.80665; %gravitational acceleration [m/s^2]
            k3 = 377600; %constant[K^2/mbar]
            k2 = 17; %constant[K/mbar]
            %
            [lat1, ~, h_ellips1, h_ortho1] = coo1.getGeodetic();
            [lat2, ~, h_ellips2, h_ortho2] = coo2.getGeodetic();
            %
            %if meteo data from raob are not available, it is possible to compute it from gnss in
            %the following way from Berberan-Santos et al. 1997 and Bai et al. 2003
            if nargin < 4
                %compute p_raob from p_gnss
                pr2 = pr1 * exp(-(g * md * (h_ortho2 - h_ortho1)) / (R * Tisa));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %compute zhd_gnss [m] at gnss receiver height
            %compute function f_gnss (f from Bevis et al. 1992)
            f1 = 1 - 0.00266 * cos(2 * lat1) - 0.00028 * (h_ellips1 / 1000); %in this formula h_ell [km] and lat_gnss [rad]
            %
            zhd1 = (2.2779 * pr1 / f1) / 1000; %zhd_gnss [mm] --> [m]
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %compute zhd_raob [m] at raob receiver height
            %compute function f_raob (f from Bevis et al. 1992)
            f2 = 1 - 0.00266 * cos (2 * lat2) - 0.00028 * (h_ellips2 / 1000); %in this formula h_ell [km] and lat_raob [rad]
            %
            zhd2 = (2.2779*pr2/f2)/1000; %zhd_raob [mm] --> [m]
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            
            % corrections
            zhd_corr = zhd2 - zhd1;
        end

        %-----------------------------------------------------------
        % IONO
        %-----------------------------------------------------------
        function [delay] = klobucharModel(lat, lon, az, el, sow, ionoparams)
            % SYNTAX
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
            delay(index,1) = Core_Utils.V_LIGHT * f(index) .* (5e-9 + a(index) .* (1 - (x(index).^2)/2 + (x(index).^4)/24));
            
            index = find(abs(x) >= 1.57);
            delay(index,1) = Core_Utils.V_LIGHT * f(index) .* 5e-9;
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
            %   k                   direction of the ray in ecef coordinates [X Y Z]
            %
            % SYNTAX
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
                iono_mf = (1-(k).^2).^(-1/2);
                iono_mf = reshape(iono_mf, input_size(1), input_size(2));
            end
            
            if nargout > 3
                k = [-sin(az_rad).*cos(el_rad) ...
                    -cos(az_rad).*cos(el_rad) ...
                    -sin(el_rad)];
                % go to global system
                k = local2globalVel2(k', lon_rad,lat_rad)';
            end
            
            lat_pp = reshape(lat_pp, input_size(1), input_size(2));
            lon_pp = reshape(lon_pp, input_size(1), input_size(2));
        end
        
        function [ZWD] = saast_wet(T, H,h)
            % Saastamoinen wet
            %
            % SYNTAX
            %   [ZWD] = saast_wet(T, H);
            %
            % INPUT:
            %   T = air temperature
            %   H = humidity
            %
            % OUTPUT:
            %   ZWD = Zenith Wet Delay
            %
            % DESCRIPTION:
            %   Zenith Wet Delay (ZWD) computation by Saastamoinen model.
            
            
            % Convert C -> K
            T = T + 273.15;
            
            % height correction must be done before
            % (keep the following line commented)
            % H = H * exp(-0.0006396 * h);
            % Convert humidity
            H = H./100;
            
            c = -37.2465 + 0.213166 * T - 2.56908 * (10^-4) * (T.^2);
            e = H .* exp(c);
            
            %ZWD (Saastamoinen model)
            ZWD = 0.0022768 * (((1255 ./ T) + 0.05) .* e);
        end
        
        function [ZHD] = saast_dry(P, h, lat)
            % Saastamoinen dry
            %
            % SYNTAX
            %   [ZHD] = saast_dry(P, h, lat);
            %
            % INPUT:
            %   P = atmospheric pressure [hPa]
            %   h = orthometric height [m]
            %   lat = latitude [deg]
            %
            % OUTPUT:
            %   ZHD = Zenith Hydrostatic Delay
            %
            % DESCRIPTION:
            %   Zenith Hydrostatic Delay (ZHD) computation by Saastamoinen model.
            
            
            %ZHD (Saastamoinen model)
            ZHD = 0.0022768 * P(:) .* (1 + 0.00266 * cosd(2*lat(:)) + 0.00000028 * h(:));
            %ZHD = 0.0022767 * P(:) ./ (1 - 0.00266 * cosd(2*lat(:)) - 0.00000028 * h(:));
        end
        
        function [cotan_term] = macmillanGrad(el)
            % formual from macmillan
            %
            % INPUT:
            %  el  elevation [rad]
            %
            % SYNTAX:
            %  [cotan_term] = Atmosphere.macmillanGrad(el)
            %
            % NOTE:
            %    angle in radians !!!
            cotan_term = cot(el);
        end
        
        function [cotan_term] = chenHerringGrad(el)
            % forumal from chen and herring
            % coefficient from chao
            %
            % INPUT:
            %  el  elevation [rad]
            %
            % SYNTAX:
            %  [cotan_term] = Atmosphere.chenHerringGrad(el)
            %
            % NOTE:
            %    angle in radians !!!
            cotan_term = 1 ./ ( sin(el).*tan(el) + 0.0032);
        end
        
    end
end
