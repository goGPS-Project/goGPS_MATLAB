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
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:        Giulio Tagliaferro
%  Contributors:      Giulio Tagliaferro, Andrea Gatti ...
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
        STD_TEMP = 291.15;
        STD_PRES = 1013.25;
        STD_HUMI = 50;
    end
    
    properties  (SetAccess = private, GetAccess = public)
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
        zhoic = struct(...
            'data', [],... % [lat x lon x zernike x time]
            'first_lat',  -87.5, ...    % first latitude
            'first_lon',  0, ...    % first longitude
            'radial_order', [], ...   % radial order 
            'azimuthal_order',[], ...  % azimuthal order
            'd_lat',      2.5, ...    % lat spacing
            'd_lon',      5, ...    % lon_spacing
            'n_lat',      71, ...    % num lat
            'n_lon',      72, ...    % num lon
            'first_time', [], ...    % times [time] of the maps
            'first_time_double', [], ...    % times [time] of the maps [seconds from GPS zero]
            'dt',         3600, ...    % time spacing
            'n_t',        [] ...   % num of epocvhs
            )
        
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
            fid = fopen(file_name,'rt');
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
            fname = Core.getState.getIonoFileName( dsa, dso);
            for i = 1 : length(fname)
                this.importIonex(fname{i});
            end
        end
        
        function importZHOIC(this,fname)
                      
            PP=load(fname,'-ascii');
            data_tmp = reshape(PP,71,72,16);
            [~,NAME,EXT] = fileparts(fname);
            time = GPS_Time([str2num(NAME(5:8)) str2num(NAME(9:10)) 1 str2num(EXT(3:4)) 0 0]);
            if isempty(this.zhoic.data)
                this.zhoic.data = data_tmp;
                this.zhoic.first_time = time;
                this.zhoic.first_time_double = time.getMatlabTime;
                this.zhoic.nt = 1; 
            else
                this.zhoic.data = cat(4,this.zhoic.data,data_tmp);
                this.zhoic.nt = this.zhoic.nt + 1; 
            end

        end
        
        function initZHOIC(this,dsa,dso)
            % initialise zhoic map importing all the necessary files
            %
            % SYNTAX
            %   initZHOIC(this, dsa, dso)
            dso = dso.getCopy();
            dsa = dsa.getCopy();
            dso.addSeconds(6*3600);
            fnp = File_Name_Processor();
            state = Core.getCurrentSettings();
            fname = fnp.dateKeyRepBatch([state.iono_dir '/IFCz${YYYY}${MM}.H${HH}'], dsa, dso);
            for i = 1 : length(fname)
                this.importZHOIC(fname{i});
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
            state = Core.getState;
            fname = state.getVMFFileName( dsa, dso);
            % import the coefficeints files
            for i = 1 : length(fname)
                this.importVMFCoeffFile(fname{i});
            end
            h_fname = state.getVMFHeightFileName();
            fid = fopen(h_fname, 'rt');
            if strcmpi(state.vmf_res,'2.5x2')
                fgetl(fid); %skip one line header
                formatSpec = [repmat([repmat(' %f',1,10) '\n'],1,14) repmat(' %f',1,5)];
                h_data = cell2mat(textscan(fid,formatSpec,91));
                h_data(:,end) = [];
            elseif strcmpi(state.vmf_res,'5x5')
                formatSpec = '%f';
                h_data = textscan(fid,formatSpec); h_data = h_data{1,1};
                h_data = reshape(h_data,72,36)';
            elseif strcmpi(state.vmf_res,'1x1')
                formatSpec = '%f';
                h_data = textscan(fid,formatSpec); h_data = h_data{1,1};
                h_data = reshape(h_data,360,180)';
            end
                
            this.vmf_coeff.ell_height = h_data;
            fclose(fid);
        end
        
        function importTidalAtmLoadHarmonics(this)
            % importing Tidal Atm and loading Harmonics
            %
            % SYNTAX
            %   importTidalAtmLoadHarmonics(this)
            fname = Core.getState.getTAtmLoadFileName();
            data = importdata(fname);
            this.atm_load_t.harmonics = permute(reshape(data(:,3:end),360 ,180, 18),[2 1 3])/1e3;
        end
        
        function importAtmLoadCoeffFile(this, filename)
            % import data of atmospehric loading file
            %
            % SYNTAX
            %   importAtmLoadCoeffFile(this, filename)
            fid = fopen([filename],'rt');
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
            state = Core.getCurrentSettings();
            fnp = File_Name_Processor;
            fid = fopen(ffile_name,'rt');
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
                    if strcmpi(state.vmf_res,'1x1')
                        n_lat = 180;
                        n_lon = 360;
                        d_lat = -1;
                        d_lon = 1;
                        first_lat = 0.5;
                        first_lon = 0.5;
                    elseif strcmpi(state.vmf_res,'2.5x2')
                        n_lat = 91;
                        n_lon = 144;
                        d_lat = -1;
                        d_lon = 1;
                        first_lat = 90;
                        first_lon = 0.0;
                    elseif strcmpi(state.vmf_res,'5x5')
                        n_lat = 36;
                        n_lon = 72;
                        d_lat = -5;
                        d_lon = 5;
                        first_lat = 2.5;
                        first_lon = 2.5;
                    end
                    data = sscanf(txt(lim(eoh + 1, 1) : lim(end,2)), '%f ');
                    ah =  reshape(data(3:6:end), n_lon, n_lat)';
                    aw =  reshape(data(4:6:end), n_lon, n_lat)';
                    zhd = reshape(data(5:6:end), n_lon, n_lat)';
                    zwd = reshape(data(6:6:end), n_lon, n_lat)';
                    
                    if isempty_obj
                        this.vmf_coeff.ah         = ah;
                        this.vmf_coeff.aw         = aw;
                        this.vmf_coeff.zhd        = zhd;
                        this.vmf_coeff.zwd        = zwd;
                        this.vmf_coeff.first_lat  = first_lat;
                        this.vmf_coeff.first_lon  = first_lon;
                        this.vmf_coeff.d_lat      = d_lat;
                        this.vmf_coeff.d_lon      = d_lon;
                        this.vmf_coeff.n_lat      = n_lat;
                        this.vmf_coeff.n_lon      = n_lon;
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
            fid = fopen(filename, 'rt');
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
        
        function [it, st, ilons, ilone, slon, ilat, slat] = getVMFIndex(this, gps_time, lat, lon)
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
            if isempty(this.ionex) || isempty(this.ionex.data)
                tec = 0;
            else
                tec = Core_Utils.linInterpLatLonTime(this.ionex.data, this.ionex.first_lat, this.ionex.d_lat, this.ionex.first_lon, this.ionex.d_lon, this.ionex.first_time_double, this.ionex.d_t, lat, lon,gps_time);
            end
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
            if nargout > 3
               [latpp, lonpp, mfpp, k] = this.getPiercePoint( lat/180*pi, lon/180*pi, h, az(:)/180*pi, el(:)/180*pi, thin_shell_height, 6371000);
            else
               [latpp, lonpp, mfpp] = this.getPiercePoint( lat/180*pi, lon/180*pi, h, az(:)/180*pi, el(:)/180*pi, thin_shell_height, 6371000);
            end
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
                    [stec, ~, ~] = this.getSTEC(lat, lon, az(t,idx_sat), el(t,idx_sat), h, t_time);
                    foi_delay(t,idx_sat) = 40.3 * 1e16 .* stec ./ this.V_LIGHT^2; % to be multipleid by wavelength^2
                end
            end
        end
        
        function [hoi_delay2_coeff, hoi_delay3_coeff, bending_coeff, ppo] = getHOIdelayCoeff(this,lat,lon, az,el,h,time)
            % get the coefficient to be multiplied foe a power of the frequency to get the hoi order ionospheric delay and bending
            % hoi_delay2 -> hoi_delay2_coeff * wavelength^3
            % hoi_delay3 -> hoi_delay3_coeff * wavelength^4
            % bending    -> bending_coeff    * wavelength^4
            
            % [1] Fritsche, M., R. Dietrich, C. Knÿfel, A. Rÿlke, S. Vey, M. Rothacher, and P. Steigenberger. Impact
            % of higher-order ionospheric terms on GPS estimates. Geophysical Research Letters, 32(23),
            % 2005. doi: 10.1029/2005GL024342.
            % [2] Odijk, Dennis. "Fast precise GPS positioning in the presence of ionospheric delays." (2002).
            
            
            hoi_delay2_coeff = zeros(size(el));
            hoi_delay3_coeff = zeros(size(el));
            bending_coeff = zeros(size(el));
            ppo = zeros([size(el)]);
            if not(isempty(this.ionex)) && not(isempty(this.ionex.data))
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
        
        function [hoi_if12] = getZHOICdelayCoeff(this,lat,lon, az,el,time)
            hoi_if12 = nan(size(az));
            lon(lon<0) = lon(lon<0)+360;
            time_d = rem(time - this.zhoic.first_time,86400);
            for i = 1 : time.length
                idx_t = round(time_d(i)/this.zhoic.dt)+1;
                idx_lat = round((lat - this.zhoic.first_lat)/this.zhoic.d_lat)+1;
                idx_lon = round((lon - this.zhoic.first_lon)/this.zhoic.d_lon)+1;
                z_coeff = squeeze(this.zhoic.data(idx_lat,idx_lon,:,idx_t));
                hoi_if12(i,:) = ifczevalX(el(i,:),az(i,:),z_coeff)/1e3;
            end
            function re = ifczevalX(ele,azi,x)
                m=length(ele);
                r    =1-sin(ele/180*pi);
                theta=      azi/180*pi;
                Z(1:m,1 )=ones(1,m);
                Z(1:m,2 )=r.*cos(theta);
                Z(1:m,3 )=r.*sin(theta);
                Z(1:m,4 )=(2*r.^2 - 1);
                Z(1:m,5 )=r.^2 .* cos(2*theta);
                Z(1:m,6 )=r.^2 .* sin(2*theta);
                Z(1:m,7 )=(3*r.^3 - 2*r) .* cos(theta);
                Z(1:m,8 )=(3*r.^3 - 2*r) .* sin(theta);
                Z(1:m,9 )= 6*r.^4 - 6*r.^2 + 1;
                Z(1:m,10)=r.^3 .* cos(3*theta);
                Z(1:m,11)=r.^3 .* sin(3*theta);
                Z(1:m,12)=(4*r.^4 - 3*r.^2) .* cos(2*theta);
                Z(1:m,13)=(4*r.^4 - 3*r.^2) .* sin(2*theta);
                Z(1:m,14)=(3*r-12*r.^3+10*r.^5).* cos(theta);
                Z(1:m,15)=(3*r-12*r.^3+10*r.^5).* sin(theta);
                Z(1:m,16)=20*r.^6 - 30*r.^4 + 12*r.^2 - 1;
                re=zeros(m,1);
                for g=1:16
                    re(1:m)=re(1:m)+Z(1:m,g)*x(g);
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
                this.log.addWarning(sprintf('No valid meteo data are present @%s\nUsing standard GPT values \n - %.1f ÿC\n - %.1f hpa\n - humidity %.1f', datestr(gps_time / 86400 + GPS_Time.GPS_ZERO, 'HH:MM'), T, P, H), 100);
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
        
        function [gmfh, gmfw] = vmf_grd(this, time, lat, lon, el, h_ell, version, interp_first)
            %angles in radians!!
            %code based on:
            %    [1]  Boehm, J., B. Werl, H. Schuh (2006),  Troposphere mapping functions for GPS and very long baseline interferometry  from European Centre for Medium-Range Weather Forecasts operational analysis data,J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
            %    [2] Kouba, Jan. "Implementation and testing of the gridded Vienna Mapping Function 1 (VMF1)." Journal of Geodesy 82.4-5 (2008): 193-205.
            %    [3] Niell, A. E. "Global mapping functions for the atmosphere delay at radio wavelengths." Journal of Geophysical Research: Solid Earth 101.B2 (1996): 3227-3246.
            %    [4] http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/vmf1_grd.m
            %
            % SYNTAX
            %   [gmfh, gmfw] = vmf(this, gps_time, lat, lon, el, h_ell, version)
            if nargin < 8
                interp_first = false;
            end
            if nargin <7
                version = 1;
            end
            if interp_first % interp first
                [ah, aw] = this.interpolateAlpha(time.getGpsTime(), lat/pi*180, lon/pi*180);
                
                
                %             h_ell_vmf = this.interpolateVMFElHeight(lat*180/pi,lon*180/pi); % get height of the station
                %             % eq (6) in [2]
                %             aw = aw - 4 * 1e-8 * (h_ell - h_ell_vmf);
                if version == 1
                    [bh, bw, ch, cw] = this.getGMFVMFBC(time,lat);
                elseif version == 3
                    [bh, bw, ch, cw] = this.getVMF3BC(time, lat);
                end
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
                    if version == 1
                        [bh, bw, ch, cw] = this.getGMFVMFBC(time, lat);
                    elseif version == 3
                        [bh, bw, ch, cw] = this.getVMF3BC(time, lat,lon);
                    end
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
        
        function [gmfh, gmfw] = herring(this, lat, el, h_el, temperature)
            % angles in radians!!
            % height in meters
            % temperature in celsius
            %
            % code based on:
            %    [3] Herring 1992 Modeling atmospheric delays in the analysis of space geodetic data (eq 8, 9)
            %
            % SYNTAX
            %   [gmfh, gmfw] = herring(this, lat, el, h_el)
            
            ah = (1.2320 + 0.0139 * cos(lat) - 0.0209 * h_el/1000 + 0.00215 * (temperature - 10)) * 1e-3;
            bh = (3.1612 - 0.1600 * cos(lat) - 0.0331 * h_el/1000 + 0.00206 * (temperature - 10)) * 1e-3;
            ch = (71.244 - 4.2930 * cos(lat) - 0.0149 * h_el/1000 - 0.00210 * (temperature - 10)) * 1e-3;
            
            n_sat = size(el,2);
           
            [gmfh] = this.mfContinuedFractionForm(repmat(ah,1,n_sat),repmat(bh,1,n_sat),repmat(ch,1,n_sat),zero2nan(el));
            
            if nargout == 2
                aw = (0.583 + 0.011 * cos(lat) - 0.052 * h_el/1000 + 0.0014 * (temperature - 10)) * 1e-3;
                bw = (1.402 - 0.102 * cos(lat) - 0.101 * h_el/1000 + 0.0020 * (temperature - 10)) * 1e-3;
                cw = (45.85 - 1.910 * cos(lat) - 1.290 * h_el/1000 + 0.0150 * (temperature - 10)) * 1e-3;
                
                [gmfw] = this.mfContinuedFractionForm(repmat(aw,1,n_sat),repmat(bw,1,n_sat),repmat(cw,1,n_sat),zero2nan(el));
            end
        end
        
        function [gmfh, gmfw] = niell(this, time, lat, el, h_ell)
            %angles in radians!!
            %code based on:
            %    [3] Niell, A. E. "Global mapping functions for the atmosphere delay at radio wavelengths." Journal of Geophysical Research: Solid Earth 101.B2 (1996): 3227-3246.
            %
            % SYNTAX
            %   [gmfh, gmfw] = niel(this, gps_time, lat, lon, zd)
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
            % wet b and c form Niell mapping function at 45ÿ lat tab 4 in [3]
            bw = 0.00146;
            cw = 0.04391;
        end
        
        function [bh, bw, ch, cw] = getVMF3BC(time, lat,lon)
            % get the coefficients b and c of the continued fraction form,
            % both for for VMF3 
            % code taken from : https://vmf.geo.tuwien.ac.at/codes/vmf3.m
            % INPUT:
            %    time - time GPS_Time object
            %    lat  - latitude in radians
            % coefficients of tab 1 in [1]
            % convert mjd to doy
            mjd = time.getMJD();
            hour = floor((mjd-floor(mjd))*24);   % get hours
            minu = floor((((mjd-floor(mjd))*24)-hour)*60);   % get minutes
            sec = (((((mjd-floor(mjd))*24)-hour)*60)-minu)*60;   % get seconds
            
            % change secs, min hour whose sec==60
            minu(sec==60) = minu(sec==60)+1;
            sec(sec==60) = 0;
            hour(minu==60) = hour(minu==60)+1;
            minu(minu==60)=0;
            
            % calc jd (yet wrong for hour==24)
            jd = mjd+2400000.5;
            
            % if hr==24, correct jd and set hour==0
            jd(hour==24)=jd(hour==24)+1;
            hour(hour==24)=0;
            
            % integer julian date
            jd_int = floor(jd+0.5);
            
            aa = jd_int+32044;
            bb = floor((4*aa+3)/146097);
            cc = aa-floor((bb*146097)/4);
            dd = floor((4*cc+3)/1461);
            ee = cc-floor((1461*dd)/4);
            mm = floor((5*ee+2)/153);
            
            day = ee-floor((153*mm+2)/5)+1;
            month = mm+3-12*floor(mm/10);
            year = bb*100+dd-4800+floor(mm/10);
            
            % first check if the specified year is leap year or not (logical output)
            leapYear = ((mod(year,4) == 0 & mod(year,100) ~= 0) | mod(year,400) == 0);
            
            days = [31 28 31 30 31 30 31 31 30 31 30 31];
            doy = sum(days(1:month-1)) + day;
%             if leapYear == 1 & month > 2
%                 doy = doy + 1;
%             end
            doy(leapYear == 1 & month > 2) =  doy(leapYear == 1 & month > 2) +1;
            doy = doy + mjd-floor(mjd);   % add decimal places
            
            
            % determine the VMF3 coefficients
            
            % Legendre functions for bh, bw, ch and cw
            anm_bh = [0.00271285863109945 -1.39197786008938e-06 1.34955672002719e-06 2.71686279717968e-07 1.56659301773925e-06;9.80476624811974e-06 -5.83922611260673e-05 -2.07307023860417e-05 1.14628726961148e-06 4.93610283608719e-06;-1.03443106534268e-05 -2.05536138785961e-06 2.09692641914244e-06 -1.55491034130965e-08 -1.89706404675801e-07;-3.00353961749658e-05 2.37284447073503e-05 2.02236885378918e-05 1.69276006349609e-06 8.72156681243892e-07;-7.99121077044035e-07 -5.39048313389504e-06 -4.21234502039861e-06 -2.70944149806894e-06 -6.80894455531746e-07;7.51439609883296e-07 3.85509708865520e-07 4.41508016098164e-08 -2.07507808307757e-08 4.95354985050743e-08;2.21790962160087e-05 -5.56986238775212e-05 -1.81287885563308e-05 -4.41076013532589e-06 4.93573223917278e-06;-4.47639989737328e-06 -2.60452893072120e-06 2.56376320011189e-06 4.41600992220479e-07 2.93437730332869e-07;8.14992682244945e-07 2.03945571424434e-07 1.11832498659806e-08 3.25756664234497e-08 3.01029040414968e-08;-7.96927680907488e-08 -3.66953150925865e-08 -6.74742632186619e-09 -1.30315731273651e-08 -2.00748924306947e-09;-2.16138375166934e-05 1.67350317962556e-05 1.93768260076821e-05 1.99595120161850e-06 -2.42463528222014e-06;5.34360283708044e-07 -3.64189022040600e-06 -2.99935375194279e-06 -2.06880962903922e-06 -9.40815692626002e-07;6.80235884441822e-07 1.33023436079845e-07 -1.80349593705226e-08 2.51276252565192e-08 -1.43240592002794e-09;-7.13790897253802e-08 7.81998506267559e-09 1.13826909570178e-09 -5.89629600214654e-09 -4.20760865522804e-09;-5.80109372399116e-09 1.13702284491976e-09 7.29046067602764e-10 -9.10468988754012e-10 -2.58814364808642e-10;1.75558618192965e-05 -2.85579168876063e-05 -1.47442190284602e-05 -6.29300414335248e-06 -5.12204538913460e-07;-1.90788558291310e-06 -1.62144845155361e-06 7.57239241641566e-07 6.93365788711348e-07 6.88855644570695e-07;2.27050351488552e-07 1.03925791277660e-07 -3.31105076632079e-09 2.88065761026675e-08 -8.00256848229136e-09;-2.77028851807614e-08 -5.96251132206930e-09 2.95987495527251e-10 -5.87644249625625e-09 -3.28803981542337e-09;-1.89918479865558e-08 3.54083436578857e-09 8.10617835854935e-10 4.99207055948336e-10 -1.52691648387663e-10;1.04022499586096e-09 -2.36437143845013e-10 -2.25110813484842e-10 -7.39850069252329e-11 7.95929405440911e-11;-3.11579421267630e-05 -3.43576336877494e-06 5.81663608263384e-06 8.31534700351802e-07 4.02619520312154e-06;6.00037066879001e-07 -1.12538760056168e-07 -3.86745332115590e-07 -3.88218746020826e-07 -6.83764967176388e-07;-9.79583981249316e-08 9.14964449851003e-08 4.77779838549237e-09 2.44283811750703e-09 -6.26361079345158e-09;-2.37742207548109e-08 -5.53336301671633e-09 -3.73625445257115e-09 -1.92304189572886e-09 -7.18681390197449e-09;-6.58203463929583e-09 9.28456148541896e-10 2.47218904311077e-10 1.10664919110218e-10 -4.20390976974043e-11;9.45857603373426e-10 -3.29683402990254e-11 -8.15440375865127e-11 -1.21615589356628e-12 -9.70713008848085e-12;1.61377382316176e-10 6.84326027598147e-12 -4.66898885683671e-12 2.31211355085535e-12 2.39195112937346e-12;2.99634365075821e-07 8.14391615472128e-06 6.70458490942443e-06 -9.92542646762000e-07 -3.04078064992750e-06;-6.52697933801393e-07 2.87255329776428e-07 -1.78227609772085e-08 2.65525429849935e-07 8.60650570551813e-08;-1.62727164011710e-07 1.09102479325892e-07 4.97827431850001e-09 7.86649963082937e-11 -6.67193813407656e-09;-2.96370000987760e-09 1.20008401576557e-09 1.75885448022883e-09 -1.74756709684384e-09 3.21963061454248e-09;-9.91101697778560e-10 7.54541713140752e-10 -2.95880967800875e-10 1.81009160501278e-10 8.31547411640954e-11;1.21268051949609e-10 -5.93572774509587e-11 -5.03295034994351e-11 3.05383430975252e-11 3.56280438509939e-11;6.92012970333794e-11 -9.02885345797597e-12 -3.44151832744880e-12 2.03164894681921e-12 -5.44852265137606e-12;5.56731263672800e-12 3.57272150106101e-12 2.25885622368678e-12 -2.44508240047675e-13 -6.83314378535235e-13;3.96883487797254e-06 -4.57100506169608e-06 -3.30208117813256e-06 3.32599719134845e-06 4.26539325549339e-06;1.10123151770973e-06 4.58046760144882e-07 1.86831972581926e-07 -1.60092770735081e-07 -5.58956114867062e-07;-3.40344900506653e-08 2.87649741373047e-08 -1.83929753066251e-08 -9.74179203885847e-09 -2.42064137485043e-09;-6.49731596932566e-09 -3.07048108404447e-09 -2.84380614669848e-09 1.55123146524283e-09 4.53694984588346e-10;5.45175793803325e-10 -3.73287624700125e-10 -1.16293122618336e-10 7.25845618602690e-11 -4.34112440021627e-11;1.89481447552805e-10 3.67431482211078e-12 -1.72180065021194e-11 1.47046319023226e-11 1.31920481414062e-11;2.10125915737167e-12 -3.08420783495975e-12 -4.87748712363020e-12 1.16363599902490e-14 1.26698255558605e-13;-8.07894928696254e-12 9.19344620512607e-13 3.26929173307443e-13 2.00438149416495e-13 -9.57035765212079e-15;1.38737151773284e-12 1.09340178371420e-13 5.15714202449053e-14 -5.92156438588931e-14 -3.29586752336143e-14;6.38137197198254e-06 4.62426300749908e-06 4.42334454191034e-06 1.15374736092349e-06 -2.61859702227253e-06;-2.25320619636149e-07 3.21907705479353e-07 -3.34834530764823e-07 -4.82132753601810e-07 -3.22410936343355e-07;3.48894515496995e-09 3.49951261408458e-08 -6.01128959281142e-09 4.78213900943443e-09 1.46012816168576e-08;-9.66682871952083e-11 3.75806627535317e-09 2.38984004956705e-09 2.07545049877203e-09 1.58573595632766e-09;1.06834370693917e-09 -4.07975055112153e-10 -2.37598937943957e-10 5.89327007480137e-11 1.18891820437634e-10;5.22433722695807e-11 6.02011995016293e-12 -7.80605402956048e-12 1.50873145627341e-11 -1.40550093106311e-12;2.13396242187279e-13 -1.71939313965536e-12 -3.57625378660975e-14 -5.01675184988446e-14 -1.07805487368797e-12;-1.24352330043311e-12 8.26105883301606e-13 4.63606970128517e-13 6.39517888984486e-14 -7.35135439920086e-14;-5.39023859065631e-13 2.54188315588243e-14 1.30933833278664e-14 6.06153473304781e-15 -4.24722717533726e-14;3.12767756884813e-14 -2.29517847871632e-15 2.53117304424948e-16 7.07504914138118e-16 -1.20089065310688e-15;2.08311178819214e-06 -1.22179185044174e-06 -2.98842190131044e-06 3.07310218974299e-06 2.27100346036619e-06;-3.94601643855452e-07 -5.44014825116083e-07 -6.16955333162507e-08 -2.31954821580670e-07 1.14010813005310e-07;6.11067575043044e-08 -3.93240193194272e-08 -1.62979132528933e-08 1.01339204652581e-08 1.97319601566071e-08;2.57770508710055e-09 1.87799543582899e-09 1.95407654714372e-09 1.15276419281270e-09 2.25397005402120e-09;7.16926338026236e-10 -3.65857693313858e-10 -1.54864067050915e-11 6.50770211276549e-11 -7.85160007413546e-12;4.90007693914221e-12 3.31649396536340e-12 4.81664871165640e-13 7.26080745617085e-12 2.30960953372164e-12;9.75489202240545e-13 -1.68967954531421e-13 7.38383391334110e-13 -3.58435515913239e-13 -3.01564710027450e-13;-3.79533601922805e-13 2.76681830946617e-13 1.21480375553803e-13 -1.57729077644850e-14 -8.87664977818700e-14;-3.96462845480288e-14 2.94155690934610e-14 6.78413205760717e-15 -4.12135802787361e-15 -1.46373307795619e-14;-8.64941937408121e-15 -1.91822620970386e-15 -8.01725413560744e-16 5.02941051180784e-16 -1.07572628474344e-15;-4.13816294742758e-15 -7.43602019785880e-17 -5.54248556346072e-17 -4.83999456005158e-17 -1.19622559730466e-16;-8.34852132750364e-07 -7.45794677612056e-06 -6.58132648865533e-06 -1.38608110346732e-06 5.32326534882584e-07;-2.75513802414150e-07 3.64713745106279e-08 -7.12385417940442e-08 -7.86206067228882e-08 2.28048393207161e-08;-4.26696415431918e-08 -4.65599668635087e-09 7.35037936327566e-09 1.17098354115804e-08 1.44594777658035e-08;1.12407689274199e-09 7.62142529563709e-10 -6.72563708415472e-10 -1.18094592485992e-10 -1.17043815733292e-09;1.76612225246125e-10 -1.01188552503192e-10 7.32546072616968e-11 1.79542821801610e-11 -2.23264859965402e-11;-9.35960722512375e-12 1.90894283812231e-12 -6.34792824525760e-13 3.98597963877826e-12 -4.47591409078971e-12;-3.34623858556099e-12 4.56384903915853e-14 2.72561108521416e-13 -3.57942733300468e-15 1.99794810657713e-13;-6.16775522568954e-14 8.25316968328823e-14 7.19845814260518e-14 -2.92415710855106e-14 -5.49570017444031e-15;-8.50728802453217e-15 8.38161600916267e-15 3.43651657459983e-15 -8.19429434115910e-16 -4.08905746461100e-15;4.39042894275548e-15 -3.69440485320477e-16 1.22249256876779e-16 -2.09359444520984e-16 -3.34211740264257e-16;-5.36054548134225e-16 3.29794204041989e-17 2.13564354374585e-17 -1.37838993720865e-18 -1.29188342867753e-17;-3.26421841529845e-17 7.38235405234126e-18 2.49291659676210e-18 8.18252735459593e-19 1.73824952279230e-20;4.67237509268208e-06 1.93611283787239e-06 9.39035455627622e-07 -5.84565118072823e-07 -1.76198705802101e-07;-3.33739157421993e-07 4.12139555299163e-07 1.58754695700856e-07 1.37448753329669e-07 1.04722936936873e-07;6.64200603076386e-09 1.45412222625734e-08 1.82498796118030e-08 2.86633517581614e-09 1.06066984548100e-09;5.25549696746655e-09 -1.33677183394083e-09 7.60804375937931e-11 -1.07918624219037e-10 8.09178898247941e-10;1.89318454110039e-10 9.23092164791765e-11 5.51434573131180e-11 3.86696392289240e-11 -1.15208165047149e-11;-1.02252706006226e-12 -7.25921015411136e-13 -1.98110126887620e-12 -2.18964868282672e-13 -7.18834476685625e-13;-2.69770025318548e-12 -2.17850340796321e-14 4.73040820865871e-13 1.57947421572149e-13 1.86925164972766e-13;1.07831718354771e-13 2.26681841611017e-14 2.56046087047783e-14 -1.14995851659554e-14 -2.27056907624485e-14;6.29825154734712e-15 8.04458225889001e-16 9.53173540411138e-16 1.16892301877735e-15 -1.04324684545047e-15;-5.57345639727027e-16 -2.93949227634932e-16 7.47621406284534e-18 -5.36416885470756e-17 -2.87213280230513e-16;1.73219775047208e-16 2.05017387523061e-17 9.08873886345587e-18 -2.86881547225742e-18 -1.25303645304992e-17;-7.30829109684568e-18 2.03711261415353e-18 7.62162636124024e-19 -7.54847922012517e-19 -8.85105098195030e-19;5.62039968280587e-18 -1.38144206573507e-19 1.68028711767211e-20 1.81223858251981e-19 -8.50245194985878e-20];
            anm_bw = [0.00136127467401223 -6.83476317823061e-07 -1.37211986707674e-06 7.02561866200582e-07 -2.16342338010651e-07;-9.53197486400299e-06 6.58703762338336e-06 2.42000663952044e-06 -6.04283463108935e-07 2.02144424676990e-07;-6.76728911259359e-06 6.03830755085583e-07 -8.72568628835897e-08 2.21750344140938e-06 1.05146032931020e-06;-3.21102832397338e-05 -7.88685357568093e-06 -2.55495673641049e-06 -1.99601934456719e-06 -4.62005252198027e-07;-7.84639263523250e-07 3.11624739733849e-06 9.02170019697389e-07 6.37066632506008e-07 -9.44485038780872e-09;2.19476873575507e-06 -2.20580510638233e-07 6.94761415598378e-07 4.80770865279717e-07 -1.34357837196401e-07;2.18469215148328e-05 -1.80674174262038e-06 -1.52754285605060e-06 -3.51212288219241e-07 2.73741237656351e-06;2.85579058479116e-06 1.57201369332361e-07 -2.80599072875081e-07 -4.91267304946072e-07 -2.11648188821805e-07;2.81729255594770e-06 3.02487362536122e-07 -1.64836481475431e-07 -2.11607615408593e-07 -6.47817762225366e-08;1.31809947620223e-07 -1.58289524114549e-07 -7.05580919885505e-08 5.56781440550867e-08 1.23403290710365e-08;-1.29252282695869e-05 -1.07247072037590e-05 -3.31109519638196e-06 2.13776673779736e-06 -1.49519398373391e-07;1.81685152305722e-06 -1.17362204417861e-06 -3.19205277136370e-08 4.09166457255416e-07 1.53286667406152e-07;1.63477723125362e-06 -2.68584775517243e-08 4.94662064805191e-09 -7.09027987928288e-08 4.44353430574937e-08;-2.13090618917978e-07 4.05836983493219e-08 2.94495876336549e-08 -1.75005469063176e-08 -3.03015988647002e-09;-2.16074435298006e-09 9.37631708987675e-09 -2.05996036369828e-08 6.97068002894092e-09 -8.90988987979604e-09;1.38047798906967e-05 2.05528261553901e-05 1.59072148872708e-05 7.34088731264443e-07 1.28226710383580e-06;7.08175753966264e-07 -9.27988276636505e-07 1.60535820026081e-07 -3.27296675122065e-07 -2.20518321170684e-07;1.90932483086199e-07 -7.44215272759193e-08 1.81330673333187e-08 4.37149649043616e-08 4.18884335594172e-08;-5.37009063880924e-08 2.22870057779431e-08 1.73740123037651e-08 -4.45137302235032e-09 9.44721910524571e-09;-6.83406949047909e-08 -1.95046676795923e-10 2.57535903049686e-09 4.82643164083020e-09 3.37657333705158e-09;3.96128688448981e-09 -6.63809403270686e-10 2.44781464212534e-10 5.92280853590699e-11 -4.78502591970721e-10;1.75859399041414e-05 -2.81238050668481e-06 -2.43670534594848e-06 3.58244562699714e-06 -1.76547446732691e-06;-1.06451311473304e-07 1.54336689617184e-06 -2.00690000442673e-07 1.38790047911880e-09 -1.62490619890017e-07;-2.72757421686155e-07 1.71139266205398e-07 -2.55080309401917e-08 -8.40793079489831e-09 -1.01129447760167e-08;2.92966025844079e-08 -2.07556718857313e-08 5.45985315647905e-09 8.76857690274150e-09 1.06785510440474e-08;-1.22059608941331e-08 6.52491630264276e-09 -1.79332492326928e-10 3.75921793745396e-10 -7.06416506254786e-10;1.63224355776652e-09 4.95586028736232e-10 -3.07879011759040e-10 -7.78354087544277e-11 1.43959047067250e-10;3.86319414653663e-10 -2.06467134617933e-10 4.37330971382694e-11 -5.00421056263711e-11 -9.40237773015723e-12;-1.23856142706451e-05 7.61047394008415e-06 -1.99104114578138e-07 6.86177748886858e-07 -1.09466747592827e-07;2.99866062403128e-07 1.87525561397390e-07 4.99374806994715e-08 4.86229763781404e-07 4.46570575517658e-07;-5.05748332368430e-07 1.95523624722285e-08 -9.17535435911345e-08 -2.56671607433547e-08 -7.11896201616653e-08;-2.66062200406494e-08 -5.40470019739274e-09 -2.29718660244954e-09 -3.73328592264404e-09 3.38748313712376e-09;5.30855327954894e-10 5.28851845648032e-10 -2.22278913745418e-10 -5.52628653064771e-11 -9.24825145219684e-10;6.03737227573716e-10 -3.52190673510919e-12 -1.30371720641414e-10 -9.12787239944822e-12 6.42187285537238e-12;1.78081862458539e-10 2.93772078656037e-12 -1.04698379945322e-11 -2.82260024833024e-11 -5.61810459067525e-12;9.35003092299580e-12 -8.23133834521577e-13 5.54878414224198e-13 -3.62943215777181e-13 2.38858933771653e-12;-1.31216096107331e-05 -5.70451670731759e-06 -5.11598683573971e-06 -4.99990779887599e-06 1.27389320221511e-07;-1.23108260369048e-06 5.53093245213587e-07 8.60093183929302e-07 2.65569700925696e-07 1.95485134805575e-07;-2.29647072638049e-07 -5.45266515081825e-08 2.85298129762263e-08 1.98167939680185e-08 5.52227340898335e-09;-2.73844745019857e-08 -4.48345173291362e-10 -1.93967347049382e-09 -1.41508853776629e-09 -1.75456962391145e-09;-2.68863184376108e-11 -2.20546981683293e-09 6.56116990576877e-10 1.27129855674922e-10 -2.32334506413213e-10;1.98303136881156e-10 6.04782006047075e-11 2.91291115431570e-11 6.18098615782757e-11 -3.82682292530379e-11;9.48294455071158e-12 -3.05873596453015e-13 5.31539408055057e-13 -7.31016438665600e-12 -1.19921002209198e-11;-2.25188050845725e-11 -3.91627574966393e-13 -6.80217235976769e-13 5.91033607278405e-13 5.02991534452191e-13;1.29532063896247e-12 1.66337285851564e-13 3.25543028344555e-13 1.89143357962363e-13 3.32288378169726e-13;-2.45864358781728e-06 4.49460524898260e-06 1.03890496648813e-06 -2.73783420376785e-06 7.12695730642593e-07;-9.27805078535168e-07 -4.97733876686731e-07 9.18680298906510e-08 -2.47200617423980e-07 6.16163630140379e-08;-1.39623661883136e-08 -1.12580495666505e-07 2.61821435950379e-08 -2.31875562002885e-08 5.72679835033659e-08;-9.52538983318497e-09 -5.40909215302433e-09 1.88698793952475e-09 -4.08127746406372e-09 1.09534895853812e-10;3.79767457525741e-09 1.11549801373366e-10 -6.45504957274111e-10 3.05477141010356e-10 1.26261210565856e-10;5.08813577945300e-11 1.43250547678637e-11 8.81616572082448e-12 2.58968878880804e-11 3.83421818249954e-11;8.95094368142044e-12 -3.26220304555971e-12 -1.28047847191896e-12 2.67562170258942e-12 2.72195031576670e-12;-6.47181697409757e-12 1.13776457455685e-12 2.84856274334969e-13 -7.63667272085395e-14 -1.34451657758826e-13;-1.25291265888343e-12 8.63500441050317e-14 -1.21307856635548e-13 5.12570529540511e-14 3.32389276976573e-14;3.73573418085813e-14 -5.37808783042784e-16 -4.23430408270850e-16 -4.75110565740493e-15 6.02553212780166e-15;8.95483987262751e-06 -3.90778212666235e-06 -1.12115019808259e-06 1.78678942093383e-06 1.46806344157962e-06;-4.59185232678613e-07 1.09497995905419e-07 1.31663977640045e-07 4.20525791073626e-08 -9.71470741607431e-08;1.63399802579572e-07 1.50909360648645e-08 -1.11480472593347e-08 -1.84000857674573e-08 7.82124614794256e-09;1.22887452385094e-08 -4.06647399822746e-10 -6.49120327585597e-10 8.63651225791194e-10 -2.73440085913102e-09;2.51748630889583e-09 4.79895880425564e-10 -2.44908073860844e-10 2.56735882664876e-10 -1.64815306286912e-10;4.85671381736718e-11 -2.51742732115131e-11 -2.60819437993179e-11 6.12728324086123e-12 2.16833310896138e-11;4.11389702320298e-12 -8.09433180989935e-13 -1.19812498226024e-12 1.46885737888520e-12 3.15807685137836e-12;-1.47614580597013e-12 4.66726413909320e-13 1.72089709006255e-13 1.13854935381418e-13 2.77741161317003e-13;-1.02257724967727e-13 1.10394382923502e-13 -3.14153505370805e-15 2.41103099110106e-14 2.13853053149771e-14;-3.19080885842786e-14 -9.53904307973447e-15 2.74542788156379e-15 2.33797859107844e-15 -2.53192474907304e-15;-5.87702222126367e-15 -1.80133850930249e-15 -3.09793125614454e-16 -1.04197538975295e-16 3.72781664701327e-16;1.86187054729085e-06 8.33098045333428e-06 3.18277735484232e-06 -7.68273797022231e-07 -1.52337222261696e-06;-5.07076646593648e-07 -8.61959553442156e-07 -3.51690005432816e-07 -4.20797082902431e-07 -3.07652993252673e-07;-7.38992472164147e-08 -8.39473083080280e-08 -2.51587083298935e-08 7.30691259725451e-09 -3.19457155958983e-08;-1.99777182012924e-09 -3.21265085916022e-09 -4.84477421865675e-10 -1.82924814205799e-09 -3.46664344655997e-10;-7.05788559634927e-11 1.21840735569025e-10 7.97347726425926e-11 1.08275679614409e-10 -1.17891254809785e-10;1.10299718947774e-11 -3.22958261390263e-11 -1.43535798209229e-11 6.87096504209595e-12 -6.64963212272352e-12;-6.47393639740084e-12 1.03156978325120e-12 -9.20099775082358e-14 -2.40150316641949e-13 1.14008812047857e-12;-1.23957846397250e-13 2.85996703969692e-13 1.91579874982553e-13 5.20597174693064e-14 -4.06741434883370e-14;-2.35479068911236e-14 1.97847338186993e-14 1.58935977518516e-15 -2.32217195254742e-15 -8.48611789490575e-15;1.03992320391626e-14 1.54017082092642e-15 1.05950035082788e-16 -1.17870898461353e-15 -1.10937420707372e-15;-1.09011948374520e-15 -6.04168007633584e-16 -9.10901998157436e-17 1.98379116989461e-16 -1.03715496658498e-16;-1.38171942108278e-16 -6.33037999097522e-17 -1.38777695011470e-17 1.94191397045401e-17 5.70055906754485e-18;1.92989406002085e-06 -3.82662130483128e-06 -4.60189561036048e-07 2.24290587856309e-06 1.40544379451550e-06;6.49033717633394e-08 2.41396114435326e-07 2.73948898223321e-07 1.10633664439332e-07 -3.19555270171075e-08;-2.91988966963297e-08 -6.03828192816571e-09 1.18462386444840e-08 1.32095545004128e-08 -5.06572721528914e-09;7.31079058474148e-09 -8.42775299751834e-10 1.10190810090667e-09 1.96592273424306e-09 -2.13135932785688e-09;7.06656405314388e-11 1.43441125783756e-10 1.46962246686924e-10 7.44592776425197e-11 -3.64331892799173e-11;-2.52393942119372e-11 1.07520964869263e-11 5.84669886072094e-12 6.52029744217103e-12 1.82947123132059e-12;-4.15669940115121e-12 -1.95963254053648e-13 2.16977822834301e-13 -2.84701408462031e-13 4.27194601040231e-13;3.07891105454129e-13 1.91523190672955e-13 1.05367297580989e-13 -5.28136363920236e-14 -3.53364110005917e-14;7.02156663274738e-15 9.52230536780849e-15 -3.41019408682733e-15 -3.59825303352899e-15 -2.62576411636150e-15;-1.75110277413804e-15 5.29265220719483e-16 4.45015980897919e-16 -3.80179856341347e-16 -4.32917763829695e-16;1.16038609651443e-16 -6.69643574373352e-17 2.65667154817303e-17 -9.76010333683956e-17 4.07312981076655e-17;5.72659246346386e-18 1.30357528108671e-18 2.49193258417535e-18 1.76247014075584e-18 7.59614374197688e-19;1.03352170833303e-17 -2.30633516638829e-18 2.84777940620193e-18 -7.72161347944693e-19 6.07028034506380e-19];
            anm_ch = [0.0571481238161787 3.35402081801137e-05 3.15988141788728e-05 -1.34477341887086e-05 -2.61831023577773e-07;5.77367395845715e-05 -0.000669057185209558 -6.51057691648904e-05 -1.61830149147091e-06 8.96771209464758e-05;-8.50773002452907e-05 -4.87106614880272e-05 4.03431160775277e-05 2.54090162741464e-06 -5.59109319864264e-06;0.00150536423187709 0.000611682258892697 0.000369730024614855 -1.95658439780282e-05 -3.46246726553700e-05;-2.32168718433966e-05 -0.000127478686553809 -9.00292451740728e-05 -6.07834315901830e-05 -1.04628419422714e-05;-1.38607250922551e-06 -3.97271603842309e-06 -8.16155320152118e-07 5.73266706046665e-07 2.00366060212696e-07;6.52491559188663e-05 -0.00112224323460183 -0.000344967958304075 -7.67282640947300e-05 0.000107907110551939;-0.000138870461448036 -7.29995695401936e-05 5.35986591445824e-05 9.03804869703890e-06 8.61370129482732e-06;-9.98524443968768e-07 -6.84966792665998e-08 1.47478021860771e-07 1.94857794008064e-06 7.17176852732910e-07;1.27066367911720e-06 1.12113289164288e-06 2.71525688515375e-07 -2.76125723009239e-07 -1.05429690305013e-07;-0.000377264999981652 0.000262691217024294 0.000183639785837590 3.93177048515576e-06 -6.66187081899168e-06;-4.93720951871921e-05 -0.000102820030405771 -5.69904376301748e-05 -3.79603438055116e-05 -3.96726017834930e-06;-2.21881958961135e-06 -1.40207117987894e-06 1.60956630798516e-07 2.06121145135022e-06 6.50944708093149e-07;2.21876332411271e-07 1.92272880430386e-07 -6.44016558013941e-09 -1.40954921332410e-07 -4.26742169137667e-07;-3.51738525149881e-08 2.89616194332516e-08 -3.40343352397886e-08 -2.89763392721812e-08 -6.40980581663785e-10;3.51240856823468e-05 -0.000725895015345786 -0.000322514037108045 -0.000106143759981636 4.08153152459337e-05;-2.36269716929413e-05 -4.20691836557932e-05 1.43926743222922e-05 2.61811210631784e-05 2.09610762194903e-05;-7.91765756673890e-07 1.64556789159745e-06 -9.43930166276555e-07 6.46641738736139e-07 -5.91509547299176e-07;3.92768838766879e-07 -1.98027731703690e-07 -5.41303590057253e-08 -4.21705797874207e-07 -6.06042329660681e-08;-1.56650141024305e-08 7.61808165752027e-08 -1.81900460250934e-08 1.30196216971675e-08 1.08616031342379e-08;-2.80964779829242e-08 -7.25951488826103e-09 -2.59789823306225e-09 -2.79271942407154e-09 4.10558774868586e-09;-0.000638227857648286 -0.000154814045363391 7.78518327501759e-05 -2.95961469342381e-05 1.15965225055757e-06;4.47833146915112e-06 1.33712284237555e-05 3.61048816552123e-06 -2.50717844073547e-06 -1.28100822021734e-05;-2.26958070007455e-06 2.57779960912242e-06 1.08395653197976e-06 1.29403393862805e-07 -1.04854652812567e-06;-3.98954043463392e-07 -2.26931182815454e-07 -1.09169545045028e-07 -1.49509536031939e-07 -3.98376793949903e-07;2.30418911071110e-08 1.23098508481555e-08 -1.71161401463708e-08 2.35829696577657e-09 1.31136164162040e-08;3.69423793101582e-09 3.49231027561927e-10 -1.18581468768647e-09 5.43180735828820e-10 5.43192337651588e-10;-1.38608847117992e-09 -1.86719145546559e-10 -8.13477384765498e-10 2.01919878240491e-10 1.00067892622287e-10;-4.35499078415956e-05 0.000450727967957804 0.000328978494268850 -3.05249478582848e-05 -3.21914834544310e-05;1.24887940973241e-05 1.34275239548403e-05 1.11275518344713e-06 7.46733554562851e-06 -2.12458664760353e-06;9.50250784948476e-07 2.34367372695203e-06 -5.43099244798980e-07 -4.35196904508734e-07 -8.31852234345897e-07;5.91775478636535e-09 -1.48970922508592e-07 2.99840061173840e-08 -1.30595933407792e-07 1.27136765045597e-07;-1.78491083554475e-08 1.76864919393085e-08 -1.96740493482011e-08 1.21096708004261e-08 2.95518703155064e-10;1.75053510088658e-09 -1.31414287871615e-09 -1.44689439791928e-09 1.14682483668460e-09 1.74488616540169e-09;1.08152964586251e-09 -3.85678162063266e-10 -2.77851016629979e-10 3.89890578625590e-11 -2.54627365853495e-10;-1.88340955578221e-10 5.19645384002867e-11 2.14131326027631e-11 1.24027770392728e-11 -9.42818962431967e-12;0.000359777729843898 -0.000111692619996219 -6.87103418744904e-05 0.000115128973879551 7.59796247722486e-05;5.23717968000879e-05 1.32279078116467e-05 -5.72277317139479e-07 -7.56326558610214e-06 -1.95749622214651e-05;1.00109213210139e-06 -2.75515216592735e-07 -1.13393194050846e-06 -4.75049734870663e-07 -3.21499480530932e-07;-2.07013716598890e-07 -7.31392258077707e-08 -3.96445714084160e-08 3.21390452929387e-08 -1.43738764991525e-08;2.03081434931767e-09 -1.35423687136122e-08 -4.47637454261816e-09 2.18409121726643e-09 -3.74845286805217e-09;3.17469255318367e-09 2.44221027314129e-10 -2.46820614760019e-10 7.55851003884434e-10 6.98980592550891e-10;9.89541493531067e-11 -2.78762878057315e-11 -2.10947962916771e-10 3.77882267360636e-11 -1.20009542671532e-12;5.01720575730940e-11 1.66470417102135e-11 -7.50624817938091e-12 9.97880221482238e-12 4.87141864438892e-12;2.53137945301589e-11 1.93030083090772e-12 -1.44708804231290e-12 -1.77837100743423e-12 -8.10068935490951e-13;0.000115735341520738 0.000116910591048350 8.36315620479475e-05 1.61095702669207e-05 -7.53084853489862e-05;-9.76879433427199e-06 9.16968438003335e-06 -8.72755127288830e-06 -1.30077933880053e-05 -9.78841937993320e-06;1.04902782517565e-07 2.14036988364936e-07 -7.19358686652888e-07 1.12529592946332e-07 7.07316352860448e-07;7.63177265285080e-08 1.22781974434290e-07 8.99971272969286e-08 5.63482239352990e-08 4.31054352285547e-08;3.29855763107355e-09 -6.95004336734441e-09 -6.52491370576354e-09 1.97749180391742e-09 3.51941791940498e-09;3.85373745846559e-10 1.65754130924183e-10 -3.31326088103057e-10 5.93256024580436e-10 1.27725220636915e-10;-1.08840956376565e-10 -4.56042860268189e-11 -4.77254322645633e-12 -2.94405398621875e-12 -3.07199979999475e-11;2.07389879095010e-11 1.51186798732451e-11 9.28139802941848e-12 5.92738269687687e-12 9.70337402306505e-13;-2.85879708060306e-12 1.92164314717053e-13 4.02664678967890e-14 5.18246319204277e-13 -7.91438726419423e-13;6.91890667590734e-13 -8.49442290988352e-14 -5.54404947212402e-15 9.71093377538790e-15 -5.33714333415971e-14;-5.06132972789792e-05 -4.28348772058883e-05 -6.90746551020305e-05 8.48380415176836e-05 7.04135614675053e-05;-1.27945598849788e-05 -1.92362865537803e-05 -2.30971771867138e-06 -8.98515975724166e-06 5.25675205004752e-06;-8.71907027470177e-07 -1.02091512861164e-06 -1.69548051683864e-07 4.87239045855761e-07 9.13163249899837e-07;-6.23651943425918e-08 6.98993315829649e-08 5.91597766733390e-08 4.36227124230661e-08 6.45321798431575e-08;-1.46315079552637e-10 -7.85142670184337e-09 1.48788168857903e-09 2.16870499912160e-09 -1.16723047065545e-09;3.31888494450352e-10 1.90931898336457e-10 -3.13671901557599e-11 2.60711798190524e-10 8.45240112207997e-11;1.36645682588537e-11 -5.68830303783976e-12 1.57518923848140e-11 -1.61935794656758e-11 -4.16568077748351e-12;9.44684950971905e-13 7.30313977131995e-12 3.14451447892684e-12 6.49029875639842e-13 -9.66911019905919e-13;-8.13097374090024e-13 5.23351897822186e-13 8.94349188113951e-14 -1.33327759673270e-13 -4.04549450989029e-13;-3.76176467005839e-14 -6.19953702289713e-14 -3.74537190139726e-14 1.71275486301958e-14 -3.81946773167132e-14;-4.81393385544160e-14 3.66084990006325e-15 3.10432030972253e-15 -4.10964475657416e-15 -6.58644244242900e-15;-7.81077363746945e-05 -0.000254773632197303 -0.000214538508009518 -3.80780934346726e-05 1.83495359193990e-05;5.89140224113144e-06 -3.17312632433258e-06 -3.81872516710791e-06 -2.27592226861647e-06 1.57044619888023e-06;-1.44272505088690e-06 -1.10236588903758e-07 2.64336813084693e-07 4.76074163332460e-07 4.28623587694570e-07;3.98889120733904e-08 -1.29638005554027e-08 -4.13668481273828e-08 1.27686793719542e-09 -3.54202962042383e-08;1.60726837551750e-09 -2.70750776726156e-09 2.79387092681070e-09 -3.01419734793998e-10 -1.29101669438296e-10;-2.55708290234943e-10 2.27878015173471e-11 -6.43063443462716e-12 1.26531554846856e-10 -1.65822147437220e-10;-3.35886470557484e-11 -3.51895009091595e-12 5.80698399963198e-12 -2.84881487149207e-12 8.91708061745902e-12;-3.12788523950588e-12 3.35366912964637e-12 2.52236848033838e-12 -8.12801050709184e-13 -2.63510394773892e-13;6.83791881183142e-14 2.41583263270381e-13 8.58807794189356e-14 -5.12528492761045e-14 -1.40961725631276e-13;-1.28585349115321e-14 -2.11049721804969e-14 5.26409596614749e-15 -4.31736582588616e-15 -1.60991602619068e-14;-9.35623261461309e-15 -3.94384886372442e-16 5.04633016896942e-16 -5.40268998456055e-16 -1.07857944298104e-15;8.79756791888023e-16 4.52529935675330e-16 1.36886341163227e-16 -1.12984402980452e-16 6.30354561057224e-18;0.000117829256884757 2.67013591698442e-05 2.57913446775250e-05 -4.40766244878807e-05 -1.60651761172523e-06;-1.87058092029105e-05 1.34371169060024e-05 5.59131416451555e-06 4.50960364635647e-06 2.87612873904633e-06;2.79835536517287e-07 8.93092708148293e-07 8.37294601021795e-07 -1.99029785860896e-08 -8.87240405168977e-08;4.95854313394905e-08 -1.44694570735912e-08 2.51662229339375e-08 -3.87086600452258e-09 2.29741919071270e-08;4.71497840986162e-09 2.47509999454076e-09 1.67323845102824e-09 8.14196768283530e-10 -3.71467396944165e-10;-1.07340743907054e-10 -8.07691657949326e-11 -5.99381660248133e-11 2.33173929639378e-12 -2.26994195544563e-11;-3.83130441984224e-11 -5.82499946138714e-12 1.43286311435124e-11 3.15150503353387e-12 5.97891025146774e-12;-5.64389191072230e-13 9.57258316335954e-13 1.12055192185939e-12 -4.42417706775420e-13 -9.93190361616481e-13;1.78188860269677e-13 7.82582024904950e-14 5.18061650118009e-14 2.13456507353387e-14 -5.26202113779510e-14;-8.18481324740893e-15 -3.71256746886786e-15 4.23508855164371e-16 -2.91292502923102e-15 -1.15454205389350e-14;6.16578691696810e-15 6.74087154080877e-16 5.71628946437034e-16 -2.05251213979975e-16 -7.25999138903781e-16;9.35481959699383e-17 6.23535830498083e-17 3.18076728802060e-18 -2.92353209354587e-17 7.65216088665263e-19;2.34173078531701e-17 -8.30342420281772e-18 -4.33602329912952e-18 1.90226281379981e-18 -7.85507922718903e-19];
            anm_cw = [0.0395329695826997 -0.000131114380761895 -0.000116331009006233 6.23548420410646e-05 5.72641113425116e-05;-0.000441837640880650 0.000701288648654908 0.000338489802858270 3.76700309908602e-05 -8.70889013574699e-06;1.30418530496887e-05 -0.000185046547597376 4.31032103066723e-05 0.000105583334124319 3.23045436993589e-05;3.68918433448519e-05 -0.000219433014681503 3.46768613485000e-06 -9.17185187163528e-05 -3.69243242456081e-05;-6.50227201116778e-06 2.07614874282187e-05 -5.09131314798362e-05 -3.08053225174359e-05 -4.18483655873918e-05;2.67879176459056e-05 -6.89303730743691e-05 2.11046783217168e-06 1.93163912538178e-05 -1.97877143887704e-06;0.000393937595007422 -0.000452948381236406 -0.000136517846073846 0.000138239247989489 0.000133175232977863;5.00214539435002e-05 3.57229726719727e-05 -9.38010547535432e-07 -3.52586798317563e-05 -7.01218677681254e-06;3.91965314099929e-05 1.02236686806489e-05 -1.95710695226022e-05 -5.93904795230695e-06 3.24339769876093e-06;6.68158778290653e-06 -8.10468752307024e-06 -9.91192994096109e-06 -1.89755520007723e-07 -3.26799467595579e-06;0.000314196817753895 -0.000296548447162009 -0.000218410153263575 -1.57318389871000e-05 4.69789570185785e-05;0.000104597721123977 -3.31000119089319e-05 5.60326793626348e-05 4.71895007710715e-05 3.57432326236664e-05;8.95483021572039e-06 1.44019305383365e-05 4.87912790492931e-06 -3.45826387853503e-06 3.23960320438157e-06;-1.35249651009930e-05 -2.49349762695977e-06 -2.51509483521132e-06 -9.14254874104858e-07 -8.57897406100890e-07;-1.68143325235195e-06 1.72073417594235e-06 1.38765993969565e-06 4.09770982137530e-07 -6.60908742097123e-07;-0.000639889366487161 0.00120194042474696 0.000753258598887703 3.87356377414663e-05 1.31231811175345e-05;2.77062763606783e-05 -9.51425270178477e-06 -6.61068056107547e-06 -1.38713669012109e-05 9.84662092961671e-06;-2.69398078539471e-06 6.50860676783123e-06 3.80855926988090e-06 -1.98076068364785e-06 1.17187335666772e-06;-2.63719028151905e-06 5.03149473656743e-07 7.38964893399716e-07 -8.38892485369078e-07 1.30943917775613e-06;-1.56634992245479e-06 -2.97026487417045e-08 5.06602801102463e-08 -4.60436007958792e-08 -1.62536449440997e-07;-2.37493912770935e-07 1.69781593069938e-08 8.35178275224265e-08 -4.83564044549811e-08 -4.96448864199318e-08;0.00134012259587597 -0.000250989369253194 -2.97647945512547e-05 -6.47889968094926e-05 8.41302130716859e-05;-0.000113287184900929 4.78918993866293e-05 -3.14572113583139e-05 -2.10518256626847e-05 -2.03933633847417e-05;-4.97413321312139e-07 3.72599822034753e-06 -3.53221588399266e-06 -1.05232048036416e-06 -2.74821498198519e-06;4.81988542428155e-06 4.21400219782474e-07 1.02814808667637e-06 4.40299068486188e-09 3.37103399036634e-09;1.10140301678818e-08 1.90257670180182e-07 -1.00831353341885e-08 1.44860642389714e-08 -5.29882089987747e-08;6.12420414245775e-08 -4.48953461152996e-09 -1.38837603709003e-08 -2.05533675904779e-08 1.49517908802329e-09;9.17090243673643e-10 -9.24878857867367e-09 -2.30856560363943e-09 -4.36348789716735e-09 -4.45808881183025e-10;-0.000424912699609112 -0.000114365438471564 -0.000403200981827193 4.19949560550194e-05 -3.02068483713739e-05;3.85435472851225e-05 -5.70726887668306e-05 4.96313706308613e-07 1.02395703617082e-05 5.85550000567006e-06;-7.38204470183331e-06 -4.56638770109511e-06 -3.94007992121367e-06 -2.16666812189101e-06 -4.55694264113194e-06;5.89841165408527e-07 1.40862905173449e-08 1.08149086563211e-07 -2.18592601537944e-07 -3.78927431428119e-07;4.85164687450468e-08 8.34273921293655e-08 1.47489605513673e-08 6.01494125001291e-08 6.43812884159484e-09;1.13055580655363e-08 3.50568765400469e-09 -5.09396162501750e-09 -1.83362063152411e-09 -4.11227251553035e-09;3.16454132867156e-09 -1.39634794131087e-09 -7.34085003895929e-10 -7.55541371271796e-10 -1.57568747643705e-10;1.27572900992112e-09 -3.51625955080441e-10 -4.84132020565098e-10 1.52427274930711e-10 1.27466120431317e-10;-0.000481655666236529 -0.000245423313903835 -0.000239499902816719 -0.000157132947351028 5.54583099258017e-05;-1.52987254785589e-05 2.78383892116245e-05 4.32299123991860e-05 1.70981319744327e-05 -1.35090841769225e-06;-8.65400907717798e-06 -6.51882656990376e-06 -2.43810171017369e-07 8.54348785752623e-07 2.98371863248143e-07;-1.68155571776752e-06 -3.53602587563318e-07 -1.00404435881759e-07 -2.14162249012859e-08 -2.42131535531526e-07;-1.08048603277187e-08 -9.78850785763030e-08 -2.32906554437417e-08 2.22003630858805e-08 -2.27230368089683e-09;-5.98864391551041e-09 7.38970926486848e-09 3.61322835311957e-09 3.70037329172919e-09 -3.41121137081362e-09;-7.33113754909726e-10 -9.08374249335220e-11 -1.78204392133739e-10 8.28618491929026e-11 -1.32966817912373e-10;-5.23340481314676e-10 1.36403528233346e-10 -7.04478837151279e-11 -6.83175201536443e-12 -2.86040864071134e-12;3.75347503578356e-11 -1.08518134138781e-11 -2.53583751744508e-12 1.00168232812303e-11 1.74929602713312e-11;-0.000686805336370570 0.000591849814585706 0.000475117378328026 -2.59339398048415e-05 3.74825110514968e-05;3.35231363034093e-05 2.38331521146909e-05 7.43545963794093e-06 -3.41430817541849e-06 7.20180957675353e-06;3.60564374432978e-07 -3.13300039589662e-06 -6.38974746108020e-07 -8.63985524672024e-07 2.43367665208655e-06;-4.09605238516094e-07 -2.51158699554904e-07 -1.29359217235188e-07 -2.27744642483133e-07 7.04065989970205e-08;6.74886341820129e-08 -1.02009407061935e-08 -3.30790296448812e-08 1.64959795655031e-08 1.40641779998855e-08;1.31706886235108e-09 -1.06243701278671e-09 -2.85573799673944e-09 3.72566568681289e-09 2.48402582003925e-09;-3.68427463251097e-11 -1.90028122983781e-10 -3.98586561768697e-11 1.14458831693287e-11 -2.27722300377854e-12;-7.90029729611056e-11 3.81213646526419e-11 4.63303426711788e-11 1.52294835905903e-11 -2.99094751490726e-12;-2.36146602045017e-11 1.03852674709985e-11 -4.47242126307100e-12 5.30884113537806e-12 1.68499023262969e-12;-3.30107358134527e-13 -4.73989085379655e-13 5.17199549822684e-13 2.34951744478255e-13 2.05931351608192e-13;0.000430215687511780 -0.000132831373000014 -3.41830835017045e-05 4.70312161436033e-06 -3.84807179340006e-05;1.66861163032403e-05 -8.10092908523550e-06 8.20658107437905e-06 6.12399025026683e-06 -1.85536495631911e-06;1.53552093641337e-06 2.19486495660361e-06 -1.07253805120137e-06 -4.72141767909137e-07 4.00744581573216e-07;2.56647305130757e-07 -8.07492046592274e-08 -2.05858469296168e-07 1.09784168930599e-07 -7.76823030181225e-08;1.77744008115031e-08 1.64134677817420e-08 4.86163044879020e-09 1.13334251800856e-08 -7.17260621115426e-09;1.61133063219326e-09 -1.85414677057024e-09 -2.13798537812651e-09 1.15255123229679e-09 2.24504700129464e-09;1.23344223096739e-10 -1.20385012169848e-10 -2.18038256346433e-12 3.23033120628279e-11 8.01179568213400e-11;-6.55745274387847e-12 1.22127104697198e-11 5.83805016355883e-12 -8.31201582509817e-12 1.90985373872656e-12;-2.89199983667265e-12 5.05962500506667e-12 1.28092925110279e-12 5.60353813743813e-13 1.76753731968770e-12;-1.61678729774956e-13 -3.92206170988615e-13 -9.04941327579237e-14 1.89847694200763e-13 4.10008676756463e-14;-1.16808369005656e-13 -9.97464591430510e-14 7.46366550245722e-15 2.53398578153179e-14 1.06510689748906e-14;-0.000113716921384790 -0.000131902722651488 -0.000162844886485788 7.90171538739454e-06 -0.000178768066961413;-2.13146535366500e-06 -3.57818705543597e-05 -1.50825855069298e-05 -2.17909259570022e-05 -8.19332236308581e-06;-2.88001138617357e-06 -2.09957465440793e-06 6.81466526687552e-08 3.58308906974448e-07 -4.18502067223724e-07;-1.10761444317605e-07 6.91773860777929e-08 8.17125372450372e-08 -2.16476237959181e-08 7.59221970502074e-08;-9.56994224818941e-09 6.64104921728432e-09 6.33077902928348e-09 2.85721181743727e-09 -6.39666681678123e-09;4.62558627839842e-10 -1.69014863754621e-09 -2.80260429599733e-10 4.27558937623863e-11 -1.66926133269027e-10;-7.23385132663753e-11 5.51961193545280e-11 3.04070791942335e-11 3.23227055919062e-12 8.47312431934829e-11;-1.61189613765486e-11 1.66868155925172e-11 1.05370341694715e-11 -4.41495859079592e-12 -2.24939051401750e-12;-8.72229568056267e-13 1.88613726203286e-12 1.21711137534390e-14 -1.13342372297867e-12 -6.87151975256052e-13;7.99311988544090e-15 4.46150979586709e-14 7.50406779454998e-14 -3.20385428942275e-14 -1.26543636054393e-14;4.80503817699514e-14 -3.35545623603729e-14 -1.18546423610485e-14 4.19419209985980e-15 -1.73525614436880e-14;-1.20464898830163e-15 -8.80752065000456e-16 -1.22214298993313e-15 1.69928513019657e-15 1.93593051311405e-16;1.68528879784841e-05 3.57144412031081e-05 -1.65999910125077e-05 5.40370336805755e-05 0.000118138122851376;-3.28151779115881e-05 1.04231790790798e-05 -2.80761862890640e-06 2.98996152515593e-06 -2.67641158709985e-06;-2.08664816151978e-06 -1.64463884697475e-06 6.79099429284834e-08 7.23955842946495e-07 -6.86378427465657e-07;-2.88205823027255e-09 2.38319699493291e-09 1.14169347509045e-07 8.12981074994402e-08 -1.56957943666988e-07;-7.09711403570189e-09 6.29470515502988e-09 3.50833306577579e-09 8.31289199649054e-09 -2.14221463168338e-09;-8.11910123910038e-10 3.34047829618955e-10 3.70619377446490e-10 3.30426088213373e-10 4.86297305597865e-11;1.98628160424161e-11 -4.98557831380098e-12 -5.90523187802174e-12 -1.27027116925122e-12 1.49982368570355e-11;2.62289263262748e-12 3.91242360693861e-12 6.56035499387192e-12 -1.17412941089401e-12 -9.40878197853394e-13;-3.37805010124487e-13 5.39454874299593e-13 -2.41569839991525e-13 -2.41572016820792e-13 -3.01983673057198e-13;-1.85034053857964e-13 4.31132161871815e-14 4.13497222026824e-15 -4.60075514595980e-14 -1.92454846400146e-14;2.96113888929854e-15 -1.11688534391626e-14 3.76275373238932e-15 -3.72593295948136e-15 1.98205490249604e-16;1.40074667864629e-15 -5.15564234798333e-16 3.56287382196512e-16 5.07242777691587e-16 -2.30405782826134e-17;2.96822530176851e-16 -4.77029898301223e-17 1.12782285532775e-16 1.58443229778573e-18 8.22141904662969e-17];
            bnm_bh = [0 0 0 0 0;0 0 0 0 0;-2.29210587053658e-06 -2.33805004374529e-06 -7.49312880102168e-07 -5.12022747852006e-07 5.88926055066172e-07;0 0 0 0 0;-4.63382754843690e-06 -2.23853015662938e-06 8.14830531656518e-07 1.15453269407116e-06 -4.53555450927571e-07;-6.92432096320778e-07 -2.98734455136141e-07 1.48085153955641e-08 1.37881746148773e-07 -6.92492118460215e-09;0 0 0 0 0;-1.91507979850310e-06 -1.83614825459598e-06 -7.46807436870647e-07 -1.28329122348007e-06 5.04937180063059e-07;-8.07527103916713e-07 2.83997840574570e-08 -6.01890498063025e-08 -2.48339507554546e-08 2.46284627824308e-08;-2.82995069303093e-07 1.38818274596408e-09 3.22731214161408e-09 2.87731153972404e-10 1.53895537278496e-08;0 0 0 0 0;-6.68210270956800e-07 -2.19104833297845e-06 1.30116691657253e-07 4.78445730433450e-07 -4.40344300914051e-07;-2.36946755740436e-07 -1.32730991878204e-07 1.83669593693860e-08 7.90218931983569e-08 -4.70161979232584e-08;1.07746083292179e-07 -4.17088637760330e-09 -1.83296035841109e-09 -5.80243971371211e-09 -2.11682361167439e-09;-5.44712355496109e-08 1.89717032256923e-09 2.27327316287804e-10 7.78400728280038e-10 8.82380487618991e-12;0 0 0 0 0;-5.61707049615673e-08 -1.09066447089585e-06 -2.25742250174119e-07 -8.64367795924377e-07 1.06411275240680e-08;2.41782935157918e-08 -3.65762298303819e-08 -6.93420659586875e-08 -3.97316214341991e-08 -2.08767816486390e-08;6.38293030383436e-08 1.11377936334470e-08 6.91424941454782e-09 1.39887159955004e-09 5.25428749022906e-09;1.09291268489958e-08 1.23935926756516e-10 3.92917259954515e-10 -1.79144682483562e-10 -9.11802874917597e-10;-4.40957607823325e-09 1.45751390560667e-10 1.24641258165301e-10 -6.45810339804674e-11 -8.92894658893326e-12;0 0 0 0 0;1.54754294162102e-08 -1.60154742388847e-06 -4.08425188394881e-07 6.18170290113531e-09 -2.58919765162122e-07;1.37130642286873e-08 -6.67813955828458e-08 -7.01410996605609e-09 3.82732572660461e-08 -2.73381870915135e-08;2.19113155379218e-08 4.11027496396868e-09 6.33816020485226e-09 -1.49242411327524e-09 -6.14224941851705e-10;6.26573021218961e-09 5.17137416480052e-10 -3.49784328298676e-10 1.13578756343208e-10 2.80414613398411e-10;1.65048133258794e-11 1.00047239417239e-10 1.05124654878499e-10 -3.03826002621926e-11 4.57155388334682e-11;6.20221691418381e-11 9.75852610098156e-12 -5.46716005756984e-12 1.31643349569537e-11 3.61618775715470e-12;0 0 0 0 0;-1.03938913012708e-06 -1.78417431315664e-07 2.86040141364439e-07 1.83508599345952e-08 -1.34452220464346e-07;-4.36557481393662e-08 7.49780206868834e-09 -8.62829428674082e-09 5.50577793039009e-09 -9.46897502333254e-09;3.43193738406672e-10 1.13545447306468e-08 1.25242388852214e-09 6.03221501959620e-10 1.57172070361180e-09;-4.73307591021391e-10 1.70855824051391e-10 -2.62470421477037e-11 2.04525835988874e-10 -1.17859695928164e-10;-3.36185995299839e-10 3.19243054562183e-11 1.17589412418126e-10 -1.35478747434514e-12 5.11192214558542e-11;3.19640547592136e-11 2.94297823804643e-12 -1.00651526276990e-11 -1.67028733953153e-12 3.03938833625503e-12;1.68928641118173e-11 -7.90032886682002e-13 -1.40899773539137e-12 7.76937592393354e-13 7.32539820298651e-13;0 0 0 0 0;2.32949756055277e-07 1.46237594908093e-07 -1.07770884952484e-07 1.26824870644476e-07 -2.36345735961108e-08;8.89572676497766e-08 7.24810004121931e-08 2.67583556180119e-08 2.48434796111361e-08 -3.55004782858686e-09;-1.00823909773603e-08 8.84433929029076e-10 -2.55502517594511e-10 -5.48034274059119e-10 -8.50241938494079e-10;1.13259819566467e-09 5.55186945221216e-10 7.63679807785295e-11 -1.70067998092043e-11 1.57081965572493e-10;-2.37748192185353e-10 2.45463764948000e-11 3.23208414802860e-11 -2.72624834520723e-12 8.14449183666500e-12;-1.54977633126025e-11 4.58754903157884e-12 -1.25864665839074e-12 2.44139868157872e-12 -1.82827441958193e-12;3.28285563794513e-12 -1.10072329225465e-12 -7.23470501810935e-13 5.85309745620389e-13 4.11317589687125e-13;4.57596974384170e-13 9.84198128213558e-14 3.34503817702830e-14 7.08431086558307e-15 2.79891177268807e-14;0 0 0 0 0;-3.67820719155580e-07 6.98497901205902e-07 1.83397388750300e-07 2.39730262495372e-07 -2.58441984368194e-07;5.17793954077994e-08 5.54614175977835e-08 1.75026214305232e-09 -2.55518450411346e-09 -6.12272723006537e-09;-7.94292648157198e-09 -1.01709107852895e-09 -1.49251241812310e-09 9.32827213605682e-10 -8.24490722043118e-10;1.36410408475679e-11 2.16390220454971e-10 1.24934806872235e-10 -6.82507825145903e-11 -4.01575177719668e-11;-1.41619917600555e-11 -1.54733230409082e-11 1.36792829351538e-11 1.11157862104733e-12 2.08548465892268e-11;-3.56521723755846e-12 4.47877185884557e-12 -6.34096209274637e-16 -1.13010624512348e-12 -2.82018136861041e-13;2.22758955943441e-12 -4.63876465559380e-13 -5.80688019272507e-13 2.45878690598655e-13 1.49997666808106e-13;-6.26833903786958e-14 2.73416335780807e-14 1.91842340758425e-14 1.67405061129010e-14 -2.45268543953704e-17;1.81972870222228e-14 5.43036245069085e-15 1.92476637107321e-15 8.78498602508626e-17 -1.42581647227657e-15;0 0 0 0 0;9.74322164613392e-07 -5.23101820582724e-07 -2.81997898176227e-07 4.54762451707384e-08 -3.34645078118827e-08;-6.75813194549663e-09 3.49744702199583e-08 -5.09170419895883e-09 5.24359476874755e-09 4.96664262534662e-09;4.53858847892396e-10 -1.49347392165963e-09 -2.00939511362154e-09 9.30987163387955e-10 9.74450200826854e-11;-4.92900885858693e-10 5.34223033225688e-12 1.08501839729368e-10 -6.43526142089173e-11 -3.11063319142619e-11;1.38469246386690e-11 -7.91180584906922e-12 2.26641656746936e-13 4.55251515177956e-12 6.05270575117769e-12;4.02247935664225e-12 1.82776657951829e-12 -1.28348801405445e-13 -2.16257301300350e-13 -5.54363979435025e-14;4.15005914461687e-13 -2.00647573581168e-13 -1.67278251942946e-13 1.30332398257985e-13 1.52742363652434e-13;6.36376500056974e-14 1.65794532815776e-14 -3.80832559052662e-15 -6.40262894005341e-16 2.42577181848072e-15;-5.55273521249151e-15 3.69725182221479e-15 2.02114207545759e-15 -4.50870833392161e-16 9.62950493696677e-17;1.00935904205024e-17 6.54751873609395e-17 -1.09138810997186e-16 -8.62396750098759e-17 -3.82788257844306e-17;0 0 0 0 0;4.21958510903678e-07 -8.30678271007705e-08 -3.47006439555247e-07 -3.36442823712421e-08 9.90739768222027e-08;2.64389033612742e-08 2.65825090066479e-09 -1.28895513428522e-08 -7.07182694980098e-10 7.10907165301180e-09;6.31203524153492e-09 -1.67038260990134e-09 1.33104703539822e-09 8.34376495185149e-10 -2.52478613522612e-10;1.18414896299279e-10 -2.57745052288455e-11 2.88295935685818e-11 -3.27782977418354e-11 -1.05705000036156e-11;-4.20826459055091e-12 -6.97430607432268e-12 -3.90660545970607e-12 -3.90449239948755e-13 -4.60384797517466e-13;-9.47668356558200e-13 6.53305025354881e-13 2.63240185434960e-13 1.40129115015734e-13 3.85788887132074e-14;2.23947810407291e-13 7.35262771548253e-15 -3.83348211931292e-14 4.20376514344176e-14 4.26445836468461e-14;-3.88008154470596e-16 2.28561424667750e-15 -8.73599966653373e-16 2.14321147947665e-15 6.38631825071920e-16;-8.62165565535721e-15 1.79742912149810e-15 1.01541125038661e-15 -7.91027655831866e-17 -4.06505132825230e-16;-2.35355054392189e-16 -6.13997759731013e-17 -2.73490528665965e-17 2.63895177155121e-17 -4.47531057245187e-18;6.01909706823530e-17 5.35520010856833e-18 -2.15530106132531e-18 -2.46778496746231e-18 -7.09947296442799e-19;0 0 0 0 0;-3.75005956318736e-07 -5.39872297906819e-07 -1.19929654883034e-07 4.52771083775007e-08 1.82790552943564e-07;7.82606642505646e-09 -1.68890832383153e-08 -8.45995188378997e-09 1.42958730598502e-09 3.21075754133531e-09;4.28818421913782e-09 -1.07501469928219e-09 8.84086350297418e-10 9.74171228764155e-10 8.59877149602304e-12;1.28983712172521e-10 -6.96375160373676e-11 -2.13481436408896e-11 1.33516375568179e-11 -1.65864626508258e-11;-4.48914384622368e-12 9.68953616831263e-13 -1.61372463422897e-12 -2.09683563440448e-12 -1.90096826314068e-12;-1.12626619779175e-13 3.34903159106509e-14 -1.21721528343657e-13 7.46246339290354e-14 3.68424909859186e-13;5.08294274367790e-14 2.83036159977090e-14 1.48074873486387e-14 -9.59633528834945e-15 -1.26231060951100e-14;-4.01464098583541e-16 1.97047929526674e-15 -5.29967950447497e-16 -3.59120406619931e-16 1.69690933982683e-16;-1.73919209873841e-15 7.52792462841274e-16 3.65589287101147e-16 -7.79247612043812e-17 -8.24599670368999e-17;-4.61555616150128e-17 4.94529746019753e-19 -1.09858157212270e-17 3.95550811124928e-18 3.23972399884100e-18;-2.27040686655766e-17 -3.27855689001215e-18 -3.30649011116861e-19 9.08748546536849e-19 8.92197599890994e-19;5.67241944733762e-18 3.84449400209976e-19 1.77668058015537e-19 2.00432838283455e-20 -2.00801461564767e-19];
            bnm_bw = [0 0 0 0 0;0 0 0 0 0;-9.56715196386889e-06 -3.68040633020420e-08 1.27846786489883e-07 1.32525487755973e-06 1.53075361125066e-06;0 0 0 0 0;-7.17682617983607e-06 2.89994188119445e-06 -2.97763578173405e-07 8.95742089134942e-07 3.44416325304006e-07;-8.02661132285210e-07 3.66738692077244e-07 -3.02880965723280e-07 3.54144282036103e-07 -1.68873066391463e-07;0 0 0 0 0;-2.89640569283461e-06 -7.83566373343614e-07 -8.36667214682577e-07 -7.41891843549121e-07 -9.23922655636489e-08;-1.06144662284862e-06 1.57709930505924e-07 1.04203025714319e-07 1.20783300488461e-07 -1.38726055821134e-07;-4.16549018672265e-07 -1.35220897698872e-07 -6.40269964829901e-08 1.63258283210837e-08 -2.57958025095959e-08;0 0 0 0 0;3.52324885892419e-06 -2.26705543513814e-07 1.53835589488292e-06 -3.75263061267433e-07 3.69384057396017e-07;-2.06569149157664e-07 -9.36260183227175e-08 -3.55985284353048e-08 -9.13671163891094e-08 6.93156256562600e-09;1.32437594740782e-07 4.44349887272663e-08 -3.38192451721674e-08 -3.97263855781102e-08 -1.93087822995800e-09;-1.29595244818942e-07 -1.40852985547683e-08 1.42587592939760e-09 7.05779876554001e-09 -1.00996269264535e-08;0 0 0 0 0;4.06960756215938e-06 -1.97898540226986e-06 7.21905857553588e-08 -1.19908881538755e-06 -5.67561861536903e-08;6.53369660286999e-08 -2.42818687866392e-07 -1.66203004559493e-08 -2.41512414151897e-08 4.45426333411018e-08;1.44650670663281e-07 8.50666367433859e-09 -4.61165612004307e-09 4.88527987491045e-09 1.06277326713172e-08;1.86770937103513e-08 -6.44197940288930e-10 -7.60456736846174e-09 -9.97186468682689e-10 8.73229752697716e-10;-1.00206566229113e-08 1.33934372663121e-09 1.41691503439220e-09 8.72352590578753e-10 -8.04561626629829e-10;0 0 0 0 0;3.07161843116618e-06 1.82962085656470e-06 1.87728623016069e-07 7.10611617623261e-07 2.26499092250481e-07;4.50766403064905e-08 -1.67752393078256e-07 2.47844723639070e-08 -3.56484348424869e-09 -1.56634836636584e-08;3.77011651881090e-08 -7.23045828480496e-09 5.22995988863761e-09 -1.03740320341306e-09 4.57839777217789e-09;8.09495635883121e-09 -3.01977244420529e-10 -2.30104544933093e-09 3.63658580939428e-10 4.39320811714867e-10;9.37087629961269e-11 1.00780920426635e-09 1.28140539913350e-10 -6.65795285522138e-12 4.71732796198631e-11;-8.88504487069155e-11 -1.63253810435461e-10 7.22669710644299e-11 5.64715132584527e-11 -1.08949308197617e-12;0 0 0 0 0;-2.64054293284174e-07 -2.37611606117256e-06 -1.83671059706264e-06 -3.12199354841993e-07 -1.05598289276114e-07;7.41706968747147e-08 -1.64359098062646e-08 -3.09750224040234e-08 -9.68640079410317e-09 -7.90399057863403e-08;-1.00254376564271e-08 1.12528248631191e-08 -2.67841549174100e-09 -2.69481819323647e-09 1.56550607475331e-09;-2.18568129350729e-09 6.26422056977450e-10 1.95007291427316e-09 3.14226463591125e-10 -3.62000388344482e-10;-9.30451291747549e-10 5.62175549482704e-11 1.01022849902012e-10 5.18675856498499e-11 5.37561696283235e-11;5.33151334468794e-11 1.07571307336725e-10 -1.31714567944652e-11 -4.17524405900018e-11 -2.16737797893502e-12;4.69916869001309e-11 -4.34516364859583e-12 -6.61054225868897e-12 -5.75845818545368e-12 -2.32180293529175e-12;0 0 0 0 0;-3.50305843086926e-06 1.76085131953403e-06 8.16661224478572e-07 4.09111042640801e-07 -9.85414469804995e-08;1.44670876127274e-07 -1.41331228923029e-08 -3.06530152369269e-08 -1.46732098927996e-08 -2.30660839364244e-08;-2.00043052422933e-08 1.72145861031776e-09 2.13714615094209e-09 1.02982676689194e-09 -1.64945224692217e-10;1.23552540016991e-09 1.42028470911613e-09 8.79622616627508e-10 -7.44465600265154e-10 -7.17124672589442e-11;-6.67749524914644e-10 -5.77722874934050e-11 3.40077806879472e-11 4.26176076541840e-11 8.23189659748212e-11;-4.62771648935992e-11 -7.24005305716782e-13 1.18233730497485e-12 5.18156973532267e-12 -1.53329687155297e-12;4.75581699468619e-12 -3.79782291469732e-12 1.33077109836853e-12 -1.02426020107120e-12 3.10385019249130e-13;1.66486090578792e-12 1.08573672403649e-12 1.26268044166279e-13 -1.23509297742757e-13 -1.81842007284038e-13;0 0 0 0 0;9.93870680202303e-08 -1.85264736035628e-06 -5.58942734710854e-07 -5.54183448316270e-07 -3.95581289689398e-08;7.88329069002365e-08 2.04810091451078e-08 3.74588851000076e-09 3.42429296613803e-08 -2.00840228416712e-08;-5.93700447329696e-10 -6.57499436973459e-10 -6.90560448220751e-09 3.56586371051089e-09 7.33310245621566e-11;-6.38101662363634e-11 4.23668020216529e-10 -2.43764895979202e-10 -9.31466610703172e-11 -3.17491457845975e-10;1.50943725382470e-11 -6.11641188685078e-11 -4.37018785685645e-11 -2.32871158949602e-11 4.19757251950526e-11;-1.18165328825853e-11 -9.91299557532438e-13 6.40908678055865e-14 2.41049422936434e-12 -8.20746054454953e-14;6.01892101914838e-12 -8.78487122873450e-13 -1.58887481332294e-12 -3.13556902469604e-13 5.14523727801645e-14;-1.50791729401891e-13 -1.45234807159695e-13 1.65302377570887e-13 -5.77094211651483e-15 9.22218953528393e-14;-1.85618902787381e-14 5.64333811864051e-14 -9.94311377945570e-15 -2.40992156199999e-15 -2.19196760659665e-14;0 0 0 0 0;-8.16252352075899e-08 1.61725487723444e-06 9.55522506715921e-07 4.02436267433511e-07 -2.80682052597712e-07;7.68684790328630e-09 -5.00940723761353e-09 -2.43640127974386e-08 -2.59119930503129e-08 3.35015169182094e-08;7.97903115186673e-09 3.73803883416618e-09 3.27888334636662e-09 1.37481300578804e-09 -1.10677168734482e-10;-1.67853012769912e-09 -1.61405252173139e-10 -1.98841576520056e-10 -1.46591506832192e-11 9.35710487804660e-11;4.08807084343221e-11 -3.74514169689568e-11 -3.03638493323910e-11 -5.02332555734577e-12 -8.03417498408344e-12;6.48922619024579e-12 1.96166891023817e-12 -1.96968755122868e-12 -5.20970156382361e-12 -1.62656885103402e-12;1.28603518902875e-12 -4.88146958435109e-13 -3.37034886991840e-13 1.37393696103000e-14 4.41398325716943e-14;1.48670014793021e-13 4.41636026364555e-14 2.06210477976005e-14 -3.43717583585390e-14 -1.21693704024213e-14;-1.67624180330244e-14 6.59317111144238e-15 2.57238525440646e-15 -3.21568425020512e-17 5.29659568026553e-15;7.85453466393227e-16 6.91252183915939e-16 -1.20540764178454e-15 -3.85803892583301e-16 3.46606994632006e-16;0 0 0 0 0;2.86710087625579e-06 -1.68179842305865e-06 -8.48306772016870e-07 -7.08798062479598e-07 -1.27469453733635e-07;2.11824305734993e-09 2.02274279084379e-08 1.61862253091554e-08 3.25597167111807e-08 3.40868964045822e-09;1.21757111431438e-08 1.68405530472906e-09 1.55379338018638e-09 -3.81467795805531e-10 2.53316405545058e-09;-9.98413758659768e-11 5.38382145421318e-10 3.92629628330704e-10 -1.43067134097778e-10 3.74959329667113e-12;-1.57270407028909e-11 -9.02797202317592e-12 8.45997059887690e-12 4.71474382524218e-12 5.41880986596427e-12;-1.20658618702054e-12 7.12940685593433e-13 1.02148613026937e-12 1.63063852348169e-13 1.74048793197708e-13;3.80559390991789e-13 1.19678271353485e-13 9.72859455604188e-14 5.42642400031729e-14 8.18796710714586e-14;-4.69629218656902e-14 5.59889038686206e-15 2.05363292795059e-15 5.38599403288686e-15 -2.68929559474202e-15;-1.88759348081742e-14 5.20975954705924e-15 -4.43585653096395e-16 5.57436617793556e-16 -3.95922805817677e-16;-9.80871456373282e-16 2.50857658461759e-17 -1.24253000050963e-16 6.00857065211394e-17 3.53799635311500e-18;2.49370713054872e-16 -1.49119714269816e-17 -3.12276052640583e-17 -2.42001662334001e-17 -1.69766504318143e-17;0 0 0 0 0;-1.69222102455713e-06 1.64277906173064e-06 5.28855114364096e-07 4.28159853268650e-07 -1.57362445882665e-07;1.67656782413678e-08 -3.77746114074055e-08 -2.21564555842165e-08 -3.37071806992217e-08 1.47454008739800e-08;1.06080499491408e-08 3.21990403709678e-09 3.87301757435359e-09 2.92241827834347e-10 -1.86619473655742e-11;1.62399669665839e-10 3.51322865845172e-10 2.67086377702958e-11 -1.31596563625491e-10 3.14164569507034e-11;-2.02180016657259e-11 2.03305178342732e-11 6.34969032565839e-12 5.99522296668787e-12 -4.46275273451008e-12;-9.88409290158885e-13 -1.47692750858224e-13 3.14655550730530e-13 -2.41857189187879e-13 4.47727504501486e-13;1.71430777754854e-13 1.73950835042486e-13 5.92323956541558e-14 8.06625710171825e-15 2.33252485755634e-14;-1.74184545690134e-15 -8.18003353124179e-16 -6.62369006497819e-16 4.16303374396147e-15 7.06513748014024e-15;-6.02936238677014e-15 1.89241084885229e-15 1.99097881944270e-17 -6.99974290696640e-16 -2.69504942597709e-17;-4.65632962602379e-16 3.70281995445114e-18 -9.04232973763345e-17 2.20847370761932e-17 7.62909453726566e-17;-6.25921477907943e-17 -2.10532795609842e-17 -1.03808073867183e-17 1.15091380049019e-18 4.66794445408388e-19;9.39427013576903e-18 9.17044662931859e-19 2.04132745117549e-18 -1.72364063154625e-19 -1.18098896532163e-18];
            bnm_ch = [0 0 0 0 0;0 0 0 0 0;3.44092035729033e-05 -1.21876825440561e-05 -1.87490665238967e-05 -2.60980336247863e-05 4.31639313264615e-06;0 0 0 0 0;-2.60125613000133e-05 1.70570295762269e-05 3.08331896996832e-05 1.66256596588688e-05 -1.07841055501996e-05;8.74011641844073e-06 -2.25874169896607e-06 6.50985196673747e-07 1.30424765493752e-06 -1.85081244549542e-07;0 0 0 0 0;3.77496505484964e-05 -1.08198973553337e-05 -1.67717574544937e-05 -3.22476096673598e-05 1.12281888201134e-05;-7.68623378647958e-07 -4.01400837153063e-06 -2.16390246700835e-06 -1.76912959937924e-06 -1.12740084951955e-06;-2.37092815818895e-06 -9.52317223759653e-07 -2.22722065579131e-07 -6.25157619772530e-08 1.86582003894639e-08;0 0 0 0 0;-6.10254317785872e-05 -2.51815503068494e-05 2.01046207874667e-05 7.21107723367308e-06 -1.30692058660457e-05;-9.60655417241537e-06 -7.31381721742373e-06 -2.52767927589636e-06 9.09039973214621e-07 -6.76454911344246e-07;-2.25743206384908e-08 2.33058746737575e-07 2.24746779293445e-07 6.78551351968876e-08 1.25076011387284e-07;-2.25744112770133e-07 -1.44429560891636e-07 -2.96810417448652e-08 -5.93858519742856e-08 -2.43210229455420e-08;0 0 0 0 0;7.45721015256308e-06 -3.81396821676410e-05 -1.41086198468687e-05 -2.28514517574713e-05 7.28638705683277e-06;-5.77517778169692e-06 -3.93061211403839e-06 -2.17369763310752e-06 -1.48060935583664e-07 -2.74200485662814e-07;4.52962035878238e-07 9.80990375495214e-07 4.67492045269286e-07 -8.31032252212116e-09 1.69426023427740e-07;7.20536791795515e-10 2.75612253452141e-09 2.47772119382536e-09 4.30621825021233e-09 -2.86498479499428e-08;-2.46253956492716e-08 -3.10300833499669e-09 8.06559148724445e-09 2.98197408430123e-10 6.32503656532846e-09;0 0 0 0 0;-6.01147094179306e-05 -3.16631758509869e-05 4.10038115100010e-06 3.55215057231403e-07 -2.23606515237408e-06;-2.85937516921923e-06 -3.67775706610630e-06 -5.06445540401637e-07 8.21776759711184e-07 -5.98690271725558e-07;7.77122595418965e-07 3.60896376754085e-07 3.88610487893381e-07 -4.39533892679537e-08 -6.26882227849174e-08;1.05759993661891e-07 2.58009912408833e-08 -1.51356049060972e-08 -1.13335813107412e-09 5.37470857850370e-10;7.99831506181984e-09 1.67423735327465e-09 2.94736760548677e-09 -1.56727133704788e-09 8.46186800849124e-10;3.07727104043851e-09 3.93584215798484e-10 3.86721562770643e-11 1.72181091277391e-10 -2.16915737920145e-10;0 0 0 0 0;-1.16335389078126e-05 -1.39864676661484e-05 2.52546278407717e-06 -8.79152625440188e-06 -8.97665132187974e-06;-3.95874550504316e-06 -1.17976262528730e-07 7.03189926369300e-07 3.38907065351535e-07 -3.67714052493558e-07;2.29082449370440e-07 5.72961531093329e-07 4.21969662578894e-08 1.24112958141431e-08 9.56404486571888e-08;1.44631865298671e-09 6.19368473895584e-09 1.67110424041236e-09 2.57979463602951e-09 -6.90806907510366e-09;1.77235802019153e-09 -8.14388846228970e-10 4.50421956523579e-09 5.67452314909707e-10 2.47610443675560e-09;4.85932343880617e-10 2.24864117422804e-10 -2.22534534468511e-10 -7.96395824973477e-11 3.12587399902493e-12;-3.20173937255409e-11 -1.29872402028088e-11 -4.24092901203818e-11 2.66570185704416e-11 -5.25164954403909e-12;0 0 0 0 0;-1.36010179191872e-05 1.77873053642413e-05 4.80988546657119e-06 3.46859608161212e-06 -1.73247520896541e-06;2.00020483116258e-06 2.43393064079673e-06 1.21478843695862e-06 1.95582820041644e-07 -3.11847995109088e-07;-8.13287218979310e-09 1.05206830238665e-08 6.54040136224164e-09 -1.96402660575990e-08 -1.40379796070732e-08;4.01291020310740e-08 2.92634301047947e-08 6.04179709273169e-09 8.61849065020545e-10 5.98065429697245e-09;-1.06149335032911e-09 -4.39748495862323e-10 8.83040310269353e-10 3.49392227277679e-10 8.57722299002622e-10;-1.25049888909390e-11 2.05203288281631e-10 1.37817670505319e-11 6.82057794430145e-11 -9.41515631694254e-11;7.47196022644130e-12 -2.51369898528782e-11 -2.12196687809200e-11 1.55282119505201e-11 9.99224438231805e-12;-7.90534019004874e-13 3.55824506982589e-12 8.00835777767281e-13 8.73460019069655e-13 1.34176126600106e-12;0 0 0 0 0;3.12855262465316e-05 1.31629386003608e-05 2.65598119437581e-06 8.68923340949135e-06 -7.51164082949678e-06;1.56870792650533e-06 1.89227301685370e-06 4.15620385341985e-07 -2.74253787880603e-07 -4.28826210119200e-07;-9.99176994565587e-08 -1.10785129426286e-07 -1.10318125091182e-07 6.22726507350764e-09 -3.39214566386250e-08;1.24872975018433e-08 1.10663206077249e-08 5.40658975901469e-09 -2.79119137105115e-09 -2.47500096192502e-09;1.11518917154060e-10 -4.21965763244849e-10 3.26786005211229e-10 1.93488254914545e-10 7.00774679999972e-10;1.50889220040757e-10 1.03130002661366e-10 -3.09481760816903e-11 -4.47656630703759e-11 -7.36245021803800e-12;-1.91144562110285e-12 -1.11355583995978e-11 -1.76207323352556e-11 8.15289793192265e-12 3.45078925412654e-12;-2.73248710476019e-12 -1.65089342283056e-13 -2.20125355220819e-13 5.32589191504356e-13 5.70008982140874e-13;8.06636928368811e-13 1.30893069976672e-13 9.72079137767479e-14 3.87410156264322e-14 -5.56410013263563e-14;0 0 0 0 0;2.02454485403216e-05 -9.77720471118669e-06 -4.35467548126223e-06 2.19599868869063e-06 -3.26670819043690e-06;-3.21839256310540e-08 8.38760368015005e-07 -5.08058835724060e-07 4.16177282491396e-08 1.53842592762120e-07;-1.57377633165313e-07 -7.86803586842404e-08 -7.40444711426898e-08 3.15259864117954e-08 5.60536231567172e-09;-3.26080428920229e-10 -3.14576780695439e-09 8.46796096612981e-10 -2.59329379174262e-09 -8.01054756588382e-10;-4.58725236153576e-11 -6.87847958546571e-11 8.18226480126754e-12 1.81082075625897e-10 1.74510532938256e-10;7.60233505328792e-11 4.76463939581321e-11 -2.47198455442033e-11 -8.83439688929965e-12 5.93967446277316e-13;-8.92919292558887e-12 -4.38524572312029e-12 -4.02709146060896e-12 4.84344426425295e-12 5.12869042781520e-12;1.91518361809952e-12 3.06846255371817e-13 -2.44830265306345e-13 7.86297493099244e-14 2.72347805801980e-13;9.09936624159538e-14 7.20650818861447e-15 2.45383991578283e-14 -4.79580974186462e-15 3.64604724046944e-14;-4.63611142770709e-14 1.73908246420636e-15 -4.41651410674801e-15 -6.61409045306922e-16 -1.60016049099639e-15;0 0 0 0 0;6.17105245892845e-06 -1.04342983738457e-05 -1.72711741097994e-05 -8.16815967888426e-07 3.42789959967593e-06;-2.44014060833825e-07 2.06991837444652e-07 -3.85805819475679e-07 1.67162359832166e-08 4.15139610402483e-07;8.18199006804020e-08 -3.20013409049159e-08 5.94000906771151e-08 2.24122167188946e-08 -1.33796186160409e-08;7.66269294674338e-11 -6.07862178874828e-10 4.95795757186248e-10 -3.07589245481422e-10 3.44456287710689e-10;-1.84076250254929e-10 -1.30985312312781e-10 -1.52547325533276e-10 -2.51000125929512e-11 -1.93924012590455e-11;-2.93307452197665e-11 2.88627386757582e-11 5.58812021182217e-12 -1.68692874069187e-13 1.80464313900575e-12;-9.59053874473003e-13 6.04803122874761e-13 -9.80015608958536e-13 1.70530372034214e-12 1.70458664160775e-12;2.80169588226043e-13 9.09573148053551e-14 2.16449186617004e-14 1.15550091496353e-13 4.97772796761321e-14;-3.04524400761371e-14 3.42845631349694e-14 2.44230630602064e-14 5.76017546103056e-16 -9.74409465961093e-15;5.98765340844291e-15 -2.63942474859535e-15 -1.80204805804437e-15 -1.84981819321183e-16 -5.85073392163660e-16;-2.37069441910133e-15 2.87429226086856e-16 -1.67055963193389e-16 2.72110684914090e-18 8.46646962667892e-17;0 0 0 0 0;-2.71386164105722e-05 -1.41834938338454e-05 -2.00777928859929e-07 5.94329804681196e-07 8.61856994375586e-06;-3.93656495458664e-08 -6.36432821807576e-07 -2.47887475106438e-07 -2.64906446204966e-08 1.10689794197004e-07;5.25319489188562e-08 9.00866357158695e-09 5.00693379572512e-08 2.47269011056404e-08 -7.27648556194598e-09;1.87207107149043e-09 -1.46428282396138e-09 -2.71812237167257e-10 8.44902265891466e-10 -5.62683870906027e-10;-1.08295119666184e-10 4.75553388543793e-11 -5.49429386495686e-11 -6.60907871731611e-11 -5.97347322824822e-11;-4.95118306815571e-12 5.31083735234970e-13 -1.93679746327378e-12 -1.61770521840510e-12 1.23276727202510e-11;6.68582682909900e-13 7.38288575160449e-13 5.47630483499201e-13 -1.00770258118914e-13 -1.65564928475981e-13;5.80963409268471e-14 6.93474288078737e-14 6.60728092794315e-15 -5.21029056725202e-15 -1.11283532854883e-16;-4.10567742688903e-15 1.62252646805882e-14 1.00774699865989e-14 -2.44793214897877e-16 -1.59283906414563e-15;1.84669506619904e-17 8.28473337813919e-17 -1.53400662078899e-16 -5.01060672199689e-17 -2.20727935766132e-16;2.65355116203636e-16 -3.70233146147684e-17 3.52689394451586e-18 -8.62215942516328e-18 9.26909361974526e-18;9.94266950643135e-17 4.17028699663441e-18 -7.65153491125819e-21 -5.62131270981041e-18 -3.03732817297438e-18];
            bnm_cw = [0 0 0 0 0;0 0 0 0 0;-0.000209104872912563 -1.41530274973540e-05 3.00318745764815e-05 -1.82864291318284e-05 -7.62965409959238e-06;0 0 0 0 0;-0.000186336519900275 0.000191256553935638 7.28356195304996e-05 3.59637869639906e-05 -2.53927226167388e-05;0.000108195343799485 -6.97050977217619e-05 -6.68037133871099e-05 2.30387653190503e-05 -1.22735483925784e-05;0 0 0 0 0;0.000119941091277039 -7.70547844186875e-05 -8.15376297964528e-05 1.06005789545203e-05 2.31177232268720e-05;-1.77494760217164e-05 -1.37061385686605e-05 -1.74805936475816e-05 -6.91745900867532e-07 -7.10231790947787e-06;-1.47564103733219e-05 2.08890785485260e-06 3.19876879447867e-06 9.43984664503715e-07 -4.90480527577521e-06;0 0 0 0 0;4.93300138389457e-05 -6.77641298460617e-05 -3.25043347246397e-05 8.33226714911921e-06 8.11499972792905e-06;-2.80449863471272e-05 -1.04367606414606e-05 1.64473584641163e-07 -3.57420965807816e-06 2.95887156564038e-06;1.88835280111533e-06 5.69125761193702e-07 -2.22757382799409e-06 -1.96699131032252e-07 -2.91861219283659e-07;-4.69918971436680e-06 -7.00778948636735e-07 2.97544157334673e-09 3.86100512544410e-07 2.30939653701027e-07;0 0 0 0 0;1.77050610394149e-05 -3.18353071311574e-05 3.04232260950316e-05 -6.26821316488169e-05 -1.75094810002378e-06;9.25605901565775e-06 -8.25179123302247e-06 6.74032752408358e-06 3.22192289084524e-06 6.09414500075259e-06;4.28233825242200e-06 2.10470570087927e-07 -4.75050074985668e-07 -4.89382663470592e-07 8.75232347469207e-07;8.50393520366934e-07 1.58764911467186e-07 -2.16267638321210e-07 -7.43341300487416e-10 1.75131729813230e-07;-2.87064111623119e-07 4.50393893102830e-08 6.63315044416690e-08 7.61199387418853e-08 -6.05694385243652e-09;0 0 0 0 0;-1.95692079507947e-05 5.15486098887851e-05 3.00852761598173e-05 1.21485028343416e-05 -6.72450521493428e-06;5.34496867088158e-06 3.90973451680699e-06 3.70148924718425e-06 5.73731499938212e-08 5.52258220288780e-07;3.39950838185315e-07 -5.63443976772634e-07 4.52082211980595e-07 -2.57094645806243e-07 -6.84885762924729e-08;2.15793276880684e-07 2.05911354090873e-07 1.33747872341142e-08 -2.07997626478952e-08 -3.69812938736019e-08;2.11952749403224e-09 4.04317822544732e-08 2.40972024883650e-09 8.56289126938059e-09 2.31035283490200e-08;-2.08402298813248e-09 -8.50243600879112e-09 2.60895410117768e-09 -6.69156841738591e-10 -5.16280278087006e-09;0 0 0 0 0;0.000124901291436683 -5.70770326719086e-05 -8.44887248105015e-05 -3.11442665354698e-05 -1.12982893252046e-05;-8.38934444233944e-06 1.56860091415414e-06 -1.77704563531825e-06 -5.70219068898717e-08 -4.30377735031244e-06;3.72965318017681e-07 6.98175439446187e-07 1.75760544807919e-08 1.59731284857151e-07 3.62363848767891e-07;-2.32148850787091e-07 -4.21888751852973e-08 8.35926113952108e-08 -2.24572480575674e-08 -6.92114100904503e-08;-2.92635642210745e-09 3.38086229163415e-09 4.72186694662901e-09 -8.32354437305758e-11 4.19673890995627e-09;-1.26452887692900e-09 1.91309690886864e-09 1.54755631983655e-09 -1.09865169400249e-09 1.83645326319994e-10;9.92539437011905e-10 -2.96318203488300e-10 1.17466020823486e-10 -5.00185957995526e-10 -8.54777591408537e-11;0 0 0 0 0;-0.000182885335404854 7.27424724520089e-05 3.05286278023427e-05 2.55324463432562e-05 -6.39859510763234e-06;-5.21449265232557e-06 -6.70572386081398e-06 -3.95473351292738e-06 -6.41023334372861e-07 -3.11616331059009e-06;2.37090789071727e-07 3.58427517014705e-07 2.55709192777007e-07 8.44593804408541e-08 9.27243162355359e-09;7.24370898432057e-08 -7.43945120337710e-09 8.61751911975683e-10 -2.34651212610623e-08 2.94052921681456e-09;-1.22127317934425e-08 -3.89758984276768e-09 4.12890383904924e-11 2.06528068002723e-09 1.73488696972270e-09;-5.44137406907620e-10 -4.81034553189921e-10 -2.56101759039694e-11 3.21880564410154e-10 -2.70195343165250e-11;1.08394225300546e-10 -7.99525492688661e-11 1.73850287030654e-10 -8.06390014426271e-11 -7.63143364291160e-13;-3.41446959267441e-11 2.72675729042792e-11 5.69674704865345e-12 -3.38402998344892e-12 -2.96732381931007e-12;0 0 0 0 0;2.91161315987250e-05 -7.24641166590735e-05 -8.58323519857884e-06 -1.14037444255820e-05 1.32244819451517e-05;1.24266748259826e-06 -4.13127038469802e-06 -8.47496394492885e-07 5.48722958754267e-07 -1.98288551821205e-06;-1.70671245196917e-08 1.36891127083540e-08 -2.80901972249870e-07 -5.45369793946222e-09 -9.58796303763498e-08;1.14115335901746e-08 2.79308166429178e-08 -1.71144803132413e-08 4.86116243565380e-09 -8.13061459952280e-09;-1.19144311035824e-09 -1.28197815211763e-09 -1.22313592972373e-09 6.23116336753674e-10 2.11527825898689e-09;4.94618645030426e-10 -1.01554483531252e-10 -3.58808808952276e-10 1.23499783028794e-10 -1.21017599361833e-10;1.33959569836451e-10 -1.87140898812283e-11 -3.04265350158941e-11 -1.42907553051431e-11 -1.09873858099638e-11;1.30277419203512e-11 -4.95312627777245e-12 2.23070215544358e-12 1.66450226016423e-12 6.26222944728474e-12;-4.40721204874728e-12 2.99575133064885e-12 -1.54917262009097e-12 8.90015664527060e-14 -1.59135267012937e-12;0 0 0 0 0;-4.17667211323160e-05 1.39005215116294e-05 1.46521361817829e-05 3.23485458024416e-05 -8.57936261085263e-06;9.48491026524450e-07 1.67749735481991e-06 6.80159475477603e-07 -1.34558044496631e-06 1.62108231492249e-06;-2.67545753355631e-07 -3.31848493018159e-08 1.05837219557465e-07 1.55587655479400e-07 -2.84996014386667e-08;-5.15113778734878e-08 8.83630725241303e-09 3.36579455982772e-09 -6.22350102096402e-09 5.03959133095369e-09;2.04635880823035e-11 -1.07923589059151e-09 -6.96482137669712e-10 -4.70238500452793e-10 -6.60277903598297e-10;-2.41897168749189e-11 1.33547763615216e-10 -5.13534673658908e-11 -8.32767177662817e-11 5.72614717082428e-11;7.55170562359940e-12 -1.57123461699055e-11 -1.48874069619124e-11 -7.10529462981252e-13 -7.99006335025107e-12;2.41883156738960e-12 2.97346980183361e-12 1.28719977731450e-12 -2.49240876894143e-12 6.71155595793198e-13;4.16995565336914e-13 -1.71584521275288e-13 -7.23064067359978e-14 2.45405880599037e-13 4.43532934905830e-13;3.56937508828997e-14 2.43012511260300e-14 -7.96090778289326e-14 -1.59548529636358e-14 8.99103763000507e-15;0 0 0 0 0;0.000117579258399489 -4.52648448635772e-05 -2.69130037097862e-05 -3.82266335794366e-05 -4.36549257701084e-06;-1.43270371215502e-06 1.21565440183855e-06 8.53701136074284e-07 1.52709810023665e-06 1.22382663462904e-06;3.06089147519664e-07 9.79084123751975e-08 7.96524661441178e-08 4.54770947973458e-08 2.22842369458882e-07;-9.94254707745127e-09 1.43251376378012e-08 1.93911753685160e-08 -6.52214645690987e-09 -1.97114016452408e-09;-9.20751919828404e-10 -9.44312829629076e-10 7.24196738163952e-11 -6.71801072324561e-11 2.33146774065873e-10;-1.43544298956410e-11 1.78464235318769e-10 7.69950023012326e-11 -4.22390057304453e-12 3.05176324574816e-11;-7.88053753973990e-12 -3.20207793051003e-12 1.01527407317625e-12 6.02788185858449e-12 1.14919530900453e-11;-1.21558899266069e-12 5.31300597882986e-13 3.44023865079264e-13 -6.22598216726224e-14 -5.47031650765402e-14;-4.15627948750943e-13 2.77620907292721e-13 -8.99784134364011e-14 1.07254247320864e-13 6.85990080564196e-14;-3.91837863922901e-14 9.74714976816180e-15 6.79982450963903e-15 -2.41420876658572e-15 -2.20889384455344e-15;9.25912068402776e-15 -4.02621719248224e-15 -2.43952036351187e-15 -1.97006876049866e-15 1.03065621527869e-16;0 0 0 0 0;-0.000103762036940193 4.38145356960292e-05 2.43406920349913e-05 7.89103527673736e-06 -1.66841465339160e-05;-1.18428449371744e-06 -1.30188721737259e-06 -1.88013557116650e-06 -1.01342046295303e-06 9.21813037802502e-07;1.51836068712460e-07 1.11362553803933e-07 1.55375052233052e-07 1.94450910788747e-09 -1.73093755828342e-08;-3.77758211813121e-09 1.23323969583610e-08 1.72510045250302e-09 -1.88609789458597e-09 1.28937597985937e-09;-1.07947760393523e-09 5.26051570105365e-10 -3.67657536332496e-11 3.16110123523840e-10 -3.24273198242170e-10;-2.00385649209820e-12 2.54703869682390e-11 4.08563622440851e-12 -4.83350348928636e-11 -3.98153443845079e-13;2.73094467727215e-12 5.08900664114903e-12 -7.66669089075134e-13 2.50015592643012e-12 4.29763262853853e-12;6.53946487537890e-13 -2.24958413781008e-13 6.74638861781238e-15 3.28537647613903e-14 2.54199700290116e-13;-1.09122051193505e-13 8.36362392931501e-14 -3.90750153912300e-14 -5.44915910741950e-14 2.43816947219217e-14;-1.41882561550134e-14 1.00455397812713e-14 2.63347255121581e-15 1.53043256823601e-15 2.49081021428095e-15;-1.17256193152654e-15 1.05648985031971e-16 1.31778372453016e-16 1.44815198666577e-16 -3.72532768618480e-16;2.66203457773766e-16 -7.67224608659658e-17 3.51487351031864e-18 4.10287131339291e-17 -6.72171711728514e-17];
            
            
            % read the respective lines from the matrices
            anm_bh_A0 = anm_bh(:,1);   anm_bh_A1 = anm_bh(:,2);   anm_bh_B1 = anm_bh(:,3);   anm_bh_A2 = anm_bh(:,4);   anm_bh_B2 = anm_bh(:,5);
            anm_bw_A0 = anm_bw(:,1);   anm_bw_A1 = anm_bw(:,2);   anm_bw_B1 = anm_bw(:,3);   anm_bw_A2 = anm_bw(:,4);   anm_bw_B2 = anm_bw(:,5);
            anm_ch_A0 = anm_ch(:,1);   anm_ch_A1 = anm_ch(:,2);   anm_ch_B1 = anm_ch(:,3);   anm_ch_A2 = anm_ch(:,4);   anm_ch_B2 = anm_ch(:,5);
            anm_cw_A0 = anm_cw(:,1);   anm_cw_A1 = anm_cw(:,2);   anm_cw_B1 = anm_cw(:,3);   anm_cw_A2 = anm_cw(:,4);   anm_cw_B2 = anm_cw(:,5);
            bnm_bh_A0 = bnm_bh(:,1);   bnm_bh_A1 = bnm_bh(:,2);   bnm_bh_B1 = bnm_bh(:,3);   bnm_bh_A2 = bnm_bh(:,4);   bnm_bh_B2 = bnm_bh(:,5);
            bnm_bw_A0 = bnm_bw(:,1);   bnm_bw_A1 = bnm_bw(:,2);   bnm_bw_B1 = bnm_bw(:,3);   bnm_bw_A2 = bnm_bw(:,4);   bnm_bw_B2 = bnm_bw(:,5);
            bnm_ch_A0 = bnm_ch(:,1);   bnm_ch_A1 = bnm_ch(:,2);   bnm_ch_B1 = bnm_ch(:,3);   bnm_ch_A2 = bnm_ch(:,4);   bnm_ch_B2 = bnm_ch(:,5);
            bnm_cw_A0 = bnm_cw(:,1);   bnm_cw_A1 = bnm_cw(:,2);   bnm_cw_B1 = bnm_cw(:,3);   bnm_cw_A2 = bnm_cw(:,4);   bnm_cw_B2 = bnm_cw(:,5);
            
            
            % conversions
            polDist = pi/2 - lat;
            
            
            % a.) calculate Legendre polynomials
            
            % degree n and order m
            nmax = 12;
            
            % unit vector
            x = sin(polDist)*cos(lon);
            y = sin(polDist)*sin(lon);
            z = cos(polDist);
            
            % Legendre polynomials
            V(1,1) = 1;
            W(1,1) = 0;
            V(2,1) = z * V(1,1);
            W(2,1) = 0;
            
            for n = 2:nmax
                V(n+1,1) = ((2*n-1) * z * V(n,1) - (n-1) * V(n-1,1)) / n;
                W(n+1,1) = 0;
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
            
            
            
            % b.) determine the coefficients bh, bw, ch and cw
            
            % initialize
            bh_A0 = 0; bh_A1 = 0; bh_B1 = 0; bh_A2 = 0; bh_B2 = 0;
            bw_A0 = 0; bw_A1 = 0; bw_B1 = 0; bw_A2 = 0; bw_B2 = 0;
            ch_A0 = 0; ch_A1 = 0; ch_B1 = 0; ch_A2 = 0; ch_B2 = 0;
            cw_A0 = 0; cw_A1 = 0; cw_B1 = 0; cw_A2 = 0; cw_B2 = 0;
            i = 0;
            
            for n = 0:nmax
                for m = 0:n
                    
                    i = i+1;
                    
                    bh_A0 = bh_A0 + (anm_bh_A0(i)*V(n+1,m+1) + bnm_bh_A0(i)*W(n+1,m+1));
                    bh_A1 = bh_A1 + (anm_bh_A1(i)*V(n+1,m+1) + bnm_bh_A1(i)*W(n+1,m+1));
                    bh_B1 = bh_B1 + (anm_bh_B1(i)*V(n+1,m+1) + bnm_bh_B1(i)*W(n+1,m+1));
                    bh_A2 = bh_A2 + (anm_bh_A2(i)*V(n+1,m+1) + bnm_bh_A2(i)*W(n+1,m+1));
                    bh_B2 = bh_B2 + (anm_bh_B2(i)*V(n+1,m+1) + bnm_bh_B2(i)*W(n+1,m+1));
                    
                    bw_A0 = bw_A0 + (anm_bw_A0(i)*V(n+1,m+1) + bnm_bw_A0(i)*W(n+1,m+1));
                    bw_A1 = bw_A1 + (anm_bw_A1(i)*V(n+1,m+1) + bnm_bw_A1(i)*W(n+1,m+1));
                    bw_B1 = bw_B1 + (anm_bw_B1(i)*V(n+1,m+1) + bnm_bw_B1(i)*W(n+1,m+1));
                    bw_A2 = bw_A2 + (anm_bw_A2(i)*V(n+1,m+1) + bnm_bw_A2(i)*W(n+1,m+1));
                    bw_B2 = bw_B2 + (anm_bw_B2(i)*V(n+1,m+1) + bnm_bw_B2(i)*W(n+1,m+1));
                    
                    ch_A0 = ch_A0 + (anm_ch_A0(i)*V(n+1,m+1) + bnm_ch_A0(i)*W(n+1,m+1));
                    ch_A1 = ch_A1 + (anm_ch_A1(i)*V(n+1,m+1) + bnm_ch_A1(i)*W(n+1,m+1));
                    ch_B1 = ch_B1 + (anm_ch_B1(i)*V(n+1,m+1) + bnm_ch_B1(i)*W(n+1,m+1));
                    ch_A2 = ch_A2 + (anm_ch_A2(i)*V(n+1,m+1) + bnm_ch_A2(i)*W(n+1,m+1));
                    ch_B2 = ch_B2 + (anm_ch_B2(i)*V(n+1,m+1) + bnm_ch_B2(i)*W(n+1,m+1));
                    
                    cw_A0 = cw_A0 + (anm_cw_A0(i)*V(n+1,m+1) + bnm_cw_A0(i)*W(n+1,m+1));
                    cw_A1 = cw_A1 + (anm_cw_A1(i)*V(n+1,m+1) + bnm_cw_A1(i)*W(n+1,m+1));
                    cw_B1 = cw_B1 + (anm_cw_B1(i)*V(n+1,m+1) + bnm_cw_B1(i)*W(n+1,m+1));
                    cw_A2 = cw_A2 + (anm_cw_A2(i)*V(n+1,m+1) + bnm_cw_A2(i)*W(n+1,m+1));
                    cw_B2 = cw_B2 + (anm_cw_B2(i)*V(n+1,m+1) + bnm_cw_B2(i)*W(n+1,m+1));
                    
                end
            end
            
            % adding the seasonal amplitudes for the specified doy to the mean values
            bh = bh_A0 + bh_A1*cos(doy/365.25*2*pi) + bh_B1*sin(doy/365.25*2*pi) + bh_A2*cos(doy/365.25*4*pi) + bh_B2*sin(doy/365.25*4*pi);
            bw = bw_A0 + bw_A1*cos(doy/365.25*2*pi) + bw_B1*sin(doy/365.25*2*pi) + bw_A2*cos(doy/365.25*4*pi) + bw_B2*sin(doy/365.25*4*pi);
            ch = ch_A0 + ch_A1*cos(doy/365.25*2*pi) + ch_B1*sin(doy/365.25*2*pi) + ch_A2*cos(doy/365.25*4*pi) + ch_B2*sin(doy/365.25*4*pi);
            cw = cw_A0 + cw_A1*cos(doy/365.25*2*pi) + cw_B1*sin(doy/365.25*2*pi) + cw_A2*cos(doy/365.25*4*pi) + cw_B2*sin(doy/365.25*4*pi);
            
            
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
            
            % set azimuth from -180 to 180
            az_rad = mod((az_rad+pi),2*pi)-pi;
            
            % latitude of the ionosphere piercing point
            lat_pp = asin(sin(lat_rad) * cos(psi_pp) + cos(lat_rad) * sin(psi_pp) .* cos(az_rad));
            
            % longitude of the ionosphere piercing point
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
            % Formula from macmillan(1995 or 1997 Atmospheric gradients and the VLBI terrestrial and celestial reference frames)
            % To be used with Hydrostatic mapping function and temperature fixed to 15 degree 
            % from Herring 1992 Modeling atmospheric delays in the analysis of space geodetic data
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
            % forumal from chen and herring (1997 Effects of atmospheric azimuthal asymmetry on the analysis of space geodetic data)
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
