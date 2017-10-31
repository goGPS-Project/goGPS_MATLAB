%   CLASS Receiver
% =========================================================================
%
% DESCRIPTION
%   Class to store receiver data (observations, and characteristics
%
% EXAMPLE
%   trg = Receiver();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Receiver

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 3
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
classdef Receiver < handle
    
    properties (SetAccess = private, GetAccess = private)
        cc = Constellation_Collector('GRECJ'); % local cc
        w_bar                                  % handle to waitbar
        state                                  % local handle of state;
        logger                                 % handle to logger
    end
    
    properties (SetAccess = public, GetAccess = public)
        file           % file rinex object
        rin_type       % rinex version format
        
        ant            % antenna number
        ant_type       % antenna type
        ant_delta_h    % antenna height from the ground [m]
        ant_delta_en   % antenna east/north offset from the ground [m]
        
        name           % marker name
        type           % marker type
        rin_obs_code   % list of types per constellation
        ph_shift       %
        
        xyz;           % approximate position of the receiver (XYZ geocentric)
        
        n_sat = 0;     % number of satellites
        n_freq = 0;    % number of stored frequencies
        n_spe = [];    % number of observations per epoch
                
        time = [];     % internal time ref of the stored epochs
        rate;          % observations rate;
        dt = 0;        % clock offset of the receiver
        
        active_ids     % rows of active satellites
        wl             % wave-lenght of each row of row_id
        f_id           % frequency number e.g. L1 -> 1,  L2 ->2, E1 -> 1, E5b -> 3 ...
        ss_id          % satellite system number
        prn            % pseudo-range number of the satellite
        go_id          % internal id for a certain satellite
        system         % char id of the satellite system corresponding to the row_id
        
        obs_validity   % validity of the row (does it contains usable values?)
        
        obs_code       % obs code for each line of the data matrix obs
        obs            % huge obbservation matrix with all the observables for all the systems / frequencies / ecc ...
        
        clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true
        
        rec2sat = struct( ...
            'avail_index', [], ...    % boolean [n_epoch x n_sat] availability of satellites
            'err_tropo',   [], ...    % double  [n_epoch x n_sat] tropo error
            'err_iono',    [], ...    % double  [n_epoch x n_sat] iono error
            'dtS',         [], ...    % double  [n_epoch x n_sat] staellite clok error at trasmission time
            'rel_clk_corr',[], ...    % double  [n_epoch x n_sat] relativistic correction at trasmission time
            'tot',         [], ...    % double  [n_epoch x n_sat] time of travel
            'az',          [], ...    % double  [n_epoch x n_sat] azimuth
            'el',          [], ...    % double  [n_epoch x n_sat] elevation
            'cs',          [], ...    % Core_Sky
            'XS_tx',       [] ...     % compute Satellite postion a t transmission time
            )
    end
    
    % ==================================================================================================================================================
    %  SETTER
    % ==================================================================================================================================================
    
    methods
        function this = Receiver(cc)
            % SYNTAX  this = Receiver(<cc>)
            this.initObs();
            this.logger = Logger.getInstance();
            this.state = Go_State.getCurrentSettings();
            if nargin == 1
                this.cc = cc;
            else
                this.cc = this.state.cc;
            end
            this.w_bar = Go_Wait_Bar.getInstance();
        end
        
        function initObs(this)
            % initialize the receiver obj
            this.file = [];             % file rinex object
            this.rin_type = 0;          % rinex version format
            
            this.ant          = 0;       % antenna number
            this.ant_type     = '';      % antenna type
            this.ant_delta_h  = 0;       % antenna height from the ground [m]
            this.ant_delta_en = [0 0];   % antenna east/north offset from the ground [m]
            
            this.name         = 'empty';  % marker name
            this.type         = '';       % marker type
            this.rin_obs_code = '';       % list of types per constellation
            this.ph_shift     = [];
            
            this.xyz          = [0 0 0];  % approximate position of the receiver (XYZ geocentric)
            
            this.n_sat = 0;               % number of satellites
            this.n_freq = 0;              % number of stored frequencies
            n_epo = 0;               % number of epochs stored
            this.n_spe = [];              % number of sat per epoch
            
            this.dt = 0;                  % clock offset of the receiver
            
            this.time = [];               % internal time ref of the stored epochs
            this.rate = 0;                % observations rate;
            
            this.active_ids = [];         % rows of active satellites
            this.wl         = [];         % wave-lenght of each row of row_id
            this.f_id       = [];         % frequency number e.g. L1 -> 1,  L2 ->2, E1 -> 1, E5b -> 3 ...
            this.ss_id      = [];         % satellite system number
            this.prn        = [];         % pseudo-range number of the satellite
            this.go_id      = [];         % internal id for a certain satellite
            this.system     = '';         % char id of the satellite system corresponding to the row_id
            
            this.obs_validity = [];       % validity of the row (does it contains usable values?)
            
            this.obs_code   = [];         % obs code for each line of the data matrix obs
            this.obs        = [];         % huge obbservation matrix with all the observables for all the systems / frequencies / ecc ...
            
            this.clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true
            
            this.initR2S();
        end
        
        function initR2S(this)
            % initialize satellite related parameters
            % SYNTAX: this.initR2S();
            this.rec2sat = struct( ...
                'avail_index', [], ...    % boolean [n_epoch x n_sat] availability of satellites
                'err_tropo',   [], ...    % double  [n_epoch x n_sat] tropo error
                'err_iono',    [], ...    % double  [n_epoch x n_sat] iono error
                'dtS',         [], ...    % double  [n_epoch x n_sat] staellite clok error at trasmission time
                'rel_clk_corr',[], ...    % double  [n_epoch x n_sat] relativistic correction at trasmission time
                'tot',         [], ...    % double  [n_epoch x n_sat] time of travel
                'az',          [], ...    % double  [n_epoch x n_sat] azimuth
                'el',          [], ...    % double  [n_epoch x n_sat] elevation
                'cs',          [], ...    % Core_Sky
                'XS_tx',       [] ...     % compute Satellite postion a t transmission time
                );
            
            %this.rec2sat.avail_index  = false(sum(this.cc.n_sat), 1);
            %this.rec2sat.avail_index(this.go_ids) = true;
            this.rec2sat.cs                   = Core_Sky.getInstance();
            this.rec2sat.tot          = NaN(this.getNumEpochs, this.cc.getNumSat);
            %  this.rec2sat.XS_tx     = NaN(n_epoch, n_pr); % --> consider what to initialize
        end
        
        function loadRinex(this, file_name)
            % SYNTAX:
            %   this.loadRinex(file_name)
            %
            % INPUT:
            %   filename = RINEX observation file(s)
            %
            % OUTPUT:
            %   pr1 = code observation (L1 carrier)
            %   ph1 = phase observation (L1 carrier)
            %   pr2 = code observation (L2 carrier)
            %   ph2 = phase observation (L2 carrier)
            %   dop1 = Doppler observation (L1 carrier)
            %   dop2 = Doppler observation (L2 carrier)
            %   snr1 = signal-to-noise ratio (L1 carrier)
            %   snr2 = signal-to-noise ratio (L2 carrier)
            %   time = receiver seconds-of-week
            %   week = GPS week
            %   date = date (year,month,day,hour,minute,second)
            %   pos = rover approximate position
            %   interval = observation time interval [s]
            %   antoff = antenna offset [m]
            %   antmod = antenna model [string]
            %   codeC1 = boolean variable to notify if the C1 code is used instead of P1
            %   marker = marker name [string]
            %
            % DESCRIPTION:
            %   Parses RINEX observation files.
            
            t0 = tic;
            
            this.logger.addMarkedMessage('Reading observations...');
            this.logger.newLine();
                        
            this.file =  File_Rinex(file_name, 9);
            
            if this.file.isValid()
                this.logger.addMessage(sprintf('Opening file %s for reading', file_name), 100);
                % open RINEX observation file
                fid = fopen(file_name,'r');
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
                
                % importing header informations
                eoh = this.file.eoh;
                this.parseRinHead(txt, lim, eoh);
                
                if (this.rin_type < 3)
                    % considering rinex 2
                    this.parseRin2Data(txt, lim, eoh);
                else
                    % considering rinex 3
                    this.parseRin3Data(txt, lim, eoh);
                end
                
                this.logger.addMessage(sprintf('Parsing completed in %.2f seconds', toc(t0)));
                this.logger.newLine();
            end
            
            % Compute the other useful status array of the receiver object
            this.updateStatus();
            this.active_ids = true(this.getNumObservables, 1);
            
            % remove empty observables
            this.remObs(~this.obs_validity)
        end
        
        function updateStatus(this)
            % Compute the other useful status array of the receiver object
            % SYNTAX this.updateStatus();
            [~, this.ss_id] = ismember(this.system, this.cc.SYS_C);
            this.ss_id = this.ss_id';
            this.n_freq = numel(unique(this.f_id));
            
            this.go_id = this.prn + reshape(this.cc.n_sat(this.ss_id),length(this.prn),1); %%% some time second vector is a colum some time is a line reshape added to uniform
            this.n_sat = numel(unique(this.go_id));
            
            % Compute number of satellite per epoch
            
            % considerig only epoch with code on the first frequency
            code_line = this.obs_code(:,1) == 'C' & this.f_id == 1;
            this.n_spe = sum(this.obs(code_line, :) ~= 0);
            % more generic approach bbut a lot slower
            %for e = 1 : this.getNumEpochs()
            %    this.n_spe(e) = numel(unique(this.go_id(this.obs(:,e) ~= 0)));
            %end
                        
            this.obs_validity = any(this.obs, 2);
                        
            this.rec2sat.avail_index  = false(sum(this.cc.n_sat), 1);
            this.rec2sat.avail_index(this.go_id) = true;
        end
        
        function remEpoch(this, id_epo)
            % remove epochs with a certain id
            % SYNTAX:   this.remObs(id_obs)
            
            this.obs(:,epo) = [];
            this.time.delId(id_opo);
            if numel(this.dt) == this.getNumObservables()
                this.dt(id_epo) = [];
            end
            
            this.obs_validity = any(this.obs, 2);
        end
        
        function remObs(this, id_obs)
            % remove observations with a certain id
            % SYNTAX:   this.remObs(id_obs)
            
            this.obs(id_obs,:) = [];            
            
            this.active_ids(id_obs) = [];
            this.wl(id_obs) = [];
            this.f_id(id_obs) = [];
            this.ss_id(id_obs) = [];
            this.prn(id_obs) = [];
            this.go_id(id_obs) = [];
            this.system(id_obs) = [];
            this.obs_validity(id_obs) = [];
            
            this.obs_code(id_obs, :) = [];
            if length(id_obs) == length(this.wl) %case constellation collector contains lees constellation than rinex (try and error fix)
                
            end
            
            this.rec2sat.avail_index  = false(sum(this.cc.n_sat), 1);
            this.rec2sat.avail_index(this.go_id) = true;
        end      
        
        function applyDtDrift(this)
            % add dt * v_light to pseudo ranges and phases
            if ~this.clock_corrected_obs
                cpp = Core_Pre_Processing;
                d_dt = cpp.diffAndPred(this.dt);
                [d_dt] = simpleFill1D(d_dt, abs(d_dt) > 1e-4);
                dt = cumsum(d_dt);
                
                dt_corr = repmat(dt' * Go_State.V_LIGHT, size(this.ph, 1), 1);
                
                this.pr = this.pr - dt_corr;
                this.ph = this.ph - dt_corr;
                this.clock_corrected_obs = true;
            end
        end
        
        function remDtDrift(this)
            % del dt * v_light to pseudo ranges and phases
            if this.clock_corrected_obs
                cpp = Core_Pre_Processing;
                d_dt = cpp.diffAndPred(this.dt);
                [d_dt] = simpleFill1D(d_dt, abs(d_dt) > 1e-4);
                dt = cumsum(d_dt);
                
                dt_corr = repmat(dt * Go_State.V_LIGHT, 1, this.n_sat, this.n_freq);
                
                this.pr = this.pr + dt_corr;
                this.ph = this.ph + dt_corr;
                this.clock_corrected_obs = false;
            end
        end
        
        function parseRinHead(this, txt, lim, eoh)
            % Parse the header of the Observation Rinex file
            % SYNTAX:
            %    this.parseRinHead(txt, nl)
            % INPUT:
            %    txt    raw txt of the RINEX
            %    lim    indexes to determine start-stop of a line in "txt"  [n_line x 2/<3>]
            %    eoh    end of header line
            
            h_std{1} = 'RINEX VERSION / TYPE';                  %  1
            h_std{2} = 'PGM / RUN BY / DATE';                   %  2
            h_std{3} = 'MARKER NAME';                           %  3
            h_std{4} = 'OBSERVER / AGENCY';                     %  4
            h_std{5} = 'REC # / TYPE / VERS';                   %  5
            h_std{6} = 'ANT # / TYPE';                          %  6
            h_std{7} = 'APPROX POSITION XYZ';                   %  7
            h_std{8} = 'ANTENNA: DELTA H/E/N';                  %  8
            h_std{9} = 'TIME OF FIRST OBS';                     %  9
            
            h_opt{1} = 'MARKER NUMBER';                         % 10
            h_opt{2} = 'INTERVAL';                              % 11
            h_opt{3} = 'TIME OF LAST OBS';                      % 12
            h_opt{4} = 'LEAP SECONDS';                          % 13
            h_opt{5} = '# OF SATELLITES';                       % 14
            h_opt{6} = 'PRN / # OF OBS';                        % 15
            
            h_rin2_only{1} = '# / TYPES OF OBSERV';             % 16
            h_rin2_only{2} = 'WAVELENGTH FACT L1/2';            % 17
            
            h_rin3_only{1} = 'MARKER TYPE';                     % 18
            h_rin3_only{2} = 'SYS / # / OBS TYPES';             % 19
            h_rin3_only{3} = 'SYS / PHASE SHIFT';               % 20
            h_rin3_only{4} = 'GLONASS SLOT / FRQ #';            % 21
            h_rin3_only{5} = 'GLONASS COD/PHS/BIS';             % 22
            
            h_opt_rin3_only{1} = 'ANTENNA: DELTA X/Y/Z';        % 23
            h_opt_rin3_only{2} = 'ANTENNA:PHASECENTER';         % 24
            h_opt_rin3_only{3} = 'ANTENNA: B.SIGHT XYZ';        % 25
            h_opt_rin3_only{4} = 'ANTENNA: ZERODIR AZI';        % 26
            h_opt_rin3_only{5} = 'ANTENNA: ZERODIR XYZ';        % 27
            h_opt_rin3_only{6} = 'CENTER OF MASS: XYZ';         % 28
            h_opt_rin3_only{7} = 'SIGNAL STRENGTH UNIT';        % 29
            h_opt_rin3_only{8} = 'RCV CLOCK OFFS APPL';         % 30
            h_opt_rin3_only{9} = 'SYS / DCBS APPLIED';          % 31
            h_opt_rin3_only{10} = 'SYS / PCVS APPLIED';         % 32
            h_opt_rin3_only{11} = 'SYS / SCALE FACTOR';         % 33
            
            head_field = {h_std{:} h_opt{:} h_rin2_only{:} h_rin3_only{:} h_opt_rin3_only{:}}';
            
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
            if dataset{2} == 'O'
                if (this.rin_type < 3)
                    if (dataset{4} ~= 'G')
                        % GPS only RINEX2 - mixed or glonass -> actually not working
                        %throw(MException('VerifyInput:InvalidObservationFile', 'RINEX2 is supported for GPS only dataset, please use a RINEX3 file '));
                    else
                        % GPS only RINEX2 -> ok
                    end
                else
                    % RINEX 3 file -> ok
                end
            else
                throw(MException('VerifyInput:InvalidObservationFile', 'This observation RINEX does not contain observations'));
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
            % 2) 'PGM / RUN BY / DATE'
            % ignoring
            % 3) 'MARKER NAME'
            fln = find(line2head == 3, 1, 'first'); % get field line
            if isempty(fln)
                this.name = 'NO_NAME';
            else
                this.name = strtrim(txt(lim(fln, 1) + (0:59)));
            end
            % 4) 'OBSERVER / AGENCY'
            % ignoring
            % 5) 'REC # / TYPE / VERS'
            % ignoring
            % 6) 'ANT # / TYPE'
            fln = find(line2head == 6, 1, 'first'); % get field line
            if isempty(fln)
                this.ant = '';
                this.ant_type = '';
            else
                this.ant = strtrim(txt(lim(fln, 1) + (0:20)));
                this.ant_type = strtrim(txt(lim(fln, 1) + (21:40)));
            end
            % 7) 'APPROX POSITION XYZ'
            fln = find(line2head == 7, 1, 'first'); % get field line
            if isempty(fln)
                this.xyz = [0 0 0];
            else
                tmp = sscanf(txt(lim(fln, 1) + (0:41)),'%f')';                                               % read value
                this.xyz = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 3), [0 0 0], tmp);          % check value integrity
            end
            % 8) 'ANTENNA: DELTA H/E/N'
            fln = find(line2head == 8, 1, 'first'); % get field line
            if isempty(fln)
                this.ant_delta_h = 0;
                this.ant_delta_en = [0 0];
            else
                tmp = sscanf(txt(lim(fln, 1) + (0:13)),'%f')';                                                % read value
                this.ant_delta_h = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 1), 0, tmp);         % check value integrity
                tmp = sscanf(txt(lim(fln, 1) + (14:41)),'%f')';                                               % read value
                this.ant_delta_en = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 2), [0 0], tmp);    % check value integrity
            end
            % 9) 'TIME OF FIRST OBS'
            % ignoring it's already in this.file.first_epoch, but the code to read it is the following
            %fln = find(line2head == 9, 1, 'first'); % get field line
            %tmp = sscanf(txt(lim(fln, 1) + (0:42)),'%f')';
            %first_epoch = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 6), this.file.first_epoch, GPS_Time(tmp));    % check value integrity
            %first_epoch.setGPS(~strcmp(txt(lim(fln, 1) + (48:50)),'GLO'));
            % 10) 'MARKER NUMBER'
            % ignoring
            % 11) INTERVAL
            fln = find(line2head == 11, 1, 'first'); % get field line
            if isempty(fln)
                this.rate = 0; % If it's zero it'll be necessary to compute it
            else
                tmp = sscanf(txt(lim(fln, 1) + (0:9)),'%f')';                                  % read value
                this.rate = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 1), 0, tmp);  % check value integrity
            end
            % 12) TIME OF LAST OBS
            % ignoring it's already in this.file.last_epoch, but the code to read it is the following
            % fln = find(line2head == 12, 1, 'first'); % get field line
            % tmp = sscanf(txt(lim(fln, 1) + (0:42)),'%f')';
            % last_epoch = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 6), this.file.first_epoch, GPS_Time(tmp));    % check value integrity
            % last_epoch.setGPS(~strcmp(txt(lim(fln, 1) + (48:50)),'GLO'));
            % 13) LEAP SECONDS
            % ignoring
            % 14) # OF SATELLITES
            fln = find(line2head == 14, 1, 'first'); % get field line
            if isempty(fln)
                this.n_sat = this.cc.getNumSat(); % If it's zero it'll be necessary to compute it
            else
                tmp = sscanf(txt(lim(fln, 1) + (0:5)),'%f')';                                  % read value
                this.n_sat = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 1), this.cc.getNumSat(), tmp);  % check value integrity
            end
            % 15) PRN / # OF OBS            % ignoring
            % 16) # / TYPES OF OBSERV
            if this.rin_type < 3
                fln = find(line2head == 16); % get field line
                rin_obs_code = [];
                if ~isempty(fln)
                    n_obs = sscanf(txt(lim(fln(1), 1) + (3:5)),'%d');
                    l = 1;
                    while l <= numel(fln)
                        n_line = ceil(n_obs / 9);
                        l_offset = 0;
                        while l_offset < n_line
                            rin_obs_code = [rin_obs_code sscanf(txt(lim(fln(l + l_offset), 1) + (6:59)),'%s')];
                            l_offset = l_offset + 1;
                        end
                        l = l + l_offset;
                    end
                    rin_obs_code = serialize([reshape(rin_obs_code, 2, numel(rin_obs_code) / 2); ' ' * ones(1, numel(rin_obs_code) / 2)])';
                end
                this.rin_obs_code = struct('g', rin_obs_code, 'r', rin_obs_code, 'e', rin_obs_code, 'j', rin_obs_code, 'c', rin_obs_code, 'i', rin_obs_code, 's', rin_obs_code);
                
            end
            % 17) WAVELENGTH FACT L1/2
            % ignoring
            % 18) MARKER TYPE
            % Assuming non geodetic type as default
            this.type = 'NON-GEODETIC';
            fln = find(line2head == 18, 1, 'first'); % get field line
            if ~isempty(fln)
                this.type = strtrim(txt(lim(fln, 1) + (0:19)));
            end
            
            % 19) SYS / # / OBS TYPES
            if this.rin_type >= 3
                fln = find(line2head == 19); % get field lines
                this.rin_obs_code = struct('g',[],'r',[],'e',[],'j',[],'c',[],'i',[],'s',[]);
                if ~isempty(fln)
                    l = 1;
                    while l <= numel(fln)
                        sys = char(txt(lim(fln(l), 1))+32);
                        n_obs = sscanf(txt(lim(fln(l), 1) + (3:5)),'%d');
                        n_line = ceil(n_obs / 13);
                        l_offset = 0;
                        while l_offset < n_line
                            this.rin_obs_code.(sys) = [this.rin_obs_code.(sys) sscanf(txt(lim(fln(l + l_offset), 1) + (6:59)),'%s')];
                            l_offset = l_offset + 1;
                        end
                        l = l + l_offset;
                    end
                end
                if ~isempty(strfind(this.rin_obs_code.c, '1'))
                    this.rin_obs_code.c(this.rin_obs_code.c == '1') = '2';
                    this.logger.addWarning('BeiDou band 1 is now defined as 2 -> Automatically converting the observation codes of the RINEX!');
                end
            end
            % 20) SYS / PHASE SHIFT
            fln = find(line2head == 20); % get field line
            if this.rin_type < 3
                this.ph_shift = struct('g', zeros(numel(this.rin_obs_code.g) / 3, 1));
            else
                this.ph_shift = struct('g',[],'r',[],'e',[],'j',[],'c',[],'i',[],'s',[]);
                for l = 1 : numel(fln)
                    sys = char(txt(lim(fln(l), 1)) + 32);
                    rin_obs_code = txt(lim(fln(l), 1) + (2:4));
                    obs_id = (strfind(this.rin_obs_code.(sys), rin_obs_code) - 1) / 3 + 1;
                    if isempty(this.ph_shift.(sys))
                        this.ph_shift.(sys) = zeros(numel(this.rin_obs_code.(sys)) / 3, 1);
                    end
                    shift = sscanf(txt(lim(fln(l), 1) + (6:14)),'%f');
                    if ~isempty(shift)
                        this.ph_shift.(sys)(obs_id) = shift;
                    end
                end
            end
            % 21) GLONASS SLOT / FRQ #
            % ignoring
            % 22) GLONASS COD/PHS/BIS
            % ignoring
            % 23) ANTENNA: DELTA X/Y/Z
            % ignoring
            % 24) ANTENNA:PHASECENTER
            % ignoring
            % 25) ANTENNA: B.SIGHT XYZ
            % ignoring
            % 26) ANTENNA: ZERODIR AZI
            % ignoring
            % 27) ANTENNA: ZERODIR XYZ
            % ignoring
            % 28) CENTER OF MASS: XYZ
            % ignoring
            % 29) SIGNAL STRENGTH UNIT
            % ignoring
            % 30) RCV CLOCK OFFS APPL
            % ignoring
            % 31) SYS / DCBS APPLIED
            % ignoring
            % 32) SYS / PCVS APPLIED
            % ignoring
            % 33) SYS / SCALE FACTOR
            % ignoring
            
        end
        
        function chooseDataTypes(this)
            % get the right attribute column to be used for a certain type/band couple
            % LEGACY????
            t_ok = 'CLDS'; % type
            
            rin_obs_col = struct('g', zeros(4, numel(this.cc.gps.F_VEC)), ...
                'r', zeros(4, size(this.cc.glo.F_VEC,2)), ...
                'e', zeros(4, numel(this.cc.gal.F_VEC)), ...
                'j', zeros(4, numel(this.cc.qzs.F_VEC)), ...
                'c', zeros(4, numel(this.cc.bds.F_VEC)), ...
                'i', zeros(4, numel(this.cc.irn.F_VEC)), ...
                's', zeros(4, numel(this.cc.sbs.F_VEC)));
            
            if this.rin_type >= 3
                
                for c = 1 : numel(this.cc.SYS_C)
                    sys_c = char(this.cc.SYS_C(c) + 32);
                    sys = char(this.cc.SYS_NAME{c} + 32);
                    
                    if ~isempty(this.rin_obs_code.g)
                        code = reshape(this.rin_obs_code.(sys_c), 3, numel(this.rin_obs_code.(sys_c)) / 3)';
                        b_ok = this.cc.(sys).CODE_RIN3_2BAND;  % band
                        a_ok = this.cc.(sys).CODE_RIN3_ATTRIB; % attribute
                        for t = 1 : numel(t_ok)
                            for b = 1 : numel(b_ok)
                                % get the observation codes with a certain type t_ok(t) and band b_ok(b)
                                obs = (code(:,1) == t_ok(t)) & (code(:,2) == b_ok(b));
                                if any(obs)
                                    % find the preferred observation among the available ones
                                    [a, id] = intersect(code(obs, 3), a_ok{b}); a = a(id);
                                    % save the id of the column in the rin_obs_col struct matrix
                                    rin_obs_col.(sys_c)(t, b) = find(obs & code(:,3) == a(1));
                                end
                            end
                        end
                    end
                end
                
            else % rinex 2
                keyboard;
                % to be done
            end
        end
        
    end
    
    % ==================================================================================================================================================
    %  GETTER
    % ==================================================================================================================================================
    methods
        function n_obs = getNumObservables(this)
            % get the number of observables stored in the object
            % SYNTAX: n_obs = this.getNumObservables()
            n_obs = size(this.obs, 1);
        end
        
        function n_epo = getNumEpochs(this)
            % get the number of epochs stored in the object
            % SYNTAX: n_obs = this.getNumEpochs()
            n_epo = size(this.obs, 2);
        end
        
        function n_pr = getNumPseudoRanges(this)
            % get the number of epochs stored in the object
            % SYNTAX: n_pr = this.getNumPseudoRanges()
            n_pr = sum(rec.obs_code(:,1) == 'C');
        end
        
        function n_sat = getNumSat(this)
            % get the number of epochs stored in the object
            % SYNTAX: n_sat = this.getNumSat()
            n_sat = numel(unique(this.go_id));
        end
        
        function pr = pr1(this, flag_valid, sys_c)
            % get p_range 1 (Legacy)
            % SYNTAX this.pr1(<flag_valid>, <sys_c>)
            switch nargin
                case 1
                    id = (this.active_ids) & (this.f_id == 1);
                case 2
                    if flag_valid
                        id = (this.active_ids) & (this.f_id == 1) & this.pr_validity;
                    else
                        id = (this.active_ids) & (this.f_id == 1);
                    end
                case 3
                    id = (this.active_ids) & (this.f_id == 1) & this.pr_validity & this.system' == sys_c;
            end
            pr = this.pr(id,:);
        end
        
        function pr = pr2(this, flag_valid, sys_c)
            % get p_range 2 (Legacy)
            % SYNTAX this.pr1(<flag_valid>, <sys_c>)
            switch nargin
                case 1
                    id = (this.active_ids) & (this.f_id == 2);
                case 2
                    if flag_valid
                        id = (this.active_ids) & (this.f_id == 2) & this.pr_validity;
                    else
                        id = (this.active_ids) & (this.f_id == 2);
                    end
                case 3
                    id = (this.active_ids) & (this.f_id == 2) & this.pr_validity & this.system' == sys_c;
            end
            pr = this.pr(id,:);
        end
        
        function [ph, wl] = ph1(this, flag_valid, sys_c)
            % get phase 1 (Legacy)
            % SYNTAX this.ph1(<flag_valid>, <sys_c>)
            switch nargin
                case 1
                    id = (this.active_ids) & (this.f_id == 1);
                case 2
                    if flag_valid
                        id = (this.active_ids) & (this.f_id == 1) & this.pr_validity;
                    else
                        id = (this.active_ids) & (this.f_id == 1);
                    end
                case 3
                    id = (this.active_ids) & (this.f_id == 1) & this.pr_validity & this.system' == sys_c;
            end
            ph = this.ph(id,:);
            wl = this.wl(id);
        end
        
        function [ph, wl] = ph2(this, flag_valid, sys_c)
            % get phase 2 (Legacy)
            % SYNTAX this.ph1(<flag_valid>, <sys_c>)
            switch nargin
                case 1
                    id = (this.active_ids) & (this.f_id == 2);
                case 2
                    if flag_valid
                        id = (this.active_ids) & (this.f_id == 2) & this.pr_validity;
                    else
                        id = (this.active_ids) & (this.f_id == 2);
                    end
                case 3
                    id = (this.active_ids) & (this.f_id == 2) & this.pr_validity & this.system' == sys_c;
            end
            ph = this.ph(id,:);
            wl = this.wl(id);
        end
        
        function [ph, wl] = getPhGps(this, flag_valid)
            % get phase 2 (Legacy)
            % SYNTAX this.ph1(<flag_valid>, <sys_c>)
            switch nargin
                case 1
                    id = (this.active_ids) & (this.system == 'G')';
                case 2
                    if flag_valid
                        id = (this.active_ids) & (this.system == 'G')' & this.pr_validity;
                    else
                        id = (this.active_ids) & (this.system == 'G')';
                    end
            end
            ph = this.ph(id,:);
            wl = this.wl(id);
        end
        function [obs, idx] = getObs(this, flag, system)
            % get observation and index corresponfing to the flag
            % SYNTAX this.getObsIdx(flag, <system>)
            if nargin > 2
                idx = this.getObsIdx(flag, system);
            else
                idx = this.getObsIdx(flag);
            end
            obs = this.obs(idx,:);
        end
        function [obs, idx] = getPrefObs(this, flag, system)
            % get observation and index corresponfing to the flag
            % SYNTAX this.getObsIdx(flag, <system>)
            idx = this.getPrefObsIdx(flag, system);
            obs = this.obs(idx,:);
        end
        function [idx] = getObsIdx(this, flag, system)
            % get observation index corresponfing to the flag
            % SYNTAX this.getObsIdx(flag, <system>)
            idx = sum(this.obs_code(:,1:length(flag)) == repmat(flag,size(this.obs_code,1),1),2) == length(flag);
            if nargin > 2
                idx = idx & [this.system == system]';
            end
            idx = idx .* [1:length(idx)]';
            idx(idx==0)=[];
        end
        function [idx] = getPrefObsIdx(this, flag, system)
            % get observation index corresponfing to the flag using best
            % channel according to the feinition in GPS_SS, GLONASS_SS
            % SYNTAX this.getObsIdx(flag, <system>)
            if length(flag)==3
                idx = sum(this.obs_code == repmat(flag,size(this.obs_code,1),1),2) == 3;
                idx = idx & [this.system == system]';
                %this.legger.addWarning(['Unnecessary Call obs_type already determined, use getObsIdx instead'])
            elseif length(flag) ==2
                flags = zeros(size(this.obs_code,1),3);
                sys_idx = [this.system == system]';
                sys = this.cc.getSys(system);
                band = find(sys.CODE_RIN3_2BAND == flag(2));
                if isempty(band)
                    this.logger.addError('Obs not found');
                end
                preferences = sys.CODE_RIN3_ATTRIB{band};
                sys_obs_code = this.obs_code(sys_idx,:);
                sz =size(sys_obs_code,1);
                complete_flag = [];
                for i = 1:length(preferences)
                    if sum(sys_obs_code == repmat([flag preferences(i)],sz,1))>0
                        complete_flag = [flag preferences(i)];
                    end
                end
                if isempty(complete_flag)
                    this.logger.addError('Obs not found');
                end
                flags = repmat(complete_flag,size(this.obs_code,1),1);
                idx = sum(this.obs_code == flags,2) == 3;
            else
                this.legger.addError(['Invalide length of obs code(' num3str(length(flag)) 'can not determine preferred observation'])
            end
            if nargin > 2
                
            end
            idx = idx .* [1:length(idx)]';
            idx(idx==0)=[];
        end
        function [obs] = getIonoFree(this, flag1, flag2, system)
            % get Iono free combination for the two selcted measurements
            % SYNTAX [obs] = this.getIonoFree(flag1, flag2, system)
            if not(flag1(1)=='C' | flag1(1)=='L' | flag2(1)=='C' | flag2(1)=='L')
                rec.logger.addWarning('Can produce IONO free combination for the selcted observation')
                return
            end
            if flag1(1)~=flag2(1)
                rec.logger.addWarning('Incompatible observation type')
                return
            end
            [obs1, idx1] = this.getPrefObs(flag1, system);
            [obs2, idx2] = this.getPrefObs(flag2, system);
            % put zeros to NaN
            obs1(obs1 == 0) = NaN;
            obs2(obs2 == 0) = NaN;
            
            
            %gte wavelenghts
            inv_wl1 = repmat(1./this.wl(idx1),1,size(obs1,2));
            inv_wl2 = repmat(1./this.wl(idx2),1,size(obs2,2));
            obs = ((inv_wl1).^2 .* obs1 - (inv_wl2).^2 .* obs2)./ ( (inv_wl1).^2 - (inv_wl2).^2 );
            
            % ste NaN to 0
            obs(isnan(obs)) = 0;
        end
    end
    
    % ==================================================================================================================================================
    %  FUNCTIONS used as utilities
    % ==================================================================================================================================================
    methods (Access = public)
        function syncPrPh(this)
            % remove all the observations that are not present for both phase and pseudo-range
            % SYNTAX: this.syncPrPh()
            sat = ~isnan(this.pr) & ~isnan(this.ph);
            this.pr(~sat) = nan;
            this.ph(~sat) = nan;
        end
        
        function syncPhFreq(this, f_to_sync)
            % remove all the observations that are not present in all the specified frequencies
            % SYNTAX: this.syncFreq(f_to_sync)
            
            go_ids = unique(this.go_id);
            id_f = false(size(this.f_id));
            for f = 1 : numel(f_to_sync)
                id_f = id_f | this.f_id == f_to_sync(f);
            end
            for s = 1 : numel(go_ids)
                sat = (this.go_id == go_ids(s)) & id_f;
                if numel(sat) == 1
                    this.ph(sat, :) = nan;
                else
                    id_ko = sum(isnan(this.ph(sat, :))) > 0;
                    this.ph(sat, id_ko) = nan;
                end
            end
        end
    end
    
    % ==================================================================================================================================================
    %  STATIC FUNCTIONS used as utilities
    % ==================================================================================================================================================
    methods (Static, Access = public)
        
        function syncronize2receivers(rec1, rec2)
            % remove all the observations that are not present for both phase and pseudo-range between two receivers
            if (rec1.n_freq == 2) && (rec2.n_freq == 2)
                sat = ~isnan(rec1.pr(:,:,1)) & ~isnan(rec1.pr(:,:,2)) & ~isnan(rec1.ph(:,:,1)) & ~isnan(rec1.ph(:,:,2)) & ...
                    ~isnan(rec2.pr(:,:,1)) & ~isnan(rec2.pr(:,:,2)) & ~isnan(rec2.ph(:,:,1)) & ~isnan(rec2.ph(:,:,2));
            else
                sat = ~isnan(rec1.pr(:,:,1)) & ~isnan(rec1.ph(:,:,1)) & ...
                    ~isnan(rec2.pr(:,:,1)) & ~isnan(rec2.ph(:,:,1));
            end
            rec1.pr(~sat) = nan;
            rec1.ph(~sat) = nan;
            rec2.pr(~sat) = nan;
            rec2.ph(~sat) = nan;
        end
        
        function [y0, pc, wl, ref] = prepareY0(trg, mst, lambda, pivot)
            % prepare y0 and pivot_correction arrays (phase only)
            % SYNTAX: [y0, pc] = prepareY0(trg, mst, lambda, pivot)
            % WARNING: y0 contains also the pivot observations and must be reduced by the pivot corrections
            %          use composeY0 to do it
            y0 = [];
            wl = [];
            pc = [];
            i = 0;
            for t = 1 : trg.n_epo
                for f = 1 : trg.n_freq
                    sat_pr = trg.p_range(t,:,f) & mst.p_range(t,:,f);
                    sat_ph = trg.phase(t,:,f) & mst.phase(t,:,f);
                    sat = sat_pr & sat_ph;
                    pc_epo = (trg.phase(t, pivot(t), f) - mst.phase(t, pivot(t), f));
                    y0_epo = ((trg.phase(t, sat, f) - mst.phase(t, sat, f)));
                    ref = median((trg.phase(t, sat, f) - mst.phase(t, sat, f)));
                    wl_epo = lambda(sat, 1);
                    
                    idx = i + (1 : numel(y0_epo))';
                    y0(idx) = y0_epo;
                    pc(idx) = pc_epo;
                    wl(idx) = wl_epo;
                    i = idx(end);
                end
            end
        end
        
        function y0 = composeY0(y0, pc, wl)
            % SYNTAX: y0 = composeY0(y0, pc, wl)
            y0 = serialize((y0 - pc) .* wl);
            y0(y0 == 0) = []; % remove pivots
        end
    end
    
    % ==================================================================================================================================================
    %  FUNCTIONS TO GET SATELLITE RELATED PARAMETER
    % ==================================================================================================================================================
    methods
        function time_tx = getTimeTx(this,sat)
            % SYNTAX:
            %   this.getTimeTx(epoch);
            %
            % INPUT:
            % OUTPUT:
            %   time_tx = transmission time
            %   time_tx =
            %
            % DESCRIPTION:
            %   Get Transmission time
            if isempty(this.rec2sat.tot)
                this.updateTOT();
            end
            time_tx = this.time - this.tot(this.rec2sat.avail_index);
            
            
        end
        function updateTOT(this)
            % SYNTAX:
            %   this.updateTOT(time_rx, dtR);
            %
            % INPUT:
            %
            % OUTPUT:
            % DESCRIPTION:
            %   Compute the signal time of travel.
            this.tot =  (this.pr1(true) - this.getErrTropo() - this.getErrIono()) / goGNSS.V_LIGHT - dtR + this.getDtS(time_rx) + this.getRelClkCorr(time_rx) - this.getGD('L1');
            
        end
        function time_of_travel = getTOT(this)
            % SYNTAX:
            %   this.getTraveltime()
            % INPUT:
            % OUTPUT:
            %   time_of_travel   = time of travel
            % DESCRIPTION:
            %   Compute the signal transmission time.
            time_of_travel = this.tot;
        end
        function dtS = getDtS(this, time_rx)
            % SYNTAX:
            %   this.getDtS(time_rx)
            %
            % INPUT:
            %   time_rx   = reception time
            %
            % OUTPUT:
            %   dtS     = satellite clock errors
            % DESCRIPTION:
            %   Compute the satellite clock error.
            dtS = zeros(size(this.avail_index));
            for s = 1 : size(dtS)
                dtS(this.avail_index(:,s)) = this.cs.clockInterpolate(this.time_rx(this.avail_index(:,s)),s);
            end
            
        end
        function group_delay = getGD(this,obs_type)
            % SYNTAX:
            %   this.getDtS(time_rx)
            %
            % INPUT:
            %  obs_type = frequency or frequency combination,possible values : L1 L2 L3
            %
            % OUTPUT:
            %   group_delay = goup delay
            % DESCRIPTION:
            %   Compute the Grup delay from DCB
            %%% !!!! ASSUME P1/P2 RECEIVER
            %%%% QZSS IRNSS ??? -> To Be investigated
            group_delay = zeros(size(this.rec2sat.avail_index))
            switch obs_type
                case 'L1'
                    dcb_factor = cc.gps.getIonoFree.alpha2;
                case 'L2'
                    dcb_factor = cc.gps.getIonoFree.alpha1;
                case 'L3'
                    dcb_factor = 0;
                otherwise
                    Logger.getInstance().addWarning(['Unknown observable ' obs_type]);
                    dcb_factor = 0;
            end
            for s = 1:size(group_delay,2)
                group_delay(this.avail_index(:,s),:) = dcb_factor*core_sky.DCB.P1P2.values(s)
            end
        end
        
        function [XS_tx_r ,XS_tx] = getXSTxRot(this, time_rx, cc)
            % SYNTAX:
            %   [XS_tx_r ,XS_tx] = this.getXSTxRot( time_rx, cc)
            %
            % INPUT:
            % time_rx = receiver time
            % cc = Constellation Collector
            % OUTPUT:
            % XS_tx = satellite position computed at trasmission time
            % XS_tx_r = Satellite postions at transimission time rotated by earth rotation occured
            % during time of travel
            % DESCRIPTION:
            %   Compute satellite positions at transmission time and rotate them by the earth rotation
            %   occured during time of travel of the signal
            [XS_tx] = this.getXSTx();
            [XS_r] = this.earthRotationCorrection(this, XS_tx, time_rx, cc);
        end
        function [XS_tx] = getXSTx(this)
            % SYNTAX:
            %   [XS_tx_frame , XS_rx_frame] = this.getXSTx()
            %
            % INPUT:
            %
            % OUTPUT:
            % XS_tx = satellite position computed at trasmission time
            % DESCRIPTION:
            % Compute satellite positions at trasmission time
            if isempty(this.tot)
                %this.updateTimeTx();
                Logger.getInstance().addError('Trasmission time still not computed')
                return
            end
            
            XS_tx  = zeros(size(this.rec2sat.avail_index));
            for s = 1 : size(XS_tx)
                idx = this.rec2sat.avail_index(:,s);
                %%% compute staeliite position a t trasmission time
                time_tx = this.time.subset(idx);
                time_tx.time_diff = time_tx.time_diff - this.rec2sat.tot(idx,s)
                [XS_tx(idx,:,:), ~] = this.rec2sat.cs.coordInterpolate(time_tx);
            end
        end
        function [XS_r] = earthRotationCorrection(this, XS)
            % SYNTAX:
            %   [XS_r] = this.earthRotationCorrection(XS)
            %
            % INPUT:
            %   XS      = positions of satellites
            %   time_rx = receiver time
            %   cc      = Constellation Collector
            % OUTPUT:
            %   XS_r    = Satellite postions rotated by earth roattion occured
            %   during time of travel
            % DESCRIPTION:
            %   Rotate the satellites position by the earth rotation
            %   occured during time of travel of the signal
            travel_time = this.getTOT();
            %%% TBD -> consider the case XS and travel_time does not match
            XS_r = zeros(size(XS));
            for s = 1 : size(XS)
                sys = this.rec2sat.cc.system(s);
                switch char(sys)
                    case 'G'
                        omegae_dot = cc.gps.ORBITAL_P.OMEGAE_DOT;
                    case 'R'
                        omegae_dot = cc.glo.ORBITAL_P.OMEGAE_DOT;
                    case 'E'
                        omegae_dot = cc.gal.ORBITAL_P.OMEGAE_DOT;
                    case 'C'
                        omegae_dot = cc.bds.ORBITAL_P.OMEGAE_DOT;
                    case 'J'
                        omegae_dot = cc.qzs.ORBITAL_P.OMEGAE_DOT;
                    case 'I'
                        omegae_dot = cc.irn.ORBITAL_P.OMEGAE_DOT;
                    otherwise
                        Logger.getInstance().addWarning('Something went wrong in satellite_positions.m\nUnrecognized Satellite system!\n');
                        omegae_dot = cc.gps.ORBITAL_P.OMEGAE_DOT;
                end
                omega_tau = omegae_dot * travel_time(this.avail_index(s,:),s);
                R3s  = cat(3,[cos(omega_tau)    sin(omega_tau)],[-sin(omega_tau)    cos(omega_tau)]); %[n_valid_epoch x 2 x 2] matrix with all travel times rotations, Z line is omitted since roattion is along Z
                XS_r(this.avail_index(s,:),s,1) = sum(R3s(:,1,:) .* XS(this.avail_index(s,:),s,1:2),3); % X
                XS_r(this.avail_index(s,:),s,2) = sum(R3s(:,2,:) .* XS(this.avail_index(s,:),s,1:2),3); % Y
                XS_r(this.avail_index(s,:),s,2) = XS(this.avail_index(s,:),s,3); % Z
            end
            
        end
        function error_tropo = getErrTropo(this, XS)
            error_iono = zeros(size(this.rec2sat.avail_index));
            for s = 1 : size(this.available_index)
                %%% compute lat lon
                [~, lat, h, lon] = cart2geod(this.XR(:,1), this.XR(:,2), this.XR(:,3));
                %%% compute az el
                if size(this.XR,2)>1
                    [az, el] = this.getAzimuthElevation(this.XR(this.avaliable_index(:,s),:) ,XS(this.avaliable_index(:,s),s,:));
                else
                    [az, el] = this.getAzimuthElevation(this.XR, XS(this.avaliable_index(:,s),s,:));
                end
                
                idx = this.available_index(:,s);
                switch this.gs.cur_settings.tropo_model
                    case 0 %no model
                        
                    case 1 %Saastamoinen with standard atmosphere
                        error_tropo(idx,s) = Atmosphere.saastamoinen_model(lat, lon, h, el);
                        
                    case 2 %Saastamoinen with GPT
                        for e = 1 : size(this.rec2sat.avail_index,1)
                            [gps_week, gps_sow, gps_dow] = this.time.getGpsWeek(e);
                            error_tropo(e,s) = Atmosphere.saastamoinen_model_GPT(lat(e), lon(e), az(e), el(e), gps_sow, this.cs.iono)
                        end
                        
                end
            end
            
        end
        function error_iono = getErrIono(this,XS)
            error_iono = zeros(size(this.rec2sat.avail_index));
            for s = 1 : size(this.available_index)
                %%% compute lat lon
                [~, lat, ~, lon] = cart2geod(this.XR(:,1), this.XR(:,2), this.XR(:,3));
                %%% compute az el
                if size(this.XR,2)>1
                    [az, el] = this.getAzimuthElevation(this.XR(this.avaliable_index(:,s),:) ,XS(this.avaliable_index(:,s),s,:));
                else
                    [az, el] = this.getAzimuthElevation(this.XR, XS(this.avaliable_index(:,s),s,:));
                end
                
                
                switch this.gs.cur_settings.tropo_model
                    case 0 %no model
                        corr = zeros(size(el));
                    case 1 %Geckle and Feen model
                        %corr = simplified_model(lat, lon, az, el, mjd);
                    case 2 %Klobuchar model
                        [week, sow] = time2weektow(zero_time + this.time_tx);
                        error_iono(idx,s) = Atmosphere.klobuchar_model(lat, lon, az, el, sow, this.cs.iono)
                        
                end
            end
        end
        function [az, el] = getAzimuthElevation(this, XS, XR)
            % SYNTAX:
            %   [az, el] = this.getAzimuthElevation(XS)
            %
            % INPUT:
            % XS = positions of satellite [n_epoch x 1]
            % XR = positions of reciever [n_epoch x 1] (optional, non static
            % case)
            % OUTPUT:
            % Az = Azimuths of satellite [n_epoch x 1]
            % El = Elevations of satellite [n_epoch x 1]
            % during time of travel
            % DESCRIPTION:
            %   Compute Azimuth and elevation of the staellite
            n_epoch = size(XS,1);
            if nargin > 2
                if size(XR,1) ~= n_epoch
                    this.log.addError('[ getAzimuthElevation ] Number of satellite positions differ from number of receiver positions');
                    return
                end
            else
                XR = repmat(this.XR(1,:),n_epoch,1);
            end
            
            az = zeros(n_epoch,1); el = zeros(n_epoch,1);
            
            [phi, lam] = cart2geod(XR(:,1), XR(:,2), XR(:,3));
            XSR = XS - XR; %%% sats orbit with origon in receiver
            
            e_unit = [-sin(lam)            cos(lam)           zeros(size(lam))       ]; % East unit vector
            n_unit = [-sin(phi).*cos(lam) -sin(phi).*sin(lam) cos(phi)]; % North unit vector
            u_unit = [ cos(phi).*cos(lam)  cos(phi).*sin(lam) sin(phi)]; % Up unit vector
            
            e = sum(e_unit .* XSR,2);
            n = sum(n_unit .* XSR,2);
            u = sum(u_unit .* XSR,2);
            
            hor_dist = sqrt( e.^2 + n.^2);
            
            zero_idx = hor_dist < 1.e-20;
            
            az(zero_idx) = 0;
            el(zero_idx) = 90;
            
            az(~zero_idx) = atan2d(e(~zero_idx),n(~zero_idx));
            el(~zero_idx) = atan2d(u(~zero_idx),hor_dist(~zero_idx));
            
            
        end
        function [dist, corr] = getRelDistance(this, XS, XR)
            % SYNTAX:
            %   [corr, distSR_corr] = this.getRelDistance(XS, XR);
            %
            % INPUT:
            % XS = positions of satellite [n_epoch x 1]
            % XR = positions of reciever [n_epoch x 1] (optional, non static
            % case)
            %
            % OUTPUT:
            %   corr = relativistic range error correction term (Shapiro delay)
            %   dist = dist
            % DESCRIPTION:
            %   Compute distance from satellite ot reciever considering
            %   (Shapiro delay) - copied from
            %   relativistic_range_error_correction.m
            n_epoch = size(XS,1);
            if nargin > 2
                if size(XR,1) ~= n_epoch
                    this.log.addError('[ getRelDistance ] Number of satellite positions differ from number of receiver positions');
                    return
                end
            else
                XR = repmat(this.XR(1,:),n_epoch,1);
            end
            
            distR = sqrt(sum(XR.^2 ,2));
            distS = sqrt(sum(XS.^2 ,2));
            
            distSR = sqrt(sum((XS-XR).^2 ,2));
            
            
            GM = 3.986005e14;
            
            
            corr = 2*GM/(goGNSS.V_LIGHT^2)*log((distR + distS + distSR)./(distR + distS - distSR));
            
            dist = distSR + corr;
            
        end
    end
    
    methods (Access = private)
        function parseRin2Data(this, txt, lim, eoh)
            % Parse the data part of a RINEX 2 file -  the header must already be parsed
            % SYNTAX: this.parseRin2Data(txt, lim, eoh)
            
            % find all the observation lines
            t_line = find([false(eoh, 1); (txt(lim(eoh+1:end,1) + 2) ~= ' ')' & (txt(lim(eoh+1:end,1) + 3) == ' ')' & lim(eoh+1:end,3) > 25]);
            n_epo = numel(t_line);
            % extract all the epoch lines
            string_time = txt(repmat(lim(t_line,1),1,25) + repmat(1:25, n_epo, 1))';
            % convert the times into a 6 col time
            date = cell2mat(textscan(string_time,'%2f %2f %2f %2f %2f %10.7f'));
            after_70 = (date(:,1) < 70); date(:, 1) = date(:, 1) + 1900 + after_70 * 100; % convert to 4 digits
            % import it as a GPS_Time obj
            this.time = GPS_Time(date, [], this.file.first_epoch.is_gps);
            this.rate = this.time.getRate();
            n_epo = numel(t_line);
            
            % get number of sat per epoch
            this.n_spe = sscanf(txt(repmat(lim(t_line,1),1,3) + repmat(29:31, n_epo, 1))', '%d');
            
            all_sat = [];
            for e = 1 : n_epo
                n_sat = this.n_spe(e);
                sat = serialize(txt(lim(t_line(e),1) + repmat((0 : ceil(this.n_spe(e) / 12) - 1)' * 69, 1, 36) + repmat(32:67, ceil(this.n_spe(e) / 12), 1))')';
                sat = sat(1:n_sat * 3);
                all_sat = [all_sat sat];
            end
            all_sat = reshape(all_sat, 3, numel(all_sat)/3)';
            
            gps_prn = unique(sscanf(all_sat(all_sat(:,1) == 'G', 2 : 3)', '%2d'));
            glo_prn = unique(sscanf(all_sat(all_sat(:,1) == 'R', 2 : 3)', '%2d'));
            gal_prn = unique(sscanf(all_sat(all_sat(:,1) == 'E', 2 : 3)', '%2d'));
            qzs_prn = unique(sscanf(all_sat(all_sat(:,1) == 'J', 2 : 3)', '%2d'));
            bds_prn = unique(sscanf(all_sat(all_sat(:,1) == 'C', 2 : 3)', '%2d'));
            irn_prn = unique(sscanf(all_sat(all_sat(:,1) == 'I', 2 : 3)', '%2d'));
            sbs_prn = unique(sscanf(all_sat(all_sat(:,1) == 'S', 2 : 3)', '%2d'));
            prn = struct('g', gps_prn', 'r', glo_prn', 'e', gal_prn', 'j', qzs_prn', 'c', bds_prn', 'i', irn_prn', 's', sbs_prn');
            
            % update the maximum number of rows to store
            n_obs = numel(prn.g) * numel(this.rin_obs_code.g) / 3 + ...
                numel(prn.r) * numel(this.rin_obs_code.r) / 3 + ...
                numel(prn.e) * numel(this.rin_obs_code.e) / 3 + ...
                numel(prn.j) * numel(this.rin_obs_code.j) / 3 + ...
                numel(prn.c) * numel(this.rin_obs_code.c) / 3 + ...
                numel(prn.i) * numel(this.rin_obs_code.i) / 3 + ...
                numel(prn.s) * numel(this.rin_obs_code.s) / 3;
            
            clear gps_prn glo_prn gal_prn qzs_prn bds_prn irn_prn sbs_prn;
            
            % order of storage
            % sat_system / obs_code / satellite
            sys_c = char(this.cc.sys_c + 32);
            n_ss = numel(sys_c); % number of satellite system
            
            % init datasets
            obs = zeros(n_obs, n_epo);
            
            this.obs_code = [];
            this.prn = [];
            this.system = [];
            this.f_id = [];
            this.wl = [];
            
            for  s = 1 : n_ss
                sys = sys_c(s);
                n_sat = numel(prn.(sys)); % number of satellite system
                this.n_sat = this.n_sat + n_sat;
                n_code = numel(this.rin_obs_code.(sys)) / 3; % number of satellite system
                % transform in n_code x 3
                obs_code = reshape(this.rin_obs_code.(sys), 3, n_code)';
                % replicate obs_code for n_sat
                obs_code = serialize(repmat(obs_code, 1, n_sat)');
                obs_code = reshape(obs_code, 3, numel(obs_code) / 3)';
                
                this.obs_code = [this.obs_code; obs_code];
                prn_ss = repmat(prn.(sys)', n_code, 1);
                this.prn = [this.prn; prn_ss];
                this.system = [this.system repmat(char(sys - 32), 1, size(obs_code, 1))];
                
                f_id = obs_code(:,2);
                ss = this.cc.(char((this.cc.SYS_NAME{s} + 32)));
                [~, f_id] = ismember(f_id, ss.CODE_RIN3_2BAND);
                
                ismember(this.system, this.cc.SYS_C);
                this.f_id = [this.f_id; f_id];
                
                if s == 2
                    wl = ss.L_VEC((max(1, f_id) - 1) * size(ss.L_VEC, 1) + ss.PRN2IDCH(min(prn_ss, ss.N_SAT))');
                    wl(prn_ss > ss.N_SAT) = NaN;
                    wl(f_id == 0) = NaN;
                else
                    wl = ss.L_VEC(max(1, f_id))';
                    wl(f_id == 0) = NaN;
                end
                this.wl = [this.wl; wl];
            end
            
            this.w_bar.createNewBar(' Parsing epochs...');
            this.w_bar.setBarLen(n_epo);
            
            n_ops = numel(this.rin_obs_code.g)/3; % number of observations per satellite
            n_lps = ceil(n_ops / 5); % number of obbservation lines per satellite
            
            mask = repmat('         0.00000',1 ,40);
            data_pos = repmat(logical([true(1, 14) false(1, 2)]),1 ,40);
            id_line  = reshape(1 : numel(mask), 80, numel(mask)/80);
            for e = 1 : n_epo % for each epoch
                n_sat = this.n_spe(e);
                sat = serialize(txt(lim(t_line(e),1) + repmat((0 : ceil(this.n_spe(e) / 12) - 1)' * 69, 1, 36) + repmat(32:67, ceil(this.n_spe(e) / 12), 1))')';
                sat = sat(~isspace(sat));
                sat = sat(1:n_sat * 3);
                sat = reshape(sat, 3, n_sat)';
                prn_e = sscanf(serialize(sat(:,2:3)'), '%02d');
                for s = 1 : size(sat, 1)
                    % line to fill with the current observation line
                    obs_line = (this.prn == prn_e(s)) & this.system' == sat(s, 1);
                    line_start = t_line(e) + ceil(n_sat / 12) + (s-1) * n_lps;
                    line = mask(1 : n_ops * 16);
                    for i = 0 : n_lps - 1
                        try
                            line(id_line(1:lim(line_start + i, 3) + 1,i+1)) = txt(lim(line_start + i, 1) : lim(line_start + i, 2));
                        catch
                            % empty last lines
                        end
                    end
                    % remove return characters
                    ck = line == ' '; line(ck) = mask(ck); % fill empty fields -> otherwise textscan ignore the empty fields
                    % try with sscanf
                    line = line(data_pos(1 : numel(line)));
                    data = sscanf(reshape(line, 14, numel(line) / 14), '%f');
                    obs(obs_line, e) = data;
                    % alternative approach with textscan
                    %data = textscan(line, '%14.3f%1d%1d');
                    %obs(obs_line(1:numel(data{1})), e) = data{1};
                end
                this.w_bar.go(e);
            end
            this.logger.newLine();
            this.obs = obs;
        end
        
        function parseRin3Data(this, txt, lim, eoh)
            % find all the observation lines
            t_line = find([false(eoh, 1); (txt(lim(eoh+1:end,1)) == '>')']);
            n_epo = numel(t_line);
            % extract all the epoch lines
            string_time = txt(repmat(lim(t_line,1),1,27) + repmat(2:28, n_epo, 1))';
            % convert the times into a 6 col time
            date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %10.7f'));
            % import it as a GPS_Time obj
            this.time = GPS_Time(date, [], this.file.first_epoch.is_gps);
            this.rate = this.time.getRate();
            n_epo = numel(t_line);
            
            % get number of observations per epoch
            this.n_spe = sscanf(txt(repmat(lim(t_line,1),1,3) + repmat(32:34, n_epo, 1))', '%d');
            d_line = find(~[true(eoh, 1); (txt(lim(eoh+1:end,1)) == '>')']);
            
            all_sat = txt(repmat(lim(d_line,1), 1, 3) + repmat(0 : 2, numel(d_line), 1));
            
            % find the data present into the file
            gps_line = d_line(txt(lim(d_line,1)) == 'G');
            glo_line = d_line(txt(lim(d_line,1)) == 'R');
            gal_line = d_line(txt(lim(d_line,1)) == 'E');
            qzs_line = d_line(txt(lim(d_line,1)) == 'J');
            bds_line = d_line(txt(lim(d_line,1)) == 'C');
            irn_line = d_line(txt(lim(d_line,1)) == 'I');
            sbs_line = d_line(txt(lim(d_line,1)) == 'S');
            % Activate only the constellation that are present in the receiver
            %this.cc.setActive([isempty(gps_line) isempty(glo_line) isempty(gal_line) isempty(qzs_line) isempty(bds_line) isempty(irn_line) isempty(sbs_line)]);
            
            gps_prn = unique(sscanf(txt(repmat(lim(gps_line,1), 1, 2) + repmat(1 : 2, numel(gps_line), 1))', '%2d'));
            glo_prn = unique(sscanf(txt(repmat(lim(glo_line,1), 1, 2) + repmat(1 : 2, numel(glo_line), 1))', '%2d'));
            gal_prn = unique(sscanf(txt(repmat(lim(gal_line,1), 1, 2) + repmat(1 : 2, numel(gal_line), 1))', '%2d'));
            qzs_prn = unique(sscanf(txt(repmat(lim(qzs_line,1), 1, 2) + repmat(1 : 2, numel(qzs_line), 1))', '%2d'));
            bds_prn = unique(sscanf(txt(repmat(lim(bds_line,1), 1, 2) + repmat(1 : 2, numel(bds_line), 1))', '%2d'));
            irn_prn = unique(sscanf(txt(repmat(lim(irn_line,1), 1, 2) + repmat(1 : 2, numel(irn_line), 1))', '%2d'));
            sbs_prn = unique(sscanf(txt(repmat(lim(sbs_line,1), 1, 2) + repmat(1 : 2, numel(sbs_line), 1))', '%2d'));
            prn = struct('g', gps_prn', 'r', glo_prn', 'e', gal_prn', 'j', qzs_prn', 'c', bds_prn', 'i', irn_prn', 's', sbs_prn');
            
            % update the maximum number of rows to store
            n_obs = numel(prn.g) * numel(this.rin_obs_code.g) / 3 + ...
                numel(prn.r) * numel(this.rin_obs_code.r) / 3 + ...
                numel(prn.e) * numel(this.rin_obs_code.e) / 3 + ...
                numel(prn.j) * numel(this.rin_obs_code.j) / 3 + ...
                numel(prn.c) * numel(this.rin_obs_code.c) / 3 + ...
                numel(prn.i) * numel(this.rin_obs_code.i) / 3 + ...
                numel(prn.s) * numel(this.rin_obs_code.s) / 3;
            
            clear gps_prn glo_prn gal_prn qzs_prn bds_prn irn_prn sbs_prn;
            
            % order of storage
            % sat_system / obs_code / satellite
            sys_c = char(this.cc.sys_c + 32);
            n_ss = numel(sys_c); % number of satellite system
            
            % init datasets
            obs = zeros(n_obs, n_epo);
            
            this.obs_code = [];
            this.prn = [];
            this.system = [];
            this.f_id = [];
            this.wl = [];
            this.n_sat = 0;
            for  s = 1 : n_ss
                sys = sys_c(s);
                n_sat = numel(prn.(sys)); % number of satellite system
                this.n_sat = this.n_sat + n_sat;
                n_code = numel(this.rin_obs_code.(sys)) / 3; % number of satellite system
                % transform in n_code x 3
                obs_code = reshape(this.rin_obs_code.(sys), 3, n_code)';
                % replicate obs_code for n_sat
                obs_code = serialize(repmat(obs_code, 1, n_sat)');
                obs_code = reshape(obs_code, 3, numel(obs_code) / 3)';
                
                this.obs_code = [this.obs_code; obs_code];
                prn_ss = repmat(prn.(sys)', n_code, 1);
                this.prn = [this.prn; prn_ss];
                this.system = [this.system repmat(char(sys - 32), 1, size(obs_code, 1))];
                
                f_id = obs_code(:,2);
                ss = this.cc.(char((this.cc.SYS_NAME{s} + 32)));
                [~, f_id] = ismember(f_id, ss.CODE_RIN3_2BAND);
                
                ismember(this.system, this.cc.SYS_C);
                this.f_id = [this.f_id; f_id];
                
                if s == 2
                    wl = ss.L_VEC((max(1, f_id) - 1) * size(ss.L_VEC, 1) + ss.PRN2IDCH(min(prn_ss, ss.N_SAT))');
                    wl(prn_ss > ss.N_SAT) = NaN;
                    wl(f_id == 0) = NaN;
                else
                    wl = ss.L_VEC(max(1, f_id))';
                    wl(f_id == 0) = NaN;
                end
                if sum(f_id == 0)
                    [~, id] = unique(double(obs_code(f_id == 0, :)) * [1 10 100]');
                    this.logger.addWarning(sprintf('These codes for the %s are not recognized, ignoring data: %s', ss.SYS_EXT_NAME, sprintf('%c%c%c ', obs_code(id, :)')));
                end
                this.wl = [this.wl; wl];
            end
            
            this.w_bar.createNewBar(' Parsing epochs...');
            this.w_bar.setBarLen(n_epo);
            
            mask = repmat('         0.00000',1 ,40);
            data_pos = repmat(logical([true(1, 14) false(1, 2)]),1 ,40);
            for e = 1 : n_epo % for each epoch
                sat = txt(repmat(lim(t_line(e) + 1 : t_line(e) + this.n_spe(e),1),1,3) + repmat(0:2, this.n_spe(e), 1));
                prn_e = sscanf(serialize(sat(:,2:3)'), '%02d');
                for s = 1 : size(sat, 1)
                    % line to fill with the current observation line
                    obs_line = find((this.prn == prn_e(s)) & this.system' == sat(s, 1));
                    line = txt(lim(t_line(e) + s, 1) + 3 : lim(t_line(e) + s, 2));
                    ck = line == ' '; line(ck) = mask(ck); % fill empty fields -> otherwise textscan ignore the empty fields
                    % try with sscanf
                    line = line(data_pos(1 : numel(line)));
                    data = sscanf(reshape(line, 14, numel(line) / 14), '%f');
                    obs(obs_line(1:size(data,1)), e) = data;
                    % alternative approach with textscan
                    %data = textscan(line, '%14.3f%1d%1d');
                    %obs(obs_line(1:numel(data{1})), e) = data{1};
                end
                this.w_bar.go(e);
            end
            this.logger.newLine();
            this.obs = obs;
            
        end
    end
    
end
