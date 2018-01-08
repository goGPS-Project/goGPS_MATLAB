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
%    |___/                    v 0.6.0 alpha 1 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, Giulio Tagliaferro ...
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
classdef Receiver < Exportable_Object
        
    properties (SetAccess = public, GetAccess = public)
        file           % file rinex object
        rin_type       % rinex version format
        rinex_ss       % flag containing the satellite system of the rinex file, G: GPS, R: GLONASS, E: Galileo, J: QZSS, C: BDS, I: IRNSS, S: SBAS payload, M: Mixed
        
        name           % marker name
        type           % marker type
                
        rid            % receiver interobservation biases
        flag_rid       % clock error for each obs code {num_obs_code}
        
        static         % static or dynamic receiver 1: static 0: dynamic
        
        n_sat = 0;     % number of satellites
        n_freq = 0;    % number of stored frequencies
        n_spe = [];    % number of observations per epoch                
    end
    
    % ==================================================================================================================================================
    %  OBSERVATIONS
    % ==================================================================================================================================================

    properties (SetAccess = public, GetAccess = public)
        active_ids     % rows of active satellites
        wl             % wave-lenght of each row of row_id
        f_id           % frequency number e.g. L1 -> 1,  L2 ->2, E1 -> 1, E5b -> 3 ...
        prn            % pseudo-range number of the satellite
        go_id          % internal id for a certain satellite
        system         % char id of the satellite system corresponding to the row_id
                
        %obs_validity   % validity of the row (does it contains usable values?)
        
        obs_code       % obs code for each line of the data matrix obs
        obs            % huge obbservation matrix with all the observables for all the systems / frequencies / ecc ...
        synt_ph;       % syntetic phases
                
        % ANTENNA ----------------------------------

        ant            % antenna number
        ant_type       % antenna type
        ant_delta_h    % antenna height from the ground [m]
        ant_delta_en   % antenna east/north offset from the ground [m]
        
        % CORRECTIONS ------------------------------
        
        rin_obs_code   % list of types per constellation
        ph_shift       % phase shift as read from RINEX files

        ocean_load_disp = [];        % ocean loading displacemnet for the station
        clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true        
        pcv = []       % phase center corrections for the receiver
        
        group_delay_status = 0; % flag to indicate if code measurement have been corrected using group delays                          (0: not corrected , 1: corrected)
        dts_delay_status   = 0; % flag to indicate if code and phase measurement have been corrected for the clock of the satellite    (0: not corrected , 1: corrected)
        sh_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for shapiro delay                 (0: not corrected , 1: corrected)
        pcv_delay_status   = 0; % flag to indicate if code and phase measurement have been corrected for pcv variations                (0: not corrected , 1: corrected)
        ol_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for ocean loading                 (0: not corrected , 1: corrected)
        pt_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for pole tides                    (0: not corrected , 1: corrected)
        pw_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for phase wind up                 (0: not corrected , 1: corrected)
        et_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for solid earth tide              (0: not corrected , 1: corrected)

        % FLAGS ------------------------------
        
        outlier_idx_ph;
        cycle_slip_idx_ph; % 1 found not repaired , -1 found repaired
        ph_idx             % idx of outlier and cycle slip in obs 
    end
    
    % ==================================================================================================================================================
    %  CELESTIAL INFORMATIONS
    % ==================================================================================================================================================

    properties (SetAccess = public, GetAccess = public)
        sat = struct( ...
            'avail_index',      [], ...    % boolean [n_epoch x n_sat] availability of satellites
            'err_tropo',        [], ...    % double  [n_epoch x n_sat] tropo error
            'err_iono',         [], ...    % double  [n_epoch x n_sat] iono error
            'solid_earth_corr', [], ...    % double  [n_epoch x n_sat] solid earth corrections
            'dtS',              [], ...    % double  [n_epoch x n_sat] staellite clok error at trasmission time
            'rel_clk_corr',     [], ...    % double  [n_epoch x n_sat] relativistic correction at trasmission time
            'tot',              [], ...    % double  [n_epoch x n_sat] time of travel
            'az',               [], ...    % double  [n_epoch x n_sat] azimuth
            'el',               [], ...    % double  [n_epoch x n_sat] elevation
            'cs',               [], ...    % Core_Sky
            'XS_tx',            [], ...    % compute Satellite postion a t transmission time
            'crx',              [], ...    % bad epochs based on crx file
            'slant_td',         [] ...    % slant total delay
            )
    end
   
    % ==================================================================================================================================================
    %  POSITION
    % ==================================================================================================================================================

    properties (SetAccess = public, GetAccess = public)
        xyz_approx     % approximate position of the receiver (XYZ geocentric)
        xyz            % position of the receiver (XYZ geocentric)
        enu            % position of the receiver (ENU local UTM)
        
        lat            % ellipsoidal latitude
        lon            % ellipsoidal longitude
        h_ellips       % ellipsoidal height
        h_ortho        % orthometric height        
    end    
    
    % ==================================================================================================================================================
    %  TIME
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        time = [];     % internal time ref of the stored epochs
        rate;          % observations rate;
        desync         % receiver clock desync (difference between nominal time and the time of the observation)
        dt_ph          % clock error for the phases generated by the de-sync process
        dt_pr          % clock error for the pseudo-ranges generated by the de-sync process
        dt_ip          % clock error correction estimated during init positioning
        dt             % reference clock error of the receiver [n_epochs x num_obs_code]
    end
    
    % ==================================================================================================================================================
    %  TROPO
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        zhd      % zenital hydrostatic delay           double   [n_epoch x 1]
        ztd      % zenital tropospheric delay          double   [n_epoch x 1]
        zwd      % zenital wet delay                   double   [n_epoch x 1]
        pwv      % precipitable water vapour           double   [n_epoch x 1]
        
        tgn      % tropospheric gradient north         double   [n_epoch x n_sat]
        tge      % tropospheric gradient east          double   [n_epoch x n_sat]
    end
    
    % ==================================================================================================================================================
    %  QUALITY INDEXES
    % ==================================================================================================================================================
    
    properties
        hdop
        khdop
        a_fix
        s_rate
    end

    % ==================================================================================================================================================
    %  USEFUL HANDLES
    % ==================================================================================================================================================
    
    properties (SetAccess = private, GetAccess = public)
        cc = Constellation_Collector('GRECJ'); % local cc
        w_bar                                  % handle to waitbar
        state                                  % local handle of state;
        log                                    % handle to log
    end

    % ==================================================================================================================================================
    %  PUBLIC METHODS
    % ==================================================================================================================================================
    
    methods
        function this = Receiver(cc, rinex_file_name)
            % SYNTAX  this = Receiver(<cc>, <rinex_file_name>)
            this.initObs();
            this.log = Logger.getInstance();
            this.state = Go_State.getCurrentSettings();
            if nargin >= 1 && ~isempty(cc)
                this.cc = cc;
            else
                this.cc = this.state.cc;
            end
            this.w_bar = Go_Wait_Bar.getInstance();
            if nargin >= 2 && ~isempty(rinex_file_name)
                this.loadRinex(rinex_file_name);
                this.loadAntModel();
            end  
        end
        
        function toString(this)
            % Display on screen information about the receiver
            % SYNTAX: this.toString();
            
            fprintf('----------------------------------------------------------------------------------\n')
            this.log.addMarkedMessage(sprintf('Receiver %s', this.name));
            fprintf('----------------------------------------------------------------------------------\n')
            this.log.addMessage(sprintf(' From     %s', this.time.first.toString()));
            this.log.addMessage(sprintf(' to       %s', this.time.last.toString()));
            this.log.newLine();
            this.log.addMessage(sprintf(' Rate of the observations [s]:            %d', this.rate));
            this.log.newLine();
            this.log.addMessage(sprintf(' Maximum number of satellites seen:       %d', this.n_sat));
            this.log.addMessage(sprintf(' Number of stored frequencies:            %d', this.n_freq));
            this.log.newLine();
            this.log.addMessage(sprintf(' Satellite System(s) seen:                "%s"', unique(this.system)));
            this.log.newLine();
            
            xyz0 = this.getAPrioriPos();
            [enu0(1), enu0(2), enu0(3)] = cart2plan(xyz0(:,1), xyz0(:,2), xyz0(:,3));
            static_dynamic = {'Dynamic', 'Static'};
            this.log.addMessage(sprintf(' %s receiver', static_dynamic{this.static + 1}));
            fprintf(' ----------------------------------------------------------\n')
            this.log.addMessage(' Receiver a-priori position:');
            this.log.addMessage(sprintf('     X = %+16.4f m        E = %+16.4f m\n     Y = %+16.4f m        N = %+16.4f m\n     Z = %+16.4f m        U = %+16.4f m', ...
                xyz0(1), enu0(1), xyz0(2), enu0(2), xyz0(3), enu0(3)));
            
            if ~isempty(this.xyz)
                enu = zero2nan(this.xyz); [enu(:, 1), enu(:, 2), enu(:, 3)] = cart2plan(zero2nan(this.xyz(:,1)), zero2nan(this.xyz(:,2)), zero2nan(this.xyz(:,3)));
                xyz_m = median(zero2nan(this.xyz), 1, 'omitnan');
                enu_m = median(enu, 1, 'omitnan');
                this.log.newLine();
                this.log.addMessage(' Receiver median position:');
                this.log.addMessage(sprintf('     X = %+16.4f m        E = %+16.4f m\n     Y = %+16.4f m        N = %+16.4f m\n     Z = %+16.4f m        U = %+16.4f m', ...
                                            xyz_m(1), enu_m(1), xyz_m(2), enu_m(2), xyz_m(3), enu_m(3)));
                
                enu = zero2nan(this.xyz); [enu(:, 1), enu(:, 2), enu(:, 3)] = cart2plan(zero2nan(this.xyz(:,1)), zero2nan(this.xyz(:,2)), zero2nan(this.xyz(:,3)));
                xyz_m = median(zero2nan(this.xyz), 1, 'omitnan');
                enu_m = median(enu, 1, 'omitnan');
                this.log.newLine();
                this.log.addMessage(' Correction of the a-priori position:');
                this.log.addMessage(sprintf('     X = %+16.4f m        E = %+16.4f m\n     Y = %+16.4f m        N = %+16.4f m\n     Z = %+16.4f m        U = %+16.4f m', ...
                                            xyz0(1) - xyz_m(1), enu0(1) - enu_m(1), xyz0(2) - xyz_m(2), enu0(2) - enu_m(2), xyz0(3) - xyz_m(3), enu0(3) - enu_m(3)));
            end
            fprintf(' ----------------------------------------------------------\n')
            
        end
        
        function reset(this)
            this.time = GPS_Time();
            this.enu = [];
            this.lat = [];
            this.lon = [];
            
            this.h_ellips = [];
            this.h_ortho = [];
            
            this.n_sat = [];
            this.hdop =  [];
            this.khdop = [];
            this.a_fix = [];
            this.s_rate = [];
            
            this.initObs;
            this.xyz = [];
            
            this.zhd  = [];
            this.zwd  = [];
            this.ztd  = [];
            this.pwv  = [];
            
            this.tgn = [];
            this.tge = [];
            
            this.sat = struct( ...
                'avail_index',      [], ...    % boolean [n_epoch x n_sat] availability of satellites
                'err_tropo',        [], ...    % double  [n_epoch x n_sat] tropo error
                'err_iono',         [], ...    % double  [n_epoch x n_sat] iono error
                'solid_earth_corr', [], ...    % double  [n_epoch x n_sat] solid earth corrections
                'dtS',              [], ...    % double  [n_epoch x n_sat] staellite clok error at trasmission time
                'rel_clk_corr',     [], ...    % double  [n_epoch x n_sat] relativistic correction at trasmission time
                'tot',              [], ...    % double  [n_epoch x n_sat] time of travel
                'az',               [], ...    % double  [n_epoch x n_sat] azimuth
                'el',               [], ...    % double  [n_epoch x n_sat] elevation
                'cs',               [], ...    % Core_Sky
                'XS_tx',            [], ...    % compute Satellite postion a t transmission time
                'crx',              [], ...    % bad epochs based on crx file
                'slant_td',         [] ...     % slant total delay
                );
            this.initR2S;
        end
        
        function initObs(this)
            % Reset the content of the receiver obj
            % SYNTAX:
            %   this.initObs;
            
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
            
            this.static       = true;     % the receivers are considered static unless specified
            
            this.n_sat = 0;               % number of satellites
            this.n_freq = 0;              % number of stored frequencies
            this.n_spe = [];              % number of sat per epoch
            
            this.rate = 0;                % observations rate;

            this.desync = 0;              % clock offset of the receiver
            this.dt_ph = 0;               % clock offset of the receiver
            this.dt_pr = 0;               % clock offset of the receiver
            this.dt_ip = 0;               % clock offset of the receiver
            this.dt = 0;                  % clock offset of the receiver
            this.flag_rid = 0;         % clock offset of the receiver
                        
            this.active_ids = [];         % rows of active satellites
            this.wl         = [];         % wave-lenght of each row of row_id
            this.f_id       = [];         % frequency number e.g. L1 -> 1,  L2 ->2, E1 -> 1, E5b -> 3 ...
            this.prn        = [];         % pseudo-range number of the satellite
            this.go_id      = [];         % internal id for a certain satellite
            this.system     = '';         % char id of the satellite system corresponding to the row_id
                        
            this.obs_code   = [];         % obs code for each line of the data matrix obs
            this.obs        = [];         % huge obbservation matrix with all the observables for all the systems / frequencies / ecc ...
            
            this.clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true
            
            this.initR2S();
        end
        
        function initR2S(this)
            % initialize satellite related parameters
            % SYNTAX: this.initR2S();
            
            this.sat.cs           = Core_Sky.getInstance();
            %this.sat.avail_index  = false(this.getNumEpochs, this.cc.getNumSat);
            %  this.sat.XS_tx     = NaN(n_epoch, n_pr); % --> consider what to initialize
        end
                
        function loadRinex(this, file_name)
            % Parses RINEX observation files.
            %
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
            
            t0 = tic;
            
            this.log.addMarkedMessage('Reading observations...');
            this.log.newLine();
            
            this.file =  File_Rinex(file_name, 9);
            
            if this.file.isValid()
                this.log.addMessage(sprintf('Opening file %s for reading', file_name), 100);
                % open RINEX observation file
                fid = fopen(file_name,'r');
                txt = fread(fid,'*char')';
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
                
                % guess rinex3 flag for incomplete flag (probably coming from rinex2 or converted rinex2 -> rinex3)
                % WARNING!! (C/A) + (P2-P1) semi codeless tracking (flag C2D) receiver not supporter (in rinex 2) convert them
                % using cc2noncc converter https://github.com/ianmartin/cc2noncc (not tested)
                
                % GPS C1 -> C1C
                idx = this.getObsIdx('C1 ','G');
                this.obs_code(idx,:) = repmat('C1C',length(idx),1);
                % GPS C2 -> C2C
                idx = this.getObsIdx('C2 ','G');
                this.obs_code(idx,:) = repmat('C2C',length(idx),1);
                % GPS C5 -> C5I
                idx = this.getObsIdx('C5 ','G');
                this.obs_code(idx,:) = repmat('C5I',length(idx),1);
                % GPS P1 -> C1W
                idx = this.getObsIdx('P1 ','G');
                this.obs_code(idx,:) = repmat('C1W',length(idx),1);
                % GPS P2 -> C2W
                idx = this.getObsIdx('P2 ','G');
                this.obs_code(idx,:) = repmat('C2W',length(idx),1);
                % GLONASS C1 -> C1C
                idx = this.getObsIdx('C1 ','R');
                this.obs_code(idx,:) = repmat('C1C',length(idx),1);
                % GLONASS C2 -> C2C
                idx = this.getObsIdx('C2 ','R');
                this.obs_code(idx,:) = repmat('C2C',length(idx),1);
                % GLONASS P1 -> C1P
                idx = this.getObsIdx('P1 ','R');
                this.obs_code(idx,:) = repmat('C1P',length(idx),1);
                % GLONASS P2 -> C2P
                idx = this.getObsIdx('P2 ','R');
                this.obs_code(idx,:) = repmat('C2P',length(idx),1);
                % other flags to be investiagated
                
                this.log.addMessage(sprintf('Parsing completed in %.2f seconds', toc(t0)));
                this.log.newLine();
            end
            
            % Compute the other useful status array of the receiver object
            this.updateStatus();
            this.active_ids = true(this.getNumObservables, 1);
            
            % remove empty observables
            this.remObs(~this.active_ids)
        end
        
        function loadAntModel(this)
            % Load and arse the antenna (ATX) file as specified into settings
            % SYNTAX: 
            %   this.loadAntModel
            filename_pcv = this.state.getAtxFile;
            fnp = File_Name_Processor();
            this.log.addMessage(sprintf('      Opening file %s for reading', fnp.getFileName(filename_pcv)));
            this.pcv = read_antenna_PCV(filename_pcv, {this.ant_type});
        end
        
        function preProcessing(this, is_static)
            % Do all operation needed in order to preprocess the data
            % remove bad observation (spare satellites or bad epochs from CRX)
            % SYNTAX:
            %   this.preProcessing();
            
            this.remBad();
            % correct for raw estimate of clock error based on the phase measure
            this.correctTimeDesync();
            % set to static or dynamic
            if (nargin == 2)
                this.static = is_static;
            else
                this.static = this.state.kf_mode == 0;
            end
            % this.TEST_smoothCodeWithDoppler(51);
            % code only solution
            this.initPositioning();
            % smooth clock estimation
            this.smoothAndApplyDt();
            this.dt_ip = this.dt;
            this.dt(:) = 0;
            % set all availability index
            this.updateAllAvailIndex();
            % update azimuth elevation
            this.updateAzimuthElevation();
            % Add a model correction for time desync -> observations are now referred to nominal time            
            this.shiftToNominal()
            
            % apply various corrections
            this.sat.cs.toCOM(); %interpolation of attitude with 15min coordinate might possibly be inaccurate switch to COM
            this.applyPCV();
            this.applyPoleTide();
            this.applyPhaseWindUpCorr();
            this.applySolidEarthTide();
            this.applyShDelay();
            this.applyOceanLoading();
            this.removeOutlierMarkCycleSlip();            
        end
        
        function TEST_smoothCodeWithDoppler(this, win_size)
            % This function has been tested in particular cases on UBLOX single frequency
            % In the future see: Optimal Doppler-aided smoothing strategy for GNSS navigation
            [pr, id_pr] = this.getPseudoRanges;
            pr_corr = Core_Pre_Processing.diffAndPred(pr + cumsum(nan2zero(this.getDoppler * this.rate)));
            for s = 1 : size(pr_corr, 2)
                pr_corr(:,s) = cumsum(nan2zero(pr_corr(:,s) - splinerMat([], movmedian(pr_corr(:,s), 3, 'omitnan'), win_size, 1e-9)));
                pr_corr(:,s) = pr_corr(:,s) - splinerMat([], pr_corr(:,s), win_size, 1e-9);
            end
            this.setPseudoRanges(this.getPseudoRanges - pr_corr, id_pr);
        end
        
        function updateStatus(this)
            % Compute the other useful status array of the receiver object
            %
            % SYNTAX:
            %   this.updateStatus();
            
            [~, ss_id] = ismember(this.system, this.cc.sys_c);
            this.n_freq = numel(unique(this.f_id));
            ss_offset = cumsum([0 this.cc.n_sat(1:end-1)]);
            this.go_id = this.prn + reshape(ss_offset(ss_id'),length(this.prn),1); %%% some time second vector is a colum some time is a line reshape added to uniform
            this.n_sat = numel(unique(this.go_id));
            
            % Compute number of satellite per epoch
            
            % considerig only epoch with code on the first frequency
            code_line = this.obs_code(:,1) == 'C' & this.f_id == 1;
            this.n_spe = sum(this.obs(code_line, :) ~= 0);
            % more generic approach bbut a lot slower
            %for e = 1 : this.getNumEpochs()
            %    this.n_spe(e) = numel(unique(this.go_id(this.obs(:,e) ~= 0)));
            %end
            
            this.active_ids = any(this.obs, 2);
            
            %this.sat.avail_index = false(this.time.length, this.cc.getNumSat());
        end
        
        function remEpoch(this, id_epo)
            % remove epochs with a certain id
            % SYNTAX:   this.remObs(id_obs)
            
            this.obs(:,epo) = 0;
            this.active_ids = any(this.obs, 2);
        end
        
        function remObs(this, id_obs)
            % remove observations with a certain id
            % SYNTAX:   this.remObs(id_obs)
            this.obs(id_obs,:) = [];
            this.obs_code(id_obs, :) = [];
            this.active_ids(id_obs) = [];
            this.wl(id_obs) = [];
            this.f_id(id_obs) = [];
            this.prn(id_obs) = [];
            this.go_id(id_obs) = [];
            this.system(id_obs) = [];
        end
        
        function remEmptyObs(this)
            %DESCRIPTION: remove empty obs lines
            empty_sat = sum(abs(this.obs),2) == 0;
            this.remObs(empty_sat);
        end
        
        function remBad(this)
            % Remove observation marked as bad in crx file and satellites
            % whose prn exceed the maximum prn (spare satellites, in maintenance, etc ..)
            % remove spare satellites
            % SYNTAX:
            %   this.remBad();
            
            this.log.addMarkedMessage('Removing satellites whose prn exceed the maximum official one')
            for s = 1 : length(this.cc.sys_c)
                sys_idx = find(this.system == this.cc.sys_c(s));
                prn = this.prn(sys_idx);
                over_max_idx = prn > this.cc.n_sat(s);
                to_remove = sys_idx(over_max_idx);
                this.obs(to_remove, :) = 0;
            end
            %remove bad epoch in crx
            this.log.addMarkedMessage('Removing observations marked as bad in Bernese .CRX file')
            [CRX, found] = load_crx(this.state.crx_dir, double(this.time.getGpsWeek), this.time.getGpsTime, this.cc);
            if found
                for s = 1 : size(CRX,1)
                    c_sat_idx = this.go_id == s;
                    this.obs(c_sat_idx,CRX(s,:)) = 0;
                end
            end
            % check empty lines
            this.remEmptyObs();
        end
        
        function nominal_time = getNominalTime(this)
            % get the nominal time aka rounded time cosidering a  constant
            % sampling rate
            % SYNTAX: nominal_time = this.getNominalTime()
            nominal_time_zero = floor(this.time.first.getMatlabTime() * 24)/24;
            rinex_time = this.time.getRefTime(nominal_time_zero);
            nominal_time = round(rinex_time * this.time.getRate) / this.time.getRate;
            ref_time = (nominal_time(1) : this.time.getRate : nominal_time(end))';
            
            % reordering observations filling empty epochs with zeros;
            nominal_time = GPS_Time(nominal_time_zero, ref_time, this.time.isGPS(), 2);
            nominal_time.toUnixTime;
        end
        
        function correctTimeDesync(this)
            %   Correction of jumps in code and phase due to dtR and time de-sync
            % SYNTAX:
            %   this.correctTimeDesync()
            this.log.addMarkedMessage('Correct for time desync');
            
            % computing nominal_time
            nominal_time_zero = floor(this.time.first.getMatlabTime() * 24)/24;
            rinex_time = this.time.getRefTime(nominal_time_zero);
            nominal_time = round(rinex_time * this.time.getRate) / this.time.getRate;
            ref_time = (nominal_time(1) : this.time.getRate : nominal_time(end))';
            
            % reordering observations filling empty epochs with zeros;
            this.time = GPS_Time(nominal_time_zero, ref_time, this.time.isGPS(), 2);
            this.time.toUnixTime;
            
            [~, id_not_empty] = intersect(ref_time, nominal_time);
            id_empty = setdiff(1 : numel(nominal_time), id_not_empty);
            
            time_desync = nan(size(nominal_time));
            time_desync(id_not_empty) = round((rinex_time - nominal_time) * 1e7) / 1e7; % the rinex time has a maximum of 7 significant decimal digits
            time_desync = simpleFill1D(time_desync,isnan(time_desync),'nearest');            
            
            this.desync = time_desync;
            clear nominal_time_zero nominal_time rinex_time
            
            this.obs(:, id_not_empty) = this.obs;
            this.obs(:, id_empty) = 0;
            this.n_spe(id_not_empty) = this.n_spe;
            this.n_spe(id_empty) = 0;
            
            % extract observations
            [ph, wl, id_ph] = this.getPhases();
            [pr, id_pr] = this.getPseudoRanges();
            
            % apply desync                        
            if any(time_desync)
                [ph_dj, dt_ph_dj] = Core_Pre_Processing.remDtJumps(ph);
                [pr_dj, dt_pr_dj] = Core_Pre_Processing.remDtJumps(pr);                
                ddt_pr = Core_Pre_Processing.diffAndPred(dt_pr_dj);                
                
                %% time_desync is a introduced by the receiver to maintain the drift of the clock into a certain range
                ddt = [0; diff(time_desync)];
                ddrifting = ddt - ddt_pr;
                drifting = cumsum(ddt - ddt_pr);
                
                % Linear interpolation of ddrifting
                jmp_reset = find(abs(ddt_pr) > 1e-7); % points where the clock is reset
                jmp_fit = setdiff(find(abs(ddrifting) > 1e-7), jmp_reset); % points where desync interpolate the clock                
                d_points = [drifting(jmp_reset); drifting(jmp_fit) - ddrifting(jmp_fit)/2];
                jmp = [jmp_reset; jmp_fit];
                drifting = interp1(jmp, d_points, (1 : numel(drifting))', 'spline');
                
                dt_ph = drifting + dt_ph_dj;
                dt_pr = drifting + dt_pr_dj;
                
                t_offset = round(mean(dt_pr(jmp) - time_desync(jmp) + ddrifting(jmp)/2) * 1e7) * 1e-7;
                dt_ph = dt_ph - t_offset;
                dt_pr = dt_pr - t_offset;
                
                ph = bsxfun(@minus, ph, dt_ph .* 299792458);
                pr = bsxfun(@minus, pr, dt_pr .* 299792458);                
            else
                [ph_dj, dt_ph_dj] = Core_Pre_Processing.remDtJumps(ph);                
                [pr_dj, dt_pr_dj] = Core_Pre_Processing.remDtJumps(pr);
                ddt_pr = Core_Pre_Processing.diffAndPred(dt_pr_dj);
                jmp_reset = find(abs(ddt_pr) > 1e-7); % points where the clock is reset
                if numel(jmp_reset) > 2
                    drifting = interp1(jmp_reset, dt_pr_dj(jmp_reset), (1 : numel(ddt_pr))', 'spline');
                else
                    drifting = 0;
                end
                dt_pr = dt_pr_dj - drifting;
                t_offset = mean(dt_pr);
                dt_ph = dt_ph_dj - drifting - t_offset;
                dt_pr = dt_pr - t_offset;

                ph = bsxfun(@minus, ph, dt_ph .* 299792458);
                pr = bsxfun(@minus, pr, dt_pr .* 299792458);                
            end
            
            if any(dt_ph_dj)
                this.log.addMessage(this.log.indent('Correcting carrier phases jumps', 6));
            else
                this.log.addMessage(this.log.indent('Correcting carrier phases for a dt drift estimated from desync interpolation', 6));
            end
            if any(dt_pr_dj)
                this.log.addMessage(this.log.indent('Correcting pseudo-ranges jumps', 6));
            else
                this.log.addMessage(this.log.indent('Correcting pseudo-ranges for a dt drift estimated from desync interpolation', 6));
            end
       
            % Saving dt into the object properties
            this.dt_ph = dt_ph;
            this.dt_pr = dt_pr;
            
            % Outlier rejection
            if (this.state.isOutlierRejectionOn())
                this.log.addMarkedMessage('Removing main outliers');
                [ph, flag_ph] = Core_Pre_Processing.flagRawObsD4(ph, ref_time - dt_ph, ref_time, 6, 5); % The minimum threshold (5 - the last parameter) is needed for low cost receiver that are applying dt corrections to the data - e.g. UBX8
                [pr, flag_pr] = Core_Pre_Processing.flagRawObsD4(pr, ref_time - dt_pr, ref_time, 6, 5); % The minimum threshold (5 - the last parameter) is needed for low cost receiver that are applying dt corrections to the data - e.g. UBX8
            end
            
            % Saving observations into the object properties
            this.setPhases(ph, wl, id_ph);
            this.setPseudoRanges(pr, id_pr);
            
            this.time.addSeconds(time_desync - this.dt_pr);
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
            this.rinex_ss = dataset{4};
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
                this.ant_type = strtrim(txt(lim(fln, 1) + (20:40)));
            end
            % 7) 'APPROX POSITION XYZ'
            fln = find(line2head == 7, 1, 'first'); % get field line
            if isempty(fln)
                this.xyz_approx = [0 0 0];
                this.xyz = [0 0 0];
            else
                tmp = sscanf(txt(lim(fln, 1) + (0:41)),'%f')';                                               % read value
                this.xyz_approx = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 3), [0 0 0], tmp);          % check value integrity
                this.xyz = this.xyz_approx;
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
                            obs_code_text = txt(lim(fln(l + l_offset), 1) + (7:59));
                            idx_code = true(length(obs_code_text),1);
                            idx_code(4:4:(floor(length(obs_code_text)/4)*4)) = false; % inedx to take only valid columns
                            obs_code_temp = obs_code_text(idx_code');
                            obs_code_temp((ceil(max(find(obs_code_temp ~= ' '))/3)*3 +1 ):end) = []; %delete empty lines at the end
                            this.rin_obs_code.(sys) = [this.rin_obs_code.(sys) obs_code_temp];
                            
                            l_offset = l_offset + 1;
                        end
                        l = l + l_offset;
                    end
                end
                if ~isempty(strfind(this.rin_obs_code.c, '1'))
                    this.rin_obs_code.c(this.rin_obs_code.c == '1') = '2';
                    this.log.addWarning('BeiDou band 1 is now defined as 2 -> Automatically converting the observation codes of the RINEX!');
                end
            end
            % 20) SYS / PHASE SHIFT
            fln = find(line2head == 20); % get field line
            if this.rin_type < 3
                this.ph_shift = struct('g', zeros(numel(this.rin_obs_code.g) / 3, 1));
            else
                this.ph_shift = struct('g',[],'r',[],'e',[],'j',[],'c',[],'i',[],'s',[]);
                for l = 1 : numel(fln)
                    if txt(lim(fln(l), 1)) ~= ' ' % ignoring phase shif only on subset of satellites
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
        
        function removeOutlierMarkCycleSlip(this)
            this.log.addMarkedMessage('Cleaning observations');
            
            % PARAMETRS
            ol_thr = 0.5; % outlier threshold
            cs_thr = 0.5; % CYCLE SLIP THR
            sa_thr = 10;  % short arc threshold
            
            %----------------------------
            % Outlier Detection
            %----------------------------
            
            this.log.addMessage(this.log.indent(sprintf('Removing observations under cut-off (%d degrees)', this.state.cut_off), 6));
            mask = this.sat.el > this.state.cut_off;
            this.obs = this.obs .* mask(:, this.go_id)';
            
            % mark all as outlier and interpolate
            % get observed values
            [ph, wl, id_ph_l] = this.getPhases;
            
            this.log.addMessage(this.log.indent('Detect outlier candidates from residual phase time derivate', 6));
            % first time derivative
            synt_ph = this.getSyntPhases;
            sensor_ph = Core_Pre_Processing.diffAndPred(ph - synt_ph);
            % subtract median (clock error)
            sensor_ph = bsxfun(@minus, sensor_ph, median(sensor_ph, 2, 'omitnan'));
            % divide for wavelenght
            sensor_ph = bsxfun(@rdivide, sensor_ph, wl');
            % outlier when they exceed 0.5 cycle
            poss_out_idx = abs(sensor_ph) > ol_thr;
            % take them off
            ph2 = ph;
            ph2(poss_out_idx) = nan;
            
            %----------------------------
            % Cycle slip detection
            %----------------------------
            
            this.log.addMessage(this.log.indent('Detect cycle slips from residual phase time derivate', 6));
            % join the nan
            sensor_ph_cs = nan(size(sensor_ph));
            for o = 1 : size(ph2,2)
                tmp_ph = ph2(:,o);
                ph_idx = not(isnan(tmp_ph));
                tmp_ph = tmp_ph(ph_idx);
                if ~isempty(tmp_ph)
                    sensor_ph_cs(ph_idx,o) = Core_Pre_Processing.diffAndPred(tmp_ph - synt_ph(ph_idx,o),1);
                end
            end
            
            % subtract median
            sensor_ph_cs2 = bsxfun(@minus, sensor_ph_cs, median(sensor_ph_cs, 2, 'omitnan'));
            % divide for wavelenght
            sensor_ph_cs2 = bsxfun(@rdivide, sensor_ph_cs2, wl');
            
            % find possible cycle slip
            % cycle slip when they exceed threhsold cycle
            poss_slip_idx = abs(sensor_ph_cs2) > cs_thr;
            
            % check if epoch before cycle slip can be restored
            poss_rest = [poss_slip_idx(2:end,:); zeros(1,size(poss_slip_idx,2))];
            poss_rest = poss_rest & poss_out_idx;
            poss_rest_line = sum(poss_rest,2);
            poss_rest_line = poss_rest_line | [false; poss_rest_line(2:end)];
            ph_rest_lines = ph(poss_rest_line,:);
            synt_ph_rest_lines = synt_ph(poss_rest_line,:);
            
            sensor_rst = Core_Pre_Processing.diffAndPred(ph_rest_lines - synt_ph_rest_lines);
            % subtract median
            sensor_rst = bsxfun(@minus, sensor_rst, median(sensor_rst, 2, 'omitnan'));
            % divide for wavelenght
            sensor_rst = bsxfun(@rdivide, sensor_rst, wl');
            for i = 1:size(sensor_rst,2)
                for c = find(poss_rest(:,i))'
                    if ~isempty(c)
                        idx = sum(poss_rest_line(1:c));
                        if abs(sensor_rst(idx,i)) < cs_thr
                            poss_out_idx(c,i) = false; %is not outlier
                            %move 1 step before the cycle slip index
                            poss_slip_idx(c+1,i) = false;
                            poss_slip_idx(c,i) = true;
                        end
                    end
                end
            end
            
            no_out_ph = ph./repmat(this.wl(id_ph_l)',size(ph,1),1);
            no_out_ph(poss_out_idx) = nan;
            this.ph_idx = find(id_ph_l);
            
            % remove too short possible arc
            to_short_idx = flagMerge(poss_slip_idx,sa_thr);
            poss_slip_idx = [diff(to_short_idx) <0 ; false(1,size(to_short_idx,2))];
            to_short_idx(poss_slip_idx) =false;
            poss_out_idx(to_short_idx) = true;
            
            n_out = sum(sum(poss_out_idx));
            this.outlier_idx_ph = sparse(poss_out_idx);
            this.cycle_slip_idx_ph = double(sparse(poss_slip_idx));
            this.log.addMessage(this.log.indent(sprintf(' - %d phase observations marked as outlier',n_out), 6));
            
            %this.removeShortArch(this.state.getMinArc);
        end
        
        function tryCycleSlipRepair(this)
            %----------------------------
            % Cycle slip repair
            %----------------------------
            
            % window used to estimate cycle slip
            % linear time
            lin_time = 900; %single diffrence is linear in 15 minutes
            max_window = 600; %maximum windows allowed (for computational reason)
            win_size = min(max_window,ceil(lin_time / this.rate/2)*2); %force even
            
            poss_out_idx = this.outlier_idx_ph;
            poss_slip_idx = this.cycle_slip_idx_ph;
            
            [ph, ~, id_ph_l] = this.getPhases();
            id_ph = find(id_ph_l);
            synt_ph = this.getSyntPhases();
            d_no_out = (ph - synt_ph )./ repmat(this.wl(id_ph_l)',size(ph,1),1);
            n_ep = this.time.length;
            poss_slip_idx = double(poss_slip_idx);
            n_cycleslip = 0;
            n_repaired = 0;
            jmps = [];
            for p = 1:size(poss_slip_idx,2)
                c_slip = find(poss_slip_idx(:,p)' == 1);
                for c = c_slip
                    n_cycleslip = n_cycleslip + 1;
                    slip_bf = find(poss_slip_idx(max(1,c - win_size /2 ): c-1, p),1,'last');
                    if isempty(slip_bf)
                        slip_bf = 0;
                    end
                    st_wind_idx = max(1,c - win_size /2 + slip_bf);
                    slip_aft = find(poss_slip_idx(c :min( c + win_size / 2,n_ep), p),1,'first');
                    if isempty(slip_aft)
                        slip_aft = 0;
                    end
                    end_wind_idx = min(c + win_size /2 - slip_bf -1,n_ep);
                    jmp_idx =  win_size /2 - slip_bf +1;
                    
                    % find best master before
                    idx_other_l = 1:size(poss_slip_idx,2) ~= p;
                    idx_win_bf = poss_slip_idx(st_wind_idx : c -1 , idx_other_l) ;
                    best_cs = min(idx_win_bf .* repmat([1:(c-st_wind_idx)]',1,sum(idx_other_l))); %find other arc with farerer cycle slip
                    poss_mst_bf =~isnan(d_no_out(st_wind_idx : c -1 , idx_other_l));
                    for j = 1 : length(best_cs)
                        poss_mst_bf(1:best_cs(j),j) = false;
                    end
                    num_us_ep_bf = sum(poss_mst_bf); % number of usable epoch befor
                    data_sorted = sort(num_us_ep_bf,'descend');
                    [~, rnk_bf] = ismember(num_us_ep_bf,data_sorted);
                    
                    % find best master after
                    idx_win_aft = poss_slip_idx(c : end_wind_idx , idx_other_l) ;
                    best_cs = min(idx_win_aft .* repmat([1:(end_wind_idx -c +1)]',1,sum(idx_other_l))); %find other arc with farerer cycle slip
                    poss_mst_aft =~isnan(d_no_out(c : end_wind_idx, idx_other_l));
                    for j = 1 : length(best_cs)
                        poss_mst_aft(1:best_cs(j),j) = false;
                    end
                    num_us_ep_aft = sum(poss_mst_aft); % number of usable epoch befor
                    data_sorted = sort(num_us_ep_aft,'descend');
                    [~, rnk_aft] = ismember(num_us_ep_aft,data_sorted);
                    
                    
                    idx_other = find(idx_other_l);
                    
                    %decide which arc to use
                    bst_1 = find(rnk_bf == 1  & rnk_aft == 1,1,'first');
                    if isempty(bst_1)
                        bst_1 = find(rnk_bf == 1 ,1,'first');
                        bst_2 = find(rnk_aft == 1,1,'first');
                        same_slope = false;
                    else
                        bst_2 = bst_1;
                        same_slope = true;
                    end
                    
                    s_diff = d_no_out(st_wind_idx : end_wind_idx, p) - [d_no_out(st_wind_idx : c -1, idx_other(bst_1)); d_no_out(c : end_wind_idx, idx_other(bst_2))];
                    jmp = nan;
                    if sum(~isnan(s_diff(1:jmp_idx-1))) > 5 & sum(~isnan(s_diff(jmp_idx:end))) > 5
                        jmp = estimateJump(s_diff, jmp_idx,same_slope,'median');
                        jmps= [jmps ; jmp];
                    end
                    
                    %{
                    ph_piv = d_no_out(st_wind_idx : end_wind_idx,p);
                    usable_s = sum(poss_mst_bf) > 0 & sum(poss_mst_aft);
                    ph_other = d_no_out(st_wind_idx : end_wind_idx,idx_other_l);
                    ph_other([~poss_mst_bf ;~poss_mst_aft]) = nan;
                    ph_other(:, ~usable_s) = [];
                    %}
                    s_diff = ph_other - repmat(ph_piv,1,size(ph_other,2));
                    
                    % repair
                    % TO DO half cycle
                    if ~isnan(jmp)
                        this.obs(id_ph(p),c:end) = nan2zero(zero2nan(this.obs(id_ph(p),c:end)) - round(jmp));
                    end
                    if abs(jmp -round(jmp)) < 0.1
                        poss_slip_idx(c, p) = - 1;
                        n_repaired = n_repaired +1;
                    else
                        poss_slip_idx(c, p) =   1;
                    end
                end
            end
            
            this.log.addMessage(this.log.indent(sprintf(' - %d of %d cycle slip repaired',n_repaired,n_cycleslip),6));
            
            this.cycle_slip_idx_ph = poss_slip_idx;
            
            % save index into object
            
        end
        
        function removeShortArch(this, min_arc)
            % removes arch shorter than
            % SYNTAX:
            %   this.removeShortArch()
            if min_arc > 1
                this.log.addMarkedMessage(sprintf('Removing arcs shorter than %d epochs', 1 + 2 * ceil((min_arc - 1)/2)));
                [ph, wl, id_ph]= this.getPhases();
                idx = ~isnan(ph);
                idx_s = flagShrink(idx, ceil((min_arc - 1)/2));
                idx_e = flagExpand(idx_s, ceil((min_arc - 1)/2));
                el_idx = xor(idx,idx_e);
                this.outlier_idx_ph(el_idx) =  true;
                this.cycle_slip_idx_ph(el_idx) = 0;
                ph(this.outlier_idx_ph) = 0;
                this.setPhases(ph, wl, id_ph);
                this.log.addMessage(this.log.indent(sprintf(' - %d observations have been removed', sum(el_idx(:))), 6));
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
        
        function s_name = getShortName(this)
            %DESCRIPTION: return the fisrt 4 character of the filename
            %corresponfing to the receiver short name
            s_name = this.file.file_name_list{1}(1:4);
        end
        
        function is_static = isStatic(this)
            % return true if the receiver is static
            % SYNTAX: is_static = this.isStatic()
            is_static = this.static;
        end
        
        function xyz = getAPrioriPos(this)
            % return apriori position
            % SYNTAX: xyz = this.getAPrioriPos()
            xyz = this.xyz_approx;
            if ~any(xyz) && ~isempty(this.xyz)
                xyz = median(this.xyz, 1);
            end
        end
        
        function xyz = getMedianPosXYZ(this)
            % return the computed median position of the receiver
            %
            % OUTPUT:
            %   xyz     geocentric coordinates
            %
            % SYNTAX: 
            %   xyz = this.getAPrioriPos()
            
            xyz = this.xyz;
            xyz = median(this.xyz, 1);
        end
        
        function [lat, lon, h_ellips, h_ortho] = getMedianPosGeodetic(this)
            % return the computed median position of the receiver
            %
            % OUTPUT:
            %   lat         latitude  [deg]
            %   lon         longitude [deg]
            %   h_ellips    ellipsoidical heigth [m]
            %   h_ortho     orthometric heigth [m]
            %
            % SYNTAX: 
            %   [lat, lon, h_ellips, h_ortho] = this.getMedianPosGeodetic();
            
            xyz = this.xyz;
            xyz = median(this.xyz, 1);
            [lat, lon, h_ellips] = cart2geod(xyz);
            if nargout == 4
                gs = Go_State.getInstance;
                gs.initGeoid();
                ondu = getOrthometricCorr(lat, lon, gs.getRefGeoid());
                h_ortho = h_ellips + ondu;
            end
            lat = lat / pi * 180;
            lon = lon / pi * 180;
        end
        
        function [mfh, mfw] = getSlantMF(this)
            % Get Mapping function for the satellite slant
            % 
            % OUTPUT:
            %   mfh: hydrostatic mapping function
            %   mfw: wet mapping function
            %
            % SYNTAX:
            %   [mfh, mfw] = this.getSlantMF()

            [lat, lon, ~, h_ortho] = this.getMedianPosGeodetic();
            [gmfh, gmfw] = gmf(this.time.first, lat./180*pi, lon./180*pi, h_ortho, (90 - this.sat.el(:))./180*pi);
            mfh = reshape(gmfh, size(this.sat.el, 1), size(this.sat.el, 2));
            mfw = reshape(gmfw, size(this.sat.el, 1), size(this.sat.el, 2));
        end
        
        function sztd = getSlantZTD(this, smooth_win_size)
            % Get the "zenithalized" total delay
            % SYNTAX:
            %   sztd = this.getSlantZenithalizedTotalDelay(<flag_smooth_data = 0>)
            [mfh, mfw] = this.getSlantMF();
            sztd = bsxfun(@plus, (zero2nan(this.sat.slant_td) - bsxfun(@times, mfh, this.zhd)) ./ mfw, this.zhd);            
            sztd(sztd <= 0) = nan;

            if nargin == 2 && smooth_win_size > 0
                t = this.time.getRefTime;
                for s = 1 : size(sztd,2)
                    id_ok = ~isnan(sztd(:, s));
                    if sum(id_ok) > 3
                        lim = getOutliers(id_ok);
                        lim = limMerge(lim, 2*smooth_win_size);
                        
                        lim = [lim(1) lim(end)];
                        for l = 1 : size(lim, 1)
                            if (lim(l, 2) - lim(l, 1) + 1) > 3
                                id_ok = lim(l, 1) : lim(l, 2);
                                sztd(id_ok, s) = splinerMat(t(id_ok), sztd(id_ok, s) - zero2nan(this.ztd(id_ok)), smooth_win_size, 0.05) + zero2nan(this.ztd(id_ok));
                            end
                        end
                    end
                end
            end
        end
        
        function [pos] = getXR(this)
            % get xyz or the same repeated by the number of epoch
            % SYNTAX: [XR] = getXR(this)
            if size(this.xyz, 1) == 1
                pos = repmat(this.xyz,this.time.length, 1);
            else
                pos = this.xyz;
            end
        end
        
        function pr = pr1(this, sys_c)
            % get p_range 1 (Legacy)
            % SYNTAX this.pr1(<flag_valid>, <sys_c>)
            pr = zeros(this.cc.n_sat_tot, size(this.obs, 2));
            
            id_pr = this.obs_code(:, 1) == 'C' & this.obs_code(:, 2) == '1';
            if (nargin == 2)
                id_pr = id_pr & (this.system == sys_c)';
            end
            pr(this.go_id(id_pr), :) = this.obs(id_pr, :);
        end
        
        function pr = pr2(this, sys_c)
            % get p_range 2 (Legacy)
            % SYNTAX this.pr1(<flag_valid>, <sys_c>)
            pr = zeros(this.cc.n_sat_tot, size(this.obs, 2));
            
            id_pr = this.obs_code(:, 1) == 'C' & this.obs_code(:, 2) == '2';
            if (nargin == 2)
                id_pr = id_pr & (this.system == sys_c)';
            end
            pr(this.go_id(id_pr), :) = this.obs(id_pr, :);
        end
        
        function [ph, wl] = ph1(this, sys_c)
            % get phase 1 (Legacy)
            % SYNTAX this.ph1(<flag_valid>, <sys_c>)
            ph = zeros(this.cc.n_sat_tot, size(this.obs, 2));
            wl = zeros(this.cc.n_sat_tot, 1);
            id_ph = this.obs_code(:, 1) == 'L' & this.obs_code(:, 2) == '1';
            if (nargin == 2)
                id_ph = id_ph & (this.system == sys_c)';
            end
            ph(this.go_id(id_ph), :) = this.obs(id_ph, :);
            wl(this.go_id(id_ph)) = this.wl(id_ph);
            ph = bsxfun(@times, zero2nan(ph), wl);
        end
        
        function [ph, wl] = ph2(this, sys_c)
            % get phase 2 (Legacy)
            % SYNTAX this.ph1(<flag_valid>, <sys_c>)
            ph = zeros(this.cc.n_sat_tot, size(this.obs, 2));
            wl = zeros(this.cc.n_sat_tot, 1);
            id_ph = this.obs_code(:, 1) == 'L' & this.obs_code(:, 2) == '2';
            if (nargin == 2)
                id_ph = id_ph & (this.system == sys_c)';
            end
            wl(this.go_id(id_ph)) = this.wl(id_ph);
            ph = bsxfun(@times, zero2nan(ph), wl);
        end
        
        function [ph, wl, id_ph] = getPhases(this, sys_c)
            % get the phases observations in meter (not cycles)
            % SYNTAX: [ph, wl, id_ph] = this.getPhases(<sys_c>)
            % SEE ALSO: setPhases getPseudoRanges setPseudoRanges
            
            id_ph = this.obs_code(:, 1) == 'L';
            if (nargin == 2)
                id_ph = id_ph & (this.system == sys_c)';
            end
            ph = this.obs(id_ph, :);
            ph(this.outlier_idx_ph') = nan;
            wl = this.wl(id_ph);
            
            ph = bsxfun(@times, zero2nan(ph), wl)';
        end
        
        function [pr, id_pr] = getPseudoRanges(this, sys_c)
            % get the pseudo ranges observations in meter (not cycles)
            % SYNTAX: [pr, id_pr] = this.getPseudoRanges(<sys_c>)
            % SEE ALSO: getPhases setPhases setPseudoRanges
            
            id_pr = this.obs_code(:, 1) == 'C';
            if (nargin == 2)
                id_pr = id_pr & (this.system == sys_c)';
            end
            pr = zero2nan(this.obs(id_pr, :)');
        end        
        
        function [dop, wl, id_dop] = getDoppler(this, sys_c)
            % get the doppler observations
            % SYNTAX: [dop, id_dop] = this.getDoppler(<sys_c>)
            % SEE ALSO: setDoppler
            
            id_dop = this.obs_code(:, 1) == 'D';
            if (nargin == 2)
                id_dop = id_dop & (this.system == sys_c)';
            end
            dop = zero2nan(this.obs(id_dop, :)');
            wl = this.wl(id_dop);
            dop = bsxfun(@times, zero2nan(dop), wl');
        end
        
        function [snr, id_snr] = getSNR(this, sys_c)
            % get the doppler observations
            % SYNTAX: [dop, id_dop] = this.getDoppler(<sys_c>)
            % SEE ALSO: setDoppler
            
            id_snr = this.obs_code(:, 1) == 'S';
            if (nargin == 2)
                id_snr = id_snr & (this.system == sys_c)';
            end
            snr = zero2nan(this.obs(id_snr, :)');
        end
        
        function [obs, idx] = getObs(this, flag, system, prn)
            % get observation and index corresponfing to the flag
            % SYNTAX this.getObsIdx(flag, <system>)
            if nargin > 3
                idx = this.getObsIdx(flag, system, prn);
            elseif nargin > 2
                idx = this.getObsIdx(flag, system);
            else
                idx = this.getObsIdx(flag);
            end
            obs = this.obs(idx,:);
        end
        
        function [idx] = getObsIdx(this, flag, system, prn)
            % get observation index corresponfing to the flag
            % SYNTAX this.getObsIdx(flag, <system>)
            idx = sum(this.obs_code(:,1:length(flag)) == repmat(flag,size(this.obs_code,1),1),2) == length(flag);
            if nargin > 2
                idx = idx & (this.system == system)';
            end
            if nargin > 3
                idx = idx & reshape(this.prn == prn, length(this.prn),1);
            end
            idx = find(idx);
            idx(idx == 0) = [];
        end
        
        function [obs, idx, snr, cycle_slips] = getPrefObsCh(this, flag, system, max_obs_type)
            % get observation index corresponfing to the flag using best
            % channel according to the feinition in GPS_SS, GLONASS_SS
            % SYNTAX this.getObsIdx(flag, <system>)
            
            cycle_slips = [];
            if length(flag)==3
                idx = sum(this.obs_code == repmat(flag,size(this.obs_code,1),1),2) == 3;
                idx = idx & [this.system == system]';
                %this.legger.addWarning(['Unnecessary Call obs_type already determined, use getObsIdx instead'])
                [obs,idx] = this.getObs(flag, system);
            elseif length(flag) >= 2
                flags = zeros(size(this.obs_code,1),3);
                sys_idx = [this.system == system]';
                sys = this.cc.getSys(system);
                band = find(sys.CODE_RIN3_2BAND == flag(2));
                if isempty(band)
                    this.log.addWarning('Obs not found',200);
                    obs = [] ; idx = [];
                    return;
                end
                preferences = sys.CODE_RIN3_ATTRIB{band}; % get preferences
                sys_obs_code = this.obs_code(sys_idx,:); % get obs code for the given system
                sz =size(sys_obs_code,1);
                complete_flags = [];
                if nargin < 4
                    max_obs_type = length(preferences);
                end
                % find the betters flag present
                for j = 1 : max_obs_type
                    for i = 1:length(preferences)
                        if sum(sum(sys_obs_code == repmat([flag preferences(i)],sz,1),2)==3)>0
                            complete_flags = [complete_flags; flag preferences(i)];
                            preferences(i) = [];
                            break
                        end
                    end
                end
                if isempty(complete_flags)
                    this.log.addWarning('Obs not found',200);
                    obs = [] ; idx = [];
                    return;
                end
                max_obs_type = size(complete_flags,1);
                idxes = [];
                prn = [];
                for j = 1 : max_obs_type
                    flags = repmat(complete_flags(j,:),size(this.obs_code,1),1);
                    idxes = [idxes  sum(this.obs_code == flags,2) == 3];
                    prn = unique( [prn; this.prn(idxes(: , end )>0)]);
                end
                
                n_opt = size(idxes,2);
                n_epochs = size(this.obs,2);
                obs = zeros(length(prn)*n_opt,n_epochs);
                snr = zeros(size(obs));
                if flag(1) == 'L'
                    cycle_slips = sparse(size(obs,1),size(obs,2));
                end
                flags = zeros(length(prn)*n_opt,3);
                for s = 1:length(prn) % for each satellite and each epoch find the best (but put them on diffrent lines)
                    sat_idx = sys_idx & this.prn==prn(s);
                    
                    tmp_obs = zeros(n_opt,n_epochs);
                    take_idx = ones(1,n_epochs)>0;
                    for i = 1 : n_opt
                        c_idx = idxes(:, i) & sat_idx;
                        snr_idx = idxCharLines(this.obs_code, ['S' complete_flags(i,2:3)]) & sat_idx;
                        if sum(c_idx)>0
                            obs((s-1)*n_opt+i,take_idx) = this.obs(c_idx,take_idx);
                            
                            snr((s-1)*n_opt+i,take_idx) = this.obs(snr_idx,take_idx);
                            if ~isempty(this.outlier_idx_ph) && flag(1) == 'L'% take off outlier
                                ph_idx = this.ph_idx == find(c_idx);
                                obs((s-1)*n_opt+i,this.outlier_idx_ph(:,ph_idx)) = 0;
                                snr((s-1)*n_opt+i,this.outlier_idx_ph(:,ph_idx)) = 0;
                                cycle_slips((s-1)*n_opt+i,:) = this.cycle_slip_idx_ph(:,ph_idx)';
                            end
                            flags((s-1)*n_opt+i,:) = this.obs_code(c_idx,:);
                        end
                        take_idx = take_idx & obs((s-1)*n_opt+i,:) == 0;
                    end
                end
                prn = reshape(repmat(prn,1,n_opt)',length(prn)*n_opt,1);
                % remove all empty lines
                empty_idx = sum(obs==0,2) == n_epochs;
                obs(empty_idx,:) = [];
                prn(empty_idx,:) = [];
                flags(empty_idx,:) = [];
                flags=char(flags);
                idx = zeros(length(prn),1);
                for i = 1:length(prn)
                    idx(i) = find(sys_idx & this.prn == prn(i) & sum(this.obs_code == repmat(flags(i,:) ,size(this.obs_code,1) ,1),2) == 3);
                end
            else
                this.log.addError(['Invalid length of obs code(' num2str(length(flag)) ') can not determine preferred observation'])
            end
        end
        
        function [obs_set] = getIonoFree(this, flag1, flag2, system, max_obs_type)
            % get Iono free combination for the two selcted measurements
            % SYNTAX [obs] = this.getIonoFree(flag1, flag2, system)
            if not(flag1(1)=='C' | flag1(1)=='L' | flag2(1)=='C' | flag2(1)=='L')
                this.log.addWarning('Can not produce IONO free combination for the selcted observation')
                return
            end
            if flag1(1)~=flag2(1)
                this.log.addWarning('Incompatible observation type')
                return
            end
            if nargin <5
                max_obs_type = 1;
            end
            [obs1, idx1, snr1, cs1] = this.getPrefObsCh(flag1, system, max_obs_type);
            [obs2, idx2, snr2, cs2] = this.getPrefObsCh(flag2, system, max_obs_type);
            
            prn1 = this.prn(idx1);
            prn2 = this.prn(idx2);
            
            common_prn = intersect(prn1, prn2);
            sset_idx1 = ismember(prn1 , common_prn);
            sset_idx2 = ismember(prn2 , common_prn);
            prn1 = prn1(sset_idx1);
            prn2 = prn2(sset_idx2);
            idx1 = idx1(sset_idx1);
            idx2 = idx2(sset_idx2);
            obs1 = obs1(sset_idx1,:);
            obs2 = obs2(sset_idx2,:);
            snr1 = snr1(sset_idx1,:);
            snr2 = snr2(sset_idx2,:);
            if flag1(1)=='L'
                cs1 = cs1(sset_idx1,:);
                cs2 = cs2(sset_idx2,:);
            end
            
            %%% find the longer idx and replicate th other one to match the
            %%% prn
            if length(idx1) > length(idx2)
                idx_tmp = zeros(size(idx1));
                obs_tmp = zeros(size(obs1));
                snr_tmp = zeros(size(snr1));
                duplicate = prn1(1:end-1) == prn1(2:end);
                idx_tmp(~duplicate) = idx2;
                obs_tmp(~duplicate,:) = obs2;
                snr_tmp(~duplicate,:) = snr2;
                idx_tmp(duplicate) = idx_tmp(find(duplicate)+1);
                obs_tmp(duplicate,:) = obs_tmp(find(duplicate)+1,:);
                snr_tmp(duplicate,:) = snr_tmp(find(duplicate)+1,:);
                idx2 = idx_tmp;
                obs2 = obs_tmp;
                snr2 = snr_tmp;
                if flag1(1)=='L'
                    cs_tmp = zeros(size(cs1));
                    cs_tmp(~duplicate,:) = cs2;
                    cs_tmp(duplicate,:) = cs_tmp(find(duplicate)+1,:);
                    cs2 = cs_tmp;
                end
            else
                idx_tmp = zeros(size(idx2));
                obs_tmp = zeros(size(obs2));
                snr_tmp = zeros(size(snr2));
                duplicate = [prn2(1:end-1) == prn2(2:end); false];
                idx_tmp(~duplicate) = idx1;
                obs_tmp(~duplicate,:) = obs1;
                snr_tmp(~duplicate,:) = snr1;
                idx_tmp(duplicate) = idx_tmp(find(duplicate)+1);
                obs_tmp(duplicate,:) = obs_tmp(find(duplicate)+1,:);
                snr_tmp(duplicate,:) = snr_tmp(find(duplicate)+1,:);
                idx1 = idx_tmp;
                obs1 = obs_tmp;
                snr1 = snr_tmp;
                if flag1(1)=='L'
                    cs_tmp = zeros(size(cs2));
                    cs_tmp(~duplicate,:) = cs1;
                    cs_tmp(duplicate,:) = cs_tmp(find(duplicate)+1,:);
                    cs1 = cs_tmp;
                end
            end
            %             obs1 = this.obs(idx1,:);
            %             obs2 = this.obs(idx2,:);
            prn = this.prn(idx1);
            
            wl1 = this.wl(idx1);
            wl2 = this.wl(idx2);
            
            if isempty(obs1)|isempty(obs2)
                obs = [];
                prn = [];
                return
            end
            
            
            % put zeros to NaN
            obs1(obs1 == 0) = NaN;
            obs2(obs2 == 0) = NaN;
            snr1(snr1 == 0) = NaN;
            snr2(snr2 == 0) = NaN;
            
            %gte wavelenghts
            n_ep = size(obs1,2);
            inv_wl1 = repmat(1./this.wl(idx1),1,n_ep); %1/wl1;
            inv_wl2 = repmat(1./this.wl(idx2),1,n_ep); % 1/wl2;%
            alpha1 = ((inv_wl1).^2 )./ ( (inv_wl1).^2 - (inv_wl2).^2 );
            alpha2 = ((inv_wl2).^2)./ ( (inv_wl1).^2 - (inv_wl2).^2 );
            obs = alpha1 .* repmat(wl1,1,n_ep) .* obs1 - alpha2.* repmat(wl2,1,n_ep) .* obs2;
            snr = sqrt((alpha1.* snr1).^2 + (-alpha2 .* snr2).^2); % snr trated as std
            wl = alpha2(:,1)./alpha1(:,2) .* wl2;%alpha2 ./ alpha1 .* repmat(wl2,1,n_ep);
            % set NaN to 0
            nanidx = isnan(obs);
            obs(nanidx) = 0;
            snr(nanidx) = 0;
            snr(isnan(snr)) = 0;
            obs_code = [repmat(system,length(idx2),1) this.obs_code(idx1,:) this.obs_code(idx2,:) repmat('I',length(idx2),1)];
            go_id = this.cc.getIndex(system, prn);
            el = this.sat.el(:,go_id);
            el(obs' ==  0) = 0;
            az = this.sat.az(:,go_id);
            az(obs' ==  0) = 0;
            obs_set = Observation_Set(this.time.getCopy(), obs' ,obs_code, wl', el, az, prn');
            obs_set.snr = snr';
            if flag1(1)=='L'
               obs_set.cycle_slip = (cs1 | cs2)';
            end
        end
        
        function [obs_set]  = getPrefIonoFree(this, obs_type, system)
            % get Preferred Iono free combination for the two selcted measurements
            % SYNTAX [obs] = this.getIonoFree(flag1, flag2, system)
            iono_pref = this.cc.getSys(system).IONO_FREE_PREF;
            is_present = zeros(size(iono_pref,1),1) < 1;
            for i = size(iono_pref,1)
                % check if there are observation for the selected channel
                if sum(iono_pref(i,1) == this.obs_code(:,2) & iono_pref(i,1) == this.obs_code(:,1)) > 0 & sum(iono_pref(i,2) == this.obs_code(:,2) & iono_pref(i,1) == this.obs_code(:,1)) > 0
                    is_present(i) = true;
                end
            end
            iono_pref = iono_pref(is_present,:);
            [obs_set]  = this.getIonoFree([obs_type iono_pref(1,1)], [obs_type iono_pref(1,2)], system);
        end
    end
    
    % ==================================================================================================================================================
    % SETTER
    % ==================================================================================================================================================
    methods
        function setStatic(this)
            % Set the internal status of the Receiver as static
            % SYNTAX: 
            %   this.setStatic()
            this.static = true;
        end
        
        function setDynamic(this)
            % Set the internal status of the Receiver as dynamic
            % SYNTAX: 
            %   this.setDynamic()
            this.static = false;
        end
        
        function setDoppler(this, dop, wl, id_dop)
            % set the snr observations
            % SYNTAX: [pr, id_pr] = this.setDoppler(<sys_c>)
            % SEE ALSO:  getDoppler
            dop = bsxfun(@rdivide, zero2nan(dop'), wl);
            this.obs(id_dop, :) = nan2zero(dop');
        end
        
        function setSNR(this, snr, id_snr)
            % set the snr observations
            % SYNTAX: [pr, id_pr] = this.setSNR(<sys_c>)
            % SEE ALSO:  getSNR
            this.obs(id_snr, :) = nan2zero(snr');
        end
        
        function setPhases(this, ph, wl, id_ph)
            % set the phases observations in meter (not cycles)
            % SYNTAX: setPhases(this, ph, wl, id_ph)
            % SEE ALSO: getPhases getPseudoRanges setPseudoRanges
            
            ph = bsxfun(@rdivide, zero2nan(ph'), wl);
            this.obs(id_ph, :) = nan2zero(ph);
        end
        
        function setPseudoRanges(this, pr, id_pr)
            % set the pseudo ranges observations in meter (not cycles)
            % SYNTAX: [pr, id_pr] = this.getPseudoRanges(<sys_c>)
            % SEE ALSO: getPhases setPhases getPseudoRanges
            this.obs(id_pr, :) = nan2zero(pr');
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
        
        
        function applyGroupDelay(this)
            if this.group_delay_status == 0
                this.GroupDelay(1);
                this.group_delay_status = 1; %applied
            end
        end
        
        function removeGroupDelay(this)
            if this.group_delay_status == 1
                this.GroupDelay(-1);
                this.group_delay_status = 0; %applied
            end
        end
        
        function applyDtSat(this)
            if this.dts_delay_status == 0
                this.DtSat(1);
                this.dts_delay_status = 1; %applied
            end
        end
        
        function removeDtSat(this)
            if this.dts_delay_status == 1
                this.DtSat(-1);
                this.dts_delay_status = 0; %applied
            end
        end
        
        function smoothAndApplyDt(this)
            % Smooth dt * c correction computed from init-positioning with a spline with base 3 * rate,
            % apply the smoothed dt to pseudo-ranges and phases
            % SYNTAX:
            %   this.smoothAndApplyDt()
            id_ko = this.dt == 0;
            lim = getOutliers(this.dt(:,1) ~= 0 & abs(Core_Pre_Processing.diffAndPred(this.dt(:,1),2)) < 1e-8);            
            dt = simpleFill1D(zero2nan(this.dt(:,1)), this.dt == 0, 'spline');
            for i = 1 : size(lim, 1)
                if lim(i,2) - lim(i,1) > 5
                    dt(lim(i,1) : lim(i,2)) = splinerMat([], dt(lim(i,1) : lim(i,2)), 3);
                end
            end
            this.dt = simpleFill1D(zero2nan(dt), id_ko, 'spline');
            this.applyDtRec(dt)
            %this.dt_pr = this.dt_pr + this.dt;
            %this.dt_ph = this.dt_ph + this.dt;
            %this.dt = zeros(size(this.dt_pr));
        end
        
        function applyDtRec(this, dt_pr, dt_ph)
            % Apply dt * c correction to pseudo-ranges and phases
            % SYNTAX:
            %   this.applyDtRec(dt)
            %   this.applyDtRec(dt_pr, dt_ph)
            
            narginchk(2,3);
            if nargin == 2
                dt_ph = dt_pr;
            end
            [pr, id_pr] = this.getPseudoRanges();
            pr = bsxfun(@minus, pr, dt_pr * 299792458);
            this.setPseudoRanges(pr, id_pr);
            [ph, wl, id_ph] = this.getPhases();
            ph = bsxfun(@minus, ph, dt_ph * 299792458);
            this.setPhases(ph, wl, id_ph);
            this.time.addSeconds( - dt_pr);
        end
        
        
        
        function [obs, prn, sys, flag] = getBestCodeObs(this)
            % INPUT:
            % OUPUT:
            %    obs = observations [n_obs x n_epoch];
            %    prn = satellite prn [n_obs x 1];
            %    sys = system [n_obs x 1];
            %    flag = flag of the observation [ n_obs x 7] iono free
            %    combination are labeled with the obs code of both observations
            % DESCRIPTION: get "best" avaliable code or code combination for the given system
            n_epoch = size(this.obs,2);
            obs = [];
            sys = [];
            prn = [];
            flag = [];
            for i=1:this.cc.getNumSat()
                sat_idx = this.getObsIdx('C',this.cc.system(i),this.cc.prn(i));
                sat_idx = sat_idx(this.active_ids(sat_idx));
                if ~isempty(sat_idx)
                    % get epoch for which iono free is possible
                    sat_obs = this.obs(sat_idx,:);
                    av_idx = sum((sat_obs > 0),1)>0;
                    freq = str2num(this.obs_code(sat_idx,2)); %#ok<ST2NM>
                    u_freq = unique(freq);
                    if length(u_freq)>1
                        sat_freq = (sat_obs > 0) .*repmat(freq,1,n_epoch);
                        u_freq_sat = zeros(length(u_freq),n_epoch);
                        for e = 1 : length(u_freq)
                            u_freq_sat(e,:) = sum(sat_freq == u_freq(e))>0;
                        end
                        iono_free = sum(u_freq_sat)>1;
                    else
                        iono_free = false(1,n_epoch);
                    end
                    freq_list = this.cc.getSys(this.cc.system(i)).CODE_RIN3_2BAND;
                    track_list = this.cc.getSys(this.cc.system(i)).CODE_RIN3_ATTRIB;
                    if sum(iono_free > 0)
                        %this.sat.avail_index(:,i) = iono_free; % epoch for which observation is present
                        % find first freq obs
                        to_fill_epoch = iono_free;
                        first_freq = [];
                        f_obs_code = [];
                        ff_idx = zeros(n_epoch,1);
                        for f = 1 :length(freq_list)
                            track_prior = track_list{f};
                            for c = 1:length(track_prior)
                                if sum(to_fill_epoch) > 0
                                    [obs_tmp,idx_tmp] = this.getObs(['C' freq_list(f) track_prior(c)],this.cc.system(i),this.cc.prn(i));
                                    %obs_tmp(obs_tmp==0) = nan;
                                    if ~isempty(obs_tmp)
                                        first_freq = [first_freq; zeros(1,n_epoch)];
                                        first_freq(end, to_fill_epoch) = obs_tmp(to_fill_epoch);
                                        ff_idx(to_fill_epoch) = size(first_freq,1);
                                        f_obs_code = [f_obs_code; this.obs_code(idx_tmp,:)];
                                        to_fill_epoch = to_fill_epoch & (obs_tmp == 0);
                                    end
                                end
                            end
                        end
                        % find second freq obs
                        to_fill_epoch = iono_free;
                        second_freq = [];
                        s_obs_code = [];
                        sf_idx = zeros(n_epoch,1);
                        for f = 1 :length(freq_list)
                            track_prior = track_list{f};
                            for c = 1:length(track_prior)
                                if sum(to_fill_epoch) > 0
                                    [obs_tmp,idx_tmp] = this.getObs(['C' freq_list(f) track_prior(c)],this.cc.system(i),this.cc.prn(i));
                                    %obs_tmp = zero2nan(obs_tmp);
                                    
                                    % check if obs has been used already as first frequency
                                    if ~isempty(obs_tmp)
                                        
                                        ff_ot_i = find(sum(f_obs_code(:,1:2) == repmat(this.obs_code(idx_tmp,1:2),size(f_obs_code,1),1),2)==2);
                                        if ~isempty(ff_ot_i)
                                            obs_tmp(sum(first_freq(ff_ot_i,:),1)>0) = 0; % delete epoch where the obs has already been used for the first frequency
                                            
                                        end
                                        
                                        if sum(obs_tmp(to_fill_epoch)>0)>0 % if there is some new observation
                                            second_freq = [second_freq; zeros(1,n_epoch)];
                                            second_freq(end, to_fill_epoch) = obs_tmp(to_fill_epoch);
                                            sf_idx(to_fill_epoch) = size(second_freq,1);
                                            s_obs_code = [s_obs_code; this.obs_code(idx_tmp,:)];
                                            to_fill_epoch = to_fill_epoch & (obs_tmp == 0);
                                        end
                                    end
                                end
                            end
                        end
                        first_freq  = zero2nan(first_freq);
                        second_freq = zero2nan(second_freq);
                        % combine the two frequencies
                        for k = 1 : size(f_obs_code,1)
                            for y = 1 : size(s_obs_code,1)
                                inv_wl1 = 1/this.wl(this.getObsIdx(f_obs_code(k,:),this.cc.system(i),this.cc.prn(i)));
                                inv_wl2 = 1/this.wl(this.getObsIdx(s_obs_code(y,:),this.cc.system(i),this.cc.prn(i)));
                                %                                 if ((inv_wl1).^2 - (inv_wl2).^2) == 0
                                %                                     keyboard
                                %                                 end
                                
                                obs_tmp = ((inv_wl1).^2 .* first_freq(k,:) - (inv_wl2).^2 .* second_freq(y,:))./ ( (inv_wl1).^2 - (inv_wl2).^2 );
                                obs_tmp(isnan(obs_tmp)) = 0;
                                if sum(obs_tmp>0)>0
                                    obs = [obs; obs_tmp];
                                    prn = [prn; this.cc.prn(i)];
                                    sys = [sys; this.cc.system(i)];
                                    flag = [flag; [f_obs_code(k,:) s_obs_code(y,:) 'I' ]];
                                end
                            end
                        end
                        %                     end
                        %                     if sum(xor(av_idx, iono_free))>0
                    else % do not mix iono free and not combined observations
                        % find best code
                        to_fill_epoch = av_idx;
                        %this.sat.avail_index(:,i) = av_idx; % epoch for which observation is present
                        for f = 1 :length(freq_list)
                            track_prior = track_list{f};
                            for c = 1:length(track_prior)
                                if sum(to_fill_epoch) > 0
                                    [obs_tmp,idx_tmp] = this.getObs(['C' num2str(freq_list(f)) track_prior(c)],this.cc.system(i),this.cc.prn(i));
                                    if ~isempty(obs_tmp)
                                        obs = [obs; zeros(1,n_epoch)];
                                        obs(end,to_fill_epoch) = obs_tmp(to_fill_epoch);
                                        prn = [prn; this.cc.prn(i)];
                                        sys = [sys; this.cc.system(i)];
                                        flag = [flag; sprintf('%-7s',this.obs_code(idx_tmp,:))];
                                        to_fill_epoch = to_fill_epoch & (obs_tmp < 0);
                                    end
                                end
                            end
                        end
                        
                    end
                    
                end
            end
            %%% Remove obs for which coordinates of satellite are non
            %%% available
            for o = 1:length(prn)
                s = this.cc.getIndex(sys(o),prn(o));
                o_idx_l = obs(o,:)>0;
                times = this.time.getSubSet(o_idx_l);
                times.addSeconds(-obs(o,o_idx_l)'/Go_State.V_LIGHT); % add roucg time of flight
                xs = this.sat.cs.coordInterpolate(times,s);
                to_remove = isnan(xs(:,1));
                o_idx = find(o_idx_l);
                to_remove = o_idx(to_remove);
                obs(o,to_remove) = 0;
                if sum(obs(o,:)) == 0 % if line has no more observation
                    obs(o,:) = [];
                    prn(o) = [];
                    sys(o) = [];
                    flag(o,:) = [];
                end
            end
        end
        
        function res = staticPPP(this)
            ls = Least_Squares_Manipulator();
            ppp_opt.tropo = true; %this.state.flag_tropo;
            ppp_opt.tropo_g = true;%this.state.flag_tropo_gradient;
            ls.setUpPPP(this, ppp_opt);
            ls.Astack2Nstack();
            % set time regularization for the troposphere and its gradients  
%             ls.setTimeRegularization(1, 0.00001);
%             ls.setTimeRegularization(2, 0.00001);
%             ls.setTimeRegularization(3, 0.00001);
            ls.setTimeRegularization(6, 1e-7*this.rate*Go_State.V_LIGHT/0.005);%
            ls.setTimeRegularization(7, this.state.std_tropo / 3600 * this.rate/0.003 * 10);
            ls.setTimeRegularization(8, this.state.std_tropo_gradient / 3600 * this.rate/0.005);
            ls.setTimeRegularization(9, this.state.std_tropo_gradient / 3600 * this.rate/0.005);
            [x, res] = ls.solve();
            
            coo    = x(1:3,1);
            clock = x(x(:,2) == 6,1);
            tropo = x(x(:,2) == 7,1);
            amb = x(x(:,2) == 5,1);
            gntropo = x(x(:,2) == 8,1);
            getropo = x(x(:,2) == 9,1);
            this.log.addMessage(sprintf('DEBUG: s02 = %f',mean(abs(res(res~=0)))));
            new_pos = this.xyz + coo';
            diff_from_rin = (new_pos  -this.xyz_approx)';
            this.log.addMessage(sprintf('DEBUG: distance from rine pos = %.3f %.3f %.3f',diff_from_rin));
            this.log.addMessage(sprintf('DEBUG: distance from rine pos enu = %.3f %.3f %.3f',global2localVel(diff_from_rin,this.xyz')));
        end
        
        function initPositioning(this, sys_c)
            % run the most appropriate init prositioning step depending on the static flag
            % calls initStaticPositioning() or initDynamicPositioning()
            % SYNTAX:
            %   this.initPositioning();
            % INPUT:
            %   sys_c = wanted system
            % Init "errors"
            this.log.addMarkedMessage('Computing position and clock errors using a code only solution')
            this.sat.err_tropo = zeros(this.time.length, this.cc.getNumSat());
            this.sat.err_iono  = zeros(this.time.length, this.cc.getNumSat());
            this.sat.solid_earth_corr  = zeros(this.time.length, this.cc.getNumSat());
            this.log.addMessage(this.log.indent('Applying satellites Differencial Code Biases (DCB)', 6))
            % if not applied apply gruop delay
            this.applyGroupDelay();
            this.log.addMessage(this.log.indent('Applying satellites clock errors and eccentricity dependent realtivistic correction', 6))
            this.applyDtSat();
            
            % get best observation for all satellites and all epochs
            this.log.addMessage(this.log.indent('Get best code combination available for each satellites and epoch', 6))
            [obs, prn, sys, flag] = this.getBestCodeObs;
            % remove unwanted system
            if nargin < 2
                sys_c = this.cc.sys_c;
            end
            sys_idx = false(size(sys));
            for s = 1:length(sys_c)
                sys_idx = sys_idx | sys == sys_c(s);
            end
            obs(~sys_idx,:) = [];
            prn(~sys_idx,:) = [];
            sys(~sys_idx,:) = [];
            flag(~sys_idx,:) = [];
            
            %this.static = 0;
            if this.isStatic()
                this.initStaticPositioning(obs, prn, sys, flag)
            else
                this.initDynamicPositioning(obs, prn, sys, flag)
            end
            
            % Apply dt from the clock estimated by initPositioning
            this.log.addMessage(this.log.indent('Smooth and apply the clock error of the receiver', 6))            
        end
        
        function initStaticPositioning(this, obs, prn, sys, flag)
            % SYNTAX:
            %   this.StaticPositioning(obs, prn, sys, flag)
            %
            % INPUT:
            % obs : observations [meters]
            % prn : prn of observations
            % sys : sys of observations
            % flag : name of obsevation [obs_code1 obs_code2 comb_code]
            %        comb_code --> Iono Free = I
            % OUTPUT:
            %
            % DESCRIPTION:
            %   Get postioning using code observables
            
            if nargin == 1
                % get best observation for all satellites and all epochs
                this.log.addMessage(this.log.indent('Get best code combination available for each satellites and epoch', 6))
                [obs, prn, sys, flag] = this.getBestCodeObs;
                % remove unwanted system
                if nargin < 2
                    sys_c = this.cc.sys_c;
                end
                sys_idx = false(size(sys));
                for s = 1:length(sys_c)
                    sys_idx = sys_idx | sys == sys_c(s);
                end
                obs(~sys_idx,:) = [];
                prn(~sys_idx,:) = [];
                sys(~sys_idx,:) = [];
                flag(~sys_idx,:) = [];
            end
            
            this.log.addMessage(this.log.indent('Starting initial static postioning',6))
            iono_free = flag(1,7) == 'I';
            % It should be this:
            % approx_pos_unknown = all(this.xyz_approx(:) == 0);
            approx_pos_unknown = true;
            opt.rid_ep = false; % do not estimate channel dipendent error at each epoch
            
            sub_sample = false;
            if  approx_pos_unknown
                this.xyz = this.xyz_approx;
                this.log.addMessage(this.log.indent('Getting coarse position on subsample of data',6))
                if sum(sum(obs,1) > 0) >= 100
                    % sub sample observations
                    sub_sample = true;
                    idx_ss = 1 : 100; % min(100, size(obs,2)) ; %(1: round(size(obs,2) / 100):size(obs,2));
                    
                    obs_ss = zeros(size(obs));
                    obs_ss(:, idx_ss) = obs(:, idx_ss);
                    prn_ss =  prn;
                    sys_ss = sys;
                    flag_ss = flag;
                    
                    % remove line that might be empty
                    empty_sat = sum(obs_ss,2) == 0;
                    obs_ss(empty_sat, :) = [];
                    prn_ss(empty_sat, :)  = [];
                    flag_ss(empty_sat, :) = [];
                    sys_ss(empty_sat, :)  = [];
                end
                cut_off = 15;
                % first estimation noatmosphere
                opt.coord_corr = 1;
                opt.max_it = 10;
                if sub_sample
                    this.codeStaticPositionig(obs_ss, prn_ss, sys_ss, flag_ss, opt);
                    [obs_ss, sys_ss, prn_ss, flag_ss] = this.removeUndCutOff(obs_ss, sys_ss, prn_ss, flag_ss, cut_off);
                else
                    this.codeStaticPositionig(obs, prn, sys, flag, opt);
                end
                
                % update atmosphere
                this.updateAzimuthElevation();
                this.updateErrTropo();
                if ~iono_free
                    this.updateErrIono();
                end
                
                % second estimation with atmosphere
                this.log.addMessage(this.log.indent('improving estimation',6))
                opt.coord_corr = 1;
                opt.max_it = 10;
                if sub_sample
                    this.codeStaticPositionig(obs_ss, prn_ss, sys_ss, flag_ss, opt);
                else
                    this.codeStaticPositionig(obs, prn, sys, flag, opt);
                end
            else
                this.xyz = this.xyz_approx;
            end
            
            if sub_sample
                % update avalibilty index
                this.updateAllAvailIndex();
            end
            
            % update Atmosphere Corrections
            this.updateAzimuthElevation();
            this.updateErrTropo('all', 1);
            if ~iono_free
                this.updateErrIono();
            end
            
            % update solid earth corrections
            this.updateSolidEarthCorr();
            % final estimation
            opt.max_it = 1;
            opt.coord_corr = 0.1;
            opt.no_pos = true;
            this.log.addMessage(this.log.indent('Get final clock error estimation to sysncronize satellite positions',6))
            this.codeStaticPositionig(obs, prn, sys, flag, opt); % get a first estimation of receiver clock offset to get correct orbit
            opt.no_pos = false;
            %%% remove obs under cu off
            [obs, sys, prn, flag] = this.removeUndCutOff(obs, sys, prn, flag, cut_off);
            this.log.addMessage(this.log.indent('Final estimation',6))
            this.codeStaticPositionig(obs, prn, sys, flag, opt);
        end
        
        function initDynamicPositioning(this, obs, prn, sys, flag)
            % SYNTAX:
            %   this.initDynamicPositioning(obs, prn, sys, flag)
            %
            % INPUT:
            % obs : observations [meters]
            % prn : prn of observations
            % sys : sys of observations
            % flag : name of obsevation [obs_code1 obs_code2 comb_code]
            %        comb_code --> Iono Free = I
            % OUTPUT:
            %
            % DESCRIPTION: get dynamic postion using code observables
            % (independent epochs, no kalman filters or regularization)
            
            if nargin == 1
                % get best observation for all satellites and all epochs
                this.log.addMessage(this.log.indent('Get best code combination available for each satellites and epoch', 6))
                [obs, prn, sys, flag] = this.getBestCodeObs;
                % remove unwanted system
                if nargin < 2
                    sys_c = this.cc.sys_c;
                end
                sys_idx = false(size(sys));
                for s = 1:length(sys_c)
                    sys_idx = sys_idx | sys == sys_c(s);
                end
                obs(~sys_idx,:) = [];
                prn(~sys_idx,:) = [];
                sys(~sys_idx,:) = [];
                flag(~sys_idx,:) = [];
            end
            
            %initialize modeled error matrix
            this.log.addMessage(this.log.indent('Starting dynamic positioning', 6))
            if isempty(this.sat.avail_index)
                this.sat.avail_index = zeros(this.time.length, this.cc.getNumSat());
            end
            this.sat.err_tropo = zeros(this.time.length, this.cc.getNumSat());
            this.sat.err_iono  = zeros(this.time.length, this.cc.getNumSat());
            this.sat.solid_earth_corr  = zeros(this.time.length, this.cc.getNumSat());
            
            %check if the combination is ionofree
            iono_free = flag(1,7) == 'I';
            n_epochs         = this.time.length;
            code_bias_flag = [sys flag];
            % compute a column with an integer that indicate whice
            % code_bias to estimate for each obs
            code_bias_ord = zeros(size(code_bias_flag,1),1);
            u_code_bias_flag = unique(cellstr(code_bias_flag),'stable');
            n_clocks = length(u_code_bias_flag);
            for i = 1 : n_clocks
                ch_idx_sat = sum([sys flag] == repmat( sprintf('%-8s',u_code_bias_flag{i}),length(sys),1),2) == 8;
                code_bias_ord(ch_idx_sat) = i;
            end
            sat = zeros(size(prn));
            for i = 1 : length(sat)
                sat(i) = this.cc.getIndex(sys(i),prn(i));
            end
            
            this.dt = zeros(n_epochs,n_clocks);
            this.flag_rid = u_code_bias_flag;
            if ~isempty(this.xyz_approx) && sum(this.xyz_approx) ~= 0
                this.xyz = repmat(this.xyz_approx, n_epochs, 1);
            else
                this.xyz = zeros(n_epochs,3);
            end
            
            XS = zeros(this.cc.getNumSat(), 3, n_epochs);
            % get Sat orbit correctiong ony for pseudorange
            for i = 1 : this.cc.getNumSat()
                c_sys = this.cc.system(i);
                c_prn = this.cc.prn(i);
                idx_sat = sys == c_sys & prn == c_prn;
                if sum(idx_sat) > 0 % if we have an obesrvation for the satellite
                    c_obs = obs(idx_sat,:);
                    
                    c_l_obs = colFirstNonZero(c_obs); %all best obs one one line
                    idx_obs = c_l_obs > 0; %epoch with obseravtion from the satellite
                    %update time of flight times
                    this.updateAvailIndex(c_l_obs,i);
                    this.updateTOT(c_l_obs,i); % update time of travel
                    
                    [XS_temp , ~] = this.getXSTxRot(i);
                    XS(i,:,idx_obs) = XS_temp';
                end
            end
            this.log.addMessage(this.log.indent('Getting very coarse postion from first valid epoch', 6))
            % get coarse postion from first valid epoch
            %find first valid epoch
            not_found = true;
            e = 0;
            while not_found && e < n_epochs
                e = e + 1;
                idx_obs = obs(:,e) > 0;
                clock_temp = unique(code_bias_ord(idx_obs));
                if sum(idx_obs) >= (3 + length(clock_temp));
                    not_found = false;
                end
            end
            
            x = [999 999 999];
            n_obs_ep = sum(idx_obs);
            v_clocks = 1 : n_clocks;
            while max(abs(x(1:3))) > 10
                clock_ep = code_bias_ord(idx_obs);
                pres_clock = (sum(repmat(clock_ep, 1, n_clocks) == repmat(v_clocks, n_obs_ep, 1), 1) > 0) .* v_clocks;
                pres_clock(pres_clock == 0) = [];
                XS_temp = XS(sat(idx_obs),:,e);
                XS_temp = XS_temp - repmat(this.xyz(e,:),sum(idx_obs),1);
                dist = sqrt(sum(XS_temp.^2,2));
                XS_temp_norm = XS_temp./repmat(dist,1,3);
                A_dcb = zeros(n_obs_ep, length(pres_clock));
                for i = 1:length(pres_clock)
                    A_dcb(clock_ep == pres_clock(i),i) = 1;
                end
                A_temp = [- XS_temp_norm A_dcb];
                y = obs(idx_obs, e) - dist;
                x = A_temp \ y;
                this.xyz(e,:) = this.xyz(e,:) +x(1:3)';
                this.dt(e,pres_clock) = x(4:end)';
                cur_xyz_est = this.xyz(e,:);
            end
            % if combination is not ionofree compute iono using coarse
            % postion estimation computed with the first valid epoch
            if ~iono_free
                for i = 1 : this.cc.getNumSat()
                    this.updateErrIono(i);
                end
            end
            residuals = zeros(size(obs));
            this.log.addMessage(this.log.indent('Getting coarse position for all epoch',6))
            for e = 1 : n_epochs
                idx_obs = obs(:,e) > 0;
                n_obs_ep = sum(idx_obs);
                v_clocks= 1:n_clocks;
                clock_ep = code_bias_ord(idx_obs);
                pres_clock = (sum(repmat(clock_ep,1,n_clocks) == repmat(v_clocks,n_obs_ep,1),1) > 0).*v_clocks;
                pres_clock(pres_clock == 0) = [];
                if   sum(idx_obs) >= (3 + length(pres_clock)) % if system is solvable
                    this.xyz(e,:) = cur_xyz_est;
                    XS_temp = XS(sat(idx_obs),:,e);
                    XS_temp = XS_temp - repmat(cur_xyz_est,sum(idx_obs),1);
                    dist = sqrt(sum(XS_temp.^2,2));
                    XS_temp_norm = XS_temp./repmat(dist,1,3);
                    A_dcb = zeros(n_obs_ep, length(pres_clock));
                    for i = 1:length(pres_clock)
                        A_dcb(clock_ep == pres_clock(i),i) = 1;
                    end
                    A_temp = [- XS_temp_norm A_dcb];
                    
                    y = obs(idx_obs, e) - dist;
                    x = A_temp \ y;
                    res = y - A_temp * x;
                    residuals(idx_obs,e) = res;
                    if not(isnan(this.xyz(e,:))) %& sum(abs(x(1:3))) < 300
                        this.xyz(e,:) = this.xyz(e,:) +x(1:3)';
                        cur_xyz_est = this.xyz(e,:);
                        this.dt(e,pres_clock) = x(4:end)'/Go_State.V_LIGHT;
                    end
                else
                    % keyboard
                end
            end
            
            %%% COMPUTE AGAIN XS BASED ON CURRENT CLOCK ESTIMATE
            %remove cut off
            this.updateAllAvailIndex();
            this.updateAzimuthElevation();
            this.updateErrTropo();
            if ~iono_free
                this.updateErrIono();
            end
            this.updateSolidEarthCorr();
            cut_off = 15;
            [obs, sys, prn, flag] = this.removeUndCutOff(obs, sys, prn, flag, cut_off);
            sat = zeros(size(prn));
            for i = 1 : length(sat)
                sat(i) = this.cc.getIndex(sys(i),prn(i));
            end
            XS = zeros(this.cc.getNumSat(), 3, n_epochs);
            dist = zeros(n_epochs, this.cc.getNumSat());
            % get Sat orbit correctiong ony for pseudorange
            for i = 1 : this.cc.getNumSat();
                c_sys = this.cc.system(i);
                c_prn = this.cc.prn(i);
                idx_sat = sys == c_sys & prn == c_prn;
                if sum(idx_sat) > 0 % if we have an observation for the satellite
                    c_obs = obs(idx_sat,:);
                    
                    c_l_obs = colFirstNonZero(c_obs); % all best obs one one line
                    idx_obs = c_l_obs > 0; % epoch with obseravtion from the satellite
                    % update time of flight times
                    this.updateAvailIndex(c_l_obs,i);
                    this.updateTOT(c_l_obs,i); % update time of travel
                    
                    [dist(:,i), XS_temp] = this.computeSyntObs('I',i); %%% consider multiple combinations (different iono corrections) on the same satellite, not handled yet
                    
                    XS(i,:,idx_obs) = rowNormalize(XS_temp)';
                end
            end
            this.log.addMessage(this.log.indent('Getting final estimation for all epoch',6))
            for e = 1 : n_epochs
                idx_obs = obs(:,e) > 0;
                n_obs_ep = sum(idx_obs);
                v_clocks= 1:n_clocks;
                clock_ep = code_bias_ord(idx_obs);
                pres_clock = (sum(repmat(clock_ep,1,n_clocks) == repmat(v_clocks,n_obs_ep,1),1) > 0).*v_clocks;
                pres_clock(pres_clock == 0) = [];
                if   sum(idx_obs) >= (3 + length(pres_clock)); % if system is solvable
                    XS_temp = XS(sat(idx_obs),:,e);
                    A_dcb = zeros(n_obs_ep, length(pres_clock));
                    for i = 1:length(pres_clock)
                        A_dcb(clock_ep == pres_clock(i),i) = 1;
                    end
                    A_temp = [- XS_temp A_dcb];
                    
                    y = obs(idx_obs, e) - dist(e, sat(idx_obs))';
                    % remove parmater with one clock only, they add no
                    % value create problem for the robust sdjustment
                    %                     one_idx = sum(A_temp(:,4:end),1) == 1;
                    %                     rm_col = 3+find(one_idx);
                    %                     rm_idx = A_temp(:,rm_col) > 0;
                    %                     y(rm_idx) = [];
                    %                     n_obs_ep = length(y);
                    %                     A_temp(rm_idx,:) = [];
                    %                     A_temp(:,rm_col) = [];
                    %                     pres_clock(rm_col-3) = [];
                    
                    b = zeros(size(y));
                    Q = speye(length(y));
                    [x, res, s02] = Least_Squares.solver(y, b, A_temp, Q);
                    %                     x = A_temp \ y;
                    %                     res = y - A_temp * x;
                    %residuals(idx_obs,e) = res;
                    %--- robust estimation ---
                    %                     for i = 1:1
                    %                         res = res/s02;
                    %                         Q = spdiags(min(res.^2,1000) ,0 ,n_obs_ep ,n_obs_ep);
                    %                         [x, res] = Least_Squares.solver(y, b, A_temp, Q);
                    %                     end
                    if max(abs(x(1:3))) <100;
                        this.xyz(e,:) = this.xyz(e,:) +x(1:3)';
                        this.dt(e,pres_clock) = x(4:end)'/Go_State.V_LIGHT;
                    else
                        this.dt(e,pres_clock) = 0;
                        this.xyz(e,:) = 0;
                    end
                end
            end
        end
        
        function [range, XS_loc] = computeSyntObs(this, obs_type, sat)
            % DESCRIPTION: get the estimate of one measurmenet based on the
            % current postion
            % INPUT:
            %   obs_type; type of obs I(ionofree) 1(first system freqeuncy) 2(second sytem frequency) 3 (third system frequency)
            % EXAMPLE:
            %   this.computeSyntObs(1,go_id)
            n_epochs = size(this.obs, 2);
            n_sat = this.cc.getNumSat();
            if isnumeric(obs_type)
                obs_type = num2str(obs_type);
            end
            if nargin < 3
                range = zeros(n_sat, n_epochs);
                for sat = 1 : n_sat
                    range(sat, :) = this.computeSyntObs(obs_type, sat);
                end
                XS_loc = [];
            else
                XS_loc = this.getXSLoc(sat);
                range = sqrt(sum(XS_loc.^2,2));
                sys = this.cc.system(sat);
                switch obs_type
                    case 'I'
                        iono_factor = 0;
                    case '1'
                        iono_factor = 1;
                    otherwise
                        wl_ref = this.cc.getGPS.F_VEC(1);
                        wl = this.cc.getSys(sys(j)).F_VEC(str2num(obs_type));
                        iono_factor= wl_ref^2/ wl^2;
                        
                end
                range = range + this.sat.err_tropo(:,sat) + iono_factor * this.sat.err_iono(:,sat);
                
                XS_loc(isnan(range),:) = [];
                %range = range';
                range = nan2zero(range)';
            end
            
        end
        
        function  updateSyntPhases(this, sys_c)
            %DESCRIPTION: update the content of the syntetic phase
            if nargin < 2
                sys_c = this.cc.sys_c;
            end
            synt_ph_obs = zero2nan(this.computeSyntCurObs( true, sys_c)');
            synt_ph_obs(this.outlier_idx_ph) = nan;
            this.synt_ph = synt_ph_obs;
        end
        
        function synt_ph = getSyntPhases(this)
            % DESCRIPTION: get current value of syntetic phase, in case not
            % present update it
            if isempty(this.synt_ph)
                this.updateSyntPhases();
            end
            synt_ph = this.synt_ph;
            synt_ph(this.outlier_idx_ph) = nan;
        end
        
        function synt_pr_obs = getSyntPrObs(this, sys_c)
            if nargin < 2
                sys_c = this.cc.sys_c;
            end
            synt_pr_obs = zero2nan(this.computeSyntCurObs(false, sys_c)');
        end
        
        function [synt_obs, xs_loc] = getSyntTwin(this, obs_set)
            % DESCRIPTION: get the syntethic twin for the observations
            % contained in obs_set
            % WARNING: the time of the observation set and the time of
            % receiver as to be the same
            synt_obs = zeros(size(obs_set.obs));
            xs_loc   = zeros(size(obs_set.obs,1),size(obs_set.obs,2),3);
            idx_ep_obs = obs_set.getTimeIdx(this.time.getSubSet(1),this.rate);
            for i = 1 : size(synt_obs,2)
                go_id = this.cc.getIndex(obs_set.obs_code(i,1),obs_set.prn(i));
                if length(obs_set.obs_code(i,:)) > 7 && obs_set.obs_code(i,8) == 'I'
                    [range, xs_loc_t] = this.computeSyntObs('I', go_id);
                else
                    f = find(this.cc.getSys(obs_set.obs_code(i,1)).CODE_RIN3_2BAND == obs_set.obs_code(i,3));
                    [range, xs_loc_t] = this.computeSyntObs(f, go_id);
                end
                
                idx_obs = obs_set.obs(:,i) ~= 0;
                idx_obs_r = idx_ep_obs(idx_obs); % <- to wihic epoch in the receiver the observation of the satellites in obesrvation set corresponds?
                idx_obs_r_l = false(size(range)); % get the logical equivalent
                idx_obs_r_l(idx_obs_r) = true;
                range_idx = range ~= 0;
                xs_idx = idx_obs_r_l(range_idx);
                synt_obs(idx_obs, i) = range(idx_obs_r);
                xs_loc(idx_obs, i,:) = permute(xs_loc_t(xs_idx,:),[1 3 2]);
            end
            
        end
                
        function timeShiftObs(this, tt)
            % DESCRIPTION: translate observations at different epoch based on linear modeling of the satellite
            % copute the sat postion at the current epoch
            XS = this.getXSLoc();
            
            % translate the time
            this.time.addSeconds(tt);
            % compute the sat postion at epoch traslated
            XS_t = this.getXSLoc();
            
            % compute the correction
            range = sqrt(sum(XS.^2,3));
            range_t = sqrt(sum(XS_t.^2,3));
            d_range = range_t - range;
            % Correct phases
            obs_idx = this.obs_code(:,1) == 'L';
            this.obs(obs_idx,:) = nan2zero(zero2nan(this.obs(obs_idx,:)) + bsxfun(@rdivide, d_range(:,this.go_id(obs_idx))', this.wl(obs_idx)));
            % Correct pseudo-ranges
            obs_idx = this.obs_code(:,1) == 'C';
            this.obs(obs_idx,:) = nan2zero(zero2nan(this.obs(obs_idx,:)) + d_range(:,this.go_id(obs_idx))');
        end
        
        function shiftToNominal(this)
            % DESCRIPTION: translate receiver observations to nominal epochs
            nominal_time = getNominalTime(this);
            tt = nominal_time - this.time;
            this.timeShiftObs(tt)            
        end
        
    end
        
    % ==================================================================================================================================================
    %  FUNCTIONS TO GET SATELLITE POSITION AT RIGHT TIME
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
            idx = this.sat.avail_index(:, sat) > 0;
            time_tx = this.time.getSubSet(idx);
            time_tx.addSeconds( - this.sat.tot(idx, sat));
        end
        
        function updateTOT(this, obs, sat)
            % SYNTAX:
            %   this.updateTOT(time_rx, dt);
            %
            % INPUT:
            %
            % OUTPUT:
            % DESCRIPTION:
            %   Update the signal time of travel based on range observations.
            % NOTE: to have satellite orbits at 1mm sysncrnonization time of
            % travel has to be know with a precision of ~100m
            %    0.001m * (3 km/s  +         2 km/s)   /       300000 km/s   ~ 100m
            %              ^                 ^                    ^
            %              |                 |                    |
            %         (sat speed) (earth rot at sat orbit) (speed of light)
            if isempty(this.sat.tot)
                this.sat.tot = zeros(size(this.sat.avail_index));
            end
            this.sat.tot(:, sat) =  ( obs' )/ goGNSS.V_LIGHT + this.dt(:, 1);
        end
        
        function updateAvailIndex(this, obs, sat)
            % DESCRIPTION: update avaliabilty of measurement on staellite
            if isempty(this.sat.avail_index)
                this.sat.avail_index = false(this.time.length, this.cc.getNumSat());
            end
            this.sat.avail_index(:,sat) = obs > 0;
        end
        
        function updateAllAvailIndex(this)
            % DESCRIPTION: update avaliabilty of measurement on all
            % satellite based on all code and phase
            if isempty(this.sat.avail_index)
                this.sat.avail_index = false(this.time.length, this.cc.getNumSat());
            end
            for s = 1 : this.cc.getNumSat()
                obs_idx = this.go_id == s & (this.obs_code(:,1) == 'C' | this.obs_code(:,1) == 'L');
                if sum(obs_idx) > 0
                    av_idx = colFirstNonZero(this.obs(obs_idx,:)) ~= 0 ;
                    this.sat.avail_index(:,s) = av_idx;
                end
            end
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
        
        function dtS = getDtS(this, sat)
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
            if nargin < 2
                dtS = zeros(size(this.sat.avail_index));
                for s = 1 : size(dtS,2)
                    dtS(this.sat.avail_index(:,s),s) = this.sat.cs.clockInterpolate(this.time.getSubSet(this.sat.avail_index(:,s)),s);
                end
            else
                idx = this.sat.avail_index(:,sat) > 0;
                dtS = this.sat.cs.clockInterpolate(this.time.getSubSet(idx), sat);
            end
        end
        
        function dtRel = getRelClkCorr(this, sat)
            % DESCRIPTION : get clock offset of the satellite due to
            % special relativity (eccntrcity term)
            idx = this.sat.avail_index(:,sat) > 0;
            [X,V] = this.sat.cs.coordInterpolate(this.time.getSubSet(idx),sat);
            dtRel = -2 * sum(conj(X) .* V, 2) / (goGNSS.V_LIGHT ^ 2); % Relativity correction (eccentricity velocity term)
        end
        
        
        function [XS_tx_r ,XS_tx] = getXSTxRot(this, sat)
            % SYNTAX:
            %   [XS_tx_r ,XS_tx] = this.getXSTxRot( sat)
            %
            % INPUT:
            % sat = go-id of the satellite
            % OUTPUT:
            % XS_tx = satellite position computed at trasmission time
            % XS_tx_r = Satellite postions at transimission time rotated by earth rotation occured
            % during time of travel
            % DESCRIPTION:
            %   Compute satellite positions at transmission time and rotate them by the earth rotation
            %   occured during time of travel of the signal
            if nargin > 1
                [XS_tx] = this.getXSTx(sat);
                [XS_tx_r]  = this.earthRotationCorrection(XS_tx, sat);
            else
                n_sat = this.cc.getNumSat();
                XS_tx_r = zeros(this.time.length, n_sat, 3);
                for i = 1 : n_sat
                    [XS_tx] = this.getXSTx(i);
                    [XS_tx_r_temp]  = this.earthRotationCorrection(XS_tx, i);
                    XS_tx_r(this.sat.avail_index(:,i) ,i ,:) = permute(XS_tx_r_temp, [1 3 2]);
                end
            end
        end
        
        function [XS_loc] = getXSLoc(this, sat)
            % SYNTAX:
            %   [XS_tx_r ,XS_tx] = this.getXSLoc( sat)
            %
            % INPUT:
            % sat = go-id of the satellite
            % OUTPUT:
            % XS_tx = satellite position computed at trasmission time
            % XS_tx_r = Satellite postions at transimission time rotated by earth rotation occured
            % during time of travel
            % DESCRIPTION:
            %   Compute satellite positions at transmission time and rotate them by the earth rotation
            %   occured during time of travel of the signal and subtract
            %   the postion term
            n_epochs = this.time.length;
            if nargin > 1
                sat_idx = this.sat.avail_index(:, sat) > 0;
                XS = this.getXSTxRot(sat);
                XS_loc = nan(n_epochs, 3);
                XS_loc(sat_idx,:) = XS;
                if size(this.xyz,1) == 1
                    XR = repmat(this.xyz, n_epochs, 1);
                else
                    XR = this.xyz;
                end
                XS_loc = XS_loc - XR;
            else
                n_sat = this.cc.getNumSat();
                XS_loc = zeros(n_epochs,n_sat,3);
                for i = 1 : n_sat
                    XS_loc(:,i ,:) = this.getXSLoc(i);
                end
            end
        end
        function [XS_tx] = getXSTx(this, sat)
            % SYNTAX:
            %   [XS_tx_frame , XS_rx_frame] = this.getXSTx()
            %
            % INPUT:
            %  obs : [1x n_epochs] pseudi range observations
            %  sta : index of the satellite
            % OUTPUT:
            % XS_tx = satellite position computed at trasmission time
            % DESCRIPTION:
            % Compute satellite positions at trasmission time
            time_tx = this.getTimeTx(sat);
            %time_tx.addSeconds(); % rel clok neglegible
            [XS_tx, ~] = this.sat.cs.coordInterpolate(time_tx,sat);
            
            
            %                 [XS_tx(idx,:,:), ~] = this.sat.cs.coordInterpolate(time_tx);
            %             XS_tx  = zeros(size(this.sat.avail_index));
            %             for s = 1 : size(XS_tx)
            %                 idx = this.sat.avail_index(:,s);
            %                 %%% compute staeliite position a t trasmission time
            %                 time_tx = this.time.subset(idx);
            %                 time_tx = time_tx.time_diff - this.sat.tot(idx,s)
            %                 [XS_tx(idx,:,:), ~] = this.sat.cs.coordInterpolate(time_tx);
            %             end
        end
        function [XS_r] = earthRotationCorrection(this, XS, sat)
            % SYNTAX:
            %   [XS_r] = this.earthRotationCorrection(XS)
            %
            % INPUT:
            %   XS      = positions of satellites
            %   time_rx = receiver time
            %   cc      = Constellation Collector
            %   sat     = satellite
            % OUTPUT:
            %   XS_r    = Satellite postions rotated by earth roattion occured
            %   during time of travel
            % DESCRIPTION:
            %   Rotate the satellites position by the earth rotation
            %   occured during time of travel of the signal
            
            %%% TBD -> consider the case XS and travel_time does not match
            XS_r = zeros(size(XS));
            
            idx = this.sat.avail_index(:,sat) > 0;
            travel_time = this.sat.tot(idx,sat);
            sys = this.cc.system(sat);
            switch char(sys)
                case 'G'
                    omegae_dot = this.cc.gps.ORBITAL_P.OMEGAE_DOT;
                case 'R'
                    omegae_dot = this.cc.glo.ORBITAL_P.OMEGAE_DOT;
                case 'E'
                    omegae_dot = this.cc.gal.ORBITAL_P.OMEGAE_DOT;
                case 'C'
                    omegae_dot = this.cc.bds.ORBITAL_P.OMEGAE_DOT;
                case 'J'
                    omegae_dot = this.cc.qzs.ORBITAL_P.OMEGAE_DOT;
                case 'I'
                    omegae_dot = this.cc.irn.ORBITAL_P.OMEGAE_DOT;
                otherwise
                    Logger.getInstance().addWarning('Something went wrong in satellite_positions.m\nUnrecognized Satellite system!\n');
                    omegae_dot = this.cc.gps.ORBITAL_P.OMEGAE_DOT;
            end
            omega_tau = omegae_dot * travel_time;
            xR  = [cos(omega_tau)    sin(omega_tau)];
            yR  = [-sin(omega_tau)    cos(omega_tau)];
            XS_r(:,1) = sum(xR .* XS(:,1:2),2); % X
            XS_r(:,2) = sum(yR .* XS(:,1:2),2); % Y
            XS_r(:,3) = XS(:,3); % Z
        end
        
        function removeAllCorrections(this)
            this.removeDtSat();
            this.removeGroupDelay();
            this.removePCV();
            this.removePoleTide();
            this.removePhaseWindUpCorr();
            this.removeSolidEarthTide();
            this.removeShDelay();
            this.removeOceanLoading();
        end
        
        function applyAllCorrections(this)
            this.applyDtSat();
            this.applyGroupDelay();
            this.applyPCV();
            this.applyPoleTide();
            this.applyPhaseWindUpCorr();
            this.applySolidEarthTide();
            this.applyShDelay();
            this.applyOceanLoading();
        end
    end
    
    % ==================================================================================================================================================
    %  FUNCTIONS TO COMPUTE APPLY AND REMOVE VARIOUS MODELED CORRECTIONS
    %  NOTE: Methods related to corrections that are applied to the observables and whose
    %  values is not stored separately (All except Iono and tropo) are
    %  structured as follow:
    %   - A : add or subtract the correction to the observations
    %   - applyA : a warpper of A to add the correction
    %   - removeA : a wrapper of A to subtract the correction
    %   - computeA : does the numerical computation of A along the los
    %  where A is the name of the correction
    % ==================================================================================================================================================
    methods
        %--------------------------------------------------------
        % Azimuth and Elevation
        % -------------------------------------------------------
        
        function updateAzimuthElevation(this, sat)
            % Upadte azimute elevation into.sat
            % SYNTAX: 
            %   this.updateAzimuthElevation(<sat>)
            
            if nargin < 2
                for i = 1 : this.cc.getNumSat()
                    if sum(this.go_id == i) > 0
                        this.updateAzimuthElevation(i);
                    end
                end
            else
                if isempty(this.sat.el)
                    this.sat.el = zeros(size(this.sat.avail_index));
                end
                if isempty(this.sat.az)
                    this.sat.az = zeros(size(this.sat.avail_index));
                end
                av_idx = this.sat.avail_index(:, sat);
                [this.sat.az(av_idx, sat), this.sat.el(av_idx, sat)] = this.computeAzimuthElevation(sat);
            end
        end
        
        function [az, el] = computeAzimuthElevation(this, sat)
            XS = this.getXSTxRot(sat);
            [az, el] = this.computeAzimuthElevationXS(XS);
        end
        function [az] = getAz(this, sat)
            %DESCRIPTION: get valid azimuth for given satellite
            az = this.sat.az(this.sat.avail_index(:,sat),sat);
        end
        function [el] = getEl(this, sat)
            %DESCRIPTION: get valid elevation for given satellite
            el = this.sat.el(this.sat.avail_index(:,sat),sat);
        end
        function [az, el] = computeAzimuthElevationXS(this, XS, XR)
            % SYNTAX:
            %   [az, el] = this.computeAzimuthElevationXS(XS)
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
                    this.log.addError('[ computeAzimuthElevationXS ] Number of satellite positions differ from number of receiver positions');
                    return
                end
            else
                XR = repmat(this.xyz(1,:),n_epoch,1);
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
            
            az(~zero_idx) = atan2d(e(~zero_idx), n(~zero_idx));
            el(~zero_idx) = atan2d(u(~zero_idx), hor_dist(~zero_idx));
        end
        
        %--------------------------------------------------------
        % Iono and Tropo
        % -------------------------------------------------------
        
        function  updateErrTropo(this, sat, flag)
            %INPUT:
            % sat : number of sat
            % flag: flag of the tropo model
            %DESCRIPTION: update the tropospheric correction
            atmo = Atmosphere();
            
            if isempty(this.sat.err_tropo)
                this.sat.err_tropo = zeros(size(this.sat.avail_index));
            end
            
            if nargin < 2 || strcmp(sat,'all')
                this.log.addMessage(this.log.indent('Updating tropospheric errors',6))
                if nargin < 3
                    flag = this.state.tropo_model;
                end
                sat = 1 : size(this.sat.avail_index,2);
            end
            
            this.sat.err_tropo(:, sat) = 0;
            
            %%% compute lat lon
            [~, lon_full, h_full, lat_full] = cart2geod(this.xyz(:,1), this.xyz(:,2), this.xyz(:,3));
            
            if nargin < 3
                flag = this.state.tropo_model;
            end
            
            if flag > 0
                geoid = Go_State.getInstance.getRefGeoid();
                if geoid.ncols > 0
                    % geoid ondulation interpolation
                    undu = getOrthometricCorr(lat_full(end), lon_full(end), geoid); % consider geoid undulation constant
                else
                    undu = [];
                end
            end
            
            for s = sat
                idx = this.sat.avail_index(:, s) > 0;
                if sum(idx)>0
                    
                    if size(this.xyz,1) > 1 % is dynamic
                        %                         XR = this.xyz(idx,:);
                        %                         [az, el] = this.computeAzimuthElevationXS(XS, XR);
                        h = h_full(idx);
                        lat = lat_full(idx);
                        lon = lon_full(idx);
                    else  % is static
                        %                         [az, el] = this.computeAzimuthElevationXS(XS);
                        h = h_full;
                        lat = lat_full;
                        lon = lon_full;
                    end
                    %%% get compute az el
                    %az = this.getAz(s);
                    el = this.getEl(s);
                    switch flag
                        case 0 % no model
                            
                        case 1 % Saastamoinen with standard atmosphere
                            if isempty(undu)
                                this.log.addWarning('Geoid not found = using undulation = 0');
                                undu = 0;
                            end
                            this.sat.err_tropo(idx, s) = atmo.saastamoinen_model(h, undu, el);
                            
                        case 2 % Saastamoinen with GPT
                            gps_time = this.time.getGpsTime();
                            lat_t = zeros(size(idx)); lon_t = zeros(size(idx)); h_t = zeros(size(idx)); el_t = zeros(size(idx));
                            lat_t(idx) = lat; lon_t(idx) = lon; h_t(idx) = h; el_t(idx) = el;
                            for e = 1 : size(idx,1)
                                if idx(e) > 0
                                    this.sat.err_tropo(e, s) = atmo.saastamoinen_model_GPT(gps_time(e), lat_t(e) / pi * 180, lon_t(e) / pi * 180, h_t(e), undu, el_t(e));
                                end
                            end
                    end
                end
            end
            
        end
        
        function updateErrIono(this, sat)
            if isempty(this.sat.err_iono)
                this.sat.err_iono = size(this.sat.avail_index);
            end
            if nargin < 2
                this.log.addMessage(this.log.indent('Updating ionospheric errors',6))
                for s = 1 : size(this.sat.avail_index,2)
                    this.updateErrIono(s);
                end
            else
                idx = this.sat.avail_index(:,sat) > 0; % epoch for which satellite is present
                if sum(idx) > 0
                    
                    XS = this.sat.cs.coordInterpolate(this.time.getSubSet(idx), sat);
                    %%% compute lat lon
                    if size(this.xyz,1) > 1
                        [~, lat, ~, lon] = cart2geod(this.xyz(idx, 1), this.xyz(idx, 2), this.xyz(idx, 3));
                    else
                        [~, lat, ~, lon] = cart2geod(this.xyz(1), this.xyz(2), this.xyz(3));
                    end
                    %%% compute az el
                    if size(this.xyz,1)>1
                        [az, el] = this.computeAzimuthElevationXS(this.xyz(idx,:) ,XS);
                    else
                        [az, el] = this.computeAzimuthElevationXS(XS);
                    end
                    
                    
                    switch this.state.iono_model
                        case 0 %no model
                            this.sat.err_iono(idx,sat) = zeros(size(el));
                        case 1 %Geckle and Feen model
                            %corr = simplified_model(lat, lon, az, el, mjd);
                        case 2 %Klobuchar model
                            [week, sow] = time2weektow(this.time.getSubSet(idx).getGpsTime());
                            if ~isempty(this.sat.cs.iono )
                                this.sat.err_iono(idx,sat) = Atmosphere.klobuchar_model(lat, lon, az, el, sow, this.sat.cs.iono);
                            else
                                this.log.addWarning('No klobuchar parameter found, iono correction not computed',100);
                            end
                            
                    end
                end
            end
        end
        
        
        %--------------------------------------------------------
        % Solid earth tide
        % -------------------------------------------------------
        
        
        function solidEarthTide(this,sgn)
            % DESCRIPTION: add or subtract ocean loading from observations
            et_corr = this.computeSolidTideCorr();
            
            for s = 1 : this.cc.getNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        if  this.obs_code(o,1) == 'L'
                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* et_corr(o_idx,s)' ./ this.wl(o);
                        else
                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* et_corr(o_idx,s)';
                        end
                    end
                end
            end
        end
        
        function applySolidEarthTide(this)
            if this.et_delay_status == 0
                this.log.addMarkedMessage('Applying Solid Earth Tide corrections');
                this.solidEarthTide(1);
                this.et_delay_status = 1; %applied
            end
        end
        
        function removeSolidEarthTide(this)
            if this.et_delay_status == 1
                this.log.addMarkedMessage('Removing Solid Earth Tide corrections');
                this.solidEarthTide(-1);
                this.et_delay_status = 0; %not applied
            end
        end
        
        function solid_earth_corr = computeSolidTideCorr(this, sat)
            
            % SYNTAX:
            %   [stidecorr] = this.getSolidTideCorr();
            %
            % INPUT:
            %
            % OUTPUT:
            %   stidecorr = solid Earth tide correction terms (along the satellite-receiver line-of-sight)
            %
            % DESCRIPTION:
            %   Computation of the solid Earth tide displacement terms.
            if nargin < 2
                sat  = 1 : this.cc.getNumSat();
            end
            solid_earth_corr = zeros(this.time.length, length(sat));
            XR = this.getXR();
            
            [~, lam, ~, phiC] = cart2geod(XR(1,1), XR(1,2), XR(1,3));
            %north (b) and radial (c) local unit vectors
            b = [-sin(phiC)*cos(lam); -sin(phiC)*sin(lam); cos(phiC)];
            c = [+cos(phiC)*cos(lam); +cos(phiC)*sin(lam); sin(phiC)];
            
            %interpolate sun moon and satellites
            time = this.time.getCopy;
            [X_sun, X_moon]  = this.sat.cs.sunMoonInterpolate(time);
            XS               = this.sat.cs.coordInterpolate(time, sat);
            %receiver geocentric position
            
            
            XR_u = rowNormalize(XR);
            
            %sun geocentric position
            X_sun_n = repmat(sqrt(sum(X_sun.^2,2)),1,3);
            X_sun_u = X_sun ./ X_sun_n;
            
            %moon geocentric position
            X_moon_n = repmat(sqrt(sum(X_moon.^2,2)),1,3);
            X_moon_u = X_moon ./ X_moon_n;
            
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
            Vsun  = repmat(sum(conj(X_sun_u) .* XR_u, 2),1,3);
            Vmoon = repmat(sum(conj(X_moon_u) .* XR_u, 2),1,3);
            r_sun2  = (GS*R^4)./(GE*X_sun_n.^3) .*(H2.*XR_u.*(1.5*Vsun.^2  - 0.5)+ 3*L2*Vsun .*(X_sun_u  - Vsun .*XR_u));
            r_moon2 = (GM*R^4)./(GE*X_moon_n.^3).*(H2.*XR_u.*(1.5*Vmoon.^2 - 0.5) + 3*L2*Vmoon.*(X_moon_u - Vmoon.*XR_u));
            r = r_sun2 + r_moon2;
            
            %nominal degree 3 Love number
            H3 = 0.292;
            %nominal degree 3 Shida number
            L3 = 0.015;
            
            %solid Earth tide displacement (degree 3)
            r_sun3  = (GS.*R^5)./(GE.*X_sun_n.^4) .*(H3*XR_u.*(2.5.*Vsun.^3  - 1.5.*Vsun)  +   L3*(7.5*Vsun.^2  - 1.5).*(X_sun_u  - Vsun .*XR_u));
            r_moon3 = (GM.*R^5)./(GE.*X_moon_n.^4).*(H3*XR_u.*(2.5.*Vmoon.^3 - 1.5.*Vmoon) +   L3*(7.5*Vmoon.^2 - 1.5).*(X_moon_u - Vmoon.*XR_u));
            r = r + r_sun3 + r_moon3;
            
            %from "conventional tide free" to "mean tide"
            radial = (-0.1206 + 0.0001*p)*p;
            north  = (-0.0252 + 0.0001*p)*sin(2*phiC);
            r = r + repmat([radial*c + north*b]',time.length,1);
            
            %displacement along the receiver-satellite line-of-sight
            
            for i  = 1 : length(sat)
                s = sat(i);
                sat_idx = this.sat.avail_index(:,s);
                az = this.getAz(s);
                el = this.getEl(s);
                LOSu = [cosd(el).*sind(az) cosd(el).*cosd(az) sind(el)];
                % oceanloadcorr(s,1) = dot(corrXYZ,LOSu);
                solid_earth_corr(sat_idx,i) = sum(conj(r(sat_idx,:)).*LOSu,2);
            end
            
        end
        
        function updateSolidEarthCorr(this, sat)
            %DESCRIPTION: upadte the correction related to solid earth
            % solid tides, ocean loading, pole tides.
            if isempty(this.sat.solid_earth_corr)
                this.log.addMessage(this.log.indent('Updating solid earth corrections',6))
                this.sat.solid_earth_corr = zeros(size(this.sat.avail_index));
            end
            if nargin < 2
                
                for s = 1 : size(this.sat.avail_index,2)
                    this.updateSolidEarthCorr(s);
                end
            else
                this.sat.solid_earth_corr(:,sat) = this.computeSolidTideCorr(sat);% + this.computeOceanLoading(sat) + this.getPoleTideCorr(sat);
            end
        end
        
        %--------------------------------------------------------
        % Ocean Loading
        % -------------------------------------------------------
        
        function oceanLoading(this,sgn)
            % DESCRIPTION: add or subtract ocean loading from observations
            ol_corr = this.computeOceanLoading();
            if isempty(ol_corr)
                this.log.addWarning('No ocean loading displacements matrix present for the receiver')
                return
            end
            
            for s = 1 : this.cc.getNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        if  this.obs_code(o,1) == 'L'
                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* ol_corr(o_idx,s)' ./ this.wl(o);
                        else
                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* ol_corr(o_idx,s)';
                        end
                    end
                end
            end
        end
        
        function applyOceanLoading(this)
            if this.ol_delay_status == 0
                this.log.addMarkedMessage('Applying Ocean Loading corrections');
                this.oceanLoading(1);
                this.ol_delay_status = 1; %applied
            end
        end
        
        function removeOceanLoading(this)
            if this.ol_delay_status == 1
                this.log.addMarkedMessage('Removing Ocean Loading corrections');
                this.oceanLoading(-1);
                this.ol_delay_status = 0; %not applied
            end
        end
        
        function [ocean_load_corr] = computeOceanLoading(this, sat) % WARNING: to be tested
            
            % SYNTAX:
            %   [oceanloadcorr] = ocean_loading_correction(time, XR, XS);
            %
            % INPUT:
            %
            % OUTPUT:
            %   oceanloadcorr = ocean loading correction terms (along the satellite-receiver line-of-sight)
            %
            % DESCRIPTION:
            %   Computation of the ocean loading displacement terms.ù
            % NOTES:
            %  Inefficent to compute them separate by staellites always call
            %  as a block
            
            if nargin < 2
                sat = 1 : this.cc.getNumSat();
            end
            
            
            ocean_load_corr = zeros(this.time.length,length(sat));
            if (isempty(this.ocean_load_disp)) || this.ocean_load_disp.available == 0
                return
            end
            
            ol_disp = this.ocean_load_disp;
            
            %terms depending on the longitude of the lunar node (see Kouba and Heroux, 2001)
            fj = 1; %(at 1-3 mm precision)
            uj = 0; %(at 1-3 mm precision)
            
            %ref: http://202.127.29.4/cddisa/data_base/IERS/Convensions/Convension_2003/SUBROUTINES/ARG.f
            tidal_waves = [1.40519E-4, 2.0,-2.0, 0.0, 0.00; ... % M2  - semidiurnal
                1.45444E-4, 0.0, 0.0, 0.0, 0.00; ... % S2  - semidiurnal
                1.37880E-4, 2.0,-3.0, 1.0, 0.00; ... % N2  - semidiurnal
                1.45842E-4, 2.0, 0.0, 0.0, 0.00; ... % K2  - semidiurnal
                0.72921E-4, 1.0, 0.0, 0.0, 0.25; ... % K1  - diurnal
                0.67598E-4, 1.0,-2.0, 0.0,-0.25; ... % O1  - diurnal
                0.72523E-4,-1.0, 0.0, 0.0,-0.25; ... % P1  - diurnal
                0.64959E-4, 1.0,-3.0, 1.0,-0.25; ... % Q1  - diurnal
                0.53234E-5, 0.0, 2.0, 0.0, 0.00; ... % Mf  - long-period
                0.26392E-5, 0.0, 1.0,-1.0, 0.00; ... % Mm  - long-period
                0.03982E-5, 2.0, 0.0, 0.0, 0.00];    % Ssa - long-period
            
            refdate = datenum([1975 1 1 0 0 0]);
            
            [week, sow] = time2weektow(this.time.getGpsTime);
            dateUTC = datevec(gps2utc(datenum(gps2date(week, sow))));
            
            %separate the fractional part of day in seconds
            fday = dateUTC(:,4)*3600 + dateUTC(:,5)*60 + dateUTC(:,6);
            dateUTC(:,4:end) = 0;
            
            %number of days since reference date (1 Jan 1975)
            days = (datenum(dateUTC) - refdate);
            
            capt = (27392.500528 + 1.000000035*days)/36525;
            
            %mean longitude of the Sun at the beginning of day
            H0 = (279.69668 + (36000.768930485 + 3.03e-4.*capt).*capt).*pi/180;
            
            %mean longitude of the Moon at the beginning of day
            S0 = (((1.9e-6*capt - 0.001133).*capt + 481267.88314137).*capt + 270.434358).*pi/180;
            
            %mean longitude of the lunar perigee at the beginning of day
            P0 = (((-1.2e-5.*capt - 0.010325).*capt + 4069.0340329577).*capt + 334.329653)*pi/180;
            
            corr = zeros(this.time.length,3);
            corrENU = zeros(this.time.length,3);
            for k = 1 : 11
                angle = tidal_waves(k,1)*fday + tidal_waves(k,2)*H0 + tidal_waves(k,3)*S0 + tidal_waves(k,4)*P0 + tidal_waves(k,5)*2*pi;
                corr  = corr + repmat(fj*ol_disp.matrix(1:3,k)',this.time.length,1).*cos(repmat(angle,1,3) + uj - repmat(ol_disp.matrix(4:6,k)'*pi/180,this.time.length,1));
            end
            corrENU(:,1) = -corr(:,1); %east
            corrENU(:,2) = -corr(:,2); %north
            corrENU(:,3) =  corr(:,3); %up
            
            % get XR
            XR = this.getXR();
            %displacement along the receiver-satellite line-of-sight
            XRcorr = local2globalPos(corrENU', XR')';
            corrXYZ = XRcorr - XR;
            
            for i  = 1 : length(sat)
                s = sat(i);
                sat_idx = this.sat.avail_index(:,s);
                az = this.getAz(s) / 180 *pi;
                el = this.getEl(s) / 180 * pi;
                LOSu = [cos(el).*sin(az) cos(el).*cos(az) sin(el)];
                % oceanloadcorr(s,1) = dot(corrXYZ,LOSu);
                ocean_load_corr(sat_idx,i) = sum(corrXYZ(sat_idx,:).*LOSu, 2);
            end
        end
        
        function importOceanLoading(this)
            %DESCRIPTION: load ocean loading displcement matrix from
            %ocean_loading.blq if satation is present
            this.ocean_load_disp = load_BLQ( this.state.getOceanFile,{this.getShortName()});
        end
        
        %--------------------------------------------------------
        % Pole Tide
        % -------------------------------------------------------
        
        
        function poleTide(this,sgn)
            % DESCRIPTION: add or subtract ocean loading from observations
            pt_corr = this.computePoleTideCorr();
            if isempty(pt_corr)
                this.log.addWarning('No ocean loading displacements matrix present for the receiver')
                return
            end
            
            for s = 1 : this.cc.getNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        if  this.obs_code(o,1) == 'L'
                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* pt_corr(o_idx,s)' ./ this.wl(o);
                        else
                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* pt_corr(o_idx,s)';
                        end
                    end
                end
            end
        end
        
        function applyPoleTide(this)
            if this.pt_delay_status == 0
                this.log.addMarkedMessage('Applying Pole Tide corrections');
                this.poleTide(1);
                this.pt_delay_status = 1; %applied
            end
        end
        
        function removePoleTide(this)
            if this.pt_delay_status == 1
                this.log.addMarkedMessage('Removing Pole Tide corrections');
                this.poleTide(-1);
                this.pt_delay_status = 0; %not applied
            end
        end
        
        
        function [pole_tide_corr] = computePoleTideCorr(this, sat)
            
            % SYNTAX:
            %   [poletidecorr] = pole_tide_correction(time, XR, XS, SP3, phiC, lam);
            %
            % INPUT:
            %
            % OUTPUT:
            %   poletidecorr = pole tide correction terms (along the satellite-receiver line-of-sight)
            %
            % DESCRIPTION:
            %   Computation of the pole tide displacement terms.
            
            if nargin < 2
                sat = 1 : this.cc.getNumSat();
            end
            
            XR = this.getXR;
            [~, lam, ~, phiC] = cart2geod(XR(1,1), XR(2,1), XR(3,1));
            
            pole_tide_corr = zeros(this.time.length,length(sat));
            erp  = this.sat.cs.erp;
            %interpolate the pole displacements
            if (~isempty(erp))
                if (length(erp.t) > 1)
                    m1 = interp1(erp.t, erp.m1, this.time.getGpsTime, 'linear', 'extrap');
                    m2 = interp1(erp.t, erp.m2, this.time.getGpsTime, 'linear', 'extrap');
                else
                    m1 = repmat(erp.m1,this.time.length,1);
                    m2 = repmat(erp.m2,this.time.length,1);
                end
                
                deltaR   = -33*sin(2*phiC).*(m1.*cos(lam) + m2.*sin(lam))*1e-3;
                deltaLam =  9* cos(  phiC).*(m1.*sin(lam) - m2.*cos(lam))*1e-3;
                deltaPhi = -9* cos(2*phiC).*(m1.*cos(lam) + m2.*sin(lam))*1e-3;
                
                corrENU(:,1) = deltaLam; %east
                corrENU(:,2) = deltaPhi; %north
                corrENU(:,3) = deltaR;   %up
                
                %displacement along the receiver-satellite line-of-sight
                XRcorr = local2globalPos(corrENU', XR')';
                corrXYZ = XRcorr - XR;
                for i  = 1 : length(sat)
                    s = sat(i);
                    sat_idx = this.sat.avail_index(:,s);
                    az = this.getAz(s) / 180 *pi;
                    el = this.getEl(s) / 180 * pi;
                    LOSu = [cos(el).*sin(az) cos(el).*cos(az) sin(el)];
                    % oceanloadcorr(s,1) = dot(corrXYZ,LOSu);
                    pole_tide_corr(sat_idx,i) = sum(corrXYZ(sat_idx,:).*LOSu, 2);
                end
            end
            
        end
        
        %--------------------------------------------------------
        % Phase Wind up
        % -------------------------------------------------------
        
        
        
        
        function phaseWindUpCorr(this,sgn)
            % DESCRIPTION: add or subtract ocean loading from observations
            ph_wind_up = this.computePhaseWindUp();
            
            for s = 1 : this.cc.getNumSat()
                obs_idx = this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* ph_wind_up(o_idx,s)'; % add or subtract???
                    end
                end
            end
        end
        
        function applyPhaseWindUpCorr(this)
            if this.pw_delay_status == 0
                this.log.addMarkedMessage('Applying Phase Wind Up corrections');
                this.phaseWindUpCorr(1);
                this.pw_delay_status = 1; %applied
            end
        end
        
        function removePhaseWindUpCorr(this)
            if this.pw_delay_status == 1
                this.log.addMarkedMessage('Removing Phase Wind Up corrections');
                this.phaseWindUpCorr(-1);
                this.pw_delay_status = 0; %not applied
            end
        end
        
        function [ph_wind_up] = computePhaseWindUp(this)
            
            %east (a) and north (b) local unit vectors
            XR = this.getXR();
            [phi, lam] = cart2geod(XR(:,1), XR(:,2), XR(:,3));
            a = [-sin(lam) cos(lam) zeros(size(lam))];
            b = [-sin(phi).*cos(lam) -sin(phi).*sin(lam) cos(phi)];
            s_time = this.time.getCopy();%SubSet(av_idx);
            s_a = a;%(av_idx,:);
            s_b = b;%(av_idx,:);
            
            
            sat = 1: this.cc.getNumSat();
            
            [x, y, z] = this.sat.cs.getSatFixFrame(s_time);
            ph_wind_up = zeros(this.time.length,length(sat));
            for s = sat
                av_idx = this.sat.avail_index(:,s);
                i_s = squeeze(x(:,s,:));
                j_s = squeeze(y(:,s,:));
                k_s = squeeze(z(:,s,:));
                
                %receiver and satellites effective dipole vectors
                Dr = s_a - k_s.*repmat(sum(k_s.*s_a,2),1,3) + cross(k_s,s_b);
                Ds = i_s - k_s.*repmat(sum(k_s.*i_s,2),1,3) - cross(k_s,j_s);
                
                %phase wind-up computation
                psi = sum(conj(k_s) .* cross(Ds, Dr),2);
                arg = sum(conj(Ds) .* Dr,2) ./ (sqrt(sum(Ds.^2,2)) .* sqrt(sum(Dr.^2,2)));
                arg(arg < -1) = -1;
                arg(arg > 1) = 1;
                dPhi = sign(psi).*acos(arg)./(2*pi);
                % find jump
                ddPhi = diff(dPhi);
                jumps = find(abs(ddPhi) > 0.9);
                jumps_sgn = sign(ddPhi(jumps));
                for t = 1 : length(jumps)
                    jump = jumps(t);
                    dPhi((jump+1):end) = dPhi((jump+1):end) - jumps_sgn(t);
                end
                ph_wind_up(av_idx,s) = dPhi(av_idx);
            end
        end
        
        %--------------------------------------------------------
        % Shapiro Delay
        % -------------------------------------------------------
        
        
        
        function shDelay(this,sgn)
            % DESCRIPTION: add or subtract shapiro delay from observations
            for s = 1 : this.cc.getNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    freqs = unique(str2num(this.obs_code(obs_idx,2)));
                    freqs = reshape(freqs,1,length(freqs));
                    for f = freqs
                        obs_idx_f = obs_idx & this.obs_code(:,2) == num2str(f);
                        sh_delays = this.computeShapirodelay(s);
                        for o = find(obs_idx_f)'
                            pcv_idx = this.obs(o, this.sat.avail_index(:, s)) ~=0; %find which correction to apply
                            o_idx = this.obs(o, :) ~=0; %find where apply corrections
                            if  this.obs_code(o,1) == 'L'
                                this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* sh_delays(pcv_idx)' ./ this.wl(o);
                            else
                                this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* sh_delays(pcv_idx)';
                            end
                        end
                    end
                end
            end
        end
        
        function applyShDelay(this)
            if this.sh_delay_status == 0
                this.log.addMarkedMessage('Applying Shapiro delay corrections');
                this.shDelay(1);
                this.sh_delay_status = 1; %applied
            end
        end
        
        function removeShDelay(this)
            if this.sh_delay_status == 1
                this.log.addMarkedMessage('Removing Shapiro delay corrections');
                this.shDelay(-1);
                this.sh_delay_status = 0; %not applied
            end
        end
        function [sh_delay] = computeShapirodelay(this, sat)
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
            XS = this.getXSTxRot(sat);
            if this.isStatic()
                XR = repmat(this.xyz,size(XS,1),1);
            else
                XR = this.xyz(this.sat.avail_index(:,sat),:);
            end
            
            distR = sqrt(sum(XR.^2 ,2));
            distS = sqrt(sum(XS.^2 ,2));
            
            distSR = sqrt(sum((XS-XR).^2 ,2));
            
            
            GM = 3.986005e14;
            
            corr = 2*GM/(goGNSS.V_LIGHT^2) * log((distR + distS + distSR)./(distR + distS - distSR)); %#ok<CPROPLC>
            
            sh_delay = 2*GM/(goGNSS.V_LIGHT^2)*log((distR + distS + distSR)./(distR + distS - distSR));
            
        end
        
        %--------------------------------------------------------
        % PCV
        % -------------------------------------------------------
        
        function PCV(this,sgn)
            % DESCRIPTION: correct measurement for PCV both of receiver
            % antenna and satellite antenna
            if ~isempty(this.pcv) || ~isempty(this.sat.cs.ant_pcv)
                this.updateAllAvailIndex();
                % getting sat - receiver vector for each epoch
                XR_sat = - this.getXSLoc();
                
                if ~isempty(this.pcv) && ~isempty(this.pcv.name)
                    %TODO correct for receiver pcv
                    f_code_history = []; % save f_code checked to print only one time the warning message
                    for s = 1 : size(this.sat.el,2)
                        sat_idx = this.sat.avail_index(:,s);
                        el = this.sat.el(sat_idx,s);
                        az = this.sat.el(sat_idx,s);
                        neu_los = [cosd(az).*cosd(el) sind(az).*cosd(el) sind(el)];
                        obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                        obs_idx = obs_idx & this.go_id == s;
                        if sum(obs_idx) > 0
                            freqs = unique(str2num(this.obs_code(obs_idx,2)));
                            for f = freqs'
                                obs_idx_f = obs_idx & this.obs_code(:,2) == num2str(f);
                                sys = this.system(obs_idx_f);
                                f_code = [sys(1) sprintf('%02d',f)];
                                
                                pco_idx = idxCharLines(this.pcv.frequency_name,f_code);
                                if sum(pco_idx)
                                    pco = this.pcv.offset(:,:,pco_idx)';
                                    pco_delays = neu_los*pco;
                                    pcv_delays = pco_delays + this.getPCV(pco_idx,el,az);
                                    for o = find(obs_idx_f)'
                                        pcv_idx = this.obs(o, this.sat.avail_index(:, s)) ~=0; %find which correction to apply
                                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                                        if  this.obs_code(o,1) == 'L'
                                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* pcv_delays(pcv_idx)' ./ this.wl(o);
                                        else
                                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* pcv_delays(pcv_idx)';
                                        end
                                    end
                                else
                                    if isempty(f_code_history) || ~sum(idxCharLines(f_code_history,f_code))
                                        this.log.addMessage(this.log.indent(sprintf('No corrections found for antenna model %s and frequency %s',this.ant_type,f_code),6));
                                        f_code_history = [f_code_history;f_code];
                                    end
                                end
                            end
                        end
                    end
                end
                if ~isempty(this.sat.cs.ant_pcv)
                    
                    % getting satellite reference frame for each epoch
                    [x, y, z] = this.sat.cs.getSatFixFrame(this.time);                    
                    % transgform into satellite reference system
                    XR_sat = cat(3,sum(XR_sat.*x,3),sum(XR_sat.*y,3),sum(XR_sat.*z,3));
                    
                    % getting az and el
                    distances = sqrt(sum(XR_sat.^2,3));
                    XR_sat_norm = XR_sat ./ repmat(distances,1,1,3);
                    
                    az = atan2(XR_sat_norm(:, :, 2),XR_sat_norm(:, :, 1)); % here azimuth is intended as angle from x axis
                    az(az<0) = az(az<0) + 2*pi;
                    el = atan2(XR_sat_norm(:, :, 3), sqrt(sum(XR_sat_norm(:, :, 1:2).^2, 3)));
                    
                    % getting pscv and applying it to the obs
                    for s = 1 : size(el,2)
                        obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                        obs_idx = obs_idx & this.go_id == s;
                        if sum(obs_idx) > 0
                            freqs = unique(str2num(this.obs_code(obs_idx,2)));
                            for f = freqs'
                                obs_idx_f = obs_idx & this.obs_code(:,2) == num2str(f);
                                az_idx = ~isnan(az(:,s));
                                az_tmp = az(az_idx,s) / pi * 180;
                                el_idx = ~isnan(el(:,s));
                                el_tmp = el(el_idx,s) / pi * 180;
                                pcv_delays = this.sat.cs.getPCV( f, s, el_tmp, az_tmp);
                                for o = find(obs_idx_f)'
                                    pcv_idx = this.obs(o, this.sat.avail_index(:, s)) ~=0; %find which correction to apply
                                    o_idx = this.obs(o, :) ~=0; %find where apply corrections
                                    if  this.obs_code(o,1) == 'L'
                                        this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* pcv_delays(pcv_idx)' ./ this.wl(o);
                                    else
                                        this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)* pcv_delays(pcv_idx)';
                                    end
                                end
                            end
                        end
                    end
                    
                end
            end
        end
        
        function pcv_delay = getPCV(this, idx, el, az) %Duplicate of method in Core_Sky consider restructuring
            % DESCRIPTION: get the pcv correction for a given satellite and a given
            % azimuth and elevations using linear or bilinear interpolation
            
             pcv_delay = zeros(size(el));

            pcv = this.pcv;
            
            %tranform el in zen
            zen = 90 - el;
            % get el idx
            zen_pcv = pcv.tablePCV_zen;
            
            min_zen = zen_pcv(1);
            max_zen = zen_pcv(end);
            d_zen = (max_zen - min_zen)/(length(zen_pcv)-1);
            zen_idx = min(max(floor((zen - min_zen)/d_zen) + 1 , 1),length(zen_pcv) - 1);
            d_f_r_el = min(max(zen_idx*d_zen - zen, 0)/ d_zen, 1) ;
            if nargin < 4 || isempty(pcv.tablePCV_azi) %no azimuth change
                pcv_val = pcv.tableNOAZI(:,:,idx); %etract the right frequency
                
                pcv_delay = d_f_r_el .* pcv_val(zen_idx)' + (1 - d_f_r_el) .* pcv_val(zen_idx + 1)';
            else
               pcv_val = pcv.tablePCV(:,:,idx); %etract the right frequency
               
               %find azimuth indexes
               az_pcv = pcv.tablePCV_azi;
                min_az = az_pcv(1);
                max_az = az_pcv(end);
                d_az = (max_az - min_az)/(length(az_pcv)-1);
                az_idx = min(max(floor((az - min_az)/d_az) + 1, 1),length(az_pcv) - 1);
                d_f_r_az = min(max(az - (az_idx-1)*d_az, 0)/d_az, 1); 
                
                %interpolate along zenital angle
                idx1 = sub2ind(size(pcv_val),az_idx,zen_idx);
                idx2 = sub2ind(size(pcv_val),az_idx,zen_idx+1);
                pcv_delay_lf =  d_f_r_el .* pcv_val(idx1) + (1 - d_f_r_el) .* pcv_val(idx2);
                idx1 = sub2ind(size(pcv_val),az_idx+1,zen_idx);
                idx2 = sub2ind(size(pcv_val),az_idx+1,zen_idx+1);
                pcv_delay1_rg = d_f_r_el .* pcv_val(idx1) + (1 - d_f_r_el) .* pcv_val(idx2);
                %interpolate alogn azimtuh
                pcv_delay = d_f_r_az .* pcv_delay_lf + (1 - d_f_r_az) .* pcv_delay1_rg;
            end
        end
        
        function applyPCV(this)
            if this.pcv_delay_status == 0
                this.log.addMarkedMessage('Applying PCV corrections');
                this.PCV(1);
                this.pcv_delay_status = 1; %applied
            end
        end
        
        function removePCV(this)
            if this.pcv_delay_status == 1
                this.log.addMarkedMessage('Removing PCV corrections');
                this.PCV(-1);
                this.pcv_delay_status = 0; %applied
            end
        end
        
    end
    
    % ==================================================================================================================================================
    %  LEGACY IMPORTER of goGPS 5.x RESULTS
    % ==================================================================================================================================================
    methods
        function legacyImportResults(this, file_prefix, run_start, run_stop)
            % Import after reset a position and tropo file (if present)
            %
            % SYNTAX:  
            %   this.legacyImportResults(file_prefix, <run_start>, <run_stop>)
            %
            % INPUT:
            %   file_name     it could include the key ${RUN} that will be substituted with a 3 digits number containing the run, from run_start to run_stop
            %   run_start     number of the first run to load
            %   run_stop      number of the last run to load
            %            
            if (nargin == 1) || isempty(file_prefix)
                [file_prefix, file_path] = uigetfile('*.txt', 'Select a _position.txt or _tropo.txt file');
                file_prefix = [file_path file_prefix];
            end
            this.reset();            
            if (length(file_prefix) > 13 && strcmp(file_prefix(end - 12 : end), '_position.txt'))
                file_prefix = file_prefix(1 : end - 13);
            end
            if (length(file_prefix) > 10 && strcmp(file_prefix(end - 9 : end), '_tropo.txt'))
                file_prefix = file_prefix(1 : end - 10);
            end
            
            if nargin < 4
                run_start = 0;
                run_stop = 0;
            end
            
            GPS_RUN = '${RUN}';
            r = 0;
            for run = run_start : run_stop
                r = r + 1;
                file_name = [strrep(file_prefix, GPS_RUN, sprintf('%03d', run)) '_position.txt'];
                this.log.addMessage(sprintf('Importing %s', file_name));
                if exist(file_name, 'file')
                    this.legacyAppendPosition(file_name);
                    
                    file_name = [strrep(file_prefix, GPS_RUN, sprintf('%03d', run)) '_tropo.txt'];
                    if exist(file_name, 'file')
                        this.log.addMessage(sprintf('Importing %s', file_name));
                        this.legacyAppendTropo(file_name)
                    else
                        this.log.addMessage(sprintf('Error loading the tropo file, it does not exists'));
                    end
                else
                    this.log.addMessage(sprintf('Error loading the position file, it does not exists'));
                end                
            end
        end
        
        function legacyAppendPosition (this, file)
            % import and append from a position file
            
            % Open position file as a string stream
            fid = fopen(file);
            txt = fread(fid,'*char')';
            txt(txt == 13) = [];
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  numel(txt)
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
            lim = [lim lim(:,2) - lim(:,1)];
            
            % corrupted lines
            ko_lines = find(lim(:, 3) ~= median(lim(:,3)));
            for l = numel(ko_lines) : -1 : 1
                txt(lim(ko_lines(l), 1) : lim(ko_lines(l), 2) + 1) = [];
            end
            
            % importing header informations
            eoh = 1; % end of header (the header of position files contains only 1 line)
            % File example:
            %    Date        GPS time           GPS week          GPS tow         Latitude        Longitude      h (ellips.)           ECEF X           ECEF Y           ECEF Z        UTM North         UTM East      h (orthom.)         UTM zone        Num. Sat.             HDOP            KHDOP      Local North       Local East          Local H    Ambiguity fix     Success rate              ZTD              ZWD              PWV
            %2017/04/03    00:00:00.000             1943        86400.000      45.80216141       9.09562643         291.5094     4398305.8406      704150.1081     4550153.9697     5072071.0952      507430.9212         244.4506             32 T               10            0.870            0.472           0.0000           0.0000           0.0000                0           0.0000          2.30187          0.04934          0.00000
            data = sscanf(txt(lim(2,1):end)','%4d/%2d/%2d    %2d:%2d:%6f             %4d %16f %16f %16f %16f %16f %16f %16f %16f %16f %16f %14d %c %16d %16f %16f %16f %16f %16f %16d %16f %16f %16f %16f\n');
            data = reshape(data, 30, numel(data)/30)';
            % import it as a GPS_Time obj
            if this.time.length() == 0
                this.time = GPS_Time(data(:,1:6));
            else
                this.time.append6ColDate(data(:,1:6));
            end
            
            lat = data(:,9);
            lon = data(:,10);
            h_ellips = data(:,11);
            h_ortho = data(:,17);
            
            xyz = data(:,12:14);
            enu = data(:,[16 15 17]);
            n_sat = data(:,20);
            hdop =  data(:,21);
            khdop = data(:,22);
            a_fix = data(:,26); 
            s_rate = data(:,27);
            
            % Append in obj
            this.xyz = [this.xyz; xyz];
            this.enu = [this.enu; enu];
            
            this.lat = [this.lat; lat];
            this.lon = [this.lon; lon];
            this.h_ellips = [this.h_ellips; h_ellips];
            this.h_ortho = [this.h_ortho; h_ortho];
            
            this.n_sat = [this.n_sat; n_sat];
            this.hdop = [this.hdop; hdop];
            this.khdop = [this.khdop; khdop];
            this.a_fix = [this.a_fix; a_fix];
            this.s_rate = [this.s_rate; s_rate];
            
            clear data;
        end
        
        function legacyAppendTropo (this, file)
            % import and append from a tropo file
            
            % Open tropo file as a string stream
            fid = fopen(file);
            txt = fread(fid,'*char')';
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  numel(txt)
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
            lim = [lim lim(:,2) - lim(:,1)];
            
            % importing header informations
            eoh = 1; % end of header (the header of tropo files contains only 1 line)
            
            % list and count satellites in view and not in view
            s = regexp(txt(lim(1,1) : lim(1,2)),'[GREJCIS]\d\d(?=-az)', 'match');
            num_sat = length(s);
            this.sat.id = reshape(cell2mat(s), 3, num_sat)';

            % extract all the epoch lines
            string_time = txt(repmat(lim(2:end,1),1,26) + repmat(0:25, size(lim,1)-1, 1))';
            % convert the times into a 6 col time
            date = cell2mat(textscan(string_time,'%4f/%2f/%2f    %2f:%2f:%6.3f'));
            
            % import it as a GPS_Time obj 
            % it should be imported from the _position file
            time = GPS_Time(date, [], 1);
            [~, id_int, id_ext] = intersect(round(this.time.getMatlabTime() * 86400 * 1e3), round(time.getMatlabTime() * 86400 * 1e3));
            n_epo = length(id_int);
            
            % extract all the ZHD lines
            string_zhd = txt(repmat(lim(2:end,1),1,17) + repmat(62:78, size(lim,1)-1, 1))';
            tmp = sscanf(string_zhd,'%f')'; clear string_zhd
            this.zhd(end + 1 : size(this.xyz, 1), 1) = nan;
            this.zhd(id_int, 1) = tmp(id_ext);
            
            % extract all the ZTD lines
            string_ztd = txt(repmat(lim(2:end,1),1,17) + repmat(78:94, size(lim,1)-1, 1))';
            tmp = sscanf(string_ztd,'%f'); clear string_ztd
            this.ztd(end + 1 : size(this.xyz, 1), 1) = nan;
            this.ztd(id_int, 1) = tmp(id_ext);
            
            % extract all the TGN lines
            string_tgn = txt(repmat(lim(2:end,1),1,17) + repmat(95:111, size(lim,1)-1, 1))';
            tmp = sscanf(string_tgn,'%f'); clear string_tgn
            this.tgn(end + 1 : size(this.xyz, 1), 1) = nan;
            this.tgn(id_int, 1) = tmp(id_ext);
            
            % extract all the TGE lines
            string_tge = txt(repmat(lim(2:end,1),1,17) + repmat(112:128, size(lim,1)-1, 1))';
            tmp = sscanf(string_tge,'%f'); clear string_tge
            this.tge(end + 1 : size(this.xyz, 1), 1) = nan;
            this.tge(id_int, 1) = tmp(id_ext);
            
            % extract all the ZWD lines
            string_zwd = txt(repmat(lim(2:end,1),1,17) + repmat(129:145, size(lim,1)-1, 1))';
            tmp = sscanf(string_zwd,'%f'); clear string_zwd
            this.zwd(end + 1 : size(this.xyz, 1), 1) = nan;
            this.zwd(id_int, 1) = tmp(id_ext);
            
            % extract all the PWV lines
            string_pwv = txt(repmat(lim(2:end,1),1,17) + repmat(146:162, size(lim,1)-1, 1))';
            tmp = sscanf(string_pwv,'%f'); clear string_pwv
            this.pwv(end + 1 : size(this.xyz, 1), 1) = nan;
            this.pwv(id_int, 1) = tmp(id_ext);
            
            %  extract all STD values if present
            slant_start = regexp(txt(lim(1,1) : lim(1,2)),'STD') - 6;
            num_sat = numel(slant_start);
            this.sat.slant_td(end + 1 : size(this.xyz, 1), 1 : num_sat) = nan;
            for s = 1 : numel(slant_start)
                tmp = sscanf(txt(bsxfun(@plus, repmat(slant_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1) - 1))', '%f');
                this.sat.slant_td(id_int, s) = tmp(id_ext);
            end
            
            % extract all azimuth and elevation lines in a matrix with 2 layers -
            % 1st is azimuth, 2nd is elevation 
            az_start = regexp(txt(lim(1,1) : lim(1,2)),'[GREJCIS]\d\d-az') - 6;
            el_start = regexp(txt(lim(1,1) : lim(1,2)),'[GREJCIS]\d\d-el') - 6;
            num_sat = numel(az_start);
            this.sat.az = [this.sat.az; zeros(n_epo, num_sat)];
            this.sat.el = [this.sat.el; zeros(n_epo, num_sat)];
            for s = 1 : num_sat
                az = sscanf(txt(bsxfun(@plus, repmat(az_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1)))', '%f');;
                el = sscanf(txt(bsxfun(@plus, repmat(el_start(s) + (0 : 15), size(lim, 1) - 1, 1), lim(2 : end, 1)))', '%f');;
                this.sat.az(id_int, s) = az(id_ext);
                this.sat.el(id_int, s) = el(id_ext);
            end
        end

    end
    
    % ==================================================================================================================================================
    %  PLOTTING FUNCTIONS
    % ==================================================================================================================================================
    
    methods (Access = public)
        function plotPositionENU(this, one_plot)
            % Plot East North Up coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionENU();
            if nargin == 1
                one_plot = false;
            end
            if size(this.xyz,1) > 1
                this.log.addMessage('Plotting positions');
                enu = zero2nan(this.xyz);
                xyz0 = this.getAPrioriPos;
                [enu0(:,1), enu0(:,2), enu0(:,3)] = cart2plan(xyz0(:,1), xyz0(:,2), xyz0(:,3));
                [enu(:,1), enu(:,2), enu(:,3)] = cart2plan(zero2nan(this.xyz(:,1)), zero2nan(this.xyz(:,2)), zero2nan(this.xyz(:,3)));
                figure;
                color_order = handle(gca).ColorOrder;
                t = this.time.getMatlabTime();
                if ~one_plot, subplot(3,1,1); end
                plot(t, 1e0 * (enu(:,1) - enu0(1)), '.-', 'MarkerSize', 5, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                ax(3) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('East [m]'); h.FontWeight = 'bold';
                grid on;
                h = title(sprintf('Receiver %s', this.name),'interpreter', 'none'); h.FontWeight = 'bold'; h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                if ~one_plot, subplot(3,1,2); end
                plot(t, 1e0 * (enu(:,2) - enu0(2)), '.-', 'MarkerSize', 5, 'LineWidth', 2, 'Color', color_order(2,:));
                ax(2) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('North [m]'); h.FontWeight = 'bold';
                grid on;
                if ~one_plot, subplot(3,1,3); end
                plot(t, 1e0 * (enu(:,3) - enu0(3)), '.-', 'MarkerSize', 5, 'LineWidth', 2, 'Color', color_order(3,:));
                ax(1) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('Up [m]'); h.FontWeight = 'bold';
                grid on;
                if one_plot
                    h = ylabel('ENU [m]'); h.FontWeight = 'bold';
                else
                    linkaxes(ax, 'x');
                end
            else
                this.log.addMessage('Plotting a single point static position is not yet supported');
            end
            grid on;
        end
        
        function plotPositionXYZ(this, one_plot)
            % Plot X Y Z coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionXYZ();
            if nargin == 1
                one_plot = false;
            end
            if size(this.xyz,1) > 1
                this.log.addMessage('Plotting positions');
                xyz0 = this.getAPrioriPos();
                
                figure;
                color_order = handle(gca).ColorOrder;
                t = this.time.getMatlabTime();
                if ~one_plot, subplot(3,1,1); end
                plot(t, 1e0 * (zero2nan(this.xyz(:,1)) - xyz0(1)), '.-', 'MarkerSize', 5, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                ax(3) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('X [m]'); h.FontWeight = 'bold';
                grid on;
                h = title(sprintf('Receiver %s', this.name),'interpreter', 'none'); h.FontWeight = 'bold'; h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                if ~one_plot, subplot(3,1,2); end
                plot(t, 1e0 * (zero2nan(this.xyz(:,2)) - xyz0(2)), '.-', 'MarkerSize', 5, 'LineWidth', 2, 'Color', color_order(2,:));
                ax(2) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('Y [m]'); h.FontWeight = 'bold';
                grid on;
                if ~one_plot, subplot(3,1,3); end
                plot(t, 1e0 * (zero2nan(this.xyz(:,3)) - xyz0(3)), '.-', 'MarkerSize', 5, 'LineWidth', 2, 'Color', color_order(3,:));
                ax(1) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); 
                grid on;
                if one_plot
                    h = ylabel('XYZ [m]'); h.FontWeight = 'bold';
                end
                linkaxes(ax, 'x');
            else
                this.log.addMessage('Plotting a single point static position is not yet supported');
            end            
        end
        
        function plotVsSynt(this)
            % Plots phases and pseudo-ranges aginst their synthesised values
            % SYNTAX: this.plotVsSynt
            
            % Phases
            [ph, ~, id_ph] = this.getPhases;
            sensor_ph = Core_Pre_Processing.diffAndPred(ph - this.getSyntPhases); sensor_ph = bsxfun(@minus, sensor_ph, median(sensor_ph, 2, 'omitnan'));
            figure; subplot(2,3,1); plot(sensor_ph); title('Phases observed vs synthesised');
            
            this.updateAzimuthElevation()
            id_ok = (~isnan(sensor_ph));
            az = this.sat.az(:,this.go_id(id_ph));
            el = this.sat.el(:,this.go_id(id_ph));
            %flag = flagExpand(abs(sensor1(id_ok)) > 0.2, 1);
            h1 = subplot(2,3,3);  polarScatter(serialize(az(id_ok))/180*pi,serialize(90-el(id_ok))/180*pi, 30, serialize(sensor_ph(id_ok)), 'filled');
            caxis([-1 1]); colormap(gat); setColorMapGat([-1 1], 0.8, [-0.15 0.15]); colorbar();
            subplot(2,3,2); scatter(serialize(az(id_ok)), serialize(el(id_ok)), 50, abs(serialize(sensor_ph(id_ok))) > 0.2, 'filled'); caxis([-1 1]);
            
            % Pseudo Ranges
            [pr, id_pr] = this.getPseudoRanges;
            sensor_pr = Core_Pre_Processing.diffAndPred(pr - this.getSyntPrObs); sensor_pr = bsxfun(@minus, sensor_pr, median(sensor_pr, 2, 'omitnan'));
            subplot(2,3,4); plot(sensor_pr); title('Pseudo-ranges observed vs synthesised');
            
            id_ok = (~isnan(sensor_pr));
            az = this.sat.az(:,this.go_id(id_pr));
            el = this.sat.el(:,this.go_id(id_pr));
            %flag = flagExpand(abs(sensor1(id_ok)) > 0.2, 1);
            h1 = subplot(2,3,6); polarScatter(serialize(az(id_ok))/180*pi,serialize(90-el(id_ok))/180*pi, 30, serialize(sensor_pr(id_ok)), 'filled');
            caxis([-20 20]); colorbar();
            subplot(2,3,5); scatter(serialize(az(id_ok)), serialize(el(id_ok)), 50, abs(serialize(sensor_pr(id_ok))) > 5, 'filled');
            caxis([-1 1]);
        end
        
        function plotDt(this)
            % Plot Clock error
            % SYNTAX: this.plotDt
            
            figure;
            t = this.time.getMatlabTime;
            plot(t, this.desync, '-k', 'LineWidth', 2);
            hold on;
            plot(t, this.dt_pr, ':', 'LineWidth', 2);
            plot(t, this.dt_ph, ':', 'LineWidth', 2);
            plot(t, this.dt_ip, '-', 'LineWidth', 2);
            plot(t, this.dt_ip + this.dt_pr, '-', 'LineWidth', 2);
            legend('desync time', 'dt pre-estimated from pseudo ranges', 'dt pre-estimated from phases', 'dt correction from LS on Code', 'dt estimated from pre-processing', 'Location', 'northeastoutside');
            xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('receiver clock error [s]'); h.FontWeight = 'bold';
            
            h = title(sprintf('dt - receiver %s', this.name),'interpreter', 'none'); h.FontWeight = 'bold'; h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
        end
        
        function plotSNR(this, sys_c)
            % Plot Signal to Noise Ration in a skyplot
            % SYNTAX: this.plotSNR(sys_c)
            
            % SNRs
            if nargin == 2
                [snr, snr_id] = this.getSNR(sys_c);
            else
                [snr, snr_id] = this.getSNR();
            end
            
            figure;
            this.updateAzimuthElevation()
            id_ok = (~isnan(snr));
            az = this.sat.az(:,this.go_id(snr_id));
            el = this.sat.el(:,this.go_id(snr_id));
            polarScatter(serialize(az(id_ok))/180*pi,serialize(90-el(id_ok))/180*pi, 45, serialize(snr(id_ok)), 'filled');
            colormap(jet);  cax = caxis(); caxis([min(cax(1), 10), max(cax(2), 55)]); setColorMap([10 55], 0.9); colorbar();
            h = title(sprintf('SNR - receiver %s', this.name),'interpreter', 'none'); h.FontWeight = 'bold'; h.Units = 'pixels'; h.Position(2) = h.Position(2) + 20; h.Units = 'data';
        end
        
        function plotDataAvailability(this, sys_c)
            % Plot all the satellite seen by the system
            % SYNTAX: this.plotDataAvailability(sys_c)
            
            if (nargin == 1)
                sys_c = 'GREJCIS';
            end
            figure;
            ss_ok = intersect(this.cc.sys_c, sys_c);
            for ss = ss_ok
                ss_id = find(this.cc.sys_c == ss);
                switch numel(ss_ok)
                    case 2
                        subplot(1,2, ss_id);
                    case 3
                        subplot(2,2, ss_id);
                    case 4
                        subplot(2,2, ss_id);
                    case 5
                        subplot(2,3, ss_id);
                    case 6
                        subplot(2,3, ss_id);
                    case 7
                        subplot(2,4, ss_id);
                end
                
                for prn = this.cc.prn(this.cc.system == ss)'
                    id_ok = find(any(this.obs((this.system == ss)' & this.prn == prn, :),1));
                    plot(id_ok, prn * ones(size(id_ok)), 's', 'Color', [0.8 0.8 0.8]);
                    hold on;
                    id_ok = find(any(this.obs((this.system == ss)' & this.prn == prn & this.obs_code(:,1) == 'C', :),1));
                    plot(id_ok, prn * ones(size(id_ok)), 'o', 'Color', [0.2 0.2 0.2]);
                    hold on;
                    id_ok = find(any(this.obs((this.system == ss)' & this.prn == prn & this.obs_code(:,1) == 'L', :),1));
                    plot(id_ok, prn * ones(size(id_ok)), '.');
                end
                prn_ss = unique(this.prn(this.system == ss));
                xlim([1 size(this.obs,2)]);
                ylim([min(prn_ss) - 1 max(prn_ss) + 1]);
                h = ylabel('PRN'); h.FontWeight = 'bold';
                ax = gca(); ax.YTick = prn_ss;
                grid on;
                h = xlabel('epoch'); h.FontWeight = 'bold';
                title(this.cc.SYS_EXT_NAME{this.cc.SYS_C == ss});
                dockAllFigures
            end
        end
               
        function plotCycleSlip(this)
            if ~isempty(this.cycle_slip_idx_ph)
                ph = this.getPhases();
                synt_ph = this.getSyntPhases();
                d_ph = ph - synt_ph;
                figure;
                plot(d_ph);
                ep = repmat([1: this.time.length]',1,size(d_ph,2));
                hold on
                scatter(ep(this.cycle_slip_idx_ph~=0),d_ph(this.cycle_slip_idx_ph~=0))
            end
        end
        
        function plotSmoothRef(this)
            ph = this.getPhases();
            synt_ph = this.getSyntPhases();
            d_ph = ph - synt_ph;
            d_ph(this.cycle_slip_idx_ph ~= 0) =nan;
            filt_length = 1000;
            idx = double(~isnan(d_ph));
            idx_s = flagShrink(idx,filt_length/2);
            idx_e = flagExpand(idx_s,filt_length/2);
            el_idx = xor(idx,idx_e);
            idx(el_idx) = 0;
            filt_arm = floor(filt_length/2);
            idx = flagShrink(idx,4);
            
            filt = (-filt_arm:filt_arm)/filt_length*5;
            filt = exp(-(filt).^2);
            filt = filt'./sum(filt);
            sats = 1 : size(ph,2);
            idx_f = zeros(size(idx));
            d_ph2= zeros(size(d_ph));
            for i = sats
                idx_f(:,i) = conv(idx(:,i),filt,'same');
                d_ph2(:,i) = d_ph(:,i) -mean(d_ph2(:,i),'omitnan');
            end
            idx1 = min(circshift(idx_f,filt_arm),circshift(idx_f,-filt_arm));
            w_ph = d_ph2 .* idx1;
            r_ph = nan(size(ph));
            
            for i = sats
                w_ph_t = w_ph(:,sats~=i);
                scale = sum(idx1(:,sats~=i),2,'omitnan');
                to_diff = sum(w_ph_t,2,'omitnan')./scale;
                x = (1:length(to_diff))';
                valid_x = ~isnan(to_diff);
                x = x(valid_x);
                xx = 0:300:length(valid_x);
                pp = spline(x,to_diff(valid_x),xx(8:end-7));
                pp = [repmat(pp(1),1,7) pp repmat(pp(end),1,7)];
                to_diff(valid_x) = to_diff(valid_x)- spline(xx,pp,x);%polyval(p,x);
                
                r_ph(:,i) = d_ph(:,i) - to_diff;
            end
            figure;
            plot(r_ph);
            ep = repmat([1: this.time.length]',1,size(r_ph,2));
            hold on
            scatter(ep(this.cycle_slip_idx_ph~=0),r_ph(this.cycle_slip_idx_ph~=0))
        end
        
        function plotAniZtdSlant(this, time_start, time_stop, show_map)
            clf;            
            t = this.time.getMatlabTime;
            
            sztd = this.getSlantZTD(900);
            if nargin >= 3
                if isa(time_start, 'GPS_Time')
                    time_start = find(this.time.getMatlabTime >= time_start.first.getMatlabTime(), 1, 'first');
                    time_stop = find(this.time.getMatlabTime <= time_stop.last.getMatlabTime(), 1, 'last');
                end
                time_start = max(1, time_start);
                time_stop = min(size(sztd,1), time_stop);
            else
                time_start = 1;
                time_stop = size(sztd,1);
            end
           
            if nargin < 4
                show_map = true;
            end
            win_size = (t(time_stop) - t(time_start)) * 86400;
            
            yl = (median(median(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan') + ([-6 6]) .* median(std(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan'));

            subplot(3,1,3);            
            plot(t, sztd,'.'); hold on;
            plot(t, this.ztd,'k', 'LineWidth', 4);
            ylim(yl);
            hl = line('XData', t(1) * [1 1],'YData', yl, 'LineWidth', 2);
            xlim(t(time_start) + [0 win_size-1] ./ 86400);
            setTimeTicks(4,'dd/mm/yyyy HH:MMPM');
            h = ylabel('ZTD [m]'); h.FontWeight = 'bold';
            grid on;
            
            % polar plot "true" Limits
            e_grid = [-1 : 0.1 : 1];
            n_grid = [-1 : 0.1 : 1];
            [ep, np] = meshgrid(e_grid, n_grid);
            fun = @(dist) exp(-((dist*1e5)/3e4).^2);
            
            ax_sky = subplot(3,1,1:2); i = time_start;
            az = (mod(this.sat.az(i,:) + 180, 360) -180) ./ 180 * pi; az(isnan(az) | isnan(sztd(i,:))) = 1e10;
            el = (90 - this.sat.el(i,:)) ./ 180 * pi; el(isnan(el) | isnan(sztd(i,:))) = 1e10;
            
            if show_map
                td = nan(size(ep));
                hm = imagesc(e_grid, n_grid, reshape(td(:), numel(n_grid), numel(e_grid))); hold on;
                hm.AlphaData = 0.5;   
                ax_sky.YDir = 'normal';
            end
            hs = polarScatter(az, el, 250, sztd(i,:), 'filled');
            xlim([-1 1]); ylim([-1 1]);
            caxis(yl); colormap(jet); colorbar;                        
            
            subplot(3,1,3); 
            for i = time_start + 1 : time_stop
                % Move scattered points
                az = (mod(this.sat.az(i,:) + 180, 360) -180) ./ 180 * pi; az(isnan(az) | isnan(sztd(i,:))) = 1e10;
                el = (90 - this.sat.el(i,:)) ./ 180 * pi; el(isnan(el) | isnan(sztd(i,:))) = 1e10;
                decl_n = el/(pi/2);
                x = sin(az) .* decl_n;
                y = cos(az) .* decl_n;
                
                id_ok = not(isnan(zero2nan(sztd(i,:))));
                if show_map
                    td = funInterp2(ep(:), np(:), x(1, id_ok)', y(1, id_ok)', sztd(i, id_ok)', fun);
                    hm.CData = reshape(td(:), numel(n_grid), numel(e_grid));
                end
                
                hs.XData = x;
                hs.YData = y;
                hs.CData = sztd(i,:);
                
                % Move time line
                hl.XData = t(i) * [1 1];
                if nargin > 4
                    xlim(t(i) + [-win_size/2 win_size/2] ./ 86400);
                end
                drawnow;
            end
            
        end
        
        function plotZtdSlant(this, time_start, time_stop, win_size)
            clf;
            t = this.time.getMatlabTime;
            
            sztd = this.getSlantZTD(900);
            if nargin >= 3
                if isa(time_start, 'GPS_Time')
                    time_start = find(this.time.getMatlabTime >= time_start.first.getMatlabTime(), 1, 'first');
                    time_stop = find(this.time.getMatlabTime <= time_stop.last.getMatlabTime(), 1, 'last');
                end
                time_start = max(1, time_start);
                time_stop = min(size(sztd,1), time_stop);
            else
                time_start = 1;
                time_stop = size(sztd,1);
            end
            
            if nargin < 4
                win_size = (t(time_stop) - t(time_start)) * 86400;
            end
            
            %yl = (median(median(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan') + ([-6 6]) .* median(std(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan'));

            plot(t, sztd,'.'); hold on;
            plot(t, this.ztd,'k', 'LineWidth', 4);
            %ylim(yl);
            xlim(t(time_start) + [0 win_size-1] ./ 86400);
            setTimeTicks(4,'dd/mm/yyyy HH:MMPM');
            h = ylabel('ZTD [m]'); h.FontWeight = 'bold';
            grid on;
            h = title(sprintf('Receiver %s ZTD', this.name),'interpreter', 'none'); h.FontWeight = 'bold'; h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
        end
    end
    
    % ==================================================================================================================================================
    %  PRIVATE FUNCTIONS
    % ==================================================================================================================================================
    
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
            all_sat(all_sat(:,1) == 32) = this.rinex_ss;
            all_sat(all_sat == 32) = '0'; % sscanf seems to misbehave with spaces
            gps_prn = unique(sscanf(all_sat(all_sat(:,1) == 'G', 2 : 3)', '%2d'));
            glo_prn = unique(sscanf(all_sat(all_sat(:,1) == 'R', 2 : 3)', '%2d'));
            gal_prn = unique(sscanf(all_sat(all_sat(:,1) == 'E', 2 : 3)', '%2d'));
            qzs_prn = unique(sscanf(all_sat(all_sat(:,1) == 'J', 2 : 3)', '%2d'));
            bds_prn = unique(sscanf(all_sat(all_sat(:,1) == 'C', 2 : 3)', '%2d'));
            irn_prn = unique(sscanf(all_sat(all_sat(:,1) == 'I', 2 : 3)', '%2d'));
            sbs_prn = unique(sscanf(all_sat(all_sat(:,1) == 'S', 2 : 3)', '%2d'));
            prn = struct('g', gps_prn', 'r', glo_prn', 'e', gal_prn', 'j', qzs_prn', 'c', bds_prn', 'i', irn_prn', 's', sbs_prn');
            
            % update the maximum number of rows to store
            n_obs = this.cc.gps.isActive * numel(prn.g) * numel(this.rin_obs_code.g) / 3 + ...
                this.cc.glo.isActive * numel(prn.r) * numel(this.rin_obs_code.r) / 3 + ...
                this.cc.gal.isActive * numel(prn.e) * numel(this.rin_obs_code.e) / 3 + ...
                this.cc.qzs.isActive * numel(prn.j) * numel(this.rin_obs_code.j) / 3 + ...
                this.cc.bds.isActive * numel(prn.c) * numel(this.rin_obs_code.c) / 3 + ...
                this.cc.irn.isActive * numel(prn.i) * numel(this.rin_obs_code.i) / 3 + ...
                this.cc.sbs.isActive * numel(prn.s) * numel(this.rin_obs_code.s) / 3;
            
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
            bad_epochs = [];
            for e = 1 : n_epo % for each epoch
                n_sat = this.n_spe(e);
                % get the list of satellites in view
                sat = serialize(txt(lim(t_line(e),1) + repmat((0 : ceil(this.n_spe(e) / 12) - 1)' * 69, 1, 36) + repmat(32:67, ceil(this.n_spe(e) / 12), 1))')';
                sat = sat(1:n_sat * 3);
                sat = reshape(sat, 3, n_sat)';
                sat(sat(:,1) == 32) = this.rinex_ss;
                sat(sat == 32) = '0';  % sscanf seems to misbehave with spaces
                prn_e = sscanf(serialize(sat(:,2:3)'), '%02d');
                if numel(prn_e) < this.n_spe(e)
                    bad_epochs = [bad_epochs; e];
                    cm = this.log.getColorMode();
                    this.log.setColorMode(false); % disable color mode for speed up
                    this.log.addWarning(sprintf('Problematic epoch found at %s\nInspect the files to detect what went wrong!\nSkipping and continue the parsing, no action taken%s', this.time.getEpoch(e).toString, char(32*ones(this.w_bar.bar_len,1))));
                    this.log.setColorMode(cm);
                else
                    for s = 1 : size(sat, 1)
                        % line to fill with the current observation line
                        obs_line = find((this.prn == prn_e(s)) & this.system' == sat(s, 1));
                        % if is empty I'm reading a constellation that is not active in the constellation collector
                        %   -> discard the unwanted satellite
                        if ~isempty(obs_line)
                            line_start = t_line(e) + ceil(n_sat / 12) + (s-1) * n_lps;
                            line = mask(1 : n_ops * 16);
                            for i = 0 : n_lps - 1
                                try
                                    line(id_line(1:lim(line_start + i, 3),i+1)) = txt(lim(line_start + i, 1) : lim(line_start + i, 2)-1);
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
                    end
                end
                this.w_bar.go(e);
            end
            obs(:, bad_epochs) = [];
            this.time.remEpoch(bad_epochs);
            this.log.newLine();
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
            n_obs = this.cc.gps.isActive * numel(prn.g) * numel(this.rin_obs_code.g) / 3 + ...
                this.cc.glo.isActive * numel(prn.r) * numel(this.rin_obs_code.r) / 3 + ...
                this.cc.gal.isActive * numel(prn.e) * numel(this.rin_obs_code.e) / 3 + ...
                this.cc.qzs.isActive * numel(prn.j) * numel(this.rin_obs_code.j) / 3 + ...
                this.cc.bds.isActive * numel(prn.c) * numel(this.rin_obs_code.c) / 3 + ...
                this.cc.irn.isActive * numel(prn.i) * numel(this.rin_obs_code.i) / 3 + ...
                this.cc.sbs.isActive * numel(prn.s) * numel(this.rin_obs_code.s) / 3;
            
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
                ss = this.cc.getSys(sys - 32);
                [~, f_id] = ismember(f_id, ss.CODE_RIN3_2BAND);
                
                ismember(this.system, this.cc.SYS_C);
                this.f_id = [this.f_id; f_id];
                
                if sys == 'r'
                    wl = ss.L_VEC((max(1, f_id) - 1) * size(ss.L_VEC, 1) + ss.PRN2IDCH(min(prn_ss, ss.N_SAT))');
                    wl(prn_ss > ss.N_SAT) = NaN;
                    wl(f_id == 0) = NaN;
                else
                    wl = ss.L_VEC(max(1, f_id))';
                    wl(f_id == 0) = NaN;
                end
                if sum(f_id == 0)
                    [~, id] = unique(double(obs_code(f_id == 0, :)) * [1 10 100]');
                    this.log.addWarning(sprintf('These codes for the %s are not recognized, ignoring data: %s', ss.SYS_EXT_NAME, sprintf('%c%c%c ', obs_code(id, :)')));
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
                    if ~isempty(obs_line)
                        line = txt(lim(t_line(e) + s, 1) + 3 : lim(t_line(e) + s, 2));
                        ck = line == ' '; line(ck) = mask(ck); % fill empty fields -> otherwise textscan ignore the empty fields
                        % try with sscanf
                        line = line(data_pos(1 : numel(line)));
                        data = sscanf(reshape(line, 14, numel(line) / 14), '%f');
                        obs(obs_line(1:size(data,1)), e) = data;
                    end
                    % alternative approach with textscan
                    %data = textscan(line, '%14.3f%1d%1d');
                    %obs(obs_line(1:numel(data{1})), e) = data{1};
                end
                this.w_bar.go(e);
            end
            this.log.newLine();
            this.obs = obs;
            
        end
        
        function synt_pr_obs = computeSyntCurObs(this, phase, sys_c)
            % DESCRIPTION: get syntetic observation for code or phase
            obs_type = {'code', 'phase'};
            this.log.addMessage(this.log.indent(sprintf('Synthesising %s observations', obs_type{phase + 1}),6));
            idx_obs = [];
            if nargin < 3
                sys_c = this.cc.sys_c;
            end
            
            for s = sys_c
                if phase
                    idx_obs = [idx_obs; this.getObsIdx('L', s)];
                else
                    idx_obs = [idx_obs; this.getObsIdx('C', s)];
                end
            end
            
            synt_pr_obs = zeros(length(idx_obs), size(this.obs,2));
            sys = this.system(idx_obs);
            prn = this.prn(idx_obs);
            sat = this.cc.getIndex(sys,prn);
            u_sat = unique(sat);
            this.updateAllAvailIndex();
            this.updateErrIono();
            this.updateErrTropo();
            for i = u_sat'
                sat_idx = find(sat == i);
                
                this.updateTOT(colFirstNonZero(this.obs(sat_idx,:)),i);
                range = this.computeSyntObs('I', i);
                for j = sat_idx'
                    o = idx_obs(j);
                    o = find(idx_obs == o);
                    c_obs_idx = idx_obs(j); % index of the observation we are currently processing
                    ep_idx = this.obs(c_obs_idx,:) ~= 0;
                    freq = this.cc.getBand(sys(j), this.obs_code(c_obs_idx,2));
                    if sys(j) == 'G' && freq == 1
                        iono_factor = 1;
                    else
                        wl_ref = this.cc.getGPS.F_VEC(1);
                        wl = this.cc.getSys(sys(j)).F_VEC(freq);
                        iono_factor= wl_ref ^ 2/ wl ^ 2;
                    end
                    synt_pr_obs(o, ep_idx) = range(ep_idx) + iono_factor * this.sat.err_iono(ep_idx,i)';
                end
                
            end
        end
        
        function DtSat(this,flag)
            % DESCRIPTION. apply clock satellite corrections for code and phase
            % IMPORTANT: if no clock is present delete the observation
            
            if isempty(this.active_ids)
                this.active_ids = false(size(this.obs,1),1);
            end
            
            for i = 1 : this.cc.getNumSat()
                prn = this.cc.prn(i);
                sys = this.cc.system(i);
                sat_idx = this.prn == prn & (this.system == sys)' & (this.obs_code(:,1) == 'C' | this.obs_code(:,1) == 'L');
                ep_idx = logical(sum(this.obs(sat_idx,:) ~= 0));
                this.updateAvailIndex(ep_idx,i);
                dts_range = ( this.getDtS(i) + this.getRelClkCorr(i) ) * goGNSS.V_LIGHT;
                for o = find(sat_idx)'
                    obs_idx_l = this.obs(o,:) ~= 0;
                    obs_idx = find(obs_idx_l);
                    dts_idx = obs_idx_l(ep_idx);
                    if this.obs_code(o,1) == 'C'
                        this.obs(o, obs_idx_l) = this.obs(o,obs_idx_l) + sign(flag) * dts_range(dts_idx)';
                    else
                        this.obs(o, obs_idx_l) = this.obs(o,obs_idx_l) + sign(flag) * dts_range(dts_idx)'./this.wl(o);
                    end
                    dts_range_2 = dts_range(dts_idx);
                    nan_idx = obs_idx(isnan(dts_range_2));
                    this.obs(o, nan_idx) = 0;
                end
            end
        end
        
        function GroupDelay(this, sgn)
            % DESCRIPTION. apply group delay corrections for code and phase
            % measurement when a value if provided from an external source
            % (Navigational file  or DCB file)
            for i = 1 : size(this.sat.cs.group_delays, 2)
                sys  = this.sat.cs.group_delays_flags(i,1);
                code = this.sat.cs.group_delays_flags(i,2:4);
                f_num = str2double(code(2));
                idx = this.getObsIdx(code, sys);
                if sum(this.sat.cs.group_delays(:,i)) ~= 0
                    if ~isempty(idx)
                        for s = 1 : size(this.sat.cs.group_delays,1)
                            sat_idx = idx((this.prn(idx) == s));
                            full_ep_idx = not(abs(this.obs(sat_idx,:)) < 0.1);
                            if this.sat.cs.group_delays(s,i) ~= 0
                                this.obs(sat_idx,full_ep_idx) = this.obs(sat_idx,full_ep_idx) + sign(sgn) * this.sat.cs.group_delays(s,i);
                            elseif ~this.cc.isRefFrequency(sys, f_num)
                                this.active_ids(idx) = false;
                                idx = this.getObsIdx(['C' code(2:end)], sys);
                                this.active_ids(sat_idx) = sgn < 0;
                            end
                        end
                    end
                else
                    % mark as bad obs frequencies that are nor reference frequencies or that have no correction
                    if ~this.cc.isRefFrequency(sys, f_num)
                        idx = this.getObsIdx(['C' code(2:end)], sys);
                        this.active_ids(idx) = sgn < 0;
                    end
                end
            end
        end
        
        function codeStaticPositionig(this, obs, prn, sys, flag, opt)
            % INPUT:
            %   opt: structure with options of the LS adjustement
            %        .coord_corr: stop if coordinate correction goes under the paramter
            %        .max_it:     maximum number of iterations
            %        .no_pos:     compute dt only
            % DESCRITION compute the postion of the receiver based on code
            % measurements
            
            if nargin < 5
                opt = struct('coord_corr', 0.1, ...
                    'max_it',  10, ...
                    'no_pos', false);
            end
            if ~isfield(opt,'no_pos')
                opt.no_pos = false;
            end
            if ~isfield(opt,'rid_ep')
                opt.rid_ep  = false;
            end
            
            n_epochs         = this.time.getLen;
            n_valid_epochs   = sum(any(obs,1));
            
            % get the type of observations to be used for the positioning
            code_bias_flag   = cellstr([sys flag]);
            u_code_bias_flag = unique(code_bias_flag);
            
            % initialize dt
            n_obs_ch         = zeros(size(u_code_bias_flag));
            n_ep_ch          = zeros(size(u_code_bias_flag));
            ch_idx_ep        = zeros(length(u_code_bias_flag),n_epochs);
            for i = 1 : length(n_obs_ch)
                ch_idx_sat = sum([sys flag] == repmat( sprintf('%-8s',u_code_bias_flag{i}), length(sys),1),2) == 8;
                n_obs_ch(i) = sum((ch_idx_sat).*sum(obs > 0, 2)); % find number of observations per channel
                ch_idx_ep(i,:) = sum(obs(ch_idx_sat,:),1) > 0;
                n_ep_ch(i) = sum(ch_idx_ep(i,:)); % find number of observations per channel
            end
            % sort the channel variables by the number of observables
            [~, b] = sort(n_obs_ch,'descend');
            u_code_bias_flag = u_code_bias_flag(b);
            n_ep_ch = n_ep_ch(b,:);
            ch_idx_ep = ch_idx_ep(b,:);
            
            % compute a column with an integer that indicate which
            % code_bias to estimate for each obs
            code_bias_ord = zeros(size(code_bias_flag,1),1);
            for i = 1 :length(u_code_bias_flag)
                ch_idx_sat = sum([sys flag] == repmat( sprintf('%-8s',u_code_bias_flag{i}),length(sys),1),2) == 8;
                code_bias_ord(ch_idx_sat) = i;
            end
            
            % get the satellite index for all obs
            sat_ids = zeros(size(prn));
            for s = 1:length(sat_ids)
                sat_ids(s) = this.cc.getIndex(sys(s),prn(s));
            end
            
            this.dt = zeros(this.time.length,1);
            this.rid = zeros(1, length(u_code_bias_flag)-1);
            this.flag_rid = u_code_bias_flag;
            
            n_tot_obs = sum(sum(obs>0));
            
            x = [999 999 999];
            
            % ls_solver.A = sparse(n_tot_obs,3+n_valid_epochs+sum(n_ep_ch(2:end))); % version with reference clock
            
            n_it = 0;
            while max(abs(x(1:3))) > opt.coord_corr && n_it < opt.max_it
                n_it = n_it + 1;
                % fill the a matrix
                XS_norm = zeros(this.cc.getNumSat(), 3, n_epochs);
                dist = zeros(n_epochs, this.cc.getNumSat());
                
                for i = 1 : this.cc.getNumSat()
                    c_sys = this.cc.system(i);
                    c_prn = this.cc.prn(i);
                    idx_sat = sys == c_sys & prn == c_prn;
                    idx_sat_i = find(idx_sat);
                    if sum(idx_sat) > 0 % if we have an obs for the satellite
                        c_obs = obs(idx_sat,:);
                        
                        c_l_obs = colFirstNonZero(c_obs); % all best obs on one line
                        idx_obs = c_l_obs > 0; % epoch with obs from the satellite
                        
                        % update time of flight times
                        this.updateAvailIndex(c_l_obs, i);
                        this.updateTOT(c_l_obs, i); % update time of travel
                        freq = flag(idx_sat_i(1), 7);
                        if freq == ' '
                            freq = flag(idx_sat_i(1), 2);
                        end
                        [dist(:,i), XS] = this.computeSyntObs(freq,i); %%% consider multiple combinations (different iono corrections) on the same satellite, not handdled yet
                        
                        XS_norm(i,:,idx_obs) = rowNormalize(XS)';
                    end
                end
                
                % preallocate normal matrix and
                n_par = 3+n_valid_epochs+length(u_code_bias_flag)-1;
                % N = spalloc(n_par,n_par,(3+length(u_code_bias_flag))*n_valid_epochs*2);
                
                num_dcb = length(u_code_bias_flag) -1;
                num_cp = 3 + num_dcb; % num constant parameter (pos + dcb)
                
                %  Normal matrix components
                %
                %     N_cp   |       N_col'
                %     -------+---------------
                %     N_col  |   N
                %            |     _
                %            |       d
                %            |         i
                %            |           a
                %            |             g
                %
                
                N_cp = zeros(num_cp);
                N_col = zeros(n_valid_epochs, num_cp); % N coord and dcb
                N_diag = zeros(n_valid_epochs, 1); % N clocks
                B = zeros(n_par,1);
                
                cur_val_ep = 0;
                %                     n_type_obs = length(obs_ep_idx_l); % number of obervation types
                %                     v_type_obs = (1 : n_type_obs)';
                v_dcb = 1: num_dcb;
                param_idx_c = [1 2 3 3+v_dcb]; % non time dependent param
                n_c_param = length(param_idx_c);
                v_c_param = 1 : n_c_param;
                % fill parts of normal matrix
                for e = 1 : n_epochs
                    obs_ep_idx_l = obs(:,e) > 0; % logical index
                    
                    if sum(obs_ep_idx_l) > 0
                        cur_val_ep = cur_val_ep +1; % current valid epoch
                        
                        num_obs_ep = sum(obs_ep_idx_l);
                        ch_obs = code_bias_ord(obs_ep_idx_l);
                        
                        A_dcb = zeros(num_obs_ep, num_dcb);
                        for i = v_dcb
                            A_dcb(ch_obs == i+1,i) = 1;
                        end
                        % construc design matrix and (y0-b) for the current epoch
                        A_ep = [-XS_norm(sat_ids(obs_ep_idx_l), : , e) A_dcb ones(num_obs_ep, 1) ];
                        y_ep = obs(obs_ep_idx_l, e) - dist(e, sat_ids(obs_ep_idx_l))';
                        % construct nomr matrix and A'*y for the
                        % current epoch
                        N_ep = A_ep'*A_ep;
                        B_ep = A_ep' * y_ep;
                        % fill the N and B matrix
                        
                        param_idx = [param_idx_c 3+num_dcb+cur_val_ep ];
                        B(param_idx) = B(param_idx) + B_ep;
                        N_diag(cur_val_ep) = N_diag(cur_val_ep) + N_ep(end,end); % sum clock N
                        for p = v_c_param
                            N_cp(param_idx_c(p),param_idx_c) = N_cp(param_idx_c(p),param_idx_c) + N_ep(p,param_idx_c);
                        end
                        N_col(cur_val_ep,param_idx_c) = N_ep(n_c_param+1,1:n_c_param);
                    end
                end
                
                if opt.no_pos
                    %[x, res] = ls_solver.solve([4:size(ls_solver.A,2)]);
                    num_param = num_dcb + n_valid_epochs;
                    N = spdiags([zeros(num_dcb,1); N_diag], 0, num_param, num_param);
                    N(1:num_dcb, 1:num_dcb) = N_cp(4:end, 4:end);
                    N((num_dcb+1):num_param, 1:num_dcb) = N_col(:, 4:end);
                    N(1:num_dcb, (num_dcb+1):num_param) = N_col(:, 4:end)';
                    B = B(4:end);
                    x = N\B;
                    x = [zeros(3,1) ; x];
                else
                    num_param = num_cp + n_valid_epochs;
                    N = spdiags([zeros(num_cp,1); N_diag], 0, num_param, num_param);
                    N(1:num_cp, 1:num_cp) = N_cp;
                    N((num_cp+1):num_param, 1:num_cp) = N_col;
                    N(1:num_cp, (num_cp+1):num_param) = N_col';
                    x = N\B;
                    %[x, res] = ls_solver.solve();
                end
                this.xyz = this.xyz + x(1:3)';
                
                if opt.rid_ep
                    for i = 1:length(u_code_bias_flag)
                        this.dt(ch_idx_ep(i,:) > 0,i) = x(((sum(n_ep_ch(1:i-1))) : (sum(n_ep_ch(1:i)) -1 ) ) + 4) / Go_State.V_LIGHT;
                    end
                else
                    this.dt(sum(obs,1) > 0,1,1) = x((4+num_dcb):end) / Go_State.V_LIGHT;
                    this.rid = x(4: (3+num_dcb)) / Go_State.V_LIGHT;
                end
            end
        end
        
        
        function [obs, sys, prn, flag] = removeUndCutOff(this, obs, sys, prn, flag, cut_off)
            % DESCRIPTION: remove obs under cut off
            for i = 1 : length(prn)
                sat = this.cc.getIndex(sys(i), prn(i)); % get go_id
                
                idx_obs = obs(i,:) ~= 0;
                this.updateAvailIndex(idx_obs, sat);
                XS = this.getXSTxRot(sat);
                if size(this.xyz,1) == 1
                    [~ , el] = this.computeAzimuthElevationXS(XS);
                else
                    [~ , el] = this.computeAzimuthElevationXS(XS,this.xyz(idx_obs,:));
                end
                
                idx_obs_f = find(idx_obs);
                el_idx = el < cut_off;
                idx_obs_f = idx_obs_f( el_idx );
                obs(i,idx_obs_f) = 0;
            end
            %%% remove possibly generated empty lines
            empty_idx = sum(obs > 0, 2) == 0;
            obs(empty_idx,:) = [];
            sys(empty_idx,:) = [];
            prn(empty_idx,:) = [];
            flag(empty_idx,:) = [];
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
        
        function [p_time] = getSyncTime(rec, obs_type, p_rate)
            % Get the common time among alle the used receivers and the target(s)
            %
            % SYNTAX: 
            %   p_time = Receiver.getSyncTime(rec, obs_type, <p_rate>);
            %
            % EXAMPLE:
            %   p_time = Receiver.getSyncTime(rec, state.obs_type, state.getProcessingRate());

            if nargin < 3
                p_rate = 1e-6;
            end
            % Do the target(s) as last
            [~, id] = sort(obs_type, 'descend');
            
            % prepare reference time
            % processing time will start with the receiver with the last first epoch
            %          and it will stop  with the receiver with the first last epoch
            p_time_zero = round(rec(1).time.first.getMatlabTime() * 24)/24; % get the reference time
            p_time_start = rec(1).time.first.getRefTime(p_time_zero);
            p_time_stop = rec(1).time.last.getRefTime(p_time_zero);
            p_rate = lcm(round(p_rate * 1e6), round(rec(1).time.getRate * 1e6)) * 1e-6;
            
            p_time = GPS_Time(); % empty initialization
            
            i = 0;
            for r = id
                if obs_type(r) > 0 % if it's not a target
                    p_time_start = max(p_time_start,  round(rec(r).time.first.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
                    p_time_stop = min(p_time_stop,  round(rec(r).time.last.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
                    p_rate = lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)) * 1e-6;
                else
                    % It's a target
                    
                    % recompute the parameters for the ref_time estimation
                    % not that in principle I can have up to num_trg_rec ref_time
                    % in case of multiple targets the reference times should be independent
                    % so here I keep the temporary rt0 rt1 r_rate var
                    % instead of ref_time_start, ref_time_stop, ref_rate
                    pt0 = max(p_time_start, round(rec(r).time.first.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
                    pt1 = min(p_time_stop, round(rec(r).time.last.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
                    pr = lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)) * 1e-6;
                    pt0 = ceil(pt0 / pr) * pr;
                    pt1 = floor(pt1 / pr) * pr;
                    
                    % return one p_time for each target
                    i = i + 1;
                    p_time(i) = GPS_Time(p_time_zero, (pt0 : pr : pt1)); %#ok<SAGROW>
                    p_time(i).toUnixTime();
                end
            end           
        end
    end

end
