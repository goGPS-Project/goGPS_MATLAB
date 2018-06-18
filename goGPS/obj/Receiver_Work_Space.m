%   CLASS Receiver_Work_Space
% =========================================================================
%
% 
%   Class to store receiver working paramters
%
% EXAMPLE
%   trg = Receiver_Work_Space();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Receiver

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea, Giulio Tagliaferro ...
%  Contributors:
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
%--------------------------------------------------------------------------
classdef Receiver_Work_Space < Receiver_Commons
    %% CONSTANTS
    properties (Constant)
        S02_IP_THR = 1e3;
    end
    
     % ==================================================================================================================================================
    %% PROPERTIES RECEIVER
    
    properties (SetAccess = public, GetAccess = public)
        
        file            % list of rinex file loaded
        rinex_file_name % list of RINEX file name
        rin_type       % rinex version format 
        
        rinex_ss       % flag containing the satellite system of the rinex file, G: GPS, R: GLONASS, E: Galileo, J: QZSS, C: BDS, I: IRNSS, S: SBAS payload, M: Mixed
        

        n_sat = 0;     % number of satellites
        n_freq = 0;    % number of stored frequencies
        n_spe = [];    % number of observations per epoch
        
        rec_settings;  % receiver settings

    end
     % ==================================================================================================================================================
     % ==================================================================================================================================================
    %% PROPERTIES OBSERVATIONS
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
        
        
        
        % CORRECTIONS ------------------------------
        
        rin_obs_code   % list of types per constellation
        ph_shift       % phase shift as read from RINEX files
        
        ocean_load_disp = [];        % ocean loading displacemnet for the station , -1 if not found
        clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true
        pcv = []                     % phase center corrections for the receiver
        
        pp_status = false;      % flag is pre-proccesed / not pre-processed
        iono_status           = 0; % flag to indicate if measurements ahave been corrected by an external ionpheric model
        group_delay_status    = 0; % flag to indicate if code measurement have been corrected using group delays                          (0: not corrected , 1: corrected)
        dts_delay_status      = 0; % flag to indicate if code and phase measurement have been corrected for the clock of the satellite    (0: not corrected , 1: corrected)
        sh_delay_status       = 0; % flag to indicate if code and phase measurement have been corrected for shapiro delay                 (0: not corrected , 1: corrected)
        pcv_delay_status      = 0; % flag to indicate if code and phase measurement have been corrected for pcv variations                (0: not corrected , 1: corrected)
        ol_delay_status       = 0; % flag to indicate if code and phase measurement have been corrected for ocean loading                 (0: not corrected , 1: corrected)
        pt_delay_status       = 0; % flag to indicate if code and phase measurement have been corrected for pole tides                    (0: not corrected , 1: corrected)
        pw_delay_status       = 0; % flag to indicate if code and phase measurement have been corrected for phase wind up                 (0: not corrected , 1: corrected)
        et_delay_status       = 0; % flag to indicate if code and phase measurement have been corrected for solid earth tide              (0: not corrected , 1: corrected)
        hoi_delay_status      = 0; % flag to indicate if code and phase measurement have been corrected for high order ionospheric effect (0: not corrected , 1: corrected)
        atm_load_delay_status = 0; % flag to indicate if code and phase measurement have been corrected for atmospheric loading           (0: not corrected , 1: corrected)
        
        % FLAGS ------------------------------
        
        ph_idx             % idx of outlier and cycle slip in obs
        
        
        if_amb;            % temporary varibale to test PPP ambiguity fixing
        est_slant;         % temporary varibale to test PPP ambiguity fixing
    end
     % ==================================================================================================================================================
    %% PROPERTIES POSITION
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        xyz_approx     % approximate position of the receiver (XYZ geocentric)
    end
     % ==================================================================================================================================================
    %% PROPERTIES CELESTIAL INFORMATIONS
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
    sat = struct( ...
            'avail_index',      [], ...    % boolean [n_epoch x n_sat] availability of satellites
            'o_out_ph',         [], ...    % outlier    processing to be saved in result
            'o_cs_ph',          [], ...    % cycle slip processing to be saved in result
            'outlier_idx_ph',   [], ...    % logical index of outliers
            'cycle_slip_idx_ph',[], ...    % logical index of cycle slips
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
            'res',              [], ...    % residual per staellite
            'slant_td',         [], ...    % slant total delay (except ionosphere delay)
            'amb_idx',          [], ...    % temporary variable to test PPP ambiguity fixing
            'amb_mat',          [], ...
            'amb',              [] ...
            )
    end
    % ==================================================================================================================================================
    %% PROPERTIES TIME
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        id_sync        % id of the epochs used for computing the solution
        rate           % observations rate;
        dt_ph          % clock error for the phases generated by the de-sync process
        dt_pr          % clock error for the pseudo-ranges generated by the de-sync process
    end
    
    % ==================================================================================================================================================
    %% PROPERTIES TROPO
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        meteo_data % meteo data object
    end
     % ==================================================================================================================================================
    %% METHODS INIT - CLEAN - RESET - REM -IMPORT
    % ==================================================================================================================================================
    
    methods
        function this = Receiver_Work_Space(cc, parent)
            % SYNTAX  this = Receiver(<cc>)
            
            this.parent = parent;
            if nargin >= 1 && ~isempty(cc)
                this.cc = cc;
            else
                this.cc = this.state.cc;
            end
            this.rec_settings = Receiver_Settings();
            this.init();
        end
        
        function init(this)
            this.log = Logger.getInstance();
            this.state = Global_Configuration.getCurrentSettings();
            this.rf = Core_Reference_Frame.getInstance();
            this.w_bar = Go_Wait_Bar.getInstance();
            this.reset();
        end
        
        function reset(this)
            this.reset@Receiver_Commons();
            this.resetObs();
            this.initR2S(); % restore handle to core_sky singleton
            this.id_sync = [];
            
            this.n_sat = [];
            
            this.pp_status = false;
            this.group_delay_status = 0; % flag to indicate if code measurement have been corrected using group delays                          (0: not corrected , 1: corrected)
            this.dts_delay_status   = 0; % flag to indicate if code and phase measurement have been corrected for the clock of the satellite    (0: not corrected , 1: corrected)
            this.sh_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for shapiro delay                 (0: not corrected , 1: corrected)
            this.pcv_delay_status   = 0; % flag to indicate if code and phase measurement have been corrected for pcv variations                (0: not corrected , 1: corrected)
            this.ol_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for ocean loading                 (0: not corrected , 1: corrected)
            this.pt_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for pole tides                    (0: not corrected , 1: corrected)
            this.pw_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for phase wind up                 (0: not corrected , 1: corrected)
            this.et_delay_status    = 0; % flag to indicate if code and phase measurement have been corrected for solid earth tide              (0: not corrected , 1: corrected)
            this.hoi_delay_status   = 0; % flag to indicate if code and phase measurement have been corrected for high order ionospheric effect          (0: not corrected , 1: corrected)
        end
        
        function resetObs(this)
            % Reset the content of the receiver obj
            % SYNTAX
            %   this.initObs;
            
            this.file = [];             % file rinex object
            this.rin_type = 0;          % rinex version format
            
            this.parent.ant          = 0;       % antenna number
            this.parent.ant_type     = '';      % antenna type
            this.parent.ant_delta_h  = 0;       % antenna height from the ground [m]
            this.parent.ant_delta_en = [0 0];   % antenna east/north offset from the ground [m]
            
            this.rin_obs_code = '';       % list of types per constellation
            this.ph_shift     = [];
            
            this.xyz          = [0 0 0];  % approximate position of the receiver (XYZ geocentric)
            
            this.parent.static       = true;     % the receivers are considered static unless specified
            
            this.n_sat = 0;               % number of satellites
            this.n_freq = 0;              % number of stored frequencies
            this.n_spe = [];              % number of sat per epoch
            
            this.rate = 0;                % observations rate;
            
            this.desync = 0;              % clock offset of the receiver
            this.dt_ph = 0;               % clock offset of the receiver
            this.dt_pr = 0;               % clock offset of the receiver
            this.dt_ip = 0;               % clock offset of the receiver
            this.dt = 0;                  % clock offset of the receiver
            this.flag_rid = 0;            % clock offset of the receiver
            
            this.active_ids = [];         % rows of active satellites
            this.wl         = [];         % wave-lenght of each row of row_id
            this.f_id       = [];         % frequency number e.g. L1 -> 1,  L2 ->2, E1 -> 1, E5b -> 3 ...
            this.prn        = [];         % pseudo-range number of the satellite
            this.go_id      = [];         % internal id for a certain satellite
            this.system     = '';         % char id of the satellite system corresponding to the row_id
            
            this.obs_code   = [];         % obs code for each line of the data matrix obs
            this.obs        = [];         % huge observation matrix with all the observables for all the systems / frequencies / ecc ...
            
            this.ph_idx = [];
            this.if_amb = [];
            this.est_slant = [];
            this.synt_ph = [];
            
            this.clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true
        end
        
        function resetWorkSpace(this)
            this.file = File_Rinex();
            this.rinex_file_name = '';
            this.rin_type = 0;
            this.time = GPS_Time();
            this.id_sync = [];
            this.n_sat = [];
            
            this.rin_obs_code = '';       % list of types per constellation
            this.ph_shift     = [];
            
            this.xyz          = [0 0 0];  % approximate position of the receiver (XYZ geocentric)
            this.enu = [];
            this.lat = [];
            this.lon = [];
            this.h_ellips = [];
            this.h_ortho = [];
            
            this.n_sat = 0;               % number of satellites
            this.n_freq = 0;              % number of stored frequencies
            this.n_spe = [];              % number of sat per epoch
            
            this.rate = 0;                % observations rate;
            
            this.desync = 0;              % clock offset of the receiver
            this.dt_ph = 0;               % clock offset of the receiver
            this.dt_pr = 0;               % clock offset of the receiver
            this.dt_ip = 0;               % clock offset of the receiver
            this.dt = 0;                  % clock offset of the receiver
            this.flag_rid = 0;            % clock offset of the receiver
            
            this.active_ids = [];         % rows of active satellites
            this.wl         = [];         % wave-lenght of each row of row_id
            this.f_id       = [];         % frequency number e.g. L1 -> 1,  L2 ->2, E1 -> 1, E5b -> 3 ...
            this.prn        = [];         % pseudo-range number of the satellite
            this.go_id      = [];         % internal id for a certain satellite
            this.system     = '';         % char id of the satellite system corresponding to the row_id
            
            this.obs_code   = [];         % obs code for each line of the data matrix obs
            this.obs        = [];         % huge observation matrix with all the observables for all the systems / frequencies / ecc ...
            
            this.ph_idx = [];
            this.if_amb = [];
            this.est_slant = [];
            this.synt_ph = [];
            
            this.clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true
            
            this.pp_status = false;
            this.iono_status = false;
            this.group_delay_status = false;
            this.dts_delay_status = false;
            this.sh_delay_status = false;
            this.pcv_delay_status = false;
            this.ol_delay_status = false;
            this.pt_delay_status = false;
            this.pw_delay_status = false;
            this.et_delay_status = false;
            this.hoi_delay_status = false;
            this.atm_load_delay_status = false;
            
            this.s0_ip =  [];
            this.s0 =  [];
            this.hdop =  [];
            this.vdop =  [];
            this.tdop =  [];
            
            this.a_fix = [];
            this.s_rate = [];
            
            this.xyz = [];
            
            this.apr_zhd  = [];
            this.zwd  = [];
            this.apr_zwd  = [];
            this.ztd  = [];
            this.pwv  = [];
            
            this.tgn = [];
            this.tge = [];
            
            this.sat.avail_index       = [];
            this.sat.outlier_idx_ph    = [];
            this.sat.cycle_slip_idx_ph = [];
            this.sat.err_tropo         = [];
            this.sat.err_iono          = [];
            this.sat.solid_earth_corr  = [];
            this.sat.dtS               = [];
            this.sat.rel_clk_corr      = [];
            this.sat.tot               = [];
            this.sat.az                = [];
            this.sat.el                = [];
            this.sat.XS_tx             = [];
            this.sat.crx               = [];
            this.sat.res               = [];
            this.sat.slant_td          = [];
            this.sat.o_cs_ph           = [];
            this.sat.o_out_ph          = [];
            this.sat.amb_mat           = [];
            this.sat.amb_idx           = [];
            this.sat.amb               = [];
        end
        
        function load(this)
            % Load the rinex file linked to this receiver
            %
            % SYNTAX
            %   this.load();
            
            if ~isempty(this.rinex_file_name)
                this.importRinex(this.rinex_file_name);
                this.importAntModel();
            end
            rf = Core_Reference_Frame.getInstance();
            coo = rf.getCoo(this.parent.getMarkerName4Ch, this.getCentralTime);
            if ~isempty(coo)
                this.xyz = coo;
            end
            
        end
        
        function appendRinex(this, rinex_file_name)
            rec = Receiver(this.cc, rinex_file_name, this.getDynamicMode);
            rec.load();
            
            % for the moment I merge receivers but then they need to run preProcessing again as a hole
            this.remAllCorrections()
            % merge the two rinex
            
        end
        
        function injestReceiver(this, rec)
            % assumption
            
        end
        
        function initTropo(this)
            % initialize tropo variables
            n_epoch = this.time.length;
            this.ztd = nan(n_epoch, 1);
            this.apr_zhd = nan(n_epoch, 1);
            this.zwd = nan(n_epoch, 1);
            this.apr_zwd = nan(n_epoch, 1);
            this.pwv = nan(n_epoch, 1);
            this.tge = nan(n_epoch, 1);
            this.tgn = nan(n_epoch, 1);
        end
               
        function initR2S(this)
            % initialize satellite related parameters
            % SYNTAX this.initR2S();
            
            this.sat.cs           = Core_Sky.getInstance();
            %this.sat.avail_index  = false(this.length, this.cc.getMaxNumSat);
            %  this.sat.XS_tx     = NaN(n_epoch, n_pr); % --> consider what to initialize
        end
        
        function subSample(this, id_sync)
            % keep epochs id_sync
            % SYNTAX   this.subSample(id_sync)
            if (nargin < 2)
                id_sync = this.id_sync;
            end
            if ~isempty(id_sync)
                this.obs(id_sync, :) = [];
            end
        end
        
        function remObs(this, id_obs)
            % remove observable with a certain id
            %
            % SYNTAX
            %   this.remObs(id_obs)
            this.obs(id_obs,:) = [];
            this.obs_code(id_obs, :) = [];
            this.active_ids(id_obs) = [];
            this.wl(id_obs) = [];
            this.f_id(id_obs) = [];
            this.prn(id_obs) = [];
            
            this.system(id_obs) = [];
            
            go_out = this.go_id(id_obs);
            this.go_id(id_obs) = [];
            go_out = setdiff(go_out, unique(this.go_id));
            
            id_out = [];
            for i = 1 : numel(id_obs)
                id_out = [id_out, find(this.ph_idx == id_obs(i))];                 %#ok<AGROW>
            end
            tmp = false(max(this.ph_idx), 1);
            tmp(this.ph_idx) = true;
            tmp(id_out) = false;
            this.ph_idx = find(tmp);
            
            % try to remove observables from other precomputed properties of the object
            if ~isempty(this.synt_ph)
                this.synt_ph(:,id_out) = [];
            end
            if ~isempty(this.sat.outlier_idx_ph)
                this.sat.outlier_idx_ph(:,id_out) = [];
            end
            if ~isempty(this.sat.cycle_slip_idx_ph)
                this.sat.cycle_slip_idx_ph(:,id_out) = [];
            end
            
            % try to remove observables from other precomputed properties of the object in sat
            
            if ~isempty(this.sat.avail_index)
                this.sat.avail_index(:, go_out) = [];
            end
            
            if ~isempty(this.sat.err_tropo)
                this.sat.err_tropo(:, go_out) = [];
            end
            if ~isempty(this.sat.err_iono)
                this.sat.err_iono(:, go_out) = [];
            end
            if ~isempty(this.sat.solid_earth_corr)
                this.sat.solid_earth_corr(:, go_out) = [];
            end
            if ~isempty(this.sat.tot)
                this.sat.tot(:, go_out) = [];
            end
            if ~isempty(this.sat.az)
                this.sat.az(:, go_out) = [];
            end
            if ~isempty(this.sat.el)
                this.sat.el(:, go_out) = [];
            end
            if ~isempty(this.sat.res)
                this.sat.res(:, go_out) = [];
            end
            if ~isempty(this.sat.slant_td)
                this.sat.slant_td(:, go_out) = [];
            end
            try
                if ~isempty(this.sat.amb)
                    this.sat.amb(:, go_out) = [];
                end
            catch
            end
        end
        
        function keepEpochs(this, good_epochs)
            % keep epochs with a certain id
            %
            % SYNTAX
            %   this.keepEpochs(good_epochs)
            bad_epochs = 1 : this.length;
            bad_epochs(good_epochs) = [];
            this.remEpochs(bad_epochs)
        end
        
        function keep(this, rate, sys_list)
            % keep epochs at a certain rate for a certain constellation
            %
            % SYNTAX
            %   this.keep(rate, sys_list)
            if nargin > 1 && ~isempty(rate)
                [~, id_sync] = GNSS_Station.getSyncTimeExpanded(this, rate);
                if ~isempty(id_sync)
                    this.keepEpochs(id_sync(~isnan(id_sync)));
                end
            end
            if nargin > 2 && ~isempty(sys_list)
                this.setActiveSys(sys_list);
                this.remSat(this.go_id(regexp(this.system, ['[^' sys_list ']'])));
            end
        end
        
        function remEpochs(this, bad_epochs)
            if ~isempty(bad_epochs)
                % remove epochs with a certain id
                %
                % SYNTAX
                %   this.remEpochs(bad_epochs)
                
                if ~isempty(this.sat.cycle_slip_idx_ph)
                    cycle_slip = this.sat.cycle_slip_idx_ph;
                    
                    tmp = false(max(bad_epochs), 1); tmp(bad_epochs) = true;
                    lim = getOutliers(tmp);
                    
                    if lim(end) == size(cycle_slip, 1)
                        lim(end,:) = [];
                    end
                    lim(:,2) = min(lim(:,2) + 1, size(cycle_slip, 1));
                    for l = 1 : size(lim, 1)
                        cycle_slip(lim(l, 2), :) = any(cycle_slip(lim(l, 1) : lim(l, 2), :));
                    end
                    cycle_slip(bad_epochs,:) = [];
                    
                    this.sat.cycle_slip_idx_ph = cycle_slip;
                end
                
                this.time.remEpoch(bad_epochs);
                
                if ~isempty(this.n_spe)
                    this.n_spe(bad_epochs) = [];
                end
                
                if ~isempty(this.synt_ph)
                    this.synt_ph(bad_epochs, :) = [];
                end
                if ~isempty(this.obs)
                    this.obs(:, bad_epochs) = [];
                end
                
                if ~isempty(this.sat.outlier_idx_ph)
                    this.sat.outlier_idx_ph(bad_epochs, :) = [];
                end
                
                if ~isempty(this.id_sync)
                    tmp = false(max(this.id_sync), 1);
                    tmp(this.id_sync) = true;
                    tmp(bad_epochs) = [];
                    this.id_sync = find(tmp);
                end
                
                if numel(this.desync) > 1
                    this.desync(bad_epochs) = [];
                end
                if numel(this.dt) > 1
                    this.dt(bad_epochs) = [];
                end
                if numel(this.dt_ph) > 1
                    this.dt_ph(bad_epochs) = [];
                end
                if numel(this.dt_pr) > 1
                    this.dt_pr(bad_epochs) = [];
                end
                if numel(this.dt_ip) > 1
                    this.dt_ip(bad_epochs) = [];
                end
                
                if ~isempty(this.apr_zhd)
                    this.apr_zhd(bad_epochs) = [];
                end
                if ~isempty(this.ztd)
                    this.ztd(bad_epochs) = [];
                end
                if ~isempty(this.zwd)
                    this.zwd(bad_epochs) = [];
                end
                if ~isempty(this.apr_zwd)
                    this.apr_zwd(bad_epochs) = [];
                end
                if ~isempty(this.pwv)
                    this.pwv(bad_epochs) = [];
                end
                if ~isempty(this.tgn)
                    this.tgn(bad_epochs) = [];
                end
                if ~isempty(this.tge)
                    this.tge(bad_epochs) = [];
                end
                
                if ~isempty(this.sat.avail_index)
                    this.sat.avail_index(bad_epochs, :) = [];
                end
                if ~isempty(this.sat.err_tropo)
                    this.sat.err_tropo(bad_epochs, :) = [];
                end
                if ~isempty(this.sat.err_iono)
                    this.sat.err_iono(bad_epochs, :) = [];
                end
                if ~isempty(this.sat.solid_earth_corr)
                    this.sat.solid_earth_corr(bad_epochs, :) = [];
                end
                if ~isempty(this.sat.tot)
                    this.sat.tot(bad_epochs, :) = [];
                end
                if ~isempty(this.sat.az)
                    this.sat.az(bad_epochs, :) = [];
                end
                if ~isempty(this.sat.el)
                    this.sat.el(bad_epochs, :) = [];
                end
                if ~isempty(this.sat.res)
                    this.sat.res(bad_epochs, :) = [];
                end
                if ~isempty(this.sat.slant_td)
                    this.sat.slant_td(bad_epochs, :) = [];
                end
            end
        end
        
        function remSat(this, go_id, prn)
            % remove satellites from receiver
            %
            % SYNTAX
            %   this.remSat(go_id)
            %   this.remSat(sys, prn);
            if nargin > 2 % interpreting as sys_c , prn
                go_id = this.getGoId(go_id, prn);
            end
            for s = 1 : numel(go_id)
                this.remObs(find(this.go_id == go_id(s))); %#ok<FNDSB>
            end
        end
        
        function remBadPrObs(this, thr)
            % remove bad pseudo-ranges 
            % and isolated observations
            %
            % INPUT
            %   thr is a threshold in meters (default = 150)
            %
            % SYNTAX
            %   this.remBadPrObs(thr)
            if nargin == 1
                thr = 150; % meters
            end
            [pr, id_pr] = this.getPseudoRanges;
            inan = isnan(pr);
            pr_fill = simpleFill1D(pr, flagExpand(~inan, 5) &  inan);
            med_pr = median(Core_Pre_Processing.diffAndPred(pr_fill,2),2,'omitnan');
            out = abs(bsxfun(@minus, Core_Pre_Processing.diffAndPred(pr_fill, 2), med_pr)) > thr;
            pr(out) = nan;
            pr = zero2nan(Core_Pre_Processing.remShortArcs(pr', 1))';
            this.setPseudoRanges(pr, id_pr);
            n_out = sum(out(:) & ~inan(:));
            this.log.addMessage(this.log.indent(sprintf(' - %d code observations marked as outlier',n_out)));
        end
        
        function [obs, sys, prn, flag] = remUndCutOff(this, obs, sys, prn, flag, cut_off)
            %  remove obs under cut off
            for i = 1 : length(prn)
                go_id = this.getGoId(sys(i), prn(i)); % get go_id
                
                idx_obs = obs(i,:) ~= 0;
                this.updateAvailIndex(idx_obs, go_id);
                XS = this.getXSTxRot(go_id);
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
        
        function remEmptyObs(this)
            %DESCRIPTION remove empty obs lines
            empty_sat = sum(abs(this.obs),2) == 0;
            this.remObs(empty_sat);
        end
        
        function remBad(this)
            % Remove observation marked as bad in crx file and satellites
            % whose prn exceed the maximum prn (spare satellites, in maintenance, etc ..)
            % remove spare satellites
            % SYNTAX
            %   this.remBad();
            
            % remove bad epoch in crx
            this.log.addMarkedMessage('Removing observations marked as bad in Bernese .CRX file')
            [CRX, found] = load_crx(this.state.crx_dir, this.time, this.cc);
            if found
                for s = 1 : size(CRX,1)
                    c_sat_idx = this.go_id == s;
                    this.obs(c_sat_idx,CRX(s,:)) = 0;
                end
            end
            
            this.log.addMarkedMessage('Removing observations for which no ephemerid or clock is present')
            nan_coord = sum(isnan(this.sat.cs.coord) | this.sat.cs.coord == 0,3) > 0;
            nan_clock = isnan(this.sat.cs.clock) | this.sat.cs.clock == 0;
            first_epoch = this.time.first;
            coord_ref_time_diff = first_epoch - this.sat.cs.time_ref_coord;
            clock_ref_time_diff = first_epoch - this.sat.cs.time_ref_clock;
            for s = 1 : this.parent.cc.getMaxNumSat()
                o_idx = this.go_id == s;
                dnancoord = diff(nan_coord(:,s));
                st_idx = find(dnancoord == 1);
                end_idx = find(dnancoord == -1);
                if ~isempty(st_idx) || ~isempty(end_idx)
                    if isempty(st_idx)
                        st_idx = 1;
                    end
                    if isempty(end_idx)
                        end_idx = length(nan_coord);
                    end
                    if end_idx(1) < st_idx(1)
                        st_idx = [1; st_idx];
                    end
                    if st_idx(end) > end_idx(end)
                        end_idx = [end_idx ; length(nan_coord)];
                    end
                    for i = 1:length(st_idx)
                        is = st_idx(i);
                        ie = end_idx(i);
                        c_rate = this.sat.cs.coord_rate;
                        bad_ep_st = min(this.time.length,max(1, floor((-coord_ref_time_diff + is * c_rate - c_rate * 10)/this.getRate())));
                        bad_ep_en = max(1,min(this.time.length, ceil((-coord_ref_time_diff + ie * c_rate + c_rate * 10)/this.getRate())));
                        this.obs(o_idx , bad_ep_st : bad_ep_en) = 0;
                    end
                else
                    if nan_coord(1,s) == 1
                        this.obs(o_idx , :) = 0;
                    end
                end
                
                % removing near empty clock -> think a better solution 
                dnanclock = diff(nan_clock(:,s));
                st_idx = find(dnanclock == 1);
                end_idx = find(dnanclock == -1);
                if ~isempty(st_idx) || ~isempty(end_idx)
                    if isempty(st_idx)
                        st_idx = 1;
                    end
                    if isempty(end_idx)
                        end_idx = length(nan_clock);
                    end
                    if end_idx(1) < st_idx(1)
                        st_idx = [1; st_idx];
                    end
                    if st_idx(end) > end_idx(end)
                        end_idx = [end_idx ; length(nan_clock)];
                    end
                    for i = 1:length(st_idx)
                        is = st_idx(i);
                        ie = end_idx(i);
                        c_rate = this.sat.cs.clock_rate;
                        bad_ep_st = min(this.time.length,max(1, floor((-clock_ref_time_diff + is*c_rate - c_rate * 1)/this.getRate())));
                        bad_ep_en = max(1,min(this.time.length, ceil((-clock_ref_time_diff + ie*c_rate + c_rate * 1)/this.getRate())));
                        this.obs(o_idx , bad_ep_st : bad_ep_en) = 0;
                    end
                else
                    if nan_clock(1,s) == 1
                        this.obs(o_idx , :) = 0;
                    end
                end
            end
            % remove moon midnight or shadow crossing epoch
            eclipsed = this.sat.cs.checkEclipseManouver(this.time);
            for i = 1: size(eclipsed,2)
                this.obs(this.go_id == i,eclipsed(:,i)~=0) = 0;
            end
            % check empty lines
            this.remEmptyObs();
            
        end
        
        function remUnderSnrThr(this, snr_thr)
            % remS observations with an SNR smaller than snr_thr
            % 
            % SYNTAX
            %   this.remUnderSnrThr(snr_thr)
            
            if (nargin == 1)
                snr_thr = this.state.getSnrThr();
            end
            
            [snr1, id_snr] = this.getObs('S1');
            
            snr1 = this.smoothSatData([],[],zero2nan(snr1'), [], 'spline', 900 / this.getRate, 10); % smoothing SNR => to be improved
            id_snr = find(id_snr);
            if ~isempty(id_snr)
                this.log.addMarkedMessage(sprintf('Removing data under %.1f dBHz', snr_thr));                
                n_rem = 0;
                for s = 1 : numel(id_snr)
                    snr_test = snr1(:, s) < snr_thr;
                    this.obs(this.go_id == (this.go_id(id_snr(s))), snr_test) = 0;
                    n_rem = n_rem + sum(snr_test(:));
                end
                this.log.addMessage(this.log.indent(sprintf(' - %d observations under SNR threshold', n_rem)));
            end
        end
        
        function remUnderCutOff(this, cut_off)
            % remS observations with an elevation lower than cut_off
            % 
            % SYNTAX
            %   this.remUnderCutOff(cut_off)
            if nargin == 1
                cut_off = this.state.getCutOff;
            end
            this.log.addMessage(this.log.indent(sprintf('Removing observations under cut-off (%d degrees)', cut_off)));
            mask = this.sat.el > cut_off;
            this.obs = this.obs .* mask(:, this.go_id)';
        end
        
        function remShortArc(this, min_arc)
            % remS arc shorter than
            % SYNTAX
            %   this.remShortArc()
            if min_arc > 1
                this.log.addMarkedMessage(sprintf('Removing arcs shorter than %d epochs', 1 + 2 * ceil((min_arc - 1)/2)));
                
                [pr, id_pr] = this.getPseudoRanges();
                idx = ~isnan(pr);
                idx_s = flagShrink(idx, ceil((min_arc - 1)/2));
                idx_e = flagExpand(idx_s, ceil((min_arc - 1)/2));
                el_idx = xor(idx,idx_e);
                pr(el_idx) = NaN;
                this.setPseudoRanges(pr, id_pr);
                this.log.addMessage(this.log.indent(sprintf(' - %d code observations have been remD', sum(el_idx(:)))));
                
                [ph, wl, id_ph] = this.getPhases();
                idx = ~isnan(ph);
                idx_s = flagShrink(idx, ceil((min_arc - 1)/2));
                idx_e = flagExpand(idx_s, ceil((min_arc - 1)/2));
                el_idx = xor(idx,idx_e);
                ph(el_idx) = NaN;
                this.setPhases(ph, wl, id_ph);
                this.log.addMessage(this.log.indent(sprintf(' - %d phase observations have been remD', sum(el_idx(:)))));
            end
        end
        
        function updateRemOutlierMarkCycleSlip(this)
            % After changing the observations Synth phases must be recomputed and
            % old outliers and cycle-slips remD before launching a new detection
            this.sat.outlier_idx_ph = false(size(this.sat.outlier_idx_ph));
            this.sat.cycle_slip_idx_ph = false(size(this.sat.cycle_slip_idx_ph));
            this.updateSyntPhases();
            this.remOutlierMarkCycleSlip();
        end
        
        function remOutlierMarkCycleSlip(this)
            this.log.addMarkedMessage('Cleaning observations');
            %% PARAMETRS
            ol_thr = 0.5; % outlier threshold
            cs_thr = 0.5; % CYCLE SLIP THR
            sa_thr = this.state.getMinArc();  % short arc threshold
            
            %----------------------------
            % Outlier Detection
            %----------------------------            
            
            % mark all as outlier and interpolate
            % get observed values
            this.sat.outlier_idx_ph = false(size(this.sat.outlier_idx_ph));
            this.sat.cycle_slip_idx_ph = false(size(this.sat.cycle_slip_idx_ph));
            [ph, wl, id_ph_l] = this.getPhases;
            
            this.log.addMessage(this.log.indent('Detect outlier candidates from residual phase time derivate'));
            % first time derivate
            synt_ph = this.getSyntPhases;
            sensor_ph = Core_Pre_Processing.diffAndPred(ph - synt_ph);
            
            % subtract median (clock error)
            sensor_ph = bsxfun(@minus, sensor_ph, median(sensor_ph, 2, 'omitnan'));
            % divide for wavelenght
            sensor_ph = bsxfun(@rdivide, sensor_ph, wl');
            
            % test sensor variance
            tmp = sensor_ph(~isnan(sensor_ph));
            tmp(abs(tmp) > 4) = [];
            std_sensor = mean(movstd(tmp(:),900));
            %%
            % if the sensor is too noisy (i.e. the a-priori position is probably not very accurate)
            % use as a sensor the time second derivate
            if std_sensor > ol_thr
                this.log.addWarning('Bad dataset, switching to second time derivative for outlier detection');
                der = 2; % use second
                % try with second time derivate
                sensor_ph = Core_Pre_Processing.diffAndPred(ph - synt_ph, der);
                % subtract median (clock error)
                sensor_ph = bsxfun(@minus, sensor_ph, median(sensor_ph, 2, 'omitnan'));
                % divide for wavelenght
                sensor_ph = bsxfun(@rdivide, sensor_ph, wl');
            else
                der = 1; % use first
            end
            % outlier when they exceed 0.5 cycle
            poss_out_idx = abs(sensor_ph) > ol_thr / der;
            % take them off
            ph2 = ph;
            ph2(poss_out_idx) = nan;
            
            %----------------------------
            % Cycle slip detection
            %----------------------------
            
            this.log.addMessage(this.log.indent('Detect cycle slips from residual phase time derivate'));
            % join the nan
            sensor_ph_cs = nan(size(sensor_ph));
            for o = 1 : size(ph2,2)
                tmp_ph = ph2(:,o);
                ph_idx = not(isnan(tmp_ph));
                tmp_ph = tmp_ph(ph_idx);
                if ~isempty(tmp_ph)
                    sensor_ph_cs(ph_idx,o) = Core_Pre_Processing.diffAndPred(tmp_ph - synt_ph(ph_idx,o), der);
                end
            end
            
            % subtract median
            sensor_ph_cs2 = bsxfun(@minus, sensor_ph_cs, median(sensor_ph_cs, 2, 'omitnan'));
            % divide for wavelenght
            sensor_ph_cs2 = bsxfun(@rdivide, sensor_ph_cs2, wl');
            
            % find possible cycle slip
            % cycle slip when they exceed threhsold cycle
            poss_slip_idx = abs(sensor_ph_cs2) > cs_thr;
            
            %check if epoch before cycle slip can be restored
            poss_rest = [poss_slip_idx(2:end,:); zeros(1,size(poss_slip_idx,2))];
            poss_rest = poss_rest & poss_out_idx;
            poss_rest_line = sum(poss_rest,2);
            if sum(poss_rest_line) > 0
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
            end
            this.ph_idx = find(id_ph_l);
            %-------------------------------------------------------
            % MELBOURNE WUBBENA based cycle slip detection
            %--------------------------------------------------------
            %             mw = this.getMelWub('1','2','G');
            %             omw = mw.obs./repmat(mw.wl,size(mw.obs,1),1);
            %             %m_omw = reshape(medfilt_mat(omw,10),size(omw));
            %             omw1 = Core_Pre_Processing.diffAndPred(zero2nan(omw),1);
            %             omw2 = Core_Pre_Processing.diffAndPred(zero2nan(omw),2);
            %             omw3 = Core_Pre_Processing.diffAndPred(zero2nan(omw),3);
            %             omw4 = Core_Pre_Processing.diffAndPred(zero2nan(omw),4);
            %             omw5 = Core_Pre_Processing.diffAndPred(zero2nan(omw),5);
            %             m_omw = (omw1 + omw2/2 + omw3/3 + omw4/4 + omw5/5) / 5;
            %             for o = 1 : length(mw.go_id)
            %                 go_id = mw.go_id(o);
            %                 idx_gi = this.go_id(id_ph_l) == go_id;
            %                 poss_slip_idx(abs(m_omw(:,o)) > 0.5,idx_gi) = 1;
            %             end
            %-------------------------------------------------------
            % SAFE CHOICE: if there is an hole put a cycle slip
            %--------------------------------------------------------
            for o = 1 : size(ph2,2)
                tmp_ph = ph2(:,o);
                ph_idx = isnan(tmp_ph);
                c_idx = [false ;diff(ph_idx) == -1];
                poss_slip_idx(c_idx,o) = 1;
            end
            
            % if majority of satellites jump set cycle slip on all
            n_obs_ep = sum(~isnan(ph2),2);
            all_but_one = (n_obs_ep - sum(poss_slip_idx,2)) < (0.7 * n_obs_ep) |  (n_obs_ep - sum(poss_slip_idx,2)) < 3;
            for c = find(all_but_one')
                poss_slip_idx(c,~isnan(ph2(c,:))) = 1;
            end
            
            % remove too short possible arc
            to_short_idx = flagMerge(poss_slip_idx,sa_thr);
            poss_slip_idx = [diff(to_short_idx) <0 ; false(1,size(to_short_idx,2))];
            to_short_idx(poss_slip_idx) = false;
            poss_out_idx(to_short_idx) = true;
            
            n_out = sum(poss_out_idx(:));
            this.sat.outlier_idx_ph = sparse(poss_out_idx);
            this.sat.cycle_slip_idx_ph = double(sparse(poss_slip_idx));
            
            % Outlier detection for some reason is not working properly -> reperform outlier detection
            sensor_ph = Core_Pre_Processing.diffAndPred(this.getPhases - this.getSyntPhases, 2);
            sensor_ph = bsxfun(@minus, sensor_ph, median(sensor_ph, 2, 'omitnan'));
            % divide for wavelenght
            sensor_ph = bsxfun(@rdivide, sensor_ph, wl');
            % outlier when they exceed 0.5 cycle
            poss_out_idx = abs(sensor_ph) > 0.5;
            this.sat.outlier_idx_ph = this.sat.outlier_idx_ph | poss_out_idx;
            
            this.sat.cycle_slip_idx_ph([false(1,size(this.sat.outlier_idx_ph,2)); (diff(this.sat.outlier_idx_ph) == -1)]) = 1;
            this.log.addMessage(this.log.indent(sprintf(' - %d phase observations marked as outlier',n_out)));
        end
        
        function cycleSlipPPPres(this)
            sensor = Core_Pre_Processing.diffAndPred(this.sat.res);
            ph_idx =  this.obs_code(:,1) == 'L';
            idx = abs(sensor) > 0.03;
            for i = 1 : size(sensor,2)
                sat_idx = this.go_id(ph_idx) == i;
                this.sat.cycle_slip_idx_ph(idx(:,i), :) = 1;
            end
        end
        
        function tryCycleSlipRepair(this)
            %----------------------------
            % Cycle slip repair
            %----------------------------
            
            % window used to estimate cycle slip
            % linear time
            lin_time = 900; %single diffrence is linear in 15 minutes
            max_window = 600; %maximum windows allowed (for computational reason)
            win_size = min(max_window,ceil(lin_time / this.getRate()/2)*2); %force even
            
            poss_out_idx = this.sat.outlier_idx_ph;
            poss_slip_idx = this.sat.cycle_slip_idx_ph;
            
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

                    s_diff = ph_other - repmat(ph_piv,1,size(ph_other,2));
                    %}
                    % repair
                    % TO DO half cycle
                    if ~isnan(jmp)
                        this.obs(id_ph(p),c:end) = nan2zero(zero2nan(this.obs(id_ph(p),c:end)) - round(jmp));
                    end
                    if false %abs(jmp -round(jmp)) < 0.01
                        poss_slip_idx(c, p) = - 1;
                        n_repaired = n_repaired +1;
                    else
                        poss_slip_idx(c, p) =   1;
                    end
                end
            end
            
            this.log.addMessage(this.log.indent(sprintf(' - %d of %d cycle slip repaired',n_repaired,n_cycleslip)));
            
            this.sat.cycle_slip_idx_ph = poss_slip_idx;
            
            % save index into object
            
        end
        
        function importRinex(this, file_name)
            % Parses RINEX observation files.
            %
            % SYNTAX
            %   this.importRinex(file_name)
            %
            % INPUT
            %   filename = RINEX observation file(s)
            %
            % OUTPUT
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
                
                % removing empty lines at end of file
                while (lim(end,1) - lim(end-1,1))  < 2
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
                % GPS L1 -> L1C
                idx = this.getObsIdx('L1 ','G');
                this.obs_code(idx,:) = repmat('L1C',length(idx),1);
                % GPS L2 -> L2W
                idx = this.getObsIdx('L2 ','G');
                this.obs_code(idx,:) = repmat('L2W',length(idx),1);
                % GPS C5 -> C5I
                idx = this.getObsIdx('C5 ','G');
                this.obs_code(idx,:) = repmat('C5I',length(idx),1);
                % GPS L5 -> L5I
                idx = this.getObsIdx('L5 ','G');
                this.obs_code(idx,:) = repmat('L5I',length(idx),1);
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
            this.remObs(~this.active_ids);
            
            this.setActiveSys(this.getAvailableSys);
        end
        
        function importAntModel(this)
            % Load and parse the antenna (ATX) file as specified into settings
            % SYNTAX
            %   this.importAntModel
            filename_pcv = this.state.getAtxFile;
            fnp = File_Name_Processor();
            this.log.addMessage(this.log.indent(sprintf('Opening file %s for reading', fnp.getFileName(filename_pcv))));
            this.pcv = Core_Utils.readAntennaPCV(filename_pcv, {this.parent.ant_type});
        end
        
        function importOceanLoading(this)
            %DESCRIPTION load ocean loading displcement matrix from
            %ocean_loading.blq if satation is present
            [this.ocean_load_disp, found] = load_BLQ( this.state.getOceanFile,{this.parent.getMarkerName4Ch});
            if not(found) && ~strcmp(this.parent.getMarkerName4Ch, this.parent.getMarkerName)
                [this.ocean_load_disp, found] = load_BLQ( this.state.getOceanFile,{this.parent.getMarkerName});
            end
            if found == 0
                this.ocean_load_disp = -1; %ocean loading parameters not found
                this.log.addWarning('Ocean loading parameters not found.');
            end
        end
        
        function importMeteoData(this)
            %DESCRIPTION load meteo data from Meteo Network object
            %and get a virtual station at receiver positions
            mn = Meteo_Network.getInstance();
            if ~isempty(mn.mds)
                this.log.addMarkedMessage('importing meteo data');
                this.meteo_data = mn.getVMS(this.parent.marker_name, this.xyz, this.getNominalTime);
            end
        end
        
        function parseRinHead(this, txt, lim, eoh)
            % Parse the header of the Observation Rinex file
            % SYNTAX
            %    this.parseRinHead(txt, nl)
            % INPUT
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
                        %throw(MException('VerifyINPUTInvalidObservationFile', 'RINEX2 is supported for GPS only dataset, please use a RINEX3 file '));
                    else
                        % GPS only RINEX2 -> ok
                    end
                else
                    % RINEX 3 file -> ok
                end
            else
                throw(MException('VerifyINPUTInvalidObservationFile', 'This observation RINEX does not contain observations'));
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
                this.parent.marker_name = 'NO_NAME';
            else
                this.parent.marker_name = strtrim(txt(lim(fln, 1) + (0:59)));
            end
            % 4) 'OBSERVER / AGENCY'
            fln = find(line2head == 4, 1, 'first'); % get field line
            if isempty(fln)
                this.parent.observer = 'goGPS';
                this.parent.agency = 'Unknown';
            else
                this.parent.observer = strtrim(txt(lim(fln, 1) + (0:19)));
                this.parent.agency = strtrim(txt(lim(fln, 1) + (20:39)));
            end
            % 5) 'REC # / TYPE / VERS'
            fln = find(line2head == 5, 1, 'first'); % get field line
            if isempty(fln)
                this.parent.number   = '000';
                this.parent.type     = 'unknown';
                this.parent.version  = '000';
            else
                this.parent.number = strtrim(txt(lim(fln, 1) + (0:19)));
                this.parent.type = strtrim(txt(lim(fln, 1) + (20:39)));
                this.parent.version = strtrim(txt(lim(fln, 1) + (40:59)));
            end
            % ignoring
            % 6) 'ANT # / TYPE'
            fln = find(line2head == 6, 1, 'first'); % get field line
            if isempty(fln)
                this.parent.ant = '';
                this.parent.ant_type = '';
            else
                this.parent.ant = strtrim(txt(lim(fln, 1) + (0:19)));
                this.parent.ant_type = strtrim(txt(lim(fln, 1) + (20:39)));
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
                this.parent.ant_delta_h = 0;
                this.parent.ant_delta_en = [0 0];
            else
                tmp = sscanf(txt(lim(fln, 1) + (0:13)),'%f')';                                                % read value
                this.parent.ant_delta_h = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 1), 0, tmp);         % check value integrity
                tmp = sscanf(txt(lim(fln, 1) + (14:41)),'%f')';                                               % read value
                this.parent.ant_delta_en = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 2), [0 0], tmp);    % check value integrity
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
                this.rin_obs_code = struct('G', rin_obs_code, 'R', rin_obs_code, 'E', rin_obs_code, 'J', rin_obs_code, 'C', rin_obs_code, 'I', rin_obs_code, 'S', rin_obs_code);
                
            end
            % 17) WAVELENGTH FACT L1/2
            % ignoring
            % 18) MARKER TYPE
            % Assuming non geodetic type as default
            this.parent.marker_type = 'NON-GEODETIC';
            fln = find(line2head == 18, 1, 'first'); % get field line
            if ~isempty(fln)
                this.parent.marker_type = strtrim(txt(lim(fln, 1) + (0:19)));
            end
            
            % 19) SYS / # / OBS TYPES
            if this.rin_type >= 3
                fln = find(line2head == 19); % get field lines
                this.rin_obs_code = struct('G',[],'R',[],'E',[],'J',[],'C',[],'I',[],'S',[]);
                if ~isempty(fln)
                    l = 1;
                    while l <= numel(fln)
                        sys = char(txt(lim(fln(l), 1)));
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
                if ~isempty(strfind(this.rin_obs_code.C, '1'))
                    this.rin_obs_code.c(this.rin_obs_code.C == '1') = '2';
                    this.log.addWarning('BeiDou band 1 is now defined as 2 -> Automatically converting the observation codes of the RINEX!');
                end
            end
            % 20) SYS / PHASE SHIFT
            fln = find(line2head == 20); % get field line
            if this.rin_type < 3
                this.ph_shift = struct('G', zeros(numel(this.rin_obs_code.G) / 3, 1));
            else
                this.ph_shift = struct('G',[],'R',[],'E',[],'J',[],'C',[],'I',[],'S',[]);
                for l = 1 : numel(fln)
                    if txt(lim(fln(l), 1)) ~= ' ' % ignoring phase shif only on subset of satellites
                        sys = char(txt(lim(fln(l), 1)));
                        
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
        
        function parseRin2Data(this, txt, lim, eoh)
            % Parse the data part of a RINEX 2 file -  the header must already be parsed
            % SYNTAX this.parseRin2Data(txt, lim, eoh)
            % remove comment line from lim
            comment_line = sum(txt(repmat(lim(1:end-2,1),1,7) + repmat(60:66, size(lim,1)-2, 1)) == repmat('COMMENT', size(lim,1)-2, 1),2) == 7;
            comment_line(1:eoh) = false;
            lim(comment_line,:) = [];
            if lim(end,3) < 32
                txt = [txt repmat(' ',1,32 - lim(end,3))];
            end
            % find all the observation lines
            t_line = find([false(eoh, 1); (txt(lim(eoh+1:end,1) + 2) ~= ' ')' & (txt(lim(eoh+1:end,1) + 3) == ' ')' & (txt(lim(eoh+1:end,1) + 28) ~= '4')' & lim(eoh+1:end,3) > 25]);
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
            prn = struct('G', gps_prn', 'R', glo_prn', 'E', gal_prn', 'J', qzs_prn', 'C', bds_prn', 'I', irn_prn', 'S', sbs_prn');
            
            % update the maximum number of rows to store
            n_obs = this.cc.gps.isActive * numel(prn.G) * numel(this.rin_obs_code.G) / 3 + ...
                this.cc.glo.isActive * numel(prn.R) * numel(this.rin_obs_code.R) / 3 + ...
                this.cc.gal.isActive * numel(prn.E) * numel(this.rin_obs_code.E) / 3 + ...
                this.cc.qzs.isActive * numel(prn.J) * numel(this.rin_obs_code.J) / 3 + ...
                this.cc.bds.isActive * numel(prn.C) * numel(this.rin_obs_code.C) / 3 + ...
                this.cc.irn.isActive * numel(prn.I) * numel(this.rin_obs_code.I) / 3 + ...
                this.cc.sbs.isActive * numel(prn.S) * numel(this.rin_obs_code.S) / 3;
            
            clear gps_prn glo_prn gal_prn qzs_prn bds_prn irn_prn sbs_prn;
            
            % order of storage
            % sat_system / obs_code / satellite
            sys_c = char(this.cc.sys_c);
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
                this.system = [this.system repmat(sys, 1, size(obs_code, 1))];
                
                f_id = obs_code(:,2);
                ss = this.cc.(char(this.cc.SYS_NAME{s} + 32));
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
            
            n_ops = numel(this.rin_obs_code.G)/3; % number of observations per satellite
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
                    this.log.addWarning(sprintf('Problematic epoch found at %s\nInspect the files to detect what went wrong!\nSkipping and continue the parsing, no action taken%s', this.time.getEpoch(e).toString, char(32*ones(this.w_bar.getBarLen(),1))));
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
            t_line = find([false(eoh, 1); (txt(lim(eoh+1:end,1)) == '>' & txt(lim(eoh+1:end,1)+31) ~= '4' )']);
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
            
            gps_prn(gps_prn > this.cc.gps.N_SAT(1)) = [];
            glo_prn(glo_prn > this.cc.glo.N_SAT(1)) = [];
            gal_prn(gal_prn > this.cc.gal.N_SAT(1)) = [];
            qzs_prn(qzs_prn > this.cc.qzs.N_SAT(1)) = [];
            bds_prn(bds_prn > this.cc.bds.N_SAT(1)) = [];
            irn_prn(irn_prn > this.cc.irn.N_SAT(1)) = [];
            sbs_prn(sbs_prn > this.cc.sbs.N_SAT(1)) = [];

            
            prn = struct('G', gps_prn', 'R', glo_prn', 'E', gal_prn', 'J', qzs_prn', 'C', bds_prn', 'I', irn_prn', 'S', sbs_prn');
            
            % update the maximum number of rows to store
            n_obs = this.cc.gps.isActive * numel(prn.G) * numel(this.rin_obs_code.G) / 3 + ...
                this.cc.glo.isActive * numel(prn.R) * numel(this.rin_obs_code.R) / 3 + ...
                this.cc.gal.isActive * numel(prn.E) * numel(this.rin_obs_code.E) / 3 + ...
                this.cc.qzs.isActive * numel(prn.J) * numel(this.rin_obs_code.J) / 3 + ...
                this.cc.bds.isActive * numel(prn.C) * numel(this.rin_obs_code.C) / 3 + ...
                this.cc.irn.isActive * numel(prn.I) * numel(this.rin_obs_code.I) / 3 + ...
                this.cc.sbs.isActive * numel(prn.S) * numel(this.rin_obs_code.S) / 3;
            
            clear gps_prn glo_prn gal_prn qzs_prn bds_prn irn_prn sbs_prn;
            
            % order of storage
            % sat_system / obs_code / satellite
            sys_c = char(this.cc.sys_c);
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
                
                n_code = numel(this.rin_obs_code.(sys)) / 3; % number of satellite system
                % transform in n_code x 3
                obs_code = reshape(this.rin_obs_code.(sys), 3, n_code)';
                % replicate obs_code for n_sat
                obs_code = serialize(repmat(obs_code, 1, n_sat)');
                obs_code = reshape(obs_code, 3, numel(obs_code) / 3)';
                
                
                prn_ss = repmat(prn.(sys)', n_code, 1);
                % discarding staellites whose number exceed the maximum ones for constellations e.g. spare satellites GLONASS
                    this.prn = [this.prn; prn_ss];
                    this.obs_code = [this.obs_code; obs_code];
                    this.n_sat = this.n_sat + n_sat;
                    this.system = [this.system repmat(sys, 1, size(obs_code, 1))];
                    
                    f_id = obs_code(:,2);
                    ss = this.cc.getSys(sys);
                    [~, f_id] = ismember(f_id, ss.CODE_RIN3_2BAND);
                    
                    ismember(this.system, this.cc.SYS_C);
                    this.f_id = [this.f_id; f_id];
                    
                    if sys == 'R'
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
    end
    % ==================================================================================================================================================
    %% METHODS GETTER - TIME
    % ==================================================================================================================================================
    
    methods
        % time
        function time = getPositionTime(this)
            % return the time of the computed positions
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getPositionTime()
            if this.isStatic()
                time = this.getCentralTime;
            else
                time = this.getTime.getCopy();
            end
        end

        % standard utility
        function toString(this)
            % Display on screen information about the receiver
            % SYNTAX this.toString();
                if ~this.isEmpty
                    fprintf('----------------------------------------------------------------------------------\n')
                    this.log.addMarkedMessage(sprintf('Receiver Work Space %s', this.parent.marker_name));
                    fprintf('----------------------------------------------------------------------------------\n')
                    this.log.addMessage(sprintf(' From     %s', this.time.first.toString()));
                    this.log.addMessage(sprintf(' to       %s', this.time.last.toString()));
                    this.log.newLine();
                    this.log.addMessage(sprintf(' Rate of the observations [s]:            %d', this.getRate()));
                    this.log.newLine();
                    this.log.addMessage(sprintf(' Maximum number of satellites seen:       %d', max(this.n_sat)));
                    this.log.addMessage(sprintf(' Number of stored frequencies:            %d', this.n_freq));
                    this.log.newLine();
                    this.log.addMessage(sprintf(' Satellite System(s) seen:                "%s"', unique(this.system)));
                    this.log.newLine();
                    
                    xyz0 = this.getAPrioriPos();
                    [enu0(1), enu0(2), enu0(3)] = cart2plan(xyz0(:,1), xyz0(:,2), xyz0(:,3));
                    static_dynamic = {'Dynamic', 'Static'};
                    this.log.addMessage(sprintf(' %s receiver', static_dynamic{this.parent.static + 1}));
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
                        this.log.newLine();
                        this.log.addMessage(sprintf('     3D distance = %+16.4f m', sqrt(sum((xyz_m - xyz0).^2))));
                    end
                    fprintf(' ----------------------------------------------------------\n')
                    this.log.addMessage(' Processing statistics:');
                    if not(isempty(this.s0_ip) || this.s0_ip == 0)
                        this.log.addMessage(sprintf('     sigma0 code positioning = %+16.4f m', this.s0_ip));
                    end
                    if not(isempty(this.s0) || this.s0 == 0)
                        this.log.addMessage(sprintf('     sigma0 PPP positioning  = %+16.4f m', this.s0));
                    end
                    fprintf(' ----------------------------------------------------------\n')
                end
        end
        
        function is_fixed = isFixed(this)
            is_fixed = this.rf.isFixed(this.parent.getMarkerName4Ch);
        end
        
        function has_apr = hasAPriori(this)
            has_apr = this.rf.hasAPriori(this.parent.getMarkerName4Ch);
        end
        
        function n_obs = getNumObservables(this)
            % get the number of observables stored in the object
            % SYNTAX n_obs = this.getNumObservables()
            n_obs = size(this.obs, 1);
        end
        
        function n_pr = getNumPseudoRanges(this)
            % get the number of epochs stored in the object
            % SYNTAX n_pr = this.getNumPseudoRanges()
            n_pr = sum(this.obs_code(:,1) == 'C');
        end
        
        function n_sat = getNumSat(this, sys_c)
            % get the number of satellites stored in the object
            %
            % SYNTAX 
            %   n_sat = getNumSat(<sys_c>)
            if nargin == 2
                n_sat = numel(unique(this.go_id( (this.system == sys_c)' & (this.obs_code(:,1) == 'C' | this.obs_code(:,1) == 'L') )));
            else
                n_sat = numel(unique(this.go_id(this.obs_code(:,1) == 'C' | this.obs_code(:,1) == 'L')));
            end
        end
        
        function n_sat = getMaxSat(sta_list, sys_c)
            % get the number of satellites stored in the object
            %
            % SYNTAX 
            %   n_sat = getNumSat(<sys_c>)
            n_sat = zeros(size(sta_list));
            
            for r = 1 : size(sta_list, 2)
                rec = sta_list(~sta_list(r).isEmpty, r);
                
                if ~isempty(rec)
                    if nargin == 2
                        n_sat(r) = max(rec.go_id( (rec.system == sys_c)' & (rec.obs_code(:,1) == 'C' | rec.obs_code(:,1) == 'L') ));
                    else
                        n_sat(r) = max(rec.go_id(rec.obs_code(:,1) == 'C' | rec.obs_code(:,1) == 'L'));
                    end
                end
            end
        end
        
        function sys_c = getAvailableSys(this)
            % get the available system stored into the object
            % SYNTAX sys_c = this.getAvailableSys()
            
            % Select only the systems present in the file
            sys_c = unique(this.system);
        end
        
        function sys_c = getActiveSys(this)
            % get the active system stored into the object
            % SYNTAX sys_c = this.getActiveSys()
            
            % Select only the systems still present in the file
            go_id = this.getActiveGoIds();
            sys_c = serialize(unique(this.getSysPrn(go_id)))';
        end
        
        function [sys_c, prn] = getSysPrn(this, go_id)
            % return sys_c and prn for a given go_id
            % SYNTAX [sys_c, prn] = this.getSysPrn(go_id)
            sys_c = char(ones(size(go_id))*32);
            prn = zeros(size(go_id));
            for i = 1:length(prn)
                tmp_id = find(this.go_id == go_id(i), 1, 'first');
                if ~isempty(tmp_id)
                    sys_c(i) = this.system(tmp_id);
                    prn(i) = this.prn(tmp_id);
                end
            end
        end        
        
        function go_id = getActiveGoIds(this)
            % return go_id present in the receiver of the active constellations
            % SYNTAX go_id = this.getActiveGoIds()
            go_id = unique(this.go_id(regexp(this.system, ['[' this.cc.sys_c ']?'])));
        end
        
        function go_id = getGoId(this, sys, prn)
            % return go_id for a given system and prn
            % SYNTAX go_id = this.getGoId(sys, prn)
            go_id = zeros(size(prn));
            if length(sys) == 1
                sys = repmat(sys,length(prn),1);
            end
            for i = 1 : length(prn)
                p = prn(i);
                s = sys(i);
                id = find((this.system == s)' & (this.prn == p), 1, 'first');
                if isempty(id)
                    go_id(i) = 0;
                else
                    go_id(i) = this.go_id(id);
                end
            end
        end
        
        function ant_id = getAntennaId(this, go_id)
            % return ant id given a go_id
            % SYNTAX ant_id = this.getAntennaId(go_id)
            [sys_c, prn] = this.getSysPrn(go_id);
            ant_id = reshape(sprintf('%c%02d', serialize([sys_c; prn]))',3, numel(go_id))';
        end
        
        function dt = getDt(this)
            dt =  this.dt(this.id_sync);
        end
        
        function dt_ip = getDtIp(this)
            dt_ip = this.dt_ip(this.id_sync);
        end
        
        function dt_pr = getDtPr(this)
            dt_pr = this.dt_pr(this.id_sync);
        end
        
        function dt_ph = getDtPh(this)
            dt_ph =  this.dt_ph(this.id_sync);
        end
        
        function dt_pp = getDtPrePro(this)
            dt_pp = this.dt_ip(this.id_sync);
        end
        
        function desync = getDesync(this)
            desync = this.desync(this.id_sync);
        end
        
        function time = getTime(this)
            % return the time stored in the object
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getTime()
            if ~isempty(this(1).id_sync)
                time = this(1).time.getEpoch(this(1).id_sync);
            else
                time = this(1).time.getCopy();
            end
        end
        
        function missing_epochs = getMissingEpochs(this)
            % return a logical array of missing (code) epochs
            %
            % SYNTAX
            %   missing_epochs = this.getMissingEpochs()
            %
            missing_epochs = all(isnan(this.getPseudoRanges)')';
        end
        
        function xyz = getAPrioriPos_mr(this)
            % return apriori position
            %
            % SYNTAX
            %   xyz = this.getAPrioriPos_mr()
            xyz = this.getAPrioriPos();
        end
        
        function xyz = getAPrioriPos(this)
            % return apriori position
            %
            % SYNTAX
            %   xyz = this.getAPrioriPos()
            xyz = this.xyz_approx; %#ok<AGROW>
            xyz = median(xyz, 1);
            if ~any(xyz) && ~isempty(this.xyz)
                xyz = median(this.getPosXYZ, 1);
            end
        end
        
         % frequencies
        
        function is_mf = isMultiFreq(this)
            is_mf = false;
            for i = 1 : this.parent.cc.getMaxNumSat()
                cur_sat_id = find(this.go_id == i, 1, 'first');
                if not(isempty(cur_sat_id))
                    sat_idx = this.getObsIdx('C',cur_sat_id);
                    sat_idx = sat_idx(this.active_ids(sat_idx));
                    if ~isempty(sat_idx)
                        % get epoch for which iono free is possible
                        freq = str2num(this.obs_code(sat_idx,2)); %#ok<ST2NM>
                        u_freq = unique(freq);
                        if length(u_freq)>1
                            is_mf = true;
                            return
                        end
                    end
                end
            end
        end
        
        function freqs = getFreqs(this, sys)
            % get presnt frequencies for system
            idx = this.system == sys;
            freq = str2num(this.obs_code(idx,2)); %#ok<ST2NM>
            freqs = unique(freq);
        end
        
        function err4el = getErr4El(this, err, go_id)
            % compute for each satellite, for each elevation the maximum of err
            % SYNTAX err4el = this.getErr4El(err, <go_id>)
            el_lim = (0 : 1 : 90)';
            n_el = size(el_lim, 1);
            n_sat = size(this.sat.el, 2);
            err4el = zeros(n_sat, n_el - 1);
            
            if nargin == 2
                go_id = (1 : n_sat)';
            end
            
            for s = unique(go_id)'
                for e = 1 : n_el - 1
                    id_el = (this.sat.el(:, s) > el_lim(e)) & (this.sat.el(:, s) <= el_lim(e + 1));
                    if sum(id_el) > 0
                        for id_s = find(go_id == s)'
                            err4el(s, e) = max([err4el(s, e); nan2zero(err(id_el, id_s))]);
                        end
                    end
                end
            end
            
            err4el = zero2nan(err4el);
        end
        
        function nominal_time = getNominalTime(this)
            % get the nominal time aka rounded time cosidering a constant
            % sampling rate
            %
            % SYNTAX
            %   nominal_time = this.getNominalTime()
            nominal_time_zero = floor(this.time.first.getMatlabTime() * 24)/24;
            rinex_time = this.time.getRefTime(nominal_time_zero);
            nominal_time = round(rinex_time / this.time.getRate) * this.time.getRate;
            ref_time = (nominal_time(1) : this.time.getRate : nominal_time(end))';
            
            % reordering observations filling empty epochs with zeros;
            nominal_time = GPS_Time(nominal_time_zero, ref_time, this.time.isGPS(), 2);
            nominal_time.toUnixTime;
        end
        
        function time_tx = getTimeTx(this, sat)
            % SYNTAX
            %   this.getTimeTx(epoch);
            %
            % INPUT
            % OUTPUT
            %   time_tx = transmission time
            %   time_tx =
            %
            % 
            %   Get Transmission time
            idx = this.sat.avail_index(:, sat) > 0;
            time_tx = this.time.getSubSet(idx);
            if isempty(this.sat.tot)
                this.updateAllTOT();
            end
            time_tx.addSeconds( - this.sat.tot(idx, sat));
        end
        
        function time_of_travel = getTOT(this)
            % SYNTAX
            %   this.getTraveltime()
            % INPUT
            % OUTPUT
            %   time_of_travel   = time of travel
            % 
            %   Compute the signal transmission time.
            time_of_travel = this.tot;
        end
        
        function dtRel = getRelClkCorr(this, sat)
            %  : get clock offset of the satellite due to
            % special relativity (eccntrcity term)
            idx = this.sat.avail_index(:,sat) > 0;
            [X,V] = this.sat.cs.coordInterpolate(this.time.getSubSet(idx),sat);
            dtRel = -2 * sum(conj(X) .* V, 2) / (Global_Configuration.V_LIGHT ^ 2); % Relativity correction (eccentricity velocity term)
        end
        
        function dtS = getDtS(this, sat)
            % SYNTAX
            %   this.getDtS(time_rx)
            %
            % INPUT
            %   time_rx   = reception time
            %
            % OUTPUT
            %   dtS     = satellite clock errors
            % 
            %   Compute the satellite clock error.
            if nargin < 2
                dtS = zeros(size(this.sat.avail_index));
                dtS(this.sat.avail_index(:,s),s) = this.sat.cs.clockInterpolate(this.time.getSubSet(this.sat.avail_index(:,s)), 1 : size(dtS,2));
            else
                idx = this.sat.avail_index(:,sat) > 0;
                if sum(idx) > 0
                    dtS = this.sat.cs.clockInterpolate(this.time.getSubSet(idx), sat);
                else
                    dtS = zeros(0,1);
                end
            end
        end
        
        function [XS_tx_r ,XS_tx] = getXSTxRot(this, go_id)
            % SYNTAX
            %   [XS_tx_r ,XS_tx] = this.getXSTxRot( sat)
            %
            % INPUT
            % sat = go-id of the satellite
            % OUTPUT
            % XS_tx = satellite position computed at trasmission time
            % XS_tx_r = Satellite postions at transimission time rotated by earth rotation occured
            % during time of travel
            % 
            %   Compute satellite positions at transmission time and rotate them by the earth rotation
            %   occured during time of travel of the signal
            if nargin > 1
                [XS_tx] = this.getXSTx(go_id);
                [XS_tx_r]  = this.earthRotationCorrection(XS_tx, go_id);
            else
                n_sat = this.parent.getMaxSat;
                XS_tx_r = zeros(this.time.length, n_sat, 3);
                for i = unique(this.go_id)'
                    [XS_tx] = this.getXSTx(i);
                    [XS_tx_r_temp]  = this.earthRotationCorrection(XS_tx, i);
                    XS_tx_r(this.sat.avail_index(:,i) ,i ,:) = permute(XS_tx_r_temp, [1 3 2]);
                end
            end
        end
        
        function [XS_loc] = getXSLoc(this, go_id)
            % SYNTAX
            %   [XS_tx_r ,XS_tx] = this.getXSLoc( sat)
            %
            % INPUT
            %   sat = go-id of the satellite
            % OUTPUT
            %   XS_tx = satellite position computed at trasmission time
            %   XS_tx_r = Satellite postions at transimission time rotated by earth rotation occured during time of travel
            % 
            %   Compute satellite positions at transmission time and rotate them by the earth rotation
            %   occured during time of travel of the signal and subtract
            %   the postion term
            n_epochs = this.time.length;
            if nargin > 1
                sat_idx = this.sat.avail_index(:, go_id) > 0;
                XS = this.getXSTxRot(go_id);
                XS_loc = nan(n_epochs, 3);
                XS_loc(sat_idx,:) = XS;
                if size(this.xyz,1) == 1
                    XR = repmat(this.xyz, n_epochs, 1);
                else
                    XR = this.xyz;
                end
                XS_loc = XS_loc - XR;
            else
                n_sat = this.parent.cc.getMaxNumSat();
                XS_loc = zeros(n_epochs, n_sat,3);
                for i = unique(this.go_id)'
                    XS_loc(:,i ,:) = this.getXSLoc(i);
                end
            end
        end
        
        function [XS_tx] = getXSTx(this, sat)
            % SYNTAX
            %   [XS_tx_frame , XS_rx_frame] = this.getXSTx()
            %
            % INPUT
            %  obs : [1x n_epochs] pseudi range observations
            %  sta : index of the satellite
            % OUTPUT
            % XS_tx = satellite position computed at trasmission time
            % 
            % Compute satellite positions at trasmission time
            time_tx = this.getTimeTx(sat);
            %time_tx.addSeconds(); % rel clok neglegible
            [XS_tx, ~] = this.sat.cs.coordInterpolate(time_tx, sat);
            
            
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
        
        function obs_code = getAvailableCode(this, sys_c)
            % get the available observation code for a specific sys
            % SYNTAX obs_code = this.getAvailableCode(sys_c)
            tmp = this.obs_code(this.system == sys_c, :);
            [~, code_ok] = unique(uint32(tmp * ([1 2^8 2^16]')));
            obs_code = tmp(code_ok, :);
        end
        
        function is_static = isStatic(this)
            % return true if the receiver is static
            % SYNTAX is_static = this.isStatic()
            is_static = this.parent.static;
        end
        
        function is_pp = isPreProcessed(this)
            is_pp = this.pp_status;
        end
        
        function [iono_mf] = getSlantIonoMF(this)
            % Get Iono mapping function for all the valid elevations
            [lat, ~, ~, h_ortho] = this.getMedianPosGeodetic_mr();
            
            atmo = Atmosphere.getInstance();
            [iono_mf] = atmo.getIonoMF(lat ./180 * pi, h_ortho, this.getEl ./180 * pi);            
        end        
        
        function [mfh, mfw] = getSlantMF(this, id_sync)
            % Get Mapping function for the satellite slant
            %
            % OUTPUT
            %   mfh: hydrostatic mapping function
            %   mfw: wet mapping function
            %
            % SYNTAX
            %   [mfh, mfw] = this.getSlantMF()
            

            n_sat = this.parent.getMaxSat;
            
            mfh = zeros(this.length, n_sat);
            mfw = zeros(this.length, n_sat);
            
            t = 1;
            
            atmo = Atmosphere.getInstance();
            [lat, lon, h_ellipse, h_ortho] = this.getMedianPosGeodetic();
            lat = median(lat);
            lon = median(lon);
            h_ortho = median(h_ortho);
            if nargin == 1                
                for s = 1 : numel(this)
                    if ~isempty(this(s))
                        if this(s).state.mapping_function == 2
                            [mfh_tmp, mfw_tmp] = atmo.vmf_grd(this(s).time, lat./180*pi, lon./180*pi, (this(s).sat.el)./180*pi,h_ellipse);
                        elseif this(s).state.mapping_function == 1
                            [mfh_tmp, mfw_tmp] = atmo.gmf(this(s).time, lat./180*pi, lon./180*pi, h_ortho, (this(s).sat.el)./180*pi);
                        end
                        
                        if isempty(this(s).id_sync)
                            mfh(t : t + this(s).length - 1, 1 : size(this(s).sat.el, 2)) = mfh_tmp;
                            mfw(t : t + this(s).length - 1, 1 : size(this(s).sat.el, 2)) = mfw_tmp;
                        else
                            mfh(t : t + this(s).length - 1, 1 : size(this(s).sat.el, 2)) = mfh_tmp(this(s).id_sync, :);
                            mfw(t : t + this(s).length - 1, 1 : size(this(s).sat.el, 2)) = mfw_tmp(this(s).id_sync, :);
                        end
                        t = t + this(s).length;
                    end
                end
            else
                for s = 1 : numel(this)
                    if ~isempty(this(s))
                        if this.state.mapping_function == 2
                            [mfh_tmp, mfw_tmp] = atmo.vmf_grd(this(s).time, lat./180*pi, lon./180*pi, (this(s).sat.el)./180*pi,h_ellipse);
                        elseif this.state.mapping_function == 1
                            [mfh_tmp, mfw_tmp] = atmo.gmf(this(s).time, lat./180*pi, lon./180*pi, h_ortho, (this(s).sat.el)./180*pi);
                        end
                        
                        
                        mfh(t : t + length(id_sync) - 1, 1 : size(this(s).sat.el, 2)) = mfh_tmp(id_sync, :);
                        mfw(t : t + length(id_sync) - 1, 1 : size(this(s).sat.el, 2)) = mfw_tmp(id_sync, :);
                        t = t + this(s).length;
                    end
                end
            end
        end
        
        function [P, T, H] = getPTH(this, use_id_sync)
            % Get Pressure temperature and humidity at the receiver location
            %
            % OUTPUT
            %   P : pressure [hPa] [mbar]
            %   T : celsius degree
            %   H : %
            %
            % SYNTAX
            %   [P, T, H] = getPTH(this,time,flag)
            flag = this.state.meteo_data;
            if nargin < 2
                use_id_sync = false;
            end
            time = this.time;
            l = time.length;
            switch flag
                case 1 % standard atmosphere
                    atmo = Atmosphere.getInstance();
                    this.updateCoordinates();
                    Pr = atmo.STD_PRES;
                    % temperature [K]
                    Tr = atmo.STD_TEMP;
                    % humidity [%]
                    Hr = atmo.STD_HUMI;
                    h = this.h_ortho;
                    if this.isStatic()
                        P = zeros(l,1);
                        T = zeros(l,1);
                        H = zeros(l,1);
                        P(:) = Pr * (1-0.0000226*h).^5.225;
                        T(:) = Tr - 0.0065*h;
                        H(:) = Hr * exp(-0.0006396*h);
                    else
                        P = Pr * (1-0.0000226*h).^5.225;
                        T = Tr - 0.0065*h;
                        H = Hr * exp(-0.0006396*h);
                    end
                case 2 % gpt
                    atmo = Atmosphere.getInstance();
                    this.updateCoordinates();
                    time = time.getGpsTime();
                    if this.isStatic()                        
                        [P,T,undu] = atmo.gpt( time, this.lat/180*pi, this.lon/180*pi, this.h_ellips, this.h_ellips - this.h_ortho);
                        H = zeros(l,1);
                        H(:) = atmo.STD_HUMI;% * exp(-0.0006396*this.h_ortho);
                    else
                        P = zeros(l,1);
                        T = zeros(l,1);
                        for i = 1 : l
                            [P(l), T(l), undu] = atmo.gpt( time(l), this.lat(min(l, numel(this.lat)))/180*pi, this.lon(min(l, numel(this.lat)))/180*pi, this.h_ellips(min(l, numel(this.lat))), this.h_ellips(min(l, numel(this.lat))) - this.h_ortho(min(l, numel(this.lat))));
                        end
                        H = atmo.STD_HUMI;%* exp(-0.0006396*this.h_ortho);
                    end
                case 3 % local meteo data
                    atmo = Atmosphere.getInstance();
                    
                    P = this.meteo_data.getPressure(time);
                    T = this.meteo_data.getTemperature(time);
                    H = this.meteo_data.getHumidity(time);
                    
                    if any(isnan(P)) || any(isnan(T))
                        this.updateCoordinates();
                        time = time.getGpsTime();
                        if this.isStatic()
                            [P,T,undu] = atmo.gpt( time, this.lat/180*pi, this.lon/180*pi, this.h_ellips, this.h_ellips - this.h_ortho);
                        else
                            P = zeros(l,1);
                            T = zeros(l,1);
                            for i = 1 : l
                                [P(l), T(l), undu] = atmo.gpt( time(l), this.lat(l)/180*pi, this.lon(l)/180*pi, this.h_ellips(l), this.h_ellips(l) - this.h_ortho(l));
                            end
                        end
                    end
                        
                    if any(isnan(H))
                        % Fall back to std values
                        H(isnan(H)) = atmo.STD_HUMI;% * exp(-0.0006396*this.h_ortho);
                    end
            end
            if use_id_sync
                P = P(this.getIdSync);
                T = T(this.getIdSync);
                H = H(this.getIdSync);
            end
        end
        
        function updateAprTropo(this)
            % update arpiori z tropo delays
            if isempty(this.tge) || all(isnan(this.tge))
                this.tge = zeros(this.length,1);
                this.tgn = zeros(this.length,1);
            end
            if this.state.zd_model == 1 % no vmf gridded value
                [P,T,H] = this.getPTH();
                updateAprZhd(this,P,T,H);
                updateAprZwd(this,P,T,H);
            else
                updateAprZhd(this);
                updateAprZwd(this);
            end
        end
        
        function updateAprZhd(this,P,T,H)
            %update zhd
            if isempty(this.lat)
                this.updateCoordinates()
            end            
            switch this.state.zd_model
                case 1 %% saastamoinen
                    if nargin < 2
                        [P, T, H] = this.getPTH();
                    end
                    h = this.h_ortho;
                    lat = this.lat;
                    this.apr_zhd = Atmosphere.saast_dry(P,h,lat);
                case 2 %% vmf gridded
                    atmo = Atmosphere.getInstance();
                    this.apr_zhd = atmo.getVmfZhd(this.time.getGpsTime, this.lat, this.lon, this.h_ellips);
            end
            
        end
        
        function updateAprZwd(this,P,T,H)
            %update zwd
            if isempty(this.lat)
                this.updateCoordinates()
            end
            switch this.state.zd_model
                case 1 %% saastamoinen
                    if nargin < 2
                        [P, T, H] = this.getPTH();
                    end
                    h = this.h_ortho;
                    this.apr_zwd = Atmosphere.saast_wet(T,H,h); 
                case 2 %% vmf gridded
                    atmo = Atmosphere.getInstance();
                    this.apr_zwd = atmo.getVmfZwd(this.time.getGpsTime, this.lat, this.lon, this.h_ellips);
            end
        end
        
        function slant_td = getSlantTD(this)
            % Get the slant total delay
            % SYNTAX
            %   slant_td = this.getSlantTD();
            if numel(this) == 1
                n_sat = size(this.sat.slant_td, 2);
            else
                n_sat = this.parent.cc.getMaxNumSat();
            end
            slant_td = zeros(this.length, n_sat);
            t = 1;
            for r = 1 : numel(this)
                slant_td(t : t + this(r).length - 1, 1 : this(r).parent.getMaxSat) = this(r).sat.slant_td(this(r).id_sync, :);
                t = t + this(r).length;
            end
        end
        
        function swtd = getSlantZWD(this, smooth_win_size, id_extract)
            % Get the "zenithalized" wet delay
            % SYNTAX
            %   sztd = this.getSlantZWD(<flag_smooth_data = 0>)
            
            if nargin < 3
                id_extract = 1 : this.getTime.length;
            end
            
            if ~isempty(this(1).zwd)
                [mfh, mfw] = this.getSlantMF();
                swtd = (zero2nan(this.getSlantTD) - bsxfun(@times, mfh, this.getAprZhd)) ./ mfw;
                swtd(swtd <= 0) = nan;
                swtd = swtd(id_extract, :);
                
                if nargin >= 2 && smooth_win_size > 0
                    t = this.getTime.getEpoch(id_extract).getRefTime;
                    for s = 1 : size(swtd,2)
                        id_ok = ~isnan(swtd(:, s));
                        if sum(id_ok) > 3
                            lim = getOutliers(id_ok);
                            lim = limMerge(lim, 2*smooth_win_size);
                            
                            %lim = [lim(1) lim(end)];
                            for l = 1 : size(lim, 1)
                                if (lim(l, 2) - lim(l, 1) + 1) > 3
                                    id_ok = lim(l, 1) : lim(l, 2);
                                    ztd = this.getZwd();
                                    swtd(id_ok, s) = splinerMat(t(id_ok), swtd(id_ok, s) - zero2nan(ztd(id_ok)), smooth_win_size, 0.05) + zero2nan(ztd(id_ok));
                                end
                            end
                        end
                    end
                end
            else
                this.log.addWarning('ZWD and slants have not been computed');
            end
        end
        
        function [pos] = getXR(this)
            % get xyz or the same repeated by the number of epoch
            % SYNTAX [XR] = getXR(this)
            if size(this.xyz, 1) == 1
                pos = repmat(this.xyz,this.time.length, 1);
            else
                pos = this.xyz;
            end
        end
        
        function wl = getWavelenght(this, id_ph)
            % get the wavelength of a specific phase observation
            % SYNTAX wl = this.getWavelenght(id_ph)
            wl = this.wl(id_ph);
        end               
        
        function obs_code = getAvailableObsCode(this)
            % Get all the obbservation codes stored into the receiver
            %
            % SYNTAX
            %   obs_code = thisgetAvailableObsCode();
            
            obs_code = this.obsNum2Code(unique(this.obsCode2Num(this.obs_code)));
        end
        
        function [ph, wl, id_ph] = getPhases(this, sys_c)
            % get the phases observations in meter (not cycles)
            % SYNTAX [ph, wl, id_ph] = this.getPhases(<sys_c>)
            % SEE ALSO: setPhases getPseudoRanges setPseudoRanges
            
            id_ph = this.obs_code(:, 1) == 'L';
            if (nargin == 2)
                id_ph = id_ph & (this.system == sys_c)';
            end
            ph = this.obs(id_ph, :);
            if not(isempty(this.sat.outlier_idx_ph))
                ph(this.sat.outlier_idx_ph') = nan;
            end
            wl = this.wl(id_ph);
            
            ph = bsxfun(@times, zero2nan(ph), wl)';
        end
        
        function [pr, id_pr] = getPseudoRanges(this, sys_c)
            % get the pseudo ranges observations in meter (not cycles)
            % SYNTAX [pr, id_pr] = this.getPseudoRanges(<sys_c>)
            % SEE ALSO: getPhases setPhases setPseudoRanges
            
            id_pr = this.obs_code(:, 1) == 'C';
            if (nargin == 2)
                id_pr = id_pr & (this.system == sys_c)';
            end
            pr = zero2nan(this.obs(id_pr, :)');
        end
        
        function [dop, wl, id_dop] = getDoppler(this, sys_c)
            % get the doppler observations
            % SYNTAX [dop, id_dop] = this.getDoppler(<sys_c>)
            % SEE ALSO: setDoppler
            
            id_dop = this.obs_code(:, 1) == 'D';
            if (nargin == 2)
                id_dop = id_dop & (this.system == sys_c)';
            end
            dop = zero2nan(this.obs(id_dop, :)');
            wl = this.wl(id_dop);
            dop = bsxfun(@times, zero2nan(dop), wl');
        end
        
        function [snr, id_snr] = getSNR(this, sys_c, freq_c)
            % get the SNR of the observations
            % SYNTAX [dop, id_dop] = this.getSNR(<sys_c>)
            if isempty(this.obs_code)
                snr = [];
                id_snr = [];
            else
                if nargin < 3                    
                    id_snr = this.obs_code(:, 1) == 'S';
                else
                    id_snr = this.obs_code(:, 1) == 'S' & this.obs_code(:, 2) == freq_c;
                end
                if (nargin >= 2) && ~isempty(sys_c)
                    id_snr = id_snr & (this.system == sys_c)';
                end
                snr = zero2nan(this.obs(id_snr, :)');
            end
        end
        
        function [quality] = getQuality(this, type)
            % get a quality index of the satellite data
            % type could be:
            %   'snr'
            %   '...'
            %
            % SYNTAX
            %  [quality] = this.getQuality()
            
            if nargin < 2
                type = 'snr';
            end
            switch lower(type)
                case 'snr'
                    quality = [];
                    for c = this.cc.SYS_C(this.cc.getActive)
                        [snr, id_snr] = this.getSNR(c,'1');
                        snr_temp = nan(size(snr,1),this.cc.getMaxNumSat(c));
                        [~, prn] = this.getSysPrn(this.go_id(id_snr));
                        snr_temp(:, prn) = snr;
                        quality = [quality snr_temp]; %#ok<AGROW>
                    end
            end
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
            obs = zero2nan(this.obs(idx,:));
        end
        
        function [idx] = getObsIdx(this, flag, system, prn)
            % get observation index corresponfing to the flag
            % SYNTAX this.getObsIdx(flag, <system>)
            %        this.getObsIdx(flag, <system>, <prn>)
            %        this.getObsIdx(flag, <go_id>)
            go_id_imp = false;
            idx = sum(this.obs_code(:,1:length(flag)) == repmat(flag,size(this.obs_code,1),1),2) == length(flag);
            if nargin > 2
                if ~ischar(system) % if is not a character is a go_id
                    [system, prn] = this.getSysPrn(system);
                    go_id_imp = true;
                end
                idx = idx & (this.system == system)';
            end
            if nargin > 3 || go_id_imp
                idx = idx & reshape(this.prn == prn, length(this.prn),1);
            end
            idx = find(idx);
            idx(idx == 0) = [];
        end
        
        function obs_set = getPrefObsSetCh(this, flag, system)
            [obs, idx, snr, cycle_slips] = this.getPrefObsCh(flag, system, 1);
            go_ids = this.getGoId(system, this.prn(idx)');
            if ~isempty(this.sat.el) && ~isempty(this.sat.az)
                el = this.sat.el(:, go_ids);
                az = this.sat.az(:, go_ids);
            else
                az = [];
                el = [];
            end
            obs_set = Observation_Set(this.time.getCopy(), obs' ,[this.system(idx)' this.obs_code(idx,:)], this.wl(idx)', el, az, this.prn(idx)');
            obs_set.cycle_slip = cycle_slips';
            obs_set.snr = snr';
            sigma = this.rec_settings.getStd(system, obs_set.obs_code(1,2:4));
            obs_set.sigma = sigma*ones(size(obs_set.prn));
            
        end
        
        function [obs, idx, snr, cycle_slips] = getPrefObsCh(this, flag, system, max_obs_type)
            % get observation index corresponfing to the flag using best
            % channel according to the definition in GPS_SS, GLONASS_SS, ...
            % SYNTAX this.getObsIdx(flag, <system>)
            obs = [];
            idx = [];
            snr = [];
            cycle_slips = [];
            if length(flag)==3
                idx = sum(this.obs_code == repmat(flag,size(this.obs_code,1),1),2) == 3;
                idx = idx & (this.system == system)';
                %this.legger.addWarning(['Unnecessary Call obs_type already determined, use getObsIdx instead'])
                [obs,idx] = this.getObs(flag, system);
            elseif length(flag) >= 2
                flags = zeros(size(this.obs_code,1),3);
                sys_idx = (this.system == system)';
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
                % find the better flag present
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
                    idxes = [idxes  (sum(this.obs_code == flags,2) == 3) & this.system' == system];
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
                for s = 1:length(prn) % for each satellite and each epoch find the best (but put them on different lines)
                    sat_idx = sys_idx & this.prn==prn(s);
                    
                    tmp_obs = zeros(n_opt,n_epochs);
                    take_idx = ones(1,n_epochs)>0;
                    for i = 1 : n_opt
                        c_idx = idxes(:, i) & sat_idx;
                        snr_idx = (idxCharLines(this.obs_code, ['S' complete_flags(i,2:3)]) | idxCharLines(this.obs_code, ['S' complete_flags(i,2) ' '])  ) & sat_idx;
                        if sum(c_idx)>0
                            obs((s-1)*n_opt+i,take_idx) = this.obs(c_idx,take_idx);
                            if sum(snr_idx)
                                snr((s-1)*n_opt+i,take_idx) = this.obs(snr_idx,take_idx);
                            else
                                snr((s-1)*n_opt+i,take_idx) = 1;
                            end
                            if ~isempty(this.sat.outlier_idx_ph) && flag(1) == 'L'% take off outlier
                                ph_idx = this.ph_idx == find(c_idx);
                                obs((s-1)*n_opt+i,this.sat.outlier_idx_ph(:,ph_idx)) = 0;
                                snr((s-1)*n_opt+i,this.sat.outlier_idx_ph(:,ph_idx)) = 0;
                                cycle_slips((s-1)*n_opt+i,:) = this.sat.cycle_slip_idx_ph(:,ph_idx)';
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
        
        function [obs_set] = getTwoFreqComb(this, flag1, flag2, fun1, fun2, know_comb)
            %INPUT flag1 : either observation code (e.g. GL1) or observation set
            %       flag1 : either observation code (e.g. GL2) or observation set
            %       fun1  : function of the two wavelegnth to be applied to
            %       compute the first frequncy coefficient
            %       fun2  : function of the two wavelegnth to be applied to
            %       compute the second frequency coefficient
            if nargin  < 6
                know_comb = 'none';
            end
            if ischar(flag1)
                system = flag1(1);
                [o1, i1, s1, cs1] = this.getPrefObsCh(flag1(2:3), system, 1);
                o1 = o1';
                s1 = s1';
                cs1 = cs1';
                [o2, i2, s2, cs2] = this.getPrefObsCh(flag2(2:3), system, 1);
                o2 = o2';
                s2 = s2';
                cs2 = cs2';
                % get prn for both frequency
                p1 =  this.prn(i1);
                p2 =  this.prn(i2);
                w1 =  this.wl(i1);
                w2 =  this.wl(i2);
                oc1 = this.obs_code(i1(1),:);
                oc2 = this.obs_code(i2(1),:);
                % form cycle to meters
                if flag1(2) == 'L'
                    o1 = o1.*repmat(w1',size(o1,1),1);
                end
                if flag2(2) == 'L'
                    o2 = o2.*repmat(w2',size(o2,1),1);
                end
                sigma1 = repmat(this.rec_settings.getStd(system,this.obs_code(i1(1),:)),size(o1,2),1);
                sigma2 = repmat(this.rec_settings.getStd(system,this.obs_code(i2(1),:)),size(o2,2),1);
            elseif strcmp(class(flag1),'Observation_Set') %observation set
                o1 = flag1.obs;
                o2 = flag2.obs;
                s1 = flag1.snr;
                s2 = flag2.snr;
                sigma1 = flag1.sigma;
                sigma2 = flag2.sigma;
                cs1 = flag1.cycle_slip;
                cs2 = flag2.cycle_slip;
                p1 =  flag1.prn;
                p2 =  flag2.prn;
                w1 =  flag1.wl;
                w2 =  flag2.wl;
                obs_code = [];
                system = flag1.obs_code(1,1);
                oc1 = '   ';
                oc2 = '   ';
            end
            common_prns = intersect(p1, p2)';
            go_id = zeros(size(common_prns))';
            obs_out = zeros(size(o1,1), length(common_prns));
            snr_out = zeros(size(o1,1), length(common_prns));
            cs_out = zeros(size(o1,1), length(common_prns));
            az = zeros(size(o1,1), length(common_prns));
            el = zeros(size(o1,1), length(common_prns));
            obs_code = repmat([system oc1 oc2],length(common_prns),1);
            wl = zeros(1,length(common_prns),1);
            sigma = zeros(1,length(common_prns));
            for p = 1 : length(common_prns)
                ii1 = p1 == common_prns(p);
                ii2 = p2 == common_prns(p);
                alpha1 = fun1(w1(ii1), w2(ii2));
                alpha2 = fun2(w1(ii1), w2(ii2));
                obs_out(:, p) = nan2zero(alpha1 * zero2nan(o1(:,ii1)) + alpha2 * zero2nan(o2(:,ii2)));
                idx_obs = obs_out(:, p) ~=0;
                % get go id for the current sat
                c_go_id = this.getGoId(system, common_prns(p));
                go_id(p) = c_go_id;
                if ~isempty(this.sat.az)
                    az(idx_obs, p) = this.sat.az(idx_obs, c_go_id);
                else
                    az = [];
                end
                if ~isempty(this.sat.el)
                    el(idx_obs, p) = this.sat.el(idx_obs, c_go_id);
                else
                    el = [];
                end
                snr_out(:, p) = nan2zero(sqrt((alpha1.* zero2nan(s1(:, ii1))).^2 + (alpha2 .* zero2nan(s2(:, ii2))).^2));
                sigma(p) = sqrt((alpha1*sigma1(ii1))^2 + (alpha2*sigma2(ii2))^2);
                if isempty(cs1) && isempty(cs2) % cycle slips only if there is at least one phase observables
                    cs_out = [];
                    wl(p) = -1;
                elseif isempty(cs1)
                    cs_out = cs2;
                    wl(p) = alpha2 * w2(ii2);
                elseif isempty(cs2)
                    cs_out = cs1;
                    wl(p) = alpha1 * w1(ii1);
                else
                    cs_out(:, p) = cs2(:, ii2) | cs1(:, ii1);
                    if strcmp(know_comb,'IF')
                        [i,j] = rat(w2(ii2)/w1(ii1),0.001);
                        j= -j;
                    elseif strcmp(know_comb,'NL')
                        i = 1;
                        j = 1;
                    elseif strcmp(know_comb,'WL')
                        i = 1;
                        j = -1;
                    else
                        i = -w2(ii2)/2;
                        j = -w1(ii1)/2; %so wavelength is -1
                    end
                    wl(p) = w1(ii1)*w2(ii2)/(i*w2(ii2) + j*w1(ii1)); % see: https://www.tekmon.gr/online-gps-tutorial/1-4-1-linear-combinations-of-simultaneous-observations-between-frequencies
                end
            end
            obs_set = Observation_Set(this.time.getCopy(), obs_out ,obs_code, wl, el, az, common_prns);
            obs_set.go_id = go_id;
            obs_set.cycle_slip = cs_out;
            obs_set.snr = snr_out;
            obs_set.sigma = sigma;
        end
        
        function [obs_set] = getPrefObsCh_os(this, flag, system)
            [o, i, s, cs] = this.getPrefObsCh(flag, system, 1);
            el = this.getEl();
            az = this.getAz();
            if ~isempty(el)
                el = el(: ,this.go_id(i));
                az = az(: ,this.go_id(i));
            end
            ne = size(o,2);
            obs_set = Observation_Set(this.time.getCopy(), o' .* repmat(this.wl(i,:)',ne,1) , [this.system(i)' this.obs_code(i,:)], this.wl(i,:), el, az, this.prn(i));
            obs_set.cycle_slip = cs';
            obs_set.snr = s';
            obs_set.sigma = this.rec_settings.getStd(system, [flag '_'])*ones(size(this.prn(i)));
        end
        
        % ---------------------------------------
        %  Two Frequency observation combination
        %  Warappers of getTwoFreqComb function
        % ---------------------------------------
        
        function [obs_set] = getIonoFree(this,flag1,flag2,system)
            fun1 = @(wl1, wl2) wl2 ^ 2 / (wl2 ^ 2 - wl1 ^ 2);
            fun2 = @(wl1, wl2) - wl1^2/(wl2 ^ 2 - wl1 ^ 2);
            [obs_set] =  this.getTwoFreqComb([system flag1], [system flag2], fun1, fun2,'IF');
            obs_set.obs_code = [obs_set.obs_code repmat('I', size(obs_set.obs_code,1), 1)];
        end
        
        function [obs_set] = getNarrowLane(this,flag1,flag2,system)
            fun1 = @(wl1,wl2) wl2/(wl2+wl1);
            fun2 = @(wl1,wl2) wl1/(wl2+wl1);
            [obs_set] =  this.getTwoFreqComb([system flag1], [system flag2], fun1, fun2,'NL');
        end
        
        function [obs_set] = getWideLane(this,flag1,flag2,system)
            fun1 = @(wl1,wl2) wl2/(wl2-wl1);
            fun2 = @(wl1,wl2) - wl1/(wl2-wl1);
            [obs_set] =  this.getTwoFreqComb([system flag1],[system flag2], fun1, fun2,'WL');
        end
        
        function [obs_set] = getGeometryFree(this,flag1,flag2,system)
            if flag1(1) == 'C' &&  flag2(1) == 'C'
                % Code geometry free is C2 - C1 -> exchange the flags
                tmp = flag1;
                flag1 = flag2;
                flag2 = tmp;
            end
            fun1 = @(wl1,wl2) 1;
            fun2 = @(wl1,wl2) -1;
            [obs_set] =  this.getTwoFreqComb([system flag1],[system flag2], fun1, fun2);
        end
        
        function [obs_set] = getMelWub(this, freq1, freq2, system)
            
            fun1 = @(wl1,wl2) 1;
            fun2 = @(wl1,wl2) -1;
            [obs_set1] = this.getWideLane(['L' freq1],['L' freq2], system); %widelane phase
            [obs_set2] = this.getNarrowLane(['C' freq1],['C' freq2], system); %narrowlane code
            [obs_set] =  this.getTwoFreqComb(obs_set1, obs_set2, fun1, fun2);
        end
        
        function [obs_set]  = getPrefIonoFree(this, obs_type, system)
            % get Preferred Iono free combination for the two selcted measurements
            % SYNTAX [obs] = this.getIonoFree(flag1, flag2, system)
            iono_pref = this.cc.getSys(system).IONO_FREE_PREF;
            is_present = false(size(iono_pref,1),1);
            for i = 1 : size(iono_pref,1)
                % check if there are observation for the selected channel
                if sum(iono_pref(i,1) == this.obs_code(:,2) & obs_type == this.obs_code(:,1)) > 0 && sum(iono_pref(i,2) == this.obs_code(:,2) & obs_type == this.obs_code(:,1)) > 0
                    is_present(i) = true;
                end
            end
            iono_pref = iono_pref(is_present,:);
            [obs_set]  = this.getIonoFree([obs_type iono_pref(1,1)], [obs_type iono_pref(1,2)], system);
        end
        
        function [obs_set]  = getSmoothIonoFreeAvg(this, obs_type, sys_c)
            % get Preferred Iono free combination for the two selcted measurements
            % SYNTAX [obs] = this.getIonoFree(flag1, flag2, system)
                        
            [ismf_l1]  = this.getSmoothIonoFree([obs_type '1'], sys_c);
            [ismf_l2]  = this.getSmoothIonoFree([obs_type '2'], sys_c);
            
            id_ph = Core_Utils.code2Char2Num(this.getAvailableObsCode) == Core_Utils.code2Char2Num('L5');
            %if any(id_ph)
            %    [ismf_l5]  = this.getSmoothIonoFree([obs_type '5'], sys_c);
            %end
            
            fun1 = @(wl1,wl2) 0.65;
            fun2 = @(wl1,wl2) 0.35;
            [obs_set] =  this.getTwoFreqComb(ismf_l1, ismf_l2, fun1, fun2);
            obs_set.iono_free = true;
            obs_set.obs_code = repmat('GL1L2AV',length(obs_set.prn),1);
        end
        
        function [obs_set]  = getSmoothIonoFree(this, obs_type, sys_c)
            % get Preferred Iono free combination for the two selcted measurements
            % SYNTAX [obs] = this.getIonoFree(flag1, flag2, system)
            
            % WARNING -> AS now it works only with 1� and 2� frequency
            
            
            [gf] = this.getGeometryFree('L1', 'L2', sys_c); %widelane phase
            
            gf.obs = smoothSatData(this, [], [], zero2nan(gf.obs), gf.cycle_slip);
            
            [obs_set1] = getPrefObsCh_os(this, obs_type, sys_c);
            ifree = this.cc.getSys(sys_c).getIonoFree();
            coeff = [ifree.alpha2 ifree.alpha1]; 
            fun1 = @(wl1, wl2) 1;
            band = this.cc.getBand(sys_c, obs_type(2));
            fun2 = @(wl1, wl2) + coeff(band);
            [obs_set] =  this.getTwoFreqComb(obs_set1, gf, fun1, fun2);
            obs_set.iono_free = true;
            synt_ph = this.getSyntTwin(obs_set);
            %%% tailored outlier detection
            sensor_ph = Core_Pre_Processing.diffAndPred(zero2nan(obs_set.obs) - zero2nan(synt_ph));
            sensor_ph = sensor_ph./(repmat(serialize(obs_set.wl)',size(sensor_ph,1),1));
%             
%             % subtract median (clock error)
             sensor_ph = bsxfun(@minus, sensor_ph, median(sensor_ph, 2, 'omitnan'));
% 
%             
%             % outlier when they exceed 0.5 cycle
            poss_out_idx = abs(sensor_ph) > 0.5;
            poss_out_idx = poss_out_idx & ~(obs_set.cycle_slip);
            obs_set.obs(poss_out_idx) = 0;
            obs_set.obs_code = repmat('GL1',length(obs_set.prn),1);
        end        
        
        function [obs_set]  = getPrefMelWub(this, system)
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
            [obs_set]  = this.getMelWub([iono_pref(1,1)], [iono_pref(1,2)], system);
        end
        
        function [rec_out_ph] = getOOutPh(this)
            % get the phase outlier to be putted in the results
            % default outiler are the ones on the iono free combination
            if ~isempty( this.sat.o_out_ph)
                rec_out_ph = this.sat.o_out_ph(this.getIdSync(),:);
            else
                rec_out_ph = zeros(length(this.getIdSync()), this.cc.getMaxNumSat());
            end
        end
        
        function [rec_cs_ph] = getOCsPh(this)
            % get the phase outlier to be putted in the results
            % default outiler are the ones on the iono free combination
            if ~isempty( this.sat.o_cs_ph)
                rec_cs_ph = this.sat.o_cs_ph(this.getIdSync(),:);
            else
                rec_cs_ph = zeros(length(this.getIdSync()), this.cc.getMaxNumSat());
            end
        end
        
        function [range, XS_loc] = getSyntObs(this, go_id_list)
            %  get the estimate of one measurmenet based on the
            % current postion
            % INPUT
            % EXAMPLE:
            %   this.getSyntObs(1,go_id)
            n_epochs = size(this.obs, 2);
            if nargin < 2
                go_id_list = 1 : this.parent.cc.getMaxNumSat();
            end
            range = zeros(length(go_id_list), n_epochs);
            XS_loc = [];
            for i = 1 : length(go_id_list)
                go_id = go_id_list(i);
                XS_loc = this.getXSLoc(go_id);
                range_tmp = sqrt(sum(XS_loc.^2,2));
                range_tmp = range_tmp + nan2zero(this.sat.err_tropo(:,go_id));
                
                XS_loc(isnan(range_tmp),:) = [];
                %range = range';
                range(i,:) = nan2zero(range_tmp)';
            end
        end
        
        function synt_pr_obs = getSyntCurObs(this, phase, sys_c)
            %  get syntetic observation for code or phase
            obs_type = {'code', 'phase'};
            this.log.addMessage(this.log.indent(sprintf('Synthesising %s observations', obs_type{phase + 1})));
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
            go_id = this.go_id(idx_obs);
            u_sat = unique(go_id);
            %this.updateAllAvailIndex();
            if ~this.isMultiFreq()
                this.updateErrIono();
            end
            if sum(this.xyz) ~=0
                this.updateErrTropo();
            end
            % for each unique go_id
            for i = u_sat'
                
                range = this.getSyntObs( i);
                
                sat_idx = find(go_id == i);
                % for all the obs with the same go_id
                for j = sat_idx'
                    o = (idx_obs == idx_obs(j));
                    c_obs_idx = idx_obs(j); % index of the observation we are currently processing
                    ep_idx = this.obs(c_obs_idx, :) ~= 0;
                    synt_pr_obs(o, ep_idx) = range(ep_idx);
                end
                
            end
        end
        
        function synt_ph = getSyntPhases(this, sys_c)
            % get current value of syntetic phase, in case not present update it
            if isempty(this.synt_ph)
                this.updateSyntPhases();
            end
            synt_ph = this.synt_ph;
            synt_ph(this.sat.outlier_idx_ph) = nan;
            if nargin == 2
                synt_ph = synt_ph(:, this.system(this.obs_code(:, 1) == 'L') == sys_c');
            end
        end                        
        
        function synt_pr_obs = getSyntPrObs(this, sys_c)
            if nargin < 2
                sys_c = this.cc.sys_c;
            end
            synt_pr_obs = zero2nan(this.getSyntCurObs(false, sys_c)');
        end
        
        function [synt_obs, xs_loc] = getSyntTwin(this, obs_set)
            %  get the syntethic twin for the observations
            % contained in obs_set
            % WARNING: the time of the observation set and the time of
            % receiver has to be the same
            synt_obs = zeros(size(obs_set.obs));
            xs_loc   = zeros(size(obs_set.obs,1), size(obs_set.obs,2),3);
            idx_ep_obs = obs_set.getTimeIdx(this.time.first, this.getRate);
            for i = 1 : size(synt_obs,2)
                go_id = obs_set.go_id(i);
                [range, xs_loc_t] = this.getSyntObs( go_id);
                idx_obs = obs_set.obs(:, i) ~= 0;
                idx_obs_r = idx_ep_obs(idx_obs); % <- to which epoch in the receiver the observation of the satellites in obesrvation set corresponds?
                idx_obs_r_l = false(size(range)); % get the logical equivalent
                idx_obs_r_l(idx_obs_r) = true;
                range_idx = range ~= 0;
                xs_idx = idx_obs_r_l(range_idx);
                synt_obs(idx_obs, i) = range(idx_obs_r);
                xs_loc(idx_obs, i, :) = permute(xs_loc_t(xs_idx, :),[1 3 2]);
            end
        end
        
        function id_sync = getIdSync(this)
            % MultiRec: works on an array of receivers
            % SYNTAX
            %  id_sync = this.getIdSync()
            id_sync = [];
            
            offset = 0;
            n_rec = numel(this);
            for r = 1 : n_rec
                id_sync = [id_sync; this(r).id_sync + offset]; %#ok<AGROW>
                offset = offset + this.length();
            end
        end
    end
    
    % ==================================================================================================================================================
    %% METHODS SETTERS
    % ==================================================================================================================================================
    methods
        function setStatic(this)
            % Set the internal status of the Receiver as static
            % SYNTAX
            %   this.setStatic()
            this.parent.static = true;
        end
        
        function setDynamic(this)
            % Set the internal status of the Receiver as dynamic
            % SYNTAX
            %   this.setDynamic()
            this.parent.static = false;
        end
        
        function setDoppler(this, dop, wl, id_dop)
            % set the snr observations
            % SYNTAX [pr, id_pr] = this.setDoppler(<sys_c>)
            % SEE ALSO:  getDoppler
            dop = bsxfun(@rdivide, zero2nan(dop'), wl);
            this.obs(id_dop, :) = nan2zero(dop');
        end
        
        function setSNR(this, snr, id_snr)
            % set the snr observations
            % SYNTAX [pr, id_pr] = this.setSNR(<sys_c>)
            % SEE ALSO:  getSNR
            this.obs(id_snr, :) = nan2zero(snr');
        end
        
        function setPhases(this, ph, wl, id_ph)
            % set the phases observations in meter (not cycles)
            % SYNTAX setPhases(this, ph, wl, id_ph)
            % SEE ALSO: getPhases getPseudoRanges setPseudoRanges
            
            ph = bsxfun(@rdivide, zero2nan(ph'), wl);
            this.obs(id_ph, :) = nan2zero(ph);
        end
        
        function resetSynthPhases(this)
            % recompute synhtesised phases
            %
            % SYNTAX 
            %   this.resetSynthPhases()
            
            this.synt_ph = [];
            this.updateSyntPhases();
        end
        
        function injectObs(this, obs, wl, f_id, obs_code, go_id)
            % Injecting observations into Receiver
            % This routine have been written for Core_SEID
            %
            % SINTEX:
            %   this.injectObs(obs, wl, f_id, obs_code, go_id);
            
            if size(obs, 2) ~= size(this.obs, 2)
                this.log.addError(sprintf('Observation injection not possible, input contains %d epochs while obs contains %d epochs', size(obs, 2), size(this.obs, 2)));
            else
                go_id = go_id(:);
                f_id = f_id(:);
                wl = wl(:);
                assert(size(obs, 1) == numel(go_id), 'Observation injection input "go_id" parameters size error');
                
                this.obs = [this.obs; obs];
                this.go_id = [this.go_id; go_id];
                this.active_ids = [this.active_ids; true(size(go_id))];
                if numel(f_id) == 1
                    this.f_id = [this.f_id; f_id * ones(size(go_id))];
                else
                    assert(size(obs, 1) == numel(f_id), 'Observation injection input "f_id" parameters size error');
                    this.f_id = [this.f_id; f_id];
                end
                if numel(wl) == 1
                    this.wl = [this.wl; wl * ones(size(go_id))];
                else
                    assert(size(obs, 1) == numel(wl), 'Observation injection input "wl" parameters size error');
                    this.wl = [this.f_id; wl];
                end
                
                [sys, prn] = this.getSysPrn(go_id);
                this.prn = [this.prn; prn];
                this.system = [this.system sys'];
                
                if size(obs_code,1) == 1
                    this.obs_code = [this.obs_code; repmat(obs_code, size(go_id, 1), 1)];
                else
                    assert(size(obs, 2) == numel(wl), 'Observation injection input "wl" parameters size error');
                    this.obs_code = [this.obs_code; obs_code];
                end
                this.resetSynthPhases();
            end
        end
        
        function setPseudoRanges(this, pr, id_pr)
            % set the pseudo ranges observations in meter (not cycles)
            % SYNTAX [pr, id_pr] = this.getPseudoRanges(<sys_c>)
            % SEE ALSO: getPhases setPhases getPseudoRanges
            this.obs(id_pr, :) = nan2zero(pr');
        end
        
        function setXYZ(this, xyz)
            %description set xyz and upadte godetic coordinates
            this.xyz = xyz;
            this.updateCoordinates();
        end
        
        function setActiveSys(this, sys_list)
            % set the active systems to bbe used for the computations
            % SYNTAX this.setActiveSys(sys_list)
            
            % Select only the systems still present in the file
            this.cc.setActive(sys_list);
        end 
    end
    
     % ==================================================================================================================================================
    %% METHODS UPDATERS
    % ==================================================================================================================================================
    methods
        function updateCoordinates(this)
            % upadte lat lon e ortometric height
            [this.lat, this.lon, this.h_ellips, this.h_ortho] = this.getMedianPosGeodetic();
        end
        
        function updateSyntPhases(this, sys_c)
            %DESCRIPTION update the content of the syntetic phase
            if nargin < 2
                sys_c = this.cc.sys_c;
            end
            synt_ph_obs = zero2nan(this.getSyntCurObs( true, sys_c)');
            synt_ph_obs(this.sat.outlier_idx_ph) = nan;
            this.synt_ph = synt_ph_obs;
        end
        
        function updateStatus(this)
            % Compute the other useful status array of the receiver object
            %
            % SYNTAX
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
            %for e = 1 : this.length()
            %    this.n_spe(e) = numel(unique(this.go_id(this.obs(:,e) ~= 0)));
            %end
            
            this.active_ids = any(this.obs, 2);
            
            %this.sat.avail_index = false(this.time.length, this.parent.cc.getMaxNumSat());
        end
        
        function updateRinObsCode(this)
            % update the RINEX observation codes contained into the receiver
            % SYNTAX  this.updateRinObsCode()
            this.rin_obs_code = struct('G',[],'R',[],'E',[],'J',[],'C',[],'I',[],'S',[]);
            sys_c = this.getActiveSys();
            for ss = sys_c
                this.rin_obs_code.(ss) = serialize(this.getAvailableCode(ss)')';
            end
        end
        
        function updateTOT(this, obs, go_id)
            % SYNTAX
            %   this.updateTOT(time_rx, dt);
            %
            % INPUT
            %
            % OUTPUT
            % 
            %   Update the signal time of travel based on range observations.
            % NOTE: to have satellite orbits at, 1mm sysncrnonization time of
            % travel has to be know with a precision of ~100m
            %    0.001m * (3 km/s  +         2 km/s)   /       300000 km/s   ~ 100m
            %              ^                 ^                    ^
            %
            %         (sat speed) (earth rot at sat orbit) (speed of light)
            if isempty(this.sat.tot)
                this.sat.tot = zeros(size(this.sat.avail_index));
            end
            if isempty(this.dt_ip) || ~any(this.dt_ip)
                this.sat.tot(:, go_id) =   nan2zero(zero2nan(obs)' / Global_Configuration.V_LIGHT + this.dt(:, 1));  %<---- check dt with all the new dts field
            else
                this.sat.tot(:, go_id) =   nan2zero(zero2nan(obs)' / Global_Configuration.V_LIGHT);  %<---- check dt with all the new dts field
            end
        end
        
        function updateAllTOT(this, synt_based)
            % upate time of travel for all satellites
            % for each receiver
            go_id = unique(this.go_id);
            for i = serialize(unique(this.go_id(this.obs_code(:,1) == 'C')))' % iterating only on sats with code observation
                
                sat_idx = (this.go_id == i) & (this.obs_code(:,1) == 'C');
                if nargin == 1 % obs based
                    c_obs = this.obs(sat_idx,:);
                    c_l_obs = colFirstNonZero(c_obs); % all best obs one each line %% CONSIDER USING ONLY 1 FREQUENCY FOR mm CONSISTENCY
                elseif synt_based % use the synteetc observation to get Time of flight
                    c_l_obs = this.getSyntObs(i);
                end
                %update time of flight times
                this.updateTOT(c_l_obs, i); % update time of travel
            end
        end
        
        function initAvailIndex(this, ep_ok)
            % initialize the avaliability index 
             this.sat.avail_index = false(this.time.length, this.parent.cc.getMaxNumSat);
             this.updateAllAvailIndex();
             if nargin == 2
                 if islogical(ep_ok)
                     this.sat.avail_index(~ep_ok,:) = false;
                 else
                     this.sat.avail_index(setdiff(1 : size(this.sat.avail_index, 1), ep_ok), :) = false;
                 end
             end
        end
        
        function updateAvailIndex(this, obs, go_gps)
            %  update avaliabilty of measurement on staellite
            this.sat.avail_index(:, go_gps) = obs > 0;
        end
        
        function updateAllAvailIndex(this)
            %  update avaliabilty of measurement on all
            % satellite based on all code and phase

            for s = unique(this.go_id)'
                obs_idx = this.go_id == s & (this.obs_code(:,1) == 'C' | this.obs_code(:,1) == 'L');
                if sum(obs_idx) > 0
                    av_idx = colFirstNonZero(this.obs(obs_idx,:)) ~= 0 ;
                    this.sat.avail_index(:,s) = av_idx;
                end
            end
        end
        
        function setAvIdx2Visibility(this)
            % set the availability index to el >0
            this.sat.avail_index = true(this.time.length, this.parent.cc.getMaxNumSat);
            this.updateAzimuthElevation;
            this.sat.avail_index = this.sat.el > 0;
        end
    end
    
    % ==================================================================================================================================================
    %% METHODS CORRECTIONS
    % ==================================================================================================================================================
    %  FUNCTIONS TO COMPUTE APPLY AND REMOVE VARIOUS MODELED CORRECTIONS
    %  NOTE: Methods related to corrections that are applied to the observables and whose
    %  values is not stored separately (All except Iono and tropo) are
    %  structured as follow:
    %   - A : add or subtract the correction to the observations
    %   - applyA : a warpper of A to add the correction
    %   - remA : a wrapper of A to subtract the correction
    %   - computeA : does the numerical computation of A along the los
    %  where A is the name of the correction
    %  NOTE ON THE SIGN:
    %   Two type of corrections are presnts: real error on the range
    %   measuremnts and correction needed to bring the a changing position
    %   at the same coordinates. The first group will be subtracted from
    %   the observation the second one added. More specifically:
    %   "REAL" CORRECTIONS: (-)
    %   - Shapiro Delay
    %   - Phase wind up
    %   - PCV
    %   CHANGE OF COORDINATES (+)
    %   + PCO
    %   + Solid Earth tides
    %   + Ocean Loading
    %   + Pole tides
    %
    % ==================================================================================================================================================
    
    methods (Access = public)
        
        %--------------------------------------------------------
        % Time
        %--------------------------------------------------------
        
        function applyDtSat(this)
            if this.dts_delay_status == 0
                this.applyDtSatFlag(1);
                this.dts_delay_status = 1; %applied
            end
        end
        
        function remDtSat(this)
            if this.dts_delay_status == 1
                this.applyDtSatFlag(-1);
                this.dts_delay_status = 0; %applied
            end
        end
        
        function applyDtSatFlag(this,flag)
            % Apply clock satellite corrections for code and phase
            % IMPORTANT: if no clock is present delete the observation
            
            if isempty(this.active_ids)
                this.active_ids = false(size(this.obs,1),1);
            end
            this.initAvailIndex();
            % for each receiver
            [~, id] = unique(this.go_id);
            for i = id'
                go_id = this.go_id(i);
                sat_idx = (this.go_id == this.go_id(i)) & (this.obs_code(:,1) == 'C' | this.obs_code(:,1) == 'L');
                ep_idx = logical(sum(this.obs(sat_idx,:) ~= 0,1));
                dts_range = ( this.getDtS(go_id) + this.getRelClkCorr(go_id) ) * Global_Configuration.V_LIGHT;
                for o = find(sat_idx)'
                    obs_idx_l = this.obs(o,:) ~= 0;
                    obs_idx = find(obs_idx_l);
                    dts_idx = obs_idx_l(ep_idx);
                    if ~isempty(dts_range(dts_idx))
                        if this.obs_code(o,1) == 'C'
                            this.obs(o, obs_idx_l) = this.obs(o,obs_idx_l) + sign(flag) * dts_range(dts_idx)';
                        else
                            this.obs(o, obs_idx_l) = this.obs(o,obs_idx_l) + sign(flag) * dts_range(dts_idx)'./this.wl(o);
                        end
                    end
                    dts_range_2 = dts_range(dts_idx);
                    nan_idx = obs_idx(isnan(dts_range_2));
                    this.obs(o, nan_idx) = 0;
                end
            end
        end
        
        function smoothAndApplyDt(this, smoothing_win, is_pr_jumping, is_ph_jumping)
            % Smooth dt * c correction computed from init-positioning with a spline with base 3 * rate,
            % apply the smoothed dt to pseudo-ranges and phases
            % SYNTAX
            %   this.smoothAndApplyDt(smoothing_win)
            if nargin == 1
                smoothing_win = 3;
            end
            if any(smoothing_win)
                this.log.addMessage(this.log.indent('Smooth and apply the clock error of the receiver'));
            else
                this.log.addMessage(this.log.indent('Apply the clock error of the receiver'));
            end
            id_ko = this.dt == 0;
            lim = getOutliers(this.dt(:,1) ~= 0 & abs(Core_Pre_Processing.diffAndPred(this.dt(:,1),2)) < 1e-7);
            % lim = [lim(1) lim(end)];
            dt_w_border = [this.dt(1); zero2nan(this.dt(:,1)); this.dt(end,1)];
            dt_pr = simpleFill1D(dt_w_border, isnan(dt_w_border), 'pchip');
            dt_pr([1 end] ) = [];
            if smoothing_win(1) > 0
                for i = 1 : size(lim, 1)
                    if lim(i,2) - lim(i,1) > 5
                        dt_pr(lim(i,1) : lim(i,2)) = splinerMat([], dt_pr(lim(i,1) : lim(i,2)), smoothing_win(1));
                    end
                end
            end
            if length(smoothing_win) == 2
                dt_pr = splinerMat([], dt_pr, smoothing_win(2));
            end
            
            this.dt = simpleFill1D(zero2nan(dt_pr), id_ko, 'pchip');
            if (~is_ph_jumping)
                id_jump = abs(diff(dt_pr)) > 1e-4; %find clock resets
                ddt = simpleFill1D(diff(dt_pr), id_jump, 'linear');
                dt_ph = dt_pr(1) + [0; cumsum(ddt)];
            else
                dt_ph = dt_pr;
            end
            this.applyDtRec(dt_pr, dt_ph)
            %this.dt_pr = this.dt_pr + this.dt;
            %this.dt_ph = this.dt_ph + this.dt;
            this.dt(:)  = 0;%zeros(size(this.dt_pr));
        end
        
        function applyDtRec(this, dt_pr, dt_ph)
            % Apply dt * c correction to pseudo-ranges and phases
            % SYNTAX
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
        
        function shiftToNominal(this)
            %  translate receiver observations to nominal epochs
            nominal_time = getNominalTime(this);
            tt = nominal_time - this.time;
            this.timeShiftObs(tt)
        end
        
        function timeShiftObs(this, seconds)
            %  translate observations at different epoch based on linear modeling of the satellite
            % copute the sat postion at the current epoch
            
            XS = this.getXSLoc();
            
            % translate the time
            this.time.addSeconds(seconds);
            % compute the sat postion at epoch traslated
            XS_t = this.getXSLoc();
            % compute the correction
            range = sqrt(sum(XS.^2,3));
            range_t = sqrt(sum(XS_t.^2,3));
            
            
            
            
            d_range = range_t - range;
            % Correct phases
            obs_idx = this.obs_code(:,1) == 'L' ;
            this.obs(obs_idx,:) = nan2zero(zero2nan(this.obs(obs_idx,:)) + bsxfun(@rdivide, d_range(:,this.go_id(obs_idx))', this.wl(obs_idx)));
            % Correct pseudo-ranges
            obs_idx = this.obs_code(:,1) == 'C';
            this.obs(obs_idx,:) = nan2zero(zero2nan(this.obs(obs_idx,:)) + d_range(:,this.go_id(obs_idx))');
        end
        
        % ----- shift to nominal at transmission time
        function shiftObs(this, seconds, go_id)
            XS = this.getXSLoc(go_id);
            % translate the time
            this.time.addSeconds(seconds);
            % compute the sat postion at epoch traslated
            XS_t = this.getXSLoc(go_id);
            this.time.addSeconds(-seconds);
            range = sqrt(sum(XS.^2,2));
            range_t = sqrt(sum(XS_t.^2,2));
            idx_sat = this.sat.avail_index(:,go_id);
            % compuet earth rotation correction
            d_earth = zeros(size(XS,1),1);
            [XS_tx_r ,XS_tx] = this.getXSTxRot(go_id);
            XS_tx_r = XS_tx_r - repmat(this.xyz,sum(idx_sat),1);
            XS_tx = XS_tx - repmat(this.xyz,sum(idx_sat),1);
            range_e = sqrt(sum(XS_tx_r.^2,2));
            range_ne = sqrt(sum(XS_tx.^2,2));
            d_earth(idx_sat) = range_e - range_ne;
            d_range = nan2zero(range_t - range - d_earth);
            % Correct phases
            obs_idx = this.obs_code(:,1) == 'L' & this.go_id == go_id ;
            this.obs(obs_idx,:) = nan2zero(zero2nan(this.obs(obs_idx,:)) + bsxfun(@rdivide, repmat(d_range',sum(obs_idx),1), this.wl(obs_idx)));
            % Correct pseudo-ranges
            obs_idx = this.obs_code(:,1) == 'C'  & this.go_id == go_id;
            this.obs(obs_idx,:) = nan2zero(zero2nan(this.obs(obs_idx,:)) + repmat(d_range',sum(obs_idx),1));
            
        end
        
        function shiftToNominalTransmission(this)
            % tranlate the obseravtion accounting for each travel time
            % so all observations are trasmitted at nominal time
            
            %firstly shift to nominal
            this.shiftToNominal();
            % then shift each sta for its time of travel
            go_ids = unique(this.go_id);
            for g = go_ids'
                tot = this.sat.tot(:,g);
                this.shiftObs(tot, g);
                this.sat.tot(:,g) = 0;
            end
        end
        
        %--------------------------------------------------------
        % Group Delay
        %--------------------------------------------------------
        
        function applyGroupDelay(this)
            if this.group_delay_status == 0
                this.log.addMarkedMessage('Adding code group delays');
                this.groupDelay(1);
                this.group_delay_status = 1; %applied
            end
        end
        
        function remGroupDelay(this)
            if this.group_delay_status == 1
                this.log.addMarkedMessage('Removing code group delays');
                this.groupDelay(-1);
                this.group_delay_status = 0; %applied
            end
        end
        
        function groupDelay(this, sgn)
            % . apply group delay corrections for code and phase
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
        
        %--------------------------------------------------------
        % Earth rotation parameters
        %--------------------------------------------------------
        
        function [XS_r] = earthRotationCorrection(this, XS, go_id)
            % SYNTAX
            %   [XS_r] = this.earthRotationCorrection(XS)
            %
            % INPUT
            %   XS      = positions of satellites
            %   time_rx = receiver time
            %   cc      = Constellation Collector
            %   sat     = satellite
            % OUTPUT
            %   XS_r    = Satellite postions rotated by earth roattion occured
            %   during time of travel
            % 
            %   Rotate the satellites position by the earth rotation
            %   occured during time of travel of the signal
            
            %%% TBD -> consider the case XS and travel_time does not match
            XS_r = zeros(size(XS));
            
            idx = this.sat.avail_index(:,go_id) > 0;
            travel_time = this.sat.tot(idx,go_id);
            sys = this.getSysPrn(go_id);
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
                    Logger.getInstance().addWarning('Something went wrong in Receiver.earthRotationCorrection() \nUnrecognized Satellite system!\n');
                    omegae_dot = this.cc.gps.ORBITAL_P.OMEGAE_DOT;
            end
            omega_tau = omegae_dot * travel_time;
            xR  = [cos(omega_tau)    sin(omega_tau)];
            yR  = [-sin(omega_tau)    cos(omega_tau)];
            XS_r(:,1) = sum(xR .* XS(:,1:2),2); % X
            XS_r(:,2) = sum(yR .* XS(:,1:2),2); % Y
            XS_r(:,3) = XS(:,3); % Z
        end
        
        %--------------------------------------------------------
        % Azimuth and Elevation
        %--------------------------------------------------------
        
        function updateAzimuthElevation(this, go_id)
            % Upadte azimute elevation into.sat
            % SYNTAX
            %   this.updateAzimuthElevation(<sat>)
            
            if nargin < 2
                for i = unique(this.go_id)'
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
                if isempty(this.sat.avail_index)
                    this.sat.avail_index = true(this.length, this.parent.cc.getMaxNumSat);
                end
                av_idx = this.sat.avail_index(:, go_id) ~= 0;
                [this.sat.az(av_idx, go_id), this.sat.el(av_idx, go_id)] = this.computeAzimuthElevation(go_id);
            end
        end
        
        function [mf, el_points] = computeEmpMF(this, el_points, show_fig)
            % Compute an empirical hisotropic and homogeneous mapping function
            %
            % SYNTAX
            %   [mf, el_points] = this.computeEmpMF(<el_points>)
            %
            % EXAMPLE
            %   [mf, el_points] = this.computeEmpMF();
            
            if nargin < 3
                show_fig = false;
            end
            
            use_median = true;
            if use_median
                nppb = 1; % number of points per bin (1 degree)
            else
                nppb = 100; % number of points per bin (1 degree)
            end
            
            el = this.getEl;
            slant_td = bsxfun(@rdivide, this.getSlantTD, this.getZtd());            
            
            % keep only valid data
            el = el(~isnan(zero2nan(slant_td)));
            slant_td = slant_td(~isnan(zero2nan(slant_td)));
            [el, id] = sort(el); slant_td = slant_td(id);
            if (nargin < 2) || isempty(el_points)
                el_points = max(0, (min(el(:)))) : 0.1 : min(90, (max(el(:))));
            end
            
            % Limit the number of points to use for spline interpolation (speedup)
            % 100 points per bin (1 degree) choosen randomly
            bin_size = hist(el, 0.5 : 0.1 : 90);
            bin_offset = [0; cumsum(bin_size(:))];
            id_ok = [];
            for b = 1 : numel(bin_size)
                if (bin_size(b) > nppb)
                    if use_median
                        bin_slant = slant_td(serialize(bin_offset(b) + (1 : bin_size(b))));
                        [~, id_tmp] = min(abs(bin_slant - median(bin_slant)));
                        id_ok = [id_ok; bin_offset(b) + id_tmp(1)];
                    else
                        id_tmp = serialize(bin_offset(b) + randperm(bin_size(b)));
                        id_ok = [id_ok; sort(id_tmp(1 : nppb))];
                    end
                else
                    id_ok = [id_ok; serialize(bin_offset(b) + (1 : bin_size(b)))];
                end
            end
            if show_fig
                figure;
                plot(el(:), slant_td(:), '.'); hold on;
                plot(el(id_ok), slant_td(id_ok), 'o');
            end
            % The mapping function at zenith is 1 => to best fit with spline I prefer to estimate
            % it with terminal value = 0
            el = el(id_ok);
            slant_td = slant_td(id_ok) - 1;
            zeinth_fill = (max(el): 0.01 : 90)';
            el = [el; zeinth_fill];
            slant_td = [slant_td; zeros(size(zeinth_fill))];
            
            [~, ~, ~, mf] = splinerMat(el, slant_td, 1, 0, el_points); % Spline base 10 degrees            
            mf = mf + 1; % readding the remD bias
            if show_fig
                plot(el_points, mf, '.','MarkerSize', 10);
            end
        end        
        
        function [az, el] = computeAzimuthElevation(this, go_id)
            % Force computation of azimuth and elevation
            %
            % SYNTAX
            %   [az, el] = this.computeAzimuthElevation(go_id)
            if (nargin < 2) || isempty(go_id)
                go_id = unique(this.go_id);
            end
            XS = this.getXSTxRot(go_id);
            [az, el] = this.computeAzimuthElevationXS(XS);
        end          
        
        function [az, el] = computeAzimuthElevationXS(this, XS, XR)
            % SYNTAX
            %   [az, el] = this.computeAzimuthElevationXS(XS)
            %
            % INPUT
            % XS = positions of satellite [n_epoch x 1]
            % XR = positions of reciever [n_epoch x 1] (optional, non static
            % case)
            % OUTPUT
            % Az = Azimuths of satellite [n_epoch x 1]
            % El = Elevations of satellite [n_epoch x 1]
            % during time of travel
            % 
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
            XSR = XS - XR; %%% sats orbit with origin in receiver
            
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
        
        function updateErrTropo(this, go_id)
            %INPUT
            % sat : number of sat
            % flag: flag of the tropo model
            %DESCRIPTION update the tropospheric correction
            
            atmo = Atmosphere.getInstance();
            
            if isempty(this.sat.err_tropo)
                this.sat.err_tropo = zeros(size(this.sat.avail_index));
            end
            
            if nargin < 2 || strcmp(go_id, 'all')
                this.log.addMessage(this.log.indent('Updating tropospheric errors'))
                
                go_id = unique(this.go_id)';
            end
            this.sat.err_tropo(:, go_id) = 0;
            
            %%% compute lat lon
            [~, ~, h_full] = this.getMedianPosGeodetic();
            if abs(h_full) > 1e4
                this.log.addWarning('Height out of reasonable value for terrestrial postioning skipping tropo update')
                return
            end
           % get meteo data
           
           
           
           % upadte the ztd zwd
           this.updateAprTropo();
           
           zwd = nan2zero(this.getZwd);
           apr_zhd = this.getAprZhd;
           cotel = zero2nan(cotd(this.sat.el(this.id_sync, :)));
           cosaz = zero2nan(cosd(this.sat.az(this.id_sync, :)));
           sinaz = zero2nan(sind(this.sat.az(this.id_sync, :)));
           [gn ,ge] = this.getGradient();
           [mfh, mfw] = this.getSlantMF();
           if any(ge(:))
               for g = go_id                   
                   this.sat.err_tropo(this.id_sync,g) = mfh(:,g) .* apr_zhd + mfw(:,g).*zwd;
               end
           else
               for g = go_id  
                   this.sat.err_tropo(this.id_sync,g) = mfh(:,g) .* apr_zhd + mfw(:,g) .* zwd + gn .* mfw(:,g) .* cotel(:,g) .* cosaz(:,g) + ge .* mfw(:,g) .* cotel(:,g) .* sinaz(:,g);
               end
           end
        end
        
        function updateErrIono(this, go_id, iono_model_override)
            if isempty(this.sat.err_iono)
                this.sat.err_iono = zeros(size(this.sat.avail_index));
            end
            this.log.addMessage(this.log.indent('Updating ionospheric errors'))
            if nargin < 2
                
                go_id  = 1 : this.parent.cc.getMaxNumSat();
            end            
          
            if nargin < 3
                iono_model_override = this.state.getIonoModel;
            end
            this.updateCoordinates();       
            switch iono_model_override
                case 1 % no model
                    this.sat.err_iono(idx,go_id) = zeros(size(el));
                case 2 % Klobuchar model
                    if ~isempty(this.sat.cs.iono )
                        for s = go_id
                            idx = this.sat.avail_idx(:,s);
                            [week, sow] = time2weektow(this.time.getSubSet(idx).getGpsTime());  
                            this.sat.err_iono(idx,go_id) = Atmosphere.klobucharModel(this.lat, this.lon, this.sat.az(idx,s), this.sat.el(idx,s), sow, this.sat.cs.iono);
                        end
                    else
                        this.log.addWarning('No klobuchar parameter found, iono correction not computed',100);
                    end
                case 3 % IONEX
                    atm = Atmosphere.getInstance();
                    this.sat.err_iono = atm.getFOIdelayCoeff(this.lat, this.lon, this.sat.az, this.sat.el, this.h_ellips, this.time);
            end
        end

        function applyIonoModel(this)
            if this.iono_status == 0 
                this.log.addMarkedMessage('Applying Ionosphere model');
                this.ionoModel(1);
                this.iono_status = 1; % applied
            end
        end
        
        function remIonoModel(this)
            if this.iono_status == 1
                this.log.addMarkedMessage('Removing Ionosphere model');
                this.ionoModel(-1);
                this.iono_status = 0; %not applied
            end
        end
        
        function ionoModel(this, sign)
            % Apply the selected iono model
            %
            % SYNTAX
            %   this.applyIonoModel()
            ph_obs = this.obs_code(:,1) == 'L';
            pr_obs = this.obs_code(:,1) == 'C';
            for i = 1 : size(this.sat.err_iono,2)
                idx_sat = this.go_id == i & (ph_obs | pr_obs);
                if sum(idx_sat) >0
                    for s = find(idx_sat)'
                        iono_delay = this.sat.err_iono(:,i) * this.wl(s).^2;
                        if ph_obs(s) > 0
                            this.obs(s,:) = nan2zero(zero2nan(this.obs(s,:)) + sign * (iono_delay'));
                        else
                            this.obs(s,:) = nan2zero(zero2nan(this.obs(s,:)) - sign * (iono_delay'));
                        end
                    end
                end
            end
        end
        
        %--------------------------------------------------------
        % Solid earth tide
        % -------------------------------------------------------
        
        function solidEarthTide(this,sgn)
            %  add or subtract ocean loading from observations
            et_corr = this.computeSolidTideCorr();
            
            for s = 1 : this.parent.cc.getMaxNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        if  this.obs_code(o,1) == 'L'
                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn)* et_corr(o_idx,s)' ./ this.wl(o);
                        else
                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn)* et_corr(o_idx,s)';
                        end
                    end
                end
            end
        end
        
        function applySolidEarthTide(this)
            if this.et_delay_status == 0 && this.state.isSolidEarth
                this.log.addMarkedMessage('Applying Solid Earth Tide corrections');
                this.solidEarthTide(1);
                this.et_delay_status = 1; %applied
            end
        end
        
        function remSolidEarthTide(this)
            if this.et_delay_status == 1
                this.log.addMarkedMessage('Removing Solid Earth Tide corrections');
                this.solidEarthTide(-1);
                this.et_delay_status = 0; %not applied
            end
        end
        
        function solid_earth_corr = computeSolidTideCorr(this, sat)
            
            % SYNTAX
            %   [stidecorr] = this.getSolidTideCorr();
            %
            % INPUT
            %
            % OUTPUT
            %   stidecorr = solid Earth tide correction terms (along the satellite-receiver line-of-sight)
            %
            % 
            %   Computation of the solid Earth tide displacement terms.
            if nargin < 2
                sat  = 1 : this.parent.cc.getMaxNumSat();
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
            GE = Galileo_SS.ORBITAL_P.GM; %Earth
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
            [XS] = this.getXSLoc();
            for i  = 1 : length(sat)
                s = sat(i);
                sat_idx = this.sat.avail_index(:,s) ~= 0;
                XSs = permute(XS(sat_idx,s,:),[1 3 2]);
                LOSu = rowNormalize(XSs);
                solid_earth_corr(sat_idx,i) = sum(conj(r(sat_idx,:)).*LOSu,2);
            end
            
        end
        
        function updateSolidEarthCorr(this, sat)
            %DESCRIPTION upadte the correction related to solid earth
            % solid tides, ocean loading, pole tides.
            if isempty(this.sat.solid_earth_corr)
                this.log.addMessage(this.log.indent('Updating solid earth corrections'))
                this.sat.solid_earth_corr = zeros(size(this.sat.avail_index));
            end
            if nargin < 2
                this.sat.solid_earth_corr = this.computeSolidTideCorr();
                %                 for s = 1 : size(this.sat.avail_index,2)
                %                     this.updateSolidEarthCorr(s);
                %                 end
            else
                this.sat.solid_earth_corr(:,sat) = this.computeSolidTideCorr(sat);% + this.computeOceanLoading(sat) + this.getPoleTideCorr(sat);
            end
        end
        
        %--------------------------------------------------------
        % Ocean Loading
        % -------------------------------------------------------
        
        function oceanLoading(this,sgn)
            %  add or subtract ocean loading from observations
            ol_corr = this.computeOceanLoading();
            if isempty(ol_corr)
                this.log.addWarning('No ocean loading displacements matrix present for the receiver')
                return
            end
            
            for s = 1 : this.parent.cc.getMaxNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        if  this.obs_code(o,1) == 'L'
                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn)* ol_corr(o_idx,s)' ./ this.wl(o);
                        else
                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn)* ol_corr(o_idx,s)';
                        end
                    end
                end
            end
        end
        
        function applyOceanLoading(this)
            if this.ol_delay_status == 0 && this.state.isOceanLoading
                this.log.addMarkedMessage('Applying Ocean Loading corrections');
                this.oceanLoading(1);
                this.ol_delay_status = 1; %applied
            end
        end
        
        function remOceanLoading(this)
            if this.ol_delay_status == 1
                this.log.addMarkedMessage('Removing Ocean Loading corrections');
                this.oceanLoading(-1);
                this.ol_delay_status = 0; %not applied
            end
        end
        
        function [ocean_load_corr] = computeOceanLoading(this, sat) % WARNING: to be tested
            
            % SYNTAX
            %   [oceanloadcorr] = ocean_loading_correction(time, XR, XS);
            %
            % INPUT
            %
            % OUTPUT
            %   oceanloadcorr = ocean loading correction terms (along the satellite-receiver line-of-sight)
            %
            % 
            %   Computation of the ocean loading displacement terms.
            % NOTES:
            %  Inefficent to compute them separate by staellites always call
            %  as a block
            
            if nargin < 2
                sat = 1 : this.parent.cc.getMaxNumSat();
            end
            
            
            ocean_load_corr = zeros(this.time.length, length(sat));
            if isempty(this.ocean_load_disp)
                this.importOceanLoading();
            end
            if isempty(this.ocean_load_disp) || ((~isstruct(this.ocean_load_disp) && this.ocean_load_disp == -1) || this.ocean_load_disp.available == 0)
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
            corrENU(:,1) = -corr(:,2); %east
            corrENU(:,2) = -corr(:,3); %north
            corrENU(:,3) =  corr(:,1); %up
            
            % get XR
            XR = this.getXR();
            %displacement along the receiver-satellite line-of-sight
            XRcorr = local2globalPos(corrENU', XR')';
            corrXYZ = XRcorr - XR;
            %displacement along the receiver-satellite line-of-sight
            [XS] = this.getXSLoc();
            for i  = 1 : length(sat)
                s = sat(i);
                sat_idx = this.sat.avail_index(:,s);
                XSs = permute(XS(sat_idx,s,:),[1 3 2]);
                LOSu = rowNormalize(XSs);
                ocean_load_corr(sat_idx,i) = sum(conj(corrXYZ(sat_idx,:)).*LOSu,2);
            end
        end
        
        %--------------------------------------------------------
        % Pole Tide
        % -------------------------------------------------------
        
        function poleTide(this,sgn)
            %  add or subtract ocean loading from observations
            pt_corr = this.computePoleTideCorr();
            if isempty(pt_corr)
                this.log.addWarning('No ocean loading displacements matrix present for the receiver')
                return
            end
            
            for s = 1 : this.parent.cc.getMaxNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        if  this.obs_code(o,1) == 'L'
                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn)* pt_corr(o_idx,s)' ./ this.wl(o);
                        else
                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn)* pt_corr(o_idx,s)';
                        end
                    end
                end
            end
        end
        
        function applyPoleTide(this)
            if this.pt_delay_status == 0 && this.state.isPoleTide
                this.log.addMarkedMessage('Applying Pole Tide corrections');
                this.poleTide(1);
                this.pt_delay_status = 1; %applied
            end
        end
        
        function remPoleTide(this)
            if this.pt_delay_status == 1
                this.log.addMarkedMessage('Removing Pole Tide corrections');
                this.poleTide(-1);
                this.pt_delay_status = 0; %not applied
            end
        end
        
        function [pole_tide_corr] = computePoleTideCorr(this, sat)
            
            % SYNTAX
            %   [poletidecorr] = pole_tide_correction(time, XR, XS, SP3, phiC, lam);
            %
            % INPUT
            %
            % OUTPUT
            %   poletidecorr = pole tide correction terms (along the satellite-receiver line-of-sight)
            %
            % 
            %   Computation of the pole tide displacement terms.
            
            if nargin < 2
                sat = 1 : this.parent.cc.getMaxNumSat();
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
                [XS] = this.getXSLoc();
                for i  = 1 : length(sat)
                    s = sat(i);
                    sat_idx = this.sat.avail_index(:,s);
                    XSs = permute(XS(sat_idx,s,:),[1 3 2]);
                    LOSu = rowNormalize(XSs);
                    pole_tide_corr(sat_idx,i) = sum(conj(corrXYZ(sat_idx,:)).*LOSu,2);
                end
            end
            
        end
        
        %--------------------------------------------------------
        % High Order Ionospheric effect + bending
        % -------------------------------------------------------
        
        function HOI(this,sgn)
            %  add or subtract ocean loading from observations
            [hoi2, hoi3, bending] = this.computeHOI();
            
            for s = 1 : this.parent.cc.getMaxNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        wl = this.wl(o);
                        wl3 = wl^3;
                        wl4 = wl^4;
                        if  this.obs_code(o,1) == 'L'
                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)*( -1/2 *hoi2(o_idx,s)'*wl3 -1/3* hoi3(o_idx,s)'*wl4 +  0*bending(o_idx,s)'*wl4) ./ this.wl(o);
                        else
                            this.obs(o,o_idx) = this.obs(o,o_idx) - sign(sgn)*( hoi2(o_idx,s)'*wl3 + hoi3(o_idx,s)'*wl4 + 0*bending(o_idx,s)'*wl4);
                        end
                    end
                end
            end
        end
        
        function applyHOI(this)
            if this.hoi_delay_status == 0 && this.state.isHOI
                this.log.addMarkedMessage('Applying High Order Ionospheric Effect');
                this.HOI(1);
                this.hoi_delay_status = 1; %applied
            end
        end
        
        function remHOI(this)
            if this.hoi_delay_status == 1
                this.log.addMarkedMessage('Removing High Order Ionospheric Effect');
                this.HOI(-1);
                this.hoi_delay_status = 0; %not applied
            end
        end
        
        function [hoi_delay2_coeff, hoi_delay3_coeff, bending_coeff] =  computeHOI(this)
            
            % SYNTAX
            %
            % INPUT
            %
            % OUTPUT
            %
            %
            % 
            %   Computation of thigh order ionospheric effect
            
            this.updateCoordinates();
            atmo = Atmosphere.getInstance();
            [hoi_delay2_coeff, hoi_delay3_coeff, bending_coeff] = atmo.getHOIdelayCoeff(this.lat,this.lon, this.sat.az,this.sat.el,this.h_ellips,this.time);
            
        end
        
        %--------------------------------------------------------
        % Atmospheric loading
        % -------------------------------------------------------
        
        function atmLoad(this,sgn)
           %  add or subtract ocean loading from observations
            al_corr = this.computeAtmLoading(); 
            if isempty(al_corr)
                this.log.addWarning('No ocean loading displacements matrix present for the receiver')
                return
            end
            
            for s = 1 : this.parent.cc.getMaxNumSat()
                obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                obs_idx = obs_idx & this.go_id == s;
                if sum(obs_idx) > 0
                    for o = find(obs_idx)'
                        o_idx = this.obs(o, :) ~=0; %find where apply corrections
                        if  this.obs_code(o,1) == 'L'
                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn)* al_corr(o_idx,s)' ./ this.wl(o);
                        else
                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn)* al_corr(o_idx,s)';
                        end
                    end
                end
            end
        end
        
        function applyAtmLoad(this)
            if this.atm_load_delay_status == 0  && this.state.isAtmLoading
                this.log.addMarkedMessage('Applying atmospheric loading effect');
                this.atmLoad(1);
                this.atm_load_delay_status = 1; %applied
            end
        end
        
        function remAtmLoad(this)
            if this.atm_load_delay_status == 1
                this.log.addMarkedMessage('Removing atmospheric loading effect');
                this.atmLoad(-1);
                this.atm_load_delay_status = 0; %not applied
            end
        end
        
        function [al_corr] =  computeAtmLoading(this, sat)
            
            % SYNTAX
            %
            % INPUT
            %
            % OUTPUT
            %
            %
            % 
            %   Computation of atmopsheric loading
            
            if nargin < 2
                sat = 1 : this.parent.cc.getMaxNumSat();
            end
            this.updateCoordinates();
            atmo = Atmosphere.getInstance();
            dsa = this.time.first.getCopy();
            dso = this.time.getSubSet(this.time.length).getCopy();
            dso.addSeconds(6*3600);
            fname = this.state.getNTAtmLoadFileName( dsa, dso);
            for i = 1 : length(fname)
                atmo.importAtmLoadCoeffFile(fname{i});
            end
            [corrXYZ] = atmo.getAtmLoadCorr(this.lat,this.lon,this.h_ellips, this.time);
            [XS] = this.getXSLoc();
            for i  = 1 : length(sat)
                s = sat(i);
                sat_idx = this.sat.avail_index(:,s);
                XSs = permute(XS(sat_idx,s,:),[1 3 2]);
                LOSu = rowNormalize(XSs);
                al_corr(sat_idx,i) = sum(conj(corrXYZ(sat_idx,:)).*LOSu,2);
            end
            
        end
        
        %--------------------------------------------------------
        % Phase Wind up
        % -------------------------------------------------------
        
        function phaseWindUpCorr(this,sgn)
            %  add or subtract ocean loading from observations
            ph_wind_up = this.computePhaseWindUp();
            
            for s = 1 : this.parent.cc.getMaxNumSat()
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
            if this.pw_delay_status == 0 && this.state.isPhaseWind
                this.log.addMarkedMessage('Applying Phase Wind Up corrections');
                this.phaseWindUpCorr(1);
                this.pw_delay_status = 1; %applied
            end
        end
        
        function remPhaseWindUpCorr(this)
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
            
            
            sat = 1: this.parent.cc.getMaxNumSat();
            
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
            %  add or subtract shapiro delay from observations
            for s = 1 : this.parent.cc.getMaxNumSat()
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
            if this.sh_delay_status == 0 && this.state.isShapiro
                this.log.addMarkedMessage('Applying Shapiro delay corrections');
                this.shDelay(1);
                this.sh_delay_status = 1; %applied
            end
        end
        
        function remShDelay(this)
            if this.sh_delay_status == 1
                this.log.addMarkedMessage('Removing Shapiro delay corrections');
                this.shDelay(-1);
                this.sh_delay_status = 0; %not applied
            end
        end
        
        function [sh_delay] = computeShapirodelay(this, sat)
            % SYNTAX
            %   [corr, distSR_corr] = this.getRelDistance(XS, XR);
            %
            % INPUT
            % XS = positions of satellite [n_epoch x 1]
            % XR = positions of reciever [n_epoch x 1] (optional, non static
            % case)
            %
            % OUTPUT
            %   corr = relativistic range error correction term (Shapiro delay)
            %   dist = dist
            % 
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
            
            %corr = 2*GM/(Global_Configuration.V_LIGHT^2) * log((distR + distS + distSR)./(distR + distS - distSR)); %#ok<CPROPLC>
            
            sh_delay = 2*GM/(Global_Configuration.V_LIGHT^2) * log((distR + distS + distSR)./(distR + distS - distSR));
            
        end
        
        %--------------------------------------------------------
        % PCV
        % -------------------------------------------------------
        
        function applyremPCV(this, sgn)
            %  correct measurement for PCV both of receiver
            % antenna and satellite antenna
            if ~isempty(this.pcv) || ~isempty(this.sat.cs.ant_pcv)
                % this.updateAllAvailIndex(); % not needed?
                % getting sat - receiver vector for each epoch
                XR_sat = - this.getXSLoc();
                
                % Receiver PCV correction
                if ~isempty(this.pcv) && ~isempty(this.pcv.name) && this.state.isRecPCV()
                    f_code_history = []; % save f_code checked to print only one time the warning message
                    for s = unique(this.go_id)'
                        sat_idx = this.sat.avail_index(:, s);
                        el = this.sat.el(sat_idx, s);
                        az = this.sat.az(sat_idx, s);
                        neu_los = [cosd(az).*cosd(el) sind(az).*cosd(el) sind(el)];
                        obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                        obs_idx = obs_idx & this.go_id == s;
                        if sum(obs_idx) > 0 && (this.pcv.n_frequency > 0)
                            freqs = unique(str2num(this.obs_code(obs_idx,2)));
                            for f = freqs'
                                obs_idx_f = obs_idx & this.obs_code(:,2) == num2str(f);
                                sys = this.system(obs_idx_f);
                                f_code = [sys(1) sprintf('%02d',f)];
                                
                                pco_idx = idxCharLines(this.pcv.frequency_name, f_code);
                                if sum(pco_idx)
                                    pco = this.pcv.offset(:, :, pco_idx)';
                                    pco_delays = neu_los * pco;
                                    pcv_delays = pco_delays - this.getPCV(pco_idx, el, az);
                                    for o = find(obs_idx_f)'
                                        pcv_idx = this.obs(o, this.sat.avail_index(:, s)) ~= 0; % find which correction to apply
                                        o_idx = this.obs(o, :) ~= 0; % find where apply corrections
                                        if  this.obs_code(o, 1) == 'L'
                                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn) * pcv_delays(pcv_idx)' ./ this.wl(o);
                                        else
                                            this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn) * pcv_delays(pcv_idx)';
                                        end
                                    end
                                else
                                    if isempty(f_code_history) || ~sum(idxCharLines(f_code_history, f_code))
                                        this.log.addMessage(this.log.indent(sprintf('No corrections found for antenna model %s on frequency %s',this.parent.ant_type, f_code)));
                                        f_code_history = [f_code_history;f_code];
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Satellite PCV correction
                if ~isempty(this.sat.cs.ant_pcv)                    
                    % getting satellite reference frame for each epoch
                    [x, y, z] = this.sat.cs.getSatFixFrame(this.time);
                    sat_ok = unique(this.go_id)';
                    % transform into satellite reference system
                    XR_sat(:, sat_ok, :) = cat(3,sum(XR_sat(:, sat_ok, :) .* x(:, sat_ok, :),3), sum(XR_sat(:, sat_ok, :) .* y(:, sat_ok, :), 3), sum(XR_sat(:, sat_ok, :) .* z(:, sat_ok, :), 3));
                    
                    % getting az and el
                    distances = sqrt(sum(XR_sat.^2,3));
                    XR_sat_norm = XR_sat ./ repmat(distances, 1, 1, 3);
                    
                    az = atan2(XR_sat_norm(:, :, 2), XR_sat_norm(:, :, 1)); % here azimuth is intended as angle from x axis
                    az(az<0) = az(az<0) + 2*pi;
                    el = atan2(XR_sat_norm(:, :, 3), sqrt(sum(XR_sat_norm(:, :, 1:2).^2, 3)));
                    
                    % getting pscv and applying it to the obs
                    for s = unique(this.go_id)'
                        obs_idx = this.obs_code(:,1) == 'C' |  this.obs_code(:,1) == 'L';
                        obs_idx = obs_idx & this.go_id == s;
                        if sum(obs_idx) > 0
                            freqs = str2num(unique(this.obs_code(obs_idx,2)));
                            for f = freqs'
                                obs_idx_f = obs_idx & this.obs_code(:,2) == num2str(f);
                                az_idx = ~isnan(az(:,s));
                                az_tmp = az(az_idx,s) / pi * 180;
                                el_idx = ~isnan(el(:,s));
                                el_tmp = el(el_idx,s) / pi * 180;
                                ant_id = this.getAntennaId(s);
                                pcv_delays = this.sat.cs.getPCV( f, ant_id, el_tmp, az_tmp);
                                for o = find(obs_idx_f)'
                                    pcv_idx = this.obs(o, az_idx) ~= 0; %find which correction to apply
                                    o_idx = this.obs(o, :) ~=0 & az_idx'; %find where apply corrections
                                    if  this.obs_code(o,1) == 'L'
                                        this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn) * pcv_delays(pcv_idx)' ./ this.wl(o); % is it a plus
                                    else
                                        this.obs(o,o_idx) = this.obs(o,o_idx) + sign(sgn) * pcv_delays(pcv_idx)';
                                    end
                                end
                            end
                        end
                    end
                    
                end
            end
        end
        
        function pcv_delay = getPCV(this, idx, el, az, method)
            % get the pcv correction for a given satellite and a given
            % azimuth and elevations using linear or bilinear interpolation
            
            if (nargin < 5)
                method = 'minterp';
            end
            pcv_delay = zeros(size(el));
            
            pcv = this.pcv;
            
            % transform el in zen
            zen = 90 - el;
            % get el idx
            zen_pcv = pcv.tablePCV_zen;
            
            min_zen = zen_pcv(1);
            max_zen = zen_pcv(end);
            d_zen = (max_zen - min_zen)/(length(zen_pcv)-1);
            zen_float = (zen - min_zen)/d_zen + 1;
            if nargin < 4 || isempty(pcv.tablePCV_azi) % no azimuth change
                pcv_val = pcv.tableNOAZI(:,:,idx); % extract the right frequency                
                %pcv_delay = d_f_r_el .* pcv_val(zen_idx)' + (1 - d_f_r_el) .* pcv_val(zen_idx + 1)';
                
                % Use polynomial interpolation to smooth PCV
                pcv_delay = Core_Utils.interp1LS(1 : numel(pcv_val), pcv_val, min(8,numel(pcv_val)), zen_float);                
            else                                                
                % find azimuth indexes
                az_pcv = pcv.tablePCV_azi;
                min_az = az_pcv(1);
                max_az = az_pcv(end);
                d_az = (max_az - min_az)/(length(az_pcv)-1);
                az_float = (mod(az, 360) - min_az)/d_az;
                az_idx = min(max(floor((mod(az,360) - min_az)/d_az) + 1, 1),length(az_pcv) - 1);
                d_f_r_az = min(max(mod(az, 360) - (az_idx-1)*d_az, 0)/d_az, 1);
                
                
                if strcmp(method, 'polyLS')
                    % Polynomial interpolation in elevation and griddedInterpolant in az
                    % filter PCV values in elevetion - slower
                    %%
                    pcv_val = pcv.tablePCV(:,:,idx); % extract the right frequency
                    step = 0.25;
                    zen_interp = (((floor(min(zen)/step)*step) : step : (ceil(max(zen)/step)*step))' - min_zen)/d_zen + 1;

                    % Interpolate with a 7th degree polinomial in elevation
                    pcv_val = Core_Utils.interp1LS(1 : numel(pcv_val), pcv_val', min(8,size(pcv_val, 2)), zen_interp);
                    
                    n_az = size(pcv_val, 2);
                    b_dx = (n_az - floor(n_az / 3)) : (n_az - 1); % dx border
                    b_sx = 2 : (n_az - floor(n_az / 3));          % sx border
                    % replicate border for azimuth periodicity
                    pcv_val = [pcv_val(:, b_dx), pcv_val, pcv_val(:, b_sx)];
                    az_val = [az_pcv(b_dx) - 360, az_pcv, az_pcv(b_sx) + 360];
                    
                    % Using scattered interpolant for azimuth to have a smooth interpolation
                    [zen_m, az_m] = ndgrid(zen_interp, az_val);
                    fun = griddedInterpolant(zen_m, az_m, pcv_val);
                    pcv_delay = fun(zen_float, mod(az, 360));
                elseif strcmp(method, 'minterp')
                    % Interpolation with griddedInterpolant
                    % leinear interpolation - faster
                    pcv_val = pcv.tablePCV(:,:,idx)'; % extract the right frequency
                    
                    n_az = size(pcv_val, 2);
                    b_dx = (n_az - floor(n_az / 3)) : (n_az - 1); % dx border
                    b_sx = 2 : (n_az - floor(n_az / 3));          % sx border
                    % replicate border for azimuth periodicity
                    pcv_val = [pcv_val(:, b_dx), pcv_val, pcv_val(:, b_sx)];
                    az_val = [az_pcv(b_dx) - 360, az_pcv, az_pcv(b_sx) + 360];
                    
                    % Using scattered interpolant for azimuth to have a smooth interpolation
                    [zen_m, az_m] = ndgrid(zen_pcv, az_val);
                    fun = griddedInterpolant(zen_m, az_m, pcv_val, 'linear');
                    pcv_delay = fun(zen, az);
                elseif strcmp(method, 'lin')
                    
                    % linear interpolation between consecutive values 
                    % bad performance
                    %%
                    zen_idx = min(max(floor((zen - min_zen)/d_zen) + 1 , 1),length(zen_pcv) - 1);
                    d_f_r_el = min(max(zen_idx*d_zen - zen, 0)/ d_zen, 1) ;
                    pcv_val = pcv.tablePCV(:,:,idx); % extract the right frequency
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
        end
        
        function applyPCV(this)
            if (this.pcv_delay_status == 0) && this.state.isRecPCV
                this.log.addMarkedMessage('Applying PCV corrections');
                this.applyremPCV(1);
                this.pcv_delay_status = 1; % applied
            end
        end
        
        function remPCV(this)
            if this.pcv_delay_status == 1
                this.log.addMarkedMessage('Removing PCV corrections');
                this.applyremPCV(-1);
                this.pcv_delay_status = 0; % applied
            end
        end
        
        %--------------------------------------------------------
        % All
        %--------------------------------------------------------
        
        function remAllCorrections(this)
            % rem all the corrections
            % SYNTAX
            %   this.remAllCorrections();
            this.remDtSat();
            this.remGroupDelay();
            this.remPCV();
            this.remPoleTide();
            this.remPhaseWindUpCorr();
            this.remSolidEarthTide();
            this.remShDelay();
            this.remOceanLoading();
            this.remHOI();
            this.remAtmLoad();
        end
        
        function applyAllCorrections(this)
            % Apply all the corrections
            % SYNTAX
            %   this.applyAllCorrections();
            this.applyDtSat();
            this.applyGroupDelay();
            this.applyPCV();
            this.applyPoleTide();
            this.applyPhaseWindUpCorr();
            this.applySolidEarthTide();
            this.applyShDelay();
            this.applyOceanLoading();
            this.applyHOI();
            this.applyAtmload();
        end
    end
    
       
    % ==================================================================================================================================================
    %% METHODS PROCESSING
    % ==================================================================================================================================================
    
    methods
        function [is_pr_jumping, is_ph_jumping] = correctTimeDesync(this)
            %   Correction of jumps in code and phase due to dtR and time de-sync
            % SYNTAX
            %   this.correctTimeDesync()
            this.log.addMarkedMessage('Correct for time desync');
            disable_dt_correction = true; % Skip to speed-up processing
            
            % computing nominal_time
            nominal_time_zero = floor(this.time.first.getMatlabTime() * 24)/24; % start of day 
            rinex_time = this.time.getRefTime(nominal_time_zero); % seconds from start of day
            nominal_time = round(rinex_time / this.time.getRate) * this.time.getRate;
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
            
            dt_ph = zeros(this.time.length, 1);
            dt_pr = zeros(this.time.length, 1);
            
            [ph_dj, dt_ph_dj, is_ph_jumping] = Core_Pre_Processing.remDtJumps(ph);
            [pr_dj, dt_pr_dj, is_pr_jumping] = Core_Pre_Processing.remDtJumps(pr);
            % apply desync
            if ~disable_dt_correction
                if any(time_desync)
                    if (numel(dt_ph_dj) > 1 && numel(dt_pr_dj) > 1)
                        id_ko = flagExpand(sum(~isnan(ph), 2) == 0, 1);
                        tmp = diff(dt_pr_dj); tmp(~id_ko) = 0; tmp = cumsum(tmp);
                        dt_ph_dj = dt_ph_dj + tmp; % use jmp estimation from pseudo-ranges
                    end
                    
                    ddt_pr = Core_Pre_Processing.diffAndPred(dt_pr_dj);
                    
                    if (max(abs(time_desync)) > 1e-4)
                        % time_desync is a introduced by the receiver to maintain the drift of the clock into a certain range
                        time_desync = simpleFill1D(time_desync, all(isnan(pr), 2), 'nearest');
                        ddt = [0; diff(time_desync)];
                        ddrifting_pr = ddt - ddt_pr;
%                         drifting_pr = cumsum(ddrifting_pr);
                        drifting_pr = time_desync - dt_pr_dj;
                        
                        % Linear interpolation of ddrifting
                        jmp_reset = find(abs(ddt_pr) > 1e-7); % points where the clock is reset
                        lim = getOutliers(any(~isnan(ph),2));
                        lim = [lim(1); lim(end)];
                        jmp_fit = setdiff([find(abs(ddrifting_pr) > 1e-7); lim(:)], jmp_reset); % points where desync interpolate the clock
                        jmp_fit(jmp_fit==numel(drifting_pr)) = [];
                        jmp_fit(jmp_fit==1) = [];
                        d_points_pr = [drifting_pr(jmp_reset); drifting_pr(jmp_fit) - ddrifting_pr(jmp_fit)/2];
                        jmp = [jmp_reset; jmp_fit];
                        if numel(d_points_pr) < 3
                            drifting_pr = zeros(size(drifting_pr));
                        else
                            %lim = interp1(jmp, d_points_pr, [1 numel(drifting_pr)]', 'linear', 'extrap');
                            %drifting_pr = interp1([1; jmp; numel(drifting_pr)], [lim(1); d_points_pr; lim(2)], 'pchip');
                            drifting_pr = interp1(jmp, d_points_pr, (1 : numel(drifting_pr))', 'pchip');
                        end
                        
                        dt_ph_dj = iif(is_pr_jumping == is_ph_jumping, dt_pr_dj, dt_ph_dj);
                        dt_ph = drifting_pr + dt_ph_dj;
                        %dt_ph = drifting_pr * double(numel(dt_ph_dj) <= 1) + detrend(drifting_pr) * double(numel(dt_ph_dj) > 1) + dt_ph_dj;
                        dt_pr = drifting_pr + dt_pr_dj;
                        if ~isempty(jmp)
                            t_offset = round(mean(dt_pr(jmp) - time_desync(jmp) + ddrifting_pr(jmp)/2) * 1e7) * 1e-7;
                            dt_ph = dt_ph - t_offset;
                            dt_pr = dt_pr - t_offset;
                        end
                    else
                        dt_ph = dt_ph_dj;
                        dt_pr = dt_pr_dj;
                    end
                    
                    ph = bsxfun(@minus, ph, dt_ph .* 299792458);
                    pr = bsxfun(@minus, pr, dt_pr .* 299792458);
                else
                    if (is_ph_jumping ~= is_pr_jumping)
                        
                        % epochs with no phases cannot estimate dt jumps
                        if (numel(dt_ph_dj) > 1 && numel(dt_pr_dj) > 1)
                            id_ko = flagExpand(sum(~isnan(ph), 2) == 0, 1);
                            tmp = diff(dt_pr_dj); tmp(~id_ko) = 0; tmp = cumsum(tmp);
                            dt_ph_dj = dt_ph_dj + tmp; % use jmp estimation from pseudo-ranges
                        end
                        
                        ddt_pr = Core_Pre_Processing.diffAndPred(dt_pr_dj);
                        jmp_reset = find(abs(ddt_pr) > 1e-7); % points where the clock is reset
                        jmp_reset(jmp_reset==numel(dt_pr_dj)) = [];
                        if numel(jmp_reset) > 2
                            drifting_pr = interp1(jmp_reset, dt_pr_dj(jmp_reset), (1 : numel(dt_pr_dj))', 'pchip');
                        else
                            drifting_pr = 0;
                        end
                        
                        dt_pr = dt_pr_dj - drifting_pr;
                        t_offset = mean(dt_pr);
                        dt_ph = dt_ph_dj - drifting_pr - t_offset;
                        dt_pr = dt_pr - t_offset;
                        
                        ph = bsxfun(@minus, ph, dt_ph .* 299792458);
                        pr = bsxfun(@minus, pr, dt_pr .* 299792458);
                    end
                end
                if any(dt_ph_dj)
                    this.log.addMessage(this.log.indent('Correcting carrier phases jumps'));
                else
                    this.log.addMessage(this.log.indent('Correcting carrier phases for a dt drift estimated from desync interpolation'));
                end
                if any(dt_pr_dj)
                    this.log.addMessage(this.log.indent('Correcting pseudo-ranges jumps'));
                else
                    this.log.addMessage(this.log.indent('Correcting pseudo-ranges for a dt drift estimated from desync interpolation'));
                end
            end
            
            % Saving dt into the object properties
            this.dt_ph = dt_ph; %#ok<PROP>
            this.dt_pr = dt_pr; %#ok<PROP>
            
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
        
        function s02 = initPositioning(this, sys_c)
            % run the most appropriate init prositioning step depending on the static flag
            % calls initStaticPositioning() or initDynamicPositioning()
            % SYNTAX
            %   this.initPositioning();
            % INPUT
            %   sys_c = wanted system
            % Init "errors"
            if this.isEmpty()
                this.log.addError('Init positioning failed: the receiver object is empty');
            else
                this.log.addMarkedMessage('Computing position and clock errors using a code only solution')
                this.sat.err_tropo = zeros(this.time.length, this.parent.cc.getMaxNumSat());
                this.sat.err_iono  = zeros(this.time.length, this.parent.cc.getMaxNumSat());
                this.sat.solid_earth_corr  = zeros(this.time.length, this.parent.cc.getMaxNumSat());
                this.log.addMessage(this.log.indent('Applying satellites Differential Code Biases (DCB)'))
                % if not applied apply group delay
                this.applyGroupDelay();
                this.log.addMessage(this.log.indent('Applying satellites clock errors and eccentricity dependent relativistic correction'))
                this.applyDtSat();
                
                %this.parent.static = 0;
                if this.isStatic()
                    %this.initStaticPositioningOld(obs, prn, sys, flag)
                    s02 = this.initStaticPositioning();
                else                    
                    s02 = this.initDynamicPositioning();
                end
            end
        end
        
        function s0 = initStaticPositioning(this)
            % SYNTAX
            %   this.StaticPositioning(obs, prn, sys, flag)
            %
            % INPUT
            % obs : observations [meters]
            % prn : prn of observations
            % sys : sys of observations
            % flag : name of observation [obs_code1 obs_code2 comb_code]
            %        comb_code --> Iono Free = I
            % OUTPUTt_glo11.
            %
            % 
            %   Get positioning using code observables
            
            
            if this.isEmpty()
                this.log.addError('Static positioning failed: the receiver object is empty');
            else
                if isempty(this.id_sync)
                    this.id_sync = 1:this.time.length;
                end
                % check if the epochs are present
                % getting the observation set that is going to be used in
                % setUPSA
                obs_set = Observation_Set();
                if this.isMultiFreq() %% case multi frequency
                    for sys_c = this.cc.sys_c
                        obs_set.merge(this.getPrefIonoFree('C', sys_c));
                    end
                else
                    for sys_c = this.cc.sys_c
                        f = this.getFreqs(sys_c);
                        obs_set.merge(this.getPrefObsSetCh(['C' num2str(f(1))], sys_c));
                    end
                end
                
                if sum(this.hasAPriori) == 0 %%% if no apriori information on the position
                    coarsePositioning(this, obs_set)
                end
                
                
                
                this.updateAllAvailIndex();
                this.updateAllTOT();
                this.updateAzimuthElevation()
                this.updateErrTropo();
                if ~this.isMultiFreq()
                    this.updateErrIono();
                end
                this.log.addMessage(this.log.indent('Improving estimation'))
                this.codeStaticPositioning(this.id_sync, 15);
                this.remBadTracking();
                
                this.updateAllTOT();
                this.log.addMessage(this.log.indent('Final estimation'))
                [~, s0] = this.codeStaticPositioning(this.id_sync, 15);
                this.log.addMessage(this.log.indent(sprintf('Estimation sigma0 %.3f m', s0) ))
                
                % final estimation of time of flight
                this.updateAllAvailIndex()
                this.updateAllTOT();
            end
        end
        
        function coarsePositioning(this, obs_set)
            % get a very coarse postioning for the receiver
            if nargin < 2
                obs_set = Observation_Set();
                if this.isMultiFreq() %% case multi frequency
                    for sys_c = this.cc.sys_c
                        obs_set.merge(this.getPrefIonoFree('C', sys_c));
                    end
                else
                    for sys_c = this.cc.sys_c
                        f = this.getFreqs(sys_c);
                        obs_set.merge(this.getPrefObsSetCh(['C' num2str(f(1))], sys_c));
                    end
                end
            end
            min_ep_thrs = 50;
            if this.time.length < min_ep_thrs
                min_ep_thrs = 1;
            end
            
            last_ep_coarse = min(100, this.time.length);
            ep_coarse = 1 : last_ep_coarse;
            while(not( sum(sum(obs_set.obs(ep_coarse,:) ~= 0, 2) > 2) > min_ep_thrs)) % checking if the selected epochs contains at least some usabele obseravables
                ep_coarse = ep_coarse + 1;
            end
            this.initAvailIndex(ep_coarse);
            this.updateAllTOT();
            
            this.log.addMessage(this.log.indent('Getting coarse position on subsample of data'))
            
            dpos = 3000; % 3 km - entry condition
            while max(abs(dpos)) > 10
                dpos = this.codeStaticPositioning(ep_coarse);
            end
            this.updateAzimuthElevation()
            this.updateErrTropo();
            if ~this.isMultiFreq()
                this.updateErrIono();
            end
            this.codeStaticPositioning(ep_coarse, 15);
        end
        
        function remBadTracking(this) %%% imprtnat check!! if ph obervation with no code are deleted elsewhere
            % requires approximate postioin and approx clock estimate
            [pr, id_pr] = this.getPseudoRanges;
            %[ph, wl, id_ph] = this.getPhases;
            sensor =  pr - this.getSyntPrObs -repmat(this.dt,1,size(pr,2))*Global_Configuration.V_LIGHT;
            bad_track = abs(sensor) > 1e4;
            bad_track = flagExpand(bad_track, 2);
            pr(bad_track) = 0;
            %ph(bad_track) = 0;
            this.setPseudoRanges(pr, id_pr);
            %this.setPhases(ph, wl, id_ph);
            
            
        end
        
        function [dpos, s02] = codeStaticPositioning(this, id_sync, cut_off, num_reweight)
            ls = Least_Squares_Manipulator();
            if nargin < 2
                if ~isempty(this.id_sync)
                    id_sync = this.id_sync;
                else
                    id_sync = 1:this.length();
                end
            end
            if nargin < 3
                cut_off = this.state.getCutOff();
            end
            if nargin < 4
                num_reweight = 0;
            end
            ls.setUpCodeSatic( this, id_sync, cut_off);
            ls.Astack2Nstack();
            [x, res, s02] = ls.solve();
            
            % REWEIGHT ON RESIDUALS -> (not well tested , uncomment to
            % enable)
            for i = 1 : num_reweight
                ls.reweightHuber();
                ls.Astack2Nstack();
                [x, res, s02] = ls.solve();
            end
            
            dpos = [x(x(:,2) == 1,1) x(x(:,2) == 2,1) x(x(:,2) == 3,1)];
            if isempty(dpos)
                dpos = [0 0 0];
            end
            this.xyz = this.getMedianPosXYZ() + dpos;
            dt = x(x(:,2) == 6,1);
            this.dt = zeros(this.time.length,1);
            this.dt(ls.true_epoch,1) = dt ./ Global_Configuration.V_LIGHT;
            isb = x(x(:,2) == 4,1);
            this.sat.res = zeros(this.length, this.parent.cc.getMaxNumSat());
            % LS does not know the max number of satellite stored
            dsz = max(id_sync) - size(res,1);
            if dsz == 0
                this.sat.res(id_sync, ls.sat_go_id) = res(id_sync, ls.sat_go_id);
            else
                if dsz > 0
                    id_sync = id_sync(1 : end - dsz);
                end
                this.sat.res(id_sync, ls.sat_go_id) = res(id_sync, ls.sat_go_id);
            end
            this.s0_ip = s02;
        end
        
        function [dpos, s0] = codeDynamicPositioning(this, id_sync, cut_off)
            ls = Least_Squares_Manipulator();
            if nargin < 2
                if ~isempty(this.id_sync)
                    id_sync = this.id_sync;
                else
                    id_sync = 1:this.length();
                end
            end
            if nargin < 3
                cut_off = this.state.getCutOff();
            end
            ls.setUpCodeDynamic( this, id_sync, cut_off);
            ls.Astack2Nstack();
            [x, res, s0] = ls.solve();
            
            % REWEIGHT ON RESIDUALS -> (not well tested , uncomment to
            % enable)
            %             ls.reweightHuber();
            %             ls.Astack2Nstack();
            %             [x, res, s02] = ls.solve();
            
            dpos = [x(x(:,2) == 1,1) x(x(:,2) == 2,1) x(x(:,2) == 3,1)];
            if size(this.xyz,1) == 1
                this.xyz = repmat(this.xyz,this.time.length,1);
            end
            this.xyz = this.xyz(id_sync,:) + dpos;
            dt = x(x(:,2) == 6,1);
            this.dt = zeros(this.time.length,1);
            this.dt(ls.true_epoch,1) = dt ./ Global_Configuration.V_LIGHT;
            isb = x(x(:,2) == 4,1);
            this.sat.res = zeros(this.length, this.parent.cc.getMaxNumSat());
            % LS does not know the max number of satellite stored
            
            
            this.sat.res = zeros(this.length, this.parent.cc.getMaxNumSat());
            dsz = max(id_sync) - size(res,1);
            if dsz == 0
                this.sat.res(id_sync, ls.sat_go_id) = res(id_sync, ls.sat_go_id);
            else
                if dsz > 0
                    id_sync = id_sync(1 : end - dsz);
                end
                this.sat.res(id_sync, ls.sat_go_id) = res(id_sync, ls.sat_go_id);
            end
        end
        
        function s02 = initDynamicPositioning(this)
            % SYNTAX
            %   this.StaticPositioning(obs, prn, sys, flag)
            %
            % INPUT
            % obs : observations [meters]
            % prn : prn of observations
            % sys : sys of observations
            % flag : name of observation [obs_code1 obs_code2 comb_code]
            %        comb_code --> Iono Free = I
            % OUTPUTt_glo11.
            %
            %
            %   Get positioning using code observables
            
            
            if this.isEmpty()
                this.log.addError('Static positioning failed: the receiver object is empty');
            else
                if isempty(this.id_sync)
                    this.id_sync = 1:this.time.length;
                end
                % check if the epochs are present
                % getting the observation set that is going to be used in
                % setUPSA
                obs_set = Observation_Set();
                if this.isMultiFreq() %% case multi frequency
                    for sys_c = this.cc.sys_c
                        obs_set.merge(this.getPrefIonoFree('C', sys_c));
                    end
                else
                    for sys_c = this.cc.sys_c
                        f = this.getFreqs(sys_c);
                        obs_set.merge(this.getPrefObsSetCh(['C' num2str(f(1))], sys_c));
                    end
                end
                
                if sum(this.hasAPriori) == 0 %%% if no apriori information on the position
                    coarsePositioning(this, obs_set)
                end
                
                
                this.updateAllAvailIndex();
                this.updateAllTOT();
                this.updateAzimuthElevation()
                this.updateErrTropo();
                if ~this.isMultiFreq()
                    this.updateErrIono();
                end
                this.log.addMessage(this.log.indent('Improving estimation'))
                this.codeDynamicPositioning(this.id_sync, 15);
                
                this.updateAllTOT();
                this.log.addMessage(this.log.indent('Final estimation'))
                [~, s02] = this.codeDynamicPositioning(this.id_sync, 15);
                this.log.addMessage(this.log.indent(sprintf('Estimation sigma02 %.3f m', s02) ))
                this.s02_ip = s02;
                
                % final estimation of time of flight
                this.updateAllAvailIndex()
                this.updateAllTOT();
                [~, s02] = this.codeDynamicPositioning(this.id_sync, 15);
            end
        end
        
        function [obs, prn, sys, flag] = getBestCodeObs(this)
            % INPUT
            % OUPUT:
            %    obs = observations [n_obs x n_epoch];
            %    prn = satellite prn [n_obs x 1];
            %    sys = system [n_obs x 1];
            %    flag = flag of the observation [ n_obs x 7] iono free
            %    combination are labeled with the obs code of both observations
            %  get "best" avaliable code or code combination for the given system
            n_epoch = size(this.obs,2);
            obs = [];
            sys = [];
            prn = [];
            flag = [];
            for i = unique(this.go_id)'
                
                sat_idx = this.getObsIdx('C',i);
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
                                inv_wl1 = 1/this.wl(this.getObsIdx(f_obs_code(k,:),i));
                                inv_wl2 = 1/this.wl(this.getObsIdx(s_obs_code(y,:),i));
                                %                                 if ((inv_wl1).^2 - (inv_wl2).^2) == 0
                                %                                     keyboard
                                %                                 end
                                
                                obs_tmp = ((inv_wl1).^2 .* first_freq(k,:) - (inv_wl2).^2 .* second_freq(y,:))./ ( (inv_wl1).^2 - (inv_wl2).^2 );
                                obs_tmp(isnan(obs_tmp)) = 0;
                                if sum(obs_tmp>0)>0
                                    obs = [obs; obs_tmp];
                                    [t_sys , t_prn] = this.getSysPrn(i);
                                    prn = [prn; t_prn];
                                    sys = [sys; t_sys];
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
            % remove obs for which coordinates of satellite are not available
            o = 1;
            while (o <= length(prn))
                s = this.cc.getIndex(sys(o),prn(o));
                o_idx_l = obs(o,:)>0;
                times = this.time.getSubSet(o_idx_l);
                times.addSeconds(-obs(o,o_idx_l)' / Global_Configuration.V_LIGHT); % add rough time of flight
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
                else
                    o = o + 1;
                end
            end
        end
        
        function computeBasicPosition(this, sys_c)
            % Compute a very simple position from code observations
            % without applying any corrections
            %
            % SYNTAX
            %   this.computeBasicPosition();
            
            if this.isEmpty()
                this.log.addError('Pre-Processing failed: the receiver object is empty');
            else
                if nargin < 2
                    sys_c = this.cc.sys_c;
                end
                this.setActiveSys(intersect(this.getActiveSys, this.sat.cs.getAvailableSys));
                this.remBad();
                % correct for raw estimate of clock error based on the phase measure
                this.correctTimeDesync();
                this.initPositioning(sys_c);
            end
        end
        
        function preProcessing(this, sys_c)
            % Do all operation needed in order to preprocess the data
            % remove bad observation (spare satellites or bad epochs from CRX)
            % SYNTAX
            %   this.preProcessing(sys_c);
            
            if this.isEmpty()
                this.log.addError('Pre-Processing failed: the receiver object is empty');
            else
                this.pp_status = false;
                if nargin < 2
                    sys_c = this.cc.sys_c;
                end
                this.setActiveSys(intersect(this.getActiveSys, this.sat.cs.getAvailableSys));
                this.remBad();
                this.remUnderSnrThr(this.state.getSnrThr());
                % correct for raw estimate of clock error based on the phase measure
                [is_pr_jumping, is_ph_jumping] = this.correctTimeDesync();
                % this.TEST_smoothCodeWithDoppler(51);
                % code only solution
                this.importMeteoData();
                this.initTropo();
                
                if (this.state.isOutlierRejectionOn())
                    this.remBadPrObs(150);
                end
                this.remShortArc(this.state.getMinArc);            

                s02 = this.initPositioning(sys_c); %#ok<*PROPLC>                
                if (min(s02) > this.S02_IP_THR)
                    this.log.addWarning(sprintf('Very BAD code solution => something is proably wrong (s02 = %.2f)', s02));
                else
                    this.remUnderCutOff();
                    this.setAvIdx2Visibility();
                    this.meteo_data = [];
                    this.importMeteoData();
                    
                    % if the clock is stable I can try to smooth more => this.smoothAndApplyDt([0 this.length/2]);
                    this.dt_ip = simpleFill1D(this.dt, this.dt == 0, 'linear') + this.getDtPr; % save init_positioning clock
                    % smooth clock estimation
                    this.smoothAndApplyDt(0, is_pr_jumping, is_ph_jumping);
                    
                    % update azimuth elevation
                    this.updateAzimuthElevation();
                    % Add a model correction for time desync -> observations are now referred to nominal time  #14
                    this.shiftToNominal()
                    this.updateAllTOT();
                    
                    % apply various corrections
                    this.sat.cs.toCOM(); %interpolation of attitude with 15min coordinate might possibly be inaccurate switch to center of mass (COM)
                    
                    this.updateAzimuthElevation();
                    ph0 = this.getPhases();
                    this.applyPCV();
                    ph1 = this.getPhases();
                    corr.pcv = ph1 - ph0;
                    this.applyPoleTide();
                    ph2 = this.getPhases();
                    corr.pt = ph2 - ph1;
                    this.applyPhaseWindUpCorr();
                    ph3 = this.getPhases();
                    corr.pwu = ph3 - ph2;
                    this.applySolidEarthTide();
                    ph4 = this.getPhases();
                    corr.set = ph4 - ph3;
                    this.applyShDelay();
                    ph5 = this.getPhases();
                    corr.shd = ph5 - ph4;
                    this.applyOceanLoading();
                    ph6 = this.getPhases();
                    corr.ocl = ph6 - ph5;
                    this.applyAtmLoad();
                    this.applyHOI();
                    
                    this.remOutlierMarkCycleSlip();
                    this.pp_status = true;
                end
            end
        end
        
        function staticPPP(this, sys_list, id_sync)
            % compute a static PPP solution
            %
            % SYNTAX
            %   this.parent.staticPPP(<sys_list>, <id_sync>)
            %
            % EXAMPLE:
            %   Use the full dataset to compute a PPP solution
            %    - this.parent.staticPPP();
            %
            %   Use just GPS + GLONASS + Galileo to compute a PPP solution
            %   using epochs from 501 to 2380
            %    - this.parent.staticPPP('GRE', 501:2380);
            %
            %   Use all the available satellite system to compute a PPP solution
            %   using epochs from 501 to 2380
            %    - this.parent.staticPPP([], 500:2380);
            
            if this.isEmpty()
                this.log.addError('staticPPP failed The receiver object is empty');
            elseif ~this.isPreProcessed()
                if ~isempty(this.s0_ip) && (this.s0_ip < Inf)
                    this.log.addError('Pre-Processing is required to compute a PPP solution');
                else
                    this.log.addError('Pre-Processing failed: skipping PPP solution');
                end
            else
                if nargin >= 2
                    if ~isempty(sys_list)
                        this.setActiveSys(sys_list);
                    end
                end
                
                if nargin < 3
                    id_sync = (1 : this.time.length())';
                end
                
                this.log.addMarkedMessage(['Computing PPP solution using: ' this.getActiveSys()]);
                this.log.addMessage(this.log.indent('Preparing the system'));
                %this.updateAllAvailIndex
                %this.updateAllTOT
                ls = Least_Squares_Manipulator();
                id_sync = ls.setUpPPP(this, id_sync);
                ls.Astack2Nstack();
                
                time = this.time.getSubSet(id_sync);
                
                if isempty(this.zwd) || all(isnan(this.zwd))
                    this.zwd = zeros(this.time.length(), 1);
                end
                if isempty(this.apr_zhd) || all(isnan(this.apr_zhd))
                    this.apr_zhd = zeros(this.time.length(),1);
                end
                if isempty(this.ztd) || all(isnan(this.ztd))
                    this.ztd = zeros(this.time.length(),1);
                end
                
                n_sat = size(this.sat.el,2);
                if isempty(this.sat.slant_td)
                    this.sat.slant_td = zeros(this.time.length(), n_sat);
                end
                
                rate = time.getRate();
                
                ls.setTimeRegularization(ls.PAR_CLK, 1e-3 * rate); % really small regularization
                ls.setTimeRegularization(ls.PAR_TROPO, (this.state.std_tropo)^2 / 3600 * rate );% this.state.std_tropo / 3600 * rate  );
                if this.state.flag_tropo_gradient
                    ls.setTimeRegularization(ls.PAR_TROPO_N, (this.state.std_tropo_gradient)^2 / 3600 * rate );%this.state.std_tropo / 3600 * rate );
                    ls.setTimeRegularization(ls.PAR_TROPO_E, (this.state.std_tropo_gradient)^2 / 3600 * rate );%this.state.std_tropo  / 3600 * rate );
                end
                this.log.addMessage(this.log.indent('Solving the system'));
                [x, res, s0] = ls.solve();
                % REWEIGHT ON RESIDUALS -> (not well tested , uncomment to
                % enable)
                % ls.reweightHuber();
                % ls.Astack2Nstack();
                % [x, res, s02] = ls.solve();
                this.id_sync = id_sync;
                
                this.sat.res = zeros(this.length, n_sat);
                this.sat.res(id_sync, ls.sat_go_id) = res(id_sync, ls.sat_go_id);
                
                %this.id_sync = unique([serialize(this.id_sync); serialize(id_sync)]);
                
                coo = [x(x(:,2) == 1,1) x(x(:,2) == 2,1) x(x(:,2) == 3,1)];
                if isempty(coo)
                    coo = [0 0 0];
                end
                
                clock = x(x(:,2) == 6,1);
                tropo = x(x(:,2) == 7,1);
                amb = x(x(:,2) == 5,1);
                
                % saving matrix of float ambiguities
                amb_mat = zeros(length(id_sync), length(ls.go_id_amb));
                id_ok = ~isnan(ls.amb_idx);
                amb_mat(id_ok) = amb(ls.amb_idx(id_ok));
                this.sat.amb_mat = nan(this.length, this.parent.cc.getMaxNumSat);
                this.sat.amb_mat(id_sync,ls.go_id_amb) = amb_mat;
                
                gntropo = x(x(:,2) == 8,1);
                getropo = x(x(:,2) == 9,1);
                this.log.addMessage(this.log.indent(sprintf('DEBUG: s02 = %f',s0)));
                this.xyz = this.xyz + coo;
                valid_ep = ls.true_epoch;
                this.dt(valid_ep, 1) = clock / Global_Configuration.V_LIGHT;
                this.sat.amb_idx = nan(this.length, this.parent.cc.getMaxNumSat);
                this.sat.amb_idx(id_sync,ls.go_id_amb) = ls.amb_idx;
                this.if_amb = amb; % to test ambiguity fixing
                this.s0 = s0;
                if s0 > 0.10
                    this.log.addWarning(sprintf('PPP solution failed, s02: %6.4f   - no update to receiver fields',s0))
                end
                if s0 < 1.50 % with over one meter of error the results are not meaningfull
                    if this.state.flag_tropo
                        zwd = this.getZwd();
                        zwd_tmp = zeros(size(this.zwd));
                        zwd_tmp(this.id_sync) = zwd;
                        this.zwd(valid_ep) = zwd_tmp(valid_ep) + tropo;
                        this.ztd(valid_ep) = this.zwd(valid_ep) + this.apr_zhd(valid_ep);
                        this.pwv = nan(size(this.zwd));
                        if ~isempty(this.meteo_data)
                            degCtoK = 273.15;
                            [~,Tall, H] = this.getPTH();
                            % weighted mean temperature of the atmosphere over Alaska (Bevis et al., 1994)
                            Tm = (Tall(valid_ep) + degCtoK)*0.72 + 70.2;
                            
                            % Askne and Nordius formula (from Bevis et al., 1994)
                            Q = (4.61524e-3*((3.739e5./Tm) + 22.1));
                            
                            % precipitable Water Vapor
                            this.pwv(valid_ep) = this.zwd(valid_ep) ./ Q * 1e3;
                        end
                        this.sat.amb = amb;
                        [mfh, mfw] = getSlantMF(this);
                        this.sat.slant_td(id_sync, :) = nan2zero(zero2nan(this.sat.res(id_sync, :)) ...
                            + zero2nan(repmat(this.zwd(id_sync, :), 1, n_sat).*mfw) ...
                            + zero2nan(repmat(this.apr_zhd(id_sync, :), 1, n_sat).*mfh));
                        this.est_slant = repmat(tropo, 1, n_sat) .*mfw .* this.sat.avail_index(id_sync, :);  % to test ambiguity fixing
                    end
                    if this.state.flag_tropo_gradient
                        if isempty(this.tgn) || all(isnan(this.tgn))
                            this.tgn = nan(this.time.length,1);
                        end
                        this.tgn(valid_ep) =  nan2zero(this.tgn(valid_ep)) + gntropo;
                        if isempty(this.tge) || all(isnan(this.tge))
                            this.tge = nan(this.time.length,1);
                        end
                        this.tge(valid_ep) = nan2zero(this.tge(valid_ep))  + getropo;
                        
                        cotel = zero2nan(cotd(this.sat.el(id_sync, :)));
                        cosaz = zero2nan(cosd(this.sat.az(id_sync, :)));
                        sinaz = zero2nan(sind(this.sat.az(id_sync, :)));
                        this.est_slant = nan2zero(zero2nan(this.est_slant) ...
                            + zero2nan(repmat(this.tgn(id_sync, :),1,n_sat) .* mfw .* cotel .* cosaz) ...
                            + zero2nan(repmat(this.tge(id_sync, :),1,n_sat) .* mfw .* cotel .* sinaz));  % to test ambiguity fixing
                    end
                    this.updateErrTropo();
                end
                this.pushResult();
            end
            
            
            
        end
        
        function dynamicPPP(this, sys_list, id_sync)
            % compute a static PPP solution
            %
            % SYNTAX
            %   this.parent.staticPPP(<sys_list>, <id_sync>)
            %
            % EXAMPLE:
            %   Use the full dataset to compute a PPP solution
            %    - this.parent.staticPPP();
            %
            %   Use just GPS + GLONASS + Galileo to compute a PPP solution
            %   using epochs from 501 to 2380
            %    - this.parent.staticPPP('GRE', 501:2380);
            %
            %   Use all the available satellite system to compute a PPP solution
            %   using epochs from 501 to 2380
            %    - this.parent.staticPPP([], 500:2380);
            
            if this.isEmpty()
                this.log.addError('staticPPP failed The receiver object is empty');
            elseif ~this.isPreProcessed()
                if (this.s02_ip < Inf)
                    this.log.addError('Pre-Processing is required to compute a PPP solution');
                else
                    this.log.addError('Pre-Processing failed: skipping PPP solution');
                end
            else
                if nargin >= 2
                    if ~isempty(sys_list)
                        this.setActiveSys(sys_list);
                    end
                end
                
                if nargin < 3
                    id_sync = (1 : this.time.length())';
                end
                
                this.log.addMarkedMessage(['Computing PPP solution using: ' this.getActiveSys()]);
                this.log.addMessage(this.log.indent('Preparing the system'));
                %this.updateAllAvailIndex
                %this.updateAllTOT
                ls = Least_Squares_Manipulator();
                id_sync = ls.setUpPPP(this, id_sync, this.state.getCutOff, true);
                ls.Astack2Nstack();
                
                time = this.time.getSubSet(id_sync);
                
                if isempty(this.zwd) || all(isnan(this.zwd))
                    this.zwd = zeros(this.time.length(), 1);
                end
                if isempty(this.apr_zhd) || all(isnan(this.apr_zhd))
                    this.apr_zhd = zeros(this.time.length(),1);
                end
                if isempty(this.ztd) || all(isnan(this.ztd))
                    this.ztd = zeros(this.time.length(),1);
                end
                
                n_sat = size(this.sat.el,2);
                if isempty(this.sat.slant_td)
                    this.sat.slant_td = zeros(this.time.length(), n_sat);
                end
                
                rate = time.getRate();
                
                %ls.setTimeRegularization(ls.PAR_CLK, 1e-3 * rate); % really small regularization
                ls.setTimeRegularization(ls.PAR_X, 1 * rate); % really small regularization
                ls.setTimeRegularization(ls.PAR_Y, 1 * rate); % really small regularization
                ls.setTimeRegularization(ls.PAR_Z, 1 * rate); % really small regularization
                ls.setTimeRegularization(ls.PAR_TROPO, (this.state.std_tropo)^2 / 3600 * rate );% this.state.std_tropo / 3600 * rate  );
                if this.state.flag_tropo_gradient
                    ls.setTimeRegularization(ls.PAR_TROPO_N, (this.state.std_tropo_gradient)^2 / 3600 * rate );%this.state.std_tropo / 3600 * rate );
                    ls.setTimeRegularization(ls.PAR_TROPO_E, (this.state.std_tropo_gradient)^2 / 3600 * rate );%this.state.std_tropo  / 3600 * rate );
                end
                this.log.addMessage(this.log.indent('Solving the system'));
                [x, res, s02] = ls.solve();
                % REWEIGHT ON RESIDUALS -> (not well tested , uncomment to
                % enable)
                % ls.reweightHuber();
                % ls.Astack2Nstack();
                % [x, res, s02] = ls.solve();
                this.id_sync = id_sync;
                
                this.sat.res = zeros(this.length, n_sat);
                this.sat.res(id_sync, ls.sat_go_id) = res(id_sync, ls.sat_go_id);
                
                %this.id_sync = unique([serialize(this.id_sync); serialize(id_sync)]);
                
                coo_x    = x(x(:,2) == 1,1);
                coo_y    = x(x(:,2) == 2,1);
                coo_z    = x(x(:,2) == 3,1);
                
                clock = x(x(:,2) == 6,1);
                tropo = x(x(:,2) == 7,1);
                amb = x(x(:,2) == 5,1);
                
                this.sat.res = zeros(this.length, n_sat);
                this.sat.res(id_sync, ls.sat_go_id) = res(id_sync, ls.sat_go_id);
                keyboard
            end
            
        end
        
        function fixPPP(this)
            % coarse estimation of amb
            this.coarseAmbEstimation();
            % apply dt
            this.smoothAndApplyDt(0);
            % remove grupo delay
            this.remGroupDelay();
            % estaimet WL
            mwb = this.getMelWub('1','2','G');
            wl_cycle = mwb.obs./repmat(mwb.wl,size(mwb.obs,1),1);
            wl_cycle = zeros(size(mwb.obs,1),this.cc.getGPS.N_SAT);
            wl_cycle(:,mwb.go_id) = mwb.obs;
            wl_cycle(:,mwb.go_id) = wl_cycle(:,mwb.go_id)./repmat(mwb.wl,size(mwb.obs,1),1);
            wsb = this.sat.cs.getWSB(this.getCentralTime());
            % take off wsb
            wl_cycle = zero2nan(wl_cycle) + repmat(wsb,size(mwb.obs,1),1);
            % get receiver wsb
            wsb_rec = median(zero2nan(wl_cycle(:)) - ceil(zero2nan(wl_cycle(:))),'omitnan'); 
            wl_cycle = wl_cycle - wsb_rec;
            keyboard
            % estimate narrow lane
            % -> amb VCV ???
            % get undifferenced
            % new estimation round with fixed ambs
        end
        
        function coarseAmbEstimation(this)
            % estimate coarse ambiguities
            %%% experimental now just for GPS L1 L2
            
            % get geometry free
            wl1 = GPS_SS.L_VEC(1);
            wl2 = GPS_SS.L_VEC(2);
            
            f1 = GPS_SS.F_VEC(1);
            f2 = GPS_SS.F_VEC(2);
            
            gamma = f1^2 / f2^2;
            
            amb_idx = this.sat.amb_idx;
            n_amb  = max(max(amb_idx));
            geom_free = zeros(this.time.length,GPS_SS.N_SAT);
            obs_set =  this.getGeometryFree('C2W','C1W','G');
            geom_free(:,obs_set.go_id) = obs_set.obs;
            iono_elong = geom_free / ( 1 -gamma); % Eq (2)
            P1_idx = idxCharLines(this.obs_code, 'C1W')'& this.system == 'G';
            P2_idx = idxCharLines(this.obs_code, 'C2W')'& this.system == 'G';
            P1 = zeros(this.time.length,GPS_SS.N_SAT);
            P2 = zeros(this.time.length,GPS_SS.N_SAT);
            P1(:,this.go_id(P1_idx)) = this.obs(P1_idx,:)';
            P2(:,this.go_id(P2_idx)) = this.obs(P2_idx,:)';
            L1 = zeros(this.time.length,GPS_SS.N_SAT);
            L2 = zeros(this.time.length,GPS_SS.N_SAT);
            L1_idx = idxCharLines(this.obs_code, 'L1C')' & this.system == 'G';
            L1(:,this.go_id(L1_idx)) = this.obs(L1_idx, :)';
            L2_idx = idxCharLines(this.obs_code, 'L2W')' & this.system == 'G';
            L2(:,this.go_id(L2_idx)) = this.obs(L2_idx, :)';
            n1_tilde = (zero2nan(P1) - 2* zero2nan(iono_elong)) / wl1 - zero2nan(L1);
            n2_tilde = (zero2nan(P2) - 2* gamma* zero2nan(iono_elong)) / wl2 - zero2nan(L2);
            
            
            l1_amb = zeros(1,n_amb);
            l2_amb = zeros(1,n_amb);
            for i =1 : n_amb;
                l1_amb(i) = median(n1_tilde(amb_idx == i),'omitnan');
                l2_amb(i) = median(n2_tilde(amb_idx == i),'omitnan');
            end
            
            l1_amb_round = round(l1_amb);
            l1_amb_round(isnan(l1_amb_round)) = 0;
            l2_amb_round = round(l2_amb);
            l2_amb_round(isnan(l2_amb_round)) = 0;
            for i =1 : n_amb;
                idx_amb = amb_idx == i;
                L1(idx_amb) = L1(idx_amb) + l1_amb_round(i);
                L2(idx_amb) = L2(idx_amb) + l2_amb_round(i);
            end
            this.obs(L1_idx, :) = L1(:,this.go_id(L1_idx))';
            this.obs(L2_idx, :) = L2(:,this.go_id(L2_idx))';
        end
        
        function interpResults(this)
            % When computing a solution with subsampling not all the epochs have an estimation
            id_rem = true(this.time.length, 1);
            id_rem(this.id_sync) = false;
            
            if ~(isempty(this.dt))
                this.dt(id_rem) = nan;
                this.dt(id_rem) = interp1(this.time.getEpoch(this.id_sync).getRefTime, this.dt(this.id_sync), this.time.getEpoch(id_rem).getRefTime, 'pchip');
            end
            if ~(isempty(this.apr_zhd))
                this.apr_zhd(id_rem) = nan;
                this.apr_zhd(id_rem) = interp1(this.time.getEpoch(this.id_sync).getRefTime, this.apr_zhd(this.id_sync), this.time.getEpoch(id_rem).getRefTime, 'pchip');
            end
            if ~(isempty(this.ztd))
                this.ztd(id_rem) = nan;
                this.ztd(id_rem) = interp1(this.time.getEpoch(this.id_sync).getRefTime, this.ztd(this.id_sync), this.time.getEpoch(id_rem).getRefTime, 'pchip');
            end
            if ~(isempty(this.zwd))
                this.zwd(id_rem) = nan;
                this.zwd(id_rem) = interp1(this.time.getEpoch(this.id_sync).getRefTime, this.zwd(this.id_sync), this.time.getEpoch(id_rem).getRefTime, 'pchip');
            end
            if ~(isempty(this.pwv))
                this.pwv(id_rem) = nan;
                this.pwv(id_rem) = interp1(this.time.getEpoch(this.id_sync).getRefTime, this.pwv(this.id_sync), this.time.getEpoch(id_rem).getRefTime, 'pchip');
            end
            if ~(isempty(this.tgn))
                this.tgn(id_rem) = nan;
                this.tgn(id_rem) = interp1(this.time.getEpoch(this.id_sync).getRefTime, this.tgn(this.id_sync), this.time.getEpoch(id_rem).getRefTime, 'pchip');
            end
            if ~(isempty(this.tge))
                this.tge(id_rem) = nan;
                this.tge(id_rem) = interp1(this.time.getEpoch(this.id_sync).getRefTime, this.tge(this.id_sync), this.time.getEpoch(id_rem).getRefTime, 'pchip');
            end
            % slant cannot be interpolated easyle, it's better to interpolate slant_std and reconvert them to slant_td
            % not yet implemented
            %if ~(isempty(this.sat.slant_td))
            %    for s = 1 : size(this.sat.slant_td, 2)
            %        this.sat.slant_td(id_rem, s) = nan;
            %        id_ok = ~id_rem & this.sat.el(:, s) > 0 & this.sat.slant_td(:, s) > 0;
            %        id_interp = id_rem & this.sat.el(:, s) > 0;
            %        if sum(id_ok) > 3
            %            this.sat.slant_td(id_interp, s) = interp1(this.time.getEpoch(id_ok).getRefTime, zero2nan(this.sat.slant_td(id_ok, s)), this.time.getEpoch(id_interp).getRefTime, 'pchip');
            %        end
            %    end
            %end
        end
    end
    
    
    % ==================================================================================================================================================
    %% METHODS IMPORT / EXPORT
    % ==================================================================================================================================================
    
    methods        
        function exportRinex3(this, file_name)
            % Export the content of the object as RINEX 3
            % SYNTAX
            %   this.exportRinex3(file_name);
            
            this.remAllCorrections();
            
            % HEADER -------------------------------------------------------------------------------------------------------------------------------------------
            txt = [];
            this.updateRinObsCode();
            
            program = 'goGPS';
            agency = 'GReD s.r.l.';
            
            date_str_UTC = local_time_to_utc(now, 30);
            [date_str_UTC, time_str_UTC] = strtok(date_str_UTC,'T');
            time_str_UTC = time_str_UTC(2:end-1);
            
            %write header
            sys_c = unique(this.system);
            txt = sprintf('%s%9.2f           OBSERVATION DATA    %-19s RINEX VERSION / TYPE\n', txt, 3.03, iif(numel(sys_c) == 1, sys_c, 'M: Mixed'));
            txt = sprintf('%s%-20s%-20s%-20sPGM / RUN BY / DATE \n', txt, program, agency, [date_str_UTC ' ' time_str_UTC ' UTC']);
            txt = sprintf('%s%-60sMARKER NAME         \n', txt, this.parent.marker_name);
            txt = sprintf('%s%-20s                                        MARKER TYPE\n', txt, this.parent.marker_type);
            txt = sprintf('%s%-20s%-40sOBSERVER / AGENCY\n', txt, this.parent.observer, this.parent.agency);
            txt = sprintf('%sRINEX FILE EXPORTED BY goGPS OPEN SOURCE SOFTWARE           COMMENT\n', txt);
            txt = sprintf('%sREFERENCE DEV SITE: https://github.com/goGPS-Project        COMMENT\n', txt);
            txt = sprintf('%s%-20s%-20s%-20sREC # / TYPE / VERS\n', txt, this.parent.number, this.parent.type, this.parent.version);
            txt = sprintf('%s%-20s%-20s                    ANT # / TYPE\n', txt, this.parent.ant, this.parent.ant_type);
            xyz = this.getAPrioriPos();
            if ~any(xyz)
                xyz = this.getMedianPosXYZ();
            end
            txt = sprintf('%s%14.4f%14.4f%14.4f                  APPROX POSITION XYZ\n', txt, xyz(1), xyz(2), xyz(3));
            txt = sprintf('%s%14.4f%14.4f%14.4f                  ANTENNA: DELTA H/E/N\n', txt, this.parent.ant_delta_h, this.parent.ant_delta_en);
            
            sys = this.getActiveSys();
            for s = 1 : length(sys)
                obs_type = this.getAvailableCode(sys(s));
                % Set order CODE PHASE DOPPLER SNR
                obs_type = obs_type([find(obs_type(:,1) == 'C') ...
                                     find(obs_type(:,1) == 'L') ...
                                     find(obs_type(:,1) == 'D') ...
                                     find(obs_type(:,1) == 'S')], :);
                                 
                % if the code is unknown use 'the least important code' non empty as obs code
                % relative to the band
                for c = find(obs_type(:,3) == ' ')'
                    ss = this.cc.getSys(sys(s));
                    band = (ss.CODE_RIN3_2BAND == obs_type(c,2));
                    code = ss.CODE_RIN3_DEFAULT_ATTRIB{band}(1);
                    obs_type(c, 3) = code;
                end                 
                n_type = length(obs_type);
                n_line = ceil(n_type / 13);
                tmp = char(32 * ones(n_line * 13, 4));
                tmp(1 : n_type, 1:3) = obs_type;
                tmp = reshape(tmp', 4*13, n_line)';
                txt = sprintf('%s%-1s   %2d %s SYS / # / OBS TYPES\n', txt, char(sys(s)), n_type, tmp(1, :));
                for l = 2 : n_line
                    txt = sprintf('%s       %s SYS / # / OBS TYPES\n', txt, tmp(l, :));
                end
            end
            %txt = sprintf('%sDBHZ                                                        SIGNAL STRENGTH UNIT\n', txt);
            txt = sprintf('%s%10.3f                                                  INTERVAL\n', txt, this.time.getRate());
            txt = sprintf('%s%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS\n', txt, this.time.first.get6ColDate());
            txt = sprintf('%sCARRIER PHASE SHIFT remD BY goGPS SOFTWARE.              COMMENT\n', txt);
            for s = 1 : length(sys)
                obs_type = this.getAvailableCode(sys(s));
                % Set order CODE PHASE DOPPLER SNR
                obs_type = obs_type(obs_type(:,1) == 'L', :);
                for c = find(obs_type(:,3) == ' ')'
                    ss = this.cc.getSys(sys(s));
                    band = (ss.CODE_RIN3_2BAND == obs_type(c,2));
                    code = ss.CODE_RIN3_DEFAULT_ATTRIB{band}(1);
                    obs_type(c, 3) = code;
                end
                for l = 1 : size(obs_type, 1)
                    txt = sprintf('%s%-1s %-3s %8.5f                                              SYS / PHASE SHIFT \n', txt, sys(s), obs_type(l, :), 0.0);
                end
            end
            if ~isempty(intersect(sys, 'R'))
                % If glonas is present
                %txt = sprintf('%s C1C    0.000 C1P    0.000 C2C    0.000 C2P    0.000        GLONASS COD/PHS/BIS\n', txt);
                txt = sprintf('%s C1C          C1P          C2C          C2P                 GLONASS COD/PHS/BIS\n', txt);
            end
            txt  = sprintf('%s                                                            END OF HEADER       \n', txt);
            
            fid = fopen(file_name, 'w');
            fprintf(fid, '%s', txt);
            txt = [];
            
            % DATA ---------------------------------------------------------------------------------------------------------------------------------------------
            
            date6col = this.time.get6ColDate();
            flag_ok = 0;
            clock_offset = 0;
            
            obs_code = Core_Utils.code3Char2Num(this.obs_code);
            for ss = sys
                % rin_obs_code tu num;
                obs_type = reshape(this.rin_obs_code.(ss)', 3, length(this.rin_obs_code.(ss))/3)';
                % Set order CODE PHASE DOPPLER SNR
                obs_type = obs_type([find(obs_type(:,1) == 'C') ...
                    find(obs_type(:,1) == 'L') ...
                    find(obs_type(:,1) == 'D') ...
                    find(obs_type(:,1) == 'S')], :);
                                
                rin_obs_code.(ss) = Core_Utils.code3Char2Num(obs_type);
                for t = 1 : numel(rin_obs_code.(ss))
                    obs_code(obs_code == rin_obs_code.(ss)(t)) = t;
                end
            end
            
            n_epo = size(this.obs, 2);
            this.w_bar.createNewBar(' Exporting Receiver as Rinex 3.03 file...');
            this.w_bar.setBarLen(n_epo);
            for e = 1 : n_epo
                id_ok = ~isnan(zero2nan(this.obs(:, e)));
                go_id = unique(this.go_id(id_ok));
                if ~isempty(id_ok)
                    %txt = sprintf('%s> %4d %02d %02d %02d %02d %10.7f  %d%3d      %15.12f\n', txt, date6col(e, :), flag_ok, numel(go_id), clock_offset);
                    txt = sprintf('%s> %4d %02d %02d %02d %02d %10.7f  %d%3d\n', txt, date6col(e, :), flag_ok, numel(go_id));
                    for ss = sys % for each satellite system in view at this epoch
                        id_ss = id_ok & this.system' == ss;
                        for prn = unique(this.prn(id_ss))' % for each satellite in view at this epoch (of the current ss)
                            % find which obs type are present
                            id_sat_obs = (id_ss & this.prn == prn);
                            str_obs = char(32 * ones(16, numel(rin_obs_code.(ss))));
                            id_rin_col = obs_code(id_sat_obs);
                            str_obs(:, id_rin_col) = reshape(sprintf('%14.3f  ', this.obs(id_sat_obs, e)),16, numel(id_rin_col));
                            txt = sprintf('%s%c%02d%s\n', txt, ss, prn, str_obs(:)');
                        end
                    end
                end
                fprintf(fid, '%s', txt);
                txt = [];
                this.w_bar.go(e);
            end
            fclose(fid);
            this.log.newLine()
            this.log.addMarkedMessage(sprintf('Receiver exported successifully into: %s', file_name));
        end
        
        function pushResult(this)
            this.parent.out.injectResult(this);
        end
    end
    
     
    % ==================================================================================================================================================
    %% METHODS UTILITIES
    % ==================================================================================================================================================
    
    methods (Access = public)        
        function data_s = smoothSatData(this, data_az, data_el, data_in, cs_mat, method, spline_base, max_gap)
            if nargin < 5
                cs_mat = [];
            end
            if nargin < 7
                spline_base = (300 / this.getRate); % 5 min
            end
            if nargin < 6 
                method = 'spline';
            end
            if nargin < 8
                max_gap = 0;
            end
            if strcmp(method,'spline')
                data_s = data_in;
                for s = 1 : size(data_s, 2)
                    if isempty(cs_mat)
                        lim = getOutliers(~isnan(data_s(:,s)));
                    else
                        lim = getOutliers(~isnan(data_s(:,s)), cs_mat(:,s));
                    end
                    if max_gap > 0
                        lim = limMerge(lim, max_gap);
                    end

                    % remove small intervals
                    %lim((lim(:,2) - lim(:,1)) < spline_base, :) = [];
                    for l = 1 : size(lim, 1)
                        arc_size = lim(l,2) - lim(l,1) + 1;
                        id_arc = lim(l,1) : lim(l,2);
                        id_arc(isnan(data_s(id_arc, s))) = [];
                        data_tmp = data_s(id_arc, s);
                        if length(data_tmp) > 3
                            data_s(id_arc, s) = splinerMat([], data_tmp, min(arc_size, spline_base));
                        end
                    end
                end
            elseif strcmp(method,'poly_quad')
                data_s = data_in;
                for s = 1 : size(data_s, 2)
                    lim = getOutliers(~isnan(data_s(:,s)), cs_mat(:,s));
                    % remove small intervals
                    lim((lim(:,2) - lim(:,1)) < spline_base, :) = [];
                    for l = 1 : size(lim, 1)
                        n_ep =  lim(l,2) - lim(l,1) +1;
                        if n_ep > 4
                            data_tmp = data_s(lim(l,1) : lim(l,2), s);
                            el_tmp = data_el(lim(l,1) : lim(l,2), s);
                            az_tmp = data_az(lim(l,1) : lim(l,2), s);
                            t = [1:n_ep]';
                            A = [ones(n_ep,1) t t.^2];
                            %                             Qxx = cholinv(A'*A);
                            par = A\data_tmp;
                            %                             par = Qxx*A'*data_tmp;
                            quad_mod = A*par;
                            
                            %                             data_coll = data_tmp - quad_mod;
                            %                             s02 = mean(data_coll.^2) / (n_ep -3);
                            %                             Cy_haty_hat = s02*A'*Qxx;
                            %                             Cvv =  0;
                            data_s(lim(l,1) : lim(l,2), s) = splinerMat(te, data_coll , 1) + quad_mod;
                            data_s(lim(l,1) : lim(l,2), s) = quad_mod;
                        end
                    end
                end
            end
        end
        
          
        function [pr, sigma_pr] = smoothCodeWithDoppler(pr, sigma_pr, pr_go_id, dp, sigma_dp, dp_go_id)
            % NOT WORKING due to iono
            % Smooth code with phase (aka estimate phase ambiguity and remove it)
            % Pay attention that the bias between phase and code is now eliminated
            % At the moment the sigma of the solution = sigma of phase
            %
            % SYNTAX
            %   [pr, sigma_ph] = smoothCodeWithDoppler(pr, sigma_pr, pr_go_id, ph, sigma_ph, ph_go_id, cs_mat)
            %
            % EXAMPLE
            %   [pr.obs, pr.sigma] = Receiver.smoothCodeWithDoppler(zero2nan(pr.obs), pr.sigma, pr.go_id, ...
            %                                                       zero2nan(dp.obs), dp.sigma, dp.go_id);
            %            

            for s = 1 : numel(dp_go_id)
                s_c = find(pr_go_id == dp_go_id(s));
                pr(isnan(dp(:,s)), s_c) = nan;
                
                lim = getOutliers(~isnan(pr(:,s)));
                for l = 1 : size(lim, 1)
                    id_arc = (lim(l,1) : lim(l,2))';
                    
                    len_a = length(id_arc);
                    A = [speye(len_a); [-speye(len_a - 1) sparse(len_a - 1, 1)] + [sparse(len_a - 1, 1) speye(len_a - 1)]];
                    Q = speye(2 * len_a - 1);
                    Q = spdiags([ones(len_a,1) * sigma_pr(s_c).^2; ones(len_a-1,1) * (2 * sigma_dp(s).^2)], 0, Q);
                    Tn = A'/Q;
                    data_tmp = [pr(id_arc, s_c); dp(id_arc, s)];
                    pr(id_arc, s) = (Tn*A)\Tn * data_tmp;
                end
            end
        end
        
        function [ph, sigma_ph] = smoothCodeWithPhase(pr, sigma_pr, pr_go_id, ph, sigma_ph, ph_go_id, cs_mat)
            % NOT WORKING when iono is present
            % Smooth code with phase (aka estimate phase ambiguity and remove it)
            % Pay attention that the bias between phase and code is now eliminated
            % At the moment the sigma of the solution = sigma of phase
            %
            % SYNTAX
            %   [ph, sigma_ph] = smoothCodeWithPhase(pr, sigma_pr, pr_go_id, ph, sigma_ph, ph_go_id, cs_mat)
            %
            % EXAMPLE
            %   [ph.obs, ph.sigma] = Receiver.smoothCodeWithPhase(zero2nan(pr.obs), pr.sigma, pr.go_id, ...
            %                                                   zero2nan(ph.obs), ph.sigma, ph.go_id, ph.cycle_slip);
            %            

            for s = 1 : numel(ph_go_id)
                s_c = find(pr_go_id == ph_go_id(s));
                pr(isnan(ph(:,s)), s_c) = nan;
                
                lim = getOutliers(~isnan(ph(:,s)), cs_mat(:,s));
                for l = 1 : size(lim, 1)
                    id_arc = (lim(l,1) : lim(l,2))';
                    
                    len_a = length(id_arc);
                    A = [speye(len_a); [-speye(len_a - 1) sparse(len_a - 1, 1)] + [sparse(len_a - 1, 1) speye(len_a - 1)]];
                    Q = speye(2 * len_a - 1);
                    Q = spdiags([ones(len_a,1) * sigma_pr(s_c).^2; ones(len_a-1,1) * (2 * sigma_ph(s).^2)], 0, Q);
                    Tn = A'/Q;
                    data_tmp = [pr(id_arc, s_c); diff(ph(id_arc, s))];
                    ph(id_arc, s) = (Tn*A)\Tn * data_tmp;
                end
            end
        end
        
        function obs_num = obsCode2Num(obs_code)
            % Convert a 3 char name into a numeric value (float)
            % SYNTAX
            %   obs_num = obsCode2Num(obs_code);            
            obs_num = Core_Utils.code3Char2Num(obs_code(:,1:3));
        end
        
        function obs_code = obsNum2Code(obs_num)
            % Convert a numeric value (float) of an obs_code into a 3 char marker
            % SYNTAX
            %   obs_code = obsNum2Code(obs_num)
            obs_code = Core_Utils.num2Code3Char(obs_num);            
        end
        
        function marker_num = markerName2Num(marker_name)
            % Convert a 4 char name into a numeric value (float)
            % SYNTAX
            %   marker_num = markerName2Num(marker_name);
            marker_num = Core_Utils.code4Char2Num(marker_name(:,1:4));
        end
        
        function marker_name = markerNum2Name(marker_num)
            % Convert a numeric value (float) of a station into a 4 char marker
            % SYNTAX
            %   marker_name = markerNum2Name(marker_num)
            marker_name = Core_Utils.num2Code4Char(marker_num);
        end
        
        function [y0, pc, wl, ref] = prepareY0(trg, mst, lambda, pivot)
            % prepare y0 and pivot_correction arrays (phase only)
            % SYNTAX [y0, pc] = prepareY0(trg, mst, lambda, pivot)
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
            % SYNTAX y0 = composeY0(y0, pc, wl)
            y0 = serialize((y0 - pc) .* wl);
            y0(y0 == 0) = []; % remove pivots
        end
        
                      
        function sync(rec, rate)
            % keep epochs at a certain rate for a certain constellation
            %
            % SYNTAX
            %   this.keep(rate, sys_list)
            if nargin > 1 && ~isempty(rate)
                [~, id_sync] = GNSS_Station.getSyncTimeExpanded(rec, rate);
            else
                [~, id_sync] = GNSS_Station.getSyncTimeExpanded(rec);
            end
            % Keep the epochs in common
            % starting from the first when all the receivers are available
            % ending whit the last when all the receivers are available
            id_sync((find(sum(isnan(id_sync), 2) == 0, 1, 'first') : find(sum(isnan(id_sync), 2) == 0, 1, 'last')), :);
            
            % keep only synced epochs
            for r = 1 : numel(rec)
                rec(r).keepEpochs(id_sync(~isnan(id_sync(:, r)), r));
            end
            
        end
        
        function [res_ph1, mean_res, var_res] = legacyGetResidualsPh1(res_bin_file_name)
            %res_code1_fix  = [];                      % double differences code residuals (fixed solution)
            %res_code2_fix  = [];                      % double differences code residuals (fixed solution)
            %res_phase1_fix = [];                      % phase differences phase residuals (fixed solution)
            %res_phase2_fix = [];                      % phase differences phase residuals (fixed solution)
            %res_code1_float  = [];                    % double differences code residuals (float solution)
            %res_code2_float  = [];                    % double differences code residuals (float solution)
            res_phase1_float = [];                     % phase differences phase residuals (float solution)
            %res_phase2_float = [];                    % phase differences phase residuals (float solution)
            %outliers_code1 = [];                      % code double difference outlier? (fixed solution)
            %outliers_code2 = [];                      % code double difference outlier? (fixed solution)
            %outliers_phase1 = [];                     % phase double difference outlier? (fixed solution)
            %outliers_phase2 = [];                     % phase double difference outlier? (fixed solution)
            % observations reading
            log = Logger.getInstance();
            log.addMessage(['Reading: ' File_Name_Processor.getFileName(res_bin_file_name)]);
            d = dir(res_bin_file_name);
            fid_sat = fopen(res_bin_file_name,'r+');       % file opening
            num_sat = fread(fid_sat, 1, 'int8');                            % read number of satellites
            num_bytes = d.bytes-1;                                          % file size (number of bytes)
            num_words = num_bytes / 8;                                      % file size (number of words)
            num_packs = num_words / (2*num_sat*6);                          % file size (number of packets)
            buf_sat = fread(fid_sat,num_words,'double');                    % file reading
            fclose(fid_sat);                                                % file closing
            %res_code1_fix    = [res_code1_fix    zeros(num_sat,num_packs)]; % observations concatenation
            %res_code2_fix    = [res_code2_fix    zeros(num_sat,num_packs)];
            %res_phase1_fix   = [res_phase1_fix   zeros(num_sat,num_packs)];
            %res_phase2_fix   = [res_phase2_fix   zeros(num_sat,num_packs)];
            %res_code1_float  = [res_code1_float  zeros(num_sat,num_packs)];
            %res_code2_float  = [res_code2_float  zeros(num_sat,num_packs)];
            res_phase1_float = [res_phase1_float zeros(num_sat,num_packs)];
            %res_phase2_float = [res_phase2_float zeros(num_sat,num_packs)];
            %outliers_code1   = [outliers_code1   zeros(num_sat,num_packs)];
            %outliers_code2   = [outliers_code2   zeros(num_sat,num_packs)];
            %outliers_phase1  = [outliers_phase1  zeros(num_sat,num_packs)];
            %outliers_phase2  = [outliers_phase2  zeros(num_sat,num_packs)];
            i = 0;                                                           % epoch counter
            for j = 0 : (2*num_sat*6) : num_words-1
                i = i+1;                                                     % epoch counter increase
                %res_code1_fix(:,i)    = buf_sat(j + [1:num_sat]);            % observations logging
                %res_code2_fix(:,i)    = buf_sat(j + [1*num_sat+1:2*num_sat]);
                %res_phase1_fix(:,i)   = buf_sat(j + [2*num_sat+1:3*num_sat]);
                %res_phase2_fix(:,i)   = buf_sat(j + [3*num_sat+1:4*num_sat]);
                %res_code1_float(:,i)  = buf_sat(j + [4*num_sat+1:5*num_sat]);
                %res_code2_float(:,i)  = buf_sat(j + [5*num_sat+1:6*num_sat]);
                res_phase1_float(:,i) = buf_sat(j + (6*num_sat+1:7*num_sat));
                %res_phase2_float(:,i) = buf_sat(j + [7*num_sat+1:8*num_sat]);
                %outliers_code1(:,i)   = buf_sat(j + [8*num_sat+1:9*num_sat]);
                %outliers_code2(:,i)   = buf_sat(j + [9*num_sat+1:10*num_sat]);
                %outliers_phase1(:,i)  = buf_sat(j + [10*num_sat+1:11*num_sat]);
                %outliers_phase2(:,i)  = buf_sat(j + [11*num_sat+1:12*num_sat]);
            end
            res_ph1 = res_phase1_float';
            mean_res = mean(mean(zero2nan(res_ph1'), 'omitnan'), 'omitnan');
            var_res = max(var(zero2nan(res_ph1'), 'omitnan'));
        end
    end        
    %% METHODS PLOTTING FUNCTIONS
    % ==================================================================================================================================================
    
    % Various debug images
    % name variant:
    %   c cartesian
    %   s scatter
    %   p polar
    %   m mixed
    methods (Access = public)
        
        function showAll(this)
            this.toString;
            this.showDataAvailability();
            this.showSNR_p();
            this.showDt();
            this.showAll@Receiver_Commons();
        end
         
        function showDt(this)
            % Plot Clock error
            %
            % SYNTAX
            %   this.plotDt
            
            rec = this;
            if ~isempty(rec)
                f = figure; f.Name = sprintf('%03d: Dt Err', f.Number); f.NumberTitle = 'off';
                t = rec.time.getMatlabTime();
                nans = zero2nan(double(~rec.getMissingEpochs()));
                plot(t, rec.getDesync .* nans, '-k', 'LineWidth', 2);
                hold on;
                plot(t, rec.getDtPr .* nans, ':', 'LineWidth', 2);
                plot(t, rec.getDtPh .* nans, ':', 'LineWidth', 2);
                plot(t, (rec.getDtIP - rec.getDtPr) .* nans, '-', 'LineWidth', 2);
                plot(t, rec.getDtPrePro .* nans, '-', 'LineWidth', 2);
                if any(rec.getDt)
                    plot(t, rec.getDt .* nans, '-', 'LineWidth', 2);
                    plot(t, rec.getTotalDt .* nans, '-', 'LineWidth', 2);
                    legend('desync time', 'dt pre-estimated from pseudo ranges', 'dt pre-estimated from phases', 'dt correction from LS on Code', 'dt estimated from pre-processing', 'residual dt from carrier phases', 'total dt', 'Location', 'NorthEastOutside');
                else
                    legend('desync time', 'dt pre-estimated from pseudo ranges', 'dt pre-estimated from phases', 'dt correction from LS on Code', 'dt estimated from pre-processing', 'Location', 'NorthEastOutside');
                end
                xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('receiver clock error [s]'); h.FontWeight = 'bold';
                h = title(sprintf('dt - receiver %s', rec.parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
            end
        end
        
        function showOutliersAndCycleSlip(this, sys_c_list)
            % Plot the outliers found
            % SYNTAX this.showOutliers()
            
            if nargin == 1
                sys_c_list = unique(this.system);
            end
            f = figure; f.Name = sprintf('%03d: CS, Outlier', f.Number); f.NumberTitle = 'off';
            ss_ok = intersect(this.cc.sys_c, sys_c_list);
            for sys_c = sys_c_list
                ss_id = find(this.cc.sys_c == sys_c);
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
                
                ep = repmat((1: this.time.length)',1, size(this.sat.outlier_idx_ph, 2));
                
                for prn = this.cc.prn(this.cc.system == sys_c)'
                    id_ok = find(any(this.obs((this.system == sys_c)' & this.prn == prn, :),1));
                    plot(id_ok, prn * ones(size(id_ok)), 's', 'Color', [0.8 0.8 0.8]);
                    hold on;
                    id_ok = find(any(this.obs((this.system == sys_c)' & this.prn == prn & this.obs_code(:,1) == 'L', :),1));
                    plot(id_ok, prn * ones(size(id_ok)), '.', 'Color', [0.7 0.7 0.7]);
                    s = find(this.go_id(this.obs_code(:,1) == 'L') == this.getGoId(sys_c, prn));
                    if any(s)
                        cs = ep(this.sat.cycle_slip_idx_ph(:, s) ~= 0);
                        plot(cs,  prn * ones(size(cs)), '.k', 'MarkerSize', 20)
                        out = ep(this.sat.outlier_idx_ph(:, s) ~= 0);
                        plot(out,  prn * ones(size(out)), '.', 'MarkerSize', 20, 'Color', [1 0.4 0]);
                    end
                end
                prn_ss = unique(this.prn(this.system == sys_c));
                xlim([1 size(this.obs,2)]);
                ylim([min(prn_ss) - 1 max(prn_ss) + 1]);
                h = ylabel('PRN'); h.FontWeight = 'bold';
                ax = gca(); ax.YTick = prn_ss;
                grid on;
                h = xlabel('epoch'); h.FontWeight = 'bold';
                h = title(sprintf('%s %s cycle-slip(b) & outlier(o)', this.cc.getSysName(sys_c), this.parent.marker_name), 'interpreter', 'none'); h.FontWeight = 'bold';
            end
        end
        
            function showOutliersAndCycleSlip_p(this, sys_c_list)
            % Plot Signal to Noise Ration in a skyplot
            % SYNTAX this.plotSNR(sys_c)
            
            % SNRs
            if nargin == 1
                sys_c_list = unique(this.system);
            end
            
            for sys_c = sys_c_list
                [~, ~, ph_id] = this.getPhases(sys_c);
                f = figure; f.Name = sprintf('%03d: CS, Out %s', f.Number, sys_c); f.NumberTitle = 'off';
                polarScatter([],[],1,[]);
                hold on;
                decl_n = (serialize(90 - this.sat.el(:, this.go_id(ph_id))) / 180*pi) / (pi/2);
                x = sin(serialize(this.sat.az(:, this.go_id(ph_id))) / 180 * pi) .* decl_n; x(serialize(this.sat.az(:, this.go_id(ph_id))) == 0) = [];
                y = cos(serialize(this.sat.az(:, this.go_id(ph_id))) / 180 * pi) .* decl_n; y(serialize(this.sat.az(:, this.go_id(ph_id))) == 0) = [];
                plot(x, y, '.', 'Color', [0.8 0.8 0.8]); % all sat
                decl_n = (serialize(90 - this.sat.el(this.id_sync(:), this.go_id(ph_id))) / 180*pi) / (pi/2);
                x = sin(serialize(this.sat.az(this.id_sync(:), this.go_id(ph_id))) / 180 * pi) .* decl_n; x(serialize(this.sat.az(this.id_sync(:), this.go_id(ph_id))) == 0) = [];
                y = cos(serialize(this.sat.az(this.id_sync(:), this.go_id(ph_id))) / 180 * pi) .* decl_n; y(serialize(this.sat.az(this.id_sync(:), this.go_id(ph_id))) == 0) = [];
                m = serialize(~isnan(zero2nan(this.obs(ph_id, this.id_sync(:)))'));
                plot(x(m), y(m), '.', 'Color', [0.4 0.4 0.4]); % sat in vew
                decl_n = (serialize(90 - this.sat.el(:, this.go_id(ph_id))) / 180*pi) / (pi/2);
                x = sin(serialize(this.sat.az(:, this.go_id(ph_id))) / 180 * pi) .* decl_n; x(serialize(this.sat.az(:, this.go_id(ph_id))) == 0) = [];
                y = cos(serialize(this.sat.az(:, this.go_id(ph_id))) / 180 * pi) .* decl_n; y(serialize(this.sat.az(:, this.go_id(ph_id))) == 0) = [];
                cut_offed = decl_n(serialize(this.sat.az(:, this.go_id(ph_id))) ~= 0) > ((90 - this.state.getCutOff) / 90);
                plot(x(cut_offed), y(cut_offed), '.', 'Color', [0.95 0.95 0.95]);
                
                for s = unique(this.go_id(ph_id))'
                    az = this.sat.az(:,s);
                    el = this.sat.el(:,s);
                    
                    id_cs = find(this.go_id(this.obs_code(:,1) == 'L') == s);
                    cs = sum(this.sat.cycle_slip_idx_ph(:, id_cs), 2) > 0;
                    out = sum(this.sat.outlier_idx_ph(:, id_cs), 2) > 0;
                    
                    decl_n = (serialize(90 - el(cs)) / 180*pi) / (pi/2);
                    x = sin(az(cs)/180*pi) .* decl_n; x(az(cs) == 0) = [];
                    y = cos(az(cs)/180*pi) .* decl_n; y(az(cs) == 0) = [];
                    plot(x, y, '.k', 'MarkerSize', 20)
                    
                    decl_n = (serialize(90 - el(out)) / 180*pi) / (pi/2);
                    x = sin(az(out)/180*pi) .* decl_n; x(az(out) == 0) = [];
                    y = cos(az(out)/180*pi) .* decl_n; y(az(out) == 0) = [];
                    plot(x, y, '.', 'MarkerSize', 20, 'Color', [1 0.4 0]);
                end
                h = title(sprintf('%s %s cycle-slip(b) & outlier(o)', this.cc.getSysName(sys_c), this.parent.marker_name), 'interpreter', 'none'); h.FontWeight = 'bold';
            end
            
        end
        
        function showDataAvailability(this, sys_c)
            % Plot all the satellite seen by the system
            % SYNTAX this.plotDataAvailability(sys_c)
            
            if (nargin == 1)
                sys_c = this.cc.sys_c;
            end
            f = figure; f.Name = sprintf('%03d: DataAvail', f.Number); f.NumberTitle = 'off';
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
                title(this.cc.getSysName(ss));
            end
        end
        
        function showObsVsSynt_m(this, sys_c)
            % Plots phases and pseudo-ranges aginst their synthesised values
            % SYNTAX
            %   this.plotVsSynt_m(<sys_c>)
            
            % Phases
            if nargin == 1
                [ph, ~, id_ph] = this.getPhases;
                sensor_ph = Core_Pre_Processing.diffAndPred(ph - this.getSyntPhases); sensor_ph = bsxfun(@minus, sensor_ph, median(sensor_ph, 2, 'omitnan'));
                f = figure; f.Name = sprintf('%03d: VsSynt', f.Number); f.NumberTitle = 'off';
            else
                [ph, ~, id_ph] = this.getPhases(sys_c);
                sensor_ph = Core_Pre_Processing.diffAndPred(ph - this.getSyntPhases(sys_c)); sensor_ph = bsxfun(@minus, sensor_ph, median(sensor_ph, 2, 'omitnan'));
                f = figure; f.Name = sprintf('%03d: VsSynt%s', f.Number, this.cc.getSysName(sys_c)); f.NumberTitle = 'off';
            end
            
            subplot(2,3,1); plot(sensor_ph); title('Phases observed vs synthesised');
            
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
        
        function f_handle = showSNR_p(this, sys_c_list, flag_smooth)
            % Plot Signal to Noise Ration in a skyplot
            % SYNTAX f_handles = this.plotSNR(sys_c)
            
            % SNRs
            if nargin < 2 || isempty(sys_c_list)
                sys_c_list = unique(this.system);
            end
            
            idx_f = 0;
            for sys_c = sys_c_list
                for b = 1 : 9 % try all the bands
                    [snr, snr_id] = this.getSNR(sys_c, num2str(b));
                    if nargin > 2 && flag_smooth
                        snr = this.smoothSatData([],[],zero2nan(snr), [], 'spline', 900 / this.getRate, 10); % smoothing SNR => to be improved
                    end
                    
                    if any(snr_id) && any(snr(:))
                        idx_f = idx_f +1;
                        f = figure; f.Name = sprintf('%03d: SNR%d %s', f.Number, b, this.cc.getSysName(sys_c)); f.NumberTitle = 'off';
                        f_handle(idx_f) = f;
                        id_ok = (~isnan(snr));
                        az = this.sat.az(:,this.go_id(snr_id));
                        el = this.sat.el(:,this.go_id(snr_id));
                        polarScatter(serialize(az(id_ok))/180*pi,serialize(90-el(id_ok))/180*pi, 45, serialize(snr(id_ok)), 'filled');
                        colormap(jet);  cax = caxis(); caxis([min(cax(1), 10), max(cax(2), 55)]); setColorMap([10 55], 0.9); colorbar();
                        h = title(sprintf('SNR%d - receiver %s - %s', b, this.parent.marker_name, this.cc.getSysExtName(sys_c)),'interpreter', 'none'); h.FontWeight = 'bold'; h.Units = 'pixels'; h.Position(2) = h.Position(2) + 20; h.Units = 'data';
                    end
                end
            end
            
            if idx_f == 0
                f_handle = [];
            end
        end
    end
end
