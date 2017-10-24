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
        obs_code       % list of types per constellation
        obs_col        % struct matrix containing the colum of an observation .constellation [4(CLDS) x n_freq]
        ph_shift       %
        
        code_C1        % boolean variable to notify if the C1 code is used instead of P1
        xyz;           % approximate position of the receiver (XYZ geocentric)
        
        n_freq = 0;    % number of stored frequencies
        n_max_obs = 0; % maximum numbber of observables
        n_epo = 0;     % number of epochs stored
        n_sat = 0;     % number of satellites
        
        dt = 0;        % clock offset of the receiver
                
        time = [];     % internal time ref of the stored epochs
        rate;          % obbservations rate;
                
        active_ids     % rows of active satellites
        wl             % wave-lenght of each row of row_id
        f_id           % frequency number e.g. L1 -> 1,  L2 ->2, E1 -> 1, E5b -> 3 ...
        prn            % pseudo-range number of the satellite
        go_id          % internal id for a certain satellite
        system         % char id of the satellite system corresponding to the row_id
                
        row_id         % observations id structure
        pr_validity    % validity of the row (does it contains values?)
        ph_validity    % validity of the row (does it contains values?)
        dop_validity   % validity of the row (does it contains values?)
        snr_validity   % validity of the row (does it contains values?)
                 
        pr             % matrix containing pseudo-range observations ordered as SFO (Satellite/Frequency, Observation)
        ph             % matrix containing phase observations ordered as SFO (Satellite/Frequency, Observation)
        snr            % matrix containing snr observations ordered as SFO (Satellite/Frequency, Observation)
        dop            % matrix containing doppler observations ordered as SFO (Satellite/Frequency, Observation)

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
            this.n_freq = 0;    % number of stored frequencies
            this.n_epo = 0;     % number of epochs stored
            this.n_sat = 0;     % number of satellites
            
            this.dt = 0;        % clock offset of the receiver
            this.time = [];     % internal time ref of the stored epochs
            this.rate;          % save rate of the observations
            
            this.active_ids = [];
            this.wl = [];
            this.f_id = [];
            this.system = [];
            
            this.pr = [];       % [n_obs x n_epo] n_obs = number of observables e.g. L1 for sat 1 2 3 + L2 for sat 1 2 3  + E1 for sat 1 2 -> 8 total obs
            this.ph = [];       % [n_obs x n_epo] n_obs = number of observables e.g. L1 for sat 1 2 3 + L2 for sat 1 2 3  + E1 for sat 1 2 -> 8 total obs
            this.snr = [];      % [n_obs x n_epo] n_obs = number of observables e.g. L1 for sat 1 2 3 + L2 for sat 1 2 3  + E1 for sat 1 2 -> 8 total obs
            this.dop = [];      % [n_obs x n_epo] n_obs = number of observables e.g. L1 for sat 1 2 3 + L2 for sat 1 2 3  + E1 for sat 1 2 -> 8 total obs
            
            this.row_id = struct('gps', struct('L1',  (nan(this.cc.gps.N_SAT,1)), ...
                                               'L2',  (nan(this.cc.gps.N_SAT,1)), ...
                                               'L5',  (nan(this.cc.gps.N_SAT,1))), ...
                                 'glo', struct('G1',  (nan(this.cc.glo.N_SAT,1)), ...
                                               'G2',  (nan(this.cc.glo.N_SAT,1)), ...
                                               'G3',  (nan(this.cc.glo.N_SAT,1))), ...
                                 'gal', struct('E1',  (nan(this.cc.gal.N_SAT,1)), ...
                                               'E5a', (nan(this.cc.gal.N_SAT,1)), ...
                                               'E5b', (nan(this.cc.gal.N_SAT,1)), ...
                                               'E5',  (nan(this.cc.gal.N_SAT,1)), ...
                                               'E6',  (nan(this.cc.gal.N_SAT,1))), ...
                                 'bds', struct('B1',  (nan(this.cc.bds.N_SAT,1)), ...
                                               'B2',  (nan(this.cc.bds.N_SAT,1)), ...
                                               'B3',  (nan(this.cc.bds.N_SAT,1))), ...
                                 'qzs', struct('J1',  (nan(this.cc.qzs.N_SAT,1)), ...
                                               'J2',  (nan(this.cc.qzs.N_SAT,1)), ...
                                               'J5',  (nan(this.cc.qzs.N_SAT,1)), ...
                                               'J6',  (nan(this.cc.qzs.N_SAT,1))), ...
                                 'irn', struct('L5',  (nan(this.cc.irn.N_SAT,1)), ...
                                               'S',   (nan(this.cc.irn.N_SAT,1))), ...
                                 'sbs', struct('L1',  (nan(this.cc.sbs.N_SAT,1)), ...
                                               'L5',  (nan(this.cc.sbs.N_SAT,1))));
            this.pr_validity = [];
            this.ph_validity = [];
            this.dop_validity = [];
            this.snr_validity = [];
            
            this.clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true
        end
        function initR2S(this)
            %%% initialize satellite related parameters
            [n_pr, n_epoch]       = size(this.pr1());
            this.rec2sat.avail_index  = not(isnan(this.pr1()'));
            this.cc.cs            = Core_Sky.getInstance();
            this.rec2sat.tot          = NaN(n_epoch, n_pr);
            %  this.rec2sat.XS_tx          = NaN(n_epoch, n_pr); % --> consider what to initialize
        end
        function update(this)
            % update the internal state of the object
            % SYNTAX this.update()
            
            active_ss = this.cc.getActive();
            active_ids = [serialize( active_ss(1) * struct2array(this.row_id.gps) .* repmat((this.cc.gps.flag_f'), this.cc.gps.N_SAT, 1)); ...
                serialize( active_ss(2) * struct2array(this.row_id.glo) .* repmat((this.cc.glo.flag_f'), this.cc.glo.N_SAT, 1)); ...
                serialize( active_ss(3) * struct2array(this.row_id.gal) .* repmat((this.cc.gal.flag_f'), this.cc.gal.N_SAT, 1)); ...
                serialize( active_ss(4) * struct2array(this.row_id.bds) .* repmat((this.cc.bds.flag_f'), this.cc.bds.N_SAT, 1)); ...
                serialize( active_ss(5) * struct2array(this.row_id.qzs) .* repmat((this.cc.qzs.flag_f'), this.cc.qzs.N_SAT, 1)); ...
                serialize( active_ss(6) * struct2array(this.row_id.irn) .* repmat((this.cc.irn.flag_f'), this.cc.irn.N_SAT, 1)); ...
                serialize( active_ss(7) * struct2array(this.row_id.sbs) .* repmat((this.cc.sbs.flag_f'), this.cc.sbs.N_SAT, 1))];
            
            wl = [serialize( repmat(this.cc.gps.L_VEC, this.cc.gps.N_SAT, 1) ); ...
                serialize( this.cc.glo.L_VEC(this.cc.glo.PRN2IDCH,:) ); ...
                serialize( repmat(this.cc.gal.L_VEC, this.cc.gal.N_SAT, 1) ); ...
                serialize( repmat(this.cc.bds.L_VEC, this.cc.bds.N_SAT, 1) ); ...
                serialize( repmat(this.cc.qzs.L_VEC, this.cc.qzs.N_SAT, 1) ); ...
                serialize( repmat(this.cc.irn.L_VEC, this.cc.irn.N_SAT, 1) ); ...
                serialize( repmat(this.cc.sbs.L_VEC, this.cc.sbs.N_SAT, 1) )];
            
            f_id = [serialize( repmat(1:size(this.cc.gps.F_VEC,2), this.cc.gps.N_SAT, 1) ); ...
                serialize( repmat(1:size(this.cc.glo.F_VEC,2), this.cc.glo.N_SAT, 1) ); ...
                serialize( repmat(1:size(this.cc.gal.F_VEC,2), this.cc.gal.N_SAT, 1) ); ...
                serialize( repmat(1:size(this.cc.bds.F_VEC,2), this.cc.bds.N_SAT, 1) ); ...
                serialize( repmat(1:size(this.cc.qzs.F_VEC,2), this.cc.qzs.N_SAT, 1) ); ...
                serialize( repmat(1:size(this.cc.irn.F_VEC,2), this.cc.irn.N_SAT, 1) ); ...
                serialize( repmat(1:size(this.cc.sbs.F_VEC,2), this.cc.sbs.N_SAT, 1) )];

            prn = [serialize( repmat(this.cc.gps.PRN, 1, size(this.cc.gps.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.glo.PRN, 1, size(this.cc.glo.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.gal.PRN, 1, size(this.cc.gal.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.bds.PRN, 1, size(this.cc.bds.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.qzs.PRN, 1, size(this.cc.qzs.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.irn.PRN, 1, size(this.cc.irn.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.sbs.PRN, 1, size(this.cc.sbs.F_VEC, 2)) )];

            go_id = [serialize( repmat(this.cc.gps.go_ids, 1, size(this.cc.gps.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.glo.go_ids, 1, size(this.cc.glo.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.gal.go_ids, 1, size(this.cc.gal.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.bds.go_ids, 1, size(this.cc.bds.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.qzs.go_ids, 1, size(this.cc.qzs.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.irn.go_ids, 1, size(this.cc.irn.F_VEC, 2)) ); ...
                serialize( repmat(this.cc.sbs.go_ids, 1, size(this.cc.sbs.F_VEC, 2)) )];

            system = [char(ones(this.cc.gps.N_SAT * size(this.cc.gps.F_VEC,2), 1) * this.cc.gps.SYS_C); ...
                char(ones(this.cc.glo.N_SAT * size(this.cc.glo.F_VEC,2), 1) * this.cc.glo.SYS_C); ...
                char(ones(this.cc.gal.N_SAT * size(this.cc.gal.F_VEC,2), 1) * this.cc.gal.SYS_C); ...
                char(ones(this.cc.bds.N_SAT * size(this.cc.bds.F_VEC,2), 1) * this.cc.bds.SYS_C); ...
                char(ones(this.cc.qzs.N_SAT * size(this.cc.qzs.F_VEC,2), 1) * this.cc.qzs.SYS_C); ...
                char(ones(this.cc.irn.N_SAT * size(this.cc.irn.F_VEC,2), 1) * this.cc.irn.SYS_C); ...
                char(ones(this.cc.sbs.N_SAT * size(this.cc.sbs.F_VEC,2), 1) * this.cc.sbs.SYS_C)]';
            
            active_ids(active_ids == 0) = NaN;
            is_active = ~isnan(active_ids);
            active_ids = active_ids(is_active);
            wl = wl(is_active);
            wl(active_ids) = wl;
            f_id = f_id(is_active);
            f_id(active_ids) = f_id;
            prn = prn(is_active);
            prn(active_ids) = prn;
            go_id = go_id(is_active);
            go_id(active_ids) = go_id;
            system = system(is_active);
            system(active_ids) = system;
            
            this.wl = wl;
            this.f_id = f_id;
            this.prn = prn;
            this.go_id = go_id;
            this.n_freq = numel(unique(f_id));
            this.system = system;
            
            this.active_ids = sort(active_ids);
            this.pr_validity = sum(nan2zero(abs(this.ph(this.active_ids, :))), 2) > 0;
            this.ph_validity = sum(nan2zero(abs(this.pr(this.active_ids, :))), 2) > 0;
            this.dop_validity = sum(nan2zero(abs(this.dop(this.active_ids, :))), 2) > 0;
            this.snr_validity = sum(nan2zero(abs(this.snr(this.active_ids, :))), 2) > 0;
        end
        
        function [active_ids, active_ss] = getIdList(this)
            % get list of active observations rows
            % SYNTAX: [active_ids, active_ss] = getIdList(this)
            active_ss = this.cc.getActive();
            if isempty(this.active_ids)
                this.update();
            end
            active_ids = this.active_ids;
        end
        
        function [full_ids] = getFullIdList(this)
            % get list of all the observations rows stored into the object
            % SYNTAX: [full_ids] = getFullIdList(this)
            full_ids = [serialize(struct2array(this.row_id.gps)); ...
                        serialize(struct2array(this.row_id.glo)); ...
                        serialize(struct2array(this.row_id.gal)); ...
                        serialize(struct2array(this.row_id.bds)); ...
                        serialize(struct2array(this.row_id.qzs)); ...
                        serialize(struct2array(this.row_id.sbs))];
            full_ids = full_ids(~isnan(full_ids));
        end
        
        function legacyImport(this, time, pos, dt, clock_corrected_obs, pr1, ph1, pr2, ph2, snr1, snr2, dop1, dop2)
            % import with OSF (Observations Satellites Frequencies)
            % SYNTAX:   
            %   this.legacyImport(this, dt, is_clock_corrected, pr1, ph1, pr2, ph2)
            % EXAMPLE:
            %   this.legacyImport(time_GPS_diff, pos_R, 0, 0, pr1_R, ph1_R, pr2_R, ph2_R, snr1_R, snr2_R, dop1_R, dop2_R);
            % REMARKS:  pr1, ph1, pr2, ph2 must have the same size and are stored according to cc (constellation collector)
            
            this.initObs();
            
            this.time = time;
            this.pos = pos;
            this.dt = dt;
            this.clock_corrected_obs = clock_corrected_obs;
            this.n_freq = 1;
                        
            if (any(pr2(:)) || any(ph2(:)))
                this.n_freq = 2;
                this.pr = zero2nan([pr1; pr2]);
                this.ph = zero2nan([ph1; ph2]);
                this.dop = zero2nan([dop1; dop2]);
                this.snr = zero2nan([snr1; snr2]);
            else
                this.pr = zero2nan(pr1);
                this.ph = zero2nan(ph1);
                this.dop = zero2nan(dop1);
                this.snr = zero2nan(snr1);
            end
            [this.n_sat, this.n_epo] = size(pr1);
            
            for s = 1 : this.cc.n_sys
                ss = this.cc.sys_c(s);
                sys_pos = find(this.cc.system == ss)';
                switch ss
                    case 'G'
                        this.row_id.gps.L1 = sys_pos;
                        if this.n_freq == 2
                            this.row_id.gps.L2 = sys_pos + numel(this.cc.system);
                        end
                    case 'R'
                        this.row_id.glo.G1 = sys_pos;
                        if this.n_freq == 2
                            this.row_id.glo.G2 = sys_pos + numel(this.cc.system);
                        end
                    case 'E'
                        this.row_id.gal.E1 = sys_pos;
                        if this.n_freq == 2
                            this.row_id.gal.E5a = sys_pos + numel(this.cc.system);
                        end
                    case 'C'
                        this.row_id.bds.B1 = sys_pos;
                        if this.n_freq == 2
                            this.row_id.bds.B2 = sys_pos + numel(this.cc.system);
                        end
                    case 'J'
                        this.row_id.qzs.J1 = sys_pos;
                        if this.n_freq == 2
                            this.row_id.qzs.J2 = sys_pos + numel(this.cc.system);
                        end
                    case 'I'
                        this.row_id.irn.L5 = sys_pos;
                        if this.n_freq == 2
                            this.row_id.irn.S = sys_pos + numel(this.cc.system);
                        end
                    case 'S'
                        this.row_id.sbs.L1 = sys_pos;
                        if this.n_freq == 2
                            this.row_id.sbs.L5 = sys_pos + numel(this.cc.system);
                        end
                end
            end
            this.update();
        end
        
        function remObs(this, id_obs)
            % remove observations with a certain id
            % SYNTAX:   this.remObs(id_obs)
            this.ph(:, id_obs) = [];
            this.pr(:, id_obs) = [];
            this.dop(:, id_obs) = [];
            this.snr(:, id_obs) = [];
            this.dt(id_obs) = [];
            this.time(id_obs) = [];
            this.n_epo = size(this.pr, 2);
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
                if strcmp(txt((lim(l,1) + 60) : lim(l,2)), h_std{1})
                    type_found = true;
                    dataset = textscan(txt(lim(1,1):lim(1,2)), '%f%c%18c%c');
                end
            end
            this.rin_type = dataset{1};
            if dataset{2} == 'O'
                if (this.rin_type < 3)
                    if (dataset{4} ~= 'G')
                        % GPS only RINEX2 - mixed or glonass -> actually not working 
                        throw(MException('VerifyInput:InvalidObservationFile', 'RINEX2 is supported for GPS only dataset, please use a RINEX3 file '));
                    else
                        % GPS only RINEX2 -> ok
                    end
                else
                    % RINEX 3 file -> ok
                end
            else
                throw(MException('VerifyInput:InvalidObservationFile', 'This observation RINEX does not contain obbservations'));
            end
                        
            % parsing ------------------------------------------------------------------------------------------------------------------------------------------
            
            % retriving the kind of header information is contained on each line
            line2head = zeros(eoh, 1);
            l = 0;
            while l < eoh
                l = l + 1;
                txt((lim(l,1) + 60) : lim(l,2));
                tmp = find(strcmp(txt((lim(l,1) + 60) : lim(l,2)), head_field));
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
            % 14) LEAP SECONDS
            % ignoring
            % 15) # OF SATELLITES
            fln = find(line2head == 15, 1, 'first'); % get field line
            if isempty(fln)
                this.n_sat = this.cc.getNumSat(); % If it's zero it'll be necessary to compute it
            else
                tmp = sscanf(txt(lim(fln, 1) + (0:5)),'%f')';                                  % read value
                this.n_sat = iif(isempty(tmp) || ~isnumeric(tmp) || (numel(tmp) ~= 1), cc.getNumSat(), tmp);  % check value integrity
            end
            % 16) PRN / # OF OBS
            if this.rin_type < 3 % if it's RINEX 2 consider gps only
                fln = find(line2head == 16); % get field line                
                this.obs_code = struct('g',[]);
                if ~isempty(fln)
                    n_obs = sscanf(txt(lim(fln(1), 1) + (0:5)),'%d');
                    for l = 1 : numel(fln)
                        n_obs = sscanf(txt(lim(fln(l), 1) + (0:5)),'%c');
                        this.obs_code.g = [this.obs_code.g sscanf(txt(lim(fln(l), 1) + (6:59)),'%s')];
                    end
                    if (n_obs > (numel(this.obs_code{1}) / 2))
                        error('In Receiver load RINEX, something bad appened: reading a rinex 2 the number of observation types is different from the number of types that have been read in the header');
                    end
                    this.obs_code.g = serialize([reshape(this.obs_code.g, 2, n_obs); 32 * ones(1, n_obs)])';                    
                end              
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
                this.obs_code = struct('g',[],'r',[],'e',[],'j',[],'c',[],'i',[],'s',[]);
                if ~isempty(fln)
                    
                    l = 1;
                    while l <= numel(fln)
                        sys = char(txt(lim(fln(l), 1))+32);
                        n_obs = sscanf(txt(lim(fln(l), 1) + (3:5)),'%d');
                        n_line = ceil(n_obs / 13);
                        l_offset = 0;
                        while l_offset < n_line
                            this.obs_code.(sys) = [this.obs_code.(sys) sscanf(txt(lim(fln(l + l_offset), 1) + (6:59)),'%s')];
                            l_offset = l_offset + 1;
                        end
                        l = l + l_offset;                    
                    end
                end
            end
            % 20) SYS / PHASE SHIFT
            fln = find(line2head == 20); % get field line
            if this.rin_type < 3
                this.ph_shift = struct('g', zeros(numel(this.obs_code.g) / 3, 1));
            else
                this.ph_shift = struct('g',[],'r',[],'e',[],'j',[],'c',[],'i',[],'s',[]);
                for l = 1 : numel(fln)                  
                    sys = char(txt(lim(fln(l), 1)) + 32);
                    obs_code = txt(lim(fln(l), 1) + (2:4));
                    obs_id = (strfind(this.obs_code.(sys), obs_code) - 1) / 3 + 1;
                    if isempty(this.ph_shift.(sys))
                        this.ph_shift.(sys) = zeros(numel(this.obs_code.(sys)) / 3, 1);
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
            
            this.chooseDataTypes();
        end
        
        function chooseDataTypes(this)
            % get the right attribute column to be used for a certain type/band couple
            t_ok = 'CLDS'; % type
            
            this.obs_col = struct('g', zeros(4, numel(this.cc.gps.F_VEC)), ...
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
                    
                    if ~isempty(this.obs_code.g)
                        code = reshape(this.obs_code.(sys_c), 3, numel(this.obs_code.(sys_c)) / 3)';
                        b_ok = this.cc.(sys).CODE_RIN3_2BAND;  % band
                        a_ok = this.cc.(sys).CODE_RIN3_ATTRIB; % attribute
                        for t = 1 : numel(t_ok)
                            for b = 1 : numel(b_ok)
                                % get the observation codes with a certain type t_ok(t) and band b_ok(b)
                                obs = (code(:,1) == t_ok(t)) & (code(:,2) == b_ok(b));
                                if any(obs)
                                    % find the preferred observation among the available ones
                                    [a, id] = intersect(code(obs, 3), a_ok{b}); a = a(id);
                                    % save the id of the column in the obs_col struct matrix
                                    this.obs_col.(sys_c)(t, b) = find(obs & code(:,3) == a(1));
                                end
                            end
                        end
                    end
                end
                
            else % rinex 2
                keyboard;
                % to be done
            end
            
            this.n_max_obs = sum(this.obs_col.g > 0, 2) * this.cc.gps.N_SAT + ...
                             sum(this.obs_col.r > 0, 2) * this.cc.glo.N_SAT + ...
                             sum(this.obs_col.e > 0, 2) * this.cc.gal.N_SAT + ...
                             sum(this.obs_col.j > 0, 2) * this.cc.qzs.N_SAT + ...
                             sum(this.obs_col.c > 0, 2) * this.cc.bds.N_SAT + ...
                             sum(this.obs_col.i > 0, 2) * this.cc.irn.N_SAT + ...
                             sum(this.obs_col.s > 0, 2) * this.cc.sbs.N_SAT;
                         
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
                        
            processing_interval = 0;
            
            % number of satellite slots for enabled constellations
            nObsTot = cc.getMaxObsSat();
                        
            %variable initialization
            this.file =  File_Rinex(file_name, 9);
            
            if this.file.isValid()
                this.logger.addMessage(sprintf('Opening file %s for reading', file_name));
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

                % importing header informations
                this.parseRinHead(txt, lim, fr.eoh);
                
                % considering rinex 3
                if (this.rin_type >= 3)
                    % find all the obbservation lines
                    t_line = find([false(eoh, 1); (txt(lim(eoh+1:end,1)) == '>')']);
                    % extract all the epoch lines
                    string_time = txt(repmat(lim(t_line,1),1,27) + repmat(2:28, this.n_epo, 1))';
                    % convert the times into a 6 col time
                    date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %10.7f'));
                    % import it as a GPS_Time obj
                    this.time = GPS_Time(date, [], this.file.first_epoch.is_gps);
                    this.n_epo = numel(t_line);
                    
                    % get number of observations per epoch
                    n_ope = sscanf(txt(repmat(lim(t_line,1),1,3) + repmat(32:34, this.n_epo, 1))', '%d');
                    for e = this.n_epo
                        txt(repmat(lim(t_line,1),1,27) + repmat(30:32, this.n_epo, 1))'                                                
                    end
                    
                    % .... init datasets
                end
                                
                pr = NaN(nObsTot, n_epo);
                ph = NaN(nObsTot, n_epo);
                dop = NaN(nObsTot,n_epo);
                snr = NaN(nObsTot,n_epo);
                date = NaN(n_epo, 6);
                pos = zeros(3,1);
                interval = 0;
                antoff = zeros(3,1);
                antmod = '';
                codeC1 = zeros(nObsTot,n_epo);
                marker = '';
                
                % to be used in get obs:
                % starting index in the total array for the various constellations
                idGPS = this.cc.getGPS().getFirstId();
                idGLONASS = this.cc.getGLONASS().getFirstId();
                idGalileo = this.cc.getGalileo().getFirstId();
                idBeiDou = this.cc.getBeiDou().getFirstId();
                idQZSS = this.cc.getQZSS().getFirstId();
                idSBAS = this.cc.getSBAS().getFirstId();
                % to speed-up get obs in MATLAB I need to precompute these values
                first_id_sys = [idGPS idGLONASS idGalileo idBeiDou idQZSS idSBAS];
                active_sys = cc.getActive();
                
                max_k = 0;                               
                
                % Let's start parsing the header
                
                %file_buf = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '');
                %file_buf = file_buf{1};
                %cur_line = 0;
                %fclose(fid);
                
                this.logger.addMessage('Start parsing (this operation could last several seconds)');
                
                %parse RINEX header
                [cur_line, obs_type, pos(:,1,f), basic_info, interval(1,1,f), sysId, antoff(:,1,f), antmod{1,1,f}, marker{1,1,f}] = RINEX_parse_hdr(file_buf, cur_line);
                
                %check the availability of basic data to parse the RINEX file
                if (basic_info == 0)
                    error(['RINEX file ' file_name ': basic data is missing in the file header'])
                end
                
                % find observation type columns
                [obsColumns, nObsTypes] = obs_type_find(obs_type, sysId);
                
                %-------------------------------------------------------------------------------
                
                k = 1;
                while (cur_line < numel(file_buf))
                    
                    %read data for the current epoch (ROVER)
                    [cur_line, date(k,:,f), num_sat, sat, sat_types] = RINEX_get_epoch(file_buf, cur_line);
                    if ~isempty(date(k,1,f))
                        if (k > n_epo)
                            %variable initialization (GPS)
                            pr1(:,k,f) = zeros(nObsTot,1);
                            pr2(:,k,f) = zeros(nObsTot,1);
                            ph1(:,k,f) = zeros(nObsTot,1);
                            ph2(:,k,f) = zeros(nObsTot,1);
                            dop1(:,k,f) = zeros(nObsTot,1);
                            dop2(:,k,f) = zeros(nObsTot,1);
                            snr1(:,k,f) = zeros(nObsTot,1);
                            snr2(:,k,f) = zeros(nObsTot,1);
                            
                            n_epo = n_epo  + 1;
                        end
                        
                        %read ROVER observations
                        [cur_line, obs] = RINEX_get_obs(file_buf, cur_line, num_sat, sat, sat_types, obsColumns, nObsTypes, cc, active_sys, first_id_sys);
                        
                        idx_P1 = obs.P1 ~= 0;
                        idx_C1 = obs.C1 ~= 0;
                        idx_codeC1 = idx_P1 - idx_C1;
                        codeC1(idx_codeC1 < 0,k,f) = 1;
                        pr1(:,k,f) = zeros(size(pr1(:,k,f)));
                        pr1(idx_P1,k,f) = obs.P1(idx_P1);
                        pr1(find(codeC1(:,k,f)),k,f) = obs.C1(find(codeC1(:,k,f))); %#ok<FNDSB>
                        
                        %         %read ROVER observations
                        %         if (~any(obs.C1) || sum(obs.P1 ~= 0) == sum(obs.C1 ~= 0))
                        %             pr1(:,k,f) = obs.P1;
                        %         else
                        %             pr1(:,k,f) = obs.C1;
                        %             codeC1(:,:,f) = 1;
                        %         end
                        pr2(:,k,f) = obs.P2;
                        ph1(:,k,f) = obs.L1;
                        ph2(:,k,f) = obs.L2;
                        dop1(:,k,f) = obs.D1;
                        dop2(:,k,f) = obs.D2;
                        snr1(:,k,f) = obs.S1;
                        snr2(:,k,f) = obs.S2;
                        k = k + 1;
                    end
                end
                
                max_k = max(max_k, k-1);
                
                time(f) = GPS_Time(date(1:k-1,:,f));
                % try to guess observation rate when not read from header
                if (interval(1,1,f) == 0)
                    interval(1,1,f) = time(f).getRate();
                end
                
                %-------------------------------------------------------------------------------
                
                %     if (processing_interval > interval(:,1,f))
                %         interval(:,1,f) = processing_interval;
                %     end
                
                if (~isempty(report) && report.opt.write == 1)
                    % extract quality parameters for report
                    j=strfind(file_name,'\');
                    if isempty(j)
                        j=strfind(file_name,'/');
                    end
                    if isempty(j)
                        j=0;
                    end
                    report.obs.filename(f)=cellstr(file_name(j(end)+1:end));
                    % create statistics on observations
                    stat_sat = ((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f))) + (ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ...
                        (pr1(:,:,f)~=0 & isfinite(pr1(:,:,f))) + (pr2(:,:,f)~=0 & isfinite(pr2(:,:,f))) + ...
                        (dop1(:,:,f)~=0 & isfinite(dop1(:,:,f))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
                    report.obs_raw.n_sat(f)=sum(sum(stat_sat,2)~=0);
                    report.obs_raw.n_epoch(f)=sum(sum(stat_sat,1)~=0);
                    report.obs_raw.n_ph1(f) = sum(sum((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f)))));
                    report.obs_raw.n_ph2(f) = sum(sum((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f)))));
                    report.obs_raw.n_pr1(f) = sum(sum((pr1(:,:,f)~=0 & isfinite(pr1(:,:,f)))));
                    report.obs_raw.n_pr2(f) = sum(sum((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))));
                    report.obs_raw.n_dop1(f) = sum(sum((dop1(:,:,f)~=0 & isfinite(dop1(:,:,f)))));
                    report.obs_raw.n_dop2(f) = sum(sum((dop2(:,:,f)~=0 & isfinite(dop2(:,:,f)))));
                    report.obs_raw.interval(f) = interval(1,1,f);
                    report.obs_raw.time_start(f)=cellstr(sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(1,1,f),date(1,2,f),date(1,3,f),date(1,4,f),date(1,5,f),date(1,6,f)));
                    report.obs_raw.time_end(f)=cellstr(sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(k-1,1,f),date(k-1,2,f),date(k-1,3,f),date(k-1,4,f),date(k-1,5,f),date(k-1,6,f)));
                    
                    stat_sat = ((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
                    if any(stat_sat(:))
                        report.obs_raw.nfreq(f)=2;
                    else
                        report.obs_raw.nfreq(f)=1;
                    end
                    report.obs_raw.n_epoch_expected(f) = time(f).getExpectedLen();
                    
                    report.obs_raw.epoch_completeness(f)=report.obs_raw.n_epoch(f)/report.obs_raw.n_epoch_expected(f)*100;
                    if report.obs_raw.nfreq(f) == 2
                        report.obs_raw.L1L2_completeness(f) = cellstr(sprintf('%6.1f', report.obs_raw.n_ph2(f)/report.obs_raw.n_ph1(f)*100));
                    else
                        report.obs_raw.L1L2_completeness(f) = cellstr(sprintf('%6s', '0'));
                    end
                end
                
                this.logger.addMessage('RINEX parsing completed');
                this.logger.newLine();
                
                % trim output (it have been pre-allocated bigger)
                pr1 = pr1(:,(1 : max_k),:);
                pr2 = pr2(:,(1 : max_k),:);
                ph1 = ph1(:,(1 : max_k),:);
                ph2 = ph2(:,(1 : max_k),:);
                dop1 = dop1(:,(1 : max_k),:);
                dop2 = dop2(:,(1 : max_k),:);
                snr1 = snr1(:,(1 : max_k),:);
                snr2 = snr2(:,(1 : max_k),:);
                date = date((1 : max_k),:,:);
                codeC1 = codeC1(:,(1 : max_k),:);
                
                this.logger.addMessage('Syncing observations if needed');
                
                %sync observations
                [time_zero, time_GPS, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, codeC1, interval] = ...
                    sync_obs(time, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, codeC1, interval, processing_interval);
                state.p_rate = interval;
                
                this.logger.addMessage('Trimming short arcs if needed');
                % remove short arcs
                min_arc = state.getMinArc();
                for f = 1 : nFiles
                    pr1(:,:,f) = remove_short_arcs(pr1(:,:,f), min_arc);
                    pr2(:,:,f) = remove_short_arcs(pr2(:,:,f), min_arc);
                    ph1(:,:,f) = remove_short_arcs(ph1(:,:,f), min_arc);
                    ph2(:,:,f) = remove_short_arcs(ph2(:,:,f), min_arc);
                end
                
                % Find when all the dataset have at least one good observations in common
                resync_flag_ok = any(pr1,1) & any(pr1,1) & any(pr1,1) & any(ph1,1);
                max_sync = find(sum(resync_flag_ok,3) == size(resync_flag_ok,3), 1, 'last');
                
                pr1 = pr1(:,(1 : max_sync),:);
                pr2 = pr2(:,(1 : max_sync),:);
                ph1 = ph1(:,(1 : max_sync),:);
                ph2 = ph2(:,(1 : max_sync),:);
                dop1 = dop1(:,(1 : max_sync),:);
                dop2 = dop2(:,(1 : max_sync),:);
                snr1 = snr1(:,(1 : max_sync),:);
                snr2 = snr2(:,(1 : max_sync),:);
                date = date((1 : max_sync),:,:);
                codeC1 = codeC1(:,(1 : max_sync),:);
                time_GPS = time_GPS(1 : max_sync);
                time = time((1 : max_sync),:,:);
                week = week((1 : max_sync),:,:);
                date = date((1 : max_sync),:,:);
                
                for f = 1 : nFiles
                    holes = find(week(:,1,f) == 0);
                    for h = holes'
                        if (h > 1)
                            time(h,:,f) = time(h-1,1,f) + interval;
                            week(h,1,f) = week(h-1,1,f);
                            date(h,:,f) = datevec(datenum(date(h-1,:,f)) + datenum([0 0 0 0 0 interval]));
                        elseif (holes(end)+1 <= length(week(:,1,f)))
                            time(h,1,f) = time(holes(end)+1,1,f) - interval*holes(end);
                            week(h,1,f) = week(holes(end)+1,1,f);
                            date(h,:,f) = datevec(datenum(date(holes(end)+1,:,f)) - datenum([0 0 0 0 0 interval*holes(end)]));
                        end
                    end
                end
                
                if (~isempty(report) && report.opt.write == 1)
                    % extract quality parameters for report
                    for f = 1 : nFiles
                        % create statistics on observations
                        stat_sat = ((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f))) + (ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ...
                            (pr1(:,:,f)~=0 & isfinite(pr1(:,:,f))) + (pr2(:,:,f)~=0 & isfinite(pr2(:,:,f))) + ...
                            (dop1(:,:,f)~=0 & isfinite(dop1(:,:,f))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
                        report.obs_sync.n_sat(f)=sum(sum(stat_sat,2)~=0);
                        report.obs_sync.n_epoch(f)=sum(sum(stat_sat,1)~=0);
                        report.obs_sync.n_ph1(f) = sum(sum((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f)))));
                        report.obs_sync.n_ph2(f) = sum(sum((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f)))));
                        report.obs_sync.n_pr1(f) = sum(sum((pr1(:,:,f)~=0 & isfinite(pr1(:,:,f)))));
                        report.obs_sync.n_pr2(f) = sum(sum((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))));
                        report.obs_sync.n_dop1(f) = sum(sum((dop1(:,:,f)~=0 & isfinite(dop1(:,:,f)))));
                        report.obs_sync.n_dop2(f) = sum(sum((dop2(:,:,f)~=0 & isfinite(dop2(:,:,f)))));
                        report.obs_sync.interval(f) = interval;
                        report.obs_sync.time_start(f)=cellstr(sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(1,1,f),date(1,2,f),date(1,3,f),date(1,4,f),date(1,5,f),date(1,6,f)));
                        report.obs_sync.time_end(f)=cellstr(sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(size(date,1),1,f),date(size(date,1),2,f),date(size(date,1),3,f),date(size(date,1),4,f),date(size(date,1),5,f),date(size(date,1),6,f)));
                        
                        stat_sat = ((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
                        if any(stat_sat(:))
                            report.obs_sync.nfreq(f)=2;
                        else
                            report.obs_sync.nfreq(f)=1;
                        end
                        report.obs_sync.n_epoch_expected(f) = length((roundmod(time(1,1,f),interval) : interval : roundmod(time(size(date,1),1,f),interval)));
                        
                        report.obs_sync.epoch_completeness(f)=report.obs_sync.n_epoch(f)/report.obs_sync.n_epoch_expected(f)*100;
                        if report.obs_sync.nfreq(f) == 2
                            report.obs_sync.L1L2_completeness(f) = cellstr(sprintf('%6.1f', report.obs_sync.n_ph2(f)/report.obs_sync.n_ph1(f)*100));
                        else
                            report.obs_sync.L1L2_completeness(f) = cellstr(sprintf('%6s', '0'));
                        end
                    end
                end
                
                this.logger.addMessage(sprintf('Parsing completed in %.2f seconds', toc(t0)));
                this.logger.newLine();
            end
        end
    end
    
    % ==================================================================================================================================================
    %  GETTER
    % ==================================================================================================================================================
    methods
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
            time_tx = this.time - this.tot(this.rec2sat.avail_index)
            
            
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
end
