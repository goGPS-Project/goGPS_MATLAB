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

    properties (Constant)
    end

    properties (SetAccess = public, GetAccess = public)
        h_antenna = 0; % Antenna height from the ground [m]

        pos = [0 0 0]; % approximate position of the receiver (XYZ geocentric)
        
        n_freq = 0;    % number of stored frequencies
        n_epo = 0;     % number of epochs stored
        n_sat = 0;     % number of satellites
        
        dt = 0;        % clock offset of the receiver
                
        time = [];     % internal time ref of the stored epochs
        
        cc = Constellation_Collector('GRECJ'); % constellation collector
        
        active_ids  % rows of active satellites
        wl          % wave-lenght of each row of row_id
        f_id        % frequency number e.g. L1 -> 1,  L2 ->2, E1 -> 1, E5b -> 3 ...
        prn         % pseudo-range number of the satellite
        go_id       % internal id for a certain satellite
        system      % char id of the satellite system corresponding to the row_id
                
        row_id         % observations id structure
        pr_validity    % validity of the row (does it contains values?)
        ph_validity    % validity of the row (does it contains values?)
        dop_validity   % validity of the row (does it contains values?)
        snr_validity   % validity of the row (does it contains values?)
                 
        pr          % matrix containing pseudo-range observations ordered as SFO (Satellite/Frequency, Observation)
        ph          % matrix containing phase observations ordered as SFO (Satellite/Frequency, Observation)
        snr         % matrix containing snr observations ordered as SFO (Satellite/Frequency, Observation)
        dop         % matrix containing doppler observations ordered as SFO (Satellite/Frequency, Observation)

        clock_corrected_obs = false; % if the obs have been corrected with dt * v_light this flag should be true
    end

    % ==================================================================================================================================================
    %  SETTER
    % ==================================================================================================================================================
    
    methods
        function this = Receiver(cc)
            % SYNTAX  this = Receiver(cc)
            this.cc = cc;
            this.initObs();
        end
        
        function initObs(this)
            % initialize the receiver obj
            this.n_freq = 0;    % number of stored frequencies
            this.n_epo = 0;     % number of epochs stored
            this.n_sat = 0;     % number of satellites
            
            this.dt = 0;        % clock offset of the receiver
            this.time = [];     % internal time ref of the stored epochs
            
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
end
