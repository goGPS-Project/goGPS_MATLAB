%   CLASS Core Block
% =========================================================================
%
% DESCRIPTION
%   Class to manage goBlock solutions
%
% EXAMPLE
%   go_block = Core_Block();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_Block
%
% Note for the future: the class uses the current obs storage of goGPS
% -> switch to objects for rover and master observations is suggested

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef Core_Block < handle

    properties (Constant, Access = private)
        FLAG_CODE_ONLY  = int8(-1);
        FLAG_CODE_PHASE = int8(0);
        FLAG_PHASE_ONLY = int8(1);
    end

    properties (Access = public)% Public Access
        log
        state

        % number of position solutions to be estimated
        n_pos = 1

        % number of position solutions to be estimated (high rate)
        n_pos_hr = 1

        % solution rate used for the high rate solution
        s_rate = 86400;

        % flag to spacify the tipe of solution code only (-1) - code/phase (0) - phase only (1)
        sol_type = Core_Block.FLAG_PHASE_ONLY

        % total number of observations that can be used
        n_obs_tot = 0

        % number of valid epochs used in goBlock
        n_epoch = 0

        % max number of valid epochs used in goBlock
        n_tot_epoch = 0

        % reference time
        time_diff

        % indexes of empty observations
        empty_epoch = []

        % matrices to keep track of the satellite configuration changes (to fill in the proper ambiguity slots)
        sat_pr_track    % satellite configuration matrix for code observations
        sat_ph_track    % satellite configuration matrix for phase observations
        pivot_track     % satellite configuration matrix for pivot tracking

        obs_track       % matrix to keep track of the obs -> epoch; PRN; flag code/phase;

        amb_prn_track   % store the prn of each ambiguity
        ref_arc         % arcs used to stabilize the solution (the system is rank deficient)
        col_ok          % columns of the design matrix to be used for the LS

        % LS variable
        A   % LS Design matrix
        Q   % LS Cofactor matrix
        y0  % LS Observations array
        b   % LS known array

        % Results
        pos0            % a-priori position
        pos             % estimated position
        pos_cov         % estimated covariance matrix of the positions
        is_fixed = 0    % 0 => float 1 => fix 2 => fix_hr

        x_float         % estimated parameter s(float solution)
        x_fix           % estimated parameters (fix solution)
        amb_fix         % ambiguities as fixed by lambda
        amb_fix_full    % full set of ambiguities
        G               % transformation matrix SD -> DD
        pos_cov_fix     % Covariance matrix of the fixed positions
        x_hr            % estimated parameters (high_rate)
        Cxx             % Covariance matrix of the parameters
        s02             % estimated variance
        v_hat           % residuals of the observarions

        phase_res       % phase residuals ([n_obs x n_amb x 2]);  first slice value, second slice weight
        id_track        % id in the design matrix of the observations in phase_res
        
        cs_factor = 0.5 % fix cs_factor * lambda jumps
    end
    
    % ==================================================================================================================================================
    %  CREATOR
    % ==================================================================================================================================================
    
    methods (Static)
        function this = Core_Block(n_epoch, n_pr_obs, n_ph_obs)
            % Core object creator initialize the structures needed for the computation:
            % EXAMPLE: go_block = Core_Block(n_epoch, n_pr_obs, n_ph_obs)
            
            this.log = Core.getLogger();
            this.state = Core.getState();
            
            % number of position solutions to be estimated
            this.n_pos = 1;
            this.n_pos_hr = 1;
            
            this.n_epoch = n_epoch;
            this.n_tot_epoch = n_epoch;
            this.time_diff = [];
            
            n_sat = this.state.cc.getNumSat();
            % matrices to keep track of the satellite configuration changes (to fill in the proper ambiguity slots)
            this.sat_pr_track = int8(zeros(n_sat, n_epoch));
            this.sat_ph_track = int8(zeros(n_sat, n_epoch));
            this.pivot_track = uint8(zeros(n_epoch, 1));
            
            % total number of observations (for matrix initialization)
            % (some satellites observations will be discarded -> this is the max size)
            if (this.state.isModePh())
                % n_obs_tot = n_pr_obs + n_ph_obs; % to use code and phase
                this.n_obs_tot = n_ph_obs;              % to only use phase
            else
                this.n_obs_tot = n_pr_obs;
            end
            
            this.obs_track = NaN(this.n_obs_tot, 3); % epoch; PRN; flag code/phase
            this.empty_epoch = [];
            
            % init LS variables
            this.y0 = NaN(this.n_obs_tot, 1);
            this.b  = NaN(this.n_obs_tot, 1);
            this.A = zeros(this.n_obs_tot, this.n_pos * 3);
            this.Q  = spalloc(this.n_obs_tot, this.n_obs_tot, 30 * this.n_obs_tot); % considering 30 satellites per epoch as mean
        end
    end
    
    % ==================================================================================================================================================
    %  PROCESSING FUNCTIONS goBlock
    % ==================================================================================================================================================
    
    methods % Public Access
        function prepare (this, ...
                time_diff,  ...
                pos_r, pos_m,  ...
                pr1_r, pr1_m, pr2_r, pr2_m, ...
                ph1_r, ph1_m, ph2_r, ph2_m,  ...
                snr_r, snr_m,  ...
                eph, sp3, iono, lambda, ant_pcv)
            
            % Fill the matrices of the LS system, this is necessary to get a solution
            %
            % SYNTAX:
            %   prepare(this, time_diff, pos_r, pos_m, pr1_r, pr1_m, pr2_r, pr2_m, ph1_r, ph1_m, ph2_r, ph2_m, snr_r, snr_m, eph, sp3, iono, lambda, ant_pcv)
            %
            % INPUT:
            %   time_diff  GPS reception time
            %   pos_r      ROVER approximate position
            %   pos_m      MASTER position
            %   pr1_r      ROVER code observations (L1 carrier)
            %   pr1_m      MASTER code observations (L1 carrier)
            %   pr2_r      ROVER code observations (L2 carrier)
            %   pr2_m      MASTER code observations (L2 carrier)
            %   ph1_r      ROVER phase observations (L1 carrier)
            %   ph1_m      MASTER phase observations (L1 carrier)
            %   ph2_r      ROVER phase observations (L2 carrier)
            %   ph2_m      MASTER phase observations (L2 carrier)
            %   snr_r      ROVER-SATELLITE signal-to-noise ratio
            %   snr_m      MASTER-SATELLITE signal-to-noise ratio
            %   eph        satellite ephemeris
            %   sp3        structure containing precise ephemeris and clock
            %   iono       ionosphere parameters
            %   lambda     wavelength matrix (depending on the enabled constellations)
            %   phase      L1 carrier (phase=1), L2 carrier (phase=2)
            %   ant_pcv antenna phase center variation
            %
            % INTERNAL INPUT:
            %   state
            %
            % INTERNAL OUTOUT (properties whos values are changed):
            %   log, pos, pos0, pos_cov, is_fixed, x_float, x_fix, x_hr, Cxx, s02, v_hat,
            %   sat_pr_track, sat_ph_track, obs_track
            %   y0, b, A, Q
            %
            % CALL:
            %   oneEpochLS
            %   addAmbiguities
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            
            this.log.addMarkedMessage('Preparing goBlock system');
            
            this.time_diff = time_diff; % reference time
            this.pos = [];          % estimated parameters
            this.pos0 = pos_r;      % a-priori position
            this.pos_cov = [];      % estimated parameters
            this.is_fixed = 0;      % flag is fixed
            
            this.x_float = [];      % estimated float parameters
            this.x_fix = [];        % estimated parameters (fix solution)
            this.x_hr = [];         % estimated parameters (fix solution + float)
            this.Cxx = [];          % Covariance matrix of the parameters
            this.s02 = [];          % estimated variance
            this.v_hat = [];        % residuals of the observarions;
            
            % up to now GPS only is available for goBLock
            frequencies = find(this.state.cc.getGPS().flag_f);
            
            % get wait bar instance
            w_bar = Go_Wait_Bar.getInstance();
            
            % variable name change for readability reasons
            n_sat = this.state.cc.getNumSat();
            
            % init epoch counter
            epoch_track = 0;
            epoch_ok = 1 : length(time_diff);
            this.pivot_track = uint8(zeros(length(time_diff), 1));
            
            % goGPS waiting bar
            w_bar.setBarLen(length(time_diff));
            w_bar.createNewBar('Building the design matrix...');
            
            p_rate = this.state.getProcessingRate();
            % init loop
            for t = epoch_ok
                eph_t = rt_find_eph (eph, time_diff(t), n_sat);
                
                [y0_epo, A_epo, b_epo, Q_epo, this.sat_pr_track(:, t), this.sat_ph_track(:, t), pivot] = this.oneEpochLS (time_diff(t), pos_r, pos_m(:,t), pr1_r(:,t), pr1_m(:,t), pr2_r(:,t), pr2_m(:,t), ph1_r(:,t), ph1_m(:,t), ph2_r(:,t), ph2_m(:,t), snr_r(:,t), snr_m(:,t), eph_t, sp3, iono, lambda, frequencies(1), p_rate, ant_pcv);
                
                if (pivot > 0)
                    n_obs = length(y0_epo);
                    
                    idx = epoch_track + (1 : n_obs)';
                    
                    this.y0( idx) = y0_epo;
                    this.b ( idx) =  b_epo;
                    this.A ( idx, :) = A_epo;
                    this.Q ( idx, idx) = Q_epo;
                    
                    this.obs_track(idx, 3) = 1;
                    this.obs_track(idx, 2) = setdiff(find(this.sat_ph_track(:,t)), pivot);
                    this.obs_track(idx, 1) = t;
                    
                    this.pivot_track(t) = pivot;
                    
                    epoch_track = epoch_track + n_obs;
                else
                    this.empty_epoch = [this.empty_epoch; t];
                end
                w_bar.goTime(t);
            end
            
            % cut the empty epochs
            this.sat_pr_track(:,this.empty_epoch) = [];
            this.sat_ph_track(:,this.empty_epoch) = [];
            this.pivot_track(this.empty_epoch) = [];
            
            this.n_epoch = length(this.pivot_track);
            
            % cut the unused lines
            id_ko = (epoch_track + 1 : size(this.A,1))';
            this.n_obs_tot = epoch_track;
            this.y0(id_ko) = [];
            this.b(id_ko) = [];
            this.A(id_ko, :) = [];
            this.Q(:, id_ko) = [];
            this.Q(id_ko, :) = [];
            this.obs_track(id_ko,:) = [];
            w_bar.close();
            
            % Add to the Design matric the columns relative to the ambiguities
            this.addAmbiguities (lambda);
        end
        
        function [pos, pos_cov] = solveFloatOld (this, full_slip_split)
            % Compute a first float solution
            %
            % METHODS CALL REQUIREMENTS:
            %  -> prepare
            %
            % SYNTAX:
            %   [pos, pos_cov] = this.solveFloat(this)
            %
            % INTERNAL INPUT:
            %   A, y0, b, Q, obs_track, amb_num, amb_prn_track, state, log
            %
            % OUTPUT:
            %   pos     coordinates of the estimated positions
            %   pos_cov         covariance of the estimated positions
            %
            % INTERNAL OUTOUT (properties whos values are changed):
            %   x_float, Cxx, s02, v_hat, pos, pos_cov, is_fixed
            %   y0,  b, A, Q
            %   obs_track, amb_prn_track, n_epoch
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.solveFloat()
            
            this.log.addMarkedMessage('Compute a float solution');
            
            if nargin == 1
                full_slip_split = this.state.getFullSlipSplit();
            end
            
            % Let's first clean short observations arcs
            [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track] = this.remShortArc(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, this.state.getMinArc());
            
            % A block is an interval of continuous phase observations
            [this.col_ok, this.ref_arc, blk_cols, blk_rows] = this.getBlockProperties();
            
            % To disable full slip split put full_slip_split = false;
            if ~full_slip_split
                blk_cols = true(size(blk_cols, 1), 1);
                blk_rows = true(size(blk_rows, 1), 1);
            end
            n_block = size(blk_cols, 2);
            
            this.log.addMarkedMessage(sprintf('Independent blocks found: %d', n_block));
            pivot_change = find(abs(diff(int8(this.pivot_track)))>0);
            row_id  = 1;
            this.ref_arc = zeros(n_block,1);
            bad_blocks = [];
            for i = 1 : n_block
                this.log.addMessage(sprintf('      Processing block %d/%d -------------------------------', i, n_block));
                
                % Extract a subset of the LS system
                y0 = this.y0(blk_rows(:, i));
                b = this.b(blk_rows(:, i));
                A = this.A(blk_rows(:, i), blk_cols(:, i));
                Q = this.Q(blk_rows(:, i), blk_rows(:, i));
                obs_track = this.obs_track(blk_rows(:, i),:);
                epoch_offset = obs_track(1,1) - 1;
                pivot_change = pivot_change - epoch_offset;
                obs_track(:,1) = obs_track(:,1) - epoch_offset;
                amb_prn_track = 4 : size(A,2);
                
                [A, y0, b, Q, obs_track, amb_prn_track] = this.remShortArc(A, y0, b, Q, obs_track, amb_prn_track, this.state.getMinArc());
                
                if size(A,1) > size(A,2)
                    
                    % Get the arc with higher quality
                    [col_ok, ref_arc] = this.getBestRefArc(y0, b, A, Q);
                    
                    if this.state.isPreCleaningOn()
                        this.log.addMessage('       - try to improve observations (risky...check the results!!!)');
                        % Try to correct cycle slips / discontinuities in the observations and increase spike variance
                        % WARNING: risky operation, do it with consciousness, check the results against disabled pre-cleaning
                        %          this feature can be used when the phase residuals show unresolved anbiguities
                        [y0, Q] = this.preCorrectObsIntAmb(y0, b, A, col_ok, Q, this.n_pos, obs_track, pivot_change); % Try to correct integer ambiguities slips (maybe missed cycle slips)
                    end
                    
                    this.log.addMessage('       - first estimation');
                    
                    % computing a first solution with float ambiguities
                    [x_float, Cxx, s02, v_hat] = this.solveLS(y0, b, A, col_ok, Q); %#ok<*ASGLU>
                    
                    this.log.addMessage('       - improve solution by outlier underweight');
                    % Improve solution by iterative increase of bad observations variance
                    [x_float, Cxx, s02, v_hat, Q_tmp] = this.improveFloatSolution(y0, b, A, col_ok, Q, v_hat, obs_track);
                    
                    % Compute phase residuals
                    [phase_res, id_track] = this.computeDataTrack( v_hat, [], A, obs_track, 1);
                    n_clean = this.state.getBlockPostCleaningLoops();
                    
                    % Try to fix missing cycle slips
                    [x_float, Cxx, s02, v_hat, y0, Q_tmp] = this.loopCorrector(y0, b, A, col_ok, Q, obs_track, amb_prn_track, n_clean);
                    
                    if this.state.isBlockForceStabilizationOn()
                        % If the system is unstable remove the arcs that are making it so
                        [A, col_ok, ref_arc, y0, b, Q, obs_track, amb_prn_track] = this.remUnstableArcs(A, col_ok, ref_arc, y0, b, Q, obs_track, amb_prn_track, Cxx);
                    end
                    
                    if (size(A,2) > 3 + 2)
                        if this.state.isOutlierRejectionOn()
                            % Delete bad observations and restore variances
                            this.log.addMessage('       - reject outliers');
                            [x_float, Cxx, s02, v_hat, y0, b, A, col_ok, Q, obs_track, amb_prn_track] = this.cleanObsHiRes(y0, b, A, col_ok, Q, v_hat, obs_track, amb_prn_track, 9, false);
                            %[~, ~, ~, ~, y0,  b, A, ~, Q, obs_track, amb_prn_track] = this.remSolitaryObs(y0, b, A, col_ok, Q, obs_track, amb_prn_track, round(this.state.getMinArc()/2));
                            [col_ok, ref_arc] = this.getBestRefArc(y0, b, A, Q);
                            [x_float, Cxx, s02, v_hat] = this.solveLS(y0, b, A, col_ok, Q);
                            [x_float, Cxx, s02, v_hat, Q_tmp] = this.improveFloatSolution(y0, b, A, col_ok, Q, v_hat, obs_track);
                            if this.state.isBlockForceStabilizationOn()
                                % If the system is unstable remove the arcs that are making it so
                                [A, ~, ref_arc, y0, b, Q, obs_track, amb_prn_track] = this.remUnstableArcs(A, col_ok, ref_arc, y0, b, Q, obs_track, amb_prn_track, Cxx);
                            end
                        else
                            Q = Q_tmp;
                        end
                    end
                end
                
                if (size(A,2) < 3 + 2) || (size(A,1) <= size(A,2))
                    % If the system is completely unstable
                    this.log.addMessage('         [ WW ] system still unstable, try to use it as it is!!!\n                (it might be stabilized by the unified solution)');
                    bad_blocks = [bad_blocks; i]; %#ok<AGROW>
                    y0 = this.y0(blk_rows(:, i));
                    b = this.b(blk_rows(:, i));
                    A = this.A(blk_rows(:, i), blk_cols(:, i));
                    Q = this.Q(blk_rows(:, i), blk_rows(:, i));
                    obs_track = this.obs_track(blk_rows(:, i),:);
                    epoch_offset = obs_track(1,1) - 1;
                    pivot_change = pivot_change - epoch_offset;
                    obs_track(:,1) = obs_track(:,1) - epoch_offset;
                    amb_prn_track = 4 : size(A,2);
                    [~, ref_arc] = this.getBestRefArc(y0, b, A, Q);
                end
                % reassemble the system
                % id on the obj matrix
                row_id_last = row_id + numel(y0) - 1;
                full_row_id = row_id : row_id_last;
                
                full_col_id = find(blk_cols(:, i));
                full_col_id = full_col_id([1 2 3 amb_prn_track]);
                
                
                % find the columns that are still used
                this.ref_arc(i) = full_col_id(ref_arc + 3);
                this.y0(full_row_id) = y0;
                this.b(full_row_id) = b;
                this.A(full_row_id, :) = 0;
                this.A(full_row_id, full_col_id) = A;
                this.Q(full_row_id, :) = 0;
                this.Q(:, full_row_id) = 0;
                this.Q(full_row_id,full_row_id) = Q;
                obs_track(:,1) = obs_track(:,1) + epoch_offset;
                this.obs_track(full_row_id, :) = obs_track;
                row_id = row_id_last + 1;
            end
            
            % remove unecessary rows (observations removed as outliers)
            this.y0(row_id_last + 1 : end) = [];
            this.b(row_id_last + 1 : end) = [];
            this.A(row_id_last + 1 : end, :) = [];
            this.Q(row_id_last + 1 : end,:) = [];
            this.Q(:, row_id_last + 1 : end) = [];
            this.obs_track(row_id_last + 1 : end, :) = [];
            
            [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track] = this.remShortArc(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, this.state.getMinArc());
            this.col_ok = setdiff(1:size(this.A, 2), this.ref_arc + 3);
            
            if full_slip_split
                [~, ~, blk_cols, ~] = this.getBlockProperties();
                if ~isempty(bad_blocks)
                    [this.col_ok, this.ref_arc] = this.getBestBlockRefArc(this.y0, this.b, this.A, this.Q, this.ref_arc, [], bad_blocks, blk_cols);
                end
            end
            
            this.log.addMarkedMessage('Compute the final float solution -------------------');
            [this.x_float, this.Cxx, this.s02, this.v_hat] = fast_least_squares_solver(this.y0, this.b, this.A(:, this.col_ok), this.Q);
            
            this.log.addMessage('       - improve solution by outlier underweight');
            [this.x_float, this.Cxx, this.s02, this.v_hat, Q_tmp] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, this.Q, this.v_hat, this.obs_track);
            [this.phase_res, this.id_track] = this.computeDataTrack();
            
            if full_slip_split
                % Refining final solution if it have been computed in blocks
                
                [this.x_float, this.Cxx, this.s02, this.v_hat, this.y0, Q_tmp] = this.loopCorrector(this.y0, this.b, this.A, this.col_ok, this.Q, this.obs_track, this.amb_prn_track, n_clean);
                [this.phase_res, this.id_track] = this.computeDataTrack();
                
                % If the system is unstable try to change the reference arc
                bad_col = find(abs(median(this.phase_res(:,:,1),'omitnan')) > 1) + 3;
                if ~isempty(bad_col)
                    [this.col_ok, this.ref_arc, bad_blocks] = this.getBestBlockRefArc(this.y0, this.b, this.A, this.Q, this.ref_arc, bad_col - 3, [], blk_cols);
                    if ~isempty(bad_blocks)
                        [this.x_float, this.Cxx, this.s02, this.v_hat, ~] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, this.Q, this.v_hat, this.obs_track);
                        [this.phase_res, this.id_track] = this.computeDataTrack();
                        bad_col = find(abs(median(this.phase_res(:,:,1),'omitnan')) > 1) + 3;
                        [this.x_float, this.Cxx, this.s02, this.v_hat, this.y0, Q_tmp] = this.loopCorrector(this.y0, this.b, this.A, this.col_ok, this.Q, this.obs_track, this.amb_prn_track, n_clean);
                        [this.phase_res, this.id_track] = this.computeDataTrack();
                    end
                end
                
                if this.state.isBlockForceStabilizationOn()
                    % If the system is still unstable remove the with median high residuals
                    while ~isempty(bad_col)
                        
                        while ~isempty(bad_col)
                            [~, ~, blk_cols, ~] = this.getBlockProperties();
                            blk_cols(1 : 3, :) = 0;
                            unstable_block = find(sum([blk_cols; -blk_cols(bad_col,:)]) < 3);
                            if ~isempty(unstable_block)
                                this.log.addWarning(sprintf('One or more block have been found unstable, removing block %s', sprintf('%d ', unstable_block)));
                                bad_col = union(bad_col, find(blk_cols(:, unstable_block)));
                                this.ref_arc(unstable_block) = [];
                            end

                            for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(bad_col < this.ref_arc(a) + 3); end
                            this.col_ok = setdiff(1 : size(this.A, 2), this.ref_arc + 3);

                            this.log.addMessage(sprintf('         [ WW ] System probably unstable, removing arcs: %s PRNs: %s', sprintf('%d ', bad_col - 3), sprintf('%d ', this.amb_prn_track(bad_col - 3))));
                            [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track] = this.remArcCol(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, bad_col);
                            [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, bad_col] = this.remShortArc(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, this.state.getMinArc());
                        end
                        this.ref_arc = setdiff(this.ref_arc, bad_col - 3);
                        for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(bad_col < this.ref_arc(a) + 3); end
                        this.col_ok = setdiff(1 : size(this.A, 2), this.ref_arc + 3);
                        [this.x_float, this.Cxx, this.s02, this.v_hat] = fast_least_squares_solver(this.y0, this.b, this.A(:, this.col_ok), this.Q);
                        [this.x_float, this.Cxx, this.s02, this.v_hat, Q_tmp] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, this.Q, this.v_hat, this.obs_track);
                        [this.phase_res, this.id_track] = this.computeDataTrack();
                        bad_col = find(abs(median(this.phase_res(:,:,1),'omitnan')) > 1) + 3;
                    end
                end
                
                % if flag_outlier
                %     this.log.addMessage('       - reject outliers');
                %     [this.phase_res, this.id_track] = this.computeDataTrack();
                %     N_inv = [];
                %     subset_out = serialize(this.id_track(abs(this.phase_res(:,:,1)) > 0.05));
                %     [x_k, s2_k, v_hat_k, Cxx_k, N_inv] = ELOBO(this.A(:, this.col_ok), this.Q, this.y0, this.b, N_inv, this.v_hat, this.x_float, this.s02, subset_out);
                %     subset_in = serialize(this.id_track(abs(this.phase_res(:,:,1)) < 10 * this.state.getMaxPhaseErrThr));
                %     this.y0 = this.y0(subset_in);
                %     this.b = this.b(subset_in);
                %     this.A = this.A(subset_in, :);
                %     Q_tmp = Q_tmp(subset_in, subset_in);
                %     this.Q = this.Q(subset_in, subset_in);
                %     this.obs_track = this.obs_track(subset_in, :);
                %     [this.x_float, this.Cxx, this.s02, this.v_hat, Q_tmp] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, Q_tmp, [], this.obs_track);
                %     [this.phase_res, this.id_track] = this.computeDataTrack();
                % end
                
                if this.state.isBlockForceStabilizationOn()
                    % If the system is unstable remove the arcs that are making it so
                    amb_var = zeros(size(this.Cxx,1) + numel(this.ref_arc), 1); amb_var(this.col_ok) = diag(this.Cxx);
                    amb_var(amb_var < 0) = 100; % negative variances means bad arcs
                    bad_col = find(amb_var(4:end) > 1) + 3;
                    
                    while ~isempty(bad_col)
                        [~, ~, blk_cols, ~] = this.getBlockProperties();
                        blk_cols(1 : 3, :) = 0;
                        unstable_block = find(sum([blk_cols; -blk_cols(bad_col,:)]) < 4);
                        if ~isempty(unstable_block)
                            this.log.addWarning(sprintf('One or more block have been found unstable, removing block %s', sprintf('%d ', unstable_block)));
                            bad_col = union(bad_col, find(blk_cols(:, unstable_block)));
                            this.ref_arc(unstable_block) = [];
                        end
                        
                        for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(bad_col < this.ref_arc(a) + 3); end
                        this.col_ok = setdiff(1 : size(this.A, 2), this.ref_arc + 3);
                        
                        this.log.addMessage(sprintf('         [ WW ] System unstable, removing arcs: %s\n                                          PRNs: %s', sprintf('%d ', bad_col - 3), sprintf('%d ', this.amb_prn_track(bad_col - 3))));
                        [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track] = this.remArcCol(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, bad_col);
                        [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, bad_col] = this.remShortArc(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, this.state.getMinArc());
                        this.ref_arc = setdiff(this.ref_arc, bad_col - 3);
                        for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(bad_col < this.ref_arc(a) + 3); end
                        this.col_ok = setdiff(1 : size(this.A, 2), this.ref_arc + 3);
                        [this.x_float, this.Cxx, this.s02, this.v_hat] = fast_least_squares_solver(this.y0, this.b, this.A(:, this.col_ok), this.Q);
                        [this.x_float, this.Cxx, this.s02, this.v_hat, Q_tmp] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, this.Q, this.v_hat, this.obs_track);
                        
                        [this.phase_res, this.id_track] = this.computeDataTrack();
                        [this.x_float, this.Cxx, this.s02, this.v_hat, this.y0, Q_tmp] = this.loopCorrector(this.y0, this.b, this.A, this.col_ok, this.Q, this.obs_track, this.amb_prn_track, n_clean);
                        [this.phase_res, this.id_track] = this.computeDataTrack();
                        
                        amb_var = zeros(size(this.Cxx,1) + numel(this.ref_arc), 1); amb_var(this.col_ok) = diag(this.Cxx);
                        amb_var(amb_var < 0) = 100; % negative variances means bad arcs
                        bad_col = find(amb_var(4:end) > 1) + 3;
                        % amb_var = zeros(size(this.Cxx,1) + numel(ref_arc)); amb_var(this.col_ok) = diag(this.Cxx);
                        % bad_col = find((amb_var(4:end) > 10) | (amb_var(4:end) > mean(amb_var(col_ok(4 : end))) + 10 * std(amb_var(col_ok(4 : end))))) + 3;
                    end
                end
            end
            
            this.Q = Q_tmp;
            
            % Compute phase residuals
            [this.phase_res, this.id_track] = this.computeDataTrack();
            
            % show residuals
            %close all; this.plotPhRes();
            
            % extract estimated position
            this.log.newLine();
            this.log.addMarkedMessage('Float solution computed, rover positions corrections:');
            d_pos = reshape(this.x_float(1:this.n_pos * 3), 3, this.n_pos);
            pos = repmat(this.pos0(:), 1, this.n_pos) + d_pos;
            this.log.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', [this.getENU(pos)' this.getDeltaENU(pos)']'));
            pos_cov = full(this.Cxx(1:this.n_pos * 3, 1:this.n_pos * 3));
            this.is_fixed = 0;
            this.pos = pos;
            this.pos_cov = pos_cov;
        end
        
        function [pos, pos_cov] = solveFloat (this, full_slip_split)
            % Compute a first float solution
            %
            % METHODS CALL REQUIREMENTS:
            %  -> prepare
            %
            % SYNTAX:
            %   [pos, pos_cov] = this.solveFloat(this)
            %
            % INTERNAL INPUT:
            %   A, y0, b, Q, obs_track, amb_num, amb_prn_track, state, log
            %
            % OUTPUT:
            %   pos     coordinates of the estimated positions
            %   pos_cov         covariance of the estimated positions
            %
            % INTERNAL OUTOUT (properties whos values are changed):
            %   x_float, Cxx, s02, v_hat, pos, pos_cov, is_fixed
            %   y0,  b, A, Q
            %   obs_track, amb_prn_track, n_epoch
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.solveFloat()
            
            this.log.addMarkedMessage('Compute a float solution');
            
            if nargin == 1
                full_slip_split = this.state.getFullSlipSplit();
            end
            
            % Let's first clean short observations arcs
            [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track] = this.remShortArc(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, this.state.getMinArc());
            
            % A block is an interval of continuous phase observations
            [this.col_ok, this.ref_arc, blk_cols, blk_rows] = this.getBlockProperties();
            
            % To disable full slip split put full_slip_split = false;
            if ~full_slip_split
                blk_cols = true(size(blk_cols, 1), 1);
                blk_rows = true(size(blk_rows, 1), 1);
            end
            n_block = size(blk_cols, 2);
            
            this.log.addMarkedMessage(sprintf('Independent blocks found: %d', n_block));
            pivot_change = find(abs(diff(int8(this.pivot_track)))>0);
            row_id  = 1;
            this.ref_arc = zeros(n_block,1);
            bad_blocks = [];
            for i = 1 : n_block
                this.log.addMessage(sprintf('      Processing block %d/%d -------------------------------', i, n_block));
                
                % Extract a subset of the LS system
                y0 = this.y0(blk_rows(:, i));
                b = this.b(blk_rows(:, i));
                A = this.A(blk_rows(:, i), blk_cols(:, i));
                Q = this.Q(blk_rows(:, i), blk_rows(:, i));
                obs_track = this.obs_track(blk_rows(:, i),:);
                epoch_offset = obs_track(1,1) - 1;
                pivot_change = pivot_change - epoch_offset;
                obs_track(:,1) = obs_track(:,1) - epoch_offset;
                amb_prn_track = 4 : size(A,2);
                

                [A, y0, b, Q, obs_track, amb_prn_track] = this.remShortArc(A, y0, b, Q, obs_track, amb_prn_track, this.state.getMinArc());
                
                if size(A,1) > size(A,2)
                    
                    % Get the arc with higher quality
                    iQ = cholinv(Q);
                    [col_ok, ref_arc] = this.getBestRefArc(y0, b, A, iQ);
                    
                    if this.state.isPreCleaningOn()
                        this.log.addMessage('       - try to improve observations (risky...check the results!!!)');
                        % Try to correct cycle slips / discontinuities in the observations and increase spike variance
                        % WARNING: risky operation, do it with consciousness, check the results against disabled pre-cleaning
                        %          this feature can be used when the phase residuals show unresolved anbiguities
                        [y0, Q] = this.preCorrectObsIntAmb(y0, b, A, Q, this.n_pos, obs_track, pivot_change); % Try to correct integer ambiguities slips (maybe missed cycle slips)
                        iQ = cholinv(Q);
                    end

                    this.log.addMessage('       - first estimation');

                    % computing a first solution with float ambiguities
                    %[x_float, Cxx, s02, v_hat] = this.solveLS(y0, b, A, col_ok, iQ); %#ok<*ASGLU>

                    this.log.addMessage('       - improve solution by outlier underweight');
                    %[x_float, Cxx, s02, v_hat, Q_tmp] = this.improveFloatSolution(y0, b, A, col_ok, Q, [], obs_track);
                    % Improve solution by iterative increase of bad observations variance
                    n_clean = this.state.getBlockPostCleaningLoops();
                    % Try to fix missing cycle slips
                    [x_float, Cxx, s02, v_hat, y0, Q_tmp, iQ_tmp, phase_res, id_track] = this.loopCorrector(y0, b, A, col_ok, Q, obs_track, amb_prn_track, n_clean);

                    if this.state.isBlockForceStabilizationOn()
                        % If the system is unstable remove the arcs that are making it so
                        [A, col_ok, ref_arc, y0, b, Q, obs_track, amb_prn_track] = this.remUnstableArcs(A, col_ok, ref_arc, y0, b, Q, obs_track, amb_prn_track, Cxx);
                        %[x_float, Cxx, s02, v_hat, Q_tmp] = this.improveFloatSolution(y0, b, A, col_ok, Q, [], obs_track);
                        [x_float, Cxx, s02, v_hat, y0, Q_tmp, iQ_tmp, phase_res, id_track] = this.loopCorrector(y0, b, A, col_ok, Q, obs_track, amb_prn_track, n_clean);
                    end

                    if (size(A,2) > 3 + 2)
                        if this.state.isOutlierRejectionOn()
                            n_obs_tmp = size(A,1);
                            % Delete bad observations and restore variances
                            this.log.addMessage('       - reject outliers');
                            [x_float, Cxx, s02, v_hat, y0, b, A, col_ok, Q, obs_track, amb_prn_track, iQ] = this.cleanObsHiRes(y0, b, A, col_ok, Q, v_hat, obs_track, amb_prn_track, 9, false);
                            %[~, ~, ~, ~, y0,  b, A, ~, Q, obs_track, amb_prn_track, iQ] = this.remSolitaryObs(y0, b, A, col_ok, Q, obs_track, amb_prn_track, round(this.state.getMinArc()/2));
                            [col_ok, ref_arc] = this.getBestRefArc(y0, b, A, iQ);
                            [x_float, Cxx, s02, v_hat, Q_tmp] = this.improveFloatSolution(y0, b, A, col_ok, Q, [], obs_track);
                            if this.state.isBlockForceStabilizationOn()
                                % If the system is unstable remove the arcs that are making it so
                                [A, ~, ref_arc, y0, b, Q, obs_track, amb_prn_track] = this.remUnstableArcs(A, col_ok, ref_arc, y0, b, Q, obs_track, amb_prn_track, Cxx);
                            end
                            
                            if size(A,1) < 0.6 * n_obs_tmp
                                this.log.addMessage('         [ WW ] too many outliers detected, use the block as it is!');
                                A = [];
                            else
                                % PLOT close all;
                                [phase_res, id_track] = this.computeDataTrack(v_hat, [], A, obs_track, 1);
                                % PLOT this.plotDataTrack(phase_res, id_track, A, amb_prn_track);% ylim([-0.5 0.5])

                                % Detect new full CS candidates -------------------------------------------

                                this.log.addMessage('       - detect new cycle slip candidates');
                                sensor = median(movstd(phase_res,3,'omitnan'),2,'omitnan');
                                detector = find(sensor > min(0.1, 9 * std(sensor,'omitnan')) | isnan(sensor));

                                % PLOT figure; plot(sensor); hold on; plot(detector, sensor(detector),'o'); dockAllFigures;

                                lim = unique([0; detector; numel(sensor)+1]);
                                lim = [lim(1 : end - 1) + 1, lim(2 : end) - 1];
                                lim = lim((lim(:,2) - lim(:,1)) > 5,:);

                                this.log.addMessage(sprintf('            %d candidate found', size(lim, 1) - 1));

                                % Strong cleaning based on the check of the residuals derivatives ---------

                                this.log.addMessage('       - stronger reject outliers');
                                [x_float, Cxx, s02, v_hat, y0, b, A, col_ok, Q, obs_track, amb_prn_track, iQ] = this.cleanObsHiResVar(y0, b, A, col_ok, Q, v_hat, obs_track, amb_prn_track, 0.997);
                                %[~, ~, ~, ~, y0,  b, A, ~, Q, obs_track, amb_prn_track, iQ] = this.remSolitaryObs(y0, b, A, col_ok, Q, obs_track, amb_prn_track, round(this.state.getMinArc()/2));
                                [col_ok, ref_arc] = this.getBestRefArc(y0, b, A, iQ);
                                [x_float, Cxx, s02, v_hat, Q_tmp] = this.improveFloatSolution(y0, b, A, col_ok, Q, [], obs_track);

                                [phase_res, id_track] = this.computeDataTrack(v_hat, [], A, obs_track, 1);
                                % PLOT this.plotDataTrack(phase_res, id_track, A, amb_prn_track);% ylim([-0.5 0.5])

                                if full_slip_split && (size(lim, 1) > 1) && (size(lim, 2) > 1) && (s02 > 0.02)

                                    for l = 1 : size(lim, 1)
                                        lim_obs = find(obs_track(:,1) >= lim(l, 1), 1, 'first');
                                        if ~(isempty(lim_obs))
                                            lim(l, 1) = find(obs_track(:,1) >= lim(l, 1), 1, 'first');
                                            lim(l, 2) = find(obs_track(:,1) <= lim(l, 2), 1, 'last');
                                        else
                                            lim(l, : ) = [];
                                        end
                                    end

                                    % Building a splitted system to compute separate and repair CS ------------

                                    this.log.addMessage('          - repair cycle slips');
                                    A_ck = spalloc(size(A, 1), 10000, sum(A(:)~=0));
                                    row = 0;
                                    s02_vec = [];
                                    v_hat3 = [];
                                    obs_track2 = [];
                                    obs_in = [];
                                    col = 3;
                                    col_ok_ck = [1 2 3];
                                    amb_prn_track_ck = [];
                                    for j = 1 : size(lim, 1)
                                        obs_id = lim(j,1) : lim(j,2);
                                        row = row(end) + (1 : numel(obs_id));
                                        arc_id = find(full(sum(abs(A(obs_id, 4 : end)))) > 0);
                                        col_id = arc_id + 3;
                                        if ~isempty(col_id)
                                            col = col(end) + (1 : numel(col_id));
                                            A_ck(row, col) = A(obs_id, col_id);
                                            A_ck(row, 1:3) = A(obs_id,1:3);
                                            amb_prn_track_ck(col-3) = amb_prn_track(col_id - 3);
                                        end
                                        obs_in = [obs_in obs_id];
                                        [col_ok_tmp] = this.getBestRefArc(y0(obs_id), b(obs_id),  A_ck(row, [1 2 3 col]), cholinv(Q(obs_id,obs_id)));
                                        col_ok_ck = [col_ok_ck (col(end) - numel(col_id) + col_ok_tmp(4:end) -3)];
                                    end
                                    A_ck = A_ck(1 : row(end), 1 : col(end));

                                    % remove the unused columns / rows
                                    A = A(obs_in, :);
                                    y0 = y0(obs_in);
                                    b = b(obs_in);
                                    Q = Q(obs_in, obs_in);
                                    iQ = cholinv(Q);
                                    obs_track = obs_track(obs_in, :);

                                    %[x, Cxx, s02, v_hat_ck] = Core_Block.solveLS(y0, b, A_ck, col_ok_ck, iQ);
                                    [x, Cxx, s02, v_hat_ck, Q_tmp] = this.improveFloatSolution(y0, b, A_ck, col_ok_ck, Q, [], obs_track);
                                    % [x, Cxx, s02, v_hat_ck, y0, Q_tmp, iQ_tmp, phase_res, id_track] = this.loopCorrector(y0, b, A_ck, col_ok_ck, Q, obs_track, amb_prn_track_ck, 4);
                                    [phase_res, id_track] = this.computeDataTrack(v_hat_ck, [], A, obs_track, 1);
                                    % PLOT this.plotDataTrack(phase_res, id_track, A, amb_prn_track);% ylim([-0.5 0.5])

                                    % Repair CS ---------------------------------------------------------------

                                    amb_correction = zeros(size(A_ck,2) - 3, 1);
                                    amb_correction(col_ok_ck(4 : end) - 3, 1) = round(x(4 : end) / this.cs_factor) * this.cs_factor;

                                    for a = 1 : size(amb_correction, 1)
                                        y0 = y0 - A_ck(:, a + 3) .* amb_correction(a, 1);
                                    end
                                    %[x, Cxx, s02, v_hat_ck, Q_tmp] = this.improveFloatSolution(y0_ck, b, A, col_ok, Q, [], obs_track);
                                    [col_ok, ref_arc] = this.getBestRefArc(y0, b, A, iQ);
                                    [x_float, Cxx, s02, v_hat, y0, Q_tmp, iQ_tmp, phase_res, id_track] = this.loopCorrector(y0, b, A, col_ok, Q, obs_track, amb_prn_track, 4);
                                    % PLOT this.plotDataTrack(phase_res, id_track, A, amb_prn_track);% ylim([-0.5 0.5])
                                end
                            end
                        end
                    end
                end
                %if (s02 > 0.02) || (size(A,2) < 3 + 2) || (size(A,1) <= size(A,2))
                if (size(A,2) < 3 + 2) || (size(A,1) <= size(A,2))
                    % If the system is completely unstable
                    if size(A,2) < 3
                        s02 = inf;
                    end
                    this.log.addMessage(sprintf('         [ WW ] system too small, try to use it unprocessed!!!\n                (it might be stabilized by the unified solution)'));

                    bad_blocks = [bad_blocks; i]; %#ok<AGROW>
                    y0 = this.y0(blk_rows(:, i));
                    b = this.b(blk_rows(:, i));
                    A = this.A(blk_rows(:, i), blk_cols(:, i));
                    Q = this.Q(blk_rows(:, i), blk_rows(:, i));
                    obs_track = this.obs_track(blk_rows(:, i),:);
                    epoch_offset = obs_track(1,1) - 1;
                    pivot_change = pivot_change - epoch_offset;
                    obs_track(:,1) = obs_track(:,1) - epoch_offset;
                    amb_prn_track = 4 : size(A,2);
                    [~, ref_arc] = this.getBestRefArc(y0, b, A, Q);
                else
                    if (s02 > 0.02)
                        % If the system is completely unstable
                        this.log.addMessage(sprintf('         [ WW ] system unstable (s02 = %.4f), try to use it as it is!!!\n                (it might be stabilized by the unified solution)', s02));
                    end
                end

                % reassemble the system
                % id on the obj matrix
                row_id_last = row_id + numel(y0) - 1;
                full_row_id = row_id : row_id_last;

                full_col_id = find(blk_cols(:, i));
                full_col_id = full_col_id([1 2 3 amb_prn_track]);


                % find the columns that are still used
                this.ref_arc(i) = full_col_id(ref_arc + 3);
                this.y0(full_row_id) = y0;
                this.b(full_row_id) = b;
                this.A(full_row_id, :) = 0;
                this.A(full_row_id, full_col_id) = A;
                this.Q(full_row_id, :) = 0;
                this.Q(:, full_row_id) = 0;
                this.Q(full_row_id,full_row_id) = Q;
                obs_track(:,1) = obs_track(:,1) + epoch_offset;
                this.obs_track(full_row_id, :) = obs_track;
                row_id = row_id_last + 1;
            end

            % remove unecessary rows (observations removed as outliers)
            this.y0(row_id_last + 1 : end) = [];
            this.b(row_id_last + 1 : end) = [];
            this.A(row_id_last + 1 : end, :) = [];
            this.Q(row_id_last + 1 : end,:) = [];
            this.Q(:, row_id_last + 1 : end) = [];
            this.obs_track(row_id_last + 1 : end, :) = [];

            [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, rem_col] = this.remShortArc(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, this.state.getMinArc());
            for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(rem_col < this.ref_arc(a) + 3); end
            this.ref_arc = unique(this.ref_arc);
            this.col_ok = setdiff(1:size(this.A, 2), this.ref_arc + 3);

            iQ = cholinv(this.Q);
            if full_slip_split
                [~, ~, blk_cols, ~] = this.getBlockProperties();
                if ~isempty(bad_blocks)
                    [this.col_ok, this.ref_arc] = this.getBestBlockRefArc(this.y0, this.b, this.A, iQ, this.ref_arc, [], bad_blocks, blk_cols);
                    this.ref_arc = unique(this.ref_arc);
                end
            else
                [this.col_ok, this.ref_arc] = this.getBestRefArc(this.y0, this.b, this.A, iQ);
            end

            this.log.addMarkedMessage('Compute the final float solution -------------------');
            this.log.addMessage('       - improve solution by outlier underweight');
            %[this.x_float, this.Cxx, this.s02, this.v_hat, Q_tmp] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, this.Q, [], this.obs_track);

            [this.x_float, this.Cxx, this.s02, this.v_hat, this.y0, Q_tmp, iQ_tmp, this.phase_res, this.id_track] = this.loopCorrector(this.y0, this.b, this.A, this.col_ok, this.Q, this.obs_track, this.amb_prn_track, n_clean);

            % If the system is unstable try to change the reference arc
            bad_col = find(abs(median(this.phase_res(:,:,1),'omitnan')) > 1) + 3;
            if ~isempty(bad_col)
                [this.col_ok, this.ref_arc, bad_blocks] = this.getBestBlockRefArc(this.y0, this.b, this.A, iQ, this.ref_arc, bad_col - 3, [], blk_cols);
                if ~isempty(bad_blocks)
                    [this.x_float, this.Cxx, this.s02, this.v_hat, this.y0, Q_tmp, iQ_tmp, this.phase_res, this.id_track] = this.loopCorrector(this.y0, this.b, this.A, this.col_ok, this.Q, this.obs_track, this.amb_prn_track, n_clean);
                    bad_col = find(abs(median(this.phase_res(:,:,1),'omitnan')) > 1) + 3;
                end
            end

            if this.state.isBlockForceStabilizationOn()
                % If the system is still unstable remove the with median high residuals
                while ~isempty(bad_col)

                    while ~isempty(bad_col)
                        [~, ~, blk_cols, ~] = this.getBlockProperties();
                        if ~full_slip_split
                            blk_cols = true(size(blk_cols, 1), 1);
                            blk_rows = true(size(blk_rows, 1), 1);
                        end
                        blk_cols(1 : 3, :) = 0;
                        unstable_block = find(sum([blk_cols; -blk_cols(bad_col,:)]) < 3);
                        if ~isempty(unstable_block)
                            this.log.addWarning(sprintf('One or more block have been found unstable, removing block %s', sprintf('%d ', unstable_block)));
                            bad_col = union(bad_col, find(blk_cols(:, unstable_block)));
                            this.ref_arc(unstable_block) = [];
                        end

                        for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(bad_col < this.ref_arc(a) + 3); end
                        this.col_ok = setdiff(1 : size(this.A, 2), this.ref_arc + 3);

                        this.log.addMessage(sprintf('         [ WW ] System unstable, removing arcs: %s PRNs: %s', sprintf('%d ', bad_col - 3), sprintf('%d ', this.amb_prn_track(bad_col - 3))));
                        [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track] = this.remArcCol(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, bad_col);
                        [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, bad_col] = this.remShortArc(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, this.state.getMinArc());
                    end
                    this.ref_arc = setdiff(this.ref_arc, bad_col - 3);
                    for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(bad_col < this.ref_arc(a) + 3); end
                    this.col_ok = setdiff(1 : size(this.A, 2), this.ref_arc + 3);
                    [this.x_float, this.Cxx, this.s02, this.v_hat, Q_tmp] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, this.Q, [], this.obs_track);
                    [this.phase_res, this.id_track] = this.computeDataTrack(this.v_hat);
                    bad_col = find(abs(median(this.phase_res(:,:,1),'omitnan')) > 1) + 3;
                end
            end

            if this.state.isBlockForceStabilizationOn()
                % If the system is unstable remove the arcs that are making it so
                amb_var = zeros(size(this.Cxx,1) + numel(this.ref_arc), 1); amb_var(this.col_ok) = diag(this.Cxx);
                amb_var(amb_var < 0) = 100; % negative variances means bad arcs
                bad_col = find(amb_var(4:end) > 1) + 3;

                while ~isempty(bad_col)
                    [~, ~, blk_cols, ~] = this.getBlockProperties();
                    if ~full_slip_split
                        blk_cols = true(size(blk_cols, 1), 1);
                        blk_rows = true(size(blk_rows, 1), 1);
                    end
                    blk_cols(1 : 3, :) = 0;
                    unstable_block = find(sum([blk_cols; -blk_cols(bad_col,:)]) < 4);
                    if ~isempty(unstable_block)
                        this.log.addWarning(sprintf('One or more block have been found unstable, removing block %s', sprintf('%d ', unstable_block)));
                        bad_col = union(bad_col, find(blk_cols(:, unstable_block)));
                        this.ref_arc(unstable_block) = [];
                    end

                    for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(bad_col < this.ref_arc(a) + 3); end
                    this.col_ok = setdiff(1 : size(this.A, 2), this.ref_arc + 3);

                    this.log.addMessage(sprintf('         [ WW ] System unstable, removing arcs: %s\n                                          PRNs: %s', sprintf('%d ', bad_col - 3), sprintf('%d ', this.amb_prn_track(bad_col - 3))));
                    [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track] = this.remArcCol(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, bad_col);
                    [this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, bad_col] = this.remShortArc(this.A, this.y0, this.b, this.Q, this.obs_track, this.amb_prn_track, this.state.getMinArc());
                    this.ref_arc = setdiff(this.ref_arc, bad_col - 3);
                    for a = 1 : numel(this.ref_arc); this.ref_arc(a) = this.ref_arc(a) - sum(bad_col < this.ref_arc(a) + 3); end
                    this.col_ok = setdiff(1 : size(this.A, 2), this.ref_arc + 3);

                    [this.x_float, this.Cxx, this.s02, this.v_hat, this.y0, Q_tmp, iQ_tmp, this.phase_res, this.id_track] = this.loopCorrector(this.y0, this.b, this.A, this.col_ok, this.Q, this.obs_track, this.amb_prn_track, n_clean);

                    amb_var = zeros(size(this.Cxx,1) + numel(this.ref_arc), 1); amb_var(this.col_ok) = diag(this.Cxx);
                    amb_var(amb_var < 0) = 100; % negative variances means bad arcs
                    bad_col = find(amb_var(4:end) > 1) + 3;
                end
            end
            if ~this.state.isBlockOneArc()
                this.Q = Q_tmp;
            end

            % show residuals
            % close all; this.plotDataTrack();

            % if I'm considering non integer cycle slips (half cycle jumps) -> correct them in the data (lambda do not use half cycles)
            % -> it should be better to change the lambda fixer to deal with half cycles
            if (this.cs_factor < 1)
                amb_correction = zeros(size(this.A,2) - 3, 1);
                amb_correction(this.col_ok(4 : end) - 3, 1) = round(this.x_float(4 : end) / this.cs_factor) * this.cs_factor;
                y0 = [this.y0];
                for a = 1 : size(amb_correction, 1)
                    y0(:) = y0(:) - this.A(:, 3 + a) .* amb_correction(a, 1);
                end
                this.y0 = y0;
                this.x_float(4:end) = this.x_float(4:end) - round(this.x_float(4 : end) / this.cs_factor) * this.cs_factor;
            end

            % extract estimated position
            this.log.newLine();
            this.log.addMarkedMessage('Float solution computed, rover positions corrections:');
            d_pos = reshape(this.x_float(1:this.n_pos * 3), 3, this.n_pos);
            pos = repmat(this.pos0(:), 1, this.n_pos) + d_pos;
            this.log.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', [this.getENU(pos)' this.getDeltaENU(pos)']'));
            pos_cov = full(this.Cxx(1:this.n_pos * 3, 1:this.n_pos * 3));
            this.is_fixed = 0;
            this.pos = pos;
            this.pos_cov = pos_cov;
        end

        function [pos, pos_cov, amb_fix, amb_cov, amb_fix_full, ref_arc, G] = solveFix (this)
            % Compute a fixed solution using LAMBDA, and the the internal object properties
            %
            % METHODS CALL REQUIREMENTS:
            %   prepare -> addAmbiguities -> solveFloat
            %
            % SYNTAX:
            %   [pos, pos_cov, amb_fix, amb_cov, amb_fix_full, ref_arc] = this.solveFix()
            %
            % INTERNAL INPUT:
            %   x_float, Cxx, A, n_pos, state, log
            %
            % OUTPUT:
            %   pos             coordinates of the estimated positions
            %   pos_cov         covariance of the estimated positions
            %   amb_fix         ambiguities as estimated by lambda (n-1 w.r.t. float solution)
            %   amb_cov         ambiguities error covariance matrix
            %   amb_fix_full    ambiguities as converted from fix to float -> to be imported as pseudo observations of the float solution
            %   ref_arc         arc used as reference in the fix solution (it's the arc that create a bias in the solution)
            %   G               transformation matrix -> float -> fix
            %
            % INTERNAL OUTOUT (properties whos values are changed):
            %   x_fix, Cxx, pos, pos_cov, is_fixed
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.addAmbiguities(lambda)
            %   go_block.solveFloat()
            %   go_block.solveFix()
            %
            % CONCRETE IMPLEMENTATION IN:
            %   solveFixPar
            %

            this.log.addMarkedMessage('Compute ambiguity fix through LAMBDA');

            [d_pos, pos_cov, is_fixed, amb_fix, amb_cov, amb_fix_full, ref_arc, G] = this.solveFixPar (this.x_float, this.Cxx, size(this.A, 2) - 3 - numel(this.ref_arc));
            this.is_fixed = is_fixed;

            pos = this.pos;
            if (is_fixed)
                % extract estimated position
                this.log.addMarkedMessage('Fixed solution computed, rover positions corrections:');
                pos = repmat(this.pos0(:), 1, this.n_pos) + repmat(d_pos(:), 1, this.n_pos);
                this.log.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', [this.getENU(pos)' this.getDeltaENU(pos)']'));
                this.pos = pos;
                this.pos_cov = pos_cov;
                this.pos_cov_fix = pos_cov;
                this.x_fix = [d_pos; amb_fix_full];
            end

        end

        function [pos] = solve(this, s_rate, full_slip_split)
            % Solve Float -> Fix -> try an estimation of positions at a different rate
            %
            % METHODS CALL REQUIREMENTS:
            %  -> prepare
            %
            % SYNTAX:
            %   [pos, pos_cov, v_hat] = this.solve(this)
            %
            % INTERNAL INPUT:
            %   full object properties
            %
            % OUTPUT:
            %   pos     coordinates of the estimated positions
            %   pos_cov         covariance of the estimated positions
            %
            % INTERNAL OUTOUT (properties whos values are changed):
            %   x_float, Cxx, s02, v_hat, pos, pos_cov, is_fixed
            %   y0,  b, A, Q
            %   obs_track, amb_prn_track, n_epoch
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   [pos, pos_cov] = go_block.solve();

            if (nargin == 1) || isempty(s_rate)
                s_rate = this.time_diff(end);
            end
            if nargin < 3
                full_slip_split = this.state.getFullSlipSplit();
            end

            % Compute the best float solution
            this.solveFloat(full_slip_split);

            if (this.state.flag_iar)
                % Solve Fix -> get a valid estimation of the integer ambiguities
                [~, ~, this.amb_fix, ~, this.amb_fix_full, ~, this.G] = this.solveFix();

                % Correct the observations for the cycle slips and merge the arcs
                if this.state.isBlockOneArc()
                    %% Block one arc
                    this.log.addMarkedMessage('Compute the float solution correcting cycle slips\n(use one arc per satellite)');
                    amb_correction = zeros(size(this.A, 2) - 3, 2);
                    amb_correction(this.col_ok(4:end) - 3, 1) = this.amb_fix_full;
                    %amb_correction(this.col_ok(4:end)-3, 2) = round(this.amb_fix_full / this.cs_factor) * this.cs_factor;
                    y0 = [this.y0 this.y0];
                    [amb_prn_track, ~, id_reorder] = unique(this.amb_prn_track);
                    n_amb = numel(amb_prn_track);
                    A = [this.A(:, 1 : 3) sparse(size(this.A, 1), n_amb)];

                    for a = 1 : size(amb_correction, 1)
                        y0(:,1) = y0(:,1) - this.A(:, 3 + a) .* amb_correction(a, 1);
                        %y0(:,2) = y0(:,2) - this.A(:, 3 + a) .* amb_correction(a, 2);
                        A(:, id_reorder(a) + 3) = A(:, id_reorder(a) + 3) + this.A(:, a + 3);
                    end

                    this.y0 = y0(:,1);
                    this.A = A;
                    this.amb_prn_track = amb_prn_track;
                    [this.col_ok, this.ref_arc] = this.getBestRefArc(y0(:,1), this.b, A, cholinv(this.Q));

                    [this.x_float, this.Cxx, this.s02, this.v_hat, Q_tmp] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, this.Q, [], this.obs_track);
                    [this.phase_res, this.id_track] = this.computeDataTrack( this.v_hat, [], this.A, this.obs_track, 1);
                    [this.x_float, this.Cxx, this.s02, this.v_hat, this.y0, Q_tmp, iQ_tmp, this.phase_res, this.id_track] = this.loopCorrector(this.y0, this.b, this.A, this.col_ok, this.Q, this.obs_track, this.amb_prn_track, this.state.getBlockPostCleaningLoops());

                    % Delete bad observations and restore variances
                    if this.state.isOutlierRejectionOn()
                        this.log.addMessage('       - reject outliers');
                        Q_tmp = this.Q;
                        [this.x_float, this.Cxx, this.s02, this.v_hat, this.y0, this.b, this.A, this.col_ok, Q_tmp, this.obs_track, this.amb_prn_track] = this.cleanObsHiRes(this.y0, this.b, this.A, this.col_ok, Q_tmp, this.v_hat, this.obs_track, this.amb_prn_track, 9, false);
                        %[~, ~, ~, ~, this.y0,  this.b, this.A, ~, Q_tmp, this.obs_track, this.amb_prn_track, iQ] = this.remSolitaryObs(this.y0, this.b, this.A, this.col_ok, Q_tmp, this.obs_track, this.amb_prn_track, round(this.state.getMinArc()/2));
                        [this.col_ok, this.ref_arc] = this.getBestRefArc(this.y0, this.b, this.A, cholinv(Q_tmp));
                        [this.x_float, this.Cxx, this.s02, this.v_hat, Q_tmp] = this.improveFloatSolution(this.y0, this.b, this.A, this.col_ok, Q_tmp, [], this.obs_track);
                    end
                    this.Q = Q_tmp;

                    [this.phase_res, this.id_track] = this.computeDataTrack(this.v_hat, [], this.A, this.obs_track, 1);

                    this.log.newLine();
                    this.log.addMarkedMessage('Second Float solution computed, rover positions corrections:');
                    d_pos = reshape(this.x_float(1:this.n_pos * 3), 3, this.n_pos);
                    pos = repmat(this.pos0(:), 1, this.n_pos) + d_pos;
                    this.log.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', [this.getENU(pos)' this.getDeltaENU(pos)']'));
                    pos_cov = full(this.Cxx(1:this.n_pos * 3, 1:this.n_pos * 3));
                    this.is_fixed = 0;
                    this.pos = pos;
                    this.pos_cov = pos_cov;
                    [~, ~, this.amb_fix, ~, this.amb_fix_full, ~, this.G] = this.solveFix();
                end

                %% HR prediction
                if nargin == 2 % if s_rate is defined
                    this.solveHighRate(s_rate);
                end
            end
            pos = this.pos;
        end

        function [pos, pos_cov, v_hat] = solveHighRate(this, s_rate, use_float, reg_factor)
            % Solve Float -> Fix -> try an estimation of positions at a different rate
            %
            % METHODS CALL REQUIREMENTS:
            %  -> prepare
            %
            % SYNTAX:
            %   [pos, pos_cov, v_hat] = this.solve(this)
            %
            % INTERNAL INPUT:
            %   full object properties
            %
            % OUTPUT:
            %   pos     coordinates of the estimated positions
            %   pos_cov         covariance of the estimated positions
            %
            % INTERNAL OUTOUT (properties whos values are changed):
            %   x_float, Cxx, s02, v_hat, pos, pos_cov, is_fixed
            %   y0,  b, A, Q
            %   obs_track, amb_prn_track, n_epoch
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   [pos, pos_cov] = go_block.solve();

            % HR prediction
            narginchk(2,4);
            if nargin == 2
                use_float = false;
            end
            if nargin < 4
                reg_factor = 0;
            end

            s_rate = s_rate(1);

            if s_rate > 0
                this.log.addMarkedMessage(sprintf('Computing high rate solution @%d seconds', s_rate));

                % Find ids of observations involved in a certain time span (solution rate)

                time_track = this.time_diff(this.obs_track(:,1));
                s_time_lim = unique([time_track(1) : s_rate : time_track(end) time_track(end)]');
                s_time_lim = [s_time_lim(1 : end-1) s_time_lim(2 : end)];
                block_id_lim = zeros(size(s_time_lim));
                n_pos_hr = size(s_time_lim, 1);  % new number of positions
                for l = 1 : n_pos_hr
                    block_id_lim(l, 1) = find(time_track >= s_time_lim(l, 1), 1, 'first');
                    block_id_lim(l, 2) = find(time_track < s_time_lim(l, 2), 1, 'last');
                end
                block_id_lim(end,2) = size(this.obs_track,1);

                % Buid the Design Matrix for the high rate estimation
                n_pos_hr = size(block_id_lim, 1);
                A_hr = sparse(size(this.A, 1), size(this.A, 2) + 3 * (n_pos_hr - this.n_pos));
                for i = 1 : n_pos_hr
                    id = block_id_lim(i,1) : block_id_lim(i,2);
                    A_hr(id, (3 * (i - 1) + 1 : 3 * i)) = this.A(id,1:3);
                end
                i = i + 1;
                A_hr(:,(3 * (i - 1) + 1) : end) = this.A(:, 4 : end);

                amb_hr_ok = [(1 : 3 * n_pos_hr) (this.col_ok(3 + 1 : end) + 3 * (n_pos_hr - this.n_pos))];
                % Solve the new system

                % Regularization matrix
                if reg_factor > 0
                    R = 2 * eye(3 * n_pos_hr) - diag(ones(3 * n_pos_hr - 3, 1), 3);
                    R(1:3,1:3) =  R(1:3,1:3) ./ 2;
                    R = sparse(R - diag(ones(3 * n_pos_hr - 3, 1), -3));
                    R = R + 1e-6*speye(3 * n_pos_hr, 3 * n_pos_hr);
                    R_full = sparse(numel(amb_hr_ok), numel(amb_hr_ok));
                    R_full(1 : 3 * n_pos_hr, 1 : 3 * n_pos_hr) = reg_factor .* R;
                    [x_float, Cxx, s02, v_hat] = this.solveRegLS(this.y0, this.b, A_hr, amb_hr_ok, cholinv(this.Q), R_full); %#ok<PROPLC>
                else
                    [x_float, Cxx, s02, v_hat] = this.solveLS(this.y0, this.b, A_hr, amb_hr_ok, cholinv(this.Q)); %#ok<PROPLC>
                end

                % Transformation matrix for the estimation of a unique position (computed as the mean)
                % T = zeros(size(x_float,1) - (n_pos_hr - this.n_pos) * 3, size(x_float,1));
                % for i = 1 : n_pos_hr
                %     T(1 : 3, (i - 1) * 3 + 1 : i * 3) = eye(3) ./ n_pos_hr;
                % end
                % T(4 : end, n_pos_hr * 3 + 1 : end) = eye(numel(x_float) - 3 * n_pos_hr);
                % x_float_mean_pos = T * x_float;
                % Cxx_mean_pos = T * Cxx * T';
                % [~, ~, ~, this.amb_fix, ~, this.amb_fix_full, ~, G] = this.solveFixPar (x_float_mean_pos, Cxx_mean_pos, size(this.A, 2) - 3 - 1);

                if this.is_fixed && ~use_float
                    % Check for positions that cannot be estimated
                    pos_nan = [];
                    if (sum(isnan(x_float(1 : 3 * n_pos_hr))) > 0)
                        this.log.addWarning('Some high rate positions have not been estimated!');
                        pos = x_float(1 : 3 * n_pos_hr);
                        pos_nan = isnan(pos);
                        pos_nan(3 * n_pos_hr + 1 : end) = false;
                        x_float(pos_nan) = [];
                        Cxx(:, pos_nan) = [];
                        Cxx(pos_nan, :) = [];
                    end

                    % Fix ambiguities
                    [d_pos, pos_cov] = this.applyFix(x_float,  Cxx, this.amb_fix, this.G);

                    if ~isempty(pos_nan)
                        pos_nan = pos_nan(1 : 3 * n_pos_hr);
                        d_pos_new = nan(3, n_pos_hr, 1);
                        d_pos_new(~pos_nan) = d_pos;
                        d_pos = d_pos_new; clear d_pos_new;
                        pos_cov_new = nan(3 * n_pos_hr, 3* n_pos_hr);
                        pos_cov_new(~pos_nan, ~pos_nan) = pos_cov;
                        pos_cov = pos_cov_new; clear pos_cov_new;
                    end
                    v_hat = this.y0 - (A_hr(:, amb_hr_ok) * [d_pos(:); this.amb_fix_full] + this.b);
                else
                    this.log.addWarning('The computed high rate solution is NOT fixed (float)!!!')
                    d_pos = reshape(x_float(1 : 3 * n_pos_hr), 3, n_pos_hr);
                    pos_cov = Cxx(1 : 3 * n_pos_hr, 1 : 3 * n_pos_hr);
                end

                this.x_hr = [d_pos(:); this.amb_fix_full];
                pos = repmat(this.pos0, 1, n_pos_hr) + d_pos;

                this.log.addMarkedMessage('High Rate solution computed, rover positions corrections (mean HR):');
                this.log.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', [mean(this.getENU(pos),1, 'omitnan')' mean(this.getDeltaENU(pos),1, 'omitnan')']'));
                this.is_fixed = this.is_fixed * 2;
                this.n_pos_hr = n_pos_hr;
                this.pos = pos;
                this.pos_cov = pos_cov;
                this.s_rate = s_rate;
            end
        end

        function [A, y0, b, Q] = applyFixSuggestion(this, amb_fix_full, amb_cov, ref_arc, A, y0, b, Q, n_pos)
            % Given a set of fixed suggestions build a new LS system adding new pseudo observations
            %
            % SYNTAX:
            %   [A, y0, b, Q] = addFixSuggestion(this, amb_fix_full, amb_cov, ref_arc, A, y0, b, Q, n_pos)

            narginchk(4, 9);
            if nargin < 9
                A = this.A;
                y0 = this.y0;
                b = this.b;
                Q = this.Q;
                n_pos = this.n_pos;
            end
            n_amb = numel(amb_fix_full) - 1;
            n_obs = size(A, 1);
            n_est = n_pos;

            % Build a Design matrix with LAMBDA fix suggestions
            idx = [1 : ref_arc - 1, ref_arc + 1 : n_amb + 1]; % indexes of the fixed ambiguities returned by lambda
            n_amb = numel(idx);
            A = [A; sparse(n_amb, size(A, 2))];               % A matrix gains n rows, 1 for each fixed ambiguity
            y0 = [y0; round(amb_fix_full(idx))];           % y0 array gains n rows with the values of the fixed ambiguities
            b = [b; zeros(n_amb, 1)];                         % b array gains n rows of zeros
            Q = [[Q sparse(size(Q, 1), n_amb)]; sparse(n_amb, size(Q, 1) + n_amb)]; % Q matrix gains n rows and cols, 1 for each fixed ambiguity

            for i = 1 : n_amb
                A(n_obs + i, 3 * n_est + idx(i)) = 1;
            end
            Q(n_obs + (1:n_amb), n_obs +  (1:n_amb)) = amb_cov;
        end

        function [A, y0, b, Q, obs_track, amb_prn_track] = remObs (this, A, y0, b, Q, obs_track, amb_prn_track, rem_obs)
            A(rem_obs,:) = [];
            y0(rem_obs) = [];
            b(rem_obs) = [];
            Q(rem_obs,:) = []; Q(:,rem_obs) = [];
            obs_track(rem_obs,:) = [];
            [A, y0, b, Q, obs_track, amb_prn_track] = this.remShortArc(A, y0, b, Q, obs_track, amb_prn_track, this.state.getMinArc());
            this.getBlockProperties();
        end
    end

    % ==================================================================================================================================================
    %  AUXILIARY FUNCTIONS goBlock
    % ==================================================================================================================================================

    methods % Public Access
        function plotPhRes (this, phase_res, amb_prn_track)
            if nargin == 1
                phase_res = this.phase_res;
                amb_prn_track = this.amb_prn_track;
            end
            this.plotDataTrack (phase_res*1e3, amb_prn_track);
            %ylim([-80 80]);
        end

        function plotDataTrack (this, data_track, amb_prn_track)
            single_plot = true;

%             if (nargin < 6)
%                 cs_plot = false;
%             end
            if nargin == 1
                data_track = this.phase_res;
                amb_prn_track = this.amb_prn_track;
%                 cs_plot = true;
            end

            x = (1 : size(data_track,1));
            if single_plot
                h = figure(); dockAllFigures(h);
            end
            for a = 1 : numel(amb_prn_track)
                if ~single_plot
                    h = figure(amb_prn_track(a));
                    h.Name = sprintf('Sat: %d', amb_prn_track(a));
                    h.NumberTitle = 'off';
                    f = subplot(4,8,amb_prn_track(a));
                end
                y = data_track(:, a, 1);
                plot(x, y,'.-', 'lineWidth', 1); hold on;
                hline = findobj(h, 'type', 'line');

                if (size(data_track, 3) == 2)
                    e = 3 * data_track(:, a, 2);
                    patchColor = min(hline(1).Color + 0.3, 1);
                    plot(x, e, x, -e, 'color', patchColor);
                    patch([x(~isnan(e)) fliplr(x(~isnan(e)))], [e(~isnan(e)); flipud(-e(~isnan(e)))], 1, ...
                        'facecolor',patchColor, ...
                        'edgecolor','none', ...
                        'facealpha', 0.1);
                end

%                 if (cs_plot)
%                     lambda_val = abs(A(id_track(~isnan(y),a), 3+a)) * 1e3 * this.cs_factor;
%                     win_size = this.state.getMinArc() + mod(1 + this.state.getMinArc(),2);
%                     win_size = round(win_size / 2) + mod(round(win_size / 2) + 1, 2);
%                     if numel(y(~isnan(y))) > win_size
%                         y(~isnan(y)) =  this.getAmbCorrection(y(~isnan(y)) ./ lambda_val, win_size) .* lambda_val;
%                         plot(x, y,':k', 'LineWidth', 1); hold on;
%                     end
%                 end
                if ~single_plot
                    f.Title.String = sprintf('Sat: %d', amb_prn_track(a));
                    ax(a) = gca; %#ok<AGROW>
                end
            end
            if ~single_plot
                linkaxes(ax);
            end
        end
    end

    % ==================================================================================================================================================
    %  GETTER FUNCTIONS goBlock
    % ==================================================================================================================================================

    methods % Public Access
        function toString(this)
            this.log.addMarkedMessage('Float solution:');
            this.log.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', ...
                [mean(this.getENU(this.getFloatPos()),1, 'omitnan')' mean(this.getDeltaENU(this.getFloatPos()),1, 'omitnan')']'));
            if this.is_fixed > 1
            this.log.addMarkedMessage('Fixed Position:');
            this.log.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', ...
                [mean(this.getENU(this.getFixPos()),1, 'omitnan')' mean(this.getDeltaENU(this.getFixPos()),1, 'omitnan')']'));
            end
            if this.is_fixed == 2
            this.log.addMarkedMessage('Position @ High Rate (mean):');
            this.log.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', ...
                [mean(this.getENU(this.getPosHR()),1, 'omitnan')' mean(this.getDeltaENU(this.getPosHR()),1, 'omitnan')']'));
            end
        end

        function [time_center, time_lim] = getTimeHR(this, s_rate)
            if nargin == 1
                s_rate = this.s_rate;
            end
            time_lim = unique([this.time_diff(1) : s_rate : this.time_diff(end) this.time_diff(end)]');
            time_lim = [time_lim(1 : end-1) time_lim(2 : end)];
            time_center = time_lim(:,1) + (time_lim(:,2) - time_lim(:,1)) / 2;
        end

        function [pos, pos_cov] = getPos(this)
            if (size(this.pos,2) > 1)
                % compute the mean of the computed positions

                % do not consider NaNs
                pos = this.pos;
                pos_cov = this.pos_cov(~isnan(pos(:)), ~isnan(pos(:)));
                pos = pos(:, ~isnan(sum(pos)));
                n_pos_hr = size(pos,2);

                % transformation matrix
                T = zeros(3, 3 * n_pos_hr);
                for i = 1 : n_pos_hr
                    T(1 : 3, (i - 1) * 3 + 1 : i * 3) = eye(3);
                end

                % mean position and covariance by LS
                [pos, pos_cov] = fast_least_squares_solver(pos(:), 0 * pos(:), T', pos_cov);
                pos = pos';
            else
                pos = this.pos';
                pos_cov = this.pos_cov;
            end
        end

        function pos = getFloatPos(this)
            pos = (repmat(this.pos0(:), 1, this.n_pos) + reshape(this.x_float(1:this.n_pos * 3), 3, this.n_pos))';
        end

        function [pos, pos_cov] = getFixPos(this)
            if ~(this.is_fixed)
                pos = this.getFloatPos();
                pos_cov = this.pos_cov;
            else
                pos = (repmat(this.pos0(:), 1, this.n_pos) + reshape(this.x_fix(1:this.n_pos * 3), 3, this.n_pos))';
                pos_cov = this.pos_cov_fix;
            end
        end

        function [pos, pos_cov] = getPosHR(this)
            if isempty(this.x_hr)
                pos = this.getFloatPos();
                pos_cov = this.pos_cov;
            else
                pos = (repmat(this.pos0(:), 1, this.n_pos_hr) + reshape(this.x_hr(1:this.n_pos_hr * 3), 3, this.n_pos_hr))';
                pos_cov = this.pos_cov;
            end
        end

        function ph_res = getPhaseResiduals(this)
            if isempty(this.phase_res)
                [this.phase_res, this.id_track] = this.computeDataTrack();
            end
            ph_res = nan(this.n_tot_epoch, this.state.cc.getNumSat());
            for i = 1 : size(this.phase_res, 2)
                ph_res(~isnan(this.phase_res(:, i, 1)), this.amb_prn_track(i)) = this.phase_res(~isnan(this.phase_res(:, i, 1)), i, 1);
            end
        end

        function [empty_epoch] = getEmptyEpochs(this)
            empty_epoch = this.empty_epoch;
        end

        function [pos_KAL, Xhat_t_t_OUT, conf_sat_OUT, Cee_OUT, pivot_OUT, nsat, fixed_amb] = getLegacyOutput(this)
            % [pos_KAL, Xhat_t_t_OUT, conf_sat_OUT, Cee_OUT, pivot_OUT, nsat, fixed_amb] = go_block.getLegacyOutput();
            pos_KAL = this.getPos()';
            [Xhat_t_t_OUT, Cee_OUT] = this.getFixPos();
            Xhat_t_t_OUT = Xhat_t_t_OUT';
            conf_sat_OUT = false(this.state.cc.getNumSat(),1);
            conf_sat_OUT(unique(this.amb_prn_track)) = true;
            pivot_OUT = this.pivot_track;
            nsat = this.state.cc.getNumSat();
            fixed_amb = this.is_fixed;
        end

        function [enu] = getENU(this, pos)
            % Coordinate transformation (UTM)
            % [enu] = this.getENU(pos)
            if nargin == 1
                pos = this.pos;
            end
            if size(pos,1) ~= 3
                pos = pos';
            end
            id_ok = find(~isnan(pos(1, :)));
            up = nan(size(pos(1, :)));
            east_utm = nan(size(pos(1, :)));
            north_utm = nan(size(pos(1, :)));
            [~, ~, up(id_ok)] = cart2geod(pos(1, id_ok), pos(2, id_ok), pos(3, id_ok));
            [east_utm(id_ok), north_utm(id_ok)] = cart2plan(pos(1, id_ok)', pos(2, id_ok)', pos(3, id_ok)');

            enu = [(east_utm(:))'; (north_utm(:))'; (up(:))']';
        end

        function [delta_enu] = getDeltaENU(this, pos)
            % Coordinate transformation (UTM)
            % [delta_enu] = this.getDeltaENU(pos)
            if nargin == 1
                pos = this.pos;
            end
            if size(pos,1) ~= 3
                pos = pos';
            end
            [~, ~, up0] = cart2geod(this.pos0(1, :), this.pos0(2, :), this.pos0(3, :));
            [east_utm0, north_utm0] = cart2plan(this.pos0(1, :)', this.pos0(2, :)', this.pos0(3, :)');

            id_ok = find(~isnan(pos(1, :)));
            up = nan(size(pos(1, :)));
            east_utm = nan(size(pos(1, :)));
            north_utm = nan(size(pos(1, :)));
            [~, ~, up(id_ok)] = cart2geod(pos(1, id_ok), pos(2, id_ok), pos(3, id_ok));
            [east_utm(id_ok), north_utm(id_ok)] = cart2plan(pos(1, id_ok)', pos(2, id_ok)', pos(3, id_ok)');

            delta_enu = [(east_utm0 - east_utm(:))'; (north_utm0 - north_utm(:))'; (up0 - up(:))']';
        end

        function [col_ok, ref_arc, block_cols, block_rows, fs_lim] = getBlockProperties(this, A, obs_track, n_pos)
            % retrive a set of arcs (columns of A) to be removed to compute a stable solution
            % SYNTAX:  [col_ok, ref_arc,  block_cols, block_rows] = this.getBlockProperties(<A>, <obs_track>, <n_pos>)
            %
            % INTERNAL OUTPUT:
            %   this.ref_arc
            %   this.col_ok
            %

            if (nargin == 1)
                A = this.A;
                obs_track = this.obs_track;
                n_pos =  this.n_pos;
            end

            % id of the full slip epochs (goBlock restart)
            id_fs = zeros(numel(this.empty_epoch), 1);
            for i = 1 : numel(this.empty_epoch)
                tmp = find(obs_track(:,1) < this.empty_epoch(i),1,'last');
                if ~isempty(tmp)
                    id_fs(i) = tmp;
                end
            end
            id_fs(id_fs == 0) = [];
            % goBlock independent computations
            fs_lim = [[1; unique(id_fs) + 1] [unique(id_fs); size(obs_track,1)]];
            fs_lim(fs_lim(:,2)-fs_lim(:,1) < 1, :) = [];
            ref_arc = nan(size(fs_lim, 1), 1);

            full_slip_split = this.state.getFullSlipSplit();
            tmp = reshape(this.obs_track(fs_lim(:), 1), size(fs_lim, 1),size(fs_lim, 2));
            tmp = tmp(2:end,1)-tmp(1:end-1,2) - 1;
            id_split = (find(tmp >= full_slip_split));

            block_cols = false(size(A, 2), size(id_split, 1) + 1);
            block_rows = false(size(A, 1), size(id_split, 1) + 1);
            b = 1;
            for i = 1 : size(fs_lim, 1)
                ref_arc(i) = find(sum(A(fs_lim(i,1) : fs_lim(i,2), 3 * n_pos + 1 : end) < 0) > 0, 1, 'last');
                block_cols(:, b) = block_cols(:, b) | [ true(1, 3 * n_pos) sum(abs(A(fs_lim(i,1) : fs_lim(i,2), 3 * n_pos + 1 : end)) ~= 0) > 0 ]';
                block_rows(fs_lim(i, 1) : fs_lim(i, 2), b) = true;
                if ismember(i, id_split)
                    b = b + 1;
                end
            end
            col_ok = setdiff(1 : size(A, 2), ref_arc + 3 * n_pos);
            if (nargin == 1)
                this.ref_arc = ref_arc;
                this.col_ok = col_ok;
            end
        end

    end

    % ==================================================================================================================================================
    %  STATIC LAUNCHERS goBlock
    % ==================================================================================================================================================

    methods (Static) % Public Access

        function [A_s, amb_prn_track] = shrinkA(A, amb_prn_track)
            [amb_prn_track, ~, id_reorder] = unique(amb_prn_track);
            n_amb = numel(amb_prn_track);
            A_s = [A(:, 1 : 3) sparse(size(A, 1), n_amb)];

            for a = 1 : size(id_reorder, 1)
                A_s(:, id_reorder(a) + 3) = A_s(:, id_reorder(a) + 3) + A(:, a + 3);
            end
        end

        function ref = getAmbCorrection(data, res, win_size)
            %ref1 = Core_Block.improveRef(data, movmedian(round(data), win_size, 'omitnan'), win_size);
            %ref2 = Core_Block.improveRef(data, movmedian(round(data+0.5) - 0.5, win_size, 'omitnan'), win_size);
            %ref = iif(sum(abs(data - ref1)) > sum(abs(data - ref2)), ref2, ref1);
            
            % reduce jumps on the first derivative
            ref1 = cumsum(round([0; diff(res)]));
            ref1 = ref1 + round(mean((res - ref1)));
            
            tmp = res - ref1;
            
            % reduce jumps on the first derivative (with a distance of "diff_dist" epochs)
            diff_dist = 2;
            ref2 = cumsum(round([zeros(diff_dist, 1); tmp(diff_dist + 1 : end) - tmp(1 : end - diff_dist)]) / diff_dist);
            
            ref = ref1 + ref2;
            
            % Correct if and only if there are jumps in the original data
            sensor = ~(abs([0; diff(data,2)]) > 0.9 & abs([0; diff(ref,2)]) > 0.5);
            lim = getOutliers(sensor);
            lim(:,1) = max(1,lim(:,1) - 1); lim(:,2) = lim(:,2) + 1;
            for i = 1 : size(lim, 1)
                idx = lim(i, 1) : lim(i, 2);
                if mod(length(idx),2); idx = idx(1 : end - 1); end
                ref(idx) = median(ref(idx));
            end
            
            % test ref on data
            ref = Core_Block.improveRef(res, ref, win_size);
            
            % compute second derivative
            d2_ref = [0; diff(data - ref,2)];
            d2_data = [0; diff(data,2)];
            
            % choose the solution that minimize the speed of change
            ref_best = abs(d2_ref) < abs(d2_data);
            ref = ref(1) + cumsum([0; diff(ref) .* ref_best]);
        end
        
        function ref = improveRef(data, ref, win_size)
            margin_win = round(win_size / 2) + mod(round(win_size / 2) + 1, 2);

            % check where the interpolation is not a plane
            lim = getOutliers(full(abs(movmedian(diff([ref(1); ref]),3)) > 0));
            lim(lim(:,2) - lim(:,1) < 7, :) = [];
            for i = 1 : size(lim, 1)
                idx = lim(i, 1) : lim(i, 2);
                idx = idx(1 : end - mod(length(idx),2));
                ref(idx) = median(ref(idx(1 : end)));
            end

            % reduce the number of jumps
            lim = unique([1; find([abs(diff(data)) > 0.4; true])]);
            %lim(diff(lim) < round(win_size / 2)) = [];
            lim = unique([1; lim; numel(ref)]);
            lim = [lim(1:end-1)+1 lim(2:end)];
            for i = 1 : size(lim, 1)
                idx = lim(i, 1) : lim(i, 2);
                if mod(length(idx),2); idx = idx(1 : end - 1); end
                ref(idx) = median(ref(idx));
            end
            ref = movmedian(ref, 2 * win_size -1);

            ref(1 : margin_win) = median(ref(1 : margin_win));
            ref(end - margin_win + 1 : end) = median(ref(end - margin_win + 1 : end));
        end

        function go_block = go (time_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, sp3, iono, lambda, ant_pcv, s_rate)
            % Separate the dataset in single blocks and solve float -> fix
            % go_block = Core_Block.go(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV, 3600);
            if nargin == 18
                s_rate = numel(time_diff);
            end
            go_block = Core_Block (numel(time_diff), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            go_block.prepare (time_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, sp3, iono, lambda, ant_pcv);
            go_block.solve(s_rate);
        end

        function go_block = goMultiHighRate(time_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, sp3, iono, lambda, ant_pcv, s_rate)
            % Separate the dataset in single blocks and solve float -> fix
            % go_block = Core_Block.goMultiHighRate(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV, 3600);
            %%
            state = Core.getCurrentSettings();
            log = Core.getLogger();

            s_rate = s_rate(1);
            idx = (unique([0 : s_rate : max(time_diff) max(time_diff)]) / state.getProcessingRate)';
            n_step = length(idx)-1;
            idx = [idx(1 : end - 1)+1 idx(2 : end)];
            idx(end, 2) = numel(time_diff());

            pos_hr = zeros(n_step, 3, 1);
            pos_cov_hr = zeros(3 * n_step, 3 * n_step);
            pos_fix = zeros(n_step, 3, 1);
            pos_float = zeros(n_step, 3, 1);

            for i = 1 : n_step
                log.addMessage(sprintf('\nProcessing HR solution: %02d/%02d ----------------------------------\n', i, n_step));
                go_block = Core_Block (numel(time_diff(idx(i, 1) : idx(i, 2))), sum(serialize(pr1_R(:,idx(i, 1) : idx(i, 2)) ~= 0)), sum(serialize(ph1_R(:,idx(i, 1) : idx(i, 2)) ~= 0)));

                go_block.prepare (time_diff(idx(i, 1) : idx(i, 2)), pos_R, pos_M, pr1_R(:,idx(i, 1) : idx(i, 2)), pr1_M(:,idx(i, 1) : idx(i, 2)), pr2_R(:,idx(i, 1) : idx(i, 2)), pr2_M(:,idx(i, 1) : idx(i, 2)), ...
                    ph1_R(:,idx(i, 1) : idx(i, 2)), ph1_M(:,idx(i, 1) : idx(i, 2)), ph2_R(:,idx(i, 1) : idx(i, 2)), ph2_M(:,idx(i, 1) : idx(i, 2)), ...
                    snr_R(:,idx(i, 1) : idx(i, 2)), snr_M(:,idx(i, 1) : idx(i, 2)),  Eph, sp3, iono, lambda, ant_pcv);

                go_block.solve([], false);

                [pos, pos_cov] = go_block.getFixPos();
                pos_hr(i,:,1) = pos;
                pos_cov_hr(1 + (i - 1) * 3 : (i * 3) , 1 + (i - 1) * 3 : (i * 3)) = pos_cov;
                pos_fix(i,:,1) = go_block.getFixPos();
                pos_float(i,:,1) = go_block.getFloatPos();
            end

            % Float and fix solutions are daily position in seamless mode => store the mean
            pos_float = pos_float';
            pos_fix = pos_fix';
            pos_fix = pos_fix(:, ~isnan(sum(pos_fix)));
            pos_float = pos_float(:, ~isnan(sum(pos_fix)));
            n_pos_hr = size(pos_hr, 1);

            % transformation matrix
            T = zeros(3, 3 * n_pos_hr);
            for i = 1 : n_pos_hr
                T(1 : 3, (i - 1) * 3 + 1 : i * 3) = eye(3) / n_pos_hr;
            end
            pos_fix = T * pos_fix(:);
            pos_float = T * pos_float(:);

            go_block.importPartialResults(time_diff, pos_hr', pos_cov_hr, pos_fix, pos_float);
        end

        function go_block = getCopyFrom(obj)
            % get a copy of a Core_Block instance
            % SYNTAX: go_block = Core_Block.getCopyFrom(obj);
            go_block = Core_Block(0,0,0);
            go_block.import(obj);
        end
    end

    % ==================================================================================================================================================
    %  STATIC public goBlock misc utilities
    % ==================================================================================================================================================

    methods (Static) % Public Access

        function [x, Cxx, s02, v_hat, P_out, N_out, Cyy] = solveLS(y0, b, A, col_ok, iQ, P, N)
            % Solve the LS problem, when P, N are provided it does not compute them
            % SYNTAX: [x, Cxx, s02, v_hat, P_out, N_out, Cyy] = Core_Block.solveLS(y0, b, A, col_ok, cholinv(Q), P, N)
            [n_obs, n_col] = size(A);

            if nargout < 5
                A = A(:, col_ok);
            end

            % least-squares solution
            if (nargin < 7) || isempty(P)
                P = A' * iQ;
                N = P * A;
            end

            if (nargout > 4)
                P_out = P;
                N_out = N;
            end

            if numel(col_ok) < size(A, 2)
                P = P(col_ok, :);
                N = N(col_ok, col_ok);
                A = A(:, col_ok);
            end

            if numel(y0) > 0
                try
                    N_inv = cholinv(N);
                catch
                    N_inv = N^-1;
                end

                Y = (y0 - b);
                L = P * Y;
                x = N_inv * L;

                % estimation of the variance of the observation error
                y_hat = A * x + b;
                v_hat = y0 - y_hat;
                s02 = (v_hat' * iQ * v_hat) / (n_obs - n_col);

                % covariance matrix
                Cxx = s02 * N_inv;

                if nargout == 7
                    Cyy = s02 * A * N_inv * A';
                end
            else
                x = []; Cxx = []; s02 = 1e100; v_hat = []; P_out = []; N_out = []; Cyy = [];
            end

        end

        function [x, Cxx, s02, v_hat, P_out, N_out, Cyy] = solveRegLS(y0, b, A, col_ok, iQ, R, P, N)
            % Solve the LS problem, when P, N are provided it does not compute them
            % SYNTAX: [x, Cxx, s02, v_hat, P_out, N_out, Cyy] = Core_Block.solveLS(y0, b, A, col_ok, cholinv(Q), P, N)
            [n_obs, n_col] = size(A);

            if nargout < 5
                A = A(:, col_ok);
            end

            % least-squares solution
            if (nargin < 8) || isempty(P)
                P = A' * iQ;
                N = P * A;
            end

            if (nargout > 4)
                P_out = P;
                N_out = N;
            end

            if numel(col_ok) < size(A, 2)
                P = P(col_ok, :);
                N = N(col_ok, col_ok);
                A = A(:, col_ok);
            end

            N = N + R;
            % small regularization to make the system more stable
            reg_f = abs(diag(N)); reg_f = min(reg_f(reg_f>0));
            N_inv = cholinv((N + 1e-6*reg_f * speye(size(N,1))));

            Y = (y0 - b);
            L = P * Y;
            x = N_inv * L;

            % estimation of the variance of the observation error
            y_hat = A * x + b;
            v_hat = y0 - y_hat;
            T = iQ * v_hat;
            s02 = (v_hat' * T) / (n_obs - n_col);

            % covariance matrix
            Cxx = s02 * N_inv;

            if nargout == 7
                Cyy = s02 * A * N_inv * A';
            end
        end

    end

    % ==================================================================================================================================================
    %  PRIVATE FUNCTIONS called by public calls goBlock
    % ==================================================================================================================================================

    methods (Access = public)

        function [A, col_ok, ref_arc, y0, b, Q, obs_track, amb_prn_track] = remArcs(this, A, col_ok, ref_arc, y0, b, Q, obs_track, amb_prn_track, Cxx)
            % If the system is unstable remove the arcs that are making it so
            amb_var = zeros(size(Cxx,1) + numel(ref_arc)); amb_var(col_ok) = diag(Cxx); amb_var(amb_var < 0) = 100;
            % Bad arcs have an high estimation variance
            bad_col = find((amb_var(4:end) > 0.5) | (amb_var(4:end) > mean(amb_var(col_ok(4 : end))) + 10 * std(amb_var(col_ok(4 : end))))) + 3;
            Q_tmp = Q;
            while ~isempty(bad_col) && (size(A,2) >= 3 + 2)
                bad_col = bad_col(1);
                this.log.addMessage(sprintf('       - System found unstable, removing arc %d - prn %d', bad_col - 3, amb_prn_track(bad_col - 3)));
                [A, y0, b, Q, obs_track, amb_prn_track] = this.remArcCol(A, y0, b, Q, obs_track, amb_prn_track, bad_col);
                [A, y0, b, Q, obs_track, amb_prn_track] = this.remShortArc(A, y0, b, Q, obs_track, amb_prn_track, this.state.getMinArc());
                [col_ok, ref_arc] = this.getBestRefArc(y0, b, A, Q);
                [~, Cxx, ~, ~, Q_tmp] = this.improveFloatSolution(y0, b, A, col_ok, Q, [], obs_track);
                amb_var = zeros(size(Cxx,1) + numel(ref_arc)); amb_var(col_ok) = diag(abs(Cxx));
                bad_col = find((amb_var(4:end) > 10) | (amb_var(4:end) > mean(amb_var(col_ok(4 : end))) + 10 * std(amb_var(col_ok(4 : end))))) + 3;
            end
            Q = Q_tmp;
        end

        function [x_float, Cxx, s02, v_hat, y0, Q, iQ, phase_res, id_track] = loopCorrector(this, y0, b, A, col_ok, Q, obs_track, amb_prn_track, n_clean)
            % Try to correct integer ambiguities (maybe missed cycle slips) on the basis of the residuals
            % SYNTAX: [x_float, Cxx, s02, v_hat, y0, Q, iQ, phase_res, id_track] = this.loopCorrector(y0, b, A, col_ok, Q, obs_track, amb_prn_track, n_clean)
            narginchk(9,9);
            [x_float, Cxx, s02, v_hat, Q_tmp, iQ] = this.improveFloatSolution(y0, b, A, col_ok, Q, [], obs_track);
            %[A_s, amb_prn_track] = this.shrinkA(A, amb_prn_track);
            [phase_res, id_track] = this.computeDataTrack( v_hat, [], A, obs_track, 1);
            iQ_tmp = iQ;
            c = 1; is_new = true;
            while (c <= n_clean) && is_new
                % Try to correct integer ambiguities (maybe missed cycle slips)
                this.log.addMessage(sprintf('       - try to fix previously undetected cycle slips %d', c));
                win_size = (this.state.getMinArc()/2 + mod(1 + this.state.getMinArc()/2, 2));
                [y0_ck, is_new_ck] = this.postCorrectIntAmb(y0, phase_res, id_track, A, amb_prn_track, win_size);

                if is_new_ck
                    % Improve solution by iterative increase of bad observations variance
                    %this.log.addMessage('          - recompute the improved solution');
                    [x_float_ck, Cxx_ck, s02_ck, v_hat_ck, Q_tmp_ck, iQ_tmp_ck] = this.improveFloatSolution(y0_ck, b, A, col_ok, Q, [], obs_track);
                    % if I'm improving the solution accept the correction
                    [phase_res_ck] = this.computeDataTrack(v_hat_ck, id_track);
                    [phase_res] = this.computeDataTrack(v_hat, id_track);
                    id_good = mean(phase_res.^2, 2, 'omitnan') > mean(phase_res_ck.^2, 2, 'omitnan');
                    if sum(id_good > 0)
                        id_good = serialize(id_track(id_good,:));
                        id_good(id_good == 0) = [];
                        is_new = true;
                        y0(id_good) = y0_ck(id_good);
                        [x_float, Cxx, s02, v_hat, Q_tmp, iQ_tmp] = this.improveFloatSolution(y0, b, A, col_ok, Q, [], obs_track);
                        [phase_res] = this.computeDataTrack(v_hat, id_track);
                    else
                        this.log.addMessage('          - cycle slips detected but rejected!');
                        is_new = false;
                    end
                else
                    this.log.addMessage('          - no cycle slips have been found!');
                    is_new = false;
                end
                c = c + 1;
            end
            Q = Q_tmp;
            iQ = iQ_tmp;
        end

        function [id_out, phase_res] = phaseCleaner(this, phase_res, id_track, thr_fix)
            % clean phase residuals and suggest outliers
            % SYNTAX: [id_out, phase_res] = phaseCleaner(this, phase_res, id_track, thr_fix);
            narginchk(3,4);
            if nargin == 3
                thr_fix = this.state.getMaxPhaseErrThr();
            end

            win_size = this.state.getMinArc() + mod(1 + this.state.getMinArc(),2);
            half_win_size = round(this.state.getMinArc() / 2) + mod(1 + round(this.state.getMinArc() / 2),2);
            id_out = [];
            for i = 1 : size(phase_res, 2)
                phase = phase_res(:,i, 1);
                phase_s = splinerMat([],movmedian(phase, win_size, 'omitnan'),half_win_size, 0.01);
                thr = 15 * median(movstd(phase-phase_s, 11, 'omitnan'), 'omitnan');
                id = find(abs(phase - phase_s) > thr | abs(phase) > thr_fix);
                phase_res(id, i, 1) = NaN;
                id_out = [id_out; id_track(id, i)]; %#ok<AGROW>
            end
        end

        function importPartialResults(this, time_diff, pos_hr, pos_cov, pos_fix, pos_float)
            this.time_diff = time_diff;
            this.is_fixed = 2;
            this.s_rate = this.state.getSolutionRate();
            this.n_pos_hr = numel(pos_hr)/3;
            this.pos = pos_hr;
            this.pos_cov = pos_cov;
            this.x_hr = pos_hr(:) - repmat(this.pos0(:), this.n_pos_hr, 1);
            this.x_fix = pos_fix - repmat(this.pos0(:), this.n_pos, 1);
            this.x_float = pos_float - repmat(this.pos0(:), this.n_pos, 1);
            this.n_pos = numel(pos_float)/3;
            this.amb_fix = [];
            this.amb_fix_full = [];
        end

        function addAmbiguities (this, lambda)
            % Add to the internal Design Matrix the columns related to the phase observations
            % (Integer ambiguities - N)
            %
            % METHODS CALL REQUIREMENTS:
            %  -> prepare
            %
            % SYNTAX:
            %   this.addAmbiguities(lambda)
            %
            % INPUT:
            %   lambda     wavelength matrix (depending on the enabled constellations)
            %
            % INTERNAL INPUT:
            %   A, amb_prn_track, n_pos, n_obs_tot, sat_ph_track
            %
            % INTERNAL OUTOUT (properties whos values are changed):
            %   A, amb_prn_track
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);

            this.log.addMarkedMessage('Set up Design Matrix to estimate integer abiguities');

            amb_num = 0;
            this.amb_prn_track = [];

            if this.state.getFullSlipSplit()
                full_slip = this.empty_epoch(:) - (0 : length(this.empty_epoch) - 1)';
                full_slip = full_slip(diff([full_slip(:)]) > 1);
            else
                full_slip = [];
            end

            if (this.state.isModeSA())
                % set the design matrix to estimate the receiver clock
                this.A(:, 3 + 1) = 1;

                % clear A avoid problems when launching addAmbiguities() more than once
                this.A = this.A(:, 1 : 3 + 1);
            else
                % clear A avoid problems when launching addAmbiguities() more than once
                this.A = this.A(:, 1 : 3);
                if (this.state.isModePh())
                    % matrix to detect re-initialized ambiguities
                    % every time there is an interruption in the satellite observations -> suppose cs
                    amb_track = [int8(zeros(size(this.sat_ph_track,1), 1)) diff(this.sat_ph_track')'];
                    amb_track(:,1) = this.sat_ph_track(:,1) ~= 0;
                    % Cycle slip at pivot change
                    % pivot_change = (abs(diff([int8(this.pivot_track(1)) int8(this.pivot_track')])) > 0);
                    % amb_track = amb_track ...
                    %     + int8(this.sat_ph_track & repmat(pivot_change,size(this.sat_ph_track,1),1)) ...
                    %     - int8(this.sat_ph_track & repmat([pivot_change(2:end) pivot_change(1)],size(this.sat_ph_track,1),1));
                    if this.state.getFullSlipSplit()
                        for i = 1 : length(full_slip)
                            amb_track(this.sat_ph_track(:,full_slip(i)) ~= 0, full_slip(i)) = 1;
                        end
                    end

                    % resize the design matrix to estimate phase ambiguities
                    %this.A = [this.A sparse(this.n_obs_tot, amb_num)];
                    rows = 0;
                    this.A = [this.A spalloc(this.n_obs_tot, 10000, this.n_obs_tot * 2)];
                    for e = 1 : this.n_epoch

                        % check if a new ambiguity column for A is needed

                        % detect new ambiguities on epoch e
                        amb_prn_new = find(amb_track(:,e) == 1);
                        if (~isempty(amb_prn_new))
                            amb_num = amb_num + length(amb_prn_new);
                            this.amb_prn_track = [this.amb_prn_track; amb_prn_new];
                        end

                        % build new columns
                        pivot_prn = this.pivot_track(e);
                        this.sat_pr_track(pivot_prn, e) = -1 * abs(this.sat_pr_track(pivot_prn, e));
                        this.sat_ph_track(pivot_prn, e) = -1 * abs(this.sat_ph_track(pivot_prn, e));

                        amb_prn_avail = find(this.sat_ph_track(:,e) == 1);
                        [~, amb_id] = intersect(flipud(this.amb_prn_track), amb_prn_avail);
                        amb_id = numel(this.amb_prn_track) - amb_id + 1;

                        pivot_id = find(this.amb_prn_track == pivot_prn, 1, 'last');
                        %sat_pr_idx = sum(this.sat_pr_track(:,e) == 1) + (1 : sum(this.sat_ph_track(:, e) == 1));
                        sat_ph_idx = sum(this.sat_pr_track(:,e) > 0) + (1 : sum(this.sat_ph_track(:, e) > 0))';
                        rows = rows(end) + sat_ph_idx;
                        this.A(rows(:) + size(this.A,1) * (3 + amb_id -1)) = -lambda(amb_prn_avail, 1);
                        this.A(rows, 3 + pivot_id) = lambda(pivot_prn, 1);
                    end
                    %[this.amb_prn_track, reorder_id] = sort(this.amb_prn_track);
                    %this.A = this.A(:, [(1 : (3))'; (3) + reorder_id]);
                    this.A = this.A(:, 1 : amb_num + 3 * this.n_pos);
                end
            end
        end

        function [d_pos, pos_cov, is_fixed, amb_fix, amb_cov, amb_fix_full, ref_arc, G] = solveFixPar (this, x_float, Cxx, amb_num)
            % Compute a fixed solution using LAMBDA, and the the internal object properties (Concrete implementation)
            %
            % METHODS CALL REQUIREMENTS:
            %   prepare -> addAmbiguities -> solveFloat
            %
            % SYNTAX:
            %   [d_pos, pos_cov, amb_fix, amb_cov, amb_fix_full, ref_arc] = this.solveFix()
            %
            % INPUT:
            %   x_float     Float solution
            %   Cxx         Covariance matrix of the solution
            %   amb_num     Ambiguities number
            %
            % INTERNAL INPUT:
            %   log
            %
            % OUTPUT:
            %   d_pos       coordinates offset of the estimated positions
            %   pos_cov         covariance of the estimated positions
            %   is_fixed        flag is fixed?
            %   amb_fix         ambiguities as estimated by lambda (n-1 w.r.t. float solution)
            %   amb_cov         ambiguities error covariance matrix
            %   amb_fix_full    ambiguities as converted from fix to float -> to be imported as pseudo observations of the float solution
            %   ref_arc         arc used as reference in the fix solution (it's the arc that create a bias in the solution)
            %            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.addAmbiguities(lambda)
            %   go_block.solveFloat()
            %   go_block.solveFix()
            %

            % Compute a first float solution
            if (nargin < 3) || (isempty(x_float))
                this.solveFloat();
                x_float = this.x_float;
                Cxx = this.Cxx;
                amb_num = size(this.A, 2) - 3 - numel(this.ref_arc);
            end

            % switch from SD to DD
            D = zeros(amb_num - 1, amb_num);
            % find best estimated arc:
            %[~, ref_arc] = sort(diag(Cxx(4 : end, 4 : end))); ref_arc = ref_arc(1);
            ref_arc = 1; % there are no differences choosing one arc or another (~1e-10 differences in the results)
            D(:, ref_arc) = 1;
            D(:, [1 : ref_arc - 1, ref_arc + 1 : end]) = -eye(amb_num - 1);
            G = zeros(3 + amb_num - 1, 3 + amb_num);
            G(1 : 3, 1 : 3) = eye(3);
            G(4 : end, 4 : end) = D;
            x = G * x_float;
            Cxx = full(G * Cxx * G');

            cov_X  = Cxx(1 : 3, 1 : 3);     % position covariance block
            cov_N  = Cxx(4 : end, 4 : end); % ambiguity covariance block
            cov_XN = Cxx(1 : 3, 4 : end);   % position-ambiguity covariance block

            is_fixed = 0;
            try
                try
                    % stabilize cov_N;
                    [U] = chol(cov_N);
                    cov_N = U'*U;
                catch ex
                    this.log.addWarning(sprintf('Phase ambiguities covariance matrix unstable - %s', ex.message));
                end

                % integer phase ambiguity solving by LAMBDA
                [d_pos_ck, amb_fix_ck, amb_cov, pos_cov, d_pos, amb_fix] = lambdafix(x(1:3), x(4:end), cov_X, cov_N, cov_XN);

                if sum(amb_fix_ck == x(4:end))
                    this.log.addWarning('LAMBDA returned a fixed solution that did not pass the ratio test\nTring to use the solution anyway');
                end

                n_cands = size(amb_fix, 2);

                amb_fix_full = zeros(amb_num, n_cands);
                amb_fix_full([1 : ref_arc - 1, ref_arc + 1 : end],:) = -(amb_fix - repmat(x_float(ref_arc + 3), 1, n_cands));
                amb_fix_full(ref_arc,:) = repmat(x_float(3 + ref_arc), 1, n_cands);

                x_new = [repmat(x(1 : 3), 1, n_cands); amb_fix_full];
                y_hat = this.A(:, this.col_ok)  * x_new + repmat(this.b, 1, n_cands);
                v_hat = repmat(this.y0, 1, n_cands) - y_hat;
                T = this.Q \ v_hat;
                s02 = diag((v_hat' * T) / (size(this.A, 1) - numel(this.col_ok)));
                [~, id_best] = sort(s02); id_best = id_best(1);

                amb_fix = amb_fix(:, id_best);
                d_pos = d_pos(:, id_best);
                amb_fix_full = amb_fix_full(:, id_best);
                is_fixed = 1;
            catch ex
                this.log.addWarning(sprintf('It was not possible to estimate integer ambiguities: a float solution will be output.\n%s',ex.message));
                pos_cov = cov_X;
                amb_cov = cov_N;
                amb_fix_full = x_float(4 : end);
                amb_fix = amb_fix_full;
                d_pos = x_float(1 : 3);
                ref_arc = 0;
            end
        end

        function [x_float, Cxx, s02, v_hat, Q, iQ] = improveFloatSolution(this, ...
                    y0, b, A, col_ok, Q, ...
                    v_hat, obs_track, thr)
            % Stabilize solution by increasing bad observations variances
            %
            % SYNTAX:
            %   [x_float, Cxx, s02, v_hat, Q, iQ] = this.improveFloatSolution(this, y0, b, A, col_ok, Q, v_hat, obs_track)

            x_float = [];
            iQ = cholinv(Q);
            if isempty(v_hat)
                [x_float, Cxx, s02, v_hat] = this.solveLS(y0, b, A, col_ok, iQ);
            end

            narginchk(7,8);
            if nargin == 8
                thr = 0;
            end

            search_for_outlier = 1;

            out_ph_old = false(size(obs_track,1),1);
            out_pr_old = false(size(obs_track,1),1);
            n_out_old = 0;
            idx_pr = find(obs_track(:,3) == -1);
            idx_ph = find(obs_track(:,3) == 1);
            thr_perc = 0.995;
            while (search_for_outlier == 1)
                % never remove more than 0.5% of data at time
                out_pr = abs(v_hat(idx_pr)) > max(0.5*thr, max(this.state.getMaxCodeErrThr() , perc(abs(v_hat(~out_pr_old)), thr_perc)));
                if isempty(out_pr)
                    out_pr = out_pr_old;
                end
                out_ph = abs(v_hat(idx_ph)) > max(0.5*thr, max(this.state.getMaxPhaseErrThr() , perc(abs(v_hat(~out_ph_old)), thr_perc)));
                if isempty(out_ph)
                    out_ph = out_ph_old;
                end
                n_out = sum(out_pr_old | out_pr | out_ph_old | out_ph);
                idx_out_pr = idx_pr(out_pr);
                idx_out_ph = idx_ph(out_ph);
                if n_out_old < n_out
                    n_out_old = n_out;
                    out_pr_old = out_pr_old | out_pr;
                    out_ph_old = out_ph_old | out_ph;
                    idx_out = [idx_out_pr idx_out_ph];
                    Q(idx_out + size(Q,1) * (idx_out - 1)) = min(1e4, max(Q(idx_out + size(Q,1) * (idx_out - 1)), (v_hat(idx_out)).^2)); % Bad observations have now their empirical error
                    %Q(idx_out + size(Q,1) * (idx_out - 1)) = min(1e4, (v_hat(idx_out)).^2); % Bad observations have now their empirical error
                    iQ = cholinv(Q);
                    [x_float, Cxx, s02, v_hat] = this.solveLS(y0, b, A, col_ok, iQ);
                else
                    if isempty(x_float)
                        [x_float, Cxx, s02, v_hat] = this.solveLS(y0, b, A, col_ok, iQ);
                    end
                    search_for_outlier = 0;
                end
                thr_perc = thr_perc * 0.95;
            end
        end

        function [x_float, Cxx, s02, v_hat, y0, b, A, col_ok, Q, obs_track, amb_prn_track, iQ]  = cleanObsHiRes (this, ...
                    y0, b, A, col_ok, Q, ...
                    v_hat, obs_track, amb_prn_track, thr, force_out)
            % Stabilize solution by removing observations eith high residuals
            % [x_float, Cxx, s02, v_hat, y0, b, A, Q, obs_track, amb_prn_track] = cleanObsHiRes(this, y0, b, A, Q, v_hat, obs_track, amb_prn_track)
            search_for_outlier = 1;

            if nargin < 10
                thr = 3;
            end
            if nargin < 11
                force_out = true;
            end
            iQ = cholinv(Q);
            [x_float, Cxx, s02, v_hat] = this.improveFloatSolution(y0, b, A, col_ok, Q, v_hat, obs_track);
            v_hat_ck = v_hat;
            A_ck = A;
            y0_ck = y0;
            b_ck = b;
            Q_ck = Q;
            obs_track_ck = obs_track;
            amb_prn_track_ck = amb_prn_track;
            col_ok_ck = col_ok;

            while (search_for_outlier == 1)
                idx_pr = find(obs_track_ck(:,3) == -1);
                idx_ph = find(obs_track_ck(:,3) == 1);
                % never remove more than 0.5% of data at time
                pr_ok = abs(v_hat_ck(idx_pr)) < max(min(this.state.getMaxCodeErrThr() , thr * sqrt(Q(idx_pr + size(Q,1) * (idx_pr - 1)))), perc(abs(v_hat_ck), 0.995));
                ph_ok = abs(v_hat_ck(idx_ph)) < max(min(this.state.getMaxPhaseErrThr() , thr * sqrt(Q(idx_ph + size(Q,1) * (idx_ph - 1)))), perc(abs(v_hat_ck), 0.995));
                idx_pr_ok = idx_pr(pr_ok);
                idx_ph_ok = idx_ph(ph_ok);
                idx_ok = [idx_pr_ok idx_ph_ok];
                if (numel(idx_ok) < numel(y0_ck))
                    n_col = size(A, 2);
                    [A_ck, y0_ck, b_ck, Q_ck, obs_track_ck, amb_prn_track_ck] = this.remShortArc(A_ck(idx_ok,:), y0_ck(idx_ok), b_ck(idx_ok), Q_ck(idx_ok, idx_ok), obs_track_ck(idx_ok,:), amb_prn_track_ck, this.state.getMinArc());
                    % If I removed an arc re-estimate the best arc
                    iQ_ck = cholinv(Q_ck);
                    if size(A_ck, 2) < n_col
                        [col_ok_ck] = this.getBestRefArc(y0_ck, b_ck, A_ck, iQ_ck);
                    end
                    if numel(y0_ck) > 0
                        %[x_float_ck, Cxx_ck, s02_ck, v_hat_ck] = this.solveLS(y0_ck, b_ck, A_ck, col_ok_ck, iQ_ck);
                        [x_float_ck, Cxx_ck, s02_ck, v_hat_ck] = this.improveFloatSolution(y0_ck, b_ck, A_ck, col_ok_ck, Q_ck, v_hat_ck, obs_track_ck);
                        % keep the cleaning only if the solution is improving
                        if force_out || (sum(diag(Cxx_ck(1:3,1:3)).^2) <= sum(diag(Cxx(1:3,1:3)).^2))
                            x_float = x_float_ck;
                            Cxx = Cxx_ck;
                            s02 = s02_ck;
                            v_hat = v_hat_ck;
                            A = A_ck;
                            y0 = y0_ck;
                            b = b_ck;
                            Q = Q_ck;
                            iQ = iQ_ck;
                            obs_track = obs_track_ck;
                            amb_prn_track = amb_prn_track_ck;
                            col_ok = col_ok_ck;
                        else
                            search_for_outlier = 0;
                        end
                    end
                else
                    search_for_outlier = 0;
                end
            end
        end

        function [x_float, Cxx, s02, v_hat, y0, b, A, col_ok, Q, obs_track, amb_prn_track, iQ]  = cleanObsHiResVar (this, ...
                y0, b, A, col_ok, Q, ...
                v_hat, obs_track, amb_prn_track, thr)
            % Stabilize solution by removing observations eith high residuals
            % [x_float, Cxx, s02, v_hat, y0, b, A, Q, obs_track, amb_prn_track] = cleanObsHiResVar(this, y0, b, A, Q, v_hat, obs_track, amb_prn_track)

            if nargin < 10
                thr = 0.997;
            end

            search_for_outlier = 1;
            iQ = cholinv(Q);
            [x_float, Cxx, s02, v_hat] = this.improveFloatSolution(y0, b, A, col_ok, Q, v_hat, obs_track);
            while (search_for_outlier == 1)
                idx_pr = find(obs_track(:,3) == -1);
                [phase_res, id_track] = this.computeDataTrack(v_hat, [], A, obs_track, 1);
                sensor = movstd([zeros(1, size(phase_res, 2)); diff(phase_res)],3,'omitnan');
                thr_var = min(this.state.getMaxPhaseErrThr() * ones(size(sensor,1),1), min(median(sensor,2, 'omitnan'), this.state.getMaxPhaseErrThr()) + perc(median(sensor,2, 'omitnan'), thr));
                idx_ph = [];
                for a = 1 : size(id_track,2)
                    idx_ph = [idx_ph; id_track(sensor(:, a) < thr_var & ~(id_track(:,a) == 0), a)];
                end
                idx_ok = [idx_pr idx_ph];

                idx_ph = find(obs_track(:,3) == 1);
                % never remove more than 0.5% of data at time
                ph_ok = abs(v_hat(idx_ph)) < max(max(this.state.getMaxPhaseErrThr() , thr * sqrt(Q(idx_ph + size(Q,1) * (idx_ph - 1)))), perc(abs(v_hat), 0.995));
                idx_ph_ok = idx_ph(ph_ok);

                idx_ok = intersect(idx_ok, idx_ph_ok);
                if (numel(idx_ok) < numel(y0))
                    col_ok_ck = col_ok;
                    n_col = size(A, 2);
                    [A_ck, y0_ck, b_ck, Q_ck, obs_track_ck, amb_prn_track_ck] = this.remShortArc(A(idx_ok,:), y0(idx_ok), b(idx_ok), Q(idx_ok, idx_ok), obs_track(idx_ok,:), amb_prn_track, this.state.getMinArc());
                    if ~isempty(A_ck)
                        % If I removed an arc re-estimate the best arc
                        iQ_ck = cholinv(Q_ck);
                        if size(A_ck, 2) < n_col
                            [col_ok_ck] = this.getBestRefArc(y0_ck, b_ck, A_ck, iQ_ck);
                            n_col = size(A_ck, 2);
                            [~, ~, ~, ~, y0_ck,  b_ck, A_ck, ~, Q_ck, obs_track_ck, amb_prn_track_ck, iQ_ck] = this.remSolitaryObs(y0_ck, b_ck, A_ck, col_ok_ck, Q_ck, obs_track_ck, amb_prn_track_ck, round(this.state.getMinArc()));
                            if size(A_ck, 2) < n_col
                                [col_ok_ck] = this.getBestRefArc(y0_ck, b_ck, A_ck, iQ_ck);
                            end
                        end

                        if numel(y0) > 0
                            %[x_float_ck, Cxx_ck, s02_ck, v_hat_ck] = this.solveLS(y0_ck, b_ck, A_ck, col_ok_ck, iQ_ck);
                            [x_float_ck, Cxx_ck, s02_ck, v_hat_ck] = this.improveFloatSolution(y0_ck, b_ck, A_ck, col_ok_ck, Q_ck, [], obs_track_ck);
                            %[x_float_ck, Cxx_ck, s02_ck, v_hat_ck, y0_ck, Q_ck] = this.loopCorrector(y0_ck, b_ck, A_ck, col_ok_ck, Q_ck, obs_track_ck, amb_prn_track_ck, 3);
                            % keep the cleaning only if the solution is improving
                            if sum(diag(Cxx_ck(1:3,1:3)).^2) <= sum(diag(Cxx(1:3,1:3)).^2)
                                x_float = x_float_ck;
                                Cxx = Cxx_ck;
                                s02 = s02_ck;
                                v_hat = v_hat_ck;
                                A = A_ck;
                                y0 = y0_ck;
                                b = b_ck;
                                Q = Q_ck;
                                iQ = iQ_ck;
                                obs_track = obs_track_ck;
                                amb_prn_track = amb_prn_track_ck;
                                col_ok = col_ok_ck;
                            else
                                search_for_outlier = 0;
                            end
                        end
                    else
                        x_float = [];
                        Cxx = [];
                        s02 = [];
                        v_hat = [];
                        A = A_ck;
                        y0 = y0_ck;
                        b = b_ck;
                        Q = Q_ck;
                        iQ = [];
                        obs_track = obs_track_ck;
                        amb_prn_track = amb_prn_track_ck;
                        col_ok = col_ok_ck;
                        search_for_outlier = 0;
                    end
                else
                    search_for_outlier = 0;
                end
            end
        end

        function [x_float, Cxx, s02, v_hat, y0, b, A, col_ok, Q, obs_track, amb_prn_track, iQ]  = remSolitaryObs (this, ...
                y0, b, A, col_ok, Q, ...
                obs_track, amb_prn_track, min_contiguous_obs)

            amb_num = numel(amb_prn_track);

            idx_pr = (obs_track(:,3) == -1); %#ok<NASGU>
            idx_ph = (obs_track(:,3) == 1);

            idx_out = [];
            for a = 1 : amb_num
                % find obs of an arc
                id_obs_ok = find(full(idx_ph & (A(:,3 + a) < 0)));
                if ~isempty(id_obs_ok)
                    % find the epochs of these obs
                    id_epoch_ok = zeros(max(obs_track(:,1)),1);
                    id_epoch_ok(obs_track(id_obs_ok)) = id_obs_ok;
                    [lim] = getOutliers(id_epoch_ok);
                    lim_ko = find(lim(:,2)-lim(:,1) + 1 < min_contiguous_obs);

                    for i = 1 : numel(lim_ko)
                        idx_out = [idx_out; id_epoch_ok(lim(lim_ko(i), 1) : lim(lim_ko(i), 2))]; %#ok<AGROW>
                    end
                end
            end
            idx_out = sort(idx_out);

            y0(idx_out) = [];
            A(idx_out,:) = [];
            Q(idx_out,:) = [];
            Q(:,idx_out) = [];
            b(idx_out) = [];
            obs_track(idx_out,:) = [];

            n_col = size(A, 2);
            [A, y0, b, Q, obs_track, amb_prn_track] = this.remShortArc(A, y0, b, Q, obs_track, amb_prn_track, this.state.getMinArc());
            if size(A, 2) < n_col
                [col_ok] = this.getBlockProperties(A, obs_track, 1);
            end
            iQ = cholinv(Q);
            [x_float, Cxx, s02, v_hat] = this.solveLS(y0, b, A, col_ok, iQ);
        end

%         function [y0, Q] = preCorrectObsIntAmb(this)
%             A = this.A;
%             Q = this.Q;
%             n_pos =  this.n_pos;
%             y0 = this.y0 - this.b;
%
%             % epoch of a change of pivot
%             pivot_change = find(abs(diff(int8(this.pivot_track)))>0);
%
%             n_amb = size(A, 2) - n_pos * 3;
%
%             % First loop find discontinuities in the observations
%             for a = 1 : n_amb
%                 idx = find(A(:, 3 + a) < 0);
%
%                 if numel(idx > 3)
%
%                     lambda_obs = abs(this.A(idx(1), a + this.n_pos * 3));
%
%                     tmp = [0; diff(y0(idx)) / lambda_obs];
%                     tmp = round(cumsum(tmp - movmedian(tmp, 3)));
%
%                     % discontinuities due to a change of pivot are wanted
%                     for i = 1 : numel(pivot_change)
%                         [~, jmp] = intersect(idx, find(this.obs_track(:,1) == pivot_change(i)));
%                         jmp = jmp + sum(this.empty_epoch < pivot_change(i));
%                         if (jmp < length(tmp)-1)
%                             tmp(jmp+1:end) = tmp(jmp+1:end) - (tmp(jmp+1) - tmp(jmp));
%                         end
%                     end
%                     tmp = movmedian(tmp,3);
%
%                     % remove the integer discontinuities in the observations
%                     y0(idx) = y0(idx) - tmp * lambda_obs;
%                 end
%             end
%
%             [~, ~, ~, v_hat] = this.improveFloatSolution(y0 + this.b, this.b, this.A, this.col_ok, this.Q, [], this.obs_track);
%
%             for a = 1 : n_amb
%                 idx = find(A(:, 3 + a) < 0);
%                 if numel(idx > 3)
%
%                     % lambda_obs = abs(this.A(idx(1), a + this.n_pos * 3));
%                     % tmp = round(v_hat(idx) / lambda_obs) * lambda_obs;
%                     % y0(idx) = y0(idx) - tmp;
%
%                     % clean obs jmp
%                     res = [0; diff(y0(idx)) - movmedian(diff(y0(idx)), 3)];
%                     id_ko = abs(res) > 0.04;
%                     id_ko = idx(id_ko(1:end-1) & id_ko(2:end));
%                     Q(id_ko + size(Q,1) * (id_ko - 1)) = Q(id_ko + size(Q,1) * (id_ko - 1)) * 4;
%                     %Q(id_ko + size(Q,1) * (id_ko - 1)) = v_hat(id_ko).^2;
%                 end
%             end
%
%             y0 = y0 + this.b;
%         end

        function [y0, Q] = preCorrectObsIntAmb(this, y0, b, A, Q, n_pos, obs_track, pivot_change)
            if nargin == 1
                y0 = this.y0 - this.b;
                b = this.b;
                A = this.A;
                Q = this.Q;
                n_pos =  this.n_pos;
                obs_track = this.obs_track;
                empty_epoch = this.empty_epoch;
            else
                empty_epoch = [];
                y0 = y0 - b;
            end

            n_amb = size(A, 2) - n_pos * 3;

            % First loop find discontinuities in the observations
            for a = 1 : n_amb
                idx = find(A(:, 3 + a) < 0);

                if numel(idx > 3)

                    lambda_obs = abs(A(idx(1), a + n_pos * 3));

                    tmp = [0; diff(y0(idx))];
                    tmp = round(cumsum(tmp - movmedian(tmp, 3)) / lambda_obs) * lambda_obs;
                    tmp = movmedian(tmp,3);
%                     tmp1 = [0; diff(y0(idx) - tmp)];
%                     tmp1(abs(tmp1) < 0.4 * lambda_obs) = 0;
%                     tmp1 = round(cumsum(tmp1) / lambda_obs) * lambda_obs;
%                     tmp = tmp + tmp1;

                    % discontinuities due to a change of pivot are wanted
                    for i = 1 : numel(pivot_change)
                        [~, jmp] = intersect(idx, find(obs_track(:,1) == pivot_change(i)));
                        jmp = jmp + sum(empty_epoch < pivot_change(i));
                        if (jmp < length(tmp)-1)
                            tmp(jmp+1:end) = tmp(jmp+1:end) - (tmp(jmp+1) - tmp(jmp));
                        end
                    end

                    % remove the integer discontinuities in the observations
                    y0(idx) = y0(idx) - tmp * lambda_obs;
                end
            end

            for a = 1 : n_amb
                idx = find(A(:, 3 + a) < 0);
                if numel(idx > 3)

                    % lambda_obs = abs(this.A(idx(1), a + this.n_pos * 3));
                    % tmp = round(v_hat(idx) / lambda_obs) * lambda_obs;
                    % y0(idx) = y0(idx) - tmp;

                    % clean obs jmp
                    res = [0; diff(y0(idx)) - movmedian(diff(y0(idx)), 3)];
                    id_ko = abs(res) > 0.04;
                    id_ko = idx(id_ko(1:end-1) & id_ko(2:end));
                    Q(id_ko + size(Q,1) * (id_ko - 1)) = Q(id_ko + size(Q,1) * (id_ko - 1)) * 4;
                    %Q(id_ko + size(Q,1) * (id_ko - 1)) = v_hat(id_ko).^2;
                end
            end

            y0 = y0 + b;
        end

        function [y0, is_new] = postCorrectIntAmb (this, ...
                    y0, phase_res, id_track, ...
                    A, amb_prn_track, win_size)
            % try to correct integer ambiguities by using a moving median on the LS residuals
            % fixed at n * lambda levels.
            % Note that the solution MUST be stable to perform this action
            %
            % SYNTAX:
            %   [y0] = postCorrectIntAmb(this, y0, phase_res, id_track, A, amb_prn_track)
            %
            % OUTPUT:
            %   y0      array with "corrected" integer ambiguities
            %   is_new  flag -> true when y0 has been changed
            if (nargin < 7)
                %%%% NCU: win_size = this.state.getMinArc() + mod(1 + this.state.getMinArc(),2);
            end

            is_new = false;
            phase_ref = nan(size(id_track));
            for a = 1 : numel(amb_prn_track)
                y = phase_res(:, a, 1);
                id_obs = id_track(~isnan(y),a);
                if ~ isempty(id_obs)
                    lambda_val = abs(A(id_obs, 3 + a)) * this.cs_factor;
                    %if (numel(y(~isnan(y))) > 1.5 * win_size) % it's probably a pivot
                        ref = this.getAmbCorrection(y0(id_obs) ./ lambda_val, y(~isnan(y)) ./ lambda_val, 1) .* lambda_val;
                        phase_ref(~isnan(y),a) = ref;
                        % If the ref correction provides a reduction of the derivative std -> use it
                        if ~(full(sum(abs(ref)) == 0)) && (std(diff(y(~isnan(y))-ref)) <= std(diff(y(~isnan(y)))))
                        %if ~(full(sum(abs(ref)) == 0))
                            figure(999); clf; plot(y(~isnan(y)),'.-'); hold on; plot(ref,'o', 'lineWidth', 2)
                            is_new = true;
                            figure(1000); clf; plot(diff(y0(id_obs)),'.-'); hold on;
                            y0(id_obs) = y0(id_obs) - ref;
                            plot(diff(y0(id_obs)),'.-'); ylim([-1 1]);
                        end
                    %end
                end
            end
        end

        function [y0, A, b, Q, sat_pr_track, sat_ph_track, pivot] =  oneEpochLS (this, ...
                    time_rx, ...
                    pos_r, pos_m, ...
                    pr1_r, pr1_m, pr2_r, pr2_m, ...
                    ph1_r, ph1_m, ph2_r, ph2_m, ...
                    snr_r, snr_m, ...
                    eph, sp3, iono, lambda, frequencies, p_rate, ant_pcv)

            % SYNTAX:
            %   prepare_dd_sys(time_rx, XR0, XM, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M, Eph, sp3, iono, lambda, phase, ant_pcv);
            %
            % INPUT:
            %   time_rx = GPS reception time
            %   XR0     = ROVER approximate position
            %   XM      = MASTER position
            %   pr1_R   = ROVER code observations (L1 carrier)
            %   pr1_M   = MASTER code observations (L1 carrier)
            %   pr2_R   = ROVER code observations (L2 carrier)
            %   pr2_M   = MASTER code observations (L2 carrier)
            %   ph1_R   = ROVER phase observations (L1 carrier)
            %   ph1_M   = MASTER phase observations (L1 carrier)
            %   ph2_R   = ROVER phase observations (L2 carrier)
            %   ph2_M   = MASTER phase observations (L2 carrier)
            %   snr_R   = ROVER-SATELLITE signal-to-noise ratio
            %   snr_M   = MASTER-SATELLITE signal-to-noise ratio
            %   Eph     = satellite ephemeris
            %   sp3     = structure containing precise ephemeris and clock
            %   iono    = ionosphere parameters
            %   lambda  = wavelength matrix (depending on the enabled constellations)
            %   frequencies  = L1 carrier (phase=1), L2 carrier (phase=2)
            %   p_rate = processing rate
            %   ant_pcv = antenna phase center variation
            %
            % DESCRIPTION:
            %   Computation of the receiver position (X,Y,Z).
            %   Relative (double difference) positioning by least squares adjustment
            %   on code and phase observations.

            cutoff = this.state.cut_off;
            snr_threshold = this.state.snr_thr;
            cond_n_threshold = this.state.cond_num_thr;

            y0 = [];
            A  = [];
            b  = [];
            Q  = [];

            %total number of satellite slots (depending on the constellations enabled)
            n_sat_tot = size(pr1_r,1);

            %topocentric coordinate initialization
            azR   = zeros(n_sat_tot,1);
            elR   = zeros(n_sat_tot,1);
            distR = zeros(n_sat_tot,1);
            azM   = zeros(n_sat_tot,1);
            elM   = zeros(n_sat_tot,1);
            distM = zeros(n_sat_tot,1);

            %--------------------------------------------------------------------------------------------
            % SATELLITE SELECTION
            %--------------------------------------------------------------------------------------------

            % Find sat in common, between master and rover
            if (length(frequencies) == 2)
                sat_pr = pr1_r & pr1_m & pr2_r & pr2_m;
                sat_ph = ph1_r & ph1_m & ph2_r & ph2_m;
            else
                if (frequencies == 1)
                    sat_pr = pr1_r & pr1_m;
                    sat_ph = ph1_r & ph1_m;
                else
                    sat_pr = pr2_r & pr2_m;
                    sat_ph = ph2_r & ph2_m;
                end
            end

            % filter satellites with no ephemeris
            if (isempty(sp3))
                eph_avail = serialize(eph(30,:) > 0);
            else
                eph_avail = sp3.avail;
            end
            sat_pr = find(sat_pr & eph_avail);
            sat_ph = find(sat_ph & eph_avail);

            min_nsat_LS = 3 + this.state.cc.getNumSys();

            flag_XR = 2;

            % satellite configuration
            sat_pr_track = int8(zeros(n_sat_tot, 1));
            sat_ph_track = int8(zeros(n_sat_tot, 1));
            pivot = 0;

            if (size(sat_pr,1) >= min_nsat_LS)

                sat_pr_old = sat_pr;

                if (frequencies == 1)
                    [pos_m, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM] ...
                        = init_positioning(time_rx, pr1_m(sat_pr),   snr_m(sat_pr),   eph, sp3, iono, [],  pos_m, [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, frequencies, p_rate,       2, 0, 0, 0);
                    if (sum(sat_pr_M) < min_nsat_LS); return; end
                    [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] ...
                        = init_positioning(time_rx, pr1_r(sat_pr_M), snr_r(sat_pr_M), eph, sp3, iono, [], pos_r, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 1, 0, 0);
                else
                    [pos_m, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM] ...
                        = init_positioning(time_rx, pr2_m(sat_pr),   snr_m(sat_pr),   eph, sp3, iono, [],  pos_m, [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, frequencies, p_rate,      2, 0, 0, 0);
                    if (sum(sat_pr_M) < min_nsat_LS); return; end
                    [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] ...
                        = init_positioning(time_rx, pr2_r(sat_pr_M), snr_r(sat_pr_M), eph, sp3, iono, [], pos_r, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 1, 0, 0);
                end

                %keep only satellites that rover and master have in common
                [sat_pr, iR, iM] = intersect(sat_pr_R, sat_pr_M);
                XS = XS(iR,:);
                sys = sys(iR);
                if (~isempty(err_tropo_R))
                    err_tropo_R = err_tropo_R(iR);
                    err_iono_R  = err_iono_R (iR);
                    err_tropo_M = err_tropo_M(iM);
                    err_iono_M  = err_iono_M (iM);
                end

                % keep only satellites that rover and master have in common both in phase and code
                [sat_pr, i_pr] = intersect(sat_pr, sat_ph);
                XS = XS(i_pr,:);
                sys = sys(i_pr);
                if (~isempty(err_tropo_R))
                    err_tropo_R = err_tropo_R(i_pr);
                    err_iono_R  = err_iono_R (i_pr);
                    err_tropo_M = err_tropo_M(i_pr);
                    err_iono_M  = err_iono_M (i_pr);
                end

                %apply cutoffs also to phase satellites
                sat_removed = setdiff(sat_pr_old, sat_pr);
                sat_ph(ismember(sat_ph,sat_removed)) = [];

                %--------------------------------------------------------------------------------------------
                % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
                %--------------------------------------------------------------------------------------------

                switch this.sol_type
                    case 0 % code and phase
                        sat_pr_track(sat_pr, 1) = 1;
                        sat_ph_track(sat_ph, 1) = 1;
                    case -1 % code only;
                        sat_pr_track(sat_pr, 1) = 1;
                    case 1 % phase only
                        sat_ph_track(sat_ph, 1) = 1;
                end

                %actual pivot
                [null_max_elR, pivot_index] = max(elR(sat_ph));
                pivot = sat_ph(pivot_index);

                %--------------------------------------------------------------------------------------------
                % PHASE CENTER VARIATIONS
                %--------------------------------------------------------------------------------------------

                %compute PCV: phase and code 1
                [~, index_ph]=intersect(sat_pr,sat_ph);

                if (~isempty(ant_pcv) && ant_pcv(2).n_frequency ~= 0) % master
                    index_master = 2;
                    PCO1_M = PCO_correction(ant_pcv(index_master), pos_r, XS, sys, 1);
                    PCV1_M = PCV_correction(ant_pcv(index_master), 90-elM(sat_pr), azM(sat_pr), sys, 1);
                    pr1_m(sat_pr) = pr1_m(sat_pr) - (PCO1_M + PCV1_M);
                    ph1_m(sat_ph)    = ph1_m(sat_ph)    - (PCO1_M(index_ph) + PCV1_M(index_ph))./lambda(sat_ph,1);

                    if (length(frequencies) == 2 || frequencies(1) == 2)
                        PCO2_M = PCO_correction(ant_pcv(index_master), pos_r, XS, sys, 2);
                        PCV2_M = PCV_correction(ant_pcv(index_master), 90-elM(sat_pr), azM(sat_pr), sys, 2);
                        pr2_m(sat_pr) = pr2_m(sat_pr) - (PCO2_M + PCV2_M);
                        ph2_m(sat_ph)    = ph2_m(sat_ph)    - (PCO2_M(index_ph) + PCV2_M(index_ph))./lambda(sat_ph,2);
                    end
                end

                if (~isempty(ant_pcv) && ant_pcv(1).n_frequency ~= 0) % rover
                    index_rover = 1;
                    PCO1_R = PCO_correction(ant_pcv(index_rover), pos_r, XS, sys, 1);
                    PCV1_R = PCV_correction(ant_pcv(index_rover), 90-elR(sat_pr), azR(sat_pr), sys, 1);
                    pr1_r(sat_pr) = pr1_r(sat_pr) - (PCO1_R + PCV1_R);
                    ph1_r(sat_ph)    = ph1_r(sat_ph)    - (PCO1_R(index_ph) + PCV1_R(index_ph))./lambda(sat_ph,1);

                    if (length(frequencies) == 2 || frequencies(1) == 2)
                        PCO1_R = PCO_correction(ant_pcv(index_rover), pos_r, XS, sys, 2);
                        PCV2_R = PCV_correction(ant_pcv(index_rover), 90-elM(sat_pr), azM(sat_pr), sys, 2);
                        pr2_r(sat_pr) = pr2_r(sat_pr) - (PCO1_R + PCV2_R);
                        ph2_r(sat_ph)    = ph2_r(sat_ph)    - (PCO1_R(index_ph) + PCV2_R(index_ph))./lambda(sat_ph,2);
                    end
                end

                %--------------------------------------------------------------------------------------------
                % PREPARE INPUT FOR LEAST SQUARES BATCH
                %--------------------------------------------------------------------------------------------

                %if at least min_nsat_LS satellites are available after the cutoffs, and if the
                % condition number in the least squares does not exceed the threshold
                if (size(sat_ph,1) >= min_nsat_LS && (isempty(cond_num) || cond_num < cond_n_threshold))

                    if (frequencies == 1)
                        [y0, A, b, Q] = this.oneEpochBuild(pos_r, pos_m, XS, pr1_r(sat_ph), ph1_r(sat_ph), snr_r(sat_ph), pr1_m(sat_ph), ph1_m(sat_ph), snr_m(sat_ph), elR(sat_ph), elM(sat_ph), err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda(sat_ph,1));
                    else
                        [y0, A, b, Q] = this.oneEpochBuild(pos_r, pos_m, XS, pr2_r(sat_ph), ph2_r(sat_ph), snr_r(sat_ph), pr2_m(sat_ph), ph2_m(sat_ph), snr_m(sat_ph), elR(sat_ph), elM(sat_ph), err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda(sat_ph,2));
                    end

                else
                    pivot = 0;
                end
            else
                pivot = 0;
            end
        end

        function [y0, A, b, Q] = oneEpochBuild ( this, ...
                   pos_r_approx, pos_m, pos_s, ...
                   pr_r, ph_r, snr_r, ...
                   pr_m, ph_m, snr_m, ...
                   el_r, el_m, ...
                   err_tropo_r, err_iono_r, ...
                   err_tropo_m, err_iono_m, pivot_id, lambda)

            % SYNTAX:
            %   [y0, A, b, Q] = this.oneEpochBuild (XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda);
            %
            % INPUT:
            %   XR_approx   = receiver approximate position (X,Y,Z)
            %   XM          = master station position (X,Y,Z)
            %   XS          = satellite position (X,Y,Z)
            %   pr_R        = receiver code observations
            %   ph_R        = receiver phase observations
            %   pr_M        = master code observations
            %   pr_M        = master phase observations
            %   snr_R       = receiversignal-to-noise ratio
            %   snr_M       = mastersignal-to-noise ratio
            %   elR         = satellite elevation (vector)
            %   elM         = satellite elevation (vector)
            %   err_tropo_R = tropospheric error
            %   err_tropo_M = tropospheric error
            %   err_iono_R  = ionospheric error
            %   err_iono_M  = ionospheric error
            %   pivot_index = index identifying the pivot satellite
            %   lambda      = vector containing GNSS wavelengths for available satellites
            %
            % OUTPUT:
            %   y0 = observation vector
            %   A = design matrix
            %   b = known term vector
            %   Q = observation covariance matrix
            %
            % DESCRIPTION:
            %   Function that prepares the input matrices for the least squares batch solution.

            % variable initialization
            global sigmaq_cod1 sigmaq_ph

            % number of observations
            n = length(pr_r);

            % approximate receiver-satellite distance
            XR_mat = pos_r_approx(:,ones(n,1))';
            XM_mat = pos_m(:,ones(n,1))';
            distR_approx = sqrt(sum((pos_s-XR_mat).^2 ,2));
            distM = sqrt(sum((pos_s-XM_mat).^2 ,2));

            % design matrix (code or phase)
            A = [((pos_r_approx(1) - pos_s(:,1)) ./ distR_approx) - ((pos_r_approx(1) - pos_s(pivot_id,1)) / distR_approx(pivot_id)), ... %column for X coordinate
                ((pos_r_approx(2) - pos_s(:,2)) ./ distR_approx) - ((pos_r_approx(2) - pos_s(pivot_id,2)) / distR_approx(pivot_id)), ... %column for Y coordinate
                ((pos_r_approx(3) - pos_s(:,3)) ./ distR_approx) - ((pos_r_approx(3) - pos_s(pivot_id,3)) / distR_approx(pivot_id))];    %column for Z coordinate

            % known term vector
            b    =     (distR_approx - distM)      - (distR_approx(pivot_id) - distM(pivot_id));       %approximate pseudorange DD
            b    = b + (err_tropo_r - err_tropo_m) - (err_tropo_r(pivot_id)  - err_tropo_m(pivot_id)); %tropospheric error DD

            if (this.sol_type <= 0)
                % known term vector
                b_pr = b + (err_iono_r  - err_iono_m)  - (err_iono_r(pivot_id)   - err_iono_m(pivot_id));  %ionoshperic error DD (code)
                % observation vector
                y0_pr = (pr_r - pr_m) - (pr_r(pivot_id) - pr_m(pivot_id));
                % remove pivot-pivot lines
                b_pr(pivot_id)    = [];
                y0_pr(pivot_id)   = [];
            end

            if (this.sol_type >= 0)
                % known term vector
                b_ph = b - (err_iono_r  - err_iono_m)  + (err_iono_r(pivot_id)   - err_iono_m(pivot_id));  %ionoshperic error DD (phase)
                % observation vector
                y0_ph = lambda.*((ph_r - ph_m) - (ph_r(pivot_id) - ph_m(pivot_id)));
                % remove pivot-pivot lines
                b_ph(pivot_id)    = [];
                y0_ph(pivot_id)   = [];
            end

            % remove pivot-pivot lines
            A(pivot_id, :) = [];

            % observation noise covariance matrix
            Q1 = cofactor_matrix(el_r, el_m, snr_r, snr_m, pivot_id);

            switch this.sol_type
                case 0 % code and phase
                    A = [A; A];
                    b = [b_pr; b_ph];
                    y0 = [y0_pr; y0_ph];
                    n = 2*n - 2;
                    Q = zeros(n);
                    Q(1:n/2,1:n/2) = sigmaq_cod1 * Q1;
                    Q(n/2+1:end,n/2+1:end) = sigmaq_ph * Q1;
                case -1 % code only;
                    b = b_pr;
                    y0 = y0_pr;
                    Q = sigmaq_cod1 * Q1;
                case 1 % phase only
                    b = b_ph;
                    y0 = y0_ph;
                    Q = sigmaq_ph * Q1;
            end
        end

        function id_track = computeIdTrack(this, A, obs_track, n_pos)
            % Get a matrix of n_obs x n_arcs with the index of the observations in y0
            % Where no observations are present it contains zeros
            % SYNTAX: id_track = this.computeIdTrack(A, obs_track, n_pos)
            if (nargin == 1)
                A = this.A;
                obs_track = this.obs_track;
                n_pos =  this.n_pos;
            end
            if ~isempty(A)
                n_arc = size(A, 2) - n_pos * 3;
                n_tot_epoch = obs_track(end,1);

                id_track = spalloc(n_tot_epoch, n_arc, round(n_tot_epoch * 2));
                for a = 1 : n_arc
                    idx = find(A(:, 3 * n_pos + a) < 0);
                    id_track(obs_track(idx, 1), a) = idx; %#ok<*SPRIX>
                end
            else
                id_track = [];
            end
        end

        function [data_track, id_track] = computeDataTrack(this, data, id_track, A, obs_track, n_pos)
            % Get a matrix of n_obs x n_arcs with the values of the observations in data
            % Where no observations are present for a certain arcs it contains zerosÿ
            %
            % INPUT:
            %   data       data to be put in the matrix [n_obs x n_set]
            %   id_track   can be empty and computed automatically contains the id of the obs to fill the matrix
            %   A          required if id_track is empty, default value from this.
            %   obs_track  required if id_track is empty, default value from this.
            %   n_pos          required if id_track is empty, default value from this.
            %
            % SYNTAX:
            %       [data_track, id_track] = this.computeDataTrack(data, id_track, A, obs_track, n_pos)
            % EXAMPLE:
            %       [data_track, id_track] = this.computeDataTrack(v_hat);
            %       [data_track, id_track] = this.computeDataTrack(v_hat, id_track);
            %       [data_track, id_track] = this.computeDataTrack(v_hat, [],  A, obs_track, n_pos);

            if nargin == 1
                data = this.v_hat;
            end
            if (nargin < 3) || isempty(id_track)
                if (nargin < 5)
                    A = this.A;
                    obs_track = this.obs_track;
                    n_pos =  this.n_pos;
                end
                id_track = this.computeIdTrack(A, obs_track, n_pos);
            end
            [n_obs, n_arc] = size(id_track);
            n_set = size(data,2);

            data_track = nan(n_obs, n_arc, n_set);
            for a = 1 : n_arc
                idx = id_track(:, a) ~= 0;
                for s = 1 : n_set
                    data_track(idx, a, s) = data(id_track(idx, a), s); %#ok<*SPRIX>
                end
            end
        end

        function go_block = getCopy(this)
            % get a copy of a Core_Block instance
            % SYNTAX: go_block = this.getCopy();
            go_block = Core_Block(0,0,0);
            go_block.import(this);
        end

        function [out_struct] = struct(this)
            % convert to struct the content of the obj
            % SYNTAX: [out_struct] = this.struct;
            out_struct.log = this.log;
            out_struct.state = this.state;
            out_struct.n_pos = this.n_pos;

            out_struct.n_pos_hr = this.n_pos_hr;
            out_struct.s_rate = this.s_rate;
            out_struct.sol_type = this.sol_type;
            out_struct.n_obs_tot = this.n_obs_tot;
            out_struct.n_epoch = this.n_epoch;
            out_struct.n_tot_epoch = this.n_tot_epoch;
            out_struct.time_diff = this.time_diff;
            out_struct.empty_epoch = this.empty_epoch;
            out_struct.sat_pr_track = this.sat_pr_track;
            out_struct.sat_ph_track = this.sat_ph_track;
            out_struct.pivot_track = this.pivot_track;

            out_struct.obs_track = this.obs_track;
            out_struct.amb_prn_track = this.amb_prn_track;
            out_struct.ref_arc = this.ref_arc;
            out_struct.col_ok = this.col_ok;

            % LS variable
            out_struct.A = this.A;
            out_struct.Q = this.Q;
            out_struct.y0 = this.y0;
            out_struct.b = this.b;

            % Results
            out_struct.pos0 = this.pos0;
            out_struct.pos = this.pos;
            out_struct.pos_cov = this.pos_cov;
            out_struct.is_fixed = this.is_fixed;

            out_struct.x_float = this.x_float;
            out_struct.x_fix = this.x_fix;
            out_struct.amb_fix = this.amb_fix;
            out_struct.amb_fix_full = this.amb_fix_full;
            out_struct.G = this.G;
            out_struct.pos_cov_fix = this.pos_cov_fix;
            out_struct.x_hr = this.x_hr;
            out_struct.Cxx = this.Cxx;
            out_struct.s02 = this.s02;
            out_struct.v_hat = this.v_hat;

            out_struct.phase_res = this.phase_res;
            out_struct.id_track = this.id_track;

            out_struct.cs_factor = this.cs_factor;
        end

        function import(this, obj)
            % import from struct or obj all the properties
            % SYNTAX: this.import(obj);
            this.log = obj.log;
            this.state = obj.state;
            this.n_pos = obj.n_pos;

            this.n_pos_hr = obj.n_pos_hr;
            this.s_rate = obj.s_rate;
            this.sol_type = obj.sol_type;
            this.n_obs_tot = obj.n_obs_tot;
            this.n_epoch = obj.n_epoch;
            this.n_tot_epoch = obj.n_tot_epoch;
            this.time_diff = obj.time_diff;
            this.empty_epoch = obj.empty_epoch;
            this.sat_pr_track = obj.sat_pr_track;
            this.sat_ph_track = obj.sat_ph_track;
            this.pivot_track = obj.pivot_track;

            this.obs_track = obj.obs_track;
            this.amb_prn_track = obj.amb_prn_track;
            this.ref_arc = obj.ref_arc;
            this.col_ok = obj.col_ok;

            % LS variable
            this.A = obj.A;
            this.Q = obj.Q;
            this.y0 = obj.y0;
            this.b = obj.b;

            % Results
            this.pos0 = obj.pos0;
            this.pos = obj.pos;
            this.pos_cov = obj.pos_cov;
            this.is_fixed = obj.is_fixed;

            this.x_float = obj.x_float;
            this.x_fix = obj.x_fix;
            this.amb_fix = obj.amb_fix;
            this.amb_fix_full = obj.amb_fix_full;
            this.G = obj.G;
            this.pos_cov_fix = obj.pos_cov_fix;
            this.x_hr = obj.x_hr;
            this.Cxx = obj.Cxx;
            this.s02 = obj.s02;
            this.v_hat = obj.v_hat;

            this.phase_res = obj.phase_res;
            this.id_track = obj.id_track;

            this.cs_factor = obj.cs_factor;
        end
        
    end

    % ==================================================================================================================================================
    %  STATIC FUNCTIONS used as utilities goBlock
    % ==================================================================================================================================================
    methods (Static, Access = private)

        function [col_ok, ref_arc] = getBestRefArc(y0, b, A, iQ)
            % compute multiple LS solutions changing the reference arc fort the LS and get the best
            % SYNTAX: id_min = refineRefArc(y0, b, A, iQ)
            amb_num = size(A,2) - 3;
            P = []; N = [];
            if (amb_num > 1)
                test = zeros(amb_num, 1);
                for i = 4 : size(A,2)
                    col_ok = setdiff(1:size(A,2), i);
                    [x, Cxx, s02, v_hat, P, N] = Core_Block.solveLS(y0, b, A, col_ok, iQ, P, N);
                    %amb_var = Cxx(4:end,4:end);
                    %amb_var(amb_var < 0) = 1e30;
                    %test(i-3) = sum(diag(amb_var));
                    test(i-3) = s02;
                end
                [~, ref_arc] = sort(test);
                ref_arc = ref_arc(1);
            else
                ref_arc = [];
            end
            col_ok = setdiff(1:size(A,2), ref_arc + 3);
        end

        function [col_ok, ref_arc, bad_blocks] = getBestBlockRefArc(y0, b, A, iQ, ref_arc, bad_arc, bad_blocks, blk_cols)
            % compute multiple LS solutions changing the reference arc fort the LS and get the best
            % SYNTAX: [col_ok, ref_arc] = this.getBestBlockRefArc(y0, b, A, iQ, ref_arc, bad_arc, blk_cols)

            if isempty(bad_blocks)
                bad_blocks = find(ismember(ref_arc, bad_arc));
            end
            P = []; N = [];
            for bb = 1 : numel(bad_blocks)
                bad_block = bad_blocks(bb);
                if ~isempty(bad_block)
                    ref_arc(bad_block) = 0;

                    amb_num = sum(blk_cols(4 : end, bad_block));
                    if (amb_num > 1)
                        test = ones(size(A, 2) - 3, 1) * 1e30;
                        block_arc = setdiff(find(blk_cols(4 : end, bad_block)), bad_arc);
                        for i = 1 : numel(block_arc)
                            col_ok = setdiff(1:size(A,2), [ref_arc + 3; block_arc(i) + 3]);
                            [~, Cxx, ~, ~, P, N] = Core_Block.solveLS(y0, b, A, col_ok, iQ, P, N);
                            amb_var = Cxx(4:end,4:end);
                            amb_var(amb_var < 0) = 1e30;
                            test(block_arc(i)) = sum(diag(amb_var));
                        end
                        [~, tmp] = sort(test);
                        tmp = tmp(1);
                    else
                        tmp = [];
                    end
                    ref_arc(bad_block) = tmp;
                end
            end
            ref_arc = unique(ref_arc);
            col_ok = setdiff(1:size(A,2), ref_arc + 3);
        end

        function [pos, pos_cov] = applyFix(x_float, Cxx, amb_fix, G)
            % Apply a fix for the ambiguities
            % SYNTAX: [pos, pos_cov] = applyFix(x_float, Cxx, amb_fix, G)
            n_amb = numel(amb_fix);
            n_pos = (size(x_float, 1) - 1 - n_amb) / 3;
            if (nargin >= 4)
                if size(G, 1) < (n_pos * 3 + n_amb)
                    G_multipos = eye((n_pos * 3 + n_amb), (n_pos * 3 + n_amb) + 1);
                    G_multipos((n_pos - 1) * 3 + 1 : end, (n_pos - 1) * 3 +1 : end) = G;
                else
                    G_multipos = G;
                end
                x_float = G_multipos * x_float;
                Cxx = full(G_multipos * Cxx * G_multipos');
            end

            cov_pos  = Cxx(1 : 3 * n_pos, 1 : 3 * n_pos);     % position covariance block
            cov_amb  = Cxx(3 * n_pos + 1 : end, 3 * n_pos + 1 : end); % ambiguity covariance block
            cov_cross = Cxx(1 : 3 * n_pos, 3 * n_pos + 1 : end);   % position-ambiguity covariance block

            try
                % stabilize cov_N;
                [U] = chol(cov_amb);
                cov_amb = U'*U;
                pos = reshape(x_float(1 : n_pos * 3) - cov_cross * cholinv(cov_amb) * (x_float(3 * n_pos + 1 : end) - amb_fix(:, 1)), 3, n_pos);
                pos_cov = cov_pos  - cov_cross * cholinv(cov_amb) * cov_cross';
            catch ex
                log = Core.getLogger();
                log.addWarning(sprintf('Phase ambiguities covariance matrix unstable - %s', ex.message));
                pos = reshape(x_float(1 : n_pos * 3) - cov_cross * (cov_amb)^-1 * (x_float(3 * n_pos + 1 : end) - amb_fix(:, 1)), 3, n_pos);
                pos_cov = cov_pos  - cov_cross * (cov_amb)^-1 * cov_cross';
            end
        end

        function [A, y0, b, Q, obs_track, amb_prn_track, rem_col] = remShortArc (A, y0, b, Q, obs_track, amb_prn_track, min_arc)
            % Remove ambiguity unkowns with arcs shorter than given threshold
            % SYNTAX: [A, y0, b, Q, obs_track, amb_prn_track, rem_col] = remShortArc (A, y0, b, Q, obs_track, amb_prn_track, min_arc)

            amb_num = numel(amb_prn_track);
            pos_num = size(A,2) - amb_num;
            rem_col = setdiff(find(sum(A~=0,1) < min_arc), 1 : pos_num);

            if (~isempty(rem_col))
                rem_obs = [];
                for r = 1 : length(rem_col)
                    rem_obs = [rem_obs; find(A(:,rem_col(r))~=0)]; %#ok<AGROW>
                end
                A(rem_obs,:) = [];
                y0(rem_obs) = [];
                b(rem_obs) = [];
                Q(rem_obs,:) = []; Q(:,rem_obs) = [];
                obs_track(rem_obs,:) = [];

                A(:,rem_col) = [];
                amb_prn_track(rem_col - pos_num) = [];
            end
        end

        function [A, y0, b, Q, obs_track, amb_prn_track] = remArcCol (A, y0, b, Q, obs_track, amb_prn_track, rem_amb)
            % remove one arc from the LS system
            % SYNTAX: [A, y0, b, Q, obs_track, amb_prn_track] = remArcCol(A, y0, b, Q, obs_track, amb_prn_track, rem_amb)
            if (~isempty(rem_amb))
                rem_obs = [];
                for r = 1 : length(rem_amb)
                    rem_obs = [rem_obs; find(A(:,rem_amb(r))~=0)]; %#ok<AGROW>
                end
                A(rem_obs,:) = [];
                y0(rem_obs) = [];
                b(rem_obs) = [];
                Q(rem_obs,:) = []; Q(:,rem_obs) = [];
                obs_track(rem_obs,:) = [];

                A(:,rem_amb) = [];
                amb_prn_track(rem_amb-3) = [];
            end
        end

        function [id_out] = findBadObs (phase_res, thr, win_size)
            % find obs far from a median solution (not currently used)
            % SYNTAX: [id_out] = findBadObs(phase_res, thr, win_size)
            id_out = [];
            for a = 1 : size(phase_res, 2)
                y = phase_res(:, a, 1);
                if nargin == 3
                    mm = medfilt_mat(y(~isnan(y)), win_size);
                else
                    mm = 0;
                end

                e_sup = y; e_sup(~isnan(y)) = mm + thr;
                e_inf = y; e_inf(~isnan(y)) = mm - thr;

                id_out = [id_out; (a-1) * size(phase_res, 1) + find(y > e_sup | y < e_inf)]; %#ok<AGROW>
            end
        end

        function [A, amb_prn_track]  = splitArcs (A, obs_track, amb_prn_track)

            amb_num = numel(amb_prn_track);

            idx_ph = (obs_track(:,3) == 1);

            for a = 1 : amb_num
                % find obs of an arc
                id_obs_ok = find(full(idx_ph & (A(:,3 + a) < 0)));
                % find the epochs of these obs
                id_epoch_ok = zeros(max(obs_track(:,1)),1);
                id_epoch_ok(obs_track(id_obs_ok)) = id_obs_ok;
                [lim] = getOutliers(id_epoch_ok);
                for i = 2 : size(lim,1)
                    new_col = zeros(size(A,1),1);
                    new_col(lim(i,1) : lim(i,2)) = A(lim(i,1) : lim(i,2), 3 + a);
                    amb_prn_track = [amb_prn_track; amb_prn_track(a)]; %#ok<AGROW>
                    A(lim(i,1) : lim(i,2), 3 + a) = 0;
                    A = [A new_col]; %#ok<AGROW>
                end
            end
        end

    end

end
