%   CLASS Core Block
% =========================================================================
%
% DESCRIPTION
%   Class to manage goBlock solutions
%
% EXAMPLE
%   go_block = Core_Block();
%
% FOR A LIST OF CONSTANTs and METHODS use doc goGNSS
%
% Note for the future: the class uses the current obs storage of goGPS
% -> switch to objects for rover and master observations is suggested

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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef Core_Block < handle

    properties % Public Access
        logger;
        state;

        % number of position solutions to be estimated
        n_pos = 1;

        % flag to spacify the tipe of solution code only (-1) - code/phase (0) - phase only (1)
        sol_type = int8(1);

        % total number of observations that can be used
        n_obs_tot = 0;

        % number of valid epochs used in goBlock
        n_epoch = 0;
        % max number of valid epochs used in goBlock
        n_tot_epoch = 0;

        % indexes of empty observations
        empty_epoch = [];

        % matrices to keep track of the satellite configuration changes (to fill in the proper ambiguity slots)
        sat_pr_track; % satellite configuration matrix for code observations
        sat_ph_track; % satellite configuration matrix for phase observations
        pivot_track; % satellite configuration matrix for pivot tracking

        obs_track; % matrix to keep track of the obs -> epoch; PRN; flag code/phase

        amb_prn_track; % store the prn of each ambiguity

        % LS variable
        A;  % LS Design matrix
        Q;  % LS Cofactor matrix
        y0; % LS Observations array
        b;  % LS known array

        % Results
        pos;          % estimated position
        pos0;         % a-priori position
        pos_cov;      % estimated covariance matrix of the positions
        is_fixed = false;  % flag is fixed

        x_float;      % estimated parameter s(float solution)
        x_fix;        % estimated parameters (fix solution)
        x_final;      % estimated parameters (fix solution + float)
        Cxx;          % Covariance matrix of the parameters
        sigma02_hat;  % estimated variance
        v_hat;        % residuals of the obbservarions;

        phase_res;    % phase residuals ([n_obs x n_amb x 2]); first slice value, second slice weight
        id_res;       % id in the design matrix of the observations in phase_res
    end

    properties (Constant, Access = private)
        FLAG_CODE_ONLY  = int8(-1);
        FLAG_CODE_PHASE = int8(0);
        FLAG_PHASE_ONLY = int8(1);
    end

    methods (Static)
        function this = Core_Block(n_epoch, n_pr_obs, n_ph_obs)
            % Core object creator initialize the structures needed for the computation:
            % EXAMPLE: go_block = Core_Block(n_epoch, n_pr_obs, n_ph_obs)
            
            this.logger = Logger.getInstance();
            this.state = Go_State.getCurrentSettings();

            % number of position solutions to be estimated
            this.n_pos = 1;

            this.n_epoch = n_epoch;
            this.n_tot_epoch = n_epoch;

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

            this.obs_track = NaN(this.n_obs_tot,3); % epoch; PRN; flag code/phase
            this.empty_epoch = [];

            % init LS variables
            this.y0 = NaN(this.n_obs_tot, 1);
            this.b  = NaN(this.n_obs_tot, 1);
            if (this.state.isModeSA)
                this.A = spalloc(this.n_obs_tot, this.n_pos * 3 + this.n_obs_tot, round(2.5 * this.n_obs_tot));
            else
                this.A = sparse(this.n_obs_tot, this.n_pos * 3);
            end
            this.Q  = sparse(this.n_obs_tot, this.n_obs_tot, 10 * this.n_obs_tot);
        end
    end

    % =========================================================================
    %  PROCESSING FUNCTIONS goBlock
    % =========================================================================

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
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);

            this.logger.addMarkedMessage('Preparing goBlock system');

            this.pos = [];          % estimated parameters
            this.pos0 = pos_r;      % a-priori position
            this.pos_cov = [];      % estimated parameters
            this.is_fixed = 0;      % flag is fixed

            this.x_float = [];      % estimated float parameters
            this.x_fix = [];        % estimated parameters (fix solution)
            this.x_final = [];      % estimated parameters (fix solution + float)
            this.Cxx = [];          % Covariance matrix of the parameters
            this.sigma02_hat = [];  % estimated variance
            this.v_hat = [];        % residuals of the obbservarions;

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
            w_bar.createNewBar('Processing...');

            % init loop
            for t = epoch_ok
                eph_t = rt_find_eph (eph, time_diff(t), n_sat);

                [y0_epo, A_epo, b_epo, Q_epo, this.sat_pr_track(:, t), this.sat_ph_track(:, t), pivot] = this.oneEpochLS (time_diff(t), pos_r, pos_m(:,t), pr1_r(:,t), pr1_m(:,t), pr2_r(:,t), pr2_m(:,t), ph1_r(:,t), ph1_m(:,t), ph2_r(:,t), ph2_m(:,t), snr_r(:,t), snr_m(:,t), eph_t, sp3, iono, lambda, frequencies(1), ant_pcv);

                if (pivot > 0)
                    n_obs = length(y0_epo);

                    idx = epoch_track + (1 : n_obs)';

                    this.y0( idx) = y0_epo;
                    this.b ( idx) =  b_epo;
                    this.A ( idx, 1:3) = A_epo;
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
            
            % Add to the Design matric the columns relative to the ambbiguities
            this.setAmbiguities (lambda);
        end

        function setAmbiguities (this, lambda)
            % Add to the internal Design Matrix the columns related to the phase observations
            % (Integer ambiguities - N)
            %
            % METHODS CALL REQUIREMENTS:
            %  -> prepare
            %
            % SYNTAX:
            %   this.setAmbiguities(lambda)
            %
            % INPUT:
            %   lambda     wavelength matrix (depending on the enabled constellations)
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.setAmbiguities(lambda)

            this.logger.addMarkedMessage('Set up Design Matrix to estimate integer abiguities');

            if (this.state.isModePh())
                sat_avail = any(this.sat_ph_track, 2);
                amb_num = sum(sat_avail);
                amb_prn = find(sat_avail);
                this.amb_prn_track = amb_prn;
                amb_idx = 1 : amb_num;
            end

            if (this.state.isModeSA())
                % set the design matrix to estimate the receiver clock
                this.A(:, 3 * this.n_pos + 1) = 1;

                % clear A avoid problems when launching setAmbiguities() more than once
                this.A = this.A(:, 1 : 3 * this.n_pos + 1);
            else
                % clear A avoid problems when launching setAmbiguities() more than once
                this.A = this.A(:, 1 : 3 * this.n_pos);
                if (this.state.isModePh())
                    % matrix to detect re-initialized ambiguities
                    % every time there os an interruption in the satellite observations -> suppose cs
                    amb_track = [int8(zeros(size(this.sat_ph_track,1), 1)) diff(this.sat_ph_track')'];

                    % resize the design matrix to estimate phase ambiguities
                    this.A = [this.A sparse(this.n_obs_tot, amb_num)];
                    rows = 0;
                    for e = 1 : this.n_epoch

                        % check if a new ambiguity column for A is needed

                        % detect new ambiguities on epoch e
                        amb_prn_new = find(amb_track(:,e) == 1);
                        if (~isempty(amb_prn_new))
                            [~, amb_idx_new] = intersect(amb_prn, amb_prn_new);

                            % add a new amb column for the same satellite only if it already had estimates previously
                            old_amb = find(any(amb_track(:,1:e-1) == -1, 2));
                            [amb_prn_new, idx] = intersect(amb_prn_new, old_amb);
                            amb_idx_new = amb_idx_new(idx);

                            if (~isempty(amb_prn_new))
                                amb_idx(amb_idx_new) = max(amb_idx) + (1 : length(amb_prn_new));
                                amb_num = amb_num + length(amb_prn_new);
                                this.amb_prn_track = [this.amb_prn_track; amb_prn_new];
                                this.A = [this.A sparse(this.n_obs_tot, numel(amb_idx_new))];
                            end
                        end

                        % build new columns
                        pivot_prn = this.pivot_track(e);
                        this.sat_pr_track(pivot_prn, e) = -1;
                        this.sat_ph_track(pivot_prn, e) = -1;

                        amb_prn_avail = find(this.sat_ph_track(:,e) == 1);
                        [~, amb_idx_avail] = intersect(amb_prn, amb_prn_avail);
                        pivot_id = amb_idx(amb_prn == pivot_prn);
                        %sat_ph_idx = sum(this.sat_pr_track(:,e) == 1) + (1 : sum(this.sat_ph_track(:, e) == 1));
                        sat_ph_idx = sum(this.sat_pr_track(:,e) == 1) + (1 : sum(this.sat_ph_track(:, e) == 1));
                        rows = rows(end) + sat_ph_idx;
                        this.A(rows + size(this.A,1) * (3 * this.n_pos + amb_idx(amb_idx_avail) -1)) = - lambda(amb_prn_avail, 1);
                        this.A(rows + size(this.A,1) * (3 * this.n_pos + pivot_id -1)) = lambda(amb_prn_avail,1);
                    end
                    %[this.amb_prn_track, reorder_id] = sort(this.amb_prn_track);
                    %this.A = this.A(:, [(1 : (3 * this.n_pos))'; (3 * this.n_pos) + reorder_id]);
                end
            end
        end

        function [pos, pos_cov] = solveFloat (this, flag_outlier)
            % Compute a first float solution using the internal object properties
            %
            % METHODS CALL REQUIREMENTS:
            %  -> prepare -> setAmbiguities
            %
            % SYNTAX:
            %   [pos, pos_cov] = this.solveFloat(this)
            %
            % OUTPUT:
            %   pos     coordinates of the estimated positions
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.setAmbiguities(lambda)
            %   go_block.solveFloat()

            this.logger.addMarkedMessage('Compute a float solution');

            if nargin < 2
                flag_outlier = this.state.flag_outlier;
            end

            amb_num = numel(this.amb_prn_track);

            if (this.state.isModePh())
                [this.A, this.y0, this.b, this.Q, this.obs_track, ~, this.amb_prn_track] = LS_short_arc_removal(this.A, this.y0, this.b, this.Q, this.obs_track, amb_num, this.amb_prn_track, this.state.getMinArc());
            end

            this.logger.addMessage('       - first estimation');
            % computing a first solution with float ambiguities
            [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat] = fast_least_squares_solver(this.y0, this.b, this.A(:, 1 : end - 1), this.Q);

            this.logger.addMessage('       - improve solution by outlier underweight');
            % Improve solution by iterative increase of bad observations variance
            if (flag_outlier)
                [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat, this.y0,  this.b, this.A, this.Q, this.obs_track, this.amb_prn_track] = this.cleanFloatSolution(this.y0, this.b, this.A, this.Q, this.v_hat, this.obs_track, this.amb_prn_track, 10);
                this.n_epoch = numel(unique(this.obs_track));
            end
            [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat] = this.improveFloatSolution(this.y0, this.b, this.A, this.Q, this.v_hat, this.obs_track);

            % If the system is unstable remove the arcs that are making it so
            bad_arc = find(diag(abs(this.Cxx)) > 10);
            if ~isempty(bad_arc)
                this.logger.addWarning('System found unstable, removing bad arcs');
                [this.A, this.y0, this.b, this.Q, this.obs_track, ~, this.amb_prn_track] = this.remAmbCol(this.A, this.y0, this.b, this.Q, this.obs_track, amb_num, this.amb_prn_track, bad_arc);
                [this.A, this.y0, this.b, this.Q, this.obs_track, ~, this.amb_prn_track] = LS_short_arc_removal(this.A, this.y0, this.b, this.Q, this.obs_track, amb_num, this.amb_prn_track, this.state.getMinArc());
                [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat] = fast_least_squares_solver(this.y0, this.b, this.A(:, 1 : end - 1), this.Q);
                [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat] = this.improveFloatSolution(this.y0, this.b, this.A, this.Q, this.v_hat, this.obs_track);
            end

            % Compute phase residuals
            [this.phase_res, this.id_res] = this.computePhRes();

            % Try to correct integer ambiguities (missed cycle slips
            this.logger.addMessage('       - try to fix previously undetected cycle slips');
            this.y0 = this.correctIntAmbiguities(this.y0, this.phase_res, this.id_res, this.A, this.amb_prn_track);

            % Improve solution by iterative increase of bad observations variance
            this.logger.addMessage('       - recompute the improved solution');

            [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat, Q] = this.improveFloatSolution(this.y0, this.b, this.A, this.Q, this.v_hat, this.obs_track); %#ok<PROPLC,*PROP>

            if (flag_outlier)
                % Delete bad observations and restore variances
                this.logger.addMessage('       - reject outliers');
                [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat, this.y0,  this.b, this.A, this.Q, this.obs_track, this.amb_prn_track] = this.cleanFloatSolution(this.y0, this.b, this.A, this.Q, this.v_hat, this.obs_track, this.amb_prn_track);
                [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat, this.y0,  this.b, this.A, this.Q, this.obs_track, this.amb_prn_track] = this.remSolitaryObs(this.y0, this.b, this.A, this.Q, this.obs_track, this.amb_prn_track, round(this.state.getMinArc()/2));
                this.n_epoch = numel(unique(this.obs_track));

                % If the system is unstable remove the arcs that are making it so
                bad_arc = find(diag(abs(this.Cxx)) > 10);
                if ~isempty(bad_arc)
                    this.logger.addWarning('System found unstable, removing bad arcs');
                    [this.A, this.y0, this.b, this.Q, this.obs_track, ~, this.amb_prn_track] = this.remAmbCol(this.A, this.y0, this.b, this.Q, this.obs_track, amb_num, this.amb_prn_track, bad_arc);
                    [this.A, this.y0, this.b, this.Q, this.obs_track, ~, this.amb_prn_track] = LS_short_arc_removal(this.A, this.y0, this.b, this.Q, this.obs_track, amb_num, this.amb_prn_track, this.state.getMinArc());
                    [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat] = fast_least_squares_solver(this.y0, this.b, this.A(:, 1 : end - 1), this.Q);
                    [this.x_float, this.Cxx, this.sigma02_hat, this.v_hat] = this.improveFloatSolution(this.y0, this.b, this.A, this.Q, this.v_hat, this.obs_track);
                end
            else
                this.Q = Q;
            end

            % Compute phase residuals
            [this.phase_res, this.id_res] = this.computePhRes();

            % show residuals
            %close all; this.plotPhRes();

            % extract estimated position
            this.logger.addMarkedMessage('Float solution computed, rover positions corrections:');
            delta_pos = reshape(this.x_float(1:this.n_pos * 3), 3, this.n_pos);
            pos = repmat(this.pos0(:), 1, this.n_pos) + delta_pos;
            this.logger.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', [mean(this.getENU(pos), 2) mean(this.getDeltaENU(pos), 2)]'));
            pos_cov = full(this.Cxx(1:this.n_pos * 3, 1:this.n_pos * 3));
            this.is_fixed = 0;
            this.pos = pos;
            this.pos_cov = pos_cov;
        end

        function [pos, pos_cov, estim_amb, amb_cov, estim_amb_float, ref_arc] = solveFix (this)
            % Compute a fixed solution using LAMBDA, and the the internal object properties
            %
            % METHODS CALL REQUIREMENTS:
            %   prepare -> setAmbiguities -> solveFloat
            %
            % SYNTAX:
            %   [pos, pos_cov, estim_amb, amb_cov, estim_amb_float, ref_arc] = this.solveFix()
            %
            % INPUT:
            %   internal structure
            %
            % OUTPUT:
            %   pos     coordinates of the estimated positions
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.setAmbiguities(lambda)
            %   go_block.solveFloat()
            %   go_block.solveFix()
            %
            % CONCRETE IMPLEMENTATION IN:
            %   solveFixPar
            %
            
            this.logger.addMarkedMessage('Compute ambiguity fix through LAMBDA');

            [delta_pos, pos_cov, is_fixed, estim_amb, amb_cov, estim_amb_float, ref_arc] = this.solveFixPar (this.x_float, this.Cxx, size(this.A, 2) - 3 * this.n_pos - 1);
            this.is_fixed = is_fixed;

            pos = this.pos;
            if (is_fixed)
                % extract estimated position
                this.logger.addMarkedMessage('Fixed solution computed, rover positions corrections:');
                pos = repmat(this.pos0(:), 1, this.n_pos) + repmat(delta_pos(:), 1, this.n_pos);
                this.logger.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', [mean(this.getENU(pos), 2) mean(this.getDeltaENU(pos), 2)]'));
                this.pos = pos;
                this.pos_cov = pos_cov;
                this.x_fix = [delta_pos; estim_amb_float];
                this.Cxx((3 * this.n_pos + 2) : end, (3 * this.n_pos + 2) : end) = amb_cov;
            end

        end

        function [pos, pos_cov, sigma02_hat, v_hat] = solve(this)
            % Solve Float -> try to extract the most stable float solution
            % Compute a first float solution using the internal object properties
            %
            % METHODS CALL REQUIREMENTS:
            %  -> prepare -> setAmbiguities
            %
            % SYNTAX:
            %   [pos, pos_cov] = this.solveFloat(this)
            %
            % OUTPUT:
            %   pos     coordinates of the estimated positions
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.setAmbiguities(lambda)
            %   go_block.solveFloat()

            this.solveFloat();

            if (this.state.flag_iar)
                % Solve Fix -> get a valid estimation of the integer ambiguities
                [~, ~, estim_amb, amb_cov, estim_amb_float, ref_arc] = this.solveFix();
                
                % Use fixed ambiguities as "suggestions" for a float solution
                if (this.is_fixed)

                    this.logger.addMarkedMessage('Using LAMBDA fix suggestions and recompute float:');
                    num_amb = numel(estim_amb);
                    num_obs = size(this.A, 1);
                    num_est = this.n_pos;

                    % Build a Design matrix with LAMBDA fix suggestions
                    idx = [1 : ref_arc - 1, ref_arc + 1 : num_amb + 1];
                    A_fix = [this.A; sparse(num_amb, size(this.A, 2))];
                    y0_fix = [this.y0; round(estim_amb_float(idx))];
                    b_fix = [this.b; zeros(num_amb, 1)];
                    Q_fix = [[this.Q sparse(size(this.Q, 1), num_amb)]; sparse(num_amb, size(this.Q, 1) + num_amb)];
                    for i = 1 : numel(estim_amb)
                        A_fix(num_obs + i, 3 * num_est + idx(i)) = 1;
                    end
                    
                    Q_fix(num_obs + (1:numel(estim_amb)), num_obs +  (1:numel(estim_amb))) = amb_cov;

                    [x_float, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0_fix, b_fix, A_fix(:, 1 : end - 1), Q_fix);
                    % [phase_res, id_res] = this.computePhRes(v_hat, A_fix, Q_fix, this.obs_track, this.amb_prn_track);
                    % close all; this.plotPhRes(this.phase_res, this.id_res, this.A, this.amb_prn_track);

                    this.sigma02_hat = sigma02_hat;
                    this.v_hat = v_hat;

                    this.logger.addMarkedMessage('Final solution computed, rover positions corrections:');
                    delta_pos = reshape(x_float(1:this.n_pos * 3), 3, this.n_pos);
                    pos = repmat(this.pos0(:), 1, this.n_pos) + delta_pos;
                    this.logger.addMessage(sprintf('       East      %12.4f   %+8.4f m\n       North     %12.4f   %+8.4f m\n       Up        %12.4f   %+8.4f m\n', [mean(this.getENU(pos), 2) mean(this.getDeltaENU(pos), 2)]'));
                    pos_cov = full(Cxx(1:this.n_pos * 3, 1:this.n_pos * 3));
                    %%
                    this.is_fixed = 2;
                    this.pos = pos;
                    this.pos_cov = pos_cov;
                    this.x_final = x_float;
                end
            end
            pos = this.pos;
        end

        function [delta_pos, pos_cov, is_fixed, estim_amb, amb_cov, estim_amb_float, ref_arc] = solveFixPar (this, x_float, Cxx, amb_num)
            % Compute a fixed solution using LAMBDA, and the float solution stored in the object structure
            %
            % SYNTAX:
            %   pos = this.solveFloat(this)
            %
            % OUTPUT:
            %   pos     coordinates of the estimated positions
            %
            % EXAMPLE:
            %   go_block = Core_Block(numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            %   go_block.prepare(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
            %   go_block.setAmbiguities(lambda)
            %   go.solveFloat()

            % Compute a first float solution
            if (nargin < 3) || (isempty(x_float))
                this.solveFloat();
                x_float = this.x_float;
                Cxx = this.Cxx;
                amb_num = size(this.A, 2) - 3 * this.n_pos - 1;
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
                    this.logger.addWarning(sprintf('Phase ambiguities covariance matrix unstable - %s', ex.message));
                end

                % integer phase ambiguity solving by LAMBDA
                [delta_pos, estim_amb, amb_cov, pos_cov] = lambdafix(x(1:3), x(4:end), cov_X, cov_N, cov_XN);

                if (estim_amb ~= x(4:end))
                    is_fixed = 1;
                end

                if (~is_fixed)
                    throw(MException('VerifyOutput:Systemunstable', 'LAMBDA did not return a fixed solution.'));
                end

                estim_amb_float = zeros(amb_num, 1);
                estim_amb_float([1 : ref_arc - 1, ref_arc + 1 : end]) = -(estim_amb - x_float(ref_arc + 3));
                estim_amb_float(ref_arc) = x_float(3 + ref_arc);
            catch ex
                this.logger.addWarning(sprintf('It was not possible to estimate integer ambiguities: a float solution will be output.\n%s',ex.message));
                pos_cov = cov_X;
                amb_cov = cov_N;
                estim_amb_float = x_float(4 : end);
                estim_amb = estim_amb_float;
                delta_pos = x_float(1 : 3);
                ref_arc = 0;
            end
        end

        function [A, y0, b, Q, obs_track, amb_num, amb_prn_track] = remObs (this, A, y0, b, Q, obs_track, amb_num, amb_prn_track, rem_obs)
            A(rem_obs,:) = [];
            y0(rem_obs) = [];
            b(rem_obs) = [];
            Q(rem_obs,:) = []; Q(:,rem_obs) = [];
            obs_track(rem_obs,:) = [];
            [A, y0, b, Q, obs_track, amb_num, amb_prn_track] = LS_short_arc_removal(A, y0, b, Q, obs_track, amb_num, amb_prn_track, this.state.getMinArc());
        end

        function plotPhRes (this, phase_res, id_res, A, amb_prn_track)
            if nargin == 1
                phase_res = this.phase_res;
                id_res = this.id_res;
                A = this.A;
                amb_prn_track = this.amb_prn_track;
            end

            x = (1 : size(phase_res,1));
            for a = 1 : numel(amb_prn_track)
                h = figure(amb_prn_track(a));
                h.Name = sprintf('Sat: %d', amb_prn_track(a));
                h.NumberTitle = 'off';
                dockAllFigures();
                y = phase_res(:, a, 1);
                e = 3*phase_res(:, a, 2);
                plot(x, y,'.-'); hold on;
                hline = findobj(h, 'type', 'line');

                patchColor = min(hline(1).Color + 0.3, 1);
                plot(x, e, x, -e, 'color', patchColor);
                patch([x(~isnan(e)) fliplr(x(~isnan(e)))], [e(~isnan(e)); flipud(-e(~isnan(e)))], 1, ...
                    'facecolor',patchColor, ...
                    'edgecolor','none', ...
                    'facealpha', 0.1);
                lambda_val = abs(A(id_res(~isnan(y),a), 3+a));
                y(~isnan(y)) = movmedian(round(y(~isnan(y)) ./ lambda_val) .* lambda_val, this.state.getMinArc());
                plot(x, y,':k', 'LineWidth', 2); hold on;
            end
        end
    end

    % =========================================================================
    %  GETTER FUNCTIONS goBlock
    % =========================================================================
    methods % Public Access
        function [pos, pos_cov] = getPos(this)
            pos = this.pos;
            pos_cov = this.pos_cov;
        end

        function pos = getFloatPos(this)
            pos = repmat(this.pos0(:), 1, this.n_pos) + reshape(this.x_float(1:this.n_pos * 3), 3, this.n_pos);
        end

        function pos = getFixPos(this)
            if isempty(this.x_fix)
                pos = this.getFloatPos();
            else
                pos = repmat(this.pos0(:), 1, this.n_pos) + reshape(this.x_fix(1:this.n_pos * 3), 3, this.n_pos);
            end
        end

        function pos = getFinalPos(this)
            if isempty(this.x_final)
                pos = this.getFloatPos();
            else
                pos = repmat(this.pos0(:), 1, this.n_pos) + reshape(this.x_final(1:this.n_pos * 3), 3, this.n_pos);
            end
        end

        function ph_res = getPhRes(this)
            ph_res = nan(this.n_tot_epoch, this.state.cc.getNumSat());
            for i = 1 : numel(this.amb_prn_track)
                ph_res(~isnan(this.phase_res(:, i, 1)), this.amb_prn_track(i)) = this.phase_res(~isnan(this.phase_res(:, i, 1)), i, 1);
            end
        end

        function [pos_KAL, Xhat_t_t_OUT, Cee_OUT, pivot_OUT, nsat, fixed_amb] = getLegacyOutput(this)
            % [pos_KAL, Xhat_t_t_OUT, Cee_OUT, pivot_OUT, nsat, fixed_amb] = go_block.getLegacyOutput();
            pos_KAL = this.getPos();
            [Xhat_t_t_OUT, Cee_OUT] = this.getPos();
            pivot_OUT = this.pivot_track;
            nsat = this.state.cc.getNumSat();
            fixed_amb = this.is_fixed;
        end

        function [A, amb_prn_track]  = splitArcs (this, A, obs_track, amb_prn_track)

            amb_num = numel(amb_prn_track);

            idx_ph = (obs_track(:,3) == 1);

            for a = 1 : amb_num
                % find obs of an arc
                id_obs_ok = find(full(idx_ph & (A(:,3 * this.n_pos + a) < 0)));
                % find the epochs of these obs
                id_epoch_ok = zeros(max(obs_track(:,1)),1);
                id_epoch_ok(obs_track(id_obs_ok)) = id_obs_ok;
                [lim] = getOutliers(id_epoch_ok);
                for i = 2 : size(lim,1)
                    new_col = zeros(size(A,1),1);
                    new_col(lim(i,1) : lim(i,2)) = A(lim(i,1) : lim(i,2), 3 * this.n_pos + a);
                    amb_prn_track = [amb_prn_track; amb_prn_track(a)];
                    A(lim(i,1) : lim(i,2), 3 * this.n_pos + a) = 0;
                    A = [A new_col];
                end
            end
        end

        function [delta_enu] = getENU(this, pos)
            if nargin == 1
                pos = this.pos;
            end
            %coordinate transformation (UTM)
            [~, ~, up] = cart2geod(pos(1, :), pos(2, :), pos(3, :));
            [east_utm, north_utm] = cart2plan(pos(1, :), pos(2, :), pos(3, :));

            delta_enu = [(east_utm(:))'; (north_utm)'; (up(:))'];
        end

        function [delta_enu] = getDeltaENU(this, pos)
            if nargin == 1
                pos = this.pos;
            end
            %coordinate transformation (UTM)
            [~, ~, up0] = cart2geod(this.pos0(1, :), this.pos0(2, :), this.pos0(3, :));
            [east_utm0, north_utm0] = cart2plan(this.pos0(1, :), this.pos0(2, :), this.pos0(3, :));

            [~, ~, up] = cart2geod(pos(1, :), pos(2, :), pos(3, :));
            [east_utm, north_utm] = cart2plan(pos(1, :), pos(2, :), pos(3, :));

            delta_enu = [(east_utm0(:) - east_utm(:))'; (north_utm0(:) - north_utm)'; (up0(:) - up(:))'];
        end
    end

    methods (Static) % Public Access

        function go_block = go (time_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, sp3, iono, lambda, ant_pcv)
            go_block = Core_Block (numel(time_diff), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
            go_block.prepare (time_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, sp3, iono, lambda, ant_pcv);
            go_block.setAmbiguities (lambda);
            go_block.solve ();
        end
    end

    methods (Access = private)

        function [x_float, Cxx, sigma02_hat, v_hat, Q] = improveFloatSolution(this, ...
                    y0, b, A, Q, ...
                    v_hat, obs_track)
            % Stabilize solution by increasing bad observations variances
            %
            % SYNTAX:
            %   [x_float, Cxx, sigma02_hat, v_hat, Q] = this.improveFloatSolution(this, y0, b, A, Q, v_hat, obs_track)

            x_float = [];
            search_for_outlier = 1;

            out_ph_old = false(size(v_hat));
            out_pr_old = false(size(v_hat));
            n_out_old = 0;
            idx_pr = find(obs_track(:,3) == -1);
            idx_ph = find(obs_track(:,3) == 1);
            while (search_for_outlier == 1)
                % never remove more than 0.5% of data at time
                out_pr = abs(v_hat(idx_pr)) > max(this.state.getMaxCodeErrThr() , perc(abs(v_hat(~out_pr_old)), 0.995));
                if isempty(out_pr)
                    out_pr = out_pr_old;
                end
                out_ph = abs(v_hat(idx_ph)) > max(this.state.getMaxPhaseErrThr() , perc(abs(v_hat(~out_ph_old)), 0.995));
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
                    Qout(idx_out + size(Q,1) * (idx_out - 1)) = v_hat(idx_out).^2; % Bad observations have now their empirical error

                    [x_float, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0, b, A(:, 1 : end - 1), Q);
                else
                    if isempty(x_float)
                        [x_float, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0, b, A(:, 1 : end - 1), Q);
                    end
                    search_for_outlier = 0;
                end
            end
        end

        function [x_float, Cxx, sigma02_hat, v_hat, y0, b, A, Q, obs_track, amb_prn_track]  = cleanFloatSolution (this, ...
                    y0, b, A, Q, ...
                    v_hat, obs_track, amb_prn_track, thr)
            % Stabilize solution by removing bad observations variances
            % [x_float, Cxx, sigma02_hat, v_hat, y0, b, A, Q, obs_track, amb_prn_track] = cleanFloatSolution(this, y0, b, A, Q, v_hat, obs_track, amb_prn_track)

            x_float = [];
            amb_num = numel(amb_prn_track);
            search_for_outlier = 1;

            if nargin < 9
                thr = 3;
            end

            out_ph_old = false(size(v_hat));
            out_pr_old = false(size(v_hat));
            while (search_for_outlier == 1)
                idx_pr = find(obs_track(:,3) == -1);
                idx_ph = find(obs_track(:,3) == 1);
                % never remove more than 0.5% of data at time
                out_pr = abs(v_hat(idx_pr)) > max(max(this.state.getMaxCodeErrThr() , thr * sqrt(Q(idx_pr + size(Q,1) * (idx_pr - 1)))), perc(abs(v_hat), 0.995));
                if isempty(out_pr)
                    out_pr = out_pr_old;
                end
                out_ph = abs(v_hat(idx_ph)) > max(max(this.state.getMaxPhaseErrThr() , thr * sqrt(Q(idx_ph + size(Q,1) * (idx_ph - 1)))), perc(abs(v_hat), 0.995));
                if isempty(out_ph)
                    out_ph = out_ph_old;
                end
                idx_out_pr = idx_pr(out_pr);
                idx_out_ph = idx_ph(out_ph);
                idx_out = [idx_out_pr idx_out_ph];
                if (~isempty(idx_out))
                    y0(idx_out) = [];
                    A(idx_out,:) = [];
                    Q(idx_out,:) = [];
                    Q(:,idx_out) = [];
                    b(idx_out) = [];
                    obs_track(idx_out,:) = [];

                    [A, y0, b, Q, obs_track, amb_num, amb_prn_track] = LS_short_arc_removal(A, y0, b, Q, obs_track, amb_num, amb_prn_track, this.state.getMinArc());

                    [x_float, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0, b, A(:, 1 : end - 1), Q);
                else
                    if isempty(x_float)
                        [x_float, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0, b, A(:, 1 : end - 1), Q);
                    end
                    search_for_outlier = 0;
                end
            end
        end

        function [x_float, Cxx, sigma02_hat, v_hat, y0, b, A, Q, obs_track, amb_prn_track]  = remSolitaryObs (this, ...
                y0, b, A, Q, ...
                obs_track, amb_prn_track, min_contiguous_obs)

            amb_num = numel(amb_prn_track);

            idx_pr = (obs_track(:,3) == -1);
            idx_ph = (obs_track(:,3) == 1);

            idx_out = [];
            for a = 1 : amb_num
                % find obs of an arc
                id_obs_ok = find(full(idx_ph & (A(:,3 * this.n_pos + a) < 0)));
                % find the epochs of these obs
                id_epoch_ok = zeros(max(obs_track(:,1)),1);
                id_epoch_ok(obs_track(id_obs_ok)) = id_obs_ok;
                [lim] = getOutliers(id_epoch_ok);
                lim_ko = find(lim(:,2)-lim(:,1) + 1 < min_contiguous_obs);

                for i = 1 : numel(lim_ko)
                    idx_out = [idx_out; id_epoch_ok(lim(lim_ko(i), 1) : lim(lim_ko(i), 2))];
                end
            end
            idx_out = sort(idx_out);

            y0(idx_out) = [];
            A(idx_out,:) = [];
            Q(idx_out,:) = [];
            Q(:,idx_out) = [];
            b(idx_out) = [];
            obs_track(idx_out,:) = [];

            [A, y0, b, Q, obs_track, amb_num, amb_prn_track] = LS_short_arc_removal(A, y0, b, Q, obs_track, amb_num, amb_prn_track, this.state.getMinArc());
            [x_float, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0, b, A(:, 1 : end - 1), Q);
        end

        function [y0] = correctIntAmbiguities (this, ...
                    y0, phase_res, id_res, ...
                    A, amb_prn_track)
            % try to correct integer ambiguities by using a moving median on the LS residuals
            % fixed at n * lambda levels.
            % Note that the solution MUST be stable to perform this action
            %
            % SYNTAX:
            %   [y0] = correctIntAmbiguities(this, y0, phase_res, id_res, A, amb_prn_track)
            %
            % OUTPUT:
            %   y0      array with "corrected" integer ambiguities
            for a = 1 : numel(amb_prn_track)
                y = phase_res(:, a, 1);
                id_obs = id_res(~isnan(y),a);
                lambda_val = abs(A(id_obs, 3 * this.n_pos + a));
                y0(id_obs) = y0(id_obs) - movmedian(round(y(~isnan(y)) ./ lambda_val) .* lambda_val, this.state.getMinArc());
            end
        end

        function [y0, A, b, Q, sat_pr_track, sat_ph_track, pivot] =  oneEpochLS (this, ...
                    time_rx, ...
                    pos_r, pos_m, ...
                    pr1_r, pr1_m, pr2_r, pr2_m, ...
                    ph1_r, ph1_m, ph2_r, ph2_m, ...
                    snr_r, snr_m, ...
                    eph, sp3, iono, lambda, frequencies, ant_pcv)

            % SYNTAX:
            %   prepare_dd_sys(time_rx, XR0, XM, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M, Eph, sp3, iono, lambda, phase, ant_pcv);
            %
            % INPUT:
            %   time_rx = GPS reception time
            %   XR0   = ROVER approximate position
            %   XM    = MASTER position
            %   pr1_R = ROVER code observations (L1 carrier)
            %   pr1_M = MASTER code observations (L1 carrier)
            %   pr2_R = ROVER code observations (L2 carrier)
            %   pr2_M = MASTER code observations (L2 carrier)
            %   ph1_R = ROVER phase observations (L1 carrier)
            %   ph1_M = MASTER phase observations (L1 carrier)
            %   ph2_R = ROVER phase observations (L2 carrier)
            %   ph2_M = MASTER phase observations (L2 carrier)
            %   snr_R = ROVER-SATELLITE signal-to-noise ratio
            %   snr_M = MASTER-SATELLITE signal-to-noise ratio
            %   Eph   = satellite ephemeris
            %   sp3   = structure containing precise ephemeris and clock
            %   iono  = ionosphere parameters
            %   lambda = wavelength matrix (depending on the enabled constellations)
            %   phase  = L1 carrier (phase=1), L2 carrier (phase=2)
            %   ant_pcv = antenna phase center variation
            %
            % DESCRIPTION:
            %   Computation of the receiver position (X,Y,Z).
            %   Relative (double difference) positioning by least squares adjustment
            %   on code and phase observations.

            cutoff = this.state.cut_off;
            snr_threshold = this.state.snr_thr;
            cond_num_threshold = this.state.cond_num_thr;

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
                eph_avail = eph(30,:);
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
                        = init_positioning(time_rx, pr1_m(sat_pr),   snr_m(sat_pr),   eph, sp3, iono, [],  pos_m, [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, frequencies,       2, 0, 0, 0); %#ok<ASGLU>
                    if (sum(sat_pr_M) < min_nsat_LS); return; end
                    [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] ...
                        = init_positioning(time_rx, pr1_r(sat_pr_M), snr_r(sat_pr_M), eph, sp3, iono, [], pos_r, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, frequencies, flag_XR, 1, 0, 0); %#ok<ASGLU>
                else
                    [pos_m, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM] ...
                        = init_positioning(time_rx, pr2_m(sat_pr),   snr_m(sat_pr),   eph, sp3, iono, [],  pos_m, [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, frequencies,       2, 0, 0, 0); %#ok<ASGLU>
                    if (sum(sat_pr_M) < min_nsat_LS); return; end
                    [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] ...
                        = init_positioning(time_rx, pr2_r(sat_pr_M), snr_r(sat_pr_M), eph, sp3, iono, [], pos_r, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, frequencies, flag_XR, 1, 0, 0); %#ok<ASGLU>
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

                %apply cutoffs also to phase satellites
                sat_removed = setdiff(sat_pr_old, sat_pr);
                sat_ph(ismember(sat_ph,sat_removed)) = [];

                % keep only satellites that rover and master have in common both in phase and code
                [sat_pr, iR, iM] = intersect(sat_pr, sat_ph);
                XS = XS(iR,:);
                sys = sys(iR);
                if (~isempty(err_tropo_R))
                    err_tropo_R = err_tropo_R(iR);
                    err_iono_R  = err_iono_R (iR);
                    err_tropo_M = err_tropo_M(iM);
                    err_iono_M  = err_iono_M (iM);
                end

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
                [null_max_elR, pivot_index] = max(elR(sat_ph)); %#ok<ASGLU>
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
                if (size(sat_ph,1) >= min_nsat_LS && (isempty(cond_num) || cond_num < cond_num_threshold))

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

        function [phase_res, id_res] = computePhRes(this, v_hat, A, Q, obs_track, amb_prn_track)
            if (nargin == 1)
                v_hat = this.v_hat;
                A = this.A;
                Q = this.Q;
                obs_track = this.obs_track;
                amb_prn_track = this.amb_prn_track;
            end
            
            % EXAMPLE: phase_res = this.computePhRes(A, Q, obs_track, amb_prn_track)
            phase_res = nan(this.n_tot_epoch, numel(amb_prn_track), 2);
            id_res = spalloc(this.n_tot_epoch, numel(amb_prn_track), round(this.n_tot_epoch * numel(amb_prn_track) * 0.5));
            %err = zeros(length(amb_prn_track), 1);
            for a = 1 : length(amb_prn_track)
                idx = find(A(:, 3 + a) < 0);
                res = v_hat(idx);
                phase_res(obs_track(idx, 1), a, 1) = res;
                id_res(obs_track(idx, 1), a) = idx; %#ok<*SPRIX>
                %err(a) = perc(movstd(res(~isnan(res)),7),0.8);
                phase_res(obs_track(idx, 1), a, 2) = sqrt(Q(idx + size(Q,1) * (idx - 1)));
            end
            %err = median(err);
        end

    end

    methods (Static, Access = private)

        function [A, y0, b, Q, sat_track, amb_num, amb_prn_track] = remAmbCol(A, y0, b, Q, sat_track, amb_num, amb_prn_track, rem_amb)
            if (~isempty(rem_amb))
                rem_obs = [];
                for r = 1 : length(rem_amb)
                    rem_obs = [rem_obs; find(A(:,rem_amb(r))~=0)]; %#ok<AGROW>
                end
                A(rem_obs,:) = [];
                y0(rem_obs) = [];
                b(rem_obs) = [];
                Q(rem_obs,:) = []; Q(:,rem_obs) = [];
                sat_track(rem_obs,:) = [];

                A(:,rem_amb) = [];
                amb_num = amb_num - length(rem_amb);
                amb_prn_track(rem_amb-3) = [];
            end
        end

        function [id_out] = findBadObs(phase_res, thr, win_size)
            id_out = [];
            for a = 1 : size(phase_res, 2)
                y = phase_res(:, a, 1);
                if nargin == 3
                    mm = movmedian(y(~isnan(y)), win_size);
                else
                    mm = 0;
                end

                e_sup = y; e_sup(~isnan(y)) = mm + thr;
                e_inf = y; e_inf(~isnan(y)) = mm - thr;

                id_out = [id_out; (a-1) * size(phase_res, 1) + find(y > e_sup | y < e_inf)]; %#ok<AGROW>

            end
        end

        function plotPhResThr(phase_res, err, amb_prn_track)
            x = (1 : size(phase_res,1));
            for a = 1 : numel(amb_prn_track)
                h = figure(amb_prn_track(a));
                h.Name = sprintf('Sat: %d', amb_prn_track(a));
                h.NumberTitle = 'off';
                dockAllFigures();
                y = phase_res(:, a, 1);
                %e = 3*phase_res(:, a, 2);
                e_sup = y; e_sup(~isnan(y)) = movmedian(y(~isnan(y)), 9) + err;
                e_inf = y; e_inf(~isnan(y)) = movmedian(y(~isnan(y)), 9) - err;

                plot(x, y,'.-'); hold on;
                hline = findobj(h, 'type', 'line');

                patchColor = min(hline(1).Color + 0.3, 1);
                plot(x, e_sup, x, e_inf, 'color', patchColor);
                patch([x(~isnan(e_sup)) fliplr(x(~isnan(e_sup)))], [e_sup(~isnan(e_sup)); flipud(e_inf(~isnan(e_sup)))], 1, ...
                    'facecolor',patchColor, ...
                    'edgecolor','none', ...
                    'facealpha', 0.1);
                hold on;
            end
        end


    end

    % =========================================================================
    %  CONSTELLATION MANAGEMENT
    % =========================================================================

    methods % Public Access
        function [ cc ] = initConstellation(obj, GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)
            % Multi-constellation set-up.
            %
            % SYNTAX:
            %   cc = initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
            %   cc = initConstellation([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]);
            %
            % INPUT:
            %   single logical array whose elements are:
            %   GPS_flag = boolean flag for enabling/disabling GPS usage
            %   GLO_flag = boolean flag for enabling/disabling GLONASS usage
            %   GAL_flag = boolean flag for enabling/disabling Galileo usage
            %   BDS_flag = boolean flag for enabling/disabling BeiDou usage
            %   QZS_flag = boolean flag for enabling/disabling QZSS usage
            %   SBS_flag = boolean flag for enabling/disabling SBAS usage (for ranging)
            %
            % OUTPUT:
            %   the results is stored within the object referenced by "cc"
            switch nargin
                case 2,  enabled_ss = GPS_flag;
                case 6,  enabled_ss = logical([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, 0]);
                case 7,  enabled_ss = logical([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]);
                otherwise, error(['Initialization of Constellation_Collector failed: ' 10 '   invalid number of parameters in the constructor call']);
            end

            obj.cc = Costellation_Collector(enabled_ss);

            cc = obj.cc;
        end

    end

    methods % Public Access (Legacy support)
    end

end
