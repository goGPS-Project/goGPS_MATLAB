%   CLASS Least_Square_Manipulator
% =========================================================================
%
% DESCRIPTION
%   Manipulate least squares system
%
% EXAMPLE
%   LSM = Least_Square_Manipulator();
%
% SEE ALSO
%   - Least_Square
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 4 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
%  Contributors:     Andrea Gatti
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
classdef Least_Squares_Manipulator < handle
    
    % ==================================================================================================================================================
    %% Parameter columns id (order)
    % ==================================================================================================================================================
    
    properties (Constant)
        PAR_X = 1;
        PAR_Y = 2;
        PAR_Z = 3;
        PAR_ISB = 4;
        PAR_AMB = 5;
        PAR_REC_CLK = 6;
        PAR_TROPO = 7;
        PAR_TROPO_N = 8;
        PAR_TROPO_E = 9;
        PAR_PCO_X = 10;
        PAR_PCO_Y = 11;
        PAR_PCO_Z = 12;
        PAR_SAT_CLK = 13;
        CLASS_NAME = {'X', 'Y', 'Z', 'ISB', 'AMB', 'REC_CLK', 'TROPO', 'TROPO_N', 'TROPO_E', 'PCO_X', 'PCO_Y', 'PCO_Z', 'SAT_CLK'};
    end
    
    properties
        A_ep         % Stacked epoch-wise design matrices [n_obs x n_param_per_epoch]
        A_idx        % index of the paramter [n_obs x n_param_per_epoch]
        A_idx_mix    % inded of the paramter for multi receiver session [n_obs x n_param_per_epoch]
        amb_idx      % index of the columns per satellite
        go_id_amb    % go ids of the amb idx
        out_idx      % logical index to tell if an observation is an outlier [ n_obs x 1]
        N_ep         % Stacked epoch-wise normal matrices [ n_param_per_epoch x n_param_per_epoch x n_obs]
        G            % hard constraints (Lagrange multiplier)
        D            % known term of the hard constraints
        y            % observations  [ n_obs x 1]
        system_split % limits of the ambiguity splits
        variance     % observation variance [ n_obs x 1]
        rw           % reweight factor
        res          % observations residuals
        epoch        % epoch of the obseravtions and of the A lines [ n_obs x 1]
        sat          % satellite of the obseravtions and of the A lines [ n_obs x 1]
        n_epochs
        param_class  % [n_param_per_epoch x 1] each paramter can be part of a class
        %                 [   1 : x
        %                     2 : y
        %                     3 : z
        %                     4 : inter channel/frequency/system biases
        %                     5 : ambiguity
        %                     6 : clock
        %                     7 : tropo
        %                     8 : tropo inclination north
        %                     9 : tropo inclination east ]
        param_flag          % 0: constant in time always same param, -1: constant in time differents param (e.g ambiguity), 1: same param changing epochwise
        time_regularization % [ param_class time_varability] simple time regularization constructed from psudo obs p_ep+1 - p_ep = 0 with given accuracy
        mean_regularization
        true_epoch          % true epoch of the epoch-wise paramters
        rec_time_idxes      % for each receiver tell which epochof the common time are used
        rate                % rate of the true epoch
        sat_go_id           % go id of the sat indexes
        receiver_id         % id of the receiver, in case of differenced observations two columns are used
        
        
        wl_amb              % wide-lane ambuiguity
        wl_fixed            % is wide-lane fixed
        
        amb_set_jmp         % cell conatinin for each receiver the jmps on all ambiguity
        
        network_solution = false;
        
        sat_jmp_idx         % satelite jmp index
        
        pos_indexs_tc = {}  % to whivh index of the sampled time the progessive index correspond
        
        apriori_info % previous knowledge about the state to be estimated (for now only ambiguity)
        x_float
        Cxx_amb
    end
    
    properties (Access = private)
        Ncc         % part of the normal matrix with costant paramters
        Nee         % diagonal part of the normal matrix with epoch wise or multi epoch wise paramters
        Nce         % cross element between constant and epoch varying paramters
        
        rf          % reference frame
        
        state       % link to the current settings
        cc          % link to the current constellation collector
        
        log         % link to the logger
    end
    
    methods
        function this = Least_Squares_Manipulator(cc)
            % Creator Brahma
            this.init();
            if nargin > 0
                this.cc = cc;
            end
        end
        
        function init(this)
            % Init the singletons (create links)
            %
            % SYNTAX
            %   this.init();
            this.state = Core.getState();
            this.rf = Core.getReferenceFrame();
            this.log = Core.getLogger();            
        end
        
        function id_sync = setUpPPP(this, rec, id_sync,  cut_off, dynamic, pos_idx)
            % Init the object for the phase stand alone positioning
            %
            % SYNTAX
            %   id_sync = this.setUpPPP(rec, id_sync,  <cut_off>, <dynamic>, <pos_idx>)
            if nargin < 4
                cut_off = [];
            end
            if nargin < 5
                dynamic = false;
            end
            if nargin < 6
                pos_idx = [];
            end
            id_sync = this.setUpSA(rec, id_sync, 'L', cut_off, '', dynamic, pos_idx);
        end
        
        function id_sync = setUpCodeSatic(this, rec, id_sync, cut_off)
            % Init the object for the code static positioning
            %
            % SYNTAX
            %   id_sync = setUpCodeSatic(this, rec, id_sync, <cut_off>)
            if nargin < 4
                cut_off = [];
            end
            id_sync = this.setUpSA(rec, id_sync, 'C', cut_off);
        end
        
        function id_sync = setUpCodeDynamic(this, rec, id_sync, cut_off)
            if nargin < 4
                cut_off = [];
            end
            id_sync = this.setUpSA(rec, id_sync, 'C', cut_off, '', true);
        end
        
        function id_sync_out = setUpSA(this, rec, id_sync_in, obs_type, cut_off, custom_obs_set, dynamic, pos_idx_vec)
            % Init the object for static positioning
            % return the id_sync of the epochs to be computed
            %
            % INPUT:
            %    rec : receiver
            %    id_sync : epoch to be used
            %    obs_type : 'C' 'L' 'CL'
            %    cut_off : cut off angle [optional]
            %
            % SYNTAX
            %   id_sync_out = this.setUpSA(rec, id_sync_in, obs_type('C'/'L'/'CL'), cut_off, custom_obs_set, <dynamic>, <pos_idx_vec>)
            if nargin < 8
                pos_idx_vec = [];
            end
            if nargin < 7
                dynamic = false;
            end
            
            % Extract the observations to be used for the solution
            phase_present = instr(obs_type, 'L');
            flag_amb_fix = this.state.flag_amb_fix;
            if nargin < 6 || isempty(custom_obs_set)
                obs_set = Observation_Set();
                if rec.isMultiFreq() && ~rec.state.isIonoExtModel %% case multi frequency
                    
                    % Using smoothed iono fromg geometry free
                    for sys_c = rec.cc.sys_c
                        for i = 1 : length(obs_type)
                            if this.state.isIonoFree || ~phase_present
                                obs_set.merge(rec.getPrefIonoFree(obs_type(i), sys_c));
                            else
                                obs_set.merge(rec.getSmoothIonoFreeAvg('L', sys_c));
                                obs_set.iono_free = true;
                                obs_set.sigma = obs_set.sigma * 1.5;
                            end
                        end
                    end
                    
                    if flag_amb_fix && phase_present
                        [this.wl_amb, this.wl_fixed, wsb_rec]  = rec.getWidelaneAmbEst();
                        f_vec = GPS_SS.F_VEC;
                        l_vec = GPS_SS.L_VEC;
                        obs_set.obs = nan2zero(obs_set.obs - (this.wl_amb(:,obs_set.go_id))*f_vec(2)^2*l_vec(2)/(f_vec(1)^2 - f_vec(2)^2));
                    end
                else
                    % Using the best combination available
                    for sys_c = rec.cc.sys_c
                        f = rec.getFreqs(sys_c);
                        for i = 1 : length(obs_type)
                            if ~isempty(f)
                                obs_set.merge(rec.getPrefObsSetCh([obs_type(i) num2str(f(1))], sys_c));
                            end
                        end
                    end
                    idx_ph = obs_set.obs_code(:, 2) == 'L';
                    if sum(idx_ph) > 0
                        obs_set.obs(:, idx_ph) = obs_set.obs(:, idx_ph) .* repmat(serialize(obs_set.wl(idx_ph))', size(obs_set.obs,1), 1);
                    end
                end
            else
                obs_set = custom_obs_set;
            end
            
            if phase_present
                n_sat = rec.cc.getMaxNumSat();
                rec.sat.o_cs_ph = zeros(rec.time.length, n_sat);
                rec.sat.o_cs_ph(:,obs_set.go_id) = obs_set.cycle_slip;
                rec.sat.o_out_ph = zeros(rec.time.length, n_sat);
                dual_freq = size(obs_set.obs_code,2) > 5;
                [~, ~, ph_idx] = rec.getPhases();
                obs_code_ph = rec.obs_code(ph_idx,:);
                go_id_ph = rec.go_id(ph_idx);
                for s = 1 : length(obs_set.go_id);
                    g = obs_set.go_id(s);
                    obs_code1 = obs_set.obs_code(s,2:4);
                    if dual_freq
                        obs_code2 = obs_set.obs_code(s,5:7);
                    else
                        obs_code2 = '   ';
                    end
                    out_idx = strLineMatch(obs_code_ph,obs_code1) & go_id_ph == g;
                    out = rec.sat.outlier_idx_ph(:,out_idx);
                    if strcmp(obs_code2,'   ')
                        out_idx = strLineMatch(obs_code_ph,obs_code2) & go_id_ph == g;
                        if any(out_idx)
                            out(:,out_idx) = out(:,out_idx) | rec.sat.outlier_idx_ph(:,out_idx);
                        end
                    end
                    rec.sat.o_out_ph(:,g) = out;
                    
                end
            end
            
            % if phase observations are present check if the computation of troposphere parameters is required
            
            if phase_present
                tropo = this.state.flag_tropo;
                tropo_g = this.state.flag_tropo_gradient;
            else
                tropo = false;
                tropo_g = false;
            end
            
            % check presence of snr data and fill the gaps if needed
            if ~isempty(obs_set.snr)
                snr_to_fill = (double(obs_set.snr ~= 0) + 2 * double(obs_set.obs ~= 0)) == 2; % obs if present but snr is not
                if sum(sum(snr_to_fill))
                    obs_set.snr = simpleFill1D(obs_set.snr, snr_to_fill);
                end
            end
            
            % remove epochs based on desired sampling
            if nargin > 2
                obs_set.keepEpochs(id_sync_in);
            end
            
            % re-apply cut off if requested
            if nargin > 4 && ~isempty(cut_off) && sum(sum(obs_set.el)) ~= 0
                obs_set.remUnderCutOff(cut_off);
            end
            
            
            % remove not valid empty epoch or with only one satellite (probably too bad conditioned)
            idx_valid_ep_l = sum(obs_set.obs ~= 0, 2) > 0;
            obs_set.setZeroEpochs(~idx_valid_ep_l);
            
            %remove too shortArc
            
            % Compute the number of ambiguities that must be estimated
            cycle_slip = obs_set.cycle_slip;
            min_arc = this.state.getMinArc;
            if phase_present && min_arc > 1
                amb_idx = obs_set.getAmbIdx();
                % amb_idx = n_coo + n_iob + amb_idx;
                amb_idx = zero2nan(amb_idx);
                
                % remove short arcs
                
                % ambiguity number for each satellite
                amb_obs_count = histcounts(serialize(amb_idx), 'Normalization', 'count', 'BinMethod', 'integers');
                assert(numel(amb_obs_count) == max(amb_idx(:))); % This should always be true
                id = 1 : numel(amb_obs_count);
                ko_amb_list = id(amb_obs_count < min_arc);
                amb_obs_count(amb_obs_count < min_arc) = [];
                for ko_amb = fliplr(ko_amb_list)
                    id_ko = amb_idx == ko_amb;
                    obs_set.remObs(id_ko,false)
                end
                
            end
            idx_valid_ep_l = sum(obs_set.obs ~= 0, 2) > 0;
            obs_set.setZeroEpochs(~idx_valid_ep_l);
            
            obs_set.remEmptyColumns();
            
            % get reference observations and satellite positions
            [synt_obs, xs_loc] = rec.getSyntTwin(obs_set);
            xs_loc = zero2nan(xs_loc);
            diff_obs = nan2zero(zero2nan(obs_set.obs) - zero2nan(synt_obs));
            
            
            % Sometime code observations may contain unreasonable values -> remove them
            if obs_type == 'C'
                % very coarse outlier detection based on diff obs
                mean_diff_obs = mean(mean(abs(diff_obs),'omitnan'),'omitnan');
                diff_obs(abs(diff_obs) > 50 * mean_diff_obs) = 0;
            end
            idx_empty_ep = sum(diff_obs~=0,2) <= 1;
            diff_obs(idx_empty_ep, :) = [];
            xs_loc(idx_empty_ep, :, :) = [];
            obs_set.remEpochs(idx_empty_ep);
            obs_set.sanitizeEmpty();
            idx_empty_col = sum(diff_obs) == 0;
            diff_obs(:, idx_empty_col) = [];
            xs_loc(:, idx_empty_col, :) = [];
            
            if phase_present
                amb_idx = obs_set.getAmbIdx();
                n_amb = max(max(amb_idx));
                this.amb_idx = amb_idx;
                amb_flag = 1;
                this.go_id_amb = obs_set.go_id;
            else
                
                n_amb = 0;
                amb_flag = 0;
                this.amb_idx = [];
            end
            
            % set up requested number of parametrs
            n_epochs = size(obs_set.obs, 1);
            this.n_epochs = n_epochs;
            n_stream = size(diff_obs, 2); % number of satellites
            n_clocks = n_epochs; % number of clock errors
            n_tropo = n_clocks; % number of epoch for ZTD estimation
            ep_p_idx = 1 : n_clocks; % indexes of epochs starting from 1 to n_epochs
            
            this.true_epoch = obs_set.getTimeIdx(rec.time.first, rec.getRate); % link between original epoch, and epochs used here
            id_sync_out = this.true_epoch;
            is_fixed = rec.isFixed() || (rec.hasGoodApriori && ~phase_present);
            n_coo_par =  ~is_fixed * 3; % number of coordinates
            
            if dynamic
                n_coo = n_coo_par * n_epochs;
            else
                if isempty(pos_idx_vec)
                    n_coo = n_coo_par;
                else
                    n_pos = length(unique(pos_idx_vec));
                    n_coo =  n_pos*3;
                end
            end
            
            % get the list  of observation codes used
            u_obs_code = cell2mat(unique(cellstr(obs_set.obs_code)));
            n_u_obs_code = size(u_obs_code, 1);
            this.log.addMessage(this.log.indent(sprintf('Setting up SA system usign %s', reshape([u_obs_code repmat(' ',n_u_obs_code,1)]', 1, n_u_obs_code * (size(u_obs_code, 2) + 1)))));
            % if multiple observations types are present inter observations biases need be compouted
            iob_idx = zeros(size(obs_set.wl));
            for c = 1 : size(u_obs_code, 1)
                idx_b = strLineMatch(obs_set.obs_code, u_obs_code(c, :));
                iob_idx(idx_b) = c - 1;
            end
            iob_p_idx = iob_idx + n_coo; % progressive index start for iob
            n_iob = size(u_obs_code, 1) - 1;
            iob_flag = double(n_iob > 0);
            
            % separte antenna phase centers
            apc_flag = this.state.isSeparateApc() & phase_present ;
            n_apc = 3 * n_iob * apc_flag;
            apc_p_idx =  zeros(length(obs_set.wl),3);
            u_sys_c = unique(obs_set.obs_code(:,1));
            idx_gps = u_sys_c == 'G';
            if sum(idx_gps) > 0 % put gps in first position if presetn
                u_sys_c(idx_gps) = [];
                u_sys_c = ['G'; u_sys_c(:)];
            end
            for i = 1 : length(u_sys_c)
                sys_idx = obs_set.obs_code(:,1) == u_sys_c(i);
                apc_p_idx(sys_idx,:) = n_coo + n_iob + repmat(max((i-2),0)*3 + (1:3), sum(sys_idx),1);
            end
            % total number of observations
            n_obs = sum(sum(diff_obs ~= 0));
            
            % Building Design matrix
            n_par = n_coo_par + iob_flag + 3 * apc_flag + amb_flag + 1 + double(tropo) + 2 * double(tropo_g); % three coordinates, 1 clock, 1 inter obs bias(can be zero), 1 amb, 3 tropo paramters
            A = zeros(n_obs, n_par); % three coordinates, 1 clock, 1 inter obs bias(can be zero), 1 amb, 3 tropo paramters
            obs = zeros(n_obs, 1);
            sat = zeros(n_obs, 1);
            
            A_idx = zeros(n_obs, n_par);
            if ~is_fixed
                if isempty(pos_idx_vec)
                    A_idx(:, 1:3) = repmat([1, 2, 3], n_obs, 1);
                end
            end
            y = zeros(n_obs, 1);
            variance = zeros(n_obs, 1);
            obs_count = 1;
            this.sat_go_id = obs_set.go_id;
            
            % Getting mapping faction values
            if tropo || tropo_g
                [~, mfw] = rec.getSlantMF(id_sync_out);
                mfw(mfw  > 60 ) = nan;
                %mfw = mfw(id_sync_out,:); % getting only the desampled values
            end
            
            for s = 1 : n_stream
                id_ok_stream = diff_obs(:, s) ~= 0; % check observation existence -> logical array for a "s" stream
                
                obs_stream = diff_obs(id_ok_stream, s);
                % snr_stream = obs_set.snr(id_ok_stream, s); % SNR is not currently used
                if tropo || tropo_g
                    el_stream = obs_set.el(id_ok_stream, s) / 180 * pi;
                    az_stream = obs_set.az(id_ok_stream, s) / 180 * pi;
                    mfw_stream = mfw(id_ok_stream, obs_set.go_id(s)); % A simpler value could be 1./sin(el_stream);
                end
                xs_loc_stream = permute(xs_loc(id_ok_stream, s, :), [1, 3, 2]);
                los_stream = rowNormalize(xs_loc_stream);
                
                n_obs_stream = length(obs_stream);
                lines_stream = obs_count + (0:(n_obs_stream - 1));
                
                %--- Observation related vectors------------
                obs(lines_stream) = ep_p_idx(id_ok_stream);
                sat(lines_stream) = s;
                y(lines_stream) = obs_stream;
                variance(lines_stream) =  obs_set.sigma(s)^2;
                % ----------- FILL IMAGE MATRIX ------------
                % ----------- coordinates ------------------
                if ~is_fixed
                    A(lines_stream, 1:3) = - los_stream;
                end
                prog_p_col = 0;
                if dynamic & ~is_fixed
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = ep_p_idx(id_ok_stream);
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = n_epochs + ep_p_idx(id_ok_stream);
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = 2*n_epochs + ep_p_idx(id_ok_stream);
                elseif ~isempty(pos_idx_vec)
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = pos_idx_vec(id_ok_stream);
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = n_pos  + pos_idx_vec(id_ok_stream);
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = 2*n_pos + pos_idx_vec(id_ok_stream);
                elseif  ~is_fixed
                    prog_p_col = prog_p_col + 3 ;
                end
                % ----------- Inster observation bias ------------------
                if n_iob > 0
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = iob_idx(s) > 0;
                    A_idx(lines_stream, prog_p_col) = max(n_coo+1, iob_p_idx(s));
                end
                % ----------- Separate antenna phase centers ------------------
                if n_apc > 0
                    if obs_set.obs_code(s,1) ~= u_sys_c(1)
                        A(lines_stream, prog_p_col + ( 1 : 3)) = - los_stream;
                    end
                    A_idx(lines_stream,prog_p_col + ( 1 : 3)) = repmat(apc_p_idx(s,:),length(lines_stream),1);
                    prog_p_col = prog_p_col +3;
                end
                % ----------- Abiguity ------------------
                if phase_present
                    prog_p_col = prog_p_col + 1;
                    if not(flag_amb_fix)
                        A(lines_stream, prog_p_col) = 1;%obs_set.wl(s);
                    else
                        A(lines_stream, prog_p_col) = 1;
                    end
                    A_idx(lines_stream, prog_p_col) = n_coo + n_iob + n_apc + amb_idx(id_ok_stream, s);
                end
                % ----------- Clock ------------------
                prog_p_col = prog_p_col + 1;
                A(lines_stream, prog_p_col) = 1;
                A_idx(lines_stream, prog_p_col) = n_coo + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                % ----------- ZTD ------------------
                if tropo
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = mfw_stream;
                    A_idx(lines_stream, prog_p_col) = n_coo + n_clocks + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                end
                % ----------- ZTD gradients ------------------
                if tropo_g
                    %cotan_term = cot(el_stream) .* mfw_stream;
                    cotan_term = 1 ./ ( sin(el_stream).*tan(el_stream) + 0.0032);
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = cos(az_stream) .* cotan_term; % noth gradient
                    A_idx(lines_stream, prog_p_col) = n_coo + 2 * n_clocks + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = sin(az_stream) .* cotan_term; % east gradient
                    A_idx(lines_stream, prog_p_col) = n_coo + 3 * n_clocks + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                end
                obs_count = obs_count + n_obs_stream;
            end
            % ---- Suppress weighting until solution is more stable/tested
            %w(:) = 1;%0.005;%this.state.std_phase;
            %---------------------
            
            %---- Set up the constraint to solve the rank deficeny problem --------------
            if phase_present
                % Ambiguity set
                %G = [zeros(1, n_coo + n_iob) (amb_obs_count) -sum(~isnan(this.amb_idx), 2)'];
                %                 if ~flag_amb_fix
                %                     G = [zeros(1, n_coo + n_iob + n_apc)  ones(1,n_amb)  -ones(1,n_clocks)]; % <- This is the right one !!!
                %                 else % in case of ambiugty fixing with cnes orbit the partial trace minimization condition gives problems
                % setting the first clock of each connected set of arc to 0
                system_jmp = find([sum(nan2zero(diff(amb_idx)),2)] == sum(~isnan(amb_idx(1 : end - 1, :)),2) | [sum(nan2zero(diff(amb_idx)),2)] == sum(~isnan(amb_idx(2 : end, :)),2));
                clock_const = zeros(1,n_clocks);
                clock_const([1]) = 1;
                G = [zeros(1, n_coo + n_iob + n_apc)  zeros(1,n_amb)  clock_const];
                for i = 1: length(system_jmp)
                    clock_const = zeros(1,n_clocks);
                    clock_const(system_jmp(i)+1) = 1;
                    G = [G ;[zeros(1, n_coo + n_iob + n_apc)  zeros(1,n_amb)  clock_const]];
                end
                %                 end
                if tropo
                    G = [G zeros(size(G,1), n_clocks)];
                end
                if tropo_g
                    G = [G zeros(size(G,1), 2*n_clocks)];
                end
                D = zeros(size(G,1),1);
                this.G = G;
                this.D = D;
            end
            this.A_ep = A;
            this.A_idx = A_idx;
            this.variance = variance;
            this.y = y;
            this.receiver_id = ones(size(y));
            this.epoch = obs;
            this.sat = sat;
            if dynamic
                this.param_flag = [1, 1, 1, -ones(iob_flag), -repmat(ones(apc_flag),1,3), -ones(amb_flag), 1, ones(tropo), ones(tropo_g), ones(tropo_g)];
            else
                this.param_flag = [zeros(1,n_coo_par) -ones(iob_flag), -repmat(ones(apc_flag),1,3), -ones(amb_flag), 1, ones(tropo), ones(tropo_g), ones(tropo_g)];
            end
            this.param_class = [this.PAR_X*ones(~is_fixed) , this.PAR_Y*ones(~is_fixed), this.PAR_Z*ones(~is_fixed), this.PAR_ISB * ones(iob_flag), this.PAR_PCO_X * ones(apc_flag), this.PAR_PCO_Y * ones(apc_flag), this.PAR_PCO_Z * ones(apc_flag), this.PAR_AMB*ones(amb_flag), this.PAR_REC_CLK, this.PAR_TROPO*ones(tropo), this.PAR_TROPO_N*ones(tropo_g), this.PAR_TROPO_E*ones(tropo_g)];
            if phase_present
                system_jmp = find([sum(nan2zero(diff(amb_idx)),2)] == sum(~isnan(amb_idx(1 : end - 1, :)),2) | [sum(nan2zero(diff(amb_idx)),2)] == sum(~isnan(amb_idx(2 : end, :)),2));
                fprintf('#### DEBUG #### \n');
                [[1; system_jmp + 1] [system_jmp; max(obs)]]
                this.system_split = [[1; system_jmp + 1] [system_jmp; max(obs)]];
            else
                this.system_split = [1 max(obs)];
            end
        end
        
        function [common_time, id_sync]  = setUpNetworkAdj(this, rec_list, coo_rate)
            % NOTE : free netwrok is set up -> soft constarint on apriori coordinates to be implemnted
            %
            % OUTPUT:
            % common_time : common gps time used for the processing
            n_rec = length(rec_list);
            obs_set_list  = Observation_Set.empty(n_rec,0);
            work_list = [rec_list.work];
            this.cc = rec_list(1).cc;
            
            %--- for each receiver get one observation set
            for i = 1 : n_rec
                obs_set_list(i) = Observation_Set();
                if work_list(i).isMultiFreq() && ~work_list(i).state.isIonoExtModel %% case multi frequency
                    for sys_c = rec_list(i).work.cc.sys_c
                        if this.state.isIonoFree
                            obs_set_list(i).merge(work_list(i).getPrefIonoFree('L', sys_c));
                        else
                            obs_set_list(i).merge(work_list(i).getSmoothIonoFreeAvg('L', sys_c));
                            obs_set_list(i).iono_free = true;
                            obs_set_list(i).sigma = obs_set_list(i).sigma * 1.5;
                        end
                    end
                else
                    for sys_c = work_list(i).cc.sys_c
                        f = work_list(i).getFreqs(sys_c);
                        if ~isempty(f)
                            obs_set_list(i).merge(work_list(i).getPrefObsSetCh(['L' num2str(f(1))], sys_c));
                        end
                    end
                    idx_ph = obs_set_list(i).obs_code(:, 2) == 'L';
                    obs_set_list(i).obs(:, idx_ph) = obs_set_list(i).obs(:, idx_ph) .* repmat(obs_set_list(i).wl(idx_ph)', size(obs_set_list(i).obs,1),1);
                end
                obs_set_list(i).sanitizeEmpty();
                
                if this.state.flag_amb_pass && this.state.getCurSession > 1
                    % remove the left buffer it is necessary only to determine if a cycle slip ha occured in the first useful epoch
                    [~,limc] = this.state.getSessionLimits(this.state.getCurSession);
                    idx_rm = obs_set_list(i).time < limc.first;
                    obs_set_list(i).remEpochs(idx_rm);
                    % remove the first epoch
                    obs_set_list(i).sanitizeEmpty();
                    if ~isempty(this.apriori_info)
                        if obs_set_list(i).time.length > 0
                            for s = 1:length(obs_set_list(i).go_id)
                                if obs_set_list(i).cycle_slip(1,s) % if there is a ccle slip remove the stored ambigutiy
                                    idx_apr_info = this.apriori_info.goids == obs_set_list(i).go_id(s) & this.apriori_info.receiver == i;
                                    if sum(idx_apr_info) > 0
                                        this.removeAprInfo(idx_apr_info);
                                    end
                                    
                                end
                            end
                        else
                            idx_apr_info = this.apriori_info.receiver == i;
                            if sum(idx_apr_info) > 0
                                this.removeAprInfo(idx_apr_info);
                            end
                        end
                    end
                end
            end
            
           
            
            
            % Sync obs_sets
            sanitized = false;
            while ~sanitized
                % remove short arcs and remove "empty" satellites
                for r = 1 : n_rec
                    obs_set_list(r).remShortArc(max(this.state.getMinArc, 1));
                    obs_set_list(r).sanitizeEmpty();
                end
                [common_time, id_sync] = obs_set_list.getSyncTimeExpanded();
                % remove satellites arcs seen only by one receiver and sigle epochs arcs
                common_obs_mat = obs_set_list.getMRObsMat(common_time, id_sync);
                
                id_rm = sum(~isnan(common_obs_mat),3) == 1;
                valid_epoch =  sum(~isnan(common_obs_mat),3) > 1;
                for i = 1 : size(id_rm,2)
                    % get one epoch arcs
                    [flag_intervals] = getOutliers(valid_epoch(:,i));
                    single_arcs = (flag_intervals(:,2) - flag_intervals(:,1))  == 0;
                    idx_sa_rm = flag_intervals(single_arcs,1);
                    if ~isempty(idx_sa_rm)
                        id_rm(idx_sa_rm, i) = true;
                    end
                end
                % rem obs in each ob set
                for r = 1 : n_rec
                    id_rm_o = false(size(obs_set_list(r).obs));
                    id_rm_o(id_sync(~isnan(id_sync(:,r)),r),:) = id_rm(~isnan(id_sync(:,r)), obs_set_list(r).go_id);
                    obs_set_list(r).remObs(id_rm_o, false);
                end
                
                
                % keep the epochs common to at least 2 receivers and that sees more than 2 staellites
                id_ok = sum(~isnan(id_sync), 2) >= 2 & sum(sum(~isnan(common_obs_mat),3) > 1,2) > 2; % NOTE: if a constraint on clock is applied the trehsold might be lowered to 1 satellite, the gain will be probabily low for consumer receiver
                
                if any(~id_ok)
                    id_sync = id_sync(id_ok, :);
                    % filter the observation sets
                    for r = 1 : n_rec
                        obs_set_list(r).keepEpochs(noNaN(id_sync(:,r)));
                    end
                elseif sum(sum(id_rm)) == 0
                    sanitized = true;
                end
            end
            
            n_sat = 0;
            for r = 1 : n_rec
                n_sat = max([n_sat; obs_set_list(r).go_id(:)]);
            end
            
             %--- for each satellite checks epochs for which all receiver-satellite observation continuity is broken
            sat_jmp_idx = true(size(id_sync, 1), n_sat);
            for r = 1 : n_rec
                sat_jmp_rec = obs_set_list(r).getArcJmpMat(id_sync(:,r));
                sat_jmp_idx(:, obs_set_list(r).go_id) = sat_jmp_idx(:, obs_set_list(r).go_id) & sat_jmp_rec;
            end
            
         

            % get the observation equation for each receiver
            A = []; Aidx = []; ep = []; sat = []; p_flag = []; p_class = []; y = []; variance = []; r = [];
            [sss_lim, ~] = this.state.getSessionLimits(this.state.getCurSession());
            st_time = sss_lim.first;
            
            for i = 1 : n_rec
                % get the poitiong idx
                if ~isempty(coo_rate) && i~=1 % first receiver do not need any sub rate since is the reference
                    [pos_idx_nh, pos_idx_tc] = Least_Squares_Manipulator.getPosIdx(obs_set_list(i).time, st_time, coo_rate);
                    this.pos_indexs_tc{end+1} = pos_idx_tc; % to be used afterwards to push back postions
                else
                    pos_idx_nh = [];
                end
                [A_rec, Aidx_rec, ep_rec, sat_rec, p_flag_rec, p_class_rec, y_rec, variance_rec, amb_set_jmp] = this.getObsEq(rec_list(i).work, obs_set_list(i), pos_idx_nh);
                A = [A ; A_rec];
                Aidx = [Aidx; Aidx_rec];
                r2c = find(~isnan(id_sync(:,i)));
                ep = [ep; r2c(ep_rec)];
                sat = [sat; sat_rec];
                p_flag = p_flag_rec;
                p_class = p_class_rec;
                y = [y; y_rec];
                variance = [variance; variance_rec];
                r = [r; ones(size(y_rec))*i];
                this.amb_set_jmp{i} = r2c(amb_set_jmp);
            end
            
            p_flag(p_flag == 0) = -1;
                        
            this.A_idx = Aidx;
            
            % get the number of used epochs for each receiver
            this.n_epochs = zeros(n_rec,1);
            for i = 1 : n_rec
                this.n_epochs(i) = length(unique(this.A_idx(r == i,p_class == this.PAR_REC_CLK)));
            end
            this.A_ep = A;
            
            this.variance = variance;
            this.y = y;
            this.epoch = ep;
            time = common_time;
            this.true_epoch = round(time.getNominalTime().getRefTime()/time.getRate) + 1; %
            this.rec_time_idxes = id_sync;
            this.rate = time.getRate;
            this.sat = sat;
            this.param_flag = p_flag;
            this.param_class = p_class;
            this.receiver_id = r;
            this.sat_jmp_idx = sat_jmp_idx;
            
            this.network_solution = true;   
            
            
            % remove ambiguity fixed previously
            if this.state.flag_amb_pass
                for i = 1 : n_rec
                end
            end
            
        end
        
        function changePosIdx(this, r_id, pos_idx)
            % Change the index of the position in the design matrix
            %
            % SYNTAX:
            %   this.changePosIdx(pos_idx)
            max_idx = max(pos_idx);
            idx_rec = this.receiver_id == r_id;
            o_max_idx = max(this.A_idx(idx_rec, 1)) - min(this.A_idx(idx_rec, 1)) +1;
            pos_idx_coo = [pos_idx pos_idx+max_idx pos_idx+max_idx*2];
            
            this.A_idx(idx_rec, 1:3) = pos_idx_coo(this.epoch(idx_rec), :);
            this.A_idx(idx_rec, 4:end) = this.A_idx(idx_rec, 4:end) + 3*(max_idx-o_max_idx);
        end
                        
        function [A, A_idx, ep, sat, p_flag, p_class, y, variance, amb_set_jmp] = getObsEq(this, rec, obs_set, pos_idx_vec)
            % get reference observations and satellite positions
            
            tropo = rec.state.flag_tropo;
            tropo_g = rec.state.flag_tropo_gradient;
            dynamic = rec.state.rec_dyn_mode;
            global_sol = false; %this.state.global_solution;
            amb_flag = true;
            phase_present = true;
            
            [synt_obs, xs_loc] = rec.getSyntTwin(obs_set);
            xs_loc = zero2nan(xs_loc);
            diff_obs = nan2zero(zero2nan(obs_set.obs) - zero2nan(synt_obs));
            % diff_obs = nan2zero(zero2nan(diff_obs) - cumsum(median(Core_Utils.diffAndPred(zero2nan(diff_obs)), 2, 'omitnan'))); % for DEBUGGING: remove receiver clock
            
            amb_idx = obs_set.getAmbIdx();
            n_amb = max(max(amb_idx));
            idx = [];
            
            id_sync_out = obs_set.getTimeIdx(rec.time.first, rec.getRate);
            
            % set up requested number of parametrs
            n_epochs = size(obs_set.obs, 1);
            n_stream = size(diff_obs, 2); % number of satellites
            n_clocks_rec = n_epochs; % number of clock errors
            n_tropo = n_clocks_rec; % number of epoch for ZTD estimation
            ep_p_idx = 1 : n_clocks_rec; % indexes of epochs starting from 1 to n_epochs
            
            n_coo_par =  ~rec.isFixed() * 3; % number of coordinates
            
            if dynamic
                n_coo = n_coo_par * n_epochs;
            else
                if isempty(pos_idx_vec)
                    n_coo = n_coo_par;
                else
                    n_pos = length(unique(pos_idx_vec));
                    n_coo =  n_pos*3;
                end
            end
            
            % get the list  of observation codes used
            u_obs_code = cell2mat(unique(cellstr(obs_set.obs_code)));
            % n_u_obs_code = size(u_obs_code, 1);
            
            % if multiple observations types are present in the observations, biases need be computed
            iob_idx = zeros(size(obs_set.wl));
            for c = 1 : size(u_obs_code, 1)
                idx_b = strLineMatch(obs_set.obs_code, u_obs_code(c, :));
                iob_idx(idx_b) = c - 1;
            end
            iob_p_idx = iob_idx + n_coo; % progressive index start for iob
            n_iob = size(u_obs_code, 1) - 1;
            iob_flag = double(n_iob > 0);
            
            % separate antenna phase centers
            apc_flag = rec.state.isSeparateApc();
            n_apc = 3 * n_iob * apc_flag;
            apc_p_idx =  zeros(length(obs_set.wl),3);
            u_sys_c = unique(obs_set.obs_code(:,1));
            idx_gps = u_sys_c == 'G';
            if sum(idx_gps) > 0 % put gps in first position if present
                u_sys_c(idx_gps) = [];
                u_sys_c = ['G'; u_sys_c(:)];
            end
            for i = 1 : length(u_sys_c)
                sys_idx = obs_set.obs_code(:,1) == u_sys_c(i);
                apc_p_idx(sys_idx,:) = n_coo + n_iob + repmat(max((i-2),0)*3 + (1:3), sum(sys_idx),1);
            end
            % total number of observations
            n_obs = sum(sum(diff_obs ~= 0));
            
            if global_sol
                n_sat_clk_par = 1;
            else
                n_sat_clk_par = 0;
            end
            
            % Building Design matrix
            n_par = n_coo_par + iob_flag + 3 * apc_flag + amb_flag + 1 + double(tropo) + 2 * double(tropo_g) + n_sat_clk_par; % three coordinates, 1 clock, 1 inter obs bias(can be zero), 1 amb, 3 tropo paramters
            A = zeros(n_obs, n_par); % three coordinates, 1 clock, 1 inter obs bias(can be zero), 1 amb, 3 tropo paramters
            ep = zeros(n_obs, 1);
            sat = zeros(n_obs, 1);
            
            A_idx = zeros(n_obs, n_par);
            if ~rec.isFixed()
                if isempty(pos_idx_vec)
                    A_idx(:, 1:3) = repmat([1, 2, 3], n_obs, 1);
                end
            end
            y = zeros(n_obs, 1);
            variance = zeros(n_obs, 1);
            obs_count = 1;
            this.sat_go_id = obs_set.go_id;
            
            % Getting mapping faction values
            if tropo || tropo_g
                [~, mfw] = rec.getSlantMF(id_sync_out);
                mfw(mfw  > 60 ) = nan;
                %mfw = mfw(id_sync_out,:); % getting only the desampled values
            end
            
%             %%% DEBUG PURPOSES: removing an empirically estimated clock
%             sensor = Core_Utils.diffAndPred(zero2nan(diff_obs));
%             sensor = median(sensor, 2, 'omitnan'); % using median to disregard cycle-slips
%             % rough estimation of clock, the median is not a good estimator
%             % but for now it could stay like this
%             clock_diff = detrend(cumsum(sensor));
%             diff_obs = nan2zero(bsxfun(@minus, zero2nan(diff_obs), clock_diff));
%             %%% END
            
           
            
            for s = 1 : n_stream
                id_ok_stream = diff_obs(:, s) ~= 0; % check observation existence -> logical array for a "s" stream
                
                obs_stream = diff_obs(id_ok_stream, s);
                % snr_stream = obs_set.snr(id_ok_stream, s); % SNR is not currently used
                if tropo || tropo_g
                    el_stream = obs_set.el(id_ok_stream, s) / 180 * pi;
                    az_stream = obs_set.az(id_ok_stream, s) / 180 * pi;
                    mfw_stream = mfw(id_ok_stream, obs_set.go_id(s)); % A simpler value could be 1./sin(el_stream);
                end
                xs_loc_stream = permute(xs_loc(id_ok_stream, s, :), [1, 3, 2]);
                los_stream = rowNormalize(xs_loc_stream);
                
                n_obs_stream = length(obs_stream);
                lines_stream = obs_count + (0:(n_obs_stream - 1));
                
                %--- Observation related vectors------------
                ep(lines_stream) = ep_p_idx(id_ok_stream);
                sat(lines_stream) = obs_set.go_id(s);
                y(lines_stream) = obs_stream;
                variance(lines_stream) =  obs_set.sigma(s)^2;
                % ----------- FILL IMAGE MATRIX ------------
                % ----------- coordinates ------------------
                if ~rec.isFixed()
                    A(lines_stream, 1:3) = - los_stream;
                end
                prog_p_col = 0;
                if dynamic && ~rec.isFixed()
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = ep_p_idx(id_ok_stream);
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = n_epochs + ep_p_idx(id_ok_stream);
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = 2*n_epochs + ep_p_idx(id_ok_stream);
                elseif ~isempty(pos_idx_vec)
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = pos_idx_vec(id_ok_stream);
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = n_pos  + pos_idx_vec(id_ok_stream);
                    prog_p_col = prog_p_col +1;
                    A_idx(lines_stream, prog_p_col) = 2*n_pos + pos_idx_vec(id_ok_stream);
                elseif  ~rec.isFixed()
                    prog_p_col = prog_p_col + 3 ;
                end
                % ----------- Inter observation bias ------------------
                if n_iob > 0
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = iob_idx(s) > 0;
                    A_idx(lines_stream, prog_p_col) = max(n_coo+1, iob_p_idx(s));
                end
                % ----------- Separate antenna phase centers ------------------
                if n_apc > 0
                    if obs_set.obs_code(s,1) ~= u_sys_c(1)
                        A(lines_stream, prog_p_col + ( 1 : 3)) = - los_stream;
                    end
                    A_idx(lines_stream,prog_p_col + ( 1 : 3)) = repmat(apc_p_idx(s,:),length(lines_stream),1);
                    prog_p_col = prog_p_col +3;
                end
                % ----------- Abiguity ------------------
                if phase_present
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = obs_set.wl(s)*this.state.getCycleSlipThr();
                    A_idx(lines_stream, prog_p_col) = n_coo + n_iob + n_apc + amb_idx(id_ok_stream, s);
                end
                % ----------- Clock ------------------
                prog_p_col = prog_p_col + 1;
                A(lines_stream, prog_p_col) = 1;
                A_idx(lines_stream, prog_p_col) = n_coo + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                % ----------- ZTD ------------------
                if tropo
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = mfw_stream;
                    A_idx(lines_stream, prog_p_col) = n_coo + n_clocks_rec + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                end
                % ----------- ZTD gradients ------------------
                if tropo_g
                    %cotan_term = cot(el_stream) .* mfw_stream;
                    cotan_term = 1 ./ ( sin(el_stream).*tan(el_stream) + 0.0032);
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = cos(az_stream) .* cotan_term; % noth gradient
                    A_idx(lines_stream, prog_p_col) = n_coo + 2 * n_clocks_rec + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                    prog_p_col = prog_p_col + 1;
                    A(lines_stream, prog_p_col) = sin(az_stream) .* cotan_term; % east gradient
                    A_idx(lines_stream, prog_p_col) = n_coo + 3 * n_clocks_rec + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                end
                obs_count = obs_count + n_obs_stream;
            end
            % ---- Suppress weighting until solution is more stable/tested
            %w(:) = 1;%0.005;%this.state.std_phase;
            %---------------------
            sat = sat;
            if dynamic
                p_flag = [1, 1, 1, -ones(iob_flag), -repmat(ones(apc_flag),1,3), -ones(amb_flag), 1, ones(tropo), ones(tropo_g), ones(tropo_g) , ones(global_sol)];
            else
                p_flag = [zeros(1,n_coo_par) -ones(iob_flag), -repmat(ones(apc_flag),1,3), -ones(amb_flag), 1, ones(tropo), ones(tropo_g), ones(tropo_g), ones(global_sol) ];
            end
            p_class = [this.PAR_X*ones(~rec.isFixed()) , this.PAR_Y*ones(~rec.isFixed()), this.PAR_Z*ones(~rec.isFixed()), this.PAR_ISB * ones(iob_flag), this.PAR_PCO_X * ones(apc_flag),...
                this.PAR_PCO_Y * ones(apc_flag), this.PAR_PCO_Z * ones(apc_flag), this.PAR_AMB*ones(amb_flag), this.PAR_REC_CLK, this.PAR_TROPO*ones(tropo), this.PAR_TROPO_N*ones(tropo_g), ...
                this.PAR_TROPO_E*ones(tropo_g), this.PAR_SAT_CLK*ones(global_sol)];
            % find the ambiguity set jmp
            amb_set_jmp = find(sum(diff(nan2zero(amb_idx)) < 0, 2) == sum(~isnan(amb_idx(1 : end - 1, :)),2) | sum(diff(nan2zero(amb_idx)) > 0,2) == sum(~isnan(amb_idx(2 : end, :)),2)) + 1;
        end
        
        function setTimeRegularization(this, param_class, time_variability)
            idx_param = this.time_regularization == param_class;
            if sum(idx_param) > 0
                this.time_regularization(idx_param, 2) = time_variability;
            else %if not prestn add it
                this.time_regularization = [this.time_regularization; [param_class, time_variability]];
            end
        end
        
        function setMeanRegularization(this, param_class, var)
            idx_param = this.time_regularization == param_class;
            if sum(idx_param) > 0
                this.mean_regularization(idx_param, 2) = var;
            else %if not present add it
                this.mean_regularization = [this.mean_regularization; [param_class, var]];
            end
        end
        
        function Astack2Nstack(this)
            %DESCRIPTION: generate N stack A'*A
            n_obs = size(this.A_ep, 1);
            this.N_ep = zeros(size(this.A_ep, 2), size(this.A_ep, 2), n_obs);
            if isempty(this.rw)
                this.rw = ones(size(this.variance));
            end
            for i = 1 : n_obs
                A_l = this.A_ep(i, :);
                
                w = 1 / this.variance(i) * this.rw(i);
                this.N_ep(:, :, i) = A_l' * w * A_l;
            end
        end
        
        function res = getResiduals(this, x)
            
            %res_l = zeros(size(this.y));
            %for o = 1 : size(this.A_ep, 1)
            %    res_l(o) = this.y(o) - this.A_ep(o, :) * x(this.A_idx(o, :), 1);
            %end
            %res_l = zeros(size(this.y));
            % speed-up of the previous lines
            if any(isnan(x))
                this.log.addError('Some parameters are NaN!');
            end
            res_l = this.y - sum(this.A_ep .* reshape(x(this.A_idx), size(this.A_idx,1), size(this.A_idx,2)),2);
            
            this.res = res_l;
            res_l(this.rw == 0) = 0;
            n_epochs = max(this.true_epoch);
            n_sat = this.cc.getNumSat;
            n_rec = max(this.receiver_id);
            if ~this.network_solution
                res = zeros(n_epochs, n_sat);
                for s = 1 : length(this.sat_go_id)
                    idx = this.sat == s;
                    ep = this.epoch(idx);
                    res(this.true_epoch(ep), this.sat_go_id(s)) = res_l(idx);
                end
            else
                res = nan(n_epochs, n_sat, n_rec);
                idx_tot = [];
                for r = 1 : n_rec
                    for s = 1 : n_sat
                        idx = this.sat == s & this.receiver_id == r;
                        ep = this.epoch(idx);
                        res(this.true_epoch(ep), s, r) = res_l(idx);
                        idx_tot = [idx_tot; find(idx)];
                    end
                end
                av_res = sum(res, 3, 'omitnan') ./  sum(~isnan(zero2nan(res)),3);
                res = res - repmat(av_res,1,1,n_rec);
                this.res(idx_tot) = res(~isnan(res));
                res = nan2zero(res);
            end
        end
        %-----------------------------------------------
        % Implemenation of M-estimators
        % Note: after reweighting the function Astackt2Nstack have to be
        % called again
        %----------------------------------------------------------------
        function weightOnResidual(this, wfun, thr, thr_propagate)
            if isempty(this.rw)
                this.rw = ones(size(this.variance));
            end
            s0 = mean(abs(this.res).*this.rw);
            res_n = this.res/s0;
            if nargin > 2
                if nargin > 3 && (thr_propagate > 0)
                    sat_err = nan(this.n_epochs, max(this.sat_go_id));
                    sat_err(this.epoch + this.sat * this.n_epochs) = this.res/s0;
                    ssat_err = Receiver_Commons.smoothSatData([],[],sat_err, [], 'spline', 30, 10); % smoothing SNR => to be improved
                    idx_ko = false(this.n_epochs, max(this.sat_go_id));
                    for s = 1 : size(idx_ko, 2)
                        idx_ko(:,s) = (movmax(abs(ssat_err(:,s)), 20) > thr_propagate) & flagExpand(abs(ssat_err(:,s)) > thr, 100);
                    end
                    idx_rw = idx_ko(this.epoch + this.sat * this.n_epochs);
                else
                    idx_rw = abs(res_n) > thr;
                end
            else
                idx_rw = true(size(res_n));
            end
            this.rw(idx_rw) =  wfun(res_n(idx_rw));
        end
        
        function reweightHuber(this)
            threshold = 2;
            wfun = @(x) threshold ./ abs(x);
            this.weightOnResidual(wfun, threshold);
        end
        
        function reweightDanish(this)
            threshold = 2;
            wfun = @(x)  exp(-x.^2 ./threshold.^2);
            this.weightOnResidual(wfun, threshold);
        end
        
        function reweightDanishWM(this)
            threshold = 2;
            wfun = @(x)  max(0.5,exp(-x.^2 ./threshold.^2));
            this.weightOnResidual(wfun, threshold);
        end
        
        function reweightHubNoThr(this)
            wfun = @(x) 1 ./ abs(x);
            this.weightOnResidual(wfun);
        end
        
        function reweightTukey(this)
            threshold = 2;
            wfun = @(x) (1 - (x ./threshold).^2).^2;
            this.weightOnResidual(wfun, threshold);
        end
        
        function snooping(this)
            threshold = 2.5;
            wfun = @(x) 0;
            this.weightOnResidual(wfun, threshold);
        end
        
        function snoopingGatt(this, thr)
            if nargin == 1
                thr = 10;
            end
            threshold_propagate = 2.5;
            wfun = @(x) 0;
            this.weightOnResidual(wfun, thr, threshold_propagate);
        end
        
        %------------------------------------------------------------------------
        function [x, res, s0, Cxx] = solve(this)
            % if N_ep if empty call A
            if isempty(this.N_ep)
                this.Astack2Nstack();
            end
            is_network = this.network_solution;
            u_r = unique(this.receiver_id);
            n_rec = length(u_r);
            idx_rec_common_l = this.param_flag == 2;
            N = [];
            B = [];
            A2N_idx_tot = [];
            for r = u_r(:)'
                idx_r_l = this.receiver_id == u_r(r);
                idx_r = find(idx_r_l);
                A_rec = this.A_idx(idx_r_l, ~idx_rec_common_l);
                idx_constant_l = this.param_flag(~idx_rec_common_l) == 0 | this.param_flag(~idx_rec_common_l) == -1;
                idx_constant = find(idx_constant_l);
                idx_non_constant = find(~idx_constant_l);
                idx_ep_wise = find(this.param_flag == 1);
                
                a_idx_const =unique(A_rec(:, idx_constant_l));
                a_idx_const(a_idx_const == 0) = [];
                a_idx_ep_wise = unique(A_rec(:, idx_ep_wise));
                a_idx_ep_wise(a_idx_ep_wise == 0) = [];
                n_constant = length(a_idx_const);
                n_class = size(this.A_ep, 2);
                n_ep_wise = length(a_idx_ep_wise);
                if isempty(n_ep_wise)
                    n_ep_wise = 0;
                end
                
                n_epochs = this.n_epochs(r);
                if  is_network
                    u_ep = unique(this.epoch(this.receiver_id == r)); % <- consider not using unique
                    r2p_ep = zeros(u_ep(end),1); % receiver to progressive epoch
                    r2p_ep(u_ep) = 1: length(u_ep);
                else
                    r2p_ep = 1:n_epochs;
                    u_ep = r2p_ep;
                end
                n_obs = size(this.A_ep, 1);
                n_ep_class = n_ep_wise / n_epochs;
                Ncc = zeros(n_constant, n_constant);
                Nce = zeros(n_ep_wise , n_constant);
                n_class_ep_wise = length(idx_non_constant);
                Ndiags = zeros(n_class_ep_wise, n_class_ep_wise, n_epochs); %permute(this.N_ep(~idx_constant_l,~idx_constant_l,:),[3,1,2]);
                B_rec = zeros(n_constant+n_ep_wise, 1);
                if isempty(this.rw)
                    this.rw = ones(size(this.variance));
                end
                %%% all costant parameters are put before in the normal matrix find the mapping between A_idx and idx in the Normal matrix
                N2A_idx = [a_idx_const; a_idx_ep_wise];
                A2N_idx = zeros(size(N2A_idx));
                A2N_idx(N2A_idx) = 1:(n_constant + n_ep_wise);
                
                for i = idx_r'
                    p_idx = this.A_idx(i, ~idx_rec_common_l);
                    p_idx(p_idx == 0) = 1;  % does not matter since terms are zeros
                    N_ep = this.N_ep(~idx_rec_common_l,~idx_rec_common_l, i);
                    A_ep = this.A_ep(i, ~idx_rec_common_l);
                    variance = this.variance(i);
                    rw = this.rw(i);
                    y = this.y(i);
                    e = r2p_ep(this.epoch(i));
                    p_c_idx = p_idx(idx_constant_l);
                    p_e_idx = p_idx(~idx_constant_l);
                    p_e_idx(p_e_idx <= 0) = 1;  % does not matter since terms are zeros
                    p_c_idx = A2N_idx(p_c_idx);
                    p_e_idx = A2N_idx(p_e_idx) - n_constant;
                    p_idx = A2N_idx(p_idx);
                    
                    % fill Ncc
                    Ncc(p_c_idx, p_c_idx) = Ncc(p_c_idx, p_c_idx) + N_ep(idx_constant, idx_constant);
                    % fill Nce
                    Nce(p_e_idx, p_c_idx) = Nce(p_e_idx, p_c_idx) + N_ep(idx_non_constant, idx_constant);
                    %fill Ndiags
                    
                    Ndiags(:, :, e) = Ndiags(:, :, e) + N_ep(idx_non_constant, idx_non_constant);
                    %fill B
                    B_rec(p_idx) = B_rec(p_idx) + A_ep' * (1 ./ variance) * rw * y;
                end
                Nee = [];
                class_ep_wise = this.param_class(idx_non_constant);
                
                diff_reg = 1./double(diff(this.true_epoch(u_ep)));
                reg_diag0 = [diff_reg; 0 ] + [0; diff_reg];
                reg_diag1 = -diff_reg ;
                Ndiags = permute(Ndiags, [3, 1, 2]);
                tik_reg = ones(n_epochs,1)/n_epochs; %%% TIkhonov on ZTD and gradients
                % build the epoch wise parameters
                for i = 1:n_ep_class
                    N_col = [];
                    for j = 1:n_ep_class
                        diag0 = Ndiags(:, i, j);
                        N_el = sparse(n_epochs, n_epochs);
                        if j == i
                            cur_class = class_ep_wise(i);
                            % Time Regularization
                            if ~isempty(this.time_regularization)
                                idx_c = this.time_regularization(:, 1) == cur_class;
                                w = 1 ./ this.time_regularization(idx_c, 2) ;
                                if sum(idx_c)
                                    diag0 = diag0 + reg_diag0 * w;
                                    diag1 = reg_diag1 * w;
                                    N_el = spdiags([0; diag1], 1, N_el);
                                    N_el = spdiags(diag1, -1, N_el);
                                end
                            end
                            % Mean zero regularization - same as tikhonov
                            if ~isempty(this.mean_regularization)
                                idx_t = this.mean_regularization(:, 1) == cur_class;
                                if sum(idx_t)
                                    w = 1 ./ this.mean_regularization(idx_t, 2) ;
                                    diag0 = diag0 + tik_reg * w;
                                end
                            end
                            
                        end
                        N_el = spdiags(diag0, 0, N_el);
                        N_col = [N_col; N_el];
                    end
                    Nee = [Nee, N_col];
                end
                N_rec = [[Ncc, Nce']; [Nce, Nee]];
                d_N = size(N,1);
                d_N_rec = size(N_rec,1);
                N = [ [N sparse(d_N, d_N_rec)]; [sparse(d_N_rec, d_N)  N_rec]];
                B = [B; B_rec];
                if r > 1
                    A2N_idx = A2N_idx +max(A2N_idx_tot);
                end
                A2N_idx_tot = [A2N_idx_tot; A2N_idx];
            end
            
            if is_network
               
                n_obs = size(this.A_idx,1);
                % create the part of the normal that considers common parameters
                
                
                a_idx_const = unique(this.A_idx(this.receiver_id == 1, idx_constant_l));
                a_idx_const(a_idx_const == 0) = [];
                a_idx_ep_wise = unique(this.A_idx(this.receiver_id == 1, idx_ep_wise));
                a_idx_ep_wise(a_idx_ep_wise == 0) = [];
                
                N2A_idx = [ a_idx_const; a_idx_ep_wise];
                %mix the receiver indexes
                par_rec_id = ones(max(max(this.A_idx(this.receiver_id == 1,:))),1);
                A_idx_not_mix = this.A_idx;
                for j = 2 : n_rec
                    rec_idx = this.receiver_id == j;
                    % update the indexes
                    par_rec_id = [par_rec_id ; j*ones(max(max(this.A_idx(this.receiver_id == j,:))),1)];
                    this.A_idx(rec_idx,:) = this.A_idx(this.receiver_id == j, :) + max(max(this.A_idx(this.receiver_id == j-1,:)));
                    
                    a_idx_const =unique(this.A_idx(rec_idx, idx_constant_l));
                    a_idx_const(a_idx_const == 0) = [];
                    a_idx_ep_wise = unique(this.A_idx(rec_idx, idx_ep_wise));
                    a_idx_ep_wise(a_idx_ep_wise == 0) = [];
                    
                    N2A_idx = [N2A_idx; a_idx_const; a_idx_ep_wise];
                end
                % get the oidx for the common parameters
                common_idx = zeros(n_obs,1);
                u_sat = unique(this.sat); % <- very lazy, once it work it has to be oprimized
                p_idx = 0;
                for i = 1 :length(u_sat)
                    sat_idx = u_sat(i) == this.sat;
                    ep_sat = this.epoch(sat_idx);
                    u_ep = unique(ep_sat);
                    [~, idx_ep] = ismember(ep_sat, u_ep);
                    common_idx(sat_idx) = idx_ep + p_idx;
                    p_idx = p_idx + max(idx_ep);
                end
                
                n_common  = max(common_idx);
                B_comm    = zeros(n_common,1);
                diag_comm = zeros(n_common,1);
                n_s_r_p = size(this.A_idx,2); % number single receiver parameters
                N_stack_comm = zeros(n_common, n_rec * n_s_r_p);
                N_stack_idx  = zeros(n_common, n_rec * n_s_r_p);
                
                for i = 1 : n_obs
                    idx_common = common_idx(i);
                    variance = this.variance(i);
                    rw = this.rw(i);
                    y = this.y(i);
                    r = this.receiver_id(i);
                    B_comm(idx_common) = B_comm(idx_common) + rw * (1 ./ variance) * y; % the entrace in the idx for the common parameters are all one
                    N_stack_comm(idx_common, (r - 1) * n_s_r_p + (1 : n_s_r_p) )= N_stack_comm(idx_common,(r - 1) * n_s_r_p +(1 : n_s_r_p)) +this.A_ep(i,:) * rw * (1 ./ variance);
                    diag_comm(idx_common) = diag_comm(idx_common) + rw * (1 ./ variance);
                    N_stack_idx(idx_common, (r - 1) * n_s_r_p + (1 : n_s_r_p))= this.A_idx(i, :);
                end
                n_par = max(max(this.A_idx));
                 x_tot = zeros(n_par,1);
                line_idx = repmat([1:n_common]',1,size(N_stack_comm,2));
                idx_full = N_stack_idx~=0;
                Nreccom = sparse(line_idx(idx_full), N_stack_idx(idx_full), N_stack_comm(idx_full), n_common,n_par);
                % add the rank deficency considtion
                %                 1) sum of the satellites clock == 0
                %                 Nreccom = [Nreccom ones(n_common,1)];
                %                 B = [B; 0];
                %                 N = [[N ; zeros(1,n_par)] zeros(n_par+1,1)];
                i_diag_comm = 1 ./ diag_comm;
                i_diag_comm = spdiags(i_diag_comm,0, n_common,n_common);
                rNcomm = Nreccom'*i_diag_comm;
                
                N = N - rNcomm*Nreccom; % N11r = N11 - N21 * inv(N22) * N12
                
                B = B - rNcomm*B_comm;  % B1r = N1 - N21 * inv(N22) * B2
                                
                % resolve the rank deficency
                % ALL paramters has a rank deficency beacause the entrance of the image matrixes are very similar and we also estimated the clock of the satellite
                % 2) remove coordinates and tropo paramters of the first receiver
                % we can do that because tropo paramters are slightly constarined in time so evan if they are non present for the first receiver the rank deficecny is avoided
                idx_rec_x = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_X));
                idx_rec_y = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_Y));
                idx_rec_z = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_Z));
                idx_rec_isb = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_ISB));
                idx_rec_t = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_TROPO));
                idx_rec_tn = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_TROPO_N));
                idx_rec_te = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_TROPO_E));
                idx_rm = [idx_rec_x; idx_rec_y; idx_rec_z; idx_rec_isb; idx_rec_t; idx_rec_tn; idx_rec_te];
                % 3 ) remove one clock per epoch for the minim receiver available
                clk_idx = this.param_class == this.PAR_REC_CLK;
                n_epochs = length(unique(this.epoch));
                idx_clk_to_rm = true(n_epochs,1);
                i = 1;
                while sum(idx_clk_to_rm) > 0
                     idx_i_c = this.receiver_id == i;
                     idx_rec_clk = unique(this.A_idx(idx_i_c, clk_idx)); 
                     idx_rec_clk_ep =  unique(this.epoch(idx_i_c));
                     idx_to_rm = idx_clk_to_rm(idx_rec_clk_ep);
                     idx_rm = [idx_rm ; idx_rec_clk(idx_to_rm)];
                     idx_clk_to_rm(idx_rec_clk_ep(idx_to_rm)) = false;
                     i = i+1;
                end
                
                idx_amb_rm = [];
                prev_info = ~isempty(this.apriori_info);
                % remove the ambiguity that are not connected
                if prev_info && this.state.getCurSession > 1%%% introduce the previous amniguity
                    conn_amb  = false(size(this.apriori_info.amb_value));
                    for i = 1:(length(this.apriori_info.amb_value))
                        % determine the mapping to new freq
                        r_id = this.apriori_info.receiver(i);
                        s_id = this.apriori_info.goids(i);
                        idx_ambs = this.receiver_id == r_id & this.sat == s_id;
                        if sum(this.epoch(idx_ambs) == 1) > 0
                            conn_amb(i) = true;
                        end
                    end
                    this.removeAprInfo(~conn_amb);
                end
                
                
                
                % 4)remove one ambiguity per satellite form the firs receiver
                if true
                    %n_jmp_sat 
                    for i = 1 :length(u_sat)
                        jmp_idx = find(diff(this.sat_jmp_idx(:,u_sat(i))) == -1) + 1;
                        if ~this.sat_jmp_idx(1,u_sat(i))
                            if ~prev_info || sum(u_sat(i) == this.apriori_info.goids) == 0
                                jmp_idx = [1; jmp_idx];
                            end
                        end
                        for j = 1: length(jmp_idx(:))
                            idx_amb_rec = [];
                            d = 1;
                            jmp = jmp_idx(j);
                            jmp2 = jmp_idx(min(j+1, length(jmp_idx)));
                            if jmp == jmp2
                                jmp2 = Inf;
                            end
                            %                         while isempty(idx_amb_rec) && d <= n_rec
                            %                            idx_amb_rec = this.A_idx(this.receiver_id == d & this.sat == u_sat(i) & this.epoch >= jmp & this.epoch < jmp2 ,this.param_class == this.PAR_AMB);
                            %                            d = d + 1;
                            %                         end
                            idx_amb_rec = this.A_idx(this.sat == u_sat(i) & this.epoch >= jmp & this.epoch < jmp2 ,this.param_class == this.PAR_AMB);
                            idx_amb_rec = Core_Utils.remBFromA(idx_amb_rec,idx_amb_rm);
                            if ~isempty(idx_amb_rec)
                                idx_amb_rec = mode(idx_amb_rec);%(1);%idx_amb_rec(min(120,length(idx_amb_rec)));
                            end
                            idx_amb_rm = [idx_amb_rm; idx_amb_rec];
                        end
                    end
                    %                 5) remove one ambiguity per each set of disjunt set of arcs of each receiver to resolve the ambiguity-receiver clock rank deficency
                    %                 first recievr does not have clocks any more so no rank defricency
                    for i = 2 : n_rec
                        jmps = this.amb_set_jmp{i}';
                        if ~prev_info
                            jmps = [1 jmps];
                        end
                        for j = 1 : length(jmps)
                            jmp = jmps(j);
                            jmp2 = jmps(min(j+1,length(jmps)));
                            if jmp == jmp2
                                jmp2 = Inf;
                            end
                            idx_amb_rec = this.A_idx(this.receiver_id == i & this.epoch >= jmp & this.epoch < jmp2,this.param_class == this.PAR_AMB);
                            g = 1;
                            %                         while sum(idx_amb_rec(g) == idx_amb_rm) > 0 && g < length(idx_amb_rec)
                            %                             g = g +1;
                            %                         end
                            idx_amb_rec = Core_Utils.remBFromA(idx_amb_rec,idx_amb_rm);
                            if ~isempty(idx_amb_rec)
                                
                                idx_amb_rec = mode(idx_amb_rec);%(1);%idx_amb_rec(min(120,length(idx_amb_rec)));
                            end
                            idx_amb_rm = [idx_amb_rm; idx_amb_rec];
                        end
                    end
                end
                idx_rm = [idx_rm ; idx_amb_rm];
                idx_rm = unique(idx_rm);
                N(idx_rm, :) = [];
                N(:, idx_rm) = [];
                B(idx_rm) = [];
                if prev_info && this.state.getCurSession > 1%%% introduce the previous amniguity
                    par_ids = zeros(size(this.apriori_info.amb_value));
                    valid_float = false(size(this.apriori_info.amb_value));
                    for i = 1:(length(this.apriori_info.amb_value))
                        % determine the mapping to new freq
                        r_id = this.apriori_info.receiver(i);
                        s_id = this.apriori_info.goids(i);
                        idx_ambs = this.receiver_id == r_id & this.sat == s_id;
                        if sum(this.epoch(idx_ambs) == 1) > 0 % if there is the ambiguity
                            par_id = this.A_idx(find(idx_ambs,1,'first'),this.param_class == this.PAR_AMB);
                            if sum(par_id == idx_rm) ==0
                                par_ids(i) = par_id;
                                if this.apriori_info.fixed(i)
                                    % put the amniguity in the result
                                    x_tot(par_id) = this.apriori_info.amb_value(i);
                                    % remove paramter from the normal matrix
                                    par_id2 = par_id - sum(idx_rm < par_id);
                                    Ni = N(par_id2,:);
                                    B = B -( Ni*this.apriori_info.amb_value(i))';
                                    N(par_id2,:) = [];
                                    N(:,par_id2) = [];
                                    B(par_id2) = [];
                                    % ------
                                    idx_rm = [idx_rm; par_id];
                                else
                                    valid_float(i) = true;
                                end
                            end
                        end
                        
                    end
                    % add the float with their vcv
                    par_ids_float = par_ids(valid_float);
                    if length(par_ids_float) > 0
                        for i = 1 : length(par_ids_float) % remove removed id
                            par_ids_float(i) = par_ids_float(i) - sum(idx_rm < par_ids_float(i));
                        end
                        idx_v_f = valid_float(this.apriori_info.fixed == 0);
                        Cambamb = this.apriori_info.Cambamb(idx_v_f,idx_v_f);
                        Napri = inv(Cambamb);
                        B(par_ids_float) = B(par_ids_float) + Napri* this.apriori_info.amb_value(valid_float);
                        N(par_ids_float,par_ids_float) = N(par_ids_float,par_ids_float) + Napri;
                    end
                end
            end
            if ~isempty(this.G)
                if ~is_network
                    G = this.G(:,N2A_idx);
                else
                    G = this.G;
                end
                N =  [[N, G']; [G, zeros(size(G,1))]];
                B = [B; this.D];
            end
            
            
            x = N \ B;
            
            x_class = zeros(size(x));
            for c = 1:length(this.param_class)
                idx_p = A2N_idx_tot(this.A_idx(:, c));
                x_class(idx_p) = this.param_class(c);
            end
            if is_network
                idx_est = true(n_par,1);
                idx_est(idx_rm) = false;
                x_tot(idx_est) = x;
                x = x_tot;
                
                idx_amb_par = find(x_class(idx_est) == this.PAR_AMB);
                n_amb = length(idx_amb_par);
            end
                %[x, flag] =  pcg(N,B,1e-9, 10000);
                
           
            
            cxx_comp = false;
            if is_network && this.state.flag_amb_pass
                % getting tht VCV matrix for the ambiuities
                %idx_amb_par = find(x_class == this.PAR_AMB);
                this.x_float = x;
                b_eye = zeros(length(B),n_amb);
                idx = sub2ind(size(b_eye),idx_amb_par,[1:n_amb]');
                b_eye(idx) = 1;
                b_eye = sparse(b_eye);
                Cxx = N\b_eye;
                Cxx = Cxx(idx_amb_par,:);
                this.Cxx_amb = Cxx;
            end
            if (this.state.flag_amb_fix && length(x(x_class == 5,1))> 0) 
                % IMPORTANT NOTE:
                % This part on ambiguity fixing use a simple integer rounding, this is done mainly for two reason:
                % 1) For long observing time  it seem sufficient in most of the cases
                % 2) The VCV of the ambiguity derived from several hours joint multi epoch adjustmnet seem far too optimistic (well known fact). 
                %    This to our understending might cause problem in the usage of more sofisticated search method in ambiguity space such as the LAMBDA method.
                % A better understanding is required
                if ~is_network
                    amb = x(x_class == 5,1);
                    amb_wl_fixed = false(size(amb));
                    amb_n1 = nan(size(amb));
                    amb_wl = nan(size(amb));
                    n_ep_wl = zeros(size(amb));
                    n_amb = max(max(this.amb_idx));
                    n_ep = size(this.wl_amb,1);
                    n_coo = max(this.A_idx(:,3));
                    for i = 1 : n_amb
                        sat = this.sat_go_id(this.sat(this.A_idx(:,4)== i+n_coo));
                        idx = n_ep*(sat(1)-1) +  this.true_epoch(this.epoch(this.A_idx(:,4)== i+n_coo));
                        amb_wl(i) = this.wl_amb(idx(1));
                        amb_wl_fixed(i)=  this.wl_fixed(idx(1));
                        n_ep_wl(i) = length(idx);
                        amb_n1(i) = amb(i)/0.1070; %(amb(i)- 0*f_vec(2)^2*l_vec(2)/(f_vec(1)^2 - f_vec(2)^2)* wl_amb);  % Blewitt 1989 eq(23)
                        
                    end
                    
                    weight = min(n_ep_wl(amb_wl_fixed),100); % <- downweight too short arc
                    weight = weight / sum(weight);
                    [~, wl_frac] = Core_Utils.getFracBias(amb_n1(amb_wl_fixed),weight);
                    
                    amb_n1_fix = round(amb_n1 - wl_frac);
                    frac_part_n1 = amb_n1 - amb_n1_fix - wl_frac;
                    
                    idx_amb = find(x_class == 5);
                    % get thc cxx of the ambiguities
%                     n_amb  = length(idx_amb);
%                     b_eye = zeros(length(B),n_amb);
%                     idx = sub2ind(size(b_eye),idx_amb,[1:n_amb]');
%                     b_eye(idx) = 1;
%                     b_eye = sparse(b_eye);
%                     Cxx_amb = N\b_eye;
%                     Cxx_amb = Cxx_amb(idx_amb,:) / 0.1070^2;
                    
                    idx_fix = abs(frac_part_n1) < 0.20 & amb_wl_fixed & n_ep_wl > 30; % fixing criteria (very rough)
                    
                    idxFix2idxFlo = 1 : length(x);
                    idxFlo2idxFix = nan(length(x),1);
                    A_fixed = false(size(this.A_idx(:,4)));
                    
                    for i = 1 : length(idx_fix)
                        Ni = N(:,idx_amb(i));
                        if idx_fix(i)
                            b_if_fix = 0.1070 * (amb_n1_fix(i));
                            B = B - Ni* ( b_if_fix);
                            A_fixed = A_fixed | this.A_idx(:,4) == idx_amb(i);
                        end
                    end
                    
                    
                    B(idx_amb(idx_fix)) = [];
                    B(end) = 0;
                    N(idx_amb(idx_fix),:) = [];
                    N(:,idx_amb(idx_fix)) = [];
                    idxFix2idxFlo(idx_amb(idx_fix)) = [];
                    idxFlo2idxFix(idx_amb(idx_fix)) = 0;
                    idxFlo2idxFix(idxFlo2idxFix ~=0) = 1:length(B);
                    xf = N \ B;
                    % ---------------- try to fix again
                    new_idx_amb = idxFlo2idxFix(idx_amb);
                    new_idx_amb(new_idx_amb == 0) = [];
                    amb_n1 = xf(new_idx_amb)/0.1070;
                    weight = min(n_ep_wl(~idx_fix),100); % <- downweight too short arc
                    weight = weight / sum(weight);
                    [~, wl_frac] = Core_Utils.getFracBias(amb_n1,weight);
                    amb_n1_fix(~idx_fix) = round(amb_n1 - wl_frac);
                    frac_part_n1(~idx_fix) = amb_n1 - amb_n1_fix(~idx_fix) - wl_frac;
                    idx_fix_old = idx_fix;
                    idx_fix(~idx_fix) = abs(frac_part_n1(~idx_fix)) < 0.20 & n_ep_wl(~idx_fix) > 30; % fixing criteria (very rough)
                    new_fix_idx = idx_fix(~idx_fix_old);
                    new2old_idx_fix = find(~idx_fix_old);
                    for i = 1 : length(new_fix_idx)
                        Ni = N(:,new_idx_amb(i));
                        if new_fix_idx(i)
                            b_if_fix = 0.1070 * (amb_n1_fix(new2old_idx_fix(i)));
                            B = B - Ni* ( b_if_fix);
                            A_fixed = A_fixed | this.A_idx(:,4) == idx_amb(new2old_idx_fix(i));
                        end
                    end
                    B(idx_amb(new_fix_idx)) = [];
                    B(end) = 0;
                    N(idx_amb(new_fix_idx),:) = [];
                    N(:,idx_amb(new_fix_idx)) = [];
                    idxFix2idxFlo = 1 : length(x);
                    idxFlo2idxFix = nan(length(x),1);
                    idxFix2idxFlo(idx_amb(idx_fix)) = [];
                    idxFlo2idxFix(idx_amb(idx_fix)) = 0;
                    idxFlo2idxFix(idxFlo2idxFix ~=0) = 1:length(B);
                    %-------------------------------------------------------------
                    xf = N \ B;
                    x_old = x;
                    x(idxFix2idxFlo) = xf;
                    x(idx_amb(idx_fix)) = amb_n1_fix(idx_fix) * 0.1070;
                    this.log.addMessage(this.log.indent(sprintf('%d of %d ambiguity fixed\n',sum(idx_fix),length(idx_fix))));
                    this.log.addMessage(this.log.indent(sprintf('%.2f %% of observation has the ambiguity fixed\n',sum(A_fixed)/length(A_fixed)*100)));
                else
                    f = 0;
                    nf_loop = 2;
                    while f < nf_loop
                        
                        idx_est2idx = 1:size(x,1); % index to convert from the actual estimated varibles to the tototal one
                        idx_est2idx = idx_est2idx(idx_est);
                       
                        
                        xe = x(idx_est);
                        amb = xe(idx_amb_par,1);
                        
                        %end
                        n_ep_amb = zeros(size(amb));
                        for i = 1 : n_amb
                            n_ep_amb(i) = sum(this.A_idx(:,4) == idx_est2idx(idx_amb_par(i)));
                        end
                        
                        if f == 0 % keep record of fixed ambiguities to show statistics at the end
                            n_amb_ini = n_amb;
                            n_obs_tot = sum(n_ep_amb);
                            n_amb_fix = 0;
                            n_obs_fix = 0;
                        end
                        
                        idx_fix = abs(fracFNI(amb)) < 0.15 & n_ep_amb > 20;
                        n_amb_fix = n_amb_fix + sum(idx_fix); % update the record of fixed ambiguities
                        n_obs_fix = n_obs_fix + sum(n_ep_amb(idx_fix)); % update the record of fixed ambiguities
                        amb_fix = round(amb);
                        
                        for i = 1 : length(idx_fix)
                            Ni = N(:,idx_amb_par(i));
                            if idx_fix(i)
                                b_if_fix = (amb_fix(i));
                                B = B - Ni* ( b_if_fix);
                                n_amb = n_amb -1;
                            end
                        end
                        
                        B(idx_amb_par(idx_fix)) = [];
                        N(idx_amb_par(idx_fix),:) = [];
                        N(:,idx_amb_par(idx_fix)) = [];
                        
                        idx_nf = true(sum(idx_est,1),1);
                        idx_nf(idx_amb_par(idx_fix)) = false;
                        xf = zeros(size(idx_nf));
                        
                        xf(idx_nf) = N \ B;
                        
                        xf(~idx_nf) = amb_fix(idx_fix);%.*wl(idx_fix);
                        
                        x(idx_est) = xf;
                        % chenge the index of estaimted values
                        n_p = length(idx_est);
                        idx_est = find(idx_est);
                        idx_est(~idx_nf) = [];
                        idx_est = num2LogIdx(idx_est,n_p);
                        if nargout > 3 && f == (nf_loop -1)
                            % getting tht VCV matrix for the ambiuities
                            b_eye = zeros(length(B),n_amb);
                            idx_t_amb_par = find(x_class(idx_est) == this.PAR_AMB);
                            idx = sub2ind(size(b_eye),idx_t_amb_par,[1:n_amb]');
                            b_eye(idx) = 1;
                            b_eye = sparse(b_eye);
                            Cxx = N\b_eye;
                            Cxx = Cxx(idx_t_amb_par,:);
                            cxx_comp = true;
                        end
                        f = f +1;
                    end
                    
                    this.log.addMessage(this.log.indent(sprintf('%d of %d ambiguity fixed\n',n_amb_fix,n_amb_ini)));
                    this.log.addMessage(this.log.indent(sprintf('%.2f %% of observation has the ambiguity fixed\n',n_obs_fix/n_obs_tot*100)));
                end
            end
            if nargout > 1
                x_res = zeros(size(x));
                x_res(N2A_idx) = x(1:end-size(this.G,1));
                if sum(isnan(x_res)) ==0
                    res = this.getResiduals(x_res);
                    s0 = mean(abs(res(res~=0)));
                    
                else
                    res = [];
                    s0 = Inf;
                end
                
            end
            x = [x, x_class];
            % restore old Idx
            if is_network
                this.A_idx_mix = this.A_idx;
                this.A_idx = A_idx_not_mix;
                x_rec = ones(size(x,1),1);
                id_rec = find(diff(x_class) < 0);
                for i = 1 : length(id_rec)
                    idx = id_rec(i)+1;
                    x_rec(idx:end) = i+1;
                end
                x = [x, x_rec];
            end
        end
        
        function reduceNormalEquation(this, keep_param)
            % reduce number of parmeters (STUB)
            N = this.N;
            B = this.B;
            N11 = N(kp_idx,kp_idx);
            N12 = N(kp_idx,rd_idx);
            N21 = N(rd_idx,kp_idx);
            N22 = N(rd_idx,rd_idx);
            RD = N12 * inv(N11);
            this.N = N11 - RD * N21;
            B1 = B(kp_idx);
            B2 = B(rd_idx);
            this.B = B1 - RD * B2;
        end
        
        function [A_full, col_type, epoch, A_small, col_type_small, epoch_small] = getDesignMatrix(this, max_ep)
            % Get the full Design matrix for debug purposes
            %   
            % SYNTAX
            %   [A_full, col_type, A_small, col_type_small] = this.getDesignMatrix(max_ep)
            %
            tic
            A_full = sparse(size(this.A_idx,1), max(this.A_idx(:)));
            A_serial_id = serialize(repmat((1 : size(this.A_idx,1))', 1, size(this.A_idx, 2)) + (this.A_idx - 1) * size(this.A_idx, 1));
            col_type = nan(max(this.A_idx(:)), 1);
            for t = 1 : size(this.A_idx, 2)
                col_type(this.A_idx(:, t)) = t;
            end
            A_full(A_serial_id) = this.A_ep(:);
            col_type = this.param_class(col_type);
            epoch = this.epoch;
            toc
            
            if nargin == 2
                id_ok = epoch < max_ep;
                epoch_small = this.epoch(id_ok);
                col_ok = sum(abs(A_full(id_ok, :)) ~= 0);
                col_type_small = col_type(col_ok > 0);
                A_small = A_full(id_ok, col_ok > 0);                
            end
        end
        
        function pos_idx = getCommonPosIdx(this)
            % get the unique position idx fro all receiver
            pos_idx = [];
            for i = 1 : length(this.pos_indexs_tc)
                pos_idx = unique([pos_idx; this.pos_indexs_tc{i}]);
            end
            if isempty(pos_idx)
                pos_idx = 1;
            end
        end
        
        function remAmb(this, p_idx, r_id, amb_value)
            % remove an ambiguity from the system
            %
            % SYNTAX:
            % this.remAmb(p_idx, r_id, amb_value)
            
            rc_idx = this.receiver_id == r_id;
            p_amb = this.param_class == this.PAR_AMB;
            o_idx = rc_idx & this.A_idx(:,p_amb)== p_idx;
            wl = this.A(o_idx,p_amb);
            wl=wl(1);
            this.A(o_idx,p_amb) = 0;
            this.A_idx(o_idx,p_amb) = 0;
            % reduce the idx by one
            rd_idx = repmat(rec_idx,1,size(this.A_idx,2)) & this.A_id > p_idx;
            this.A_idx(rd_idx) = this.A_idx(rd_idx) - 1;
            this.obs(o_idx) = this.obs(o_idx) + amb_value * wl;
            
        end
        
         function removeAprInfo(this,idx)
            % remove ambigutiy form apriori info
            %
            % SYNTAX:
            %    this.removeAprInfo(idx)
            this.apriori_info.amb_value(idx) = [];
            %this.apriori_info.freqs(idx) = [];
            this.apriori_info.goids(idx) = [];
            this.apriori_info.receiver(idx) = [];
            idx_f = idx(~this.apriori_info.fixed);
            this.apriori_info.Cambamb(idx_f,:) = [];
            this.apriori_info.Cambamb(:,idx_f) = [];
            this.apriori_info.fixed(idx) = [];
                        
        end
    end
    
    methods (Static)
        function [resulting_gps_time, idxes, sat_jmp_idx] = intersectObsSet(obs_set_lst)
            % get the times common to at leas 2 receiver between observation set and return the index of each observation set
            %
            % SYNTAX:
            % [resulting_gps_time, idxes] = intersectObsSet(obs_set_lst)
            
            % OUTPUT:
            % resulting_gps_time = common time in gpstime (double)
            % idxe = double [n_time, n_obs_set] indx to which the commontimes correnspond in each ob_set_time
            % sat_jmp_idx = [n_time, n_sat] index of all satellite jumping
            
            % get min time, max time, and min common rate
            n_rec = length(obs_set_lst);
            min_gps_time = Inf;
            max_gps_time = 0;
            rate = zeros(n_rec,1);
            n_sat = 0;
            for i = 1 : n_rec
                % note: excluding receivers which samples are shifted w.r.t. the start of the minute e.g -> sampling rate of 30 s starting at 5s and 35s of the minute
                tmp_rate = obs_set_lst(i).time.getRate();
                if median(fracFNI(obs_set_lst(i).time.getSecond()/ tmp_rate)*tmp_rate) > 1 % if more than 50 percent of the samples are shifted of more than one second fromstart of the minute skip the receiver
                    rate(i) = -1;
                    this.log.addWarning(sprintf('Receiver %d sample are shifted from start of the minute, skipping',i));
                else
                    min_gps_time = min(min_gps_time, obs_set_lst(i).time.getNominalTime().first.getGpsTime);
                    max_gps_time = max(max_gps_time, obs_set_lst(i).time.getNominalTime().last.getGpsTime);
                    rate(i) = tmp_rate;
                    n_sat = max([n_sat; obs_set_lst(i).go_id(:)]);
                end
            end
            rate(rate < 0 ) = [];
            loop = true;
            while (loop) % find the minimum rate common to at least 2 receivers
                min_rate = min(rate);
                rate_idx = rate <= min_rate;
                if sum(rate_idx) > 1
                    loop = false;
                else
                    rate(rate_idx) = [];
                end
            end
            % for each receiver check which time of the common time are present
            resulting_gps_time = min_gps_time : min_rate : max_gps_time;
            idxes = zeros(length(resulting_gps_time), n_rec);
            % number of reciever seeing a satellite per epoch
            presence_mat = zeros(length(resulting_gps_time),n_sat);
            for i = 1 : n_rec
                [idx_is, idx_pos] = ismembertol(obs_set_lst(i).time.getNominalTime().getGpsTime(),resulting_gps_time); % tolleranc to 1 ms double check cause is already nominal rtime
                idx_pos = idx_pos(idx_pos > 0);

                obs_set_lst(i).keepEpochs(idx_is)
                pos_vec = 1 : obs_set_lst(i).time.length;
                idxes(idx_pos,i) = pos_vec;
                
                for j = 1 : n_sat
                    presence_idx = obs_set_lst(i).obs(idx_is,obs_set_lst(i).go_id == j) ~= 0;
                    presence_idx = idx_pos(presence_idx);
                    presence_mat(presence_idx,j) = presence_mat(presence_idx,j) + 1;
                end
                
            end
            
            % remove not useful observations
            for i = 1 : n_rec
                idx_rm = false(size(obs_set_lst(i).obs));
                for j = 1 : n_sat
                    % get observation observed only from one satellite
                    idx = presence_mat(:,j) == 1;
                    idx = idxes(idx,i);
                    idx(idx==0) = [];
                    idx_rm(idx,obs_set_lst(i).go_id == j) = true;
                    % get one epoch arcs
                    [flag_intervals] = getOutliers(presence_mat(:,j)>1);
                    single_arcs = (flag_intervals(:,2) - flag_intervals(:,1))  == 0;
                    idx = idxes(flag_intervals(single_arcs,1));
                    idx(idx==0) = [];
                    idx_rm(idx, obs_set_lst(i).go_id == j) = true;
                end
                obs_set_lst(i).remObs(idx_rm, false);
                %obs_set_lst(i).remShortArc();
            end
            
            
            % -- removing epoch for which no satellite is seen by at least teo receivers
            idx_rem = sum(presence_mat > 1, 2) == 0;
            idxes(idx_rem,:) = [];
            resulting_gps_time(idx_rem) = [];
            
            % --- for each satellite checks epochs for which all receiver-satellite observation continuity is broken
            sat_jmp_idx =  true(length(resulting_gps_time),n_sat);
            for i = 1 : n_sat
                for j = 1 : n_rec
                    goid_idx = obs_set_lst(j).go_id == i;
                    for k = find(goid_idx)'
                        idx_rec = obs_set_lst(j).obs(:,k) == 0 | obs_set_lst(j).cycle_slip(:,k);
                        [idx_is, idx_pos] = ismembertol(obs_set_lst(j).time.getNominalTime().getGpsTime(), resulting_gps_time); % tolleranc to 1 ms double check cause is already nominal rtime
                        idx_pos = idx_pos(idx_pos > 0);
                        sat_jmp_idx(idx_pos(idx_is),i) =  sat_jmp_idx(idx_pos(idx_is),i) & idx_rec(idx_is);
                    end
                end
            end
        end
        
        function [pos_idx_nh, pos_idx_tc] = getPosIdx(time, st_time, coo_rate)
            % given a time and the sampling rate return the position index referring to the given sampling rate, the first index is porgressive, the seocond id time consistent
            sec_from_sod = time.getRefTime(st_time.getMatlabTime);
            pos_idx_tc = max(1,ceil((sec_from_sod - 0.002) / coo_rate));
            u_pos = unique(pos_idx_tc);
            pos_idx_nh = pos_idx_tc;
            for i = 1 : length(u_pos)
                pos_idx_nh(pos_idx_tc==u_pos(i)) = i;
            end
            pos_idx_tc = unique(pos_idx_tc);
        end
    end
end
