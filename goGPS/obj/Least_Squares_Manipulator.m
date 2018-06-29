%   CLASS Least_Square_Manipulator
% =========================================================================
%
% DESCRIPTION
%   Efficently manipulate sparse least squares system
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
%    |___/                    v 0.6.0 alpha 3 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
classdef Least_Squares_Manipulator < handle
    
    properties (Constant)
        PAR_X = 1;
        PAR_Y = 2;
        PAR_Z = 3;
        PAR_ISB = 4;
        PAR_AMB = 5;
        PAR_CLK = 6;
        PAR_TROPO = 7;
        PAR_TROPO_N = 8;
        PAR_TROPO_E = 9;
    end
    
    properties
        A_ep % Stacked epochwise design matrices [n_obs x n_param_per_epoch]
        A_idx % index of the paramter [n_obs x n_param_per_epoch]
        amb_idx % index of the columns per satellite
        go_id_amb % go ids of the amb idx
        out_idx % index to tell if observation is outlier [ n_obs x 1]
        N_ep  % Stacked epochwise normal matrices [ n_param_per_epoch x n_param_per_epoch x n_obs]
        G % hard constraints (Lagrange multiplier)
        D % known term of the hard constraints
        y % observations  [ n_obs x 1]
        system_split % limits of the ambiguity splits
        variance % observation variance [ n_obs x 1]
        rw % reweight factor
        res % observations residuals
        epoch % epoch of the obseravtions and of the A lines [ n_obs x 1]
        sat % satellite of the obseravtions and of the A lines [ n_obs x 1]
        n_epochs
        param_class % [n_param_per_epoch x 1] each paramter can be part of a class
        %   [ 1 : x
        %     2 : y
        %     3 : z
        %     4 : inter channel/frequency/system biases,
        %     5 : ambiguity,
        %     6 : clock
        %     7 : tropo
        %     8 : tropo inclination north
        %     9 : tropo inclination east ]
        param_flag % 0: constant in time always same param, -1: constant in time differents param (e.g ambiguity), 1: same param changing epochwise
        time_regularization %[ param_class time_varability] simple time regularization constructed from psudo obs p_ep+1 - p_ep = 0 with given accuracy
        mean_regularization
        true_epoch % true epoch of the epochwise paramters
        sat_go_id  % go id of the satindeces
        receiver_id % id of the receiver , in case of differenced obsercation two columns are used
        
        wl_amb    % widelane ambuiguity
        wl_fixed  % is widelane fixed
    end
    
    properties(Access = private)
        Ncc % part of the normal matrix with costant paramters
        Nee % diagonal part of the normal matrix with epoch wise or multi epoch wise paramters
        Nce % cross element between constant and epoch varying paramters
        state % current settings
        rf % reference frame
        
        log % logger
    end
    
    methods
        function this = Least_Squares_Manipulator()
            this.state = Global_Configuration.getCurrentSettings();
            this.rf = Core_Reference_Frame.getInstance();
            this.log = Logger.getInstance();
        end
        
        function id_sync = setUpPPP(this, rec, id_sync,  cut_off, dynamic, pos_idx)
            if nargin < 4
                cut_off = [];
            end
            if nargin < 5
                dynamic = false;
            end
            if nargin < 6
                pos_idx = [];
            end
            id_sync = this.setUpSA(rec, id_sync, 'L', cut_off,'',dynamic,pos_idx);
        end
        
        function id_sync = setUpCodeSatic(this, rec, id_sync, cut_off)
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
            % return the id_sync of the epochs to be computed
            % get double frequency iono_free for all the systems
            % INPUT:
            %    rec : receiver
            %    id_sync : epoch to be used
            %    obs_type : 'C' 'L' 'CL'
            %    cut_off : cut off angle [optional]
            %
            if nargin < 8
                pos_idx_vec = [];
            end
            if nargin < 7
                dynamic = false;
            end
            % extract the observations to be used for the solution
            phase_present = instr(obs_type, 'L');
            flag_amb_fix = this.state.flag_amb_fix;
            if nargin < 6 || isempty(custom_obs_set)
                obs_set = Observation_Set();
                if rec.isMultiFreq() && ~rec.state.isIonoExtModel %% case multi frequency
                    if not(flag_amb_fix) || ~phase_present
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
                    else
                         [obs_set, this.wl_amb, this.wl_fixed]  = rec.getIonoFreeWidelaneFixed();
                    end
                else
                    for sys_c = rec.cc.sys_c
                        f = rec.getFreqs(sys_c);
                        for i = 1 : length(obs_type)
                            obs_set.merge(rec.getPrefObsSetCh([obs_type(i) num2str(f(1))], sys_c));
                        end
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
                    out_idx = idxCharLines(obs_code_ph,obs_code1) & go_id_ph == g;
                    out = rec.sat.outlier_idx_ph(:,out_idx);
                    if strcmp(obs_code2,'   ') 
                         out_idx = idxCharLines(obs_code_ph,obs_code2) & go_id_ph == g;
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
            
            this.true_epoch = obs_set.getTimeIdx(rec.time.first, rec.getRate); % link between original epoch, and epochs used here
            
            % remove not valid empty epoch or with only one satellite (probably too bad conditioned)
            idx_valid_ep_l = sum(diff_obs ~= 0, 2) > 0;
            diff_obs(~idx_valid_ep_l, :) = [];
            xs_loc(~idx_valid_ep_l, :, :) = [];
            
            %id_sync_out(~idx_valid_ep_l) = [];
            
            this.true_epoch(~idx_valid_ep_l) = [];
            id_sync_out = this.true_epoch;
            
            % removing possible empty column (sat)
            idx_valid_stream = sum(diff_obs, 1) ~= 0;
            diff_obs(:, ~idx_valid_stream) = [];
            xs_loc(:, ~idx_valid_stream, :) = [];
            
            % removing non valid epochs also from obs_set
            obs_set.remEpochs(~idx_valid_ep_l);
            obs_set.sanitizeEmpty();
            
            % set up requested number of parametrs
            n_epochs = size(obs_set.obs, 1);
            this.n_epochs = n_epochs;
            n_stream = size(diff_obs, 2); % number of satellites
            
            n_clocks = n_epochs; % number of clock errors
            n_tropo = n_clocks; % number of epoch for ZTD estimation
            ep_p_idx = 1 : n_clocks; % indexes of epochs starting from 1 to n_epochs
            
            % Compute the number of ambiguities that must be computed
            cycle_slip = obs_set.cycle_slip;
            cycle_slip(diff_obs == 0) = 0;
            if phase_present
                amb_idx = ones(size(cycle_slip));
                for s = 1:n_stream
                    if s > 1
                        amb_idx(:, s) = amb_idx(:, s) + amb_idx(n_epochs, s-1);
                    end
                    cs = find(cycle_slip(:, s) > 0)';
                    for c = cs
                        % check if cycle slip is not marked at first epoch of the stream
                        if c ~= find(diff_obs(:, s) ~= 0, 1, 'first')
                            amb_idx(c:end, s) = amb_idx(c:end, s) + 1;
                        end
                    end
                end
                % amb_idx = n_coo + n_iob + amb_idx;
                amb_idx = zero2nan(amb_idx .* (diff_obs ~= 0));
                
                % remove short arcs
                min_arc = this.state.getMinArc;
                % ambiguity number for each satellite
                amb_obs_count = histcounts(serialize(amb_idx), 'Normalization', 'count', 'BinMethod', 'integers');
                assert(numel(amb_obs_count) == max(amb_idx(:))); % This should always be true
                id = 1 : numel(amb_obs_count);
                ko_amb_list = id(amb_obs_count < min_arc);
                amb_obs_count(amb_obs_count < min_arc) = [];
                for ko_amb = fliplr(ko_amb_list)
                    id_ko = amb_idx == ko_amb;
                    diff_obs(id_ko) = 0;
                    amb_idx(id_ko) = nan;
                    amb_idx(amb_idx > ko_amb) = amb_idx(amb_idx > ko_amb) - 1;
                end
                obs_set.remObs(diff_obs == 0);
                
                % I need to refilter and recompute some things...
                % remove not valid empty epoch or with only one satellite (probably too bad conditioned)
                idx_valid_ep_l = sum(diff_obs ~= 0, 2) > 0;
                diff_obs(~idx_valid_ep_l, :) = [];
                xs_loc(~idx_valid_ep_l, :, :) = [];
                id_sync_out(~idx_valid_ep_l) = [];
                amb_idx(~idx_valid_ep_l, :) = [];
                
                this.true_epoch(~idx_valid_ep_l) = [];
                
                % removing possible empty column
                idx_valid_stream = sum(diff_obs, 1) ~= 0;
                diff_obs(:, ~idx_valid_stream) = [];
                xs_loc(:, ~idx_valid_stream, :) = [];
                amb_idx(:, ~idx_valid_stream, :) = [];
                
                % removing non valid epochs also from obs_set
                obs_set.sanitizeEmpty();
                
                n_epochs = size(obs_set.obs, 1);
                this.n_epochs = n_epochs;
                n_stream = size(diff_obs, 2); % number of satellites
                n_clocks = n_epochs; % number of clock errors
                
                % Store amb_idx
                n_amb = max(max(amb_idx));
                amb_idx = Core_Utils.remEmptyAmbIdx(amb_idx,n_amb);
                amb_flag = 1;
                this.amb_idx = amb_idx;
                this.go_id_amb = obs_set.go_id;
                
                % get ambiguity wl
                wl_amb = zeros(size(amb_obs_count));
                for s = 1 : size(amb_idx, 2)
                    wl_amb(unique(amb_idx(~isnan(amb_idx(:, s)),s))) = obs_set.wl(s);
                end
            else
                n_amb = 0;
                amb_flag = 0;
                this.amb_idx = [];
            end
            
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
            % if multiple observations types are present inter observations biases need be compouted
            iob_idx = zeros(size(obs_set.wl));
            for c = 1 : size(u_obs_code, 1)
                idx_b = idxCharLines(obs_set.obs_code, u_obs_code(c, :));
                iob_idx(idx_b) = c - 1;
            end
            iob_p_idx = iob_idx + n_coo; % progressive index start for iob
            n_iob = size(u_obs_code, 1) - 1;
            iob_flag = double(n_iob > 0);
            
            % total number of observations
            n_obs = sum(sum(diff_obs ~= 0));
            
            % Building Design matrix
            n_par = n_coo_par + iob_flag + amb_flag + double(tropo) + 2 * double(tropo_g); % three coordinates, 1 clock, 1 inter obs bias(can be zero), 1 amb, 3 tropo paramters
            A = zeros(n_obs, n_par); % three coordinates, 1 clock, 1 inter obs bias(can be zero), 1 amb, 3 tropo paramters
            obs = zeros(n_obs, 1);
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
                if ~rec.isFixed()
                    A(lines_stream, 1:3) = - los_stream;
                end
                if dynamic & ~rec.isFixed()
                    A_idx(lines_stream, 1) = ep_p_idx(id_ok_stream);
                    A_idx(lines_stream, 2) = n_epochs + ep_p_idx(id_ok_stream);
                    A_idx(lines_stream, 3) = 2*n_epochs + ep_p_idx(id_ok_stream);
                elseif ~isempty(pos_idx_vec)
                    A_idx(lines_stream, 1) = pos_idx_vec(id_ok_stream);
                    A_idx(lines_stream, 2) = n_pos  + pos_idx_vec(id_ok_stream);
                    A_idx(lines_stream, 3) = 2*n_pos + pos_idx_vec(id_ok_stream);
                end
                % ----------- Inster observation bias ------------------
                if n_iob > 0
                    A(lines_stream, n_coo_par + 1) = iob_idx(s) > 0;
                    A_idx(lines_stream, n_coo_par + 1) = max(n_coo+1, iob_p_idx(s));
                end
                % ----------- Abiguity ------------------
                if phase_present
                    amb_offset = n_coo_par + iob_flag + 1;
                    if not(flag_amb_fix)
                        A(lines_stream, amb_offset) = 1;%obs_set.wl(s);
                    else
                        A(lines_stream, amb_offset) = 1;
                    end
                    A_idx(lines_stream, amb_offset) = n_coo + n_iob + amb_idx(id_ok_stream, s);
                end
                % ----------- Clock ------------------
                A(lines_stream, n_coo_par+iob_flag+amb_flag + 1) = 1;
                A_idx(lines_stream, n_coo_par+iob_flag+amb_flag + 1) = n_coo + n_iob + n_amb + ep_p_idx(id_ok_stream);
                % ----------- ZTD ------------------
                if tropo
                    A(lines_stream, n_coo_par+iob_flag+amb_flag + 2) = mfw_stream;
                    A_idx(lines_stream, n_coo_par+iob_flag+amb_flag + 2) = n_coo + n_clocks + n_iob + n_amb + ep_p_idx(id_ok_stream);
                end
                % ----------- ZTD gradients ------------------
                if tropo_g
                    cotan_term = cot(el_stream) .* mfw_stream;
                    A(lines_stream, n_coo_par+iob_flag+amb_flag + 3) = cos(az_stream) .* cotan_term; % noth gradient
                    A(lines_stream, n_coo_par+iob_flag+amb_flag + 4) = sin(az_stream) .* cotan_term; % east gradient
                    
                    A_idx(lines_stream, n_coo_par+iob_flag+amb_flag + 3) = n_coo + 2 * n_clocks + n_iob + n_amb + ep_p_idx(id_ok_stream);
                    A_idx(lines_stream, n_coo_par+iob_flag+amb_flag + 4) = n_coo + 3 * n_clocks + n_iob + n_amb + ep_p_idx(id_ok_stream);
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
                G = [zeros(1, n_coo + n_iob) ones(1,n_amb) -ones(1,n_clocks)]; % <- This is the right one !!!
                if tropo
                    G = [G zeros(1, n_clocks)];
                end
                if tropo_g
                    G = [G zeros(1, 2*n_clocks)];
                end
                D = 0;
                this.G = G;
                this.D = D;
            end
            this.A_ep = A;
            this.A_idx = A_idx;
            this.variance = variance;
            this.y = y;
            this.epoch = obs;
            this.sat = sat;
            if dynamic
                this.param_flag = [1, 1, 1, -ones(iob_flag), -ones(amb_flag), 1, ones(tropo), ones(tropo_g), ones(tropo_g)];
            else
                this.param_flag = [zeros(1,n_coo_par) -ones(iob_flag), -ones(amb_flag), 1, ones(tropo), ones(tropo_g), ones(tropo_g)];
            end
            this.param_class = [1:n_coo_par, 4 * ones(iob_flag), 5*ones(amb_flag), 6, 7*ones(tropo), 8*ones(tropo_g), 9*ones(tropo_g)];
            if phase_present
                system_jmp = find([sum(nan2zero(diff(amb_idx)),2)] == sum(~isnan(amb_idx(1 : end - 1, :)),2) | [sum(nan2zero(diff(amb_idx)),2)] == sum(~isnan(amb_idx(2 : end, :)),2));
                fprintf('#### DEBUG #### \n');
                [[1; system_jmp + 1] [system_jmp; max(obs)]]
                this.system_split = [[1; system_jmp + 1] [system_jmp; max(obs)]];
            else
                this.system_split = [1 max(obs)];
            end
        end
        
        function setUpSDNetworkAdj(this, rec_list, datum_definition, rate, obs_type, sys_c)
            % set up a network adjustement for the receiver in the list using sigle difference ( rec - rec)  observations
            % SYNTAX: ls_manipulator = Network.setUpNetworkAdj(rec_list)
            % INPUT:
            %      rec_list = list of receivers (cell)
            %      datum_definition = struct with field .type [ F (Free network) | C Costrained ( either by lagrange multiplier or pseudo observations)]
            %                                           .station matrix,
            %                                           case free netwrok [st_nums]
            %                                           case constarined network [ st_num (1 hard 2 soft) sigma_x sigma_y sigma_z] if hard sigma ccan be left to 0 since the will be ignored
            %      rate : rate of observations in seconds
            %      obs_type: I : iono free, 1 first frequency, 2 second frequency, 3 third frequency
            % NOTE: if one wants to set up a single baseline processing he can simpy set up and hard contrsaint on the coordinates of the master see (setUpBaselineAdj)
            % IMPRTANT: when specifing the rate it is assumed for instance that 30 seconds observations are sampled in the vivicinity if 00 and 30 seconds of the minute, if
            %           this is not the case the method will not work
            if nargin < 5
                sys_c = this.state.cc.getAvailableSys();
            end
            % check consistency onf obs_type:
            if ~isempty(strfind(obs_type,'I')) && length(obs_type > 1)
                obs_type = 'I';
                this.log.addWarning('Network Adjustement: When using Iono Free combination no other frequency or combination is allowed');
            end
            
            obs_sets  = {};
            % rates of all recievers this serve also as an indicatpr of if the receiver is used in fact it remains nan of the receiver is not used
            rates = nan(length(rec_list),1);
            
            min_time = GPS_Time(datenum([9999 0 0 0 0 0])); % Warning! to be updated in year 9999
            % get Observations
            for i = 1 : length(rec_list)
                % chek if sampling is compatible
                rate_ratio = max(rate, rec_list(i).getRate()) / min(rate, rec_list(i).getRate());
                st_time_ratio = rec_list(i).time.first.getSecond() / rec_list(i).getRate();
                n_skip_epochs = [];
                if abs( rate_ratio - round(rate_ratio)) > 10e-9
                    this.log.addWarning(sprintf('Receiver %s not compatible with selected rate, skipping',rec_list(i).marker_name));
                    obs_sets{end+1} = [];
                elseif abs( st_time_ratio - round(st_time_ratio)) > 10e-3
                    this.log.addWarning(sprintf('Receiver %s has starting time not multiple of sampling rate. Not supported skipping',rec_list(i).marker_name));
                    obs_sets{end+1} = [];
                else
                    if abs(st_time_ratio - round(st_time_ratio)) > 1e-3 % chek if we have to skip epochs
                        skip_time = ceil(st_time_ratio)*rate - second;
                        n_skip_epochs = skip_time / rec_list(i).getRate();
                    end
                    if rec_list(i).time.first < min_time % detemine the staring time to form the single diffferences
                        second = rec_list(i).time.first.getSecond();
                        st_time_ratio = second / rate;
                        if abs(st_time_ratio - round(st_time_ratio)) > 1e-3 % chek if we have to skip epochs
                            min_time = rec_list(i).time.first ;
                            mid_time.addSeconds(skip_time); %skip epochs and round the time
                            
                        else % round to an epoch muliple of rate
                            round_time =round(st_time_ratio)*rate - second;
                            min_time = rec_list(i).time.first ;
                            mid_time.addSeconds(round_time);
                        end
                    end
                    obs_set = Observation_Set();
                    for j = 1 : length(obs_type) % get the observiont set
                        if obs_type(j) == 'I'
                            for sys_c = rec.cc.sys_c
                                for i = 1 : length(obs_type)
                                    obs_set.merge(rec_list(i).getPrefIonoFree(obs_type(i), sys_c));
                                end
                            end
                        else
                            for sys_c = rec.cc.sys_c
                                f = rec.getFreqs(sys_c);
                                for i = 1 : length(obs_type)
                                    obs_set.merge(rec_list(i).getPrefObsSetCh([obs_type(i) obs_type(j)], sys_c));
                                end
                            end
                        end
                    end
                    if ~isempty(n_skip_epochs) % remove epochs to homogenize start times
                        obs_set.removeEpochs([1 : n_skip_epochs])
                    end
                    obs_sets{end+1} = obs_set; % add to the list
                    rates(i) = rec_list(i).getRate();
                    
                end
            end
            % get offset from min epochs and estiamted num of valid epochs
            
            start_off_set = nan(length(rec_list),1);
            est_num_ep = nan(length(rec_list),1);
            for i = 1 : length(rec_list)
                if ~ismepty(obs_sets{i})
                    start_off_set(i) = round((obs_sets{i}.time.first - min_time)/rate); % off ste from min epochs
                    est_num_ep = round(obs_sets{i}.time.length / rate);
                end
            end
            %estimate roughly an upper bound for the numebr of observations, so that memeory can be preallocated
            % num_pochs x num_sat/2 x num_valid_receiver
            num_sat = this.state.cc.getNumSat(sys_c);
            est_n_obs = max(est_num_ep,'omitnan') * num_sat/2 * sum(~isnan(rates));
            n_par = n_coo + iob_flag + amb_flag + double(tropo) + 2 * double(tropo_g); % three coordinates, 1 clock, 1 inter obs bias(can be zero), 1 amb, 3 tropo paramters
            A = zeros(est_n_obs, n_par);
            A_idx = zeros(est_n_obs, n_par);
            
            cond_stop = true;
            while cond_stop
                for i = 1 : num_sat
                    
                end
                
            end
            
            
        end
        
        function setUpBaselineAdj(this, receivers, master)
            % set up baseline adjustment
            % SYNTAX: ls_manipulator = Network.setUpBaselineAdj(receivers, master)
            % INPUT: receivers = {rec1 , rec2}
            %        master = 1  first is master , 2 second is master
            % NOTE: wrapper to setUPNetworkAdj function
            datum_definition.type = 'C';
            datum_definition.station = [master 1 0 0 0]; % hard constarint on master
            ls_manipulator = setUpSDNetworkAdj(receivers, datum_definition);
            
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
            else %if not prestn add it
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
            res_l = zeros(size(this.y));
            for o = 1:size(this.A_ep, 1)
                res_l(o) = this.y(o) - this.A_ep(o, :) * x(this.A_idx(o, :), 1);
            end
            this.res = res_l;
            res_l(this.rw == 0) = 0;
            n_epochs = max(this.true_epoch);
            n_sat = max(this.sat_go_id);
            res = zeros(n_epochs, n_sat);
            for i = 1:length(this.sat_go_id)
                idx = this.sat == i;
                ep = this.epoch(idx);
                res(this.true_epoch(ep), this.sat_go_id(i)) = res_l(idx);
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
            wfun = @(x) - exp(x.^2 ./threshold.^2);
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
            idx_constant_l = this.param_flag == 0 | this.param_flag == -1;
            idx_constant = find(idx_constant_l);
            idx_non_constant = find(~idx_constant_l);
            a_idx_const =unique(this.A_idx(:, idx_constant_l));
            a_idx_ep_wise = unique(this.A_idx(:, ~idx_constant_l));
            n_constant = length(a_idx_const);
            n_class = size(this.A_ep, 2);
            n_ep_wise = max(max(this.A_idx(:, ~idx_constant_l))) - n_constant;
            if isempty(n_ep_wise)
                n_ep_wise = 0;
            end
            n_epochs = this.n_epochs;
            n_obs = size(this.A_ep, 1);
            n_ep_class = n_ep_wise / n_epochs;
            Ncc = zeros(n_constant, n_constant);
            Nce = zeros(n_ep_wise, n_constant);
            n_class_ep_wise = length(idx_non_constant);
            Ndiags = zeros(n_class_ep_wise, n_class_ep_wise, n_epochs); %permute(this.N_ep(~idx_constant_l,~idx_constant_l,:),[3,1,2]);
            B = zeros(n_constant+n_ep_wise, 1);
            if isempty(this.rw)
                this.rw = ones(size(this.variance));
            end
            %%% all costant parameters are put before in the normal matrix find the mapping between A_idx and idx in the Normal matrix
            N2A_idx = [a_idx_const; a_idx_ep_wise];
            A2N_idx = zeros(size(N2A_idx));
            A2N_idx(N2A_idx) = 1:(n_constant+n_ep_wise);
            
            for i = 1 : n_obs
                p_idx = this.A_idx(i, :);
                p_idx(p_idx == 0) = 1;  % does not matter since terms are zeros
                N_ep = this.N_ep(:, :, i);
                A_ep = this.A_ep(i, :);
                variance = this.variance(i);
                rw = this.rw(i);
                y = this.y(i);
                e = this.epoch(i);
                p_c_idx = p_idx(idx_constant_l);
                p_e_idx = p_idx(~idx_constant_l);
                p_e_idx(p_e_idx <= 0) = 1;  % does not matter since terms are zeros
                p_c_idx = A2N_idx(p_c_idx);
                p_e_idx = A2N_idx(p_e_idx)-n_constant;
                p_idx =A2N_idx(p_idx);
                
                % fill Ncc
                Ncc(p_c_idx, p_c_idx) = Ncc(p_c_idx, p_c_idx) + N_ep(idx_constant, idx_constant);
                % fill Nce
                Nce(p_e_idx, p_c_idx) = Nce(p_e_idx, p_c_idx) + N_ep(idx_non_constant, idx_constant);
                %fill Ndiags
                
                Ndiags(:, :, e) = Ndiags(:, :, e) + N_ep(idx_non_constant, idx_non_constant);
                %fill B
                B(p_idx) = B(p_idx) + A_ep' * (1 ./ variance) * rw * y;
            end
            Nee = [];
            class_ep_wise = this.param_class(idx_non_constant);
            
            rate = median(diff(this.true_epoch));
            diff_reg = 1./double(diff(this.true_epoch));
            reg_diag0 = [diff_reg; 0 ] + [0; diff_reg];
            reg_diag1 = -diff_reg ;
            Ndiags = permute(Ndiags, [3, 1, 2]);
            tik_reg = ones(n_epochs,1)/n_epochs; %%% TIkhonov on ZTD and gradients
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
            N = [[Ncc, Nce']; [Nce, Nee]];
            if ~ isempty(this.G)
                G = this.G(N2A_idx);
                N =  [[N, G']; [G, zeros(size(G,1))]];
                B = [B; this.D];
            end
            if nargout > 3
                %inverse by partitioning, taken from:
                % Mikhail, Edward M., and Friedrich E. Ackermann. "Observations and least squares." (1976). pp 447
                %{
                Ncc = sparse(Ncc);
                Nce = sparse(Nce);
                invNcc = (Ncc)^(-1);
                invNee = (Nee)^(-1);
                a22ia21 = invNee * Nce;
                invN11 = (Ncc - (Nce') * a22ia21)^(-1);
                invN12 = -invN11 * (a22ia21');
                invN21 = invN12';
                invN22 = invNee - a22ia21 * invN12;
                Cxx = [[invN11; invN21], [invN12; invN22]];
                %}
                Cxx = inv(N);
                x = Cxx * B;
                
            else
                x = N \ B;
                %[x, flag] =  pcg(N,B,1e-9, 10000);
            end
            x_class = zeros(size(x));
            for c = 1:length(this.param_class)
                idx_p = A2N_idx(this.A_idx(:, c));
                x_class(idx_p) = this.param_class(c);
            end
            if (this.state.flag_amb_fix && length(x(x_class == 5,1))> 0)
                amb = x(x_class == 5,1);              
                amb_nl = amb / (0.1070/2);
                %amb_nl = amb * (f_vec(1) + f_vec(2))/f_vec(1);
                frac_part = nan(size(amb));
                idx_fix = false(size(amb));
                amb_nl_fix = frac_part;
                n_amb = max(max(this.amb_idx));
                n_ep = size(this.wl_amb,1);
                n_coo = max(this.A_idx(:,3));
                for i = 1 : n_amb
                    sat = this.sat_go_id(this.sat(this.A_idx(:,4)== i+n_coo));
                    idx = n_ep*(sat(1)-1) +  this.true_epoch(this.epoch(this.A_idx(:,4)== i+n_coo));
                    wl_amb =this.wl_amb(idx);
                    wl_amb = wl_amb(1);
                    if mod(wl_amb,2) == 0
                        amb_nl_fix(i) = Core_Utils.round_even(amb_nl(i)); %nearest even
                    else
                        amb_nl_fix(i) = Core_Utils.round_odd(amb_nl(i)); %nearest odd
                    end
                    %amb_nl_fix(i) = round(amb_nl(i));
                    frac_part(i) = (amb_nl(i) - amb_nl_fix(i))/2; 
                    idx_fix(i) = abs(frac_part(i)) < 0.15 & this.wl_fixed(idx(1));
                end
                
                idx_amb = find(x_class == 5);
               
                idxFix2idxFlo = 1 : length(x);
                idxFlo2idxFix = nan(length(x),1);
                A_fixed = false(size(this.A_idx(:,4)));
                for i = 1 : length(idx_fix)
                    if idx_fix(i)
                        Ni = N(:,idx_amb(i));
                        B = B - Ni * amb_nl_fix(i) * (0.1070 / 2);
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
                % ---------------- consider second round
                x_old = x;
                x(idxFix2idxFlo) = xf;
                x(idx_amb(idx_fix)) = amb_nl_fix(idx_fix) * (0.1070 / 2);
                this.log.addMessage(this.log.indent(sprintf('%d of %d ambiguity fixed\n',sum(idx_fix),length(idx_fix))));
                this.log.addMessage(this.log.indent(sprintf('%.2f %% of observation has the ambiguity fixed\n',sum(A_fixed)/length(A_fixed)*100)));
                
            end
            if nargout > 1
                x_res = zeros(size(x));
                x_res(N2A_idx) = x(1:end-size(this.G,1));
                res = this.getResiduals(x_res);                
                s0 = mean(abs(res(res~=0)));
                if nargout > 3
                    Cxx = s0 * Cxx;
                end
            end
            x = [x, x_class];
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
    end
end
