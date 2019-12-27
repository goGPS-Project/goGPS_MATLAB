%   CLASS LS_Manipulator
% =========================================================================
%
% DESCRIPTION
%   Manipulate least squares system
%
% EXAMPLE
%   LSM = LS_Manipulator();
%
% SEE ALSO
%   - Least_Square
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
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
classdef LS_Manipulator < Exportable
    
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
        PAR_TROPO_V = 10;
        PAR_PCO_X = 10;
        PAR_PCO_Y = 11;
        PAR_PCO_Z = 12;
        PAR_SAT_CLK = 13;
        PAR_ANT_MP = 14;
        CLASS_NAME = {'X', 'Y', 'Z', 'ISB', 'AMB', 'REC_CLK', 'TROPO', 'TROPO_N', 'TROPO_E', 'PCO_X', 'PCO_Y', 'PCO_Z', 'SAT_CLK', 'ANT_MP'};
    end
    
    properties
        A_ep         % Stacked epoch-wise design matrices [n_obs x n_param_per_epoch]
        A_idx        % index of the paramter [n_obs x n_param_per_epoch]
        A_idx_mix    % inded of the paramter for multi receiver session [n_obs x n_param_per_epoch]
        amb_idx      % index of the columns per satellite
        tropo_idx    % index of tropo_par;
        tropo_g_idx  % index of tropo gradients par
        tropo_time_start % index of tiem start
        go_id_amb    % go ids of the amb idx
        phase_idx       % index of the phase measurements
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
        param_flag          % 0: constant in time always same param, -1: constant in time differents param (e.g ambiguity), 1: same param changing epochwise
        time_regularization % [ param_class time_varability] simple time regularization constructed from psudo obs p_ep+1 - p_ep = 0 with given accuracy
        mean_regularization
        true_epoch          % true epoch of the epoch-wise paramters
        time                % true eoch in gpstime
        rec_time_idxes      % for each receiver tell which epochof the common time are used
        rate                % rate of the true epoch
        sat_go_id           % go id of the sat indexes
        receiver_id         % id of the receiver, in case of differenced observations two columns are used
        
        
        wl_amb              % wide-lane ambuiguity
        wl_fixed            % is wide-lane fixed
        
        amb_set_jmp         % cell containing for each receiver the jmps on all ambiguity
        
        network_solution = false;
        
        sat_jmp_idx         % satellite jmp index
        
        pos_indexs_tc = {}  % to which index of the sampled time the progessive index correspond
        central_coo         % index of the central coordinate
        
        apriori_info % previous knowledge about the state to be estimated (for now only ambiguity)
        x_float
        Cxx_amb
        
        is_tropo_decorrel % are tropo paramter decorrelated enough
        is_coo_decorrel % are coo paramter decorrelated enough
        
        ant_mp_est = false;  % estimate antenna multipath
        
        dist_matr = [];
        distance_regularization = [];
        
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
        function this = LS_Manipulator(cc)
            % Creator Brahma
            this.init();
            if nargin < 1
                cc = Core.getState.getConstellationCollector;
            end
            this.cc = cc;
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
        
        function id_sync = setUpPPP(this, rec, sys_list, id_sync,  cut_off, dynamic, pos_idx, tropo_idx_vec)
            % Init the object for the phase stand alone positioning
            %
            % SYNTAX
            %   id_sync = this.setUpPPP(rec, sys_list, id_sync,  <cut_off>, <dynamic>, <pos_idx>)                  

            if nargin < 5
                cut_off = [];
            end
            if nargin < 6
                dynamic = false;
            end
            if nargin < 7
                pos_idx = [];
            end
            if nargin < 8
                tropo_idx_vec = [];
            end
            id_sync = this.setUpSA(rec, sys_list, id_sync, 'L', cut_off, '', dynamic, pos_idx,tropo_idx_vec);
        end
        
        function id_sync = setUpCodeSatic(this, rec, sys_list, id_sync, cut_off)
            % Init the object for the code static positioning
            %
            % SYNTAX
            %   id_sync = setUpCodeSatic(this, rec, sys_c, id_sync, <cut_off>)
            if nargin < 5
                cut_off = [];
            end
            id_sync = this.setUpSA(rec, sys_list, id_sync, 'C', cut_off);
        end
        
        function id_sync = setUpCodeDynamic(this, rec, sys_list, id_sync, cut_off)
            if nargin < 5
                cut_off = [];
            end
            id_sync = this.setUpSA(rec, sys_list, id_sync, 'C', cut_off, '', true);
        end
        
        function id_sync_out = setUpSA(this, rec, sys_list, id_sync_in, obs_type, cut_off, custom_obs_set, dynamic, pos_idx_vec, tropo_rate)
            % Init the object for static positioning
            % return the id_sync of the epochs to be computed
            %
            % INPUT:
            %    rec : receiver
            %    sys_list : list of constellations to be used
            %    id_sync : epoch to be used
            %    obs_type : 'C' 'L' 'CL'
            %    cut_off : cut off angle [optional]
            %
            % SYNTAX
            %   id_sync_out = this.setUpSA(rec, sys_list, id_sync_in, obs_type('C'/'L'/'CL'), cut_off, custom_obs_set, <dynamic>, <pos_idx_vec>)

            if isempty(sys_list)
                sys_list = rec.getActiveSys;
            end
            if nargin < 10
                tropo_rate = [];
            end
            if nargin < 9
                pos_idx_vec = [];
            end
            if nargin < 8
                dynamic = false;
            end   
                        
            % Extract the observations to be used for the solution
            phase_present = instr(obs_type, 'L');
            if phase_present && isempty(rec.findObservableByFlag('L'))
                this.log.addError(sprintf('Receiver %s does not contains phase observations for computing the requested solution', rec.parent.getMarkerName()));
                id_sync_out = [];
            else
                flag_amb_fix = this.state.getAmbFixPPP();
                if nargin < 7 || isempty(custom_obs_set)
                    obs_set = Observation_Set();
                    if rec.isMultiFreq() && ~rec.state.isIonoExtModel %% case multi frequency
                        
                        % Using smoothed iono from geometry free                        
                        for sys_c = sys_list
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
                            obs_set.wl = ones(size(obs_set.wl))*Core_Utils.V_LIGHT / (f_vec(1) + f_vec(2)); % <- set wavelength as narrow lane   
                        end
                    else
                        % Using the best combination available
                        for sys_list = rec.getActiveSys
                            f = rec.getFreqs(sys_list);
                            for i = 1 : length(obs_type)
                                if ~isempty(f)
                                    c_o_s = rec.getPrefObsSetCh([obs_type(i) num2str(f(1))], sys_list);
                                    c_o_s.sigma = c_o_s.sigma + rec.getResidualIonoError;
                                    obs_set.merge(c_o_s);
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
                if not(phase_present)
                    obs_set.wl(:) = -1;
                end
                
                if phase_present
                    cc = Core.getState.getConstellationCollector;
                    n_sat = cc.getMaxNumSat();
                    rec.sat.cycle_slip = false(rec.time.length, n_sat);
                    rec.sat.cycle_slip(:,obs_set.go_id) = obs_set.cycle_slip;
                    rec.sat.cycle_slip = sparse(rec.sat.cycle_slip);
                    rec.sat.outliers = sparse(false(rec.time.length, n_sat));
                    dual_freq = size(obs_set.obs_code,2) > 5;
                    [~, ~, ph_idx] = rec.getPhases();
                    obs_code_ph = rec.obs_code(ph_idx,:);
                    go_id_ph = rec.go_id(ph_idx);
                    ph_idx = find(obs_set.obs_code(:,2) == 'L');
                    o_ph_goi = obs_set.go_id(ph_idx);
                    for s = 1 : length(o_ph_goi)
                        g = o_ph_goi(s);
                        obs_code1 = obs_set.obs_code(ph_idx(s),2:4);
                        if dual_freq
                            obs_code2 = obs_set.obs_code(ph_idx(s),5:7);
                        else
                            obs_code2 = '   ';
                        end
                        out_idx = strLineMatch(obs_code_ph,obs_code1) & go_id_ph == g;
                        out = rec.sat.outliers_ph_by_ph(:,out_idx);
                        if strcmp(obs_code2,'   ')
                            out_idx = strLineMatch(obs_code_ph,obs_code2) & go_id_ph == g;
                            if any(out_idx)
                                out(:,out_idx) = out(:,out_idx) | rec.sat.outliers_ph_by_ph(:,out_idx);
                            end
                        end
                        rec.sat.outliers(:,g) = out;
                    end
                end
                % check presence of snr data and fill the gaps if needed
                if ~isempty(obs_set.snr)
                    snr_to_fill = (double(obs_set.snr ~= 0) + 2 * double(obs_set.obs ~= 0)) == 2; % obs if present but snr is not
                    if sum(sum(snr_to_fill))
                        obs_set.snr = simpleFill1D(obs_set.snr, snr_to_fill);
                    end
                end
                
                % remove epochs based on desired sampling
                if nargin > 3
                    obs_set.keepEpochs(id_sync_in);
                end
                
                % re-apply cut off if requested
                if nargin > 5 && ~isempty(cut_off) && sum(sum(obs_set.el)) ~= 0
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
                    
                    % remove short arcs
                    
                    % ambiguity number for each satellite
                    amb_obs_count = histcounts(serialize(amb_idx), 'Normalization', 'count', 'BinMethod', 'integers');
                    if any(amb_idx(:) == 0)
                        amb_obs_count(1) = [];
                    end
                    assert(numel(amb_obs_count) == max(amb_idx(:))); % This should always be true
                    id = 1 : numel(amb_obs_count);
                    ko_amb_list = id(amb_obs_count < min_arc);
                    amb_obs_count(amb_obs_count < min_arc) = [];
                    for ko_amb = fliplr(ko_amb_list)
                        id_ko = amb_idx == ko_amb;
                        obs_set.remObs(id_ko,false)
                    end
                    
                end
                
                % Remove all the observables for an epoch with less than this.state.getMinNSat()
                epoch_ko = sum(obs_set.obs ~= 0, 2) < this.state.getMinNSat() & sum(obs_set.obs ~= 0, 2) > 0;
                obs_set.obs(epoch_ko, :) = 0;
                
                idx_valid_ep_l = sum(obs_set.obs ~= 0, 2) > 0;
                obs_set.setZeroEpochs(~idx_valid_ep_l);
                
                obs_set.remEmptyColumns();
                % if phase observations are present check if the computation of troposphere parameters is required
                
                if phase_present
                    tropo = this.state.flag_tropo;
                    tropo_g = this.state.flag_tropo_gradient;
                else
                    tropo = false;
                    tropo_g = false;
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
                idx_empty_ep = sum(diff_obs ~= 0,2) <= 1;
                obs_set.remEpochs(idx_empty_ep);
                obs_set.sanitizeEmpty();
                if dynamic
                    pos_idx_vec = 1:obs_set.time.length;
                end
                [A, A_idx, ep, sat, p_flag, p_class, y, variance, amb_set_jmp, id_sync_out] = this.getObsEq( rec, obs_set, pos_idx_vec, tropo_rate);
                this.true_epoch = id_sync_out;
                this.A_ep = A;
                this.A_idx = A_idx;
                this.variance = variance;
                this.y = y;
                this.receiver_id = ones(size(y));
                this.epoch = ep;
                %             sat2progsat = zeors(max(obs_set.go_id));
                %             sat2progsat(obs_set.go_id) = 1: length(unique(obs_set.go_id));
                %             this.sat = sat2progsat(sat);
                this.sat = sat;
                if dynamic
                    p_flag(p_class == this.PAR_X | p_class == this.PAR_Y  | p_class == this.PAR_Z ) = 1;
                end
                this.param_flag = p_flag;
                this.param_class = p_class;
                this.amb_idx = obs_set.getAmbIdx();
                this.go_id_amb = obs_set.go_id;
                this.sat_go_id = obs_set.go_id;
                n_epochs = size(obs_set.obs, 1);
                this.n_epochs = n_epochs;
                phase_s = find(obs_set.wl ~=  -1);
                this.phase_idx = phase_s;
                
                
                %---- Set up the constraint to solve the rank deficeny problem --------------
                if phase_present
                    % Ambiguity set
                    %G = [zeros(1, n_coo + n_iob) (amb_obs_count) -sum(~isnan(this.amb_idx), 2)'];
                    %                 if ~flag_amb_fix
                    %                     G = [zeros(1, n_coo + n_iob + n_apc)  ones(1,n_amb)  -ones(1,n_clocks)]; % <- This is the right one !!!
                    %                 else % in case of ambiugty fixing with cnes orbit the partial trace minimization condition gives problems
                    % setting the first clock of each connected set of arc to 0
                    %system_jmp = find([sum(nan2zero(diff(amb_idx)),2)] == sum(~isnan(amb_idx(1 : end - 1, :)),2) | [sum(nan2zero(diff(amb_idx)),2)] == sum(~isnan(amb_idx(2 : end, :)),2));
                    %amb_set_jmp_bnd = [0; amb_set_jmp; this.n_epochs];
                    amb_set_jmp_bnd = [0; amb_set_jmp; this.n_epochs];
                    %                 clock_const = zeros(1,n_clocks);
                    %                 amb_const = zeros(1,n_amb);
                    %                 amb_const(1) = 1;
                    G = [];% [zeros(1, n_coo + n_iob + n_apc)  amb_const  clock_const];
                    amb_par = this.A_idx(:,this.param_class == this.PAR_AMB);
                    min_amb = min(amb_par);
                    max_amb = max(amb_par);
                    max_par = max(this.A_idx(:,end));
                    n_amb = max_amb - min_amb +1;
                    amb_idx = this.amb_idx;
                    clock_const = zeros(1,this.n_epochs);
                    for i = 1: (length(amb_set_jmp_bnd)-1)
                        %clock_const(amb_set_jmp_bnd(i)+1) = 1;
                        amb_const = zeros(1,n_amb);
                        amb_idx_const = noZero(amb_idx((amb_set_jmp_bnd(i)+1):amb_set_jmp_bnd(i+1),:));
                        amb_idx_const = mode(amb_idx_const);
                        if amb_idx_const ~= 0
                            amb_const(amb_idx_const) = 1;
                            G = [G ;[zeros(1, min_amb-1) amb_const clock_const]];
                        end
                    end
                    
                    G = [G zeros(size(G,1),max_par - (max_amb + this.n_epochs) +1 )];
                    
                    D = zeros(size(G,1),1);
                    this.G = G;
                    this.D = D;
                end
                if phase_present
                    % fprintf('#### DEBUG #### \n');
                    lim = [[1; amb_set_jmp + 1] [amb_set_jmp; max(ep)]];
                    Core.getLogger.addMessage(Core.getLogger.indent(sprintf('Dataset is splitted in this way:\n%s', sprintf('%6d -> %6d\n', lim'))))
                    this.system_split = [[1; amb_set_jmp + 1] [amb_set_jmp; max(ep)]];
                else
                    this.system_split = [1 max(ep)];
                end
            end
        end
        
        function [common_time, id_sync]  = setUpNetworkAdj(this, rec_list, coo_rate, wl_struct, frequency_idx)
            % NOTE : free netwrok is set up -> soft constarint on apriori coordinates to be implemnted
            %
            % OUTPUT:
            % common_time : common gps time used for the processing
            n_rec = length(rec_list);
            obs_set_list  = Observation_Set.empty(n_rec,0);
            work_list = [rec_list.work];
            
            %--- for each receiver get one observation set
            for i = 1 : n_rec
                obs_set_list(i) = Observation_Set();
                if ~isempty(wl_struct)%% case multi frequency
                    for sys_c = this.cc.sys_c
                        if this.state.isIonoFree
                            if this.state.getAmbFixNET > 1 % use the same tracking used in the computation of the widelane
                                trcks = wl_struct.combination_codes(wl_struct.combination_codes(:,1) == sys_c,2:end);
                                
                                temp_o_set = work_list(i).getIonoFree(['L',trcks(1,1:2)],['L',trcks(1,3:4)],sys_c);
                            else
                                temp_o_set = work_list(i).getPrefIonoFree('L', sys_c);
                            end
                        else
                            temp_o_set = work_list(i).getSmoothIonoFreeAvg('L', sys_c);
                            temp_o_set.iono_free = true;
                            temp_o_set.sigma = obs_set_list(i).sigma * 1.5;
                        end
                        if this.state.getAmbFixNET > 1
                            f_vec = this.cc.getSys(sys_c).F_VEC;
                            l_vec = this.cc.getSys(sys_c).L_VEC;
                            rnx3_bnd = wl_struct.combination_codes(wl_struct.combination_codes(:,1) == sys_c,[2 4]);
                            id_b1 = find(this.cc.getSys(sys_c).CODE_RIN3_2BAND == rnx3_bnd(1,1));
                            id_b2 = find(this.cc.getSys(sys_c).CODE_RIN3_2BAND == rnx3_bnd(1,2));
                            temp_o_set.obs = nan2zero(zero2nan(temp_o_set.obs) - (nan2zero(wl_struct.amb_mats{i}(:,temp_o_set.go_id)))*f_vec(id_b2)^2*l_vec(2)/(f_vec(id_b1)^2 - f_vec(id_b2)^2));
                            temp_o_set.wl = ones(size(temp_o_set.wl))*Core_Utils.V_LIGHT / (f_vec(id_b1) + f_vec(id_b2)); % <- set wavelength as narrow lane                            
                        end
                        obs_set_list(i).merge(temp_o_set);
                    end
                    
                else
                    for sys_c = intersect(this.cc.getAvailableSys, work_list(i).getAvailableSys)
                        f = work_list(i).getFreqs(sys_c);
                        if ~isempty(f)
                            obs_set_list(i).merge(work_list(i).getPrefObsSetCh(['L' num2str(f(frequency_idx))], sys_c));
                        end
                    end
                    idx_ph = obs_set_list(i).obs_code(:, 2) == 'L';
                    obs_set_list(i).obs(:, idx_ph) = obs_set_list(i).obs(:, idx_ph) .* repmat(obs_set_list(i).wl(idx_ph)', size(obs_set_list(i).obs,1),1);
                end
                obs_set_list(i).sanitizeEmpty();
                
                if this.state.flag_amb_pass && this.state.getCurSession > 1
                    % remove the left buffer it is necessary only to determine if a cycle slip ha occured in the first useful epoch
                    [~,limc] = this.state.getSessionLimits();
                    idx_rm = obs_set_list(i).time < limc.first;
                    obs_set_list(i).remEpochs(idx_rm);
                    % remove the first epoch
                    obs_set_list(i).sanitizeEmpty();
                    if ~isempty(this.apriori_info)
                        if obs_set_list(i).time.length > 0
                            for s = 1:length(obs_set_list(i).go_id)
                                if obs_set_list(i).cycle_slip(1,s) % if there is a cycle slip remove the stored ambigutiy
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
            
            n_obs = 0;
            for r = 1 : n_rec
                n_obs = n_obs + numel(obs_set_list(r).obs);
            end
                                
            if n_obs == 0
                this.log.addError('Network solution failed: no observations are flagged for usage');
                common_time = [];
                id_sync = [];
                return
            end
            
            % Sync obs_sets
            sanitized = false;
            id_sync = nan;
            while ~sanitized && ~isempty(id_sync)
                % remove short arcs and remove "empty" satellites
                [~, id_sync] = obs_set_list.getSyncTimeExpanded();
                if ~isempty(id_sync)
                    id_ko = sum(~isnan(id_sync),2) < 2;
                    for r = 1 : n_rec
                        lid_ko_rec = false(obs_set_list(r).time.length, 1);
                        lid_ko_rec(noNaN(id_sync(id_ko, r))) = true;
                        obs_set_list(r).remEpochs(lid_ko_rec);
                        obs_set_list(r).remShortArc(max(this.state.getMinArc, 1));
                        obs_set_list(r).sanitizeEmpty();
                    end
                    
                    [common_time, id_sync] = obs_set_list.getSyncTimeExpanded();
                    if ~isempty(id_sync)
                        
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
                end
            end
            
            % Exit condition
            if isempty(id_sync)
                this.log.addError('Network solution failed: no observations are flagged for usage');
                common_time = [];
                id_sync = [];
                return
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
            [sss_lim, ~] = this.state.getSessionLimits();
            st_time = sss_lim.first;
            this.tropo_time_start = common_time.first;
            for i = 1 : n_rec
                % get the position idx
                if ~isempty(coo_rate) % first receiver do not need any sub rate since is the reference
                    [pos_idx_nh, pos_idx_tc] = LS_Manipulator.getPosIdx(obs_set_list(i).time, st_time, coo_rate);
                    this.pos_indexs_tc{end+1} = pos_idx_tc; % to be used afterwards to push back postions
                else
                    if this.state.isSepCooAtBoundaries
                        [ss_lim_ext, ss_lim_int] = this.state.getSessionLimits();
                        pos_idx_nh = ones(obs_set_list(i).time.length,1);
                        idx_bf  = obs_set_list(i).time < ss_lim_int.first;
                        idx_aft = obs_set_list(i).time > ss_lim_int.last;
                        
                        if sum(idx_bf) > 0
                            pos_idx_nh = pos_idx_nh +1;
                            pos_idx_nh(idx_bf) = 1;
                            central_coo = 2;
                        else
                            central_coo = 1;
                        end
                        this.pos_indexs_tc{end+1} = unique([ones(sum(idx_bf) > 0) central_coo ones(sum(idx_aft) > 0)+central_coo])'; % to be used afterwards to push back postions
                        this.central_coo(end+1) = central_coo;
                        pos_idx_nh(idx_aft) = max(pos_idx_nh) + 1;
                    else
                        pos_idx_nh = [];
                    end
                end
                order_tropo = this.state.spline_tropo_order;
                order_tropo_g = this.state.spline_tropo_gradient_order;
                tropo_rate = [this.state.spline_rate_tropo*double(order_tropo>0)  this.state.spline_rate_tropo_gradient*double(order_tropo_g>0)];
                [A_rec, Aidx_rec, ep_rec, sat_rec, p_flag_rec, p_class_rec, y_rec, variance_rec, amb_set_jmp] = this.getObsEq(rec_list(i).work, obs_set_list(i), pos_idx_nh, tropo_rate);
                A = [A ; A_rec];
                Aidx = [Aidx; Aidx_rec];
                r2c = find(~isnan(id_sync(:,i)));
                ep = [ep; r2c(ep_rec)];
                sat = [sat; obs_set_list(i).go_id(sat_rec)];
                p_flag = p_flag_rec;
                p_class = p_class_rec;
                y = [y; y_rec];
                variance = [variance; variance_rec];
                r = [r; ones(size(y_rec))*i];
                this.amb_set_jmp{i} = r2c(amb_set_jmp);
                this.amb_idx{i} = obs_set_list(i).getAmbIdx();
                this.go_id_amb{i} = obs_set_list(i).go_id;
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
            this.time = common_time;
            this.rec_time_idxes = id_sync;
            this.rate = time.getRate;
            this.sat = sat;
            this.param_flag = p_flag;
            this.param_class = p_class;
            this.receiver_id = r;
            this.sat_jmp_idx = sat_jmp_idx;
            
            this.network_solution = true;
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
        
        
        function [A, A_idx, ep, sat, p_flag, p_class, y, variance, amb_set_jmp, id_sync_out] = getObsEq(this, rec, obs_set, pos_idx_vec, tropo_rate)
            % get reference observations and satellite positions
            if nargin < 5
                tropo_rate = [];
            end
            tropo = rec.state.flag_tropo && obs_set.hasPhase();
            tropo_g = rec.state.flag_tropo_gradient  && obs_set.hasPhase();
            dynamic = rec.state.rec_dyn_mode;
            phase_present =  obs_set.hasPhase();
            
            [synt_obs, xs_loc] = rec.getSyntTwin(obs_set);
            xs_loc = zero2nan(xs_loc);
            diff_obs = nan2zero(zero2nan(obs_set.obs) - zero2nan(synt_obs));
            % diff_obs = nan2zero(zero2nan(diff_obs) - cumsum(median(Core_Utils.diffAndPred(zero2nan(diff_obs)), 2, 'omitnan'))); % for DEBUGGING: remove receiver clock
            
            if obs_set.hasPhase()
                amb_idx = obs_set.getAmbIdx();
                n_amb = double(max(max(amb_idx)));
                
                amb_flag = 1;
            else
                n_amb = 0;
                amb_flag = 0;
            end
            
            id_sync_out = obs_set.getTimeIdx(rec.time.first, rec.getRate);
            
            % set up requested number of parametrs
            n_epochs = size(obs_set.obs, 1);
            n_stream = size(diff_obs, 2); % number of satellites
            n_clocks = n_epochs; % number of clock errors
            n_tropo = n_clocks; % number of epoch for ZTD estimation
            ep_p_idx = 1 : n_clocks; % indexes of epochs starting from 1 to n_epochs
            
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
            this.log.addMessage(this.log.indent(sprintf('Setting up SA system using %s', reshape([u_obs_code repmat(' ',n_u_obs_code,1)]', 1, n_u_obs_code * (size(u_obs_code, 2) + 1)))));
            % if multiple observations types are present inter observations biases need be computed
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
            
            % Building Design matrix
            order_tropo = this.state.spline_tropo_order;
            order_tropo_g = this.state.spline_tropo_gradient_order;
            tropo_v_g = false && obs_set.hasPhase(); 
            n_par = n_coo_par + iob_flag + 3 * apc_flag + amb_flag + 1 + double(tropo) + double(order_tropo > 0 & tropo)*(order_tropo -1) + 2 * double(tropo_g) + 2*double(order_tropo_g > 0&tropo_g)*(order_tropo_g)+ 16*double(this.ant_mp_est) + double(tropo_v_g); % three coordinates, 1 clock, 1 inter obs bias(can be zero), 1 amb, 3 tropo paramters
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
            
            % Getting mapping faction values
            if tropo || tropo_g
                [~, mfw] = rec.getSlantMF(id_sync_out);
                mfw(mfw  > 60 ) = nan;
                %mfw = mfw(id_sync_out,:); % getting only the desampled values
            end
            phase_s = find(obs_set.wl ~=  -1);
            islinspline = ~isempty(tropo_rate);
            
            if ~isempty(tropo_rate) && sum(abs(tropo_rate)) ~= 0
                if isempty(this.tropo_time_start)
                    
                    delta_tropo_time_sart = 0;
                else
                    delta_tropo_time_sart = round((obs_set.time.first  - this.tropo_time_start)/rec.time.getRate)*rec.time.getRate; % network case make the spline coherent between recievers
                end
            end
            
            if isempty(tropo_rate) || tropo_rate (1) == 0
                n_tropo = n_clocks;
                
            else
                
                %get_tropo_indexes
                edges = 0:tropo_rate(1)/rec.time.getRate:rec.time.length;
                if edges(end) ~= rec.time.length
                    edges = [edges rec.time.length];
                end
                
                tropo_idx = discretize(id_sync_out-1+delta_tropo_time_sart,edges);
                if isempty(this.tropo_time_start);
                    tropo_idx = tropo_idx - tropo_idx(1) +1;
                end
                this.tropo_idx = tropo_idx;
                n_tropo = max(tropo_idx) + order_tropo;
                tropo_dt = rem(id_sync_out-1+delta_tropo_time_sart,tropo_rate(1)/rec.time.getRate)/(tropo_rate(1)/rec.time.getRate);
            end
            if isempty(tropo_rate) || tropo_rate(2) == 0
                n_tropo_g = n_clocks;
            else
                %get_tropo gradients _indexes
                edges = 0:tropo_rate(2)/rec.time.getRate:rec.time.length;
                if edges(end) ~= rec.time.length
                    edges = [edges rec.time.length];
                end
                tropo_g_idx = discretize(id_sync_out-1+delta_tropo_time_sart,edges);
                if isempty(this.tropo_time_start)
                    tropo_g_idx = tropo_g_idx - tropo_g_idx(1) +1;
                end
                this.tropo_g_idx = tropo_g_idx;
                delta_tropo_time_sart = 0;
                n_tropo_g = max(tropo_g_idx) + order_tropo_g;
                tropo_g_dt = rem(id_sync_out-1 + delta_tropo_time_sart,tropo_rate(2)/rec.time.getRate)/(tropo_rate(2)/rec.time.getRate);
            end
            for s = 1 : n_stream
                id_ok_stream = diff_obs(:, s) ~= 0; % check observation existence -> logical array for a "s" stream
                
                obs_stream = diff_obs(id_ok_stream, s);
                % snr_stream = obs_set.snr(id_ok_stream, s); % SNR is not currently used
                el_stream = obs_set.el(id_ok_stream, s) / 180 * pi;
                if tropo || tropo_g
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
                if (this.state.getWeigthingStrategy == 1) || ~any(el_stream)
                    variance(lines_stream) =  obs_set.sigma(s)^2;
                elseif this.state.getWeigthingStrategy == 2
                    variance(lines_stream) =  (obs_set.sigma(s)./3./sin(el_stream)).^2;
                else
                    variance(lines_stream) =  obs_set.sigma(s)^2;
                end
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
                    if sum(phase_s == s) > 0
                        
%                         if this.state.flag_ppp_amb_fix
%                             A(lines_stream, prog_p_col) = 1;
%                         else
                            A(lines_stream, prog_p_col) = obs_set.wl(s);
%                         end
                        A_idx(lines_stream, prog_p_col) = n_coo + n_iob + n_apc + amb_idx(id_ok_stream, phase_s == s);
                    else
                        A_idx(lines_stream, prog_p_col) = 0;
                        A(lines_stream, prog_p_col) = 0;
                    end
                end
                % ----------- Clock ------------------
                prog_p_col = prog_p_col + 1;
                A(lines_stream, prog_p_col) = 1;
                A_idx(lines_stream, prog_p_col) = n_coo + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                % ----------- ZTD ------------------
                if tropo
                    
                    if isempty(tropo_rate)  || tropo_rate(1) == 0
                        prog_p_col = prog_p_col + 1;
                        A(lines_stream, prog_p_col) = mfw_stream;
                        A_idx(lines_stream, prog_p_col) = n_coo + n_clocks + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                    else
                        spline_v = Core_Utils.spline(tropo_dt(id_ok_stream),order_tropo);
                        for o = 1 : (order_tropo + 1)
                            prog_p_col = prog_p_col + 1;
                            A(lines_stream, prog_p_col) = mfw_stream.*spline_v(:,o);
                            A_idx(lines_stream, prog_p_col) = n_coo + n_clocks + n_iob + n_apc + n_amb + tropo_idx(id_ok_stream) + o-1;
                        end
                    end
                end
                % ----------- ZTD gradients ------------------
                if tropo_g
                    %cotan_term = cot(el_stream) .* mfw_stream;
                    cotan_term = 1 ./ ( sin(el_stream).*tan(el_stream) + 0.0032);
                    
                    if isempty(tropo_rate) || tropo_rate(2) == 0
                        prog_p_col = prog_p_col + 1;
                        A(lines_stream, prog_p_col) = cos(az_stream) .* cotan_term; % noth gradient
                        A_idx(lines_stream, prog_p_col) = n_coo + n_clocks + n_tropo + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                        prog_p_col = prog_p_col + 1;
                        A(lines_stream, prog_p_col) = sin(az_stream) .* cotan_term; % east gradient
                        A_idx(lines_stream, prog_p_col) = n_coo + n_clocks + 2*n_tropo + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                    else
                        spline_v = Core_Utils.spline(tropo_g_dt(id_ok_stream),order_tropo_g);
                        for o = 1 : (order_tropo_g + 1)
                            prog_p_col = prog_p_col + 1;
                            A(lines_stream, prog_p_col) = cos(az_stream) .* cotan_term .*spline_v(:,o); % noth gradient spli
                            A_idx(lines_stream, prog_p_col) = n_coo + n_clocks + n_tropo + n_iob + n_apc + n_amb + tropo_g_idx(id_ok_stream) + o-1;
                        end
                        for o = 1 : (order_tropo_g + 1)
                            prog_p_col = prog_p_col + 1;
                            A(lines_stream, prog_p_col) = sin(az_stream) .* cotan_term.*spline_v(:,o); % east gradient
                            A_idx(lines_stream, prog_p_col) = n_coo + n_clocks + n_tropo + n_tropo_g + n_iob + n_apc + n_amb  + tropo_g_idx(id_ok_stream) + o -1;
                        end
                    end
                    
                end
                if tropo_v_g
                       prog_p_col = prog_p_col + 1;
                        A(lines_stream, prog_p_col) = mfw_stream*rec.h_ellips;
                        A_idx(lines_stream, prog_p_col) =n_coo + n_clocks + n_tropo + 2*n_tropo_g + n_iob + n_apc + n_amb + ep_p_idx(id_ok_stream);
                end
                if this.ant_mp_est
                    prog_p_col = prog_p_col + 1;
                    n_el = 7;
                    n_az = 28;
                    [idx, val] = Core_Utils.hemisphereCubicSpline(n_az,n_el,az_stream,el_stream);
                    A(lines_stream, prog_p_col+(0:16)) = sin(az_stream) .* cotan_term.*spline_v(:,o); % east gradient
                    A_idx(lines_stream, prog_p_col+(0:16)) = n_coo + n_clocks + n_tropo + n_tropo_g + n_iob + n_apc + n_amb  + tropo_g_idx(id_ok_stream) + o -1;
                end
                obs_count = obs_count + n_obs_stream;
            end
            % ---- Suppress weighting until solution is more stable/tested
            %w(:) = 1;%0.005;%this.state.std_phase;
            %---------------------
            
            
            ep = obs;
            e_spline_mat_t = ones(1,double(tropo)*(order_tropo+1));
            e_spline_mat_tg = ones(1,double(tropo_g)*(order_tropo_g+1));
            if dynamic
                p_flag = [1, 1, 1, -ones(iob_flag), -repmat(ones(apc_flag),1,3), -ones(amb_flag), 1, ones(tropo), ones(tropo_g), ones(tropo_g), ones(tropo_v_g)];
            else
                p_flag = [zeros(1,n_coo_par) -ones(iob_flag), -repmat(ones(apc_flag),1,3), -ones(amb_flag), 1, (1 -2 * double(order_tropo > 0))*e_spline_mat_t, (1 -2 * double(order_tropo_g > 0))*e_spline_mat_tg, (1 -2 * double(order_tropo_g > 0))*e_spline_mat_tg,];
            end
           
            p_class = [this.PAR_X*ones(~is_fixed) , this.PAR_Y*ones(~is_fixed), this.PAR_Z*ones(~is_fixed), this.PAR_ISB * ones(iob_flag), this.PAR_PCO_X * ones(apc_flag), this.PAR_PCO_Y * ones(apc_flag), this.PAR_PCO_Z * ones(apc_flag),...
                this.PAR_AMB*ones(amb_flag), this.PAR_REC_CLK, this.PAR_TROPO*e_spline_mat_t, this.PAR_TROPO_N*e_spline_mat_tg, this.PAR_TROPO_E*e_spline_mat_tg , this.PAR_TROPO_V * ones(tropo_v_g) ];
            if obs_set.hasPhase()
                amb_set_jmp = find(sum(diff(int32(amb_idx)) < 0, 2) == sum((amb_idx(1 : end - 1, :)) > 0, 2) | sum(diff(int32(amb_idx)) > 0, 2) == sum((amb_idx(2 : end, :)) > 0,2)) + 1;
            else
                amb_set_jmp = [];
            end
            
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
        
        function [res, av_res] = getResiduals(this, x)
            %res_l = zeros(size(this.y));
            %for o = 1 : size(this.A_ep, 1)
            %    res_l(o) = this.y(o) - this.A_ep(o, :) * x(this.A_idx(o, :), 1);
            %end
            %res_l = zeros(size(this.y));
            % speed-up of the previous lines
            av_res = [];
            if any(isnan(x))
                this.log.addError('Some parameters are NaN!');
            end
            aidx_temp = this.A_idx;
            aidx_temp(aidx_temp == 0) = 1;
            res_l = this.y - sum(this.A_ep .* reshape(x(aidx_temp), size(this.A_idx,1), size(this.A_idx,2)),2);
            
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
        function flag_recompute = weightOnResidual(this, wfun, thr, thr_propagate, flag_clean_margin)
            if isempty(this.rw)
                this.rw = ones(size(this.variance));
            end
            s0 = mean(abs(this.res).*this.rw);
            res_n = this.res/s0;
            if nargin > 2
                % propagate outlier flag ( snooping gatt) -----------------------------------------------------------------------------------------------------
                if nargin > 3 && (thr_propagate > 0)
                    sat_err = nan(this.n_epochs, max(this.sat_go_id));
                    sat_err(this.epoch + (this.sat_go_id(this.sat) - 1) * this.n_epochs) = res_n;
                    ssat_err = Receiver_Commons.smoothSatData([],[],sat_err, [], 'spline', 30, 10);
                    
                    idx_ko = Core_Utils.snoopGatt(ssat_err, thr, thr_propagate);
                    
                    % Check derivate for big fluctuation 2 * local std + 5 mm
                    % sensor = Core_Utils.diffAndPred(movmedian(sat_err * s0, 5));
                    % for f = 1 : size(sensor,2)
                    %    % Ideally this should add CS
                    %    tmp_ko = flagExpand(abs(sensor(:,f)) > 2 * movstd(sensor(:,f),13) + 0.005, 1);
                    %    idx_ko(tmp_ko, f) = true;
                    % end

                    idx_rw = idx_ko(this.epoch + (this.sat_go_id(this.sat) - 1) * this.n_epochs);
                    
                    if nargin > 4 && flag_clean_margin
                        % Find begin and end of an arc, if it's out of thr flag it till comes down under threshold --------------------------------------------
                        ot = abs(ssat_err) > thr_propagate;
                        for s = 1 : size(sat_err, 2)
                            
                            % Beginning of the arc
                            tmp = ot(:, s) + ~isnan(sat_err(:, s));
                            
                            id_start_bad = find(diff([0; tmp]) == 2);
                            for i = 1 : numel(id_start_bad)
                                id_stop_bad = find(ot(id_start_bad(i) : end, s) == 0, 1, 'first'); % find when the arc is now under thr
                                
                                % rem over threshold elements
                                if isempty(id_stop_bad)
                                    idx_ko(id_start_bad(i) : end, s) = true;
                                else
                                    idx_ko(id_start_bad(i) + (0 : (id_stop_bad - 1)), s) = true;
                                end
                            end
                            
                            % End of the arc (flip method)
                            ot(:, s) = flipud(ot(:, s));
                            tmp = ot(:, s) + flipud(~isnan(sat_err(:, s)));
                            
                            id_start_bad = find(diff([0; tmp]) == 2);
                            for i = 1 : numel(id_start_bad)
                                id_stop_bad = find(ot(id_start_bad(i) : end, s) == 0, 1, 'first'); % find when the arc is now under thr
                                
                                % rem over threshold elements
                                if isempty(id_stop_bad)
                                    idx_ko(size(idx_ko, 1) + 1 - (id_start_bad(i) : size(idx_ko, 1)), s) = true;
                                else
                                    idx_ko(size(idx_ko, 1) + 1 - (id_start_bad(i) + (0 : (id_stop_bad - 1))), s) = true;
                                end
                            end                                                        
                        end % ---------------------------------------------------------------------------------------------------------------------------------                     
                        idx_rw = idx_rw | idx_ko(this.epoch + (this.sat_go_id(this.sat) - 1) * this.n_epochs);
                    end
                else
                    idx_rw = abs(res_n) > thr;
                end
            else
                idx_rw = true(size(res_n));
            end
            this.rw(idx_rw) =  wfun(res_n(idx_rw));
            flag_recompute = any(idx_rw);
            if sum(this.rw(idx_rw) < 1e-3) > 0 % observation with weight less than 1e-3 are removed from the adjustment otherwise parameters that depend only on them may suffer numerical instability
                flag_recompute = true;
                this.remObs(this.rw < 1e-3);
            end
        end
        
        function flag_recompute = remOverThr(this, thr)
            % Remove all the observations with residuals over a certain threshold
            %
            % SYNTAX:
            %   flag_recompute = this.remOverThr(thr)
            if isempty(this.rw)
                this.rw = ones(size(this.variance));
            end
            res_n = this.res;
            idx_rw = abs(res_n) > thr;
            this.rw(idx_rw) =  0;
            flag_recompute = any(idx_rw);
            if sum(this.rw(idx_rw) < 1e-3) > 0 % observation with weight less than 1e-3 are removed from the adjustment otherwise parameters that depend only on them may suffer numerical instability
                flag_recompute = true;
                this.remObs(this.rw < 1e-3);
            end
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
        
        function flag_recompute = snoopingGatt(this, thr)
            % Outlier detection on first tredshold + remove till under second threshold 
            % that is threshold_propagate = 2 sigma
            %
            % SYNTAX:
            %   flag_recompute = this.snoopingGatt(thr)
            if nargin == 1
                thr = 10;
            end
            threshold_propagate = 2.5;
            wfun = @(x) 0;
            flag_recompute = this.weightOnResidual(wfun, thr, threshold_propagate);
        end
        
        function flag_recompute = snoopingGatt2(this, thr)
            % Outlier detection on first tredshold + remove till under second threshold
            % + remove beginnin and end of arcs over threshold_propagate (2 sigma)
            % SYNTAX:
            %   flag_recompute = this.snoopingGatt2(thr)
            if nargin == 1
                thr = 10;
            end
            threshold_propagate = 2;
            wfun = @(x) 0;
            flag_recompute = this.weightOnResidual(wfun, thr, threshold_propagate, true);
        end
        
        %------------------------------------------------------------------------
        
        function [x, res, s0, Cxx, l_fixed, av_res] = solve(this)
            av_res = [];
            l_fixed = 0;
            Cxx = [];
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
            N2A_idx = [];
            A2N_idx_tot = [];
            x_rec = [];
            for r = u_r(:)'
                idx_r_l = this.receiver_id == u_r(r);
                idx_r = find(idx_r_l);
                A_rec = this.A_idx(idx_r_l, ~idx_rec_common_l);
                idx_constant_l = this.param_flag(~idx_rec_common_l) == 0 | this.param_flag(~idx_rec_common_l) == -1;
                idx_constant = find(idx_constant_l);
                idx_non_constant = find(~idx_constant_l);
                idx_ep_wise = find(this.param_flag == 1);
                u_class = unique(this.param_class);
                n_param_class = zeros(1,length(u_class));
                for u = 1 : length(u_class)
                    n_param_class(u) = max(max(this.A_idx(:,this.param_class == u_class(u)))) - min(min(this.A_idx(:,this.param_class == u_class(u)))) + 1;
                end
                a_idx_const = unique(A_rec(:, idx_constant_l));
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
                % constraint over rate paramters
                if ~isempty(this.time_regularization)
                    cc_class = intersect(this.time_regularization(:, 1),this.param_class(idx_constant));
                    if ~isempty(cc_class)
                        Ncc_reg = sparse(n_constant, n_constant);
                        diag_cc0 = zeros(n_constant,1);
                        diag_cc1 = zeros(n_constant-1,1);
                        idx_const_u = false(size(n_param_class));
                        for b = 1 : length(u_class)
                            idx_const_u(b) = unique(idx_constant_l(this.param_class == u_class(b)));
                        end
                        n_param_class_const_cum = cumsum(n_param_class(idx_const_u));
                        for c = cc_class'
                            idx_c = this.time_regularization(:, 1) == c;
                            w = 1 ./ this.time_regularization(idx_c, 2) ;
                            reg = ones(n_param_class(u_class == c)-1,1)*w;
                            idx_bg = n_param_class_const_cum(find(u_class(idx_const_u) == c)-1)+1;
                            idx_en = n_param_class_const_cum(u_class(idx_const_u) == c)-1;
                            diag_cc1(idx_bg:idx_en) = -reg;
                            reg0 = [0; reg] + [reg; 0];
                            diag_cc0(idx_bg:(idx_en+1)) = reg0;
                        end
                        Ncc_reg = spdiags( diag_cc0,0,Ncc_reg);
                        Ncc_reg = spdiags([0; diag_cc1],1,Ncc_reg);
                        Ncc_reg = spdiags(diag_cc1,-1,Ncc_reg);
                        Ncc = Ncc + Ncc_reg;
                    end
                    
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
                x_rec = [x_rec; r*ones(size(N_rec,1),1)];
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
                
                idx_net_comm = find(this.param_class == this.PAR_TROPO_V);
                
                for j = 2 : n_rec
                    rec_idx = this.receiver_id == j;
                    % update the indexes
                    par_rec_id = [par_rec_id ; j*ones(max(max(this.A_idx(this.receiver_id == j,:))),1)];
                    this.A_idx(rec_idx,this.param_class ~= this.PAR_TROPO_V) = this.A_idx(this.receiver_id == j, this.param_class ~= this.PAR_TROPO_V) + max(max(this.A_idx(this.receiver_id == j-1,this.param_class ~= this.PAR_TROPO_V)));
                    
                    a_idx_const =unique(this.A_idx(rec_idx, idx_constant_l));
                    a_idx_const(a_idx_const == 0) = [];
                    if ~isempty(idx_net_comm)
                        a_idx_ep_wise = unique(this.A_idx(rec_idx, idx_ep_wise(idx_ep_wise ~= idx_net_comm)));
                    else
                        a_idx_ep_wise = unique(this.A_idx(rec_idx, idx_ep_wise));
                    end
                    a_idx_ep_wise(a_idx_ep_wise == 0) = [];
                    
                    N2A_idx = [N2A_idx; a_idx_const; a_idx_ep_wise];
                end
                if ~isempty(idx_net_comm)
                    this.A_idx(:,this.param_class == this.PAR_TROPO_V) = this.A_idx(:,this.param_class == this.PAR_TROPO_V) - min(this.A_idx(:,this.param_class == this.PAR_TROPO_V)) + max(max(this.A_idx(:,this.param_class ~= this.PAR_TROPO_V)));
                    N2A_idx = [N2A_idx; unique(this.A_idx(:,this.param_class == this.PAR_TROPO_V))];
                    if ~isempty(this.distance_regularization) && Core.isGReD
                        N = GReD_Utility.regularizeTropoDist(this,N,u_ep, x_rec,this.dist_matr,this.distance_regularization.fun_tropo);
                        N = GReD_Utility.regularizeGradientsDist(this,N,u_ep,x_rec, this.dist_matr,this.distance_regularization.fun_gradients);
                    end
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
                    N_stack_idx(idx_common, (r - 1) * n_s_r_p + (1 : n_s_r_p))= A2N_idx_tot(this.A_idx(i, :));
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
                % ALL paramters has a rank deficency because the entrance of the image matrixes are very similar and we also estimated the clock of the satellite
                % 2) remove coordinates and tropo paramters of the first receiver
                % we can do that because tropo paramters are slightly constarined in time so evan if they are non present for the first receiver the rank deficecny is avoided

                idx_rec_isb = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_ISB));
                idx_rec_t = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_TROPO));
                idx_rec_tn = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_TROPO_N));
                idx_rec_te = unique(this.A_idx(this.receiver_id == 1,this.param_class == this.PAR_TROPO_E));
                idx_rm = [idx_rec_isb;];%
                if ~this.is_coo_decorrel
                    % strong regularize the mean of the coordinates to zero
                   % for each time span of coordinate find the index of
                    % the paranter
                    idx_rec_x = zeros(max([serializeCell(this.pos_indexs_tc);1])+1,n_rec);
                    idx_rec_y = zeros(max([serializeCell(this.pos_indexs_tc);1])+1,n_rec);
                    idx_rec_z = zeros(max([serializeCell(this.pos_indexs_tc);1])+1,n_rec);
                    for rr = 1:n_rec
                        idx_rec_xtmp =  unique(this.A_idx(this.receiver_id == rr,this.param_class == this.PAR_X));
                        idx_rec_ytmp =  unique(this.A_idx(this.receiver_id == rr,this.param_class == this.PAR_Y));
                        idx_rec_ztmp =  unique(this.A_idx(this.receiver_id == rr,this.param_class == this.PAR_Z));
                        if ~isempty(this.pos_indexs_tc)
                            idx_rec_x(this.pos_indexs_tc{rr},rr) = idx_rec_xtmp;
                            idx_rec_y(this.pos_indexs_tc{rr},rr) = idx_rec_ytmp;
                            idx_rec_z(this.pos_indexs_tc{rr},rr) = idx_rec_ztmp;
                        else
                            idx_rec_x(1,rr) = idx_rec_xtmp;
                            idx_rec_y(1,rr) = idx_rec_ytmp;
                            idx_rec_z(1,rr) = idx_rec_ztmp;
                        end
                    end
                    % set the mean regularization
                    for ss = 1 : size(idx_rec_x,1)
                        if sum(idx_rec_x(ss,:)) > 0
                            idx_rec_xtmp = noZero(idx_rec_x(ss,:));
                            N(idx_rec_xtmp,idx_rec_xtmp) = N(idx_rec_xtmp,idx_rec_xtmp) + 10*ones(length(idx_rec_xtmp));
                            idx_rec_ytmp = noZero(idx_rec_y(ss,:));
                            N(idx_rec_ytmp,idx_rec_ytmp) = N(idx_rec_ytmp,idx_rec_ytmp) + 10*ones(length(idx_rec_ytmp));
                            idx_rec_ztmp = noZero(idx_rec_z(ss,:));
                            N(idx_rec_ztmp,idx_rec_ztmp) = N(idx_rec_ztmp,idx_rec_ztmp) + 10*ones(length(idx_rec_ztmp));
                        end
                    end
                    
                end
                if ~this.is_tropo_decorrel
                    idx_rm = [idx_rm ; idx_rec_t; idx_rec_tn; idx_rec_te];
                end
                % 3 ) remove one clock per epoch for the minim receiver available
                clk_idx = this.param_class == this.PAR_REC_CLK;
                n_epochs = length(unique(this.epoch));
                idx_clk_to_rm = true(n_epochs,1);
                i = 1;
                while sum(idx_clk_to_rm) > 0 && i
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
                    for s = 1 :length(u_sat)
                        jmp_idx = find(diff(this.sat_jmp_idx(:,u_sat(s))) == -1);
                        if ~this.sat_jmp_idx(1,u_sat(s))
                            if ~prev_info || sum(u_sat(s) == this.apriori_info.goids) == 0
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
                            %                            idx_amb_rec = this.A_idx(this.receiver_id == d & this.sat == u_sat(s) & this.epoch >= jmp & this.epoch < jmp2 ,this.param_class == this.PAR_AMB);
                            %                            d = d + 1;
                            %                         end
                            r = 1;
                            not_found = true;
                            while r <= n_rec && not_found % always tak off ftom the first reciever, thi should avoid very nasty case of overcobstarininf
                                idx_amb_rec = this.A_idx(this.receiver_id == r & this.sat == u_sat(s) & this.epoch >= jmp & this.epoch < jmp2 ,this.param_class == this.PAR_AMB);
                                idx_amb_rec = Core_Utils.remBFromA(idx_amb_rec,idx_amb_rm);
                                if ~isempty(idx_amb_rec)
                                    idx_amb_rec = mode(idx_amb_rec);%(1);%idx_amb_rec(min(120,length(idx_amb_rec)));
                                    not_found = false;
                                end
                                idx_amb_rm = [idx_amb_rm; idx_amb_rec];
                                r = r + 1;
                            end
                        end
                    end
                end
                %                 5) remove one ambiguity per each set of disjunt set of arcs of each receiver to resolve the ambiguity-receiver clock rank deficency
                %                 first recievr does not have clocks any more so no rank defricency
                if true
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
                            
                            
                            idx_clock_rec = this.A_idx(this.receiver_id == i & this.epoch >= jmp & this.epoch < jmp2,this.param_class == this.PAR_REC_CLK);
                            %if isempty(intersect(idx_rm,idx_clock_rec))  % chekc i clck has been removed for the receiver in that ambiguity set
                            idx_space = this.receiver_id == i & this.epoch >= jmp & this.epoch < jmp2;
                            idx_amb_rec = this.A_idx(idx_space, this.param_class == this.PAR_AMB);
                            if sum(this.param_class == this.PAR_ISB) > 0
                                 idx_isb = this.A_idx(idx_space, this.param_class == this.PAR_ISB) - double(this.A_ep(idx_space, this.param_class == this.PAR_ISB) == 0); % very not nice trick to have also a differnt index when the value is 0
                            end
                            %if isempty(intersect(idx_amb_rm,idx_amb_rec))
                            %if
                            %                                     g = 1;
                            %                         while sum(idx_amb_rec(g) == idx_amb_rm) > 0 && g < length(idx_amb_rec)
                            %                             g = g +1;
                            %                         end
                            %idx_amb_rec = Core_Utils.remBFromA(idx_amb_rec,idx_amb_rm);
                            if ~isempty(idx_amb_rec)
                                
                                
                                 if sum(this.param_class == this.PAR_ISB) > 0
                                     idx_amb_rect = idx_amb_rec;
                                     idx_amb_rec = [];
                                     u_idx_isb = unique(idx_isb)';
                                     for uii = u_idx_isb
                                         idx_amb_rectt = Core_Utils.remBFromA(idx_amb_rect(idx_isb == uii),idx_amb_rm);
                                         idx_amb_rec = [idx_amb_rec; mode(idx_amb_rectt)];
                                     end
                                 else
                                     idx_amb_rec = mode(idx_amb_rec);%(1);%idx_amb_rec(min(120,length(idx_amb_rec)));
                                 end
                            end
                            idx_amb_rm = [idx_amb_rm; idx_amb_rec];
                            %end
                            %end
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
            
            warning off
            x = N \ B;
            warning on
            
            x_class = zeros(size(x));
            for c = 1:length(this.param_class)
                idx_pars = this.A_idx(:, c);
                idx_p = A2N_idx_tot(idx_pars(idx_pars~=0));
                x_class(idx_p) = this.param_class(c);
            end
            if is_network
                idx_est = true(n_par,1);
                idx_est(idx_rm) = false;
                x_tot(idx_est) = x;
                x = x_tot;
                
                idx_amb_par = find(x_class(idx_est) == this.PAR_AMB);
                n_amb = length(idx_amb_par);
                %                 x_rec = ones(size(x,1),1);
                %                 id_rec = find(diff(x_class) < 0);
                %                 for i = 1 : length(id_rec)
                %                     idx = id_rec(i)+1;
                %                     x_rec(idx:end) = i+1;
                %                 end
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
            if ~is_network
                if (this.state.getAmbFixPPP && ~isempty(x(x_class == this.PAR_AMB,1)))
                    amb = x(x_class == this.PAR_AMB,1);
                    amb_wl_fixed = false(size(amb));
                    amb_n1 = nan(size(amb));
                    amb_wl = nan(size(amb));
                    n_ep_wl = zeros(size(amb));
                    n_amb = max(max(this.amb_idx));
                    n_ep = size(this.wl_amb,1);
                    if sum(this.param_class == this.PAR_X) >0
                    n_coo = max(this.A_idx(:,3));
                    else
                        n_coo = 0;
                    end
                    for i = 1 : n_amb
                        sat = this.sat_go_id(this.sat(this.A_idx(:,this.param_class == this.PAR_AMB) == i+n_coo));
                        idx = n_ep*(sat(1)-1) +  this.true_epoch(this.epoch(this.A_idx(:,this.param_class == this.PAR_AMB)== i+n_coo));
                        if ~isempty( this.wl_fixed) % case local single frequency PPP
                            amb_wl(i) = this.wl_amb(idx(1));
                            amb_wl_fixed(i)=  this.wl_fixed(idx(1));
                        end
                        n_ep_wl(i) = length(idx);
                        amb_n1(i) = amb(i); %(amb(i)- 0*f_vec(2)^2*l_vec(2)/(f_vec(1)^2 - f_vec(2)^2)* wl_amb);  % Blewitt 1989 eq(23)
                        
                    end
                    
                    weight = min(n_ep_wl(amb_wl_fixed),100); % <- downweight too short arc
                    weight = weight / sum(weight);
                    
                    
                    idx_amb = find(x_class == this.PAR_AMB);
                    % get thc cxx of the ambiguities
                    n_amb  = length(idx_amb);
                    b_eye = zeros(length(B),n_amb);
                    idx = sub2ind(size(b_eye),idx_amb,[1:n_amb]');
                    b_eye(idx) = 1;
                    b_eye = sparse(b_eye);
                    Cxx_amb = N\b_eye;
                    Cxx_amb = Cxx_amb(idx_amb,:);
                    idx_constarined = abs(amb_n1) < 1e-5;
                    l_fixed = false(size(amb_n1,1),1);
                    amb_fixed = zeros(size(amb_n1,1),1);
                    % ILS shrinking, method 1
                    [af, is_fixed,lf] = Fixer.fix(amb_n1(~idx_constarined), Cxx_amb(~idx_constarined,~idx_constarined), 'sequential_best_integer_equivariant');
                    
                    if is_fixed
                        amb_fixed(~idx_constarined,:) = af;
                        l_fixed(~idx_constarined,:) = lf;
                        % FIXED!!!!
                        idx_est = true(size(x,1),1);
                        amb_fix = amb_fixed(:, 1);
                        
                        for i = 1 : n_amb
                            Ni = N(:, idx_amb(i));
                            if l_fixed(i)
                                B = B - Ni * amb_fix(i);
                            end
                        end
                        
                        % remove fixed ambiguity from B and N
                        B(idx_amb(l_fixed(:,1))) = [];
                        N(idx_amb(l_fixed(:,1)), :) = [];
                        N(:, idx_amb(l_fixed(:,1))) = [];
                        
                        % recompute the solution
                        idx_nf = true(sum(idx_est,1),1);
                        idx_nf(idx_amb(l_fixed(:,1))) = false;
                        xf = zeros(size(idx_nf));
                        
                        xf(idx_nf) = N \ B;
                        xf(~idx_nf) = amb_fix(l_fixed(:,1));
                        
                        x(idx_est) = xf;
                    else
                        this.log.addWarning('The ambiguities cannot be fixed!!!');
                    end
                    
                end
            else
                if ((this.state.getAmbFixNET > 1) && ~isempty(x(x_class == 5,1)))
                    
                    % Ambiguity fixing
                    
                    % Get ambiguity array
                    
                    xe = x(idx_est);
                    % get some information about the ambiguty (receiver,
                    % satellite, and number of epoch used in the
                    % compesation)
                    amb = xe(idx_amb_par, 1);
                    amb_rec = x_rec(idx_est);
                    amb_rec = amb_rec(idx_amb_par);
                    n_ep_amb = zeros(size(amb));
                    amb_sat = zeros(size(amb));
                    idx_amb_par_tot = find(idx_est); idx_amb_par_tot = idx_amb_par_tot(idx_amb_par);% <- amb indices in the A matrix
                    for i = 1 : n_amb
                        idx = this.A_idx(:,this.param_class == this.PAR_AMB)== idx_amb_par_tot(i);
                        n_ep_amb(i) = sum(idx);
                        idx_first = find(idx,1,'first');
                        amb_sat(i) = this.sat(idx_first);
                    end
                    
                    % Getting tht VCV matrix for the ambiguities
                    b_eye = zeros(length(B), n_amb);
                    idx_t_amb_par = find(x_class(idx_est) == this.PAR_AMB);
                    idx = sub2ind(size(b_eye),idx_t_amb_par, [1:n_amb]');
                    b_eye(idx) = 1;
                    b_eye = sparse(b_eye);
                    C_amb_amb = N \ b_eye;
                    C_amb_amb = C_amb_amb(idx_t_amb_par, :);
                    C_amb_amb = (C_amb_amb + C_amb_amb') ./ 2; % make it symmetric (sometimes it is not due to precion loss)
                    Cxx = C_amb_amb;
                    present_sat = unique(this.sat);
                    a_id = this.cc.getAntennaId(1:max(present_sat)); % sat in network are equals to go_id
                    s2syc_c = a_id(:,1);
                    sys_c = unique(s2syc_c);

%                     if ~isempty(this.wl_amb) && (n_rec > 2 || numel(sys_c > 1));
%                         % estimate narrowlanes phase delays and remove them
%                         % from abiguity vector and from observations
%                         for r = 2 : n_rec
%                             for j = 1:numel(sys_c)
%                                 s = sys_c(j);
%                             id_amb_r = amb_rec == r & s2syc_c(amb_sat) == sys_c(j);
%                             if sum(id_amb_r) > 0
%                                 weigth = min(n_ep_amb(id_amb_r),100)./100;
%                                 weigth = weigth./sum(weigth);
%                                 [~, frac_bias] = Core_Utils.getFracBias(amb(id_amb_r), weigth);
%                                 amb(id_amb_r) = amb(id_amb_r) - frac_bias;
%                             end
%                             end
%                         end
%                     end
                    
                    
                    % ILS shrinking, method 1
                    [amb_fixed, is_fixed, l_fixed] = Fixer.fix(amb, C_amb_amb, Main_Settings.NET_AMB_FIX_FIXER_APPROACH{this.state.getAmbFixNET});
                    
                    if is_fixed
                        % FIXED!!!!
                        amb_fix = amb_fixed(:, 1);
                        
                        for i = 1 : n_amb
                            Ni = N(:, idx_amb_par(i));
                            if l_fixed(i)
                                B = B - Ni * amb_fix(i);
                            end
                        end
                        
                        % remove fixed ambiguity from B and N
                        B(idx_amb_par(l_fixed(:,1))) = [];
                        N(idx_amb_par(l_fixed(:,1)), :) = [];
                        N(:, idx_amb_par(l_fixed(:,1))) = [];
                        
                        % recompute the solution
                        idx_nf = true(sum(idx_est,1),1);
                        idx_nf(idx_amb_par(l_fixed(:,1))) = false;
                        xf = zeros(size(idx_nf));
                        
                        xf(idx_nf) = N \ B;
                        xf(~idx_nf) = amb_fix(l_fixed(:,1));
                        x(idx_est) = xf;
                    else
                        this.log.addWarning('The ambiguities cannot be fixed!!!');
                    end
                end
            end
            
            if nargout > 1
                x_res = zeros(size(x));
                if not(isempty(x_res))
                    x_res(N2A_idx) = x(1:end-size(this.G,1));
                end
                if sum(isnan(x_res)) ==0
                    [res, av_res] = this.getResiduals(x_res);
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
                x = [x, x_rec];
            end
        end
        
        function [x, res, s0, Cxx, l_fixed, av_res] = solvePPP(this)
            av_res = [];
            l_fixed = 0;
            Cxx = [];
            % if N_ep if empty call A
            if isempty(this.N_ep)
                this.Astack2Nstack();
            end
            is_network = this.network_solution;
            u_r = unique(this.receiver_id);
            N = [];
            B = [];
            N2A_idx = [];
            A2N_idx_tot = [];
            x_rec = [];
            for r = u_r(:)'
                idx_r_l = this.receiver_id == u_r(r);
                idx_r = find(idx_r_l);
                A_rec = this.A_idx(idx_r_l, :);
                idx_constant_l = this.param_flag(:) == 0 | this.param_flag(:) == -1;
                idx_constant = find(idx_constant_l);
                idx_non_constant = find(~idx_constant_l);
                idx_ep_wise = this.param_flag == 1;
                u_class = unique(this.param_class);
                n_param_class = zeros(1,length(u_class));
                for u = 1 : length(u_class)
                    n_param_class(u) = max(max(this.A_idx(:,this.param_class == u_class(u)))) - min(min(this.A_idx(:,this.param_class == u_class(u)))) + 1;
                end
                a_idx_const = unique(A_rec(:, idx_constant_l));
                a_idx_const(a_idx_const == 0) = [];
                a_idx_ep_wise = unique(A_rec(:, idx_ep_wise));
                a_idx_ep_wise(a_idx_ep_wise == 0) = [];
                n_constant = length(a_idx_const);
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
                    p_c_idx = A2N_idx(this.A_idx(i, idx_constant_l));
                    p_e_idx = A2N_idx(max(1, this.A_idx(i, ~idx_constant_l))) - n_constant;
                    
                    % fill Ncc
                    Ncc(p_c_idx, p_c_idx) = Ncc(p_c_idx, p_c_idx) + this.N_ep(idx_constant, idx_constant, i);
                    % fill Nce
                    Nce(p_e_idx, p_c_idx) = Nce(p_e_idx, p_c_idx) + this.N_ep(idx_non_constant, idx_constant, i);
                    %fill Ndiags
                    
                    Ndiags(:, :, r2p_ep(this.epoch(i))) = Ndiags(:, :, r2p_ep(this.epoch(i))) + this.N_ep(idx_non_constant, idx_non_constant, i);
                    %fill B
                    b_id = A2N_idx(this.A_idx(i,:));
                    B_rec(b_id) = B_rec(b_id) + this.A_ep(i, :)' * (1 ./ this.variance(i)) * this.rw(i) * this.y(i);
                end
                % constraint over rate paramters
                if ~isempty(this.time_regularization)
                    cc_class = intersect(this.time_regularization(:, 1),this.param_class(idx_constant));
                    if ~isempty(cc_class)
                        Ncc_reg = sparse(n_constant, n_constant);
                        diag_cc0 = zeros(n_constant,1);
                        diag_cc1 = zeros(n_constant-1,1);
                        idx_const_u = false(size(n_param_class));
                        for b = 1 : length(u_class)
                            idx_const_u(b) = unique(idx_constant_l(this.param_class == u_class(b)));
                        end
                        n_param_class_const_cum = cumsum(n_param_class(idx_const_u));
                        for c = cc_class'
                            idx_c = this.time_regularization(:, 1) == c;
                            w = 1 ./ this.time_regularization(idx_c, 2) ;
                            reg = ones(n_param_class(u_class == c)-1,1)*w;
                            idx_bg = n_param_class_const_cum(find(u_class(idx_const_u) == c)-1)+1;
                            idx_en = n_param_class_const_cum(u_class(idx_const_u) == c)-1;
                            diag_cc1(idx_bg:idx_en) = -reg;
                            reg0 = [0; reg] + [reg; 0];
                            diag_cc0(idx_bg:(idx_en+1)) = reg0;
                        end
                        Ncc_reg = spdiags( diag_cc0,0,Ncc_reg);
                        Ncc_reg = spdiags([0; diag_cc1],1,Ncc_reg);
                        Ncc_reg = spdiags(diag_cc1,-1,Ncc_reg);
                        Ncc = Ncc + Ncc_reg;
                    end
                    
                end
                
                Nee = [];
                class_ep_wise = this.param_class(idx_non_constant);
                
                diff_reg = 1./double(diff(this.true_epoch(u_ep)));
                reg_diag0 = [diff_reg; 0 ] + [0; diff_reg];
                reg_diag1 = -diff_reg ;
                Ndiags = permute(Ndiags, [3, 1, 2]);
                tik_reg = ones(n_epochs,1)/n_epochs; %%% Tikhonov on ZTD and gradients
                % build the epoch wise parameters
                for i = 1 : n_ep_class
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
                x_rec = [x_rec; r*ones(size(N_rec,1),1)];
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
            
            warning off
            if (size(N, 1) > 1e5)
                log = Core.getLogger();
                log.addMessage(log.indent(sprintf('The full system is sparse %d x %d\nplease wait for the solver to end...', size(N, 1), size(N, 2))));
                ts = tic();
            end
            x = N \ B;
            if (size(N, 1) > 1e5)
                toc(ts)
            end
            warning on
            
            x_class = zeros(size(x));
            for c = 1:length(this.param_class)
                idx_pars = this.A_idx(:, c);
                idx_p = A2N_idx_tot(idx_pars(idx_pars~=0));
                x_class(idx_p) = this.param_class(c);
            end
            
            % Ambiguity fixing
            if (this.state.getAmbFixPPP && ~isempty(x(x_class == this.PAR_AMB,1)))
                amb = x(x_class == this.PAR_AMB,1);
                amb_wl_fixed = false(size(amb));
                amb_n1 = nan(size(amb));
                amb_wl = nan(size(amb));
                n_ep_wl = zeros(size(amb));
                n_amb = max(max(this.amb_idx));
                n_ep = size(this.wl_amb,1);
                if sum(this.param_class == this.PAR_X) >0
                    n_coo = max(this.A_idx(:,3));
                else
                    n_coo = 0;
                end
                for i = 1 : n_amb
                    sat = this.sat_go_id(this.sat(this.A_idx(:,this.param_class == this.PAR_AMB) == i+n_coo));
                    idx = n_ep*(sat(1)-1) +  this.true_epoch(this.epoch(this.A_idx(:,this.param_class == this.PAR_AMB)== i+n_coo));
                    if ~isempty( this.wl_fixed) % case local single frequency PPP
                        amb_wl(i) = this.wl_amb(idx(1));
                        amb_wl_fixed(i)=  this.wl_fixed(idx(1));
                    end
                    n_ep_wl(i) = length(idx);
                    amb_n1(i) = amb(i); %(amb(i)- 0*f_vec(2)^2*l_vec(2)/(f_vec(1)^2 - f_vec(2)^2)* wl_amb);  % Blewitt 1989 eq(23)
                    
                end
                
                weight = min(n_ep_wl(amb_wl_fixed),100); % <- downweight too short arc
                weight = weight / sum(weight);
                
                
                idx_amb = find(x_class == this.PAR_AMB);
                % get thc cxx of the ambiguities
                n_amb  = length(idx_amb);
                b_eye = zeros(length(B),n_amb);
                idx = sub2ind(size(b_eye),idx_amb,[1:n_amb]');
                b_eye(idx) = 1;
                b_eye = sparse(b_eye);
                Cxx_amb = N\b_eye;
                Cxx_amb = Cxx_amb(idx_amb,:);
                idx_constarined = abs(amb_n1) < 1e-5;
                l_fixed = false(size(amb_n1,1),1);
                amb_fixed = zeros(size(amb_n1,1),1);
                % ILS shrinking, method 1
                [af, is_fixed,lf] = Fixer.fix(amb_n1(~idx_constarined), Cxx_amb(~idx_constarined,~idx_constarined), 'sequential_best_integer_equivariant');
                
                if is_fixed
                    amb_fixed(~idx_constarined,:) = af;
                    l_fixed(~idx_constarined,:) = lf;
                    % FIXED!!!!
                    idx_est = true(size(x,1),1);
                    amb_fix = amb_fixed(:, 1);
                    
                    for i = 1 : n_amb
                        Ni = N(:, idx_amb(i));
                        if l_fixed(i)
                            B = B - Ni * amb_fix(i);
                        end
                    end
                    
                    % remove fixed ambiguity from B and N
                    B(idx_amb(l_fixed(:,1))) = [];
                    N(idx_amb(l_fixed(:,1)), :) = [];
                    N(:, idx_amb(l_fixed(:,1))) = [];
                    
                    % recompute the solution
                    idx_nf = true(sum(idx_est,1),1);
                    idx_nf(idx_amb(l_fixed(:,1))) = false;
                    xf = zeros(size(idx_nf));
                    
                    xf(idx_nf) = N \ B;
                    xf(~idx_nf) = amb_fix(l_fixed(:,1));
                    
                    x(idx_est) = xf;
                else
                    this.log.addWarning('The ambiguities cannot be fixed!!!');
                end
            end
            
            if nargout > 1
                x_res = zeros(size(x));
                if not(isempty(x_res))
                    x_res(N2A_idx) = x(1:end-size(this.G,1));
                end
                if sum(isnan(x_res)) ==0
                    [res,av_res] = this.getResiduals(x_res);
                    s0 = mean(abs(res(res~=0)));
                else
                    res = [];
                    s0 = Inf;
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
        
        function remObs(this, idx_obs)
            % remove observations from the system
            %
            % SYNTAX:
            %   this.remObs(idx_obs)
            param_to_el = unique(this.A_idx(idx_obs,:));
            epoch_to_el = unique(this.epoch(idx_obs));
            
            if sum(this.param_class == this.PAR_AMB) > 0
                inearInd = sub2ind(size(this.amb_idx), this.epoch(idx_obs), this.sat(idx_obs));
                amb_to_el = unique(this.amb_idx(inearInd));
                this.amb_idx(inearInd) = 0;
            end
            
            % remove lines from the design matrix
            this.A_idx(idx_obs,:) = [];
            this.A_ep(idx_obs,:) = [];
            this.N_ep(:,:,idx_obs) = [];
            this.variance(idx_obs) = [];
            this.epoch(idx_obs) = [];
            this.sat(idx_obs) = [];
            this.receiver_id(idx_obs) = [];
            this.rw(idx_obs) = [];
            this.res(idx_obs) = [];
            this.y(idx_obs) = [];
                        
            param_actual = unique(this.A_idx);
            %change paramter indexes
            el_count = 0;
            for i = 1: length(param_to_el)
                if sum(param_to_el(i) == param_actual) == 0
                    idx_maj = this.A_idx > param_to_el(i);
                    this.A_idx(idx_maj) = this.A_idx(idx_maj) - 1;
                    param_actual(param_actual > param_to_el(i)) = param_actual(param_actual > param_to_el(i)) -1;
                    if sum(this.param_class == this.PAR_AMB) > 0
                        % this matrix exist only when ambiguities are present
                        this.G(:,param_to_el(i)) = [];
                    end
                    param_to_el = param_to_el -1;
                end
            end
            for i = 1: length(epoch_to_el)
                if sum(epoch_to_el(i) == this.epoch) ==0
                    idx_maj = this.epoch > epoch_to_el(i);
                    this.epoch(idx_maj) = this.epoch(idx_maj) - 1;
                    idx_maj = this.system_split > epoch_to_el(i);
                    this.system_split(idx_maj) = this.system_split(idx_maj) - 1;
                    this.true_epoch(epoch_to_el(i)) = [];
                    this.amb_idx(epoch_to_el(i) - el_count,:) = [];
                    epoch_to_el = epoch_to_el -1;
                end
            end
            this.n_epochs = length(unique(this.epoch));
            % If I have ambiguities in the receiver
            if sum(this.param_class == this.PAR_AMB) > 0                
                amb_actual = unique(this.amb_idx);
                %change paramter indexes
                for i = 1: length(amb_to_el)
                    if sum(amb_to_el(i) == amb_actual) ==0
                        idx_maj = this.amb_idx > amb_to_el(i);
                        this.amb_idx(idx_maj) = this.amb_idx(idx_maj) - 1;
                        amb_actual(amb_actual > amb_to_el(i)) = amb_actual(amb_actual > amb_to_el(i)) -1;
                        amb_to_el = amb_to_el -1;
                    end
                end
            end
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
            % given a time and the sampling rate return the position index referring to the given sampling rate, the first index is prgressive, the seocond id time consistent
            sec_from_sod = time.getRefTime(st_time.getMatlabTime);
            pos_idx_tc = max(1,ceil((sec_from_sod - 0.002) / coo_rate));
            u_pos = unique(pos_idx_tc);
            pos_idx_nh = pos_idx_tc;
            for i = 1 : length(u_pos)
                pos_idx_nh(pos_idx_tc==u_pos(i)) = i;
            end
            pos_idx_tc = unique(pos_idx_tc);
            pos_idx_tc = pos_idx_tc - (pos_idx_tc(1) - 1);
        end
    end
end
