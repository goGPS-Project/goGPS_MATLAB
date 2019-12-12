%   CLASS Network
% =========================================================================
%
% DESCRIPTION
%   Class to manage network processing of receivers
%
% EXAMPLE
%   Net = Network();
%
% SEE ALSO
% FOR A LIST OF CONSTANTs and METHODS use doc Network

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 4 ION
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
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
classdef Network < handle
    properties
        rec_list
        state
        cc
        net_id           % id of the receiver in core, this uniquely identify the network
        
        common_time      % gps_time
        rec_time_indexes % indexes
        coo              % [n_coo x n_rec x n_sol] receiver coordinates
        coo_vcv          % [6 x n_rec x n_sol] receiver coordinates
        coo_rate         % rate of the coordinate solution in seconds
        clock            % [n_epoch x n_rec] reciever clock
        ztd              % [n_epoch x n_rec] reciever ZTD
        ztd_gn           % [n_epoch x n_rec] reciever ZTD gradients north
        ztd_ge           % [n_epoch x n_rec] reciever ZTD gradients east
        amb              % {n_rec} recievers ambiguity
        log
        pos_indexs_tc    % index for subpositions
        central_coo      % index of the central coordinate
        id_ref
        wl_mats          % widelane matrices
        wl_comb_codes    % codes of the widelanes (e.g. G12 B27) these are the rinex3 codes
        tropo_idx        % index of the splien tropo
        tropo_g_idx      % index fo the spline tropo gradient
        sat_wb           % staellite widelane biases
        rec_eb           % electronic bias of the receiver
        sat_eb           % electronic bias of the satellites
        
        apriori_info     % field to keep apriori info [ambiguity, tropo, ...] to be used in the adjustment
        is_tropo_decorrel % are station apart enough to estimate differents tropo?
        is_coo_decorrel   % are station decorrelated enough
    end
    
    methods
        function this = Network(rec_list, net_id)
            if nargin < 2
                net_id = 1 : length(rec_list);
            end
            this.net_id = net_id;
            this.rec_list = rec_list;
            this.init();
        end
        
        function init(this)
            this.state = Core.getState;
            this.cc = this.state.cc;
            this.log = Core.getLogger();
        end
        
        function reset(this)
            % clear the object keeping only its id and apriori info and the receivers
            %
            % SYNTAX
            %    this.reset()
            this.common_time = [];
            this.rec_time_indexes = [];
            this.coo = [];
            this.coo_vcv = [];
            this.coo_rate = [];
            this.clock = [];
            this.ztd = [];
            this.ztd_gn = [];
            this.ztd_ge = [];
            this.amb = [];
            this.pos_indexs_tc = [];
            this.id_ref = [];
            this.wl_mats  = [];
            this.wl_comb_codes = [];
        end
        
        function adjust(this, id_ref, coo_rate, reduce_iono, export_clk, frequency_idx, free_net)
            % Adjust the GNSS network
            %
            % INPUT
            %     id_ref : [1,n_rec]  receivers numeric index to be choosen as reference, their value mean will be set to zero
            %   coo_rate : rate of the solution
            %reduce_iono : reduce for ionospheric delay
            %
            % SYNATAX
            %    this. adjustNetwork(id_ref, <coo_rate>, <reduce_iono>)
            if nargin < 3
                coo_rate = [];
            end
            if nargin < 4
                reduce_iono = false;
            end
            if nargin < 5
                export_clk = false;
            end
            
            if nargin < 7
                free_net = false;
            end
            
            
            
            if nargin < 2 || any(isnan(id_ref)) || isempty(id_ref)
                lid_ref = true(size(this.net_id));
            else
                % convert to logical
                lid_ref = false(numel(this.rec_list),1);
                [~, idx_ref] = intersect(this.net_id, id_ref);
                lid_ref(idx_ref) = true;
            end
            is_empty_recs = ~this.rec_list.hasPhases_mr;
            n_valid_rec = sum(~is_empty_recs);
            n_valid_ref = sum(~is_empty_recs(lid_ref));
            if n_valid_ref < numel(id_ref)
                this.log.addError('One or more reference stations for Network solution are missing! Skipping NET');
            elseif (n_valid_rec < 2)
                this.log.addError('Not enough receivers (< 2), skipping network solution');
            else
                % if iono reduction is requested take off single frequency
                % receiver
                if reduce_iono
                    r = 1;
                    while (r <= length(this.rec_list))
                        if ~this.rec_list(r).work.isMultiFreq
                            this.log.addWarning(sprintf('Receiver %s is not multi frequency, removing it from network processing.',this.rec_list(r).getMarkerName4Ch));
                            this.rec_list(r) = [];
                            id_ref(id_ref == r) = [];
                            id_ref(id_ref > r) = id_ref(id_ref > r) -1;
                        else
                            r = r +1;
                        end
                    end
                end
                
                % set up the the network adjustment
                if nargin < 2 || any(isnan(id_ref)) || isempty(id_ref)
                    lid_ref = true(size(this.net_id));
                else
                    % convert to logical
                    lid_ref = false(numel(this.rec_list),1);
                    [~, id_ref] = intersect(this.net_id, id_ref);
                    lid_ref(id_ref) = true;
                end
                l_fixed = 0; % nothing is fixed
                
                e = find(is_empty_recs);
                if ~isempty(e)
                    this.rec_list(e) = [];
                    lid_ref(e) = [];
                end
                id_ref = find(lid_ref);
                this.id_ref = id_ref;
                
                if this.state.getReweightNET() < 2
                    n_clean = 0;
                else
                    this.log.addMessage(this.log.indent('Network solution performing 4 loops of outlier detection on the residuals'), 2);
                    n_clean = 3;
                end
                
                while n_clean >= 0
                    ls = LS_Manipulator();
                    
                    if this.state.flag_amb_pass && this.state.getCurSession > 1 && ~isempty(this.apriori_info)
                        f_time = GPS_Time.now();
                        f_time.addSeconds(1e9);
                        rate = 0;
                        for i = 1 : length(this.rec_list)
                            [~,limc] = this.state.getSessionLimits();
                            idx_fv = find(this.rec_list(i).work.time >= limc.first,1,'first');
                            f_time = min(this.rec_list(i).work.time.getEpoch(idx_fv), f_time);
                            rate = max(rate,this.rec_list(i).work.time.getRate);
                        end
                        t_dist = f_time - this.apriori_info.epoch;
                        if t_dist > -0.02 && t_dist <= (rate + 0.02)
                            ls.apriori_info = this.apriori_info;
                        end
                    end
                    wl_struct = [];
                    if reduce_iono
                        if this.state.getAmbFixNET() > 1
                            this.estimateWB();
                            ls.wl_amb = this.wl_mats;
                        end
                        wl_struct = struct();
                        wl_struct.amb_mats = this.wl_mats;
                        wl_struct.combination_codes = this.wl_comb_codes;
                    end
                    [this.common_time, this.rec_time_indexes]  = ls.setUpNetworkAdj(this.rec_list, coo_rate, wl_struct,frequency_idx);
                    
                    
                    this.is_tropo_decorrel = this.state.isReferenceTropoEnabled;
                    this.is_coo_decorrel = free_net;
                    if isempty(this.rec_time_indexes)
                        return
                    else
                        n_time = this.common_time.length;
                        n_rec = length(this.rec_list);
                        rate = this.common_time.getRate();
                        ls.setTimeRegularization(ls.PAR_REC_CLK, this.state.std_clock* rate); % really small regularization
                        if this.state.flag_tropo
                            if this.state.spline_tropo_order == 0
                                ls.setTimeRegularization(ls.PAR_TROPO, (this.state.std_tropo)^2 / 3600 * rate );% this.state.std_tropo / 3600 * rate  );
                            else
                                ls.setTimeRegularization(ls.PAR_TROPO, (this.state.std_tropo)^2 / 3600 * this.state.spline_rate_tropo );% this.state.std_tropo / 3600 * rate  );
                            end
                        end
                        if this.state.flag_tropo_gradient
                            if this.state.spline_tropo_gradient_order == 0
                                ls.setTimeRegularization(ls.PAR_TROPO_N, (this.state.std_tropo_gradient)^2 / 3600 * rate );%this.state.std_tropo / 3600 * rate );
                                ls.setTimeRegularization(ls.PAR_TROPO_E, (this.state.std_tropo_gradient)^2 / 3600 * rate );%this.state.std_tropo  / 3600 * rate );
                            else
                                ls.setTimeRegularization(ls.PAR_TROPO_N, (this.state.std_tropo_gradient)^2 / 3600 * this.state.spline_rate_tropo_gradient ); %this.state.std_tropo / 3600 * rate );
                                ls.setTimeRegularization(ls.PAR_TROPO_E, (this.state.std_tropo_gradient)^2 / 3600 *  this.state.spline_rate_tropo_gradient); %this.state.std_tropo  / 3600 * rate );
                            end
                        end
                        if Core.isGReD
                            dist_matrix = zeros(n_rec);
                            for i = 1:n_rec
                                for j = i+1:n_rec
                                    coo1 = Coordinates(); coo1.xyz = this.rec_list(i).work.xyz;
                                    coo2 = Coordinates(); coo2.xyz = this.rec_list(j).work.xyz;
                                    dist_matrix(i,j) = coo1.ellDistanceTo(coo2);
                                end
                            end
                            state = Core.getCurrentSettings();
                            ls.dist_matr = dist_matrix;
                            ls.distance_regularization.fun_tropo = @(d)  state.tropo_spatial_reg_sigma.*2.^(-d./state.tropo_spatial_reg_d_distance);
                            ls.distance_regularization.fun_gradients = @(d)  state.tropo_gradient_spatial_reg_sigma.*2.^(-d./state.tropo_gradient_spatial_reg_d_distance);
                            
                        end
                        ls.is_tropo_decorrel = this.is_tropo_decorrel;
                        ls.is_coo_decorrel = this.is_coo_decorrel;
                        [x, res, s0, Cxx, l_fixed, av_res] = ls.solve;
                        this.tropo_idx = ls.tropo_idx;
                        this.tropo_g_idx = ls.tropo_g_idx;
                        %[x, res] = ls.solve;
                        %res = res(any(res(:,:,2)'), :, :);
                        if Core.isGReD && export_clk
                            GReD_Utility.substituteClK(av_res, ls.time);
                            if Core.getCurrentCore.state.net_amb_fix_approach > 1
                                GReD_Utility.substituteWSB(this.sat_wb', this.common_time.getCentralTime);
                            end
                        end
                        
                        % cleaning -----------------------------------------------------------------------------------------------------------------------------
                        %s0 = mean(abs(res(res~=0)));
                        this.log.addMessage(this.log.indent(sprintf('Network solution computed,  s0 = %.4f', s0)));
                        %figure; plot(res(:,:,2)); ylim([-0.05 0.05]); dockAllFigures;
                        if s0 < 0.0025 && ... % low sigma0 of the computation
                                ~(any(abs(mean(zero2nan(reshape(res, size(res,1), size(res,2) * size(res,3), 1)), 'omitnan')) > 2e-3) && this.state.getReweightNET() == 3) && ... % mean of residuals all below 5 mm
                                all(abs(res(:)) < this.state.getMaxPhaseErrThr) && ...  % no residuals above thr level
                                all(std(zero2nan(reshape(permute(res(:,:,:),[1 3 2]), size(res,1) * size(res,3), size(res,2))),1,2,'omitnan') < 9e-3) % low dispersion of the residuals
                            
                            % I can be satisfied
                            n_clean = -1;
                        end
                        out_found = 0;
                        while (out_found == 0) && (n_clean >= 0)
                            for r = 1 : length(this.rec_list)
                                tmp_work = this.rec_list(r).work;
                                
                                % Index of receiver obs present in res
                                id_sync_res = round((round(tmp_work.time.getRefTime(this.common_time.first.getMatlabTime) * 1e5) * 1e-5 /  this.common_time.getRate) + 1);
                                id_sync_rec = (1 : numel(id_sync_res))';
                                id_ok = id_sync_res > 0 & id_sync_res <= size(res,1);
                                id_sync_res = id_sync_res(id_ok);
                                id_sync_rec = id_sync_rec(id_ok);
                                
                                go_id = tmp_work.go_id(this.rec_list(r).work.findObservableByFlag('L'));
                                
                                res_rec = zeros(size(tmp_work.sat.outliers_ph_by_ph));
                                res_rec(id_sync_rec, :) = res(id_sync_res, go_id, r);
                                
                                if s0 > 1
                                    % The solution is really bad, something bad happened
                                    % Extreme experiment to recover a valid positioning
                                    id_ko = abs(mean(zero2nan(res_rec), 'omitnan')) > 1;
                                    id_ko = sparse(repmat(id_ko, size(res_rec, 1), 1));
                                else
                                    % Search for outliers until found or n_clean == 1 (last clean)
                                    [id_ko] = tmp_work.search4outliers(res_rec, n_clean + 1);
                                    
                                    % At the first loops  (3-2) just start filtering bad satellites
                                    % bad satellites are satellites with std > 15mm or std > max(10mm, 90% other rec)
                                    if n_clean > 1
                                        bad_sat = std(res_rec, 'omitnan') > max(0.01, perc(noNaN(std(res_rec, 'omitnan')), 0.9)) | max(abs(res_rec)) > 0.15;
                                        if any(bad_sat) && any(any(id_ko(:, bad_sat)))
                                            id_ko(:, ~bad_sat) = false;
                                        end
                                    end
                                    
                                    % Add satellites with very bad res mean
                                    if n_clean < 3 && this.state.getReweightNET() == 3
                                        % bad satellites are satellites with median > 5 mm or median > std(other rec median)
                                        % this check is also performend depending on the std of the residuals
                                        % on bad receivers should be less stringent
                                        mean_sat_res = median(zero2nan(res_rec), 'omitnan');
                                        bad_sat = abs(mean_sat_res) > max(2e-3, 2*std(mean_sat_res, 'omitnan')) & ...
                                            abs(mean_sat_res) > std(zero2nan(res_rec(:)), 'omitnan') & ...
                                            (abs(mean_sat_res) == max(abs(mean_sat_res))); % flag bad as the only the satellite
                                        if any(bad_sat) % && any(any(id_ko(:, bad_sat)))
                                            %  checking where the mean is above threshold 2mm
                                            mean_lev = movmean(zero2nan(res_rec(:, bad_sat)),10, 'omitnan');
                                            id_ko(:, bad_sat) = abs(mean_lev) > 1e-3;
                                        end
                                    end
                                end
                                if any(id_ko(:))
                                    tmp_work.addOutliers(id_ko, true);
                                    out_found = out_found + 1;
                                end
                            end
                            n_clean = n_clean - 1;
                        end
                        % end of cleaning ----------------------------------------------------------------------------------------------------------------------
                        
                    end
                end
                
                if s0 < 0.05 || true
                    % initialize array for results
                    this.initOut(ls);
                    this.addAdjValues(x);
                    this.changeReferenceFrame(id_ref);
                    this.addAprValues();
                    if this.state.flag_coo_rate
                        % save old apriori values to be used later
                        coo_old = zeros(n_rec,3);
                        for i = 1 : n_rec
                            coo_old(i,:) =  this.rec_list(i).work.xyz;
                        end
                    end
                    this.pushBackInReceiver(s0, res, ls, l_fixed);
                    %%% from widelane l1 to l1 l2
                    if this.state.getAmbFixNET > 1 && false
                        this.pushBackAmbiguities(x(x(:,2) == ls.PAR_AMB,1),wl_struct,ls.amb_idx,ls.go_id_amb,ls.rec_time_idxes);
                    end
                else
                    this.log.addWarning(sprintf('s0 ( %.4f) too high! not updating the results',s0));
                end
                
                % pass ambiguity
                if this.state.flag_amb_pass
                    if s0 < 0.01
                        max_ep = max(ls.epoch);
                        id_amb = find(x(:,2) == ls.PAR_AMB);
                        x_float = ls.x_float;
                        this.apriori_info.amb_value = [];
                        this.apriori_info.freqs = [];
                        this.apriori_info.goids = [];
                        this.apriori_info.epoch = [];
                        this.apriori_info.receiver = [];
                        this.apriori_info.fixed = [];
                        this.apriori_info.Cambamb = [];
                        keep_id = [];
                        nnf = 1;
                        for i =1:length(id_amb)
                            a = id_amb(i);
                            a_idx = ls.A_idx_mix(:,ls.param_class == ls.PAR_AMB) == a;
                            ep_amb = ls.epoch(a_idx);
                            is_fixed = abs(fracFNI(x_float(a,1))) < eps(x_float(a,1));
                            if sum(ep_amb == max_ep) > 0  && (is_fixed || ls.Cxx_amb(nnf,nnf) > 0)% if ambugity belongs to alst epochs and values of the vcv matrix is acceptble
                                sat = ls.sat(a_idx);
                                this.apriori_info.amb_value = [this.apriori_info.amb_value; x_float(a,1)];
                                
                                this.apriori_info.goids = [this.apriori_info.goids; sat(1)];
                                this.apriori_info.epoch = this.common_time.last;
                                this.apriori_info.receiver = [this.apriori_info.receiver; x(a,3)];
                                this.apriori_info.fixed = [this.apriori_info.fixed; is_fixed ];
                                if ~this.apriori_info.fixed(end)
                                    keep_id = [keep_id nnf];
                                end
                            end
                            if ~is_fixed
                                nnf = nnf +1;
                            end
                        end
                        this.apriori_info.Cambamb = ls.Cxx_amb(keep_id,keep_id);
                    end
                end
                
                % additional coordinate rate
                if this.state.flag_coo_rate
                    [sss_lim] = this.state.getSessionLimits();
                    st_time = sss_lim.first;
                    for i = 1 : 3
                        if this.state.coo_rates(i) ~= 0
                            this.coo_rate = this.state.coo_rates(i);
                            ls.pos_indexs_tc = {};
                            for j = 1 : n_rec
                                [pos_idx_nh, pos_idx_tc] = LS_Manipulator.getPosIdx(this.common_time.getEpoch(~isnan(this.rec_time_indexes(:,j))), st_time, this.state.coo_rates(i));
                                ls.pos_indexs_tc{end+1} = pos_idx_tc; % to be used afterwards to push back postions
                                r2c = ~isnan(this.rec_time_indexes(:,j));
                                pos_idx = zeros(n_time,1);
                                pos_idx(r2c) = pos_idx_nh;
                                ls.changePosIdx(j, pos_idx);
                            end
                            [x, res] = ls.solve;
                            s0 = mean(abs(res(res~=0)));
                            this.log.addMessage(this.log.indent(sprintf('Network solution computed,  s0 = %.4f',s0)));
                            % intilaize array for results
                            this.initOut(ls);
                            if s0 < 0.02
                                this.addAdjValues(x);
                                this.changeReferenceFrame(id_ref);
                                for k = 1 : n_rec
                                    % for all paramter take the apriori in the receiver and sum the netwrok estimated correction
                                    n_coo_set = size(this.coo,3);
                                    this.coo(k,:,:) = this.coo(k,:,:) + repmat(coo_old(k,:),1,1,n_coo_set);
                                end
                            else
                                this.log.addWarning(sprintf('s0 ( %.4f) for sub coo solution ( %f s)too high! not updating the results ', s0, this.state.coo_rates(i)));
                            end
                            time_coo = st_time.getCopy;
                            time_coo.addSeconds([0 : this.state.coo_rates(i) :  (sss_lim.last - st_time)]);
                            this.pushBackSubCooInReceiver(time_coo, this.state.coo_rates(i));
                        end
                    end
                end
                for i = 1 : n_rec
                    %this.rec_list(i).work.pushResult();
                    this.rec_list(i).work.updateErrTropo();
                end
            end
            
        end
        
        function adjustNew(this, id_ref, coo_rate, reduce_iono, export_clk, free_net)
            % Adjust the GNSS network
            %
            % INPUT
            %     id_ref : [1,n_rec]  receivers numeric index to be choosen as reference, their value mean will be set to zero
            %   coo_rate : rate of the solution
            %reduce_iono : reduce for ionospheric delay
            %
            % SYNATAX
            %    this. adjustNetwork(id_ref, <coo_rate>, <reduce_iono>)
            if nargin < 3
                coo_rate = [];
            end
            if nargin < 4
                reduce_iono = false;
            end
            if nargin < 5
                export_clk = false;
            end
            
            if nargin < 6
                free_net = false;
            end
            
            
            
            if nargin < 2 || any(isnan(id_ref)) || isempty(id_ref)
                lid_ref = true(size(this.net_id));
            else
                % convert to logical
                lid_ref = false(numel(this.rec_list),1);
                [~, idx_ref] = intersect(this.net_id, id_ref);
                lid_ref(idx_ref) = true;
            end
            is_empty_recs = this.rec_list.isEmpty_mr;
            n_valid_rec = sum(~is_empty_recs);
            n_valid_ref = sum(~is_empty_recs(lid_ref));
            if n_valid_ref < numel(id_ref)
                this.log.addError('One or more reference stations for Network solution are missing! Skipping NET');
            elseif (n_valid_rec < 2)
                this.log.addError('Not enough receivers (< 2), skipping network solution');
            else
                % if iono reduction is requested take off single frequency
                % receiver
                if reduce_iono
                    r = 1;
                    while (r <= length(this.rec_list))
                        if ~this.rec_list(r).work.isMultiFreq
                            this.log.addWarning(sprintf('Receiver %s is not multi frequency, removing it from network processing.',this.rec_list(r).getMarkerName4Ch));
                            this.rec_list(r) = [];
                            id_ref(id_ref == r) = [];
                            id_ref(id_ref > r) = id_ref(id_ref > r) -1;
                        else
                            r = r +1;
                        end
                    end
                end
                
                % set up the the network adjustment
                if nargin < 2 || any(isnan(id_ref)) || isempty(id_ref)
                    lid_ref = true(size(this.net_id));
                else
                    % convert to logical
                    lid_ref = false(numel(this.rec_list),1);
                    [~, id_ref] = intersect(this.net_id, id_ref);
                    lid_ref(id_ref) = true;
                end
                l_fixed = 0; % nothing is fixed
                is_empty_recs = this.rec_list.isEmpty_mr;
                
                e = find(is_empty_recs);
                if ~isempty(e)
                    this.rec_list(e) = [];
                    lid_ref(e) = [];
                end
                id_ref = find(lid_ref);
                this.id_ref = id_ref;
                
                if this.state.getReweightNET() < 2
                    n_clean = 0;
                else
                    this.log.addMessage(this.log.indent('Network solution performing 4 loops of outlier detection on the residuals'), 2);
                    n_clean = 3;
                end
                
                ls = LS_Manipulator_new();
                param_selection = [ls.PAR_REC_X ;
                    ls.PAR_REC_Y;
                    ls.PAR_REC_Z;
                    %ls.PAR_REC_PPB;
                    ls.PAR_REC_EB;
                    %ls.PAR_SAT_PPB;
                    ls.PAR_SAT_EBFR;
                    ls.PAR_SAT_EB;
                    ls.PAR_AMB;
                    
                    ls.PAR_REC_CLK_PR;
                    ls.PAR_SAT_CLK_PR;
                    ls.PAR_REC_CLK_PH;
                    ls.PAR_SAT_CLK_PH;
                    ];
                 if reduce_iono
                     param_selection = [param_selection;
                     ls.PAR_IONO;];
                 end
                if this.state.flag_tropo
                    param_selection = [param_selection;
                        ls.PAR_TROPO;];
                end
                %                 param_selection = [param_selection;
                %                         ls.PAR_TROPO_V;];
                if this.state.flag_tropo_gradient
                    param_selection = [param_selection;
                        ls.PAR_TROPO_N;
                        ls.PAR_TROPO_E;];
                end
                
                
                glonass_r_sum = 0;
                for i = 1 : length(this.rec_list)
                    if sum(this.rec_list(i).work.system == 'R') > 0
                        glonass_r_sum = glonass_r_sum + 1;
                    end
                end
                if glonass_r_sum && false
                    param_selection = [param_selection;
                        ls.PAR_REC_EB_LIN;];
                end
                parametrization = LS_Parametrization();
                
                if ~reduce_iono
                    parametrization.iono(2) = LS_Parametrization.ALL_REC;
                end
                
                if this.state.spline_tropo_order == 1
                    parametrization.tropo(1)  = LS_Parametrization.SPLINE_LIN;
                    parametrization.tropo_opt.spline_rate = this.state.spline_rate_tropo;
                elseif this.state.spline_tropo_order == 3
                    parametrization.tropo(1)  = LS_Parametrization.SPLINE_CUB;
                    parametrization.tropo_opt.spline_rate = this.state.spline_rate_tropo;
                end
                
                if this.state.spline_tropo_gradient_order == 1
                    parametrization.tropo_e(1)  = LS_Parametrization.SPLINE_LIN;
                    parametrization.tropo_n(1)  = LS_Parametrization.SPLINE_LIN;
                    parametrization.tropo_e_opt.spline_rate = this.state.spline_rate_tropo_gradient;
                    parametrization.tropo_n_opt.spline_rate = this.state.spline_rate_tropo_gradient;
                elseif this.state.spline_tropo_gradient_order == 3
                    parametrization.tropo_e(1)  = LS_Parametrization.SPLINE_CUB;
                    parametrization.tropo_n(1)  = LS_Parametrization.SPLINE_CUB;
                    parametrization.tropo_e_opt.spline_rate = this.state.spline_rate_tropo_gradient;
                    parametrization.tropo_n_opt.spline_rate = this.state.spline_rate_tropo_gradient;
                end
                
                ls.setUpNET(this.rec_list, coo_rate, '???', param_selection, parametrization);
               
                if this.state.flag_free_net_tropo
                    ls.free_tropo = true;
                end
                
                this.is_tropo_decorrel = this.state.isReferenceTropoEnabled;
                this.is_coo_decorrel = free_net;
                
                
                if this.state.flag_tropo
                    ls.timeRegularization(ls.PAR_TROPO, (this.state.std_tropo)^2 / 3600);
                    
                end
                if this.state.flag_tropo_gradient
                    ls.timeRegularization(ls.PAR_TROPO_N, (this.state.std_tropo_gradient /10)^2 / 3600);
                    ls.timeRegularization(ls.PAR_TROPO_E, (this.state.std_tropo_gradient /10)^2 / 3600);
                end
                
                
                if Core.isGReD
                    % distance regularization to be set up
                    
                end
                this.common_time = ls.unique_time;
                ls.solve(Core.getState.net_amb_fix_approach >1);
%                 idx_fix = ls.class_par == ls.PAR_AMB;
%                 idx_fix(idx_fix) = abs(fracFNI(ls.x(idx_fix))) < 1e-9; % fixed ambiguoty
%                 ls.removeEstParam(idx_fix);
                ls.simpleSnoop(Core.getState.getMaxPhaseErrThr, Core.getState.getMaxCodeErrThr);
                ls.solve(Core.getState.net_amb_fix_approach >1);                
                s0 = mean(abs(ls.res(ls.phase_obs > 0 & ~ls.outlier_obs)));
                if s0 < 0.05
                    % initialize array for results
                    this.initOutNew(ls);
                    this.addAdjValuesNew(ls);
                    this.changeReferenceFrame(id_ref);
                    this.addAprValues();
                    this.pushBackInReceiverNew(ls);
                    %%% from widelane l1 to l1 l2
                    if this.state.getAmbFixNET > 1 && false
                        this.pushBackAmbiguities(x(x(:,2) == ls.PAR_AMB,1),wl_struct,ls.amb_idx,ls.go_id_amb,ls.rec_time_idxes);
                    end
                else
                    this.log.addWarning(sprintf('s0 ( %.4f) too high! not updating the results',s0));
                end
            end
            
        end
        
        
        
        function wl_struct = estimateWL(this)
            % Estimate widelane
            %
            % SYNTAX:
            %    wl_struct = estimateWL(this)
            this.estimateWB();
            wl_struct = struct();
            wl_struct.amb_mats = this.wl_mats;
            wl_struct.combination_codes = this.wl_comb_codes;
            
        end
        
        function phase2pseudoranges(this)
            % reosolve all DD ambiguity with respect to one receiver,
            % concept taken from AMBIZAP
            %
            % SYNTAX:
            % this.phase2pseudaranges
            
            % distance matrix
            n_r = length(this.rec_list);
            DM = zeros(n_r);
            s = [];
            t = [];
            weight = [];
            for i = 1 : n_r
                this.rec_list(i).work.updateCoordinates;
                for j = i+1:n_r
                    DM(i,j) = sphericalDistance(this.rec_list(i).work.lat,this.rec_list(i).work.lon,this.rec_list(j).work.lat,this.rec_list(j).work.lon);
                    s = [s i];
                    t = [t j];
                    weight = [weight DM(i,j)];
                end
            end
            G = graph(s,t,weight);
            [T,pred] = minspantree(G);
            T = table2array(T.Edges);
            % minimum spanning tree
            nodes_level = [T(1,1)];
            nodes_level_next = [];
            n_expl = 1;
            while n_expl < n_r
                for n = nodes_level
                    % find a_b l braches
                    node_brnch = [];
                    id_b1 = T(:,1) == n;
                    id_b2 = T(:,2) == n;
                    node_brnch = [node_brnch T(id_b1,2)'];
                    node_brnch = [node_brnch T(id_b2,1)'];
                    for b = node_brnch
                        % do baseline processing
                        this.log.addMessage(sprintf('Finxig ambiguities betwewn %s and %s',this.rec_list(n).getMarkerName4Ch, this.rec_list(b).getMarkerName4Ch));
                        net_tmp =  Network(this.rec_list([n b]));
                        net_tmp.adjust([],[],true);
                        % substsitute the ambiguities
                        n_expl = n_expl +1;
                    end
                    T(id_b1,:) = [];
                    T(id_b2,:) = [];
                    
                    nodes_level_next= [nodes_level_next node_brnch];
                end
                nodes_level = nodes_level_next;
                nodes_level_next = [];
            end
            
            
        end
        
        function initOut(this,ls)
            n_time = this.common_time.length;
            n_rec = length(this.rec_list);
            n_set_coo = length(ls.getCommonPosIdx);
            if this.state.isSepCooAtBoundaries
                n_set_coo = 1;
            end
            this.pos_indexs_tc = ls.pos_indexs_tc;
            this.central_coo = ls.central_coo;
            this.clock = zeros(n_time, n_rec);
            this.coo = nan(n_rec, 3, n_set_coo);
            this.coo_vcv = nan(n_rec, 6, n_set_coo);
            this.ztd = nan(n_time, n_rec);
            this.ztd_gn = nan(n_time, n_rec);
            this.ztd_ge = nan(n_time, n_rec);
        end
        
        function initOutNew(this,ls)
            n_time = ls.unique_time.length;
            n_rec = length(this.rec_list);
            n_set_coo = length(unique(ls.time_par(ls.class_par == ls.PAR_REC_X,1)));
            if this.state.isSepCooAtBoundaries
                n_set_coo = 1;
            end
            this.clock = zeros(n_time, n_rec);
            this.coo = nan(n_rec, 3, n_set_coo);
            this.coo_vcv = nan(n_rec, 6, n_set_coo);
            this.ztd = nan(n_time, n_rec);
            this.ztd_gn = nan(n_time, n_rec);
            this.ztd_ge = nan(n_time, n_rec);
        end
        
        function addAdjValues(this, x)
            n_rec = length(this.rec_list);
            % --- fill the correction values in the network
            for i = 1 : n_rec
                % if all value in the receiver are set to nan initilaize them to zero
                if sum(isnan(this.rec_list(i).work.ztd)) == length(this.rec_list(i).work.ztd)
                    this.rec_list(i).work.ztd(:) = 0;
                    this.rec_list(i).work.tge(:) = 0;
                    this.rec_list(i).work.tgn(:) = 0;
                end
                % for all paramter take the apriori in the receiver and sum the netwrok estimated correction
                idx_rec = x(:,3) == i;
                if i > 0 % coordiantes are always zero on first receiver
                    coo = [x(x(:,2) == 1 & idx_rec,1) x(x(:,2) == 2 & idx_rec,1) x(x(:,2) == 3 & idx_rec,1)];
                    if isempty(coo)
                        coo = [ 0 0 0];
                    end
                    if ~isempty(this.pos_indexs_tc) && ~this.state.isSepCooAtBoundaries
                        this.coo(i,:,this.pos_indexs_tc{i}) = nan2zero(this.coo(i,:,this.pos_indexs_tc{i})) + permute(coo, [3 2 1]);
                    else
                        if this.state.isSepCooAtBoundaries && numel(coo) > 3
                            this.coo(i,:) = nan2zero(this.coo(i,:)) + coo(this.central_coo(i),:);
                        else
                            this.coo(i,:) = nan2zero(this.coo(i,:)) + coo;
                        end
                    end
                else
                    this.coo(i,:) = nan2zero(this.coo(i,:));
                end
                clk = x(x(:,2) == LS_Manipulator.PAR_REC_CLK & idx_rec,1);
                this.clock(~isnan(this.rec_time_indexes(:,i)),i) = nan2zero(this.clock(~isnan(this.rec_time_indexes(:,i)),i)) + clk;
                
                if this.state.flag_tropo
                    if this.state.spline_rate_tropo ~= 0 && this.state.spline_tropo_order > 0
                        tropo_dt = rem(this.common_time.getRefTime - this.common_time.getRate,this.state.spline_rate_tropo)/(this.state.spline_rate_tropo);
                        tropo_idx = max(1, ceil((this.common_time.getRefTime(this.common_time.first.getMatlabTime))/this.state.spline_rate_tropo+eps));
                        
                        spline_base = Core_Utils.spline(tropo_dt, this.state.spline_tropo_order);
                        tropo = x(x(:,2) == LS_Manipulator.PAR_TROPO & idx_rec,1);
                        if max(tropo_idx) > (numel(tropo) - this.state.spline_tropo_order) % if the receiver time was shorter than the common one
                            tropo = [tropo; tropo(end) * ones(numel(tropo)-max(tropo_idx),1)];
                        end
                        ztd = sum(spline_base.*tropo(repmat(tropo_idx,1,this.state.spline_tropo_order+1)+repmat((0:this.state.spline_tropo_order),numel(tropo_idx),1)),2);
                        this.ztd(:,i) = nan2zero(this.ztd(:,i))  + ztd;
                    else
                        ztd = x(x(:,2) == LS_Manipulator.PAR_TROPO & idx_rec,1);
                        this.ztd(~isnan(this.rec_time_indexes(:,i)),i) = nan2zero(this.ztd(~isnan(this.rec_time_indexes(:,i)),i))  + ztd;
                    end
                end
                
                if this.state.flag_tropo_gradient
                    if this.state.spline_rate_tropo_gradient ~= 0 && this.state.spline_tropo_gradient_order > 0
                        tropo_dt = rem(this.common_time.getRefTime - this.common_time.getRate,this.state.spline_rate_tropo_gradient)/(this.state.spline_rate_tropo_gradient);
                        tropo_g_idx = max(1, ceil((this.common_time.getRefTime(this.common_time.first.getMatlabTime))/this.state.spline_rate_tropo_gradient + eps));
                        
                        spline_base = Core_Utils.spline(tropo_dt,this.state.spline_tropo_gradient_order);
                        gntropo = x(x(:,2) == LS_Manipulator.PAR_TROPO_N & idx_rec,1);
                        if max(tropo_g_idx) > numel(gntropo) - this.state.spline_tropo_gradient_order % if the receiver time was shorter than the common one
                            gntropo = [gntropo; gntropo(end)*ones(numel(gntropo)-max(tropo_g_idx),1)];
                            getropo = [getropo; getropo(end)*ones(numel(getropo)-max(tropo_g_idx),1)];
                        end
                        gntropo = sum(spline_base.*gntropo(repmat(tropo_g_idx,1,this.state.spline_tropo_gradient_order+1)+repmat((0:this.state.spline_tropo_gradient_order),numel(tropo_g_idx),1)),2);
                        this.ztd_gn(:,i) = nan2zero(this.ztd_gn(:,i))  + gntropo;
                        getropo = x(x(:,2) == LS_Manipulator.PAR_TROPO_E & idx_rec,1);
                        getropo = sum(spline_base.*getropo(repmat(tropo_g_idx,1,this.state.spline_tropo_gradient_order+1)+repmat((0:this.state.spline_tropo_gradient_order),numel(tropo_g_idx),1)),2);
                        this.ztd_ge(:,i) = nan2zero(this.ztd_ge(:,i))  + getropo;
                    else
                        gn = x(x(:,2) == LS_Manipulator.PAR_TROPO_N & idx_rec,1);
                        this.ztd_gn(~isnan(this.rec_time_indexes(:,i)),i) = nan2zero(this.ztd_gn(~isnan(this.rec_time_indexes(:,i)),i)) + gn;
                        
                        ge = x(x(:,2) == LS_Manipulator.PAR_TROPO_E & idx_rec,1);
                        this.ztd_ge(~isnan(this.rec_time_indexes(:,i)),i) = nan2zero(this.ztd_ge(~isnan(this.rec_time_indexes(:,i)),i)) + ge;
                    end
                end
            end
        end
        
        function addAdjValuesNew(this, ls)
            n_rec = length(this.rec_list);
            % --- fill the correction values in the network
            rec_vcv = ls.rec_par([find(ls.class_par == ls.PAR_REC_X); find(ls.class_par == ls.PAR_REC_Y); find(ls.class_par == ls.PAR_REC_Z)]);
            rec_vcv(rec_vcv == this.id_ref) = [];

            for i = 1 : n_rec
                % if all value in the receiver are set to nan initilaize them to zero
                if sum(isnan(this.rec_list(i).work.ztd)) == length(this.rec_list(i).work.ztd)
                    this.rec_list(i).work.ztd(:) = 0;
                    this.rec_list(i).work.tge(:) = 0;
                    this.rec_list(i).work.tgn(:) = 0;
                end
                % for all paramter take the apriori in the receiver and sum the netwrok estimated correction
                idx_rec = ls.rec_par == i;
                [~, int_lim] = this.state.getSessionLimits();
                
                if i > 0 % coordiantes are always zero on first receiver
                    coo_vcv = [];
                    if this.state.isSepCooAtBoundaries
                        % Push coordinates from LS object to rec
                        idx_x = ls.class_par == ls.PAR_REC_X & idx_rec;
                        if sum(idx_x) > 0
                            x_coo = ls.x(idx_x);
                            [x_coo_time1, x_coo_time2] = ls.getTimePar(idx_x);
                            idx_save = x_coo_time1 - int_lim.first > -5e-2 & x_coo_time2 - int_lim.last < 5e-2;
                            cox = mean(x_coo(idx_save));
                        else
                            cox = 0;
                        end
                        
                        idx_y = ls.class_par == ls.PAR_REC_Y & idx_rec;
                        if sum(idx_y) > 0
                            y_coo = ls.x(idx_y);
                            [y_coo_time1, y_coo_time2] = ls.getTimePar(idx_y);
                            idx_save = y_coo_time1 - int_lim.first > -5e-2 & y_coo_time2 - int_lim.last < 5e-2;
                            coy = mean(y_coo(idx_save));
                        end
                        
                        idx_z = ls.class_par == ls.PAR_REC_Z & idx_rec;
                        if sum(idx_z) > 0
                            z_coo = ls.x(idx_z);
                            [z_coo_time1, z_coo_time2] = ls.getTimePar(idx_z);
                            idx_save = z_coo_time1 - int_lim.first > -5e-2 & z_coo_time2 - int_lim.last < 5e-2;
                            coz = mean(z_coo(idx_save));
                            
                        else
                            coz = 0;
                        end
                        
                        coo = [cox coy coz];
                    else
                        if i ~= this.id_ref
                            coo_vcv = ls.coo_vcv(rec_vcv == i,rec_vcv == i);
                            if ~isempty(coo_vcv)
                                coo_vcv = [coo_vcv(1,1) (coo_vcv(1,2) + coo_vcv(2,1))/2  (coo_vcv(1,3) + coo_vcv(3,1))/2 coo_vcv(2,2) (coo_vcv(2,3) + coo_vcv(3,2))/2 coo_vcv(3,3)];
                            else
                                coo_vcv = zeros(1,6);
                            end
                        end
                        coo = [mean(ls.x( ls.class_par == ls.PAR_REC_X & idx_rec)) mean(ls.x(ls.class_par == ls.PAR_REC_Y & idx_rec)) mean(ls.x(ls.class_par == ls.PAR_REC_Z & idx_rec))];
                    end
                    
                    if isempty(coo)
                        coo = [ 0 0 0];
                    end
                    
                    this.coo(i,:) = nan2zero(this.coo(i,:)) + coo;
                    if ~isempty(coo_vcv)
                        this.coo_vcv(i,:) = coo_vcv;
                    end
                    
                else
                    this.coo(i,:) = nan2zero(this.coo(i,:));
                end
                idx_clk = ls.class_par == LS_Manipulator_new.PAR_REC_CLK & idx_rec;
                clk = ls.x(idx_clk);
                time_clk = ls.time_par(idx_clk);
                [~,idx_time_clk] = ismember(time_clk, this.common_time.getNominalTime.getRefTime(ls.time_min.getMatlabTime));
                this.clock(idx_time_clk,i) = nan2zero(this.clock(idx_time_clk,i)) + clk;
                
                if this.state.flag_tropo
                    if this.state.spline_rate_tropo ~= 0 && this.state.spline_tropo_order > 0
                        idx_trp = ls.class_par == LS_Manipulator_new.PAR_TROPO & idx_rec;
                        tropo = ls.x(idx_trp);
                        tropo_dt = rem(this.common_time.getNominalTime(ls.obs_rate) - ls.getTimePar(idx_trp).minimum, this.state.spline_rate_tropo)/ this.state.spline_rate_tropo;
                        tropo_idx = floor((this.common_time.getNominalTime(ls.obs_rate) - ls.getTimePar(idx_trp).minimum)/this.state.spline_rate_tropo);
                        [~,tropo_idx] = ismember(tropo_idx*this.state.spline_rate_tropo, ls.getTimePar(idx_trp).getNominalTime(ls.obs_rate).getRefTime(ls.getTimePar(idx_trp).minimum.getMatlabTime));
                        valid_ep = tropo_idx ~=0;
                        spline_base = Core_Utils.spline(tropo_dt,this.state.spline_tropo_order);
                        
                        ztd =sum(spline_base .* tropo(repmat(tropo_idx(valid_ep), 1, this.state.spline_tropo_order + 1) + repmat((0 : this.state.spline_tropo_order), numel(tropo_idx(valid_ep)), 1)), 2);
                        
                        ztd = sum(spline_base.*tropo(repmat(tropo_idx,1,this.state.spline_tropo_order+1)+repmat((0:this.state.spline_tropo_order),numel(tropo_idx),1)),2);
                        this.ztd(:,i) = nan2zero(this.ztd(:,i))  + ztd;
                    else
                        idx_trp = ls.class_par == LS_Manipulator_new.PAR_TROPO & idx_rec;
                        tropo = ls.x(idx_trp);
                        time_tropo = ls.time_par(idx_trp);
                        [~,idx_time_tropo] = ismember(time_tropo, this.common_time.getNominalTime(ls.obs_rate).getRefTime(ls.time_min.getMatlabTime));
                        this.ztd(idx_time_tropo,i) = nan2zero(this.clock(idx_time_tropo,i)) + tropo;
                    end
                end
                
                if this.state.flag_tropo_gradient
                    if this.state.spline_rate_tropo_gradient ~= 0 && this.state.spline_tropo_gradient_order > 0
                        idx_trp_n = ls.class_par == LS_Manipulator_new.PAR_TROPO & idx_rec;
                        tropo_n = ls.x(idx_trp_n);
                        idx_trp_e = ls.class_par == LS_Manipulator_new.PAR_TROPO & idx_rec;
                        
                        tropo_e = ls.x(idx_trp_e);
                        tropo_dt = rem(this.common_time.getNominalTime(ls.obs_rate) - ls.getTimePar(idx_trp_n).minimum, this.state.spline_rate_tropo)/ this.state.spline_rate_tropo;
                        tropo_idx = floor((this.common_time.getNominalTime(ls.obs_rate) - ls.getTimePar(idx_trp_n).minimum)/this.state.spline_rate_tropo);
                        [~,tropo_idx] = ismember(tropo_idx*this.state.spline_rate_tropo, ls.getTimePar(idx_trp_n).getNominalTime(ls.obs_rate).getRefTime(ls.getTimePar(idx_trp_n).minimum.getMatlabTime));
                        valid_ep = tropo_idx ~=0;
                        spline_base = Core_Utils.spline(tropo_dt,this.state.spline_tropo_order);
                        
                        tropo_n =sum(spline_base .* tropo_n(repmat(tropo_idx(valid_ep), 1, this.state.spline_tropo_order + 1) + repmat((0 : this.state.spline_tropo_order), numel(tropo_idx(valid_ep)), 1)), 2);
                        tropo_e =sum(spline_base .* tropo_e(repmat(tropo_idx(valid_ep), 1, this.state.spline_tropo_order + 1) + repmat((0 : this.state.spline_tropo_order), numel(tropo_idx(valid_ep)), 1)), 2);
                        this.ztd_gn(:,i) = nan2zero(this.ztd_gn(:,i))  + tropo_n;
                        this.ztd_ge(:,i) = nan2zero(this.ztd_ge(:,i))  + tropo_e;
                    else
                        idx_tropo_n = ls.class_par == LS_Manipulator_new.PAR_TROPO_N & idx_rec;
                        tropo_n = ls.x(idx_tropo_n);
                        time_tropo_n = ls.time_par(idx_tropo_n);
                        [~,idx_time_tropo_n] = ismember(time_tropo_n, this.common_time.getNominalTime(ls.obs_rate).getRefTime(ls.time_min.getMatlabTime));
                        this.ztd_gn(idx_time_tropo_n,i) = nan2zero(this.clock(idx_time_tropo_n,i)) + tropo_n;
                        
                        idx_tropo_e = ls.class_par == LS_Manipulator_new.PAR_TROPO_E & idx_rec;
                        tropo_e = ls.x(idx_tropo_e);
                        time_tropo_e = ls.time_par(idx_tropo_e);
                        [~,idx_time_tropo_e] = ismember(time_tropo_e, this.common_time.getNominalTime(ls.obs_rate).getRefTime(ls.time_min.getMatlabTime));
                        this.ztd_ge(idx_time_tropo_e,i) = nan2zero(this.clock(idx_time_tropo_e,i)) + tropo_n;
                    end
                end
                
                
            end
        end
        
        function changeReferenceFrame(this, id_ref)
            n_rec = length(this.rec_list);
            n_time = this.common_time.length;
            % ALL OF THIS MAKES NO SENSE TO ME (Andrea). Now it should ;) (Giulio)
            %--- transform the result in the desired free network
            
            if ~isnan(id_ref(1))
                
                S = zeros(n_rec);
                S(:, id_ref) = - 1 / numel(id_ref);
                S = S + eye(n_rec);  % < - this should be an S trasform but i am not sure
                % it is the paramter itself  the mean of the reference paramter
                % it is in matrix form so it can be used in the future for variance covariance matrix of the coordinates
                
                % Applying the S transform I obtain the corrections with respect to the reference
                for i = 1 : size(this.coo,3)
                    this.coo(:,1,i) = S * this.coo(:,1,i);
                    this.coo(:,2,i) = S * this.coo(:,2,i);
                    this.coo(:,3,i) = S * this.coo(:,3,i);
                end
                
                % apply the S transform to the epochwise parameters
                for i = 1 : n_time
                    id_present = ~isnan(this.clock(i,:));
                    id_ref_t = intersect(id_ref, find(id_present));
                    if isempty(id_ref_t)
                        S = nan;
                    else
                        n_rec_t = sum(id_present);
                        S = zeros(n_rec_t);
                        S(:,id_ref) = - 1 / numel(id_ref_t);
                        S = S + eye(n_rec_t);
                    end
                    % clock
                    this.clock(i,:) = (S*this.clock(i,:)')';
                    % ztd
                    if ~this.is_tropo_decorrel
                        if this.state.flag_tropo
                            this.ztd(i,:) = (S*this.ztd(i,:)')';
                        end
                        % gradients
                        if this.state.flag_tropo_gradient
                            this.ztd_gn(i,:) = (S*this.ztd_gn(i,:)')';
                            this.ztd_ge(i,:) = (S*this.ztd_ge(i,:)')';
                        end
                    end
                end
            end
        end
        
        function addAprValues(this)
            n_rec = length(this.rec_list);
            % --- add the apriori values
            for i = 1 : n_rec
                % if all value in the receiver are set to nan initilaize them to zero
                if sum(isnan(this.rec_list(i).work.ztd)) == length(this.rec_list(i).work.ztd)
                    this.rec_list(i).work.ztd(:) = 0;
                    this.rec_list(i).work.tge(:) = 0;
                    this.rec_list(i).work.tgn(:) = 0;
                end
                % for all paramter take the apriori in the receiver and sum the netwrok estimated correction
                n_coo_set = size(this.coo,3);
                this.coo(i,:,:) = this.coo(i,:,:) + repmat(this.rec_list(i).work.xyz,1,1,n_coo_set);
                %
                [idx_is, idx_pos] = ismembertol(this.rec_list(i).work.getTime.getGpsTime(), this.common_time.getGpsTime, 0.002, 'DataScale', 1);
                idx_pos = idx_pos(idx_pos > 0);
                clk_rec = this.rec_list(i).work.getDt();
                this.clock(idx_pos,i) = this.clock(idx_pos,i) + clk_rec(idx_is);
                
                if this.state.flag_tropo
                    ztd_rec = this.rec_list(i).work.getZtd();
                    ztd_rec_apr = this.rec_list(i).work.getZwd() + this.rec_list(i).work.getAprZhd();
                    ztd_rec(ztd_rec == 0) = ztd_rec_apr(ztd_rec == 0);
                    this.ztd(idx_pos,i) = this.ztd(idx_pos,i) + ztd_rec(idx_is);
                end
                
                if this.state.flag_tropo_gradient
                    [gn_rec, ge_rec] = this.rec_list(i).work.getGradient();
                    
                    this.ztd_gn(idx_pos,i) = this.ztd_gn(idx_pos,i) + gn_rec(idx_is);
                    
                    this.ztd_ge(idx_pos,i) = this.ztd_ge(idx_pos,i) + ge_rec(idx_is);
                end
            end
        end
        
        function estimateWB(this)
            % estimate widelane satellite bias e and widelane receiver bias
            % for a network
            % NOTE:
            %    the bias willd dpend on what code bias have been applied
            %    to observations so if code biases are changed
            %
            % SYNTAX:
            %      this.estimateWB()
            
            % firstly we define which widelane we seek to fix: once one
            % independet set of widelane has been fixed all widelane can be
            % fixed
            % for all constallation we seek to fix the widelana and
            % extrawidelane in the case of galileo and qzss are to be studied
            % NOTE: frequency order specified in [const]_SS class
            wide_laneM = {};
            wide_laneM{1} = [ 1 -1 0; %% GPS
                0  1 -1;];
            wide_laneM{2} = [ 1 -1 0; %% GLONASS
                0  1 -1;];
            wide_laneM{3} = [ 1 -1  0  0 0; %% GALILEO
                0  1 0  0 -1;
                0  1  -1 0 0;
                0  0  1 -1 0; ] ;
            wide_laneM{5} = [ 1 -1 0; %% BEIDOU
                0  1 -1;];
            wide_laneM{4} = [1 -1  0  0; % QZSS
                0  1 -1  0;
                0  0  1 -1] ;
            
            wide_laneM{6} = [1 -1]; %% IRNSS
            wide_laneM{7} = [1-1]; %% SBAS
            n_r = length(this.rec_list);
            sys_cs = this.cc.getAvailableSys;
            sys_cs(sys_cs == 'R') = []; %gloNASS excluded fro now from ambiguty fixing
            n_sat_tot = this.cc.getMaxNumSat();
            % intilaize a matrix to store the widelanes
            this.wl_mats = {};
            for r = 1 : n_r
                this.wl_mats{r} = nan(this.rec_list(r).work.time.length,n_sat_tot);
            end
            this.sat_wb = zeros(this.cc.getMaxNumSat(),1);
            for sys_c = sys_cs
                sys_idx = this.cc.SYS_C == sys_c;
                WM = wide_laneM{sys_idx}; %% get the right widelane matrix
                n_sat = this.cc.getNumSat(sys_c);
                for w = 1 %: size(WM,2)
                    sys_SS = this.cc.getSys(sys_c);
                    b1 = sys_SS.CODE_RIN3_2BAND(WM(w,:) == 1);
                    b2 = sys_SS.CODE_RIN3_2BAND(WM(w,:) == -1);
                    has_frs = false(n_r,1);
                    coderin3attr1 = sys_SS.CODE_RIN3_ATTRIB{WM(w,:) == 1};
                    coderin3attr2 = sys_SS.CODE_RIN3_ATTRIB{WM(w,:) == -1};
                    b1code_aval = false(n_r,n_sat, length(coderin3attr1));
                    b2code_aval = false(n_r,n_sat, length(coderin3attr2));
                    % matrix to know if the code have been aligned
                    aligned1 =  false(n_r,n_sat, length(coderin3attr2));
                    aligned2 =  false(n_r,n_sat, length(coderin3attr2));
                    for r = 1 : n_r
                        % check wether both frequency are available
                        rec = this.rec_list(r).work;
                        has_fr = sum(strLineMatch(rec.obs_code(:,1:2),['L' b1])) > 0 & sum(strLineMatch(rec.obs_code(:,1:2),['L' b2])) & sum(strLineMatch(rec.obs_code(:,1:2),['C' b1]))& sum(strLineMatch(rec.obs_code(:,1:2),['C' b2]));
                        has_frs(r) = has_fr;
                        if has_fr
                            % get all code for each reciever for each
                            % satellite for both frequency ordered by best noise (NOTE: there is
                            % a repeptiton here since frequency are checked
                            % multiple times, once everything works it is
                            % going to be inporved)
                            for s = 1:n_sat
                                idx_s = rec.obs_code(:,1) == 'C' & rec.prn == s & rec.system' == sys_c;
                                lid_sb1 = idx_s & rec.obs_code(:,2) == b1;
                                lid_sb2 = idx_s & rec.obs_code(:,2) == b2;
                                id_sb1 = find(lid_sb1);
                                id_sb2 = find(lid_sb2);
                                track_code_b1 = rec.obs_code(lid_sb1,3);
                                track_code_b2 = rec.obs_code(lid_sb2,3);
                                for t = 1 : length(track_code_b1)
                                    b1code_aval(r,s,coderin3attr1 == track_code_b1(t)) = true;
                                    aligned1(r,s,coderin3attr1 == track_code_b1(t)) = rec.aligned(id_sb1(t));
                                end
                                for t = 1 : length(track_code_b2)
                                    b2code_aval(r,s,coderin3attr2 == track_code_b2(t)) = true;
                                    aligned2(r,s,coderin3attr2 == track_code_b2(t)) = rec.aligned(id_sb2(t));
                                end
                            end
                        end
                    end
                    % remove satellites and receiver that doe not have the
                    % frequency
                    full_rec_lid = sum(sum(b1code_aval,3)> 0 & sum(b2code_aval,3)>0,2) > 0;
                    full_sat_lid = sum(sum(b1code_aval,3)> 0 & sum(b2code_aval,3)>0,1) > 0;
                    b1code_aval(~full_rec_lid,:,:) = [];
                    b2code_aval(~full_rec_lid,:,:) = [];
                    b1code_aval(:,~full_sat_lid,:) = [];
                    b2code_aval(:,~full_sat_lid,:) = [];
                    aligned1(~full_rec_lid,:,:) = [];
                    aligned2(~full_rec_lid,:,:) = [];
                    aligned1(:,~full_sat_lid,:) = [];
                    aligned2(:,~full_sat_lid,:) = [];
                    % remove code that are not observed by more than one
                    % receiver
                    nt2cb1 = sum(sum(squeeze(sum(b1code_aval,1)>1)),1)> 0;
                    nt2cb2 = sum(sum(squeeze(sum(b2code_aval,1)>1)),1)> 0;
                    selected_code_1 = find(nt2cb1, 1, 'first'); % best code avaliable to more than two receiver
                    selected_code_2 = find(nt2cb2, 1, 'first'); % best code avaliable to more than two receiver
                    track_1 = coderin3attr1(selected_code_1);
                    track_2 = coderin3attr2(selected_code_2);
                    if sum(nt2cb1) == 0 ||  sum(nt2cb1) == 0
                        this.log.addWarning(sprintf('No common code on the ferequency %s and freqeuncy %s for system %s',b1,b2,sys_c));
                    else
                        this.log.addMessage(this.log.indent(sprintf('Estimating %s%s%s widelane using tracking %s on frequency %s and tracking %s on frequency %s',sys_c,b1,b2,track_1,b1,track_2,b2)));
                        b1code_aval(:,:,~nt2cb1) = [];
                        b2code_aval(:,:,~nt2cb2) = [];
                        rec_with_no_cod1 = sum(sum(aligned1(:,:,:),3),2) == 0 & sum(b1code_aval(:,:,1),2) == 0;
                        rec_with_no_cod2 = sum(sum(aligned2(:,:,:),3),2) == 0 & sum(b2code_aval(:,:,1),2) == 0;
                        rec_excluded = rec_with_no_cod1 | rec_with_no_cod2;
                        if sum(rec_excluded) > 0
                            for r = serialize(find(rec_excluded))'
                                log_str = sprintf('Receiver %s excluded from widelane fixing because does not have:', this.rec_list(r).getMarkerName4Ch);
                                if rec_with_no_cod1(r)
                                    log_str = [log_str sprintf(' code %s on band %s', track_1,b1)];
                                end
                                
                                if rec_with_no_cod1(r) & rec_with_no_cod2(r)
                                    log_str = [log_str ','];
                                end
                                
                                if rec_with_no_cod2(r)
                                    log_str = [log_str sprintf(' code %s on band %s', track_2,b2)];
                                end
                                this.log.addMessage(this.log.indent(log_str));
                            end
                        end
                        %% get the melbourne wubbena combination for the selected combination
                        mel_wubs = [];
                        i = 1;
                        usable_rec = find(~rec_excluded)';
                        for r = usable_rec
                            fun1 = @(wl1,wl2) 1;
                            fun2 = @(wl1,wl2) -1;
                            [ph_wl] = this.rec_list(r).work.getWideLane(['L' b1 ],['L' b2 ], sys_c); %widelane phase
                            this.wl_comb_codes = [this.wl_comb_codes; [ph_wl.obs_code(1,[1 3 4 6 7])]]; % we have to make sure that the same tracking is also used when forming the phase ionofree combinations afterwards,
                            [pr_nl] = this.rec_list(r).work.getNarrowLane(['C' b1],['C'  b2], sys_c); %narrowlane code
                            [mw] =  this.rec_list(r).work.getTwoFreqComb(ph_wl, pr_nl, fun1, fun2);
                            mel_wubs = [mel_wubs mw];
                            i = i + 1;
                        end
                        
                        [rec_wb, sat_wb, wl_mats, go_ids] = this.estimateWideLaneAndBias(mel_wubs);
                        this.sat_wb(go_ids) = sat_wb;
                        for r = 1 :length(usable_rec)
                            this.wl_mats{r}(:,go_ids) = wl_mats{r};
                        end
                    end
                    % CONSIDER TO DO: save staellite wb to be used by other
                    % recieivers
                end
                
                
                
            end
        end
        
        function pushBackInReceiver(this, s0, res, ls, l_fixed)
            % Save in work the results computed by the network object
            %
            % INPUT
            %   s0          sigma of the solution
            %   res         all the residuals
            %   ls          Least Squares solver object
            %   l_fixed     array of flag for the fixed ambiguities
            %
            % SYNTAX
            %    this = pushBackInReceiver(s0, res, l_fixed)
            
            if nargin < 5
                l_fixed = 0;
            end
            n_rec = length(this.rec_list);
            
            % --- push back the results in the receivers
            for i = 1 : n_rec
                this.rec_list(i).work.xyz = this.coo(i,:);
                idx_res_av = ~isnan(this.clock(:, i));
                [idx_is, idx_pos] = ismembertol(this.common_time.getEpoch(idx_res_av).getGpsTime(), this.rec_list(i).work.time.getGpsTime, 0.002, 'DataScale', 1);
                idx_pos = idx_pos(idx_pos > 0);
                clk = this.clock(idx_res_av, i);
                this.rec_list(i).work.dt(idx_pos) = clk(idx_is) ./ Core_Utils.V_LIGHT;
                if this.state.flag_tropo
                    ztd = this.ztd(idx_res_av, i);
                    this.rec_list(i).work.ztd(idx_pos) = ztd(idx_is);
                    %zhd = this.rec_list(i).work.getAprZhd();
                    this.rec_list(i).work.zwd(idx_pos) = ztd(idx_is) - this.rec_list(i).work.apr_zhd(idx_pos);
                end
                if this.state.flag_tropo_gradient
                    gn = this.ztd_gn(idx_res_av, i);
                    this.rec_list(i).work.tgn(idx_pos) = gn(idx_is);
                    ge = this.ztd_ge(idx_res_av, i);
                    this.rec_list(i).work.tge(idx_pos) = ge(idx_is);
                end
                % sigma of the session
                this.rec_list(i).work.quality_info.s0 = s0;
                this.rec_list(i).work.quality_info.n_epochs = ls.n_epochs(i);
                this.rec_list(i).work.quality_info.n_obs = size(ls.epoch, 1);
                this.rec_list(i).work.quality_info.n_sat = length(unique(ls.sat));
                this.rec_list(i).work.quality_info.n_sat_max = max(hist(unique(ls.epoch * 1000 + ls.sat), max(ls.epoch)));
                this.rec_list(i).work.quality_info.fixing_ratio = (sum(l_fixed(:,1)) / size(l_fixed, 1)) * 100;
                
                % residual
                this.rec_list(i).work.sat.res(:) = 0;
                this.rec_list(i).work.sat.res(idx_pos, :) = res(idx_is, :, i);
                
            end
            
            
            
        end
        
        function pushBackInReceiverNew(this,ls)
            % Save in work the results computed by the network object
            %
            % INPUT
            %   s0          sigma of the solution
            %   res         all the residuals
            %   ls          Least Squares solver object
            %   l_fixed     array of flag for the fixed ambiguities
            %
            % SYNTAX
            %    this = pushBackInReceiver(s0, res, l_fixed)
            
            if nargin < 3
                l_fixed = 0;
            end
            n_rec = length(this.rec_list);
            
            % --- push back the results in the receivers
            for i = 1 : n_rec
                this.rec_list(i).work.xyz = this.coo(i,:);
                this.rec_list(i).work.xyz_vcv = this.coo_vcv(i,:);
                idx_res_av = ~isnan(this.clock(:, i));
                [idx_is, idx_pos] = ismembertol(this.common_time.getEpoch(idx_res_av).getGpsTime(), this.rec_list(i).work.time.getGpsTime, 0.002, 'DataScale', 1);
                idx_pos = idx_pos(idx_pos > 0);
                clk = this.clock(idx_res_av, i);
                this.rec_list(i).work.dt(idx_pos) = clk(idx_is) ./ Core_Utils.V_LIGHT;
                if this.state.flag_tropo
                    ztd = this.ztd(idx_res_av, i);
                    this.rec_list(i).work.ztd(idx_pos) = ztd(idx_is);
                    %zhd = this.rec_list(i).work.getAprZhd();
                    this.rec_list(i).work.zwd(idx_pos) = ztd(idx_is) - this.rec_list(i).work.apr_zhd(idx_pos);
                end
                if this.state.flag_tropo_gradient
                    gn = this.ztd_gn(idx_res_av, i);
                    this.rec_list(i).work.tgn(idx_pos) = gn(idx_is);
                    ge = this.ztd_ge(idx_res_av, i);
                    this.rec_list(i).work.tge(idx_pos) = ge(idx_is);
                end
                s0 = mean(abs(ls.res(ls.phase_obs > 0 & ~ls.outlier_obs)));
                % sigma of the session
                this.rec_list(i).work.quality_info.s0 = s0;
                this.rec_list(i).work.quality_info.n_epochs = length(unique(ls.time_par(ls.rec_par == i)));
                this.rec_list(i).work.quality_info.n_obs = sum(ls.receiver_obs == i);
                this.rec_list(i).work.quality_info.n_sat = length(unique(ls.satellite_obs(ls.receiver_obs == i)));
                this.rec_list(i).work.quality_info.n_sat_max = max(hist(unique(ls.time_obs.getEpoch(ls.receiver_obs == i).getNominalTime().getRefTime(ls.time_obs.minimum.getMatlabTime) * 1000 + double(ls.satellite_obs(ls.receiver_obs == i))), this.rec_list(i).work.quality_info.n_epochs ));
                this.rec_list(i).work.quality_info.fixing_ratio = (sum(l_fixed(:,1)) / size(l_fixed, 1)) * 100; %TBD
                
                % residual
                idx_rec = find( ls.receiver_obs == i);
                %                 % save phase residuals
                idx_ph = find(this.rec_list(i).work.obs_code(:,1) == 'L');
                this.rec_list(i).work.sat.res_ph_by_ph = nan(this.rec_list(i).work.time.length, length(idx_ph));
                for j = 1 : length(idx_ph)
                    ip = idx_ph(j);
                    id_code = Core_Utils.findAinB({[this.rec_list(i).work.system(ip) this.rec_list(i).work.obs_code(ip,:)]}, ls.unique_obs_codes);
                    idx_res = idx_rec(ls.obs_codes_id_obs(idx_rec) == id_code & ls.satellite_obs(idx_rec) == this.rec_list(i).work.go_id(ip));
                    if any(idx_res)
                        [~,idx_time] = ismember(ls.ref_time_obs(idx_res),this.rec_list(i).work.time.getNominalTime.getRefTime(this.rec_list(i).work.time.first.getMatlabTime));
                        this.rec_list(i).work.sat.res_ph_by_ph(idx_time,j) = ls.res(idx_res);
                    end
                end
                % save phase residuals
                idx_pr = find(this.rec_list(i).work.obs_code(:,1) == 'C');
                this.rec_list(i).work.sat.res_pr_by_pr = nan(this.rec_list(i).work.time.length, length(idx_ph));
                for j = 1 : length(idx_pr)
                    ip = idx_pr(j);
                    id_code = Core_Utils.findAinB({[this.rec_list(i).work.system(ip) this.rec_list(i).work.obs_code(ip,:)]}, ls.unique_obs_codes);
                    idx_res = idx_rec(ls.obs_codes_id_obs(idx_rec) == id_code & ls.satellite_obs(idx_rec) == this.rec_list(i).work.go_id(ip));
                    if any(idx_res)
                        [~,idx_time] = ismember(ls.ref_time_obs(idx_res),this.rec_list(i).work.time.getNominalTime.getRefTime(this.rec_list(i).work.time.first.getMatlabTime));
                        this.rec_list(i).work.sat.res_pr_by_pr(idx_time,j) = ls.res(idx_res);
                    end
                end
                % push back electronic bias
                if sum(ls.class_par == LS_Manipulator_new.PAR_REC_EB) > 0
                    idx_eb = find(ls.class_par == LS_Manipulator_new.PAR_REC_EB & ls.rec_par == i);
                    for ii = idx_eb'
                        o_code = ls.unique_obs_codes{ls.obs_codes_id_par(ii)};
                        data = ls.x(ii);
                        this.rec_list(i).work.tracking_bias{ii} = Electronic_Bias(o_code,data);
                    end
                end
                % push back ambiguities
                if sum(ls.class_par == LS_Manipulator_new.PAR_AMB) > 0
                    idx_fix = ls.class_par == ls.PAR_AMB & ls.rec_par == i;
                    idx_fix(idx_fix) = abs(fracFNI(ls.x(idx_fix))) < 1e-9; % fixed ambiguoty
                    idx_fix = find(idx_fix);
                    [ ph,wl,id_ph ] = this.rec_list(i).work.getPhases();
                    amb_mat = Core_Utils.getAmbIdx(this.rec_list(i).work.sat.cycle_slip_ph_by_ph, ph);
                    rec_time_ref = this.rec_list(i).work.time.getRefTime(ls.time_min.getMatlabTime);
                    for amb = idx_fix'
                        % get the index in phases and add them
                        o_code = ls.unique_obs_codes{ls.obs_codes_id_par(amb)};
                        sat = ls.sat_par(amb);
                        time_amb = ls.time_par(amb, :);
                        col_idx = strLineMatch(this.rec_list(i).work.obs_code(id_ph,:),o_code(2:end)) & this.rec_list(i).work.go_id(id_ph) == sat;
                        row_idx = rec_time_ref >= time_amb(1) & rec_time_ref < time_amb(2);
                        a_id = amb_mat(row_idx,col_idx);
                        a_id = a_id(1);
                        ph(amb_mat == a_id) = ph(amb_mat == a_id) - ls.x(amb)*wl(col_idx);
                    end
                    this.rec_list(i).work.setPhases(ph,wl,id_ph );
                end
            end
            if sum(ls.param_class == LS_Manipulator_new.PAR_SAT_EB) > 0
                cs = Core.getCoreSky();
                n_sat = max(ls.sat_par);
                for s = 1 : n_sat
                    idx_eb = find(ls.class_par == LS_Manipulator_new.PAR_SAT_EB & ls.sat_par == s);
                    for ii = idx_eb'
                        if sum(ls.class_par == LS_Manipulator_new.PAR_SAT_EBFR) > 0
                            wl_id = ls.wl_id_par(ii);
                            eb_fr_id = find(ls.class_par == LS_Manipulator_new.PAR_SAT_EBFR & ls.wl_id_par == wl_id & ls.sat_par == s);
                            o_code = ls.unique_obs_codes{ls.obs_codes_id_par(ii)};
                            data = ls.x(eb_fr_id) + ls.x(ii);
                            time_data = ls.getTimePar(eb_fr_id);
                            time_data_min = time_data.minimum;
                            ref_time = time_data.getRefTime(time_data_min.getMatlabTime);
                            % time_apr has to be sampled regualarly
                            time_data_final = 0 : time_data.getRate : (time_data.maximum - time_data_min);
                            [~,iii] = ismembertol(ref_time,time_data_final, 1e-7);
                            iii = Core_Utils.ordinal2logical(iii,length(time_data_final));
                            data_final = zeros(size(time_data_final));
                            data_final(iii) = data;
                            not_fill = find(~iii);
                            iii = find(iii);
                            for nn = not_fill'
                                [~,i_dist ]= min(abs(nn - iii));
                                data_final(nn) = data_final(iii(i_dist));
                            end
                            prm  = ls.ls_parametrization.getParametrization(LS_Manipulator_new.PAR_SAT_EBFR);
                            cs.tracking_bias{s}{ii} = Electronic_Bias(o_code, data_final, time_data_final, prm);
                        else
                            o_code = ls.unique_obs_codes{ls.obs_codes_id_par(ii)};
                            data = ls.x(ii);
                            cs.tracking_bias{s}{ii} = Electronic_Bias(o_code, data);
                        end
                    end
                end
            end
        end
        
        function pushBackAmbiguities(this, x_l1, wl_struct, amb_idx, go_id_ambs,rec_time_indexes)
            % push back in the reciever the reconstructed ambiguites
            n_a_prec = 0;
            for i = 1:length(amb_idx)
                % create a mat containing the l1 and the l2 amniguty
                amb_idx_rec = nan(size(wl_struct.amb_mats{i}));
                amb_idx_rec(rec_time_indexes(:,i),go_id_ambs{i}) = amb_idx{i};
                l1_amb_mat =nan(size(amb_idx_rec));
                n_a_r = max(max(noNaN(amb_idx{i})));
                for a = 1:n_a_r
                    l1_amb_mat(amb_idx_rec == a) = x_l1(a + n_a_prec);
                end
                n_a_prec = n_a_r;
                l2_amb_mat = l1_amb_mat - wl_struct.amb_mats{i};
                % get the measuremnts
                [ph, wl,lid_ph] = this.rec_list(i).work.getPhases();
                id_ph = find(lid_ph);
                cs_slip = this.rec_list(i).work.sat.cycle_slip_ph_by_ph;
                for j = 1 %: size(this.wl_comb_codes,1) %for now onluy single constellatio
                    sys_c = this.wl_comb_codes(j,1);
                    freq_used = this.wl_comb_codes(j,[2 4]);
                    ff= 1;
                    % for each frequency
                    for f = freq_used
                        % get the index of the frquency in the phases
                        lid_f = strLineMatch(this.rec_list(i).work.obs_code(lid_ph,2:3),this.wl_comb_codes(j,1 +(ff-1)*2+(1:2)));
                        id_f = find(lid_f);
                        c_wl = wl(id_f(1));
                        % get the Abx index for the phase measuremetn if
                        % the selected frequency
                        amb_idx_f = Core_Utils.getAmbIdx(cs_slip(:,lid_f),nan2zero(ph(:,lid_f)));
                        for a = unique(noNaN(amb_idx_rec))'
                            sat = find(sum(amb_idx_rec == a)>0); % sta
                            ep_net = find(amb_idx_rec(:,sat) == a); % ep of the network
                            ep = this.rec_time_indexes(ep_net,i); % epoch of the recievrr
                            col_cur_f = this.rec_list(i).work.go_id(id_ph(lid_f)) == sat; % col of the cuurent frequency pahses
                            a_f = noNaN(amb_idx_f(ep,col_cur_f));
                            if length(a_f) > 0
                                a_f = a_f(1);
                                if ff == 1
                                    amb_term = l1_amb_mat(amb_idx_rec == a);
                                    amb_term = amb_term(1);
                                else
                                    amb_term = l2_amb_mat(amb_idx_rec == a);
                                    amb_term = amb_term(1);
                                end
                                ph(amb_idx_f(:,col_cur_f) == a_f,id_f(col_cur_f)) = ph(amb_idx_f(:,col_cur_f) == a_f,id_f(col_cur_f)) - amb_term*c_wl;
                                amb_idx_f(amb_idx_f == a_f) = nan;
                            end
                        end
                        ph_temp = ph(:,lid_f);
                        ph_temp(~isnan(amb_idx_f)) = 0; % if not fixed take off
                        ph(:,lid_f) = ph_temp;
                        this.rec_list(i).work.sat.cycle_slip_ph_by_ph(:,lid_f) = false; % remove cycle slips
                        ff = ff +1;
                        
                    end
                end
                this.rec_list(i).work.setPhases(ph,wl,id_ph);
            end
        end
        
        
        function pushBackSubCooInReceiver(this, time, rate)
            n_rec = length(this.rec_list);
            
            % --- push back the results in the receivers
            for i = 1 : n_rec
                coo = struct();
                coo.coo = Coordinates.fromXYZ(permute(this.coo(i,:,:),[3 2 1]));
                coo.time = time.getCopy();
                coo.rate = rate;
                if isempty( this.rec_list(i).work.add_coo)
                    this.rec_list(i).work.add_coo = coo;
                else
                    this.rec_list(i).work.add_coo(end + 1) = coo;
                end
            end
        end
        
        function alignCodeObservables(this)
            % estimate DCB from a station that contains multiple tracking on the same frequency and apply them to all the reciver of the netwrok
            %
            % SYNTAX:
            %    this.alignCodeObservables()
            n_sat_tot = this.cc.getMaxNumSat();
            sys_cs = this.cc.getAvailableSys;
            n_r = length(this.rec_list);
            
            for sys_c = sys_cs
                n_sat = this.cc.getNumSat(sys_c);
                sys_SS = this.cc.getSys(sys_c);
                bands = sys_SS.CODE_RIN3_2BAND;
                for b = 1 : length(bands)
                    band = bands(b);
                    coderin3attr = sys_SS.CODE_RIN3_ATTRIB{b};
                    % fill a matrix taht contains for each receiver the
                    % channels tracked for all satellites
                    code_aval = false(n_r,n_sat, length(coderin3attr));
                    best_code_rec_idx = nan(n_r,1); % index of the best tracking availiable for the receievr
                    for r = 1 : n_r
                        % check wether both frequency are available
                        rec = this.rec_list(r).work;
                        has_fr = sum(strLineMatch(rec.obs_code(rec.system == sys_c,1:2),['C' band])) > 0;
                        has_frs(r) = has_fr;
                        if has_fr
                            % fill the availability matrix for the reciever
                            for s = 1:n_sat
                                idx_s = rec.obs_code(:,1) == 'C' & rec.prn == s & rec.system' == sys_c;
                                if sum(idx_s) > 0
                                    idx_sb = idx_s & rec.obs_code(:,2) == band;
                                    track_code = rec.obs_code(idx_sb,3);
                                    
                                    for t = 1 : length(track_code)
                                        code_aval(r,s,coderin3attr == track_code(t)) = true;
                                    end
                                    best_code_rec_idx(r) = min(best_code_rec_idx(r),find(code_aval(r,s,:),1,'first'));
                                end
                            end
                        end
                    end
                    if sum(has_fr) >0
                        % Find the reference receiver
                        cod_aval_no_sat = sum(permute(code_aval,[1 3 2]),3) > 0;
                        n_code_rec = sum(cod_aval_no_sat,2);
                        [ref_code_id] = min(best_code_rec_idx);
                        n_code_rec(best_code_rec_idx ~= ref_code_id) = 0;
                        [~,b_rec_idx] = max(n_code_rec); % the reciever with the best code that tracks the more channel
                        present_code = sum(sum(permute(code_aval,[3 2 1]),3) > 0,2) >0; % code available at least at one receiver
                        aligned_code  = false(size(present_code));
                        aligned_code(ref_code_id) = true; %code aligned
                        not_alignable_codes = false;
                        this.log.addMessage(sprintf('Aligning all %s psudorange observations to tarcking %s',[sys_c band],coderin3attr(ref_code_id)));
                        
                        while sum(~aligned_code(present_code)) > 0 && ~not_alignable_codes
                            t_best_rec = find(cod_aval_no_sat(b_rec_idx,:));
                            trck_ref = coderin3attr(t_best_rec(1));
                            for t = 2 : length(t_best_rec);
                                trck = coderin3attr(t_best_rec(t));
                                if ~aligned_code(t_best_rec(t)) % if still not aligned
                                    for ss = 1 : n_sat
                                        [obs_ref] = this.rec_list(b_rec_idx).work.getObs(['C' band trck_ref], sys_c, ss);
                                        [obs_cur] = this.rec_list(b_rec_idx).work.getObs(['C' band trck], sys_c, ss);
                                        if numel(obs_ref) > 0 & numel(obs_cur) >0 % sat missing
                                            dcb = strongMean((zero2nan(obs_cur) - zero2nan(obs_ref))');
                                            
                                            for rr = 1 : n_r
                                                [obs2align, idx_obs] = this.rec_list(rr).work.getObs(['C' band trck], sys_c, ss);
                                                obs2align = obs2align - dcb;
                                                this.rec_list(rr).work.setObs(obs2align, idx_obs)
                                                this.rec_list(rr).work.aligned(idx_obs) = true;
                                            end
                                        end
                                    end
                                    aligned_code(t_best_rec(t)) = true;
                                end
                            end
                            if sum(~aligned_code(present_code)) > 0
                                not_alignable_codes = true;
                                for tt = find(~aligned_code&present_code)
                                    % find if a not aligned code is linkable to
                                    % the other tracking trough an other
                                    % receievr
                                    for tt2 = find(aligned_code)'
                                        if sum(cod_aval_no_sat(:,tt).*cod_aval_no_sat(:,tt2)) > 0 % if at least a receiver has both codes
                                            b_rec_idx = find(cod_aval_no_sat(:,tt).*cod_aval_no_sat(:,tt2),1,'first');
                                            ref_code_id = tt2;
                                            not_alignable_codes = false;
                                            break
                                        end
                                    end
                                    if ~not_alignable_codes
                                        break
                                    end
                                end
                                if not_alignable_codes
                                    this.log.addMessage(sprintf('Tracking %s on banf %s are not linkable to others tracking so can not be aligned',coderin3attr(~aligned_code& present_code),[sys_c band]));
                                end
                            end
                        end
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                end
            end
        end
        
        function exportCrd(this, file_prefix)
            % export the current value of the coordinate to a bernese CRD file, if multiple session are available for the stations save only the last coordinates
            %
            % SYNTAX:
            % this.exportCrd(this)
            if nargin < 2 || isempty(file_prefix)
                %[~,file_prefix] =fileparts( rec.state.getHomeDir);
                file_prefix = [this.state.getOutPrefix '_'];
            end
            if ndims(this.coo) < 3
                st_time  = this.common_time.first;
                en_time = this.common_time.last;
                coo = this.coo;
            else
                st_time = this.state.getSessionLimits.first;
                st_time.addSeconds(this.coo_rate * (size(this.coo,3) - 1));
                en_time = this.state.getSessionLimits.first;
                en_time.addSeconds(this.coo_rate * size(this.coo,3));
                coo = this.coo(:,:,end);
            end
            [~,~,sod_s] = st_time.getDOY();
            [~,~,sod_f] = en_time.getDOY();
            if sum(sod_f == '0') == 5
                sod_f = '86400';
            end
            [~,doy] = st_time.getDOY();
            fpath  = sprintf('%s/%s%02d%03d.%05d-%05d.CRD',this.state.getOutDir, file_prefix, st_time.getYY, doy, sod_s,sod_f);
            fid = fopen(fpath,'w');
            now_time = GPS_Time.now();
            fprintf(fid, ['                                                                 ' upper(now_time.toString('dd-mmm-yy HH:MM')) ' \n']);
            
            fprintf(fid, ['--------------------------------------------------------------------------------\n']);
            fprintf(fid, ['LOCAL GEODETIC DATUM: WGS - 84          EPOCH: ' st_time.toString('dd-mm-yy HH:MM:SS') '\n\n']);
            
            fprintf(fid,'NUM  STATION NAME           X (M)          Y (M)          Z (M)     FLAG\n\n');
            n_rec = length(this.rec_list);
            for i = 1 : n_rec
                fprintf(fid,sprintf('%3d  %s              %13.5f  %13.5f  %13.5f    %s\n', i, upper(this.rec_list(i).getMarkerName4Ch), coo(i,:), iif(sum(this.id_ref == i) > 0, 'F', 'P')));
            end
            fclose(fid);
        end
        
    end
    methods (Static)
        function [rec_wb, sat_wb, wl_amb_mat, go_ids] = estimateWideLaneAndBias(mel_wub_mat)
            % estiamet the wodelanes and their biases
            %
            % SYNTAX:
            %    [rec_wb, sat_wb, wl_amb_mat] = Network.estimateWideLaneAndBias(mel_wub_cell)
            
            % get the common sat idx
            go_ids = [];
            n_rec = length(mel_wub_mat);
            for r = 1 : n_rec
                go_ids = [go_ids mel_wub_mat(r).go_id];
            end
            [occurreces, go_ids ]=hist(go_ids,unique(go_ids));
            go_ids(occurreces <2) = [];
            % estimate the widelane satelite bias form the first available
            % receiver
            n_sat =  length(go_ids);
            sat_wb = nan(1,n_sat);
            sat_wb_r_id = zeros(1, n_sat); % using which receiver the bias has been calculated
            for s = 1 : length(go_ids)
                go_id = go_ids(s);
                found = false;
                r = 1;
                while ~found && r <= n_rec
                    idx_s = mel_wub_mat(r).go_id == go_id;
                    if sum(idx_s) >0
                        found = true;
                        sat_wb(s) = Core_Utils.estimateFracBias(zero2nan(mel_wub_mat(r).getObsCy(find(idx_s))), mel_wub_mat(r).cycle_slip(:,idx_s));
                        sat_wb_r_id(s) = r;
                    end
                    r = r +1;
                end
            end
            rec_wb = zeros(1, n_rec);
            for r = 1 : n_rec
                % eliminate the satellites that has been used to compute
                % the satellite bias
                to_eliminate_sat = sat_wb_r_id == r;
                [~,rec_gi_idx, idx_satwb] = intersect(mel_wub_mat(r).go_id, go_ids(~to_eliminate_sat)); %idx of the go id in the obervation set
                if ~isempty(rec_gi_idx)
                    sat_wb_r = sat_wb(~to_eliminate_sat);
                    sat_wb_r = sat_wb_r(idx_satwb);
                    rec_wb(r) = Core_Utils.estimateFracBias(zero2nan(mel_wub_mat(r).getObsCy(rec_gi_idx)) - repmat(sat_wb_r,mel_wub_mat(r).time.length,1), mel_wub_mat(r).cycle_slip(:,rec_gi_idx));
                    % update the satellite bias knpwing the receiver one
                    if sum(to_eliminate_sat) > 0
                        for s = find(to_eliminate_sat)
                            go_id = go_ids(s);
                            idx_s = mel_wub_mat(r).go_id == go_id;
                            sat_wb(s) = sat_wb(s) - rec_wb(r);
                        end
                    end
                end
            end
            wl_amb_mat = {};
            % final estamiation of the widelane
            for r = 1 : n_rec
                cy = mel_wub_mat(r).getObsCy;
                [~,rec_gi_idx, goids_idx2] = intersect(go_ids, mel_wub_mat(r).go_id ); %idx of the go id in the obervation set
                cy = zero2nan(cy(:,goids_idx2)) - rec_wb(r) - repmat(sat_wb(rec_gi_idx),mel_wub_mat(r).time.length,1);
                wl_amb_mat{r} = nan(mel_wub_mat(r).time.length,length(go_ids));
                [~,wl_amb_mat{r}(:,rec_gi_idx)]  = Core_Utils.estimateFracBias(cy, mel_wub_mat(r).cycle_slip(:,goids_idx2));
            end
            % OPTIONAL ->refine the bias estimation with a LS adjustemtn
        end
        
    end
end
