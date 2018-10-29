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
%    |___/                    v 0.6.0 alpha 4 - nightly
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
classdef Network < handle
    properties
        rec_list
        state
        cc
        net_id           % id of the receiver in core, this uniquely identify the network
        
        common_time      % gps_time
        rec_time_indexes % indexes
        coo              % [n_coo x n_rec x n_sol] receiver coordinates
        coo_rate         % rate of the coordinate solution in seconds
        clock            % [n_epoch x n_rec] reciever clock
        ztd              % [n_epoch x n_rec] reciever ZTD
        ztd_gn           % [n_epoch x n_rec] reciever ZTD gradients north
        ztd_ge           % [n_epoch x n_rec] reciever ZTD gradients east
        amb              % {n_rec} recievers ambiguity
        log
        pos_indexs_tc       % index for subpositions
        id_ref
        
        apriori_info     % field to keep apriori info [ambiguity, tropo, ...] to be used in the adjustment
    end
    methods
        function this = Network(rec_list, net_id)
            if nargin < 2
                net_id = [];
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
            this.coo_rate = [];
            this.clock = [];
            this.ztd = [];
            this.ztd_gn = [];
            this.ztd_ge = [];
            this.amb = [];
            this.pos_indexs_tc = [];
            this.id_ref = [];
        end
        
        function adjust(this, id_ref, coo_rate)
            % Adjust the GNSS network
            %
            % INPUT
            %     id_ref : [1,n_rec] boolean, receivers to be choosen as reference, their value mean will be set to zero
            %
            % SYNATAX
            %    this. adjustNetwork(id_ref)
            if nargin < 3
                coo_rate = [];
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

            is_empty_recs = this.rec_list.isEmptyWork_mr;
            if sum(~is_empty_recs) > 1
                e = find(is_empty_recs);
                if ~isempty(e)
                    this.rec_list(e) = [];
                    lid_ref(e) = [];
                end
                id_ref = find(lid_ref);
                this.id_ref = id_ref;
                
                if this.state.getReweight() == 1
                    n_clean = 0;
                else
                    this.log.addMessage(this.log.indent('Network reweight perform only a simple outlier detection on the residuals'), 2);
                    n_clean = 3;
                end

                while n_clean >= 0
                    ls = LS_Manipulator(this.rec_list(1).cc);
                    
                    if this.state.flag_amb_pass && this.state.getCurSession > 1 && ~isempty(this.apriori_info)
                        f_time = GPS_Time.now();
                        f_time.addSeconds(1e9);
                        rate = 0;
                        for i = 1 : length(this.rec_list)
                            [~,limc] = this.state.getSessionLimits(this.state.getCurSession);
                            idx_fv = find(this.rec_list(i).work.time >= limc.first,1,'first');
                            f_time = min(this.rec_list(i).work.time.getEpoch(idx_fv), f_time);
                            rate = max(rate,this.rec_list(i).work.time.getRate);
                        end
                        t_dist = f_time - this.apriori_info.epoch;
                        if t_dist > -0.02 && t_dist <= (rate + 0.02)
                            ls.apriori_info = this.apriori_info;
                        end
                    end
                    
                    [this.common_time, this.rec_time_indexes]  = ls.setUpNetworkAdj(this.rec_list, coo_rate);
                    if isempty(this.rec_time_indexes)
                        return
                    else
                        n_time = this.common_time.length;
                        n_rec = length(this.rec_list);
                        if this.state.flag_tropo
                            ls.setTimeRegularization(ls.PAR_TROPO, (this.state.std_tropo)^2 / 3600 * ls.rate );
                        end
                        if this.state.flag_tropo_gradient
                            ls.setTimeRegularization(ls.PAR_TROPO_N, (this.state.std_tropo_gradient)^2 / 3600 * ls.rate );
                            ls.setTimeRegularization(ls.PAR_TROPO_E, (this.state.std_tropo_gradient)^2 / 3600 * ls.rate );
                        end
                        [x, res, s0, Cxx, l_fixed] = ls.solve;
                        %[x, res] = ls.solve;
                        %res = res(any(res(:,:,2)'), :, :);
                        
                        % cleaning -----------------------------------------------------------------------------------------------------------------------------
                        %s0 = mean(abs(res(res~=0)));
                        this.log.addMessage(this.log.indent(sprintf('Network solution computed,  s0 = %.4f', s0)));
                        %figure; plot(res(:,:,2)); ylim([-0.05 0.05]); dockAllFigures;
                        if s0 < 0.0025 && ... % low sigma0 of the computation
                           all(abs(res(:)) < this.state.pp_max_phase_err_thr) && ...  % no residuals above thr level
                           all(std(zero2nan(reshape(permute(res(:,:,:),[1 3 2]), size(res,1) * size(res,3), size(res,2))),1,2,'omitnan') < 9e-3) % low dispersion of the residuals
                                
                            % I can be satisfied
                            n_clean = -1;
                        end
                         if n_clean > 0
                            out_found = 0;
                            for r = 1 : length(this.rec_list)
                                tmp_work = this.rec_list(r).work;
                                
                                % Index of receiver obs present in res
                                id_sync_res = round((round(tmp_work.time.getRefTime(this.common_time.first.getMatlabTime) * 1e5) * 1e-5 /  this.common_time.getRate) + 1);
                                id_sync_rec = (1 : numel(id_sync_res))';
                                id_ok = id_sync_res > 0 & id_sync_res <= size(res,1);
                                id_sync_res = id_sync_res(id_ok);
                                id_sync_rec = id_sync_rec(id_ok);
                                
                                go_id = tmp_work.go_id(this.rec_list(r).work.findObservableByFlag('L'));
                                
                                res_rec = zeros(size(tmp_work.sat.outlier_idx_ph));
                                res_rec(id_sync_rec, :) = res(id_sync_res, go_id, r);
                                                                
                                if s0 > 1
                                    % The solution is really bad, something bad happened
                                    % Extreme experiment to recover a valid positioning
                                    id_ko = abs(mean(zero2nan(res_rec), 'omitnan')) > 1;
                                    id_ko = sparse(repmat(id_ko, size(res_rec, 1), 1));
                                else
                                    id_ko = sparse(tmp_work.search4outliers(res_rec, n_clean+1));
                                    % At the first loop start filtering bad satellites
                                    if n_clean > 1
                                        bad_sat = std(res_rec, 'omitnan') > max(0.01, perc(noNaN(std(res_rec, 'omitnan')), 0.9)) | max(abs(res_rec)) > 0.15;
                                        if any(bad_sat) && any(any(id_ko(:, bad_sat)))
                                            id_ko(:, ~bad_sat) = false;
                                        end
                                    end                                    
                                end
                                if any(id_ko(:))
                                    tmp_work.addOutliers(id_ko, true);
                                    out_found = out_found + 1;
                                end    
                            end
                            if out_found == 0
                                % no need to iterate
                                n_clean = -1;
                            end
                        end
                        n_clean = n_clean - 1;
                        % end of cleaning ----------------------------------------------------------------------------------------------------------------------

                    end
                end
                
                if s0 < 0.02
                    % intilaize array for results
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
                    [sss_lim] = this.state.getSessionLimits(this.state.getCurSession());
                    st_time = sss_lim.first;
                    for i = 1 : 3
                        if this.state.coo_rates(i) ~= 0
                            this.coo_rate = this.state.coo_rates(i);
                            ls.pos_indexs_tc = {};
                            for j = 2 : n_rec
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
            else
                this.log.addWarning('Not enough receivers (< 2), skipping network solution');
            end
            
        end
        
        function initOut(this,ls)
            n_time = this.common_time.length;
            n_rec = length(this.rec_list);
            n_set_coo = length(ls.getCommonPosIdx);
            this.pos_indexs_tc = ls.pos_indexs_tc;
            this.clock = zeros(n_time, n_rec);
            this.coo = nan(n_rec, 3, n_set_coo);
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
                if i > 1 % coordiantes are always zero on first receiver
                    coo = [x(x(:,2) == 1 & idx_rec,1) x(x(:,2) == 2 & idx_rec,1) x(x(:,2) == 3 & idx_rec,1)];
                    if ~isempty(this.pos_indexs_tc)
                        this.coo(i,:,this.pos_indexs_tc{i-1}) = nan2zero(this.coo(i,:,this.pos_indexs_tc{i-1})) + permute(coo, [3 2 1]);
                    else
                        this.coo(i,:) = nan2zero(this.coo(i,:)) + coo;
                    end
                else
                    this.coo(i,:) = nan2zero(this.coo(i,:));
                end
                clk = x(x(:,2) == LS_Manipulator.PAR_REC_CLK & idx_rec,1);
                this.clock(~isnan(this.rec_time_indexes(:,i)),i) = nan2zero(this.clock(~isnan(this.rec_time_indexes(:,i)),i)) + clk;
                
                if this.state.flag_tropo
                    ztd = x(x(:,2) == LS_Manipulator.PAR_TROPO & idx_rec,1);
                    this.ztd(~isnan(this.rec_time_indexes(:,i)),i) = nan2zero(this.ztd(~isnan(this.rec_time_indexes(:,i)),i))  + ztd;
                end
                
                if this.state.flag_tropo_gradient
                    gn = x(x(:,2) == LS_Manipulator.PAR_REC_CLK & idx_rec,1);
                    this.ztd_gn(~isnan(this.rec_time_indexes(:,i)),i) = nan2zero(this.ztd_gn(~isnan(this.rec_time_indexes(:,i)),i)) + gn;
                    
                    ge = x(x(:,2) == LS_Manipulator.PAR_REC_CLK & idx_rec,1);
                    this.ztd_ge(~isnan(this.rec_time_indexes(:,i)),i) = nan2zero(this.ztd_ge(~isnan(this.rec_time_indexes(:,i)),i)) + ge;
                end
            end
        end
        
        function changeReferenceFrame(this, id_ref)
            n_rec = length(this.rec_list); 
            n_time = this.common_time.length;
            % ALL OF THIS MAKES NO SENSE TO ME (Andrea). Now it should ;) (Giulio)
            %--- transform the result in the desired free network
            
            if ~isnan(id_ref(1)) && ~(length(id_ref) == 1 && id_ref == 1)
                
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
            for sys_c = this.cc.getAvailableSys
                sys_idx = this.cc.SYS_C == sys_c;
                WM = wide_laneM{sys_idx}; %% get the right widelane matrix
                n_sat = this.cc.getNumSat(sys_c);
                for w = 1 : size(WM,2)
                    sys_SS = this.cc.getSys(sys_c);
                    b1 = sys_SS.CODE_RIN3_2BAND(WM(w,:) == 1);
                    b2 = sys_SS.CODE_RIN3_2BAND(WM(w,:) == -1);
                    has_frs = false(n_r,1);
                    coderin3attr1 = sys_SS.CODE_RIN3_ATTRIB{WM(w,:) == 1};
                    coderin3attr2 = sys_SS.CODE_RIN3_ATTRIB{WM(w,:) == -1};
                    b1code_aval = false(n_r,n_sat, length(coderin3attr1));
                    b2code_aval = false(n_r,n_sat, length(coderin3attr2));
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
                                idx_sb1 = idx_s & rec.obs_code(:,2) == b1;
                                idx_sb2 = idx_s & rec.obs_code(:,2) == b2;
                                track_code_b1 = rec.obs_code(idx_sb1,3);
                                track_code_b2 = rec.obs_code(idx_sb2,3);
                                for t = 1 : length(track_code_b1)
                                    b1code_aval(r,s,coderin3attr1 == track_code_b1(t)) = true;
                                end
                                for t = 1 : length(track_code_b2)
                                    b2code_aval(r,s,coderin3attr2 == track_code_b2(t)) = true;
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
                    % remove code that are not observed by more than one
                    % receiver
                    nt2cb1 = sum(sum(squeeze(sum(b1code_aval,1)>1)),1)> 0; 
                    nt2cb2 = sum(sum(squeeze(sum(b2code_aval,1)>1)),1)> 0; 
                    b1code_aval(:,:,~nt2cb1) = [];
                    b2code_aval(:,:,~nt2cb2) = [];
                    % Now set up the least squares system
                    
                    
                    % now we have to align all the codes so when we average
                    % them we do not have different bias
                    
                      
                        % if the best code for different receiver is not
                        % the same check wether the differential code
                        % bias can be estimated from the other receiver
                        
                        % if they can be estimated apply the bias 
                        
                        % solve the widelane cycle slip
                        
                        % the esttimate the fractional bias 
                        
                        % we are all very happy
                        end
                        
                    for r = 1 : length(this.rec_list)
                        % get the phase widelane
                        this.rec_list(r).getLin
                    end
                end
        end
        
      
            
  
       
        
        function pushBackInReceiver(this, s0, res)
            n_rec = length(this.rec_list);
            
            % --- push back the results in the receivers
            for i = 1 : n_rec
                this.rec_list(i).work.xyz = this.coo(i,:);
                idx_res_av = ~isnan(this.clock(:, i));
                [idx_is, idx_pos] = ismembertol(this.common_time.getEpoch(idx_res_av).getGpsTime(), this.rec_list(i).work.time.getGpsTime, 0.002, 'DataScale', 1);
                idx_pos = idx_pos(idx_pos > 0);
                clk = this.clock(idx_res_av, i);
                this.rec_list(i).work.dt(idx_pos) = clk(idx_is);
                if this.state.flag_tropo
                    ztd = this.ztd(idx_res_av, i);
                    this.rec_list(i).work.ztd(idx_pos) = ztd(idx_is);
                    zhd = this.rec_list(i).work.getAprZhd();
                    this.rec_list(i).work.zwd(idx_pos) = ztd(idx_is) - zhd(idx_pos);
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
                st_time = this.state.getSessionLimits(this.state.getCurSession).first;
                st_time.addSeconds(this.coo_rate * (size(this.coo,3) - 1));
                en_time = this.state.getSessionLimits(this.state.getCurSession).first;
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
end
