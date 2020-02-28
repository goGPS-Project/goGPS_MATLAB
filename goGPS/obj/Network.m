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
%    |___/                    v 1.0b7
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
        
        apriori_info      % field to keep apriori info [ambiguity, tropo, ...] to be used in the adjustment
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
            works = [this.rec_list.work];
            is_empty_recs = ~works.hasRangeObs_mr;
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
                l_fixed = 0; % nothing is Fise
                is_empty_recs = ~this.rec_list.hasRangeObs_mr;
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
                phase = true;
                
                parametrization = LS_Parametrization();
                
                if phase
                    param_selection =  [ls.PAR_AMB];
                else
                    param_selection =  [];
                end
                state = this.state;
                if this.state.flag_coo_net
                    param_selection = [param_selection;
                        LS_Manipulator_new.PAR_REC_X;
                        LS_Manipulator_new.PAR_REC_Y;
                        LS_Manipulator_new.PAR_REC_Z;
                        ];
                    
                    % time parametrization coordinates
                    parametrization.rec_x(1) = state.tparam_coo_net;
                    parametrization.rec_y(1) = state.tparam_coo_net;
                    parametrization.rec_z(1) = state.tparam_coo_net;
                    
                    
                    parametrization.rec_x(4) = state.fparam_coo_net;
                    parametrization.rec_y(4) = state.fparam_coo_net;
                    parametrization.rec_z(4) = state.fparam_coo_net;
                end
                
                param_selection = [param_selection;
                    ls.PAR_IONO;];
                
                if this.state.flag_ztd_net
                    param_selection = [param_selection;
                        ls.PAR_TROPO;];
                end
                
                if this.state.flag_grad_net
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
                if state.flag_rec_clock_net
                    if state.flag_phpr_rec_clock_net
                        param_selection =  [param_selection;
                            LS_Manipulator_new.PAR_REC_CLK_PR;
                            LS_Manipulator_new.PAR_REC_CLK_PH;
                            ];
                    else
                        param_selection =  [param_selection;
                            LS_Manipulator_new.PAR_REC_CLK;
                            LS_Manipulator_new.PAR_REC_PPB;
                            ];
                    end
                end
                if state.flag_sat_clock_net
                    if state.flag_phpr_sat_clock_net
                        param_selection =  [param_selection;
                            LS_Manipulator_new.PAR_SAT_CLK_PR;
                            LS_Manipulator_new.PAR_SAT_CLK_PH;
                            ];
                    else
                        param_selection =  [param_selection;
                            LS_Manipulator_new.PAR_SAT_CLK;
                            LS_Manipulator_new.PAR_SAT_PPB;
                            ];
                    end
                end
                
                if state.flag_rec_trkbias_net
                    param_selection =  [param_selection;
                        LS_Manipulator_new.PAR_REC_EB;];
                    parametrization.setTimeParametrization(LS_Manipulator_new.PAR_REC_EB, state.tparam_rec_trkbias_net );
                    if state.tparam_rec_trkbias_net > 1 && state.rate_rec_trkbias_net > 0
                        parametrization.setRate(LS_Manipulator_new.PAR_REC_EB, state.rate_rec_trkbias_net );
                    end
                    
                end
                
                if state.flag_rec_ifbias_net
                    param_selection =  [param_selection;
                        LS_Manipulator_new.PAR_REC_EBFR;];
                    parametrization.setTimeParametrization(LS_Manipulator_new.PAR_REC_EBFR, state.tparam_rec_ifbias_net );
                    if state.tparam_rec_ifbias_net > 1 && state.rate_rec_ifbias_net > 0
                        parametrization.setRate(LS_Manipulator_new.PAR_REC_EBFR, state.rate_rec_ifbias_net );
                    end
                end
                
                if state.flag_sat_trkbias_net
                    param_selection =  [param_selection;
                        LS_Manipulator_new.PAR_SAT_EB;];
                    parametrization.setTimeParametrization(LS_Manipulator_new.PAR_SAT_EB, state.tparam_sat_trkbias_net );
                    if state.tparam_sat_trkbias_net > 1 && state.rate_sat_trkbias_net > 0
                        parametrization.setRate(LS_Manipulator_new.PAR_SAT_EB, state.rate_sat_trkbias_net );
                    end
                    
                end
                
                if state.flag_sat_ifbias_net
                    param_selection =  [param_selection;
                        LS_Manipulator_new.PAR_SAT_EBFR;];
                    parametrization.setTimeParametrization(LS_Manipulator_new.PAR_SAT_EBFR, state.tparam_sat_ifbias_net );
                    if state.tparam_sat_ifbias_net > 1 && state.rate_sat_ifbias_net > 0
                        parametrization.setRate(LS_Manipulator_new.PAR_SAT_EBFR, state.rate_sat_ifbias_net );
                    end
                end
                
                
              
                
                if ~this.state.flag_iono_net
                    parametrization.iono(2) = LS_Parametrization.ALL_REC;
                end
                 if state.tparam_ztd_net == 1
                    parametrization.tropo(1) = parametrization.EP_WISE;
                elseif state.tparam_ztd_net == 2
                    parametrization.tropo(1) = parametrization.SPLINE_LIN;
                    parametrization.tropo_opt.spline_rate = state.rate_ztd_net;
                elseif state.tparam_ztd_net == 3
                    parametrization.tropo(1) = parametrization.SPLINE_CUB;
                    parametrization.tropo_opt.spline_rate = state.rate_ztd_net;
                end
                                
                % Use spline for estimating ZTD gradients
                if state.tparam_grad_net == 1
                    parametrization.tropo_n(1) = parametrization.EP_WISE;
                    parametrization.tropo_e(1) = parametrization.EP_WISE;
                elseif state.tparam_grad_net == 2
                    parametrization.tropo_n(1) = parametrization.SPLINE_LIN;
                    parametrization.tropo_n_opt.spline_rate = state.rate_grad_net;
                    
                    parametrization.tropo_e(1) = parametrization.SPLINE_LIN;
                    parametrization.tropo_e_opt.spline_rate = state.rate_grad_net;
                elseif state.tparam_grad_net == 3
                    parametrization.tropo_n(1) = parametrization.SPLINE_CUB;
                    parametrization.tropo_n_opt.spline_rate = state.rate_grad_net;
                    
                    parametrization.tropo_e(1) = parametrization.SPLINE_CUB;
                    parametrization.tropo_e_opt.spline_rate = state.rate_grad_net;
                end
                ls.setUpNET(this.rec_list, coo_rate, '???', param_selection, parametrization);
               
                if this.state.flag_free_net_tropo
                    ls.free_tropo = true;
                end
                
                this.is_tropo_decorrel = this.state.isReferenceTropoEnabled;
                this.is_coo_decorrel = free_net;
                
                
                 if state.areg_ztd_net > 0
                   ls.absValRegularization(ls.PAR_TROPO, (state.areg_ztd_net)^2);
                end
                if state.areg_grad_net > 0
                    ls.absValRegularization(ls.PAR_TROPO_N, (state.areg_grad_net)^2);
                    ls.absValRegularization(ls.PAR_TROPO_E, (state.areg_grad_net)^2);
                end
                if state.areg_rec_clock_net > 0
                    if  state.flag_phpr_rec_clock_net
                        ls.absValRegularization(ls.PAR_REC_CLK_PH, (state.areg_rec_clock_net)^2);
                        ls.absValRegularization(ls.PAR_REC_CLK_PR, (state.areg_rec_clock_net)^2);
                    else
                        ls.absValRegularization(ls.PAR_REC_CLK, (state.areg_rec_clock_net)^2);
                    end
                end
                if state.areg_rec_ifbias_net > 0
                    ls.absValRegularization(ls.PAR_REC_EBFR, (state.areg_rec_ifbias_net)^2);
                end
                if state.areg_rec_trkbias_net > 0
                    ls.absValRegularization(ls.PAR_REC_EB, (state.areg_rec_trkbias_net)^2);
                end
                
                if state.dreg_ztd_net > 0
                    ls.timeRegularization(ls.PAR_TROPO, (state.dreg_ztd_net)^2/ 3600);
                end
                if state.dreg_grad_net > 0
                    ls.timeRegularization(ls.PAR_TROPO_N, (state.dreg_grad_net)^2/ 3600);
                    ls.timeRegularization(ls.PAR_TROPO_E, (state.dreg_grad_net)^2/ 3600);
                end
                
                if state.dreg_rec_ifbias_net > 0
                    ls.timeRegularization(ls.PAR_REC_EBFR, (state.dreg_rec_ifbias_net)^2/ 3600);
                end
                if state.dreg_rec_trkbias_net > 0
                    ls.timeRegularization(ls.PAR_REC_EB, (state.dreg_rec_trkbias_net)^2/ 3600);
                end
                
                if state.areg_sat_clock_net > 0
                    if  state.flag_phpr_sat_clock_net
                        ls.absValRegularization(ls.PAR_SAT_CLK_PH, (state.areg_sat_clock_net)^2);
                        ls.absValRegularization(ls.PAR_SAT_CLK_PR, (state.areg_sat_clock_net)^2);
                    else
                        ls.absValRegularization(ls.PAR_SAT_CLK, (state.areg_sat_clock_net)^2);
                    end
                end
                if state.areg_sat_ifbias_net > 0
                    ls.absValRegularization(ls.PAR_SAT_EBFR, (state.areg_sat_ifbias_net)^2);
                end
                if state.areg_sat_trkbias_net > 0
                    ls.absValRegularization(ls.PAR_SAT_EB, (state.areg_sat_trkbias_net)^2);
                end
                
                if state.dreg_sat_ifbias_net > 0
                    ls.timeRegularization(ls.PAR_SAT_EBFR, (state.dreg_sat_ifbias_net)^2/ 3600);
                end
                if state.dreg_sat_trkbias_net > 0
                    ls.timeRegularization(ls.PAR_SAT_EB, (state.dreg_sat_trkbias_net)^2/ 3600);
                end
                
                
                if Core.isGReD
                    % distance regularization to be set up
                    
                end
                this.common_time = ls.unique_time;
                ls.solve(Core.getState.net_amb_fix_approach >1);
%                 idx_fix = ls.class_par == ls.PAR_AMB;
%                 idx_fix(idx_fix) = abs(fracFNI(ls.x(idx_fix))) < 1e-9; % fixed ambiguoty
%                 ls.removeEstParam(idx_fix);
                ls.snoopGatt(Core.getState.getMaxPhaseErrThr, Core.getState.getMaxCodeErrThr);
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
            n_set_coo = size(unique(ls.time_par(ls.class_par == ls.PAR_REC_X & ls.rec_par == 1,1)),1);
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
                        if i ~= this.id_ref && false
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
                    
                    this.coo(i,1:size(coo,2)) = nan2zero(this.coo(i,1:size(coo,2))) + coo;
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
                state = this.state;
                if this.state.flag_ztd_net
                    if state.tparam_ztd_net > 1
                        if state.tparam_ztd_net == 2
                            spline_order = 1;
                        elseif state.tparam_ztd_net == 3
                            spline_order = 3;
                        end
                        idx_trp = ls.class_par == LS_Manipulator_new.PAR_TROPO & idx_rec;
                        if sum(idx_trp) > 0
                            tropo = ls.x(idx_trp);
                            tropo_dt = rem(this.common_time.getNominalTime(ls.obs_rate) - ls.getTimePar(idx_trp).minimum, this.state.rate_ztd_net)/ this.state.rate_ztd_net;
                            tropo_idx = floor((this.common_time.getNominalTime(ls.obs_rate) - ls.getTimePar(idx_trp).minimum)/this.state.rate_ztd_net);
                            [~,tropo_idx] = ismember(tropo_idx*this.state.rate_ztd_net, ls.getTimePar(idx_trp).getNominalTime(ls.obs_rate).getRefTime(ls.getTimePar(idx_trp).minimum.getMatlabTime));
                            valid_ep = tropo_idx ~=0 & tropo_idx <= (length(tropo)-3);
                            spline_base = Core_Utils.spline(tropo_dt(valid_ep),spline_order);
                            
                            ztd =sum(spline_base .* tropo(repmat(tropo_idx(valid_ep), 1, spline_order + 1) + repmat((0 : spline_order), numel(tropo_idx(valid_ep)), 1)), 2);
                            
                            this.ztd(valid_ep,i) = nan2zero(this.ztd(valid_ep,i))  + ztd;
                        end
                    else
                        idx_trp = ls.class_par == LS_Manipulator_new.PAR_TROPO & idx_rec;
                        tropo = ls.x(idx_trp);
                        time_tropo = ls.time_par(idx_trp);
                        [~,idx_time_tropo] = ismember(time_tropo, this.common_time.getNominalTime(ls.obs_rate).getRefTime(ls.time_min.getMatlabTime));
                        this.ztd(idx_time_tropo,i) = nan2zero(this.clock(idx_time_tropo,i)) + tropo;
                    end
                end
                
                if this.state.flag_grad_net
                    if state.tparam_grad_net > 1
                           if state.tparam_grad_net == 2
                            spline_order = 1;
                            elseif state.tparam_grad_net == 3
                                spline_order = 3;
                            end 
                        idx_trp_n = ls.class_par == LS_Manipulator_new.PAR_TROPO_E & idx_rec;
                        if sum(idx_trp_n) > 0
                        tropo_n = ls.x(idx_trp_n);
                        idx_trp_e = ls.class_par == LS_Manipulator_new.PAR_TROPO_N & idx_rec;
                        
                        tropo_e = ls.x(idx_trp_e);
                        tropo_dt = rem(this.common_time.getNominalTime(ls.obs_rate) - ls.getTimePar(idx_trp_n).minimum, this.state.rate_grad_net)/ this.state.rate_grad_net;
                        tropo_idx = floor((this.common_time.getNominalTime(ls.obs_rate) - ls.getTimePar(idx_trp_n).minimum)/this.state.rate_grad_net);
                        [~,tropo_idx] = ismember(tropo_idx*this.state.rate_grad_net, ls.getTimePar(idx_trp_n).getNominalTime(ls.obs_rate).getRefTime(ls.getTimePar(idx_trp_n).minimum.getMatlabTime));
                            valid_ep = tropo_idx ~=0 & tropo_idx <= (length(tropo_n)-3);
                        spline_base = Core_Utils.spline(tropo_dt(valid_ep),spline_order);
                        
                        tropo_n =sum(spline_base .* tropo_n(repmat(tropo_idx(valid_ep), 1, spline_order + 1) + repmat((0 : spline_order), numel(tropo_idx(valid_ep)), 1)), 2);
                        tropo_e =sum(spline_base .* tropo_e(repmat(tropo_idx(valid_ep), 1, spline_order + 1) + repmat((0 : spline_order), numel(tropo_idx(valid_ep)), 1)), 2);
                        this.ztd_gn(valid_ep,i) = nan2zero(this.ztd_gn(valid_ep,i))  + tropo_n;
                        this.ztd_ge(valid_ep,i) = nan2zero(this.ztd_ge(valid_ep,i))  + tropo_e;
                        end
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
                        if this.state.flag_ztd_net
                            this.ztd(i,:) = (S*this.ztd(i,:)')';
                        end
                        % gradients
                        if this.state.flag_grad_net
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
                if ~isempty(idx_pos)
                clk_rec = this.rec_list(i).work.getDt();
                this.clock(idx_pos,i) = this.clock(idx_pos,i) + clk_rec(idx_is);
                
                if this.state.flag_ztd_net
                    ztd_rec = this.rec_list(i).work.getZtd();
                    ztd_rec_apr = this.rec_list(i).work.getZwd() + this.rec_list(i).work.getAprZhd();
                    ztd_rec(ztd_rec == 0) = ztd_rec_apr(ztd_rec == 0);
                    this.ztd(idx_pos,i) = this.ztd(idx_pos,i) + ztd_rec(idx_is);
                end
                
                if this.state.flag_grad_net
                    [gn_rec, ge_rec] = this.rec_list(i).work.getGradient();
                    
                    this.ztd_gn(idx_pos,i) = this.ztd_gn(idx_pos,i) + gn_rec(idx_is);
                    
                    this.ztd_ge(idx_pos,i) = this.ztd_ge(idx_pos,i) + ge_rec(idx_is);
                end
                end
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
                if this.state.flag_ztd_net
                    ztd = this.ztd(idx_res_av, i);
                    this.rec_list(i).work.ztd(idx_pos) = ztd(idx_is);
                    %zhd = this.rec_list(i).work.getAprZhd();
                    this.rec_list(i).work.zwd(idx_pos) = ztd(idx_is) - this.rec_list(i).work.apr_zhd(idx_pos);
                end
                if this.state.flag_grad_net
                    gn = this.ztd_gn(idx_res_av, i);
                    this.rec_list(i).work.tgn(idx_pos) = gn(idx_is);
                    ge = this.ztd_ge(idx_res_av, i);
                    this.rec_list(i).work.tge(idx_pos) = ge(idx_is);
                end
                s0 = mean(abs(ls.res(ls.phase_obs > 0 & ~ls.outlier_obs)));
                % sigma of the session
                this.rec_list(i).work.quality_info.s0 = s0;
                this.rec_list(i).work.quality_info.n_epochs = length(unique(ls.time_par(ls.rec_par == i & ~ls.out_par)));
                idx_obs = ls.receiver_obs == i & ~ls.outlier_obs;
                this.rec_list(i).work.quality_info.n_obs = sum(idx_obs);
                this.rec_list(i).work.quality_info.n_sat = length(unique(ls.satellite_obs(idx_obs)));
                this.rec_list(i).work.quality_info.n_sat_max = max(hist(unique(ls.time_obs.getEpoch(idx_obs).getNominalTime(ls.obs_rate).getRefTime(ls.time_obs.minimum.getMatlabTime) * 1000 + double(ls.satellite_obs(idx_obs))), this.rec_list(i).work.quality_info.n_epochs ));
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
                        
                        [~,idx_time] = ismember(ls.ref_time_obs(idx_res),this.rec_list(i).work.time.getNominalTime.getRefTime(ls.time_min.getMatlabTime));     
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
                        [~,idx_time] = ismember(ls.ref_time_obs(idx_res),this.rec_list(i).work.time.getNominalTime.getRefTime(ls.time_min.getMatlabTime));
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
                        if ~isempty(a_id)
                        a_id = a_id(1);
                        ph(amb_mat == a_id) = ph(amb_mat == a_id) - ls.x(amb)*wl(col_idx);
                        end
                    end
                    this.rec_list(i).work.setPhases(ph,wl,id_ph );
                end
            end
            if sum(ls.param_class == LS_Manipulator_new.PAR_SAT_EB) > 0 && false
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
            fid = fopen(fpath,'Wb');
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
        
    end
end
