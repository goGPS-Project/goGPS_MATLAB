%   CLASS Observation_set
% =========================================================================
%
% DESCRIPTION
%   observation set to be used for calalculations
%
% REFERENCES
%


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:        Giulio Tagliaferro
%  Contributors:      Giulio Tagliaferro, Andrea Gatti
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

classdef Observation_Set < handle
    properties
        time     % GPS_Time [n_epochs x 1]
        obs      % the actual observatiions [n_epochs x n_obs_type]
        obs_code % the code for each observation [n_obs_type x max_n_char_obs_code] , general obs_code structure [sys xxx yyy comb_type], exmaple [GC1WC2WI]
                 % comb_type : I -> iono free
                 %             W -> widelane
                 %             N -> narrowlane
                 %             G -> geometry free
                 %             M -> melbourne wubbena
        wl       % wavelength of the observations (only for phase measurements)
        el       % elevation 
        az       % azimuth
        prn      % prn
        snr      % signla to naoise ration
        cycle_slip % cycle slip index [ n_epochs x n_obs_type] (sparse) 1 cycle slip 0 no cycle slip        
        go_id    % go_ids of the observations
        sigma    % teoretical precision of the measurements [m]
        iono_free%
    end
    
    methods
        
        function this = Observation_Set(time, obs, obs_code, wl, el, az, prn)
            % SYNTAX: 
            %   obs_set = Observation_Set(time, obs, obs_code, wl, el, az, prn)
            if nargin == 0
            return
            end
            this.time = time;
            this.obs = obs;
            this.obs_code = obs_code;
            this.wl = wl;
            this.el = el;
            this.az = az;
            this.prn = prn;
            cc = Core.getCurrentSettings().getConstellationCollector;
            this.go_id = cc.getIndex(this.obs_code(:,1), this.prn);
                        
        end
        
        function merge(this, obs_set)
            % Merge observation with same time stamps
            % 
            % SYNTAX:
            %   this.merge(obs_set)
            %
            %cases empty object
            if isempty(this.time) %case void set
                this.time = obs_set.time;
            end
            %case optional field empty
            if isempty(this.snr)
                this.snr = nan(size(this.obs));
            end
            if isempty(this.cycle_slip)
                this.cycle_slip = sparse(size(this.obs,1),size(this.obs,2));
            end
            %merge
            this.obs = [this.obs obs_set.obs];
            this.obs_code = [this.obs_code ;obs_set.obs_code];
            this.wl = [this.wl(:); obs_set.wl(:)];
            this.el = [this.el obs_set.el];
            this.az = [this.az obs_set.az];
            this.prn = [this.prn(:); obs_set.prn(:)];
            this.go_id = [this.go_id(:); obs_set.go_id(:)];
            this.sigma = [this.sigma(:); obs_set.sigma(:)];
            
            if isempty(obs_set.snr)
                snr2 = nan(size(obs_set.obs));
            else
                snr2 = obs_set.snr;
            end
            this.snr = [this.snr snr2];
            
            if isempty(obs_set.cycle_slip)
                cycle_slip2 = sparse(size(obs_set.obs,1),size(obs_set.obs,2));
            else
                 cycle_slip2 = obs_set.cycle_slip;
            end
            this.cycle_slip = [this.cycle_slip cycle_slip2];
        end
                
        function remUnderCutOff(this, cut_off)
            % Remove observations under the selcted cut off
            %
            % SYNTAX : 
            %   this.remUnderCutOff(cut_off);
            %
            idx = this.el < cut_off;
            this.remObs(idx);
        end
        
        function remUnderSnrThr(this, snr_thr)
            % Remove observations under the selcted snr threshold
            %
            % SYNTAX : 
            %   this.remUnderSnrThr(snr_thr);
            %
            idx = this.snr < snr_thr;
            this.remObs(idx);
        end
        
        function remObs(this, idx,sanitize)
            % Remove the observations identified by the index
            % idx, remove all corrsponding paramaters ( snr el az cycle
            % slip) then sanitize the object for empty row or columns
            %
            % SYNTAX : 
            %      this.remObs(idx)
            %
            if nargin < 3
                sanitize = true;
            end
            if ~isempty(this.obs)
                this.obs(idx) = 0;
            end
            if ~isempty(this.snr)
                this.snr(idx) = 0;
            end
            if ~isempty(this.el)
                this.el(idx) = 0;
            end
            if ~isempty(this.az)
                this.az(idx) = 0;
            end
            if ~isempty(this.cycle_slip) && sum(this.cycle_slip(idx)) > 0
                for s = 1 : size(idx, 2)
                    idx2 = find(idx(1 : (end - 1),s));
                    for i = idx2(:)'
                        this.cycle_slip(i+1, s) = this.cycle_slip(i+1, s) | this.cycle_slip(i, s); %<= move CS
                    end
                end
                this.cycle_slip(idx) = 0;
            end
            if sanitize
                this.sanitizeEmpty();
            end
        end
        
        function sanitizeEmpty(this)
            %Remove empty lines and rows in obs
            %
            % SYNTAX
            %     this.sanitizeEmpty()
            %
            
            %remove empty lines
            idx_e = sum(this.obs,2) == 0;
            if sum(idx_e)
                this.remEpochs(idx_e);
            end
            %remove empty column
            this.remEmptyColumns();
        end
        
        function remEpochs(this, lid_rem)
            % Remove obseravtions at desidered epoch
            %
            % SYNTAX
            %     this.remEpochs(idx)
            %
            
            if sum(lid_rem) > 0
                this.obs(lid_rem,:) = [];
                if ~isempty(this.el)
                    this.el(lid_rem,:) = [];
                end
                if ~isempty(this.az)
                    this.az(lid_rem,:) = [];
                end
                this.snr(lid_rem,:) = [];
                
                if ~isempty(this.cycle_slip)
                    lim = getFlagsLimits(lid_rem);
                    if lim(end) == numel(lid_rem)
                        lim(end,:) = [];
                    end
                    lim(:,2) = lim(:,2) + 1;
                    for l = 1 : size(lim, 1)
                        this.cycle_slip(lim(l, 2), :) = any(this.cycle_slip(lim(l, 1) : lim(l, 2), :));
                    end
                    this.cycle_slip(lid_rem,:) = [];
                end
                
                this.time = this.time.getSubSet(~lid_rem);
            end
        end
        
        function setZeroEpochs(this, idx_rem)
            % Remove obseravtions at desidered epoch
            %
            % SYNTAX
            %     this.remEpochs(idx)
            %
            
            if sum(idx_rem) > 0
                this.obs(idx_rem,:) = 0;
                if ~isempty(this.el)
                    this.el(idx_rem,:) = 0;
                end
                if ~isempty(this.az)
                    this.az(idx_rem,:) = 0;
                end
                this.snr(idx_rem,:) = 0;
                
                lim = getFlagsLimits(idx_rem);
                if lim(end) == numel(idx_rem)
                    lim(end,:) = [];
                end
                lim(:,2) = lim(:,2) + 1;
                for l = 1 : size(lim, 1)
                    this.cycle_slip(lim(l, 2), :) = any(this.cycle_slip(lim(l, 1) : lim(l, 2), :));
                end
                this.cycle_slip(idx_rem,:) = 0;
            end
        end
        
        function keepEpochs(this, idx)
            % Remove all onservations not at epochs than are not idx
            %
            % SYNTAX
            %     this.keepEpochs(idx)
            %
            idx_rem = true(this.time.length, 1);
            idx_rem(idx) = false;
            this.remEpochs(idx_rem);
        end
        
        function idx = getTimeIdx(this, time_st, rate)
            % Using a start time and a rate return at which integer
            % multiple o rate the observations are closer
            %
            % SYNTAX
            %     idx = this.getTimeIdx(this,time_st, rate)
            %
            if nargin == 3
                idx = round((this.time - time_st) / rate) +1;
            else % assume time_st is the full time
                rate = time_st.getRate();
                time_rec = time_st.getNominalTime(rate).getRefTime(round(time_st.first.getMatlabTime));
                time_obs = this.time.getNominalTime(rate).getRefTime(round(time_st.first.getMatlabTime));
                [~, idx] = intersect(round(time_rec * rate), round(time_obs * rate));                
            end
        end
        
        function has_phase = hasPhase(this)
            % tell if there are pahse measuremt
            %
            % SYNTAX:
            %    has_phase = hasPhase(this)
            has_phase = sum(this.wl ~= -1)>0;
        end
        
        function obs_cy = getObsCy(this, idx)
            % Get the observation in cycle
            % (Instead of meters)
            %
            % SYNTAX:
            % obs_cy = this.getObsCy(<idx>)
            if nargin < 2
                idx = 1:length(this.wl);
            end
            obs_cy = this.obs(:,idx) ./ repmat(this.wl(idx), size(this.obs,1),1);
        end
        
        function amb_mat = getRoundedAmb(this, idx)
            % Get the rounded cycle
            %
            % SYNTAX:
            % obs_cy = this.getObsCy(<idx>)
            if nargin < 2
                idx = 1:length(this.wl);
            end
            obs_cy = this.getObsCy(idx);
            amb_idx = this.getAmbIdx();
            
            amb_mat = nan(size(obs_cy));
            for i = 1 : max(amb_idx(:,end))
                idx_amb = amb_idx == i;
                amb = round(mean(obs_cy(idx_amb),'omitnan'));
                amb_mat(idx_amb) = amb;
            end
        end
        
        function [p_time, id_sync] = getSyncTimeExpanded(obs_set_list, p_rate)
            % Get the common time among all the observation_sets
            %
            % SYNTAX
            %   [p_time, id_sync] = getSyncTimeExpanded(obs_set_list, <p_rate>, <use_pos_time>);
            %
            % EXAMPLE:
            %   [p_time, id_sync] = obs_set_list.getSyncTimeExpanded(30);
            
            
            if nargin < 2 || isempty(p_rate)
                p_rate = 1e-6;
                
                for r = 1 : numel(obs_set_list)
                    if obs_set_list(r).time.isEmpty
                        p_rate = inf;
                    else
                        p_rate = lcm(round(p_rate * 1e6), round(obs_set_list(r).time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
                    end
                end
            end
            
            % prepare reference time
            % processing time will start with the obs_set with the last first epoch
            %          and it will stop  with the obs_set with the first last epoch
            
            % first_id_ok = find(~obs_set_list.isEmpty_mr, 1, 'first');
            first_id_ok = 1;
            if ~isempty(first_id_ok)
                p_time_zero = round(obs_set_list(first_id_ok).time.first.getMatlabTime() * 24)/24; % get the reference time
            end
            
            % Get all the common epochs
            t = [];
            for r = 1 : numel(obs_set_list)
                rec_rate = min(1, obs_set_list(r).time.getRate);
                t = [t; round(obs_set_list(r).time.getRefTime(p_time_zero) / rec_rate) * rec_rate];
                % p_rate = lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
            end
            t = unique(t);
            
            % If p_rate is specified use it
            if nargin > 1
                if ~isempty(t)
                    t = intersect(t, (t(1) : p_rate : t(end) + p_rate)');
                end
            end
            
            % Create reference time
            p_time = GPS_Time(p_time_zero, t);
            id_sync = nan(p_time.length(), numel(obs_set_list));
            
            % Get intersected times
            for r = 1 : numel(obs_set_list)
                rec_rate = min(1, obs_set_list(r).time.getRate);
                [~, id1, id2] = intersect(t, round(obs_set_list(r).time.getRefTime(p_time_zero) / rec_rate) * rec_rate);
                id_sync(id1, r) = id2;
            end
        end
        
        function [obs_c_matrix] = getMRObsMat(obs_set_list, p_time, id_sync)
            % get the common obesrvation matrix for all eleent of the reciever list,
            % first direction time, second direction satellites. third direction receiver
            %
            % SYNTAX:
            %   [obs_c_matrix] = getMRObsMat(obs_set_list, p_time, id_sync)
            if nargin < 2
                [p_time, id_sync] = getSyncTimeExpanded(obs_set_list);
            end
            % finde the unique vectors of go id
            n_r = length(obs_set_list);
            go_ids = [];
            for i = 1 : n_r
                go_ids = [go_ids; obs_set_list(i).go_id(:)];
            end
            n_s = max(go_ids);
            %initialize the matrix
            obs_c_matrix = nan(p_time.length, n_s, n_r);
            % fill it
            for i = 1 : n_r
                id_sync_t = id_sync(:,i);
                id_sync_t(isnan(id_sync_t)) = [];
                obs_c_matrix(~isnan(id_sync(:,i)), obs_set_list(i).go_id, i) = obs_set_list(i).obs(id_sync_t,:);
            end
            % tranform zeros to nan
            obs_c_matrix = zero2nan(obs_c_matrix);
        end
        
        function removeColumn(this, idx_col)
            % Remove colums from observations
            %
            % SYNTAX:
            %   this.removeColumn(idx)
            %
            this.obs(:,idx_col) = [];
            if ~isempty(this.el)
                this.el(:,idx_col) = [];
            end
            if ~isempty(this.az)
                this.az(:,idx_col) = [];
            end
            if ~isempty(this.cycle_slip)
                this.cycle_slip(:,idx_col) = [];
            end
            this.snr(:,idx_col) = [];
            this.obs_code(idx_col,:) = [];
            this.wl(idx_col) = [];
            this.prn(idx_col) = [];
            this.go_id(idx_col) = [];
            this.sigma(idx_col) = [];
        end
        
        function remEmptyColumns(this)
           % Remove empty colums from observations
            %
            % SYNTAX:
            %   this.removeEmptyColumn()
            %
            if ~isempty(this.obs)
                idx_rem_c = not(sum(nan2zero(this.obs)) ~= 0);
                if sum(idx_rem_c)
                    this.removeColumn(idx_rem_c);
                end
            end
        end
        
        function amb_idx = getAmbIdx(this)
            % get matrix of same dimesion of the observation showing the ambiguity index of the obsarvation
            %
            % SYNTAX:
            % this.getAmbIdx()
            amb_idx = Core_Utils.getAmbIdx(this.cycle_slip(:,this.wl ~= -1), this.obs(:,this.wl ~= -1));
        end
        
        function n_obs = getNumObs(this)
            % get totoal number of observations
            %
            % SYNTAX:
            %  this.getNumObs()
            n_obs = 0;
            for r = 1 : length(this)
                n_obs = n_obs + sum(sum(this.obs~=0 & ~isnan(this.obs))); 
            end
        end
        
        function arc_jmp_mat = getArcJmpMat(this, id_comm)
            % get a matrix of the dimension of observation set that has true value if the ambiguity is jumping and false value if the ambiguity is not jumping
            %
            % SYNTAX:
            % arc_jmp_mat = this.getArcJmpMat(<id_comm>)
            if nargin < 2
                arc_jmp_mat = false(size(this.obs));
                cycle_slip = this.cycle_slip;
                obs = this.obs;
            else
                % rfer everything to the common index
                arc_jmp_mat = false(length(id_comm),size(this.obs,2));
                cycle_slip = false(length(id_comm),size(this.obs,2));
                obs = false(length(id_comm),size(this.obs,2));
                id_sync = id_comm(~isnan(id_comm));
                cycle_slip(~isnan(id_comm),:) = this.cycle_slip(id_sync,:);
                obs(~isnan(id_comm),:) = this.obs(id_sync,:);
            end
            ne = size(obs,1);
            for s = 1 : size(obs,2)
                css = find(cycle_slip(:,s));
                for cs = css'
                    % marks as jump the epoch of the cycle slip and all the epochs before the cycle slip that have no obersvatiobn
                    arc_jmp_mat(cs,s) = true;
                    it = cs -1;
                    while it > 0 && obs(it,s) ==0
                        arc_jmp_mat(it,s) = true;
                        it = it -1;
                    end
                    
                end
                % mark as jmp all the epoch after the last valid one
                it = ne;
                while it > 0 && obs(it,s) == 0
                        arc_jmp_mat(it,s) = true;
                        it = it -1;
                end
            end
        end
        
        function remShortArc(this, min_arc_len)
            % Remove short arcs
            %
            % SYNTAX:
            % this.remShortArc(min_arc_len)
            if nargin < 2
                min_arc_len = 2;
            end
            cs = full(this.cycle_slip);
            cs(this.obs == 0) = true;
            n_ep = size(this.obs,1);
            id_rm = false(size(this.obs));
            for i = 1 : size(this.obs, 2)
                % remove single arcs
                for j = 2:(n_ep)
                    if cs(j, i) && cs(j-1, i) 
                        id_rm(j-1, i) = true;
                    end
                end
                if cs(j, i) & cs(j-1, i)
                    id_rm(j, i) = true;
                end
                % remove arcs less than desidered
                [lim] = getFlagsLimits(this.obs(:, i)~=0 & ~this.cycle_slip(:, i));
                if ~isempty(lim)
                    single_arcs = find((lim(:,2) - lim(:,1)) < (min_arc_len - 1));
                    for s = 1 : numel(single_arcs)
                        id_rm(lim(single_arcs(s),1) : lim(single_arcs(s),2), i) = true;
                    end
                end

            end
            id_rm(this.obs == 0) = false; %do not remove observation that are not there
            this.remObs(id_rm);
        end
        
        function is_empty = isEmpty(this)
            % Check if the object is empty
            %
            % SYNTAX:
            %     is_empty = this.isEmpty()
            is_empty =  isempty(this.obs);
        end
    end
end
