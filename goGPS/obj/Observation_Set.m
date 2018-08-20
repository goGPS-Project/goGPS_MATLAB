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
%    |___/                    v 0.6.0 alpha 4 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
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
            cc = Global_Configuration.getCurrentSettings().getConstellationCollector;
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
        
        function remEpochs(this, idx_rem)
            % Remove obseravtions at desidered epoch
            %
            % SYNTAX
            %     this.remEpochs(idx)
            %
            
            if sum(idx_rem) > 0
                this.obs(idx_rem,:) = [];
                if ~isempty(this.el)
                    this.el(idx_rem,:) = [];
                end
                if ~isempty(this.az)
                    this.az(idx_rem,:) = [];
                end
                this.snr(idx_rem,:) = [];
                
                lim = getOutliers(idx_rem);
                if lim(end) == numel(idx_rem)
                    lim(end,:) = [];
                end
                lim(:,2) = lim(:,2) + 1;
                for l = 1 : size(lim, 1)
                    this.cycle_slip(lim(l, 2), :) = any(this.cycle_slip(lim(l, 1) : lim(l, 2), :));
                end
                this.cycle_slip(idx_rem,:) = [];
                                
                this.time = this.time.getSubSet(~idx_rem);
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
                
                lim = getOutliers(idx_rem);
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
            idx = round((this.time -time_st)/rate) +1;
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
                    p_rate = lcm(round(p_rate * 1e6), round(obs_set_list(r).time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
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
                t = intersect(t, (t(1) : p_rate : t(end) + p_rate)');
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
                [p_time, id_sync] = getSyncTimeExpanded(obs_set_list)
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
                obs_c_matrix(id_sync_t, obs_set_list(i).go_id, i) = obs_set_list(i).obs;
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
            this.cycle_slip(:,idx_col) = [];
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
            
            amb_idx = ones(size(this.cycle_slip));
            n_epochs = size(amb_idx,1);
            n_stream = size(amb_idx,2);
            for s = 1:n_stream
                if s > 1
                    amb_idx(:, s) = amb_idx(:, s) + amb_idx(n_epochs, s-1);
                end
                cs = find(this.cycle_slip(:, s) > 0)';
                for c = cs
                    amb_idx(c:end, s) = amb_idx(c:end, s) + 1;
                end
            end
            amb_idx = zero2nan(amb_idx .* (this.obs ~= 0));
            amb_idx = Core_Utils.remEmptyAmbIdx(amb_idx);
        end
        
        function remShortArc(this, min_arc_len)
            % Remove short arcs
            %
            % SYNTAX:
            % this.remShortArc(min_arc_len)
            if nargin < 2
                min_arc_len = 2;
            end
            amb_idx = this.getAmbIdx();
            for i = 1 : size(this.obs, 2)
                [lim] = getOutliers(isnan(amb_idx(:,i)));
                single_arcs = find((lim(:,2) - lim(:,1)) < (min_arc_len - 1));
                for s = 1 : numel(single_arcs)
                    this.obs(lim(single_arcs(s),1) : lim(single_arcs(s),2), i) = 0;
                end
                
                [lim] = getOutliers(this.obs(:,i)~=0);
                single_arcs = find((lim(:,2) - lim(:,1)) < (min_arc_len - 1));
                for s = 1 : numel(single_arcs)
                    this.obs(lim(single_arcs(s),1) : lim(single_arcs(s),2), i) = 0;
                end
            end
        end
    end
end
