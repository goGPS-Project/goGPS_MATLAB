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
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 1 - nightly
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

classdef Observation_Set < handle
    properties
        obs
        obs_code
        wl
        el
        az
        prn
        snr
        cycle_slip
        time
        go_id
    end
    methods 
        function this = Observation_Set(time, obs, obs_code, wl, el, az, prn)
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
            cc = Go_State.getCurrentSettings().getConstellationCollector;
            this.go_id = cc.getIndex(this.obs_code(:,1), this.prn);
            
        end
        function merge(this, obs_set)
            % DESCRIPTION: merge observation with same time stamps
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
            this.wl = [this.wl obs_set.wl];
            this.el = [this.el obs_set.el];
            this.az = [this.az obs_set.az];
            this.prn = [this.prn obs_set.prn];
            this.go_id = [this.go_id obs_set.go_id];
            
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
        function obs = getUnderCutOff(this, cut_off)
            obs = this.obs;
            obs(this.el < cut_off) = 0;
        end
        function remUnderCutOff(this, cut_off)
            idx = this.el < cut_off;
            this.remObs(idx);
        end
        function remObs(this,idx)
            this.obs(idx) = 0;
            this.snr(idx) = 0;
            this.el(idx) = 0;
            this.az(idx) = 0;
            idx2 = find(idx);
            idx2(idx2==(size(idx,1) * size(idx,2))) = [];
            for i = idx2'
                this.cycle_slip(i+1) = this.cycle_slip(i+1) | this.cycle_slip(i);
            end
            this.cycle_slip(idx) = 0;
            this.sanitizeEmpty();
        end
        function sanitizeEmpty(this)
            %remove empty lines
            idx_e = sum(this.obs,2) == 0;
            this.remEpochs(idx_e);
            %remove empty column
            idx_rem_c = not(sum(nan2zero(this.obs)) ~= 0);
            this.removeColumn(idx_rem_c);
        end
        function remEpochs(this, idx_rem)
            this.obs(idx_rem,:) = [];
            this.snr(idx_rem,:) = [];
            this.el(idx_rem,:) = [];
            this.az(idx_rem,:) = [];
            idx_rem_n = find(idx_rem);
            idx_rem_n(idx_rem_n == this.time.length) =[];
            eps = find(idx_rem);
            for e = 1: eps
                this.cycle_slip(e+1,:) = this.cycle_slip(e+1,:) |this.cycle_slip(e,:) ;% bring cycle slips to next epochs
            end
            this.cycle_slip(idx_rem,:) = [];
            this.time = this.time.getSubSet(~idx_rem);

        end
        function idx = getTimeIdx(this,time_st, rate)
            idx = round((this.time -time_st)/rate) +1;
        end
        function removeColumn(this, idx_col)
            this.obs(:,idx_col) = [];
            this.az(:,idx_col) = [];
            this.el(:,idx_col) = [];
            this.cycle_slip(:,idx_col) = [];
            this.snr(:,idx_col) = [];
            this.obs_code(idx_col,:) = [];
            this.wl(idx_col) = [];
            this.prn(idx_col) = [];
            this.go_id(idx_col) = [];
        end
        function plotCycleSlip(this, rec)
            if ~isempty(this.cycle_slip)
                synt_ph = rec.getSyntTwin(this);
                obs = this.obs - synt_ph;
                figure;
                plot(zero2nan(obs));
                ep = repmat([1: this.time.length]',1,size(obs,2));
                hold on
                scatter(ep(this.cycle_slip~=0),obs(this.cycle_slip~=0))
            end
        end
    end
end