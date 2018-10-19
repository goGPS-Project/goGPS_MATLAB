%   CLASS Core_SEID
% =========================================================================
%
% DESCRIPTION
%   class to compute a SEID processing
%   Satellite specific Epoch differenced Ionospheric Delay model
%
% EXAMPLE
%   seid = Core_SEID();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_SEID


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

classdef Core_SEID < handle

    properties % Public Access
        log;
        state;
    end

    methods (Static)
        function this = Core_SEID()
            % Core object creator
            this.log = Logger.getInstance();
            this.state = Core.getState();
        end

        % Concrete implementation.  See Singleton superclass.
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_core_seid_
            unique_instance_core_seid_ = [];

            if isempty(unique_instance_core_seid_)
                this = Core();
                unique_instance_core_seid_ = this;
            else
                this = unique_instance_core_seid_;
            end
        end
    end

    % =========================================================================
    %  MAIN STATIC FUNCTION
    % =========================================================================

    methods (Static) % Public Access
        function getSyntL2(ref, trg)
            % Compute L2 synthetic observations to be injected into the receiver
            %
            % SYNTAX
            %   Core_SEID.getSyntL2(ref, trg)
            %%
            trg = trg(~trg.isEmpty_mr);
            ref = ref(~ref.isEmpty_mr);
            if ~isempty(trg) && ~isempty(ref)                
                rec(1:numel(ref)) = ref;
                rec(numel(ref) + (1 : numel(trg))) = trg;
                obs_type(1:numel(ref)) = 2;
                obs_type(numel(ref) + (1 : numel(trg))) = 0;
                [p_time, id_sync] = Receiver_Commons.getSyncTimeTR(rec, obs_type);
                log = Logger.getInstance();
                
                log.addMarkedMessage('Starting SEID processing')
                log.addMessage(log.indent('Getting Geometry free from reference receivers'));
                systems = unique(ref(1).system);
                for r = 1 : numel(ref)
                    phase_gf(r) = Observation_Set();
                    code_gf(r) = Observation_Set();
                    for sys = systems
                        phase_gf(r).merge(ref(r).getGeometryFree('L1','L2',sys));
                        code_gf(r).merge(ref(r).getGeometryFree('C1','C2',sys));
                    end
                    phase_gf(r).obs(phase_gf(r).cycle_slip > 0) = 0;
                    phase_gf(r).obs = ref(r).smoothSatData([], [], zero2nan(phase_gf(r).obs), phase_gf(r).cycle_slip, [], 300 / phase_gf(r).time.getRate); 
                    
                    %[phase_gf(r).obs, phase_gf(r).sigma] = Receiver.smoothCodeWithPhase(zero2nan(code_gf(r).obs), code_gf(r).sigma, code_gf(r).go_id, ...
                    %zero2nan(phase_gf(r).obs), phase_gf(r).sigma, phase_gf(r).go_id, phase_gf(r).cycle_slip);
                    
                    % Smoothing iono
                    %phase_gf(r).obs = ref(r).smoothSatData([], [], zero2nan(phase_gf(r).obs), phase_gf(r).cycle_slip, 'spline', (300 / ref(r).getRate));
                    %code_gf(r).obs = ref(r).smoothSatData([], [], zero2nan(code_gf(r).obs), [], 'spline', (300 / ref(r).getRate), (300 / ref(r).getRate));
                    %code_gf(r).obs = phase_gf(r).obs;
                    
                    [lat, lon, ~, h_ortho] = rec(r).getMedianPosGeodetic;
                    [pierce_point(r).lat, pierce_point(r).lon, pierce_point(r).mf] = Atmosphere.getPiercePoint(lat / 180 * pi, lon / 180 * pi, h_ortho, code_gf(r).az / 180 * pi, zero2nan(code_gf(r).el / 180 * pi), 350*1e3);
                end
                
                max_sat = 0;
                for r = 1 : numel(ref)
                    max_sat = max(max_sat, max(code_gf(r).go_id));
                end
                
                min_gap = 3;
                
                % Extract syncronized C4 L4 diff
                for t = 1 : numel(trg)
                    log.addMessage(log.indent(sprintf('Computing interpolated geometry free for target %d / %d', t, numel(trg))));
                    
                    max_sat_trg = max(max_sat, max(trg(t).go_id));
                    
                    ph_gf = nan(size(id_sync{t}, 1), max_sat_trg, numel(ref));
                    pr_gf = nan(size(id_sync{t}, 1), max_sat_trg, numel(ref));
                    for r = 1 : numel(ref)
                        ph_gf(:, phase_gf(r).go_id, r) = zero2nan(phase_gf(r).obs(id_sync{t}(:,r), :));
                        % Import CS and outliers from receivers
                        %                     for s = 1 : numel(ref(r).ph_idx)
                        %                         %ph_gf(find(ref(r).outlier_idx_ph(id_sync{t}(:,r),s)), ref(r).go_id(ref(r).ph_idx(s)), r) = nan;
                        %                         %ph_gf(find(ref(r).cycle_slip_idx_ph(id_sync{t}(:,r),s)), ref(r).go_id(ref(r).ph_idx(s)), r) = nan;
                        %
                        %                         % fill small gaps
                        %                         lim = getOutliers(isnan(ph_gf(:, ref(r).go_id(s), r)));
                        %                         lim(lim(:,2) - lim(:,1) > min_gap,:) = [];
                        %                         idx = false(size(ph_gf, 1), 1);
                        %                         for l = 1 : size(lim, 1)
                        %                             idx(lim(l, 1) : lim(l, 2)) = true;
                        %                         end
                        %                         if sum(idx) > 0
                        %                             ph_gf(:,ref(r).go_id(s), r) = simpleFill1D(ph_gf(:, ref(r).go_id(s), r), idx);
                        %                         end
                        %                     end
                        pr_gf(:, code_gf(r).go_id, r) = zero2nan(code_gf(r).obs(id_sync{t}(:,r), :));
                    end
                    ph_gf_diff = diff(ph_gf);
                    
                    % ph_gf pr_gf ph_gf_diff have max_sat satellite data stored
                    % pierce_point(r).lat could have different size receiver by receiver
                    % indexes convarsion is id = phase_gf(r).go_id
                    
                    % % DEBUG: plot
                    % hold off;
                    % for r = 1 : numel(ref)
                    %     % Id of non nan values
                    %     id_ok = reshape(~isnan(pierce_point(r).lon(:)) & ~isnan(serialize(pr_gf(:, phase_gf(r).go_id, r))), size(pierce_point(r).lat, 1), size(pierce_point(r).lat, 2));
                    %     tmp = ph_gf_diff(:, phase_gf(r).go_id, r);
                    %     id_ok(end, :) = false;
                    %     lat_lim = minMax(pierce_point(r).lat / pi * 180); lat_lim(1) = lat_lim(1) - 0.5; lat_lim(2) = lat_lim(2) + 0.5;
                    %     lon_lim = minMax(pierce_point(r).lon / pi * 180); lon_lim(1) = lon_lim(1) - 0.5; lon_lim(2) = lon_lim(2) + 0.5;
                    %     prettyScatter(tmp(id_ok(2 : end, :)), pierce_point(r).lat(id_ok) / pi * 180, pierce_point(r).lon(id_ok) / pi * 180, lat_lim(1), lat_lim(2), lon_lim(1), lon_lim(2), '10m'); hold on; colormap(jet);
                    % end
                    ph1 = [];
                    id_ph = [];
                    for sys = systems
                        [ph1_t, id_ph_t] = trg(t).getObs('L1',sys);
                        ph1 = [ph1; ph1_t];
                        id_ph = [id_ph; id_ph_t];
                    end
                    [lat, lon, ~, h_ortho] = trg(t).getMedianPosGeodetic;
                    ph1_goid = trg(t).go_id(id_ph)';
                    trg_go_id = unique(ph1_goid);
                    [lat_pp, lon_pp, iono_mf] = Atmosphere.getPiercePoint(lat / 180 * pi, lon / 180 * pi, h_ortho, trg(t).sat.az(:, trg_go_id) / 180 * pi, zero2nan(trg(t).sat.el(:, trg_go_id) / 180 * pi), 350*1e3);
                    
                    % It is necessary to better sync satellites in view
                    % this part of the code needs to be improved
                    trg_pr_gf = nan(trg(t).length, max(trg_go_id));
                    trg_ph_gf = nan(trg(t).length, max(trg_go_id));
                    for s = 1 : numel(trg_go_id)
                        lat_sat = nan(size(id_sync{t},1), numel(ref));
                        lon_sat = nan(size(id_sync{t},1), numel(ref));
                        for r = 1 : numel(ref)
                            id_sat = unique(code_gf(r).go_id) == trg_go_id(s);
                            if sum(id_sat) == 1
                                lat_sat(:, r) = pierce_point(r).lat(id_sync{t}(:,r), id_sat);
                                lon_sat(:, r) = pierce_point(r).lon(id_sync{t}(:,r), id_sat);
                            end
                        end
                        trg_pr_gf(id_sync{t}(:,t + numel(ref)), trg_go_id(s)) = Core_SEID.satDataInterp(lat_sat, lon_sat, squeeze(pr_gf(:,trg_go_id(s),:)), lat_pp(id_sync{t}(:,t + numel(ref)), s), lon_pp(id_sync{t}(:,t + numel(ref)), s));
                        %trg_ph_gf(id_sync{t}(:,t + numel(ref)), trg_go_id(s)) = Core_SEID.satDataInterp(lat_sat, lon_sat, squeeze(ph_gf(:,trg_go_id(s),:)), lat_pp(id_sync{t}(:,t + numel(ref)), s), lon_pp(id_sync{t}(:,t + numel(ref)), s));
                        trg_ph_gf(id_sync{t}(2 : end, t + numel(ref)), trg_go_id(s)) = Core_SEID.satDataInterp(lat_sat(2 : end, :), lon_sat(2 : end, :), squeeze(ph_gf_diff(:,trg_go_id(s),:)),  lat_pp(id_sync{t}(2 : end,t + numel(ref)), s), lon_pp(id_sync{t}(2 : end,t + numel(ref)), s));
                    end
                    
                    % Interpolate the diff (derivate) of L4, now rebuild L4 by cumsum (integral)
                    
                    trg_ph_gf(abs(trg_ph_gf) > 0.5) = nan; % remove outliers
                    inan = isnan(trg_ph_gf);

                    %% experimental 
%                     for i = trg(t).go_id(id_ph)
%                         trg_ph_gf(:,i) = simpleFill1D(trg_ph_gf(:,i),isnan(trg_ph_gf(:,i)),'linear');
%                     end
                    trg_ph_gf(abs(trg_ph_gf) > 0.5) = nan; % remove outliers
                    inan = isnan(trg_ph_gf);
                    trg_ph_gf = cumsum(nan2zero(trg_ph_gf));
                    
                    trg_ph_gf(inan) = nan;
                    
                    wl1 = trg(t).state.getConstellationCollector().gps.L_VEC(1);
                    wl2 = trg(t).state.getConstellationCollector().gps.L_VEC(2);
                    ph2 = nan(size(ph1));
                    % L1 * wl1 - L2 * wl2 = gf
                    % L2 = (L1 * wl1 - gf) / wl2;
                    ph2 = (ph1 * wl1 - trg_ph_gf(:, trg(t).go_id(id_ph))') / wl2;
                    
                    [~, ~, ~, flag] = trg(t).getBestCodeObs();
                    pr1 = [];
                    id_pr = [];
                    for sys = systems
                        [pr1_t, id_pr_t] = trg(t).getObs(flag(1,1:3),sys);
                        pr1 = [pr1; pr1_t];
                        id_pr = [id_pr; id_pr_t];
                    end
                    pr1_goid = trg(t).go_id(id_pr);
                    % C2 - C1 = gf
                    % C2 = C1 + gf
                    pr2 = pr1 + trg_pr_gf(:, trg(t).go_id(id_pr))';
                    
                    %compute ~L2
                    %fix_til_L2(PRN,idx_diff_L4) = (L1{target_sta}(PRN,idx_diff_L4)*lambda(PRN,1) - satel(PRN).til_L4(idx_diff_L4))/lambda(PRN,2);
                    
                    %compute ~P2
                    %fix_til_P2(PRN,idx_diff_L4) = P1{target_sta}(PRN,idx_diff_L4) + satel(PRN).til_P4(idx_diff_L4);
                    
                    % Remove the L2 stored in the object
                    id_ph = [];
                    for sys = systems
                        [~, id_ph_t] = trg(t).getObs('L2',sys);
                        id_ph = [id_ph; id_ph_t];
                    end
                    if ~isempty(id_ph)
                        log.addMessage(log.indent(sprintf('Removing L2 observations already present in the target receiver %d / %d', t, numel(trg))));
                        trg(t).remObs(id_ph);
                    end
                    id_pr = [];
                    for sys = systems
                        [~, id_pr_t] = trg(t).getObs('C2',sys);
                        id_pr = [id_pr; id_pr_t];
                    end
                    if ~isempty(id_pr)
                        log.addMessage(log.indent(sprintf('Removing C2 observations already present in the target receiver %d / %d', t, numel(trg))));
                        trg(t).remObs(id_pr);
                    end
                    
                    % Inject the new synthesised phase
                    log.addMessage(log.indent(sprintf('Injecting SEID L2 into target receiver %d / %d', t, numel(trg))));
                    trg(t).injectObs(nan2zero(pr2), wl2, 2, 'C2F', pr1_goid);
                    trg(t).injectObs(nan2zero(ph2), wl2, 2, 'L2F', ph1_goid);
                    %trg(t).injectObs(nan2zero(ref(1).getObs('C2')), wl2, 2, 'C2 ', trg_go_id);
                    %trg(t).injectObs(nan2zero(ref(1).getObs('L2')), wl2, 2, 'L2 ', trg_go_id);
                    
                    trg(t).keepEpochs(id_sync{t}(:,t + numel(ref)));
                    trg(t).updateDetectOutlierMarkCycleSlip();
                end
                
                log.addMarkedMessage('Syncing times, computing reference time');
            end
        end               
        
        function remIono(ref, trg)
            % Compute Ionosphere from reference receivers and remove it from the observations in target
            %
            % SYNTAX
            %   Core_SEID.remIono(ref, trg)
            %%
            trg = trg(~trg.isEmpty_mr);
            if ~isempty(trg)
                rec(1:numel(ref)) = ref;
                rec(numel(ref) + (1 : numel(trg))) = trg;
                obs_type(1:numel(ref)) = 2;
                obs_type(numel(ref) + (1 : numel(trg))) = 0;
                [p_time, id_sync] = Receiver_Commons.getSyncTimeTR(rec, obs_type);
                log = Logger.getInstance();
                
                log.addMarkedMessage('Starting SEID processing')
                log.addMessage(log.indent('Getting Geometry free from reference receivers'));
                
                sys_c = 'G';
                for r = 1 : numel(ref)
                    % comine code and phase
                    iono_ref(r) = ref(r).getPrefGeometryFree('L',sys_c);
                    gf_pr = ref(r).getPrefGeometryFree('C',sys_c);
                    idx_nan = iono_ref(r).obs == 0;

                    el = iono_ref(r).el / 180 * pi;
                    az = iono_ref(r).az / 180 * pi;
                    [iono_ref(r).obs] = ref(r).ionoCodePhaseSmt(zero2nan(gf_pr.obs), gf_pr.sigma.^2, zero2nan(iono_ref(r).obs), iono_ref(r).sigma.^2, iono_ref(r).getAmbIdx(), 0.01, el);
                    
                    % Apply mapping function
                    [lat, lon, ~, h_ortho] = ref(r).getMedianPosGeodetic;
                    [pierce_point(r).lat, pierce_point(r).lon, iono_ref(r).mf] = Atmosphere.getPiercePoint(lat / 180 * pi, lon / 180 * pi, h_ortho, az, zero2nan(el), 350*1e3);
                    iono_ref(r).obs(idx_nan) = nan;
                    %iono_ref(r).obs = ref(r).applyMF(iono_ref(r).obs, iono_ref(r).mf, 1);
                    %iono_ref(r).obs(idx_nan) = nan;

                    iono_ref(r).obs = ref(r).smoothSatData([], [], zero2nan(iono_ref(r).obs), false(size(iono_ref(r).cycle_slip)), [], 300 / iono_ref(r).time.getRate); % <== supposing no more cycle slips
                end
                
                max_sat = 0;
                for r = 1 : numel(ref)
                    max_sat = max(max_sat, max(iono_ref(r).go_id));
                end
                                
                % Extract syncronized C4 L4 diff
                for t = 1 : numel(trg)
                    log.addMessage(log.indent(sprintf('Computing interpolated geometry free for target %d / %d', t, numel(trg))));
                    
                    iono_sync = nan(size(id_sync{t}, 1), max_sat, numel(ref));
                    for r = 1 : numel(ref)
                        iono_sync(:, iono_ref(r).go_id, r) = zero2nan(iono_ref(r).obs(id_sync{t}(:,r), :));% - mean(iono_ref(r).obs(id_sync{t}(:,r), :), 1, 'omitnan');
                    end
                    
                    % DIFF: 
                    iono_diff = diff(nan2zero(iono_sync));
                    
                    % ph_gf pr_gf ph_gf_diff have max_sat satellite data stored
                    % pierce_point(r).lat could have different size receiver by receiver
                    % indexes convarsion is id = phase_gf(r).go_id
                    
                    % % DEBUG: plot
                    % hold off;
                    % for r = 1 : numel(ref)
                    %     % Id of non nan values
                    %     id_ok = reshape(~isnan(pierce_point(r).lon(:)) & ~isnan(serialize(pr_gf(:, phase_gf(r).go_id, r))), size(pierce_point(r).lat, 1), size(pierce_point(r).lat, 2));
                    %     tmp = ph_gf_diff(:, phase_gf(r).go_id, r);
                    %     id_ok(end, :) = false;
                    %     lat_lim = minMax(pierce_point(r).lat / pi * 180); lat_lim(1) = lat_lim(1) - 0.5; lat_lim(2) = lat_lim(2) + 0.5;
                    %     lon_lim = minMax(pierce_point(r).lon / pi * 180); lon_lim(1) = lon_lim(1) - 0.5; lon_lim(2) = lon_lim(2) + 0.5;
                    %     prettyScatter(tmp(id_ok(2 : end, :)), pierce_point(r).lat(id_ok) / pi * 180, pierce_point(r).lon(id_ok) / pi * 180, lat_lim(1), lat_lim(2), lon_lim(1), lon_lim(2), '10m'); hold on; colormap(jet);
                    % end
                    
                    [ph1, id_ph] = trg(t).getObs('L1','G');
                    [lat, lon, ~, h_ortho] = trg(t).getMedianPosGeodetic;
                    ph1_goid = trg(t).go_id(id_ph)';
                    trg_go_id = unique(ph1_goid);
                    [lat_pp, lon_pp, iono_mf] = Atmosphere.getPiercePoint(lat / 180 * pi, lon / 180 * pi, h_ortho, trg(t).sat.az(:, trg_go_id) / 180 * pi, zero2nan(trg(t).sat.el(:, trg_go_id) / 180 * pi), 350*1e3);
                                        
                    % It is necessary to better sync satellites in view
                    iono_trg = nan(trg(t).length, max(trg_go_id));
                    for s = 1 : numel(trg_go_id)
                        lat_sat = nan(size(id_sync{t},1), numel(ref));
                        lon_sat = nan(size(id_sync{t},1), numel(ref));
                        for r = 1 : numel(ref)
                            id_sat = unique(iono_ref(r).go_id) == trg_go_id(s);
                            if sum(id_sat) == 1
                                lat_sat(:, r) = pierce_point(r).lat(id_sync{t}(:,r), id_sat);
                                lon_sat(:, r) = pierce_point(r).lon(id_sync{t}(:,r), id_sat);
                            end
                        end
                        % DIFF: 
                        iono_trg(id_sync{t}(2 : end, t + numel(ref)), trg_go_id(s)) = Core_SEID.satDataInterp(lat_sat(2 : end, :), lon_sat(2 : end, :), squeeze(iono_diff(:,trg_go_id(s),:)),  lat_pp(id_sync{t}(2 : end,t + numel(ref)), s), lon_pp(id_sync{t}(2 : end,t + numel(ref)), s));
                        %iono_trg(id_sync{t}(:, t + numel(ref)), trg_go_id(s)) = Core_SEID.satDataInterp(lat_sat, lon_sat, squeeze(iono_sync(:, trg_go_id(s), :)), lat_pp(id_sync{t}(:, t + numel(ref)), s), lon_pp(id_sync{t}(:, t + numel(ref)), s));
                    end
                    
                    % Interpolate the diff (derivate) of L4, now rebuild L4 by cumsum (integral)
                    iono_trg(abs(iono_trg) > 0.5) = nan; % remove outliers
                    inan = isnan(iono_trg);
                    iono_trg = cumsum(nan2zero(iono_trg));
                    iono_trg(inan) = nan;
                    
                    %iono_trg(:, trg_go_id) = zero2nan(trg(t).applyMF(iono_trg(:, trg_go_id), 1 ./ zero2nan(iono_mf), 1));
                        
                    % Interpolate the diff (derivate) of L4, now rebuild L4 by cumsum (integral)
                    
                    wl1 = trg(t).state.getConstellationCollector().gps.L_VEC(1);
                    wl2 = trg(t).state.getConstellationCollector().gps.L_VEC(2);
                    ph2 = nan(size(ph1));
                    % L1 * wl1 - L2 * wl2 = gf
                    % L2 = (L1 * wl1 - gf) / wl2;
                    ph2 = (ph1 * wl1 - iono_trg(:, trg(t).go_id(id_ph))') / wl2;
                    
                    [~, ~, ~, flag] = trg(t).getBestCodeObs();
                    [pr1, id_pr] = trg(t).getObs(flag(1,1:3),'G');
                    pr1_goid = trg(t).go_id(id_pr);
                    % C2 - C1 = gf
                    % C2 = C1 + gf
                    pr2 = pr1 + iono_trg(:, trg(t).go_id(id_pr))';
                    
                    % Remove the L2 stored in the object
                    [ph_old, id_ph] = trg(t).getObs('L2','G');
                    if ~isempty(id_ph)
                        log.addMessage(log.indent(sprintf('Removing L2 observations already present in the target receiver %d / %d', t, numel(trg))));
                        trg(t).remObs(id_ph);
                    end
                    [pr_old, id_pr] = trg(t).getObs('C2','G');
                    if ~isempty(id_pr)
                        log.addMessage(log.indent(sprintf('Removing C2 observations already present in the target receiver %d / %d', t, numel(trg))));
                        trg(t).remObs(id_pr);
                    end
                    
                    % Inject the new synthesised phase
                    log.addMessage(log.indent(sprintf('Injecting SEID L2 into target receiver %d / %d', t, numel(trg))));
                    trg(t).injectObs(nan2zero(pr2), wl2, 2, 'C2F', pr1_goid);
                    trg(t).injectObs(nan2zero(ph2), wl2, 2, 'L2F', ph1_goid);
                    
                    trg(t).keepEpochs(id_sync{t}(:,t + numel(ref)));
                    trg(t).updateDetectOutlierMarkCycleSlip();
                end       
                log.addMarkedMessage('Syncing times, computing reference time');
            end            
        end
    end

    methods (Static)
        function data_q = satDataInterp(lat_in, lon_in, data_in, lat_q, lon_q, method)
            % Interpolate the data given for n points "_in" to query point "_q"
            %
            % INPUT:
            %   lat_in  [n_epoch x n_ref]  latitude [deg]
            %   lon_in  [n_epoch x n_ref]  longitude [deg]
            %   data_in [n_epoch x n_ref]  data [deg]
            %   lat_q   [n_epoch x 1]      query latitude
            %   lon_q   [n_epoch x 1]      query longitude
            %
            % SYNTAX:
            %   data_q = Core_SEID.satDataInterp(lat_in, lon_in, data_in, lat_q, lon_q);
            
            if nargin < 6
                method = 'distance';
            end
            
            switch method
                case {'plane'}
                    % consider a plane of interpolation, the coordinates of the pierce point must be considered as 
                    % spherical distances in phi/lam from the interpolation point... the reference (center) is the interpolation point
                    % ...to be done
                    % To be tested: regularized solution for plane parameter estimation in time
                   
                    data_q = nan(size(lat_q));
                    for i = 1 : size(lat_in, 1)
                        id_ok = ~isnan(data_in(i, :)) & ~isnan(lat_in(i, :)) & ~isnan(lon_in(i, :));
                        if sum(id_ok) == size(lat_in, 2) % require all the stations to interpolate ionosphere
                             A = ones(sum(id_ok), 3);
                             A(:,2) = lon_in(i, :) .* cos(lat_in(i, :));
                             A(:,3) = lat_in(i, :);
                             data_q(i) = [1 (lon_q(i) .* cos(lat_q(i))) lat_q(i)] * ((A'*A)\A' * data_in(i, id_ok)');
                             % Test diffent reference frame for interpolation
                             %data_q(i) = [1 (lon_q(i)) lat_q(i)] * ((A'*A)\A' * data_in(i, id_ok)');
                             % mybe a 3D Kriging is better
                        end
                    end
                case  {'distance'}
                    % very simple interpolation by spherical distance
                    data_q = nan(size(lat_q));
                    id_ok = ~isnan(data_in) & ~isnan(lat_in) & ~isnan(lon_in);
                    %ep_ok = find(sum(id_ok, 2) == size(lat_in, 2)); % require all the stations to interpolate ionosphere
                    ep_ok = find(sum(id_ok, 2) > 0); % require all the stations to interpolate ionosphere
                    for i = ep_ok'
                        d = sphericalDistance(lat_in(i, :) ./ pi * 180, lon_in(i, :) ./ pi * 180, lat_q(i) ./ pi * 180, lon_q(i) ./ pi * 180);
                        w = (1 ./ (d(id_ok(i,:)) + eps))';
                        data_q(i) = (data_in(i, id_ok(i,:)) * w) ./ sum(w);
                    end

            end
        end
        
    end

end
