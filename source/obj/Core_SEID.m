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
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, Giulio Tagliaferro
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
    
    properties (Constant)
        R_TER = 6378137 % from WGS84
    end

    properties % Public Access
        state;
    end

    methods (Static)
        function this = Core_SEID()
            % Core object creator
            this.state = Core.getState();
        end      
    end

    % =========================================================================
    %  MAIN STATIC FUNCTION
    % =========================================================================

    methods (Static) % Public Access
        function thin_shell_height = getThinShellHeight()
            % Get Thin Shell Height from Core_Athmosphere
            % 
            % SYNTAX
            %   thin_shell_height = getThinShellHeight(this)
            
            thin_shell_height = Core.getAtmosphere.getThinShellHeight();            
        end
            
        function getSyntL2(ref, trg, method, flag_mf)
            % Compute L2 synthetic observations to be injected into the receiver
            %
            % NOTE
            %   Interpolation method based on the plane is unstable 
            %   when data are missing or bad
            %   
            % INPUT
            %   ref      reference receivers
            %   trg      target receivers
            %   method   distance (default) | plane
            %   flag_mf  use ionospheric mapping functions for the interpolation
            %
            % SYNTAX
            %   Core_SEID.getSyntL2(ref, trg)
            %%
            if nargin < 3 || isempty(method)
                method = 'distance';
            end
            if nargin < 4 || isempty(flag_mf)
                flag_mf = false;
            end

            trg = trg(~trg.isEmpty_mr);
            ref = ref(~ref.isEmpty_mr);
            if ~isempty(trg) && ~isempty(ref)                
                % rec = [ref trg];
                rec(1:numel(ref)) = ref;
                rec(numel(ref) + (1 : numel(trg))) = trg;
                % Obs type == 2 means ref
                obs_type(1:numel(ref)) = 2;
                % Obs type == 0 means trg
                obs_type(numel(ref) + (1 : numel(trg))) = 0;
                
                % check if one of the target receiver does have too few epochs
                if ~strcmp(method, 'ls')
                    [rec, ref, obs_type] = Core_SEID.removeRecWithFewObs(ref,rec,obs_type);
                end
                log = Core.getLogger();
               
                % get synced time
                [p_time, id_sync] = Receiver_Commons.getSyncTimeTR(rec, obs_type);
                
                if strcmp(method, 'ls')
                    log.addMarkedMessage('Starting SID processing')
                else
                    log.addMarkedMessage('Starting SEID processing')
                end
                log.addMessage(log.indent('Getting Geometry free from reference receivers'));
                systems = unique(ref(1).system);
                for r = 1 : numel(ref)
                    ph_ref_gf(r) = Observation_Set();
                    pr_ref_gf(r) = Observation_Set();
                    for sys = systems
                        frqs =  num2str(ref(r).getFreqs(sys));
                        if ~isempty(frqs)
                            
                            % Extract geometry free for each reference receiver
                            ph_gf_t = ref(r).getGeometryFree(['L' frqs(1)],['L' frqs(2)],sys);
                            pr_gf_t = ref(r).getGeometryFree(['C' frqs(1)],['C' frqs(2)],sys);
                            
                            % Get wavelength factor
                            wl1_pr = zeros(length(pr_gf_t.go_id),1);
                            wl2_pr = zeros(length(pr_gf_t.go_id),1);
                            for i = 1 : length(pr_gf_t.go_id)
                                l_vec = Core.getConstellationCollector.getWavelength(pr_gf_t.go_id(i));
                                wl1_pr(i) = l_vec(Core.getConstellationCollector.getSys(sys).CODE_RIN3_2BAND == frqs(1));
                                wl2_pr(i) = l_vec(Core.getConstellationCollector.getSys(sys).CODE_RIN3_2BAND == frqs(2));
                            end
                            iono_factor_pr = wl2_pr.^2 - wl1_pr.^2;
                            
                            wl1_ph = zeros(length(ph_gf_t.go_id),1);
                            wl2_ph = zeros(length(ph_gf_t.go_id),1);
                            for i = 1 : length(ph_gf_t.go_id)
                                l_vec = Core.getConstellationCollector.getWavelength(ph_gf_t.go_id(i));
                                wl1_ph(i) = l_vec(Core.getConstellationCollector.getSys(sys).CODE_RIN3_2BAND == frqs(1));
                                wl2_ph(i) = l_vec(Core.getConstellationCollector.getSys(sys).CODE_RIN3_2BAND == frqs(2));
                            end
                            iono_factor_ph = wl2_ph.^2 - wl1_ph.^2;
                            
                            ph_gf_t.obs = ph_gf_t.obs ./ repmat(iono_factor_ph', size(ph_gf_t.obs,1),1);
                            pr_gf_t.obs = pr_gf_t.obs ./ repmat(iono_factor_pr', size(pr_gf_t.obs,1),1);
                            
                            ph_ref_gf(r).merge(ph_gf_t);
                            pr_ref_gf(r).merge(pr_gf_t);
                        end
                    end
                    ph_ref_gf(r).obs(ph_ref_gf(r).cycle_slip > 0) = 0;
                    % Smoothing 
                    %ph_ref_gf(r).obs = ref(r).smoothSatData([], [], zero2nan(ph_ref_gf(r).obs), ph_ref_gf(r).cycle_slip, [], 300 / ph_ref_gf(r).time.getRate); 
                                        
                    % Smooth code with phase (remove biases)
                    %[ph_ref_gf(r).obs, ph_ref_gf(r).sigma] = Receiver_Work_Space.smoothCodeWithPhase(zero2nan(pr_ref_gf(r).obs), pr_ref_gf(r).sigma, pr_ref_gf(r).go_id, ...
                    %                                                                      zero2nan(ph_ref_gf(r).obs), ph_ref_gf(r).sigma, ph_ref_gf(r).go_id, ph_ref_gf(r).cycle_slip);
                    
                    % Smoothing with splines
                    %ph_ref_gf(r).obs = ref(r).smoothSatData([], [], zero2nan(ph_ref_gf(r).obs), ph_ref_gf(r).cycle_slip, 'spline', (300 / ref(r).getRate));
                    %pr_ref_gf(r).obs = ref(r).smoothSatData([], [], zero2nan(pr_ref_gf(r).obs), [], 'spline', (300 / ref(r).getRate), (300 / ref(r).getRate));
                    %pr_ref_gf(r).obs = ph_ref_gf(r).obs;
                    
                    % Get pierce point
                    [lat, lon, ~, h_ortho] = rec(r).getMedianPosGeodetic;
                    n_max_sat = max([ph_ref_gf(r).go_id; pr_ref_gf(r).go_id]);
                    az = nan(size(pr_ref_gf(r).az, 1), n_max_sat);
                    el = nan(size(pr_ref_gf(r).az, 1), n_max_sat);
                    for s = 1 : numel(pr_ref_gf(r).go_id)
                        go_id = pr_ref_gf(r).go_id(s);
                        az(:, go_id) = pr_ref_gf(r).az(:, s);
                        el(:, go_id) = pr_ref_gf(r).el(:, s);
                    end
                    for s = 1 : numel(ph_ref_gf(r).go_id)
                        go_id = ph_ref_gf(r).go_id(s);
                        az(isnan(az(:,go_id)) == 0, go_id) = ph_ref_gf(r).az(isnan(az(:,go_id)) == 0, s);
                        el(isnan(el(:,go_id)) == 0, go_id) = ph_ref_gf(r).el(isnan(el(:,go_id)) == 0, s);
                    end

                    [pierce_point(r).lat, pierce_point(r).lon, pierce_point(r).mf] = deal(nan(size(pr_ref_gf(r).az)));
                    [pierce_point(r).lat(:, pr_ref_gf(r).go_id), pierce_point(r).lon(:, pr_ref_gf(r).go_id), pierce_point(r).mf(:, pr_ref_gf(r).go_id)] = Atmosphere.getPiercePoint(lat / 180 * pi, lon / 180 * pi, h_ortho, zero2nan(az(:, pr_ref_gf(r).go_id) / 180 * pi), zero2nan(el(:, pr_ref_gf(r).go_id) / 180 * pi), Core_SEID.getThinShellHeight());
                end
                
                max_sat = 0;
                for r = 1 : numel(ref)
                    max_sat = max(max_sat, max(pr_ref_gf(r).go_id));
                end
                
                min_gap = 3;
                
                % Extract syncronized C4 L4 diff
                for t = 1 : numel(trg)
                    if any(any(id_sync{t}))
                        log.addMessage(log.indent(sprintf('Computing interpolated geometry free for target %d / %d', t, numel(trg))));
                                                
                        % Create a unique matrix containing all the geometry free of the ref receivers ---------------------------------------------------------
                        max_sat_trg = max(max_sat, max(trg(t).go_id));
                        ph_gf = nan(size(id_sync{t}, 1), max_sat_trg, numel(ref));
                        pr_gf = nan(size(id_sync{t}, 1), max_sat_trg, numel(ref));
                        sigma_ph_gf = nan(max_sat_trg, numel(ref));
                        sigma_pr_gf = nan(max_sat_trg, numel(ref));
                        for r = 1 : numel(ref)
                            id_ok = find(~isnan(id_sync{t}(:,r)));
                            id_ok_ref = id_sync{t}(id_ok,r);
                            ph_gf(id_ok, ph_ref_gf(r).go_id, r) = zero2nan(ph_ref_gf(r).obs(id_ok_ref, :));
                            % Import CS and outliers from receivers
                            for s = 1 : numel(ref(r).ph_idx)
                                if not(isempty(ref(r).sat.outliers_ph_by_ph))
                                    ph_gf(id_ok(find(ref(r).sat.outliers_ph_by_ph(id_ok_ref,s))), ref(r).go_id(ref(r).ph_idx(s)), r) = nan;
                                end
                                if not(isempty(ref(r).sat.outliers_ph_by_ph))
                                    ph_gf(id_ok(find(ref(r).sat.cycle_slip_ph_by_ph(id_ok_ref,s))), ref(r).go_id(ref(r).ph_idx(s)), r) = nan;
                                end
                                %
                                %                         % fill small gaps
                                %                         lim = getFlagsLimits(isnan(ph_gf(:, ref(r).go_id(s), r)));
                                %                         lim(lim(:,2) - lim(:,1) > min_gap,:) = [];
                                %                         idx = false(size(ph_gf, 1), 1);
                                %                         for l = 1 : size(lim, 1)
                                %                             idx(lim(l, 1) : lim(l, 2)) = true;
                                %                         end
                                %                         if sum(idx) > 0
                                %                             ph_gf(:,ref(r).go_id(s), r) = simpleFill1D(ph_gf(:, ref(r).go_id(s), r), idx);
                                %                         end
                            end
                            pr_gf(id_ok, pr_ref_gf(r).go_id, r) = zero2nan(pr_ref_gf(r).obs(id_ok_ref, :));
                            sigma_pr_gf(pr_ref_gf(r).go_id, r) = pr_ref_gf(r).sigma(:);
                            sigma_ph_gf(ph_ref_gf(r).go_id, r) = ph_ref_gf(r).sigma(:);
                        end
                        pr_gf_diff = diff(pr_gf);
                        ph_gf_diff = diff(ph_gf);
                        
                        % Extract L1 from target ---------------------------------------------------------------------------------------------------------------
                        % ph_gf pr_gf ph_gf_diff have max_sat satellite data stored
                        % pierce_point(r).lat could have different size receiver by receiver
                        % indexes conversion is id = ph_ref_gf(r).go_id
                        
                        ph1_trg = [];
                        id_ph = [];
                        for sys = systems
                            frqs = num2str(trg(t).getFreqs(sys));
                            if ~isempty(frqs)
                                [ph1_t, id_ph_t] = trg(t).getObs(['L' frqs(1) '?'],sys);
                                ph1_trg = [ph1_trg; ph1_t];
                                id_ph = [id_ph; id_ph_t];
                            end
                        end
                        [lat, lon, ~, h_ortho] = trg(t).getMedianPosGeodetic;
                        ph1_goid = trg(t).go_id(id_ph)';
                        trg_go_id = unique(ph1_goid);
                        [lat_pp, lon_pp, iono_mf] = Atmosphere.getPiercePoint(lat / 180 * pi, lon / 180 * pi, h_ortho, trg(t).sat.az(:, trg_go_id) / 180 * pi, zero2nan(trg(t).sat.el(:, trg_go_id) / 180 * pi), Core_SEID.getThinShellHeight());
                        
                        % Create a unique matrix containing all the pierce_points/mf of the ref receivers (sat by sat) ------- AND INTERPOLATE!!!!! ------------
                        % It is necessary to better sync satellites in view
                        % this part of the code needs to be improved
                        max_n_sat = max(trg(t).go_id);
                        if strcmp(method, 'ls')
                            trg_gf = nan(trg(t).time.length, max_n_sat);
                        else
                            trg_pr_gf = nan(trg(t).time.length, max_n_sat);
                            trg_ph_gf = nan(trg(t).time.length, max_n_sat);
                        end
                        wl1 = zeros(max_n_sat,1);
                        wl2 = zeros(max_n_sat,1);
                        band = char(zeros(max_n_sat,1));
                        
                        for s = 1 : numel(trg_go_id)
                            mf_sat = nan(size(id_sync{t},1), numel(ref));
                            lat_sat = nan(size(id_sync{t},1), numel(ref));
                            lon_sat = nan(size(id_sync{t},1), numel(ref));
                            for r = 1 : numel(ref)
                                id_sat = trg_go_id(s);
                                if id_sat <= size(pierce_point(r).lat,2)
                                    id_ok = (~isnan(id_sync{t}(:,r)));
                                    id_ok_ref = id_sync{t}(id_ok,r);
                                    mf_sat(id_ok, r) = pierce_point(r).mf(id_ok_ref, id_sat);
                                    lat_sat(id_ok, r) = pierce_point(r).lat(id_ok_ref, id_sat);
                                    lon_sat(id_ok, r) = pierce_point(r).lon(id_ok_ref, id_sat);
                                end
                            end
                            
                            sys = Core.getConstellationCollector.getSysPrn(trg_go_id(s));
                            frqs = num2str(trg(t).getFreqs(sys));
                            l_vec = Core.getConstellationCollector.getWavelength(trg_go_id(s));                            
                            wl1(trg_go_id(s)) = l_vec(Core.getConstellationCollector.getSys(sys).CODE_RIN3_2BAND == frqs(1)); % <- we have that
                            wl2(trg_go_id(s)) = l_vec(find(Core.getConstellationCollector.getSys(sys).CODE_RIN3_2BAND ~= frqs(1),1,'first'));  % <- we want to generate that
                            band(trg_go_id(s)) = Core.getConstellationCollector.getSys(sys).CODE_RIN3_2BAND(find(Core.getConstellationCollector.getSys(sys).CODE_RIN3_2BAND ~= frqs(1),1,'first'));
                            iono_factor = wl2(trg_go_id(s))^2 - wl1(trg_go_id(s))^2;
                            
                            % Apply MF (of the ref)
                            if flag_mf
                                for r = 1 : numel(ref)
                                    pr_gf_diff(:,trg_go_id(s), r) = pr_gf_diff(:,trg_go_id(s), r) ./ mf_sat(2:end, r);
                                    ph_gf_diff(:,trg_go_id(s), r) = ph_gf_diff(:,trg_go_id(s), r) ./ mf_sat(2:end, r);
                                end
                            end
                            
                            % CALL INTERPOLATION!!!
                            
                            if strcmp(method, 'ls')
                                % Least Squares approach
                                data_pr = squeeze(pr_gf(:,trg_go_id(s),:));
                                data_ph = squeeze(ph_gf(:,trg_go_id(s),:));
                                sigma_pr = sigma_pr_gf(trg_go_id(s),:);
                                sigma_ph = sigma_ph_gf(trg_go_id(s),:);
                                trg_gf(id_sync{t}(:,       t + numel(ref)), trg_go_id(s)) = iono_factor * Core_SEID.satDataInterpLS(lat_sat, lon_sat, data_pr, data_ph, sigma_pr, sigma_ph, lat_pp(id_sync{t}(:, t + numel(ref)), s), lon_pp(id_sync{t}(:, t + numel(ref)), s));
                            else
                                % SEID approach based on derivate
                                trg_pr_gf(id_sync{t}(2 : end, t + numel(ref)), trg_go_id(s)) = iono_factor * Core_SEID.satDataInterp(lat_sat(2 : end, :), lon_sat(2 : end, :), squeeze(pr_gf_diff(:,trg_go_id(s),:)), lat_pp(id_sync{t}(2 : end,t + numel(ref)), s), lon_pp(id_sync{t}(2 : end, t + numel(ref)), s), method);
                                trg_ph_gf(id_sync{t}(2 : end, t + numel(ref)), trg_go_id(s)) = iono_factor * Core_SEID.satDataInterp(lat_sat(2 : end, :), lon_sat(2 : end, :), squeeze(ph_gf_diff(:,trg_go_id(s),:)), lat_pp(id_sync{t}(2 : end,t + numel(ref)), s), lon_pp(id_sync{t}(2 : end, t + numel(ref)), s), method);
                            end
                            
                            % Apply (inverse) MF of the (trg)
                            if flag_mf
                                if strcmp(method, 'ls')
                                    trg_gf(:, trg_go_id(s)) = trg_gf(:, trg_go_id(s)) .* iono_mf(:, s);
                                else
                                    trg_pr_gf(:, trg_go_id(s)) = trg_pr_gf(:, trg_go_id(s)) .* iono_mf(:, s);
                                    trg_ph_gf(:, trg_go_id(s)) = trg_ph_gf(:, trg_go_id(s)) .* iono_mf(:, s);
                                end
                            end
                        end
                        
                        % Interpolate the diff (derivate) of L4, now rebuild L4 by cumsum (integral)
                        
                        if strcmp(method, 'ls')
                            id_ko = movstd(Core_Utils.diffAndPred(trg_gf(:,s), 2), 10) > 0.01;
                            trg_gf(id_ko) = nan;
                            trg_ph_gf = trg_gf;
                            trg_pr_gf = trg_gf;
                        else                            
                            trg_ph_gf(abs(trg_ph_gf) > 0.5) = nan; % remove outliers
                            
                            %% experimental (fill gaps)
                            %                     inan = isnan(trg_ph_gf);
                            %                     for i = trg(t).go_id(id_ph)
                            %                         trg_ph_gf(:,i) = simpleFill1D(trg_ph_gf(:,i),isnan(trg_ph_gf(:,i)),'linear');
                            %                     end
                            %                     trg_ph_gf(abs(trg_ph_gf) > 0.5) = nan; % remove outliers
                            
                            % aesthetic approach
                            corr_diff = trg_ph_gf;
                            for s = 1 : size(corr_diff, 2)
                                lim = getFlagsLimits(~isnan(corr_diff(:,s)));
                                for l = 1 : size(lim, 1)
                                    id_arc = lim(l,1) : lim(l,2);
                                    trg_ph_gf(id_arc, s) = cumsum(trg_ph_gf(id_arc, s));
                                end
                            end
                            
                            corr_diff = trg_pr_gf;
                            for s = 1 : size(corr_diff, 2)
                                lim = getFlagsLimits(~isnan(corr_diff(:,s)));
                                for l = 1 : size(lim, 1)
                                    id_arc = lim(l,1) : lim(l,2);
                                    trg_pr_gf(id_arc, s) = cumsum(trg_pr_gf(id_arc, s));
                                end
                            end
                            % cumsum approach (faster but ugly)
                            %inan = isnan(trg_ph_gf);
                            %trg_ph_gf = cumsum(nan2zero(trg_ph_gf));
                            %trg_ph_gf(inan) = nan;
                        end
                        %for s = 1:32
                        %    figure(100+s); clf; dockAllFigures; plot(movstd(Core_Utils.diffAndPred(trg_ph_gf(:,s), 2),10)); hold on; plot(movstd(Core_Utils.diffAndPred(trg_gf(:,s), 2),10));
                        %end
                        %                     wl1 = trg(t).state.getConstellationCollector().gps.L_VEC(1);
                        %                     wl2 = trg(t).state.getConstellationCollector().gps.L_VEC(2);
                        ph2 = nan(size(ph1_trg));
                        % L1 * wl1 - L2 * wl2 = gf
                        % L2 = (L1 * wl1 - gf) / wl2;
                        ph2 = (ph1_trg .* repmat(wl1(trg(t).go_id(id_ph)),1,size(ph1_trg,2)) - trg_ph_gf(:, trg(t).go_id(id_ph))') ./ repmat(wl2(trg(t).go_id(id_ph)),1,size(ph1_trg,2));
                        
                        [~, ~, ~, flag] = trg(t).getBestCodeObs();
                        pr1 = [];
                        id_pr = [];
                        for sys = systems
                            frqs = num2str(trg(t).getFreqs(sys));
                            if ~isempty(frqs)
                                [pr1_t, id_pr_t] = trg(t).getObs(['C' frqs(1) '?'],sys);
                                pr1 = [pr1; pr1_t];
                                id_pr = [id_pr; id_pr_t];
                            end
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
                        for sys = systems
                            bnd = unique(band(Core.getConstellationCollector.getSysPrn(1:max_n_sat) == sys));
                            bnd(bnd == 0) = [];
                            [~, id_ph] = trg(t).getObs(['L' bnd],sys);
                            if ~isempty(id_ph)
                                log.addMessage(log.indent(sprintf('Removing %sL%s observations already present in the target receiver %d / %d',sys,bnd, t, numel(trg))));
                                trg(t).remObs(id_ph);
                            end
                        end
                        
                        for sys = systems
                            bnd = unique(band(Core.getConstellationCollector.getSysPrn(1:max_n_sat) == sys));
                            bnd(bnd == 0) = [];
                            [~, id_pr] = trg(t).getObs(['C' bnd],sys);
                            if ~isempty(id_pr)
                                log.addMessage(log.indent(sprintf('Removing %sC%s observations already present in the target receiver %d / %d',sys,bnd, t, numel(trg))));
                                trg(t).remObs(id_pr);
                            end
                        end
                        
                        % Inject the new synthesised phase
                        sys_o = Core.getConstellationCollector.getSysPrn(trg_go_id);
                        for sys = unique(sys_o)
                            sys_pr = Core.getConstellationCollector.getSysPrn(pr1_goid);
                            sys_ph = Core.getConstellationCollector.getSysPrn(ph1_goid);
                            
                            
                            lid_pr = sys_pr == sys;
                            lid_ph = sys_ph == sys;
                            
                            bnd = unique(band(pr1_goid(lid_pr)));
                             bnd(bnd == 0) = [];
                            log.addMessage(log.indent(sprintf('Injecting SEID L%s into target receiver %d / %d', bnd,t, numel(trg))));
                            trg(t).injectObs(nan2zero(pr2(lid_pr,:)), wl2(pr1_goid(lid_pr))', str2num(bnd),[ 'C' bnd 'F'], pr1_goid(lid_pr));
                            trg(t).injectObs(nan2zero(ph2(lid_ph,:)), wl2(ph1_goid(lid_ph))', str2num(bnd),[ 'L' bnd 'F'], ph1_goid(lid_ph));
                        end
                        %trg(t).injectObs(nan2zero(ref(1).getObs('C2')), wl2, 2, 'C2 ', trg_go_id);
                        %trg(t).injectObs(nan2zero(ref(1).getObs('L2')), wl2, 2, 'L2 ', trg_go_id);
                        
                        trg(t).keepEpochs(id_sync{t}(:,t + numel(ref)));
                        trg(t).updateDetectOutlierMarkCycleSlip();
                    end
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
            
            spline_smooth_time = 900;
            ref = ref(~ref.isEmpty_mr);
            trg = trg(~trg.isEmpty_mr);
            if ~isempty(trg)
                rec(1:numel(ref)) = ref;
                rec(numel(ref) + (1 : numel(trg))) = trg;
                obs_type(1:numel(ref)) = 2;
                obs_type(numel(ref) + (1 : numel(trg))) = 0;
                %%% check if one of the target receiver does have too few
                %%% epoch
                [rec, ref, obs_type] = Core_SEID.removeRecWithFewObs(ref,rec,obs_type);
                
                [p_time, id_sync] = Receiver_Commons.getSyncTimeTR(rec, obs_type);
                log = Core.getLogger();
                
                log.addMarkedMessage('Starting REMIONO processing')
                log.addMessage(log.indent('Getting Geometry free from reference receivers'));
                
                sys_c = 'G';
                for r = 1 : numel(ref)
                    % combine code and phase
                    ph_ref_gf(r) = ref(r).getPrefGeometryFree('L',sys_c);
                    pr_ref_gf = ref(r).getPrefGeometryFree('C',sys_c);
                    
                    % sync pr with ph
                    [~, id_ko] = setdiff(pr_ref_gf.go_id, ph_ref_gf(r).go_id);
                    if ~isempty(id_ko)
                        % remove pseudoranges not present as phases 
                        pr_ref_gf.obs(:,id_ko) = [];
                        pr_ref_gf.obs_code(id_ko, :) = [];
                        pr_ref_gf.wl(id_ko) = [];
                        pr_ref_gf.el(:,id_ko) = [];
                        pr_ref_gf.az(:,id_ko) = [];
                        pr_ref_gf.prn(id_ko) = [];
                        if ~isempty(pr_ref_gf.snr)
                            pr_ref_gf.snr(:,id_ko) = [];
                        end
                        if ~isempty(pr_ref_gf.cycle_slip)
                            pr_ref_gf.cycle_slip(:,id_ko) = [];
                        end
                        pr_ref_gf.go_id(id_ko) = [];
                        pr_ref_gf.sigma(id_ko) = [];
                    end
                    
                    idx_nan = ph_ref_gf(r).obs == 0;

                    el = ph_ref_gf(r).el / 180 * pi;
                    az = ph_ref_gf(r).az / 180 * pi;
                    [ph_ref_gf(r).obs] = ref(r).ionoCodePhaseSmt(zero2nan(pr_ref_gf.obs), pr_ref_gf.sigma.^2, zero2nan(ph_ref_gf(r).obs), ph_ref_gf(r).sigma.^2, ph_ref_gf(r).getAmbIdx(), 0.01, el);
                    
                    % Apply mapping function
                    [lat, lon, ~, h_ortho] = ref(r).getMedianPosGeodetic;
                    [pierce_point(r).lat, pierce_point(r).lon, pierce_point(r).mf] = Atmosphere.getPiercePoint(lat / 180 * pi, lon / 180 * pi, h_ortho, az, zero2nan(el), Core_SEID.getThinShellHeight());
                    ph_ref_gf(r).obs(idx_nan) = nan;
                    %iono_ref(r).obs = ref(r).applyMF(iono_ref(r).obs, iono_ref(r).mf, 1);
                    %iono_ref(r).obs(idx_nan) = nan;

                    %iono_ref(r).obs = ref(r).smoothSatData([], [], zero2nan(iono_ref(r).obs), false(size(iono_ref(r).cycle_slip)), [], 300 / iono_ref(r).time.getRate); % <== supposing no more cycle slips
                end
                
                max_sat = 0;
                for r = 1 : numel(ref)
                    max_sat = max(max_sat, max(ph_ref_gf(r).go_id));
                end
                                
                % Extract syncronized C4 L4 diff
                for t = 1 : numel(trg)
                    if any(any(id_sync{t}))
                        % trimming the target receiver to mach the id_sync of the reference stations
                        log.addMessage(log.indent(sprintf('Keeping only the epochs in common between "%s" and the reference stations', trg(t).parent.getMarkerName4Ch) ));
                        trg(t).keepEpochs(id_sync{t}(:,numel(ref) + t));
                        id_sync{t}(:,numel(ref) + t) = 1 : size(id_sync{t}, 1);
                        
                        log.addMessage(log.indent(sprintf('Computing interpolated geometry free for target %d / %d', t, numel(trg))));
                        
                        iono_sync = nan(size(id_sync{t}, 1), max_sat, numel(ref));
                        for r = 1 : numel(ref)
                            id_ok = find(~isnan(id_sync{t}(:,r)));
                            id_ok_ref = id_sync{t}(id_ok,r);
                            iono_sync(:, ph_ref_gf(r).go_id, r) = zero2nan(ph_ref_gf(r).obs(id_sync{t}(:,r), :));% - mean(iono_ref(r).obs(id_sync{t}(:,r), :), 1, 'omitnan');
                            for s = 1 : numel(ref(r).ph_idx)
                                iono_sync(id_ok(find(ref(r).sat.outliers_ph_by_ph(id_ok_ref,s))), ref(r).go_id(ref(r).ph_idx(s)), r) = nan;
                                iono_sync(id_ok(find(ref(r).sat.cycle_slip_ph_by_ph(id_ok_ref,s))), ref(r).go_id(ref(r).ph_idx(s)), r) = nan;
                                %
                                %                         % fill small gaps
                                %                         lim = getFlagsLimits(isnan(ph_gf(:, ref(r).go_id(s), r)));
                                %                         lim(lim(:,2) - lim(:,1) > min_gap,:) = [];
                                %                         idx = false(size(ph_gf, 1), 1);
                                %                         for l = 1 : size(lim, 1)
                                %                             idx(lim(l, 1) : lim(l, 2)) = true;
                                %                         end
                                %                         if sum(idx) > 0
                                %                             ph_gf(:,ref(r).go_id(s), r) = simpleFill1D(ph_gf(:, ref(r).go_id(s), r), idx);
                                %                         end
                            end
                        end
                        
                        % DIFF:
                        % Simple approach ---------------------------------------------
                        % iono_diff = diff(nan2zero(iono_sync)); % deprecate ("loses one epoch")
                        % Complex approach --------------------------------------------
                        % Compute median iono (to reduce the signal)
                        iono_diff = nan(size(iono_sync));
                        for r = 1: numel(ph_ref_gf)
                            tmp = Core_Utils.diffAndPred(zero2nan(iono_sync(:,:,r)));
                            mm_tmp = movmedian(tmp, 3);
                            id_ko = abs(tmp - mm_tmp) > 0.03;
                            tmp(id_ko) = mm_tmp(id_ko);
                            iono_diff(:,:,r) = tmp;
                        end
                        
                        % Reduce the iono signal by a smoothed median iono diff
                        med_iono_diff = median(iono_diff, 3, 'omitnan');
                        med_iono_diff = Receiver_Commons.smoothSatData([], [], zero2nan(med_iono_diff), false(size(med_iono_diff)), [], spline_smooth_time / p_time(t).getRate, 6); % <== supposing no more cycle slips
                        for r = 1: numel(ph_ref_gf)
                            d_tmp = iono_diff(:,:,r) - med_iono_diff;
                            d_tmp(abs(d_tmp) > 0.015) = nan;
                            iono_diff(:,:,r) = d_tmp;
                        end
                        % -------------------------------------------------------------
                        
                        % ph_gf pr_gf ph_gf_diff have max_sat satellite data stored
                        % pierce_point(r).lat could have different size receiver by receiver
                        % indexes conversion is id = ph_ref_gf(r).go_id
                        
                        [ph1, id_ph] = trg(t).getObs('L1','G');
                        [lat, lon, ~, h_ortho] = trg(t).getMedianPosGeodetic;
                        ph1_goid = trg(t).go_id(id_ph)';
                        trg_go_id = unique(ph1_goid);
                        [lat_pp, lon_pp, iono_mf] = Atmosphere.getPiercePoint(lat / 180 * pi, lon / 180 * pi, h_ortho, trg(t).sat.az(:, trg_go_id) / 180 * pi, zero2nan(trg(t).sat.el(:, trg_go_id) / 180 * pi), 350*1e3);
                        
                        % It is necessary to better sync satellites in view
                        iono_trg = nan(trg(t).time.length, max_n_sat);
                        for s = 1 : numel(trg_go_id)
                            lat_sat = nan(size(id_sync{t},1), numel(ref));
                            lon_sat = nan(size(id_sync{t},1), numel(ref));
                            for r = 1 : numel(ref)
                                id_sat = unique(ph_ref_gf(r).go_id) == trg_go_id(s);
                                if sum(id_sat) == 1
                                    lat_sat(:, r) = pierce_point(r).lat(id_sync{t}(:,r), id_sat);
                                    lon_sat(:, r) = pierce_point(r).lon(id_sync{t}(:,r), id_sat);
                                end
                            end
                            % DIFF:
                            tmp = Core_SEID.satDataInterp(lat_sat(:, :), lon_sat(:, :), squeeze(iono_diff(:,trg_go_id(s),:)),  lat_pp(id_sync{t}(:,t + numel(ref)), s), lon_pp(id_sync{t}(:,t + numel(ref)), s));
                            tmp(abs(tmp) > 0.1) = nan; % remove outliers
                            iono_trg(id_sync{t}(:, t + numel(ref)), trg_go_id(s)) = tmp;
                        end
                        iono_trg = Receiver_Commons.smoothSatData([], [], zero2nan(iono_trg), false(size(iono_trg)), [], spline_smooth_time / p_time(t).getRate, 6); % <== supposing no more cycle slips
                        iono_trg(id_sync{t}(:, t + numel(ref)), :) = iono_trg(id_sync{t}(:, t + numel(ref)), :) + med_iono_diff;
                        
                        % Interpolate the diff (derivate) of L4, now rebuild L4 by cumsum (integral)
                        inan = isnan(iono_trg);
                        iono_trg = cumsum(nan2zero(iono_trg));
                        iono_trg(inan) = nan;
                        
                        %% SEID approach => Synth L2
                        seid_approach = false;
                        if seid_approach
                            wl1 = trg(t).state.getConstellationCollector().gps.L_VEC(1);
                            wl2 = trg(t).state.getConstellationCollector().gps.L_VEC(2);
                            ph2 = nan(size(ph1));
                            % L1 * wl1 - L2 * wl2 = gf
                            % L2 = (L1 * wl1 - gf) / wl2;
                            ph2 = (ph1 * wl1 - iono_trg(:, trg(t).go_id(id_ph))') / wl2;
                            
                            % Also remove Iono
                            ph1 = ph1 - iono_trg(:, trg(t).go_id(id_ph))' * (wl1 / (wl1^2 - wl2^2));
                            trg(t).obs(id_ph, :) = nan2zero(ph1);
                            ph2 = ph2 - iono_trg(:, trg(t).go_id(id_ph))' * (wl2 / (wl1^2 - wl2^2));
                            
                            [~, ~, ~, flag] = trg(t).getBestCodeObs();
                            [pr1, id_pr] = trg(t).getObs(flag(1,1:3),'G');
                            pr1_goid = trg(t).go_id(id_pr);
                            % C2 - C1 = gf
                            % C2 = C1 + gf
                            pr2 = pr1 + iono_trg(:, trg(t).go_id(id_pr))';
                            
                            % Also remove Iono
                            pr1 = pr1 + iono_trg(:, trg(t).go_id(id_pr))' * (wl1^2 / (wl1^2 - wl2^2));
                            trg(t).obs(id_pr, :) = nan2zero(pr1);
                            pr2 = pr2 + iono_trg(:, trg(t).go_id(id_pr))' * (wl2^2 / (wl1^2 - wl2^2));
                            
                            % Remove the all the other frequency stored into the receiver
                            % Bacause now L1 are inconsistent with the other bands
                            % A different approach could be to correct all the other bands
                            for f = 2:8
                                band = num2str(f);
                                [~, id_ph] = trg(t).getObs(['L' band],'G');
                                if ~isempty(id_ph)
                                    log.addMessage(log.indent(sprintf('Removing L%d observations already present in the target receiver %d / %d', f, t, numel(trg))));
                                    trg(t).remObs(id_ph);
                                end
                                [~, id_pr] = trg(t).getObs(['C' band],'G');
                                if ~isempty(id_pr)
                                    log.addMessage(log.indent(sprintf('Removing C%d observations already present in the target receiver %d / %d', f, t, numel(trg))));
                                    trg(t).remObs(id_pr);
                                end
                            end
                            
                            % Inject the new synthesised phase
                            log.addMessage(log.indent(sprintf('Injecting SEID L2 into target receiver %d / %d', t, numel(trg))));
                            trg(t).injectObs(nan2zero(pr2), wl2, 2, 'C2F', pr1_goid);
                            trg(t).injectObs(nan2zero(ph2), wl2, 2, 'L2F', ph1_goid);
                        else %% REMIONO approach => rem I
                            trg(t).setDefaultRIE('rem_iono');
                            
                            wl1 = trg(t).state.getConstellationCollector().gps.L_VEC(1);
                            wl2 = trg(t).state.getConstellationCollector().gps.L_VEC(2);
                            
                            alpha = wl1^2 / (wl1^2 - wl2^2);
                            
                            % Correct phase 1 for iono delay
                            [~, ~, ~, flag] = trg(t).getBestCodeObs();
                            [pr1, id_pr] = trg(t).getObs(flag(1,1:3),'G');
                            pr1_goid = trg(t).go_id(id_pr);
                            pr1 = pr1 + iono_trg(:, trg(t).go_id(id_pr))' * (wl1^2 / (wl1^2 - wl2^2));
                            %pr2 = pr2 + iono_trg(:, trg(t).go_id(id_pr))' * (wl2^2 / (wl1^2 - wl2^2));
                            %ph1 = (ph1 * wl1 - iono_trg(:, trg(t).go_id(id_ph))' * (wl1^2 / (wl1^2 - wl2^2))) / wl1;
                            ph1 = ph1 - iono_trg(:, trg(t).go_id(id_ph))' * (wl1 / (wl1^2 - wl2^2));
                            %ph2 = ph2 - iono_trg(:, trg(t).go_id(id_ph))' * (wl2 / (wl1^2 - wl2^2));
                            
                            % Inject the new synthesised phase
                            log.addMessage(log.indent(sprintf('Injecting iono reduced observations into target receiver %d / %d', t, numel(trg))));
                            trg(t).obs(id_pr, :) = nan2zero(pr1);
                            trg(t).obs(id_ph, :) = nan2zero(ph1);
                            
                            % Remove the all the other frequency stored into the receiver
                            % Bacause now L1 are inconsistent with the other bands
                            % A different approach could be to correct all the other bands
                            for f = 2:8
                                band = num2str(f);
                                [~, id_ph] = trg(t).getObs(['L' band],'G');
                                if ~isempty(id_ph)
                                    log.addMessage(log.indent(sprintf('Removing L%d observations already present in the target receiver %d / %d', f, t, numel(trg))));
                                    trg(t).remObs(id_ph);
                                end
                                [~, id_pr] = trg(t).getObs(['C' band],'G');
                                if ~isempty(id_pr)
                                    log.addMessage(log.indent(sprintf('Removing C%d observations already present in the target receiver %d / %d', f, t, numel(trg))));
                                    trg(t).remObs(id_pr);
                                end
                            end
                        end
                        trg(t).keepEpochs(id_sync{t}(:,t + numel(ref)));
                        trg(t).updateDetectOutlierMarkCycleSlip();
                    end
                end
                log.addMarkedMessage('Syncing times, computing reference time');
            end            
        end
        
        
        function [rec, ref, obs_type] = removeRecWithFewObs(ref,rec,obs_type)
            % remove refrence with few observations
            %
            % SYNTAX:
            %   [rec,ref] = removeRecWithFewObs(ref,rec)
            log = Core.getLogger();
            n_ep = zeros(numel(ref),1);
            for i = 1 : numel(ref)
                n_ep(i) = rec(i).time.length * rec(i).time.getRate;
            end
            idx_poss_rem = n_ep < 0.7*median(n_ep);
            if sum(~idx_poss_rem) > 3
                idx_poss_rem = find(idx_poss_rem);
                
                if any(idx_poss_rem)
                    recs_str = [];
                    for j = 1 : length(idx_poss_rem)
                        recs_str = [recs_str ' ' rec(idx_poss_rem(j)).parent.getMarkerName4Ch];
                        log.addWarning(['Receveivers' recs_str ' has to many empty epoch, the will not be used as reference']);
                    end
                end
                rec(idx_poss_rem) = [];
                ref(idx_poss_rem) = [];
                obs_type(idx_poss_rem) = [];
                
            end
        end
    end

    methods (Static)
        function data_q = satDataInterp(lat_in, lon_in, data_in, lat_q, lon_q, method)
            % Interpolate the data given for n points "_in" to query point "_q"
            %
            %
            % INPUT:
            %   lat_in     [n_epoch x n_ref]  latitude [deg]
            %   lon_in     [n_epoch x n_ref]  longitude [deg]
            %   data_in    [n_epoch x n_ref]  data [deg]
            %   lat_q      [n_epoch x 1]      query latitude
            %   lon_q      [n_epoch x 1]      query longitude
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
                             % maybe a 3D Kriging is better
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
        
        function data_q = satDataInterpLS(lat_in, lon_in, data_pr_in, data_ph_in, sigma_pr_in, sigma_ph_in, lat_q, lon_q)
            % Interpolate the data given for n points "_in" to query point "_q"
            %
            % NOTE:
            %   For methods ('plane' and 'distance') one between data_pr_in and data_ph_in 
            %   should be empty, they work one data type at a time
            %   For method 'mq' both the set can be passed to the interpolator
            %  
            %
            % INPUT:
            %   lat_in     [n_epoch x n_ref]  latitude [deg]
            %   lon_in     [n_epoch x n_ref]  longitude [deg]
            %   data_pr_in [n_epoch x n_ref]  pseudo-range data [deg]
            %   data_ph_in [n_epoch x n_ref]  carrier phase data [deg]
            %   lat_q      [n_epoch x 1]      query latitude
            %   lon_q      [n_epoch x 1]      query longitude
            %
            % SYNTAX:
            %   data_q = Core_SEID.satDataInterp(lat_in, lon_in, data_in, lat_q, lon_q);
            
            %% prepare input data and sigma
            data_pr = data_pr_in;
            data_ph = data_ph_in;
            sigma_pr = sigma_pr_in;
            sigma_ph = sigma_ph_in;
            if isempty(data_ph_in)
               sigma_ph = sigma_pr * nan;
               data_ph = data_pr_in * nan;
            end
            if isempty(data_pr_in)
               sigma_pr = sigma_ph * nan;
               data_pr = data_ph_in * nan;
            end
            
            n_ref = size(data_ph,2);
            
            % prepare weights
            w = zeros(size(data_pr,1), n_ref);
            id_ok = find(any(data_ph, 2) | any(data_pr, 2));
            for i = id_ok'
                d = sphericalDistance(lat_in(i, :) ./ pi * 180, lon_in(i, :) ./ pi * 180, lat_q(i) ./ pi * 180, lon_q(i) ./ pi * 180);
                w(i,:) = (1 ./ (d + eps))';
                w(i,:) = w(i,:) ./ sum(nan2zero(w(i,:)));
            end
            
            % filter data with no lat/lon pierce point
            data_ph(isnan(zero2nan(w))) = nan;
            % do not estimate anything it without phase
            data_pr(isnan(data_ph)) = nan; 
            
            % dimensions of data_pr and data_ph should be the same
            n_obs = sum(~isnan(data_pr(:))) + sum(~isnan(data_ph(:)));
            n_bia = size(data_pr, 2);
            for r = 1 : size(data_ph, 2)
                lim = getFlagsLimits(~isnan(data_ph(:,r)));
                n_bia = n_bia + size(lim, 1);
            end
            
            id_ok = find(any(data_ph, 2) | any(data_pr, 2));
            offset = cumsum(any(data_ph, 2) | any(data_pr, 2));
            
            % prepare design matrix
            A = sparse([], [], [], n_obs, n_bia + numel(id_ok), n_obs * n_bia * 2);
            n = 0;
            id_bia = size(data_pr, 2);
            y0 = zeros(n_obs,1);
            sigma = zeros(n_obs,1);
            low_cost_effects = 5; % Consider pseudoranges 5 times worse than expected
            for r = 1 : size(data_pr, 2)
                % pseudo-ranges - one bias per receiver
                id_r = offset(~isnan(data_pr(:,r)));
                if numel(id_r) > 0
                    n = n(end) + (1 : numel(id_r));
                    A(n, r) = 1;
                    A(n(:) + (n_bia + id_r(:) - 1) * n_obs) = 1;
                    y0(n) = data_pr(~isnan(data_pr(:,r)), r);
                    sigma(n) = (low_cost_effects * sigma_pr(r)) ./ w(~isnan(data_pr(:,r)), r);
                end
                
                % phases - one bias per arc
                lim = getFlagsLimits(~isnan(data_ph(:,r)));                
                for l = 1 : size(lim, 1)
                    id_bia = id_bia + 1;
                    id_r = (1 + 0:(lim(l,2) - lim(l,1) + 1))' + offset(lim(l,1)) - 1;
                    if numel(id_r) > 0
                        n = n(end) + (1 : numel(id_r));
                        A(n, id_bia) = 1;
                        A(n(:) + (n_bia + id_r(:) - 1) * n_obs) = 1;
                        y0(n) = data_ph(lim(l,1) : lim(l,2), r);
                        sigma(n) = sigma_ph(r) ./ w(lim(l,1) : lim(l,2), r);
                    end
                end
            end
            invQ = spdiags(1 ./ (sigma.^2), 0, numel(sigma), numel(sigma));
            
            % System solve
            N = (A' * invQ * A);
            b = ((A' * invQ) * y0);
            
            % Lagrange multiplier: set the mean to the mean of the gf of code
            N = [N sparse([ones(n_ref, 1); zeros(size(N, 1) - n_ref, 1)])];
            N = [N; sparse([ones(1, n_ref) zeros(1, size(N, 2) - n_ref)])];
            b = [b; 0];
            
            warning off;
            x = N\b;
            warning on;
            
            % gettin the result
            data_q = nan(size(data_pr,1), 1);
            data_q(id_ok) = x((n_bia +1) : end - 1);
            
            data_q = zero2nan(data_q);
            % figure(100); imagesc(A)
            % fnum = 103;
            % figure(fnum); clf; dockAllFigures
            % plot(data_pr,'.-'); hold on;
            % plot(data_ph,'.-'); hold on;
            % plot(zero2nan(data_q), 'm', 'LineWidth', 3);
        end
        
    end

end
