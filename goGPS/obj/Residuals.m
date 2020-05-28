classdef Residuals < Exportable
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 1.0b7
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2020 Andrea Gatti, Giulio Tagliaferro, Eugenio Realini
    %  Written by:       Andrea Gatti
    %  Contributors:     ...
    %  A list of all the historical goGPS contributors is in CREDITS.nfo
    %--------------------------------------------------------------------------
    %
    %   This program is free software: you can redistribute it and/or modify
    %   it under the terms of the GNU General Public License as published by
    %   the Free Software Foundation, either version 3 of the License, or
    %   (at your option) any later version.
    %
    %   This program is distributed in the hope that it will be useful,
    %   but WITHOUT ANY WARRANTY; without even the implied warranty of
    %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %   GNU General Public License for more details.
    %
    %   You should have received a copy of the GNU General Public License
    %   along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %--------------------------------------------------------------------------
    % 01100111 01101111 01000111 01010000 01010011
    %--------------------------------------------------------------------------
    
    %% CONSTANTS
    properties (Constant)
        RES_TYPE = {'0: no residuals', '1: PREPRO', '2: U1 engine', '3: U2 engine'};
    end
    
    properties
        type        % 0,1,2,3 see RES_TYPE
        time        % time as GPS_Time                    GPS_Time [1 x 1] stores n_epoch
        value       % matrix of residuals
        
        prn         % prn of the satellite (1 x col of pr/ph)
        obs_code    % type of tracking of the column (e.g. GL1C, GL1CL2WI, ...)
        
        rec_coo     % <optional> Coordinates of the receiver
    end
    
    methods
        % Creator
        function this = Residuals()
            % Object creator
            this.reset();
        end
        
    end
    
    % =========================================================================
    %%  METHODS OBJ MANAGEMENT
    % =========================================================================
    
    methods % Public Access
        function is_empty = isEmpty(this)
            % Return the empty status
            %
            % SYNTAX
            %   is_empty = this.isEmpty()
            
            is_empty = this.type == 0 || this.time.isEmpty;
        end
        
        function reset(this)
            % Reset the stored residuals
            % Empty the object
            %
            % SYNTAX
            %   this.reset
            
            this.type = 0;
            this.time = GPS_Time();
            this.value = [];
            
            this.prn = [];
            this.obs_code = '';
            this.rec_coo = Coordinates();
        end
        
        function import(this, type, time, value, prn, obs_code, rec_coo)
            % Import new residuals (and delete the previous content)
            %
            % INPUT
            %   type        % 0,1,2,3 see RES_TYPE
            %   time        % time as GPS_Time                        GPS_Time [1 x 1] stores n_epoch
            %   value       % matrix of residuals
            %
            %   prn         % prn of the satellite (1 x col of pr/ph)
            %   obs_code    % type of tracking of the column (e.g. GL1C, GL1CL2WI, ...)
            %
            %   rec_coo     % <optional> Coordinates of the receiver
            %
            % SYNTAX
            %   this.import(type, time, value, prn, obs_code, rec_coo)
            this.init(type, time, value, prn, obs_code, rec_coo)
        end
        
        function append(this, type, time, value, prn, obs_code, rec_coo)
            % Append new residuals to the one already stored
            %
            % INPUT
            %   type        % 0,1,2,3 see RES_TYPE
            %   time        % time as GPS_Time                        GPS_Time [1 x 1] stores n_epoch
            %   value       % matrix of residuals
            %
            %   prn         % prn of the satellite (1 x col of pr/ph)
            %   obs_code    % type of tracking of the column (e.g. GL1C, GL1CL2WI, ...)
            %
            %   rec_coo     % <optional> Coordinates of the receiver
            %
            % SYNTAX
            %   this.import(type, time, value, prn, obs_code, rec_coo)
            
            res = Residuals();
            res.init(type, time, value, prn, obs_code, rec_coo);
            
            this.injest(res);
        end
        
        function injest(this, res)
            % import and append from a residual object file
            
            % 1) Remove the old overlapped residuals
            % Find the first epoch to remove from the existing res
            if ~res.time.isEmpty
                if isempty(this.time)
                    id_start = this.time.length + 1;
                    id_stop = 0;
                    
                else
                    id_start = find(this.time.getNominalTime.getRefTime(res.time.getNominalTime.first.getMatlabTime) > 0, 1, 'first');
                    id_stop = find(this.time.getNominalTime.getRefTime(res.time.getNominalTime.last.getMatlabTime) <= 0, 1, 'last');
                    
                    if isempty(id_start)
                        id_start = this.time.length + 1;
                    end
                    if isempty(id_stop)
                        id_stop = 0;
                    end
                end
                % Find the last epoch to remove from the existing res
                id_ko = id_start : id_stop;
                if ~isempty(id_ko)
                    this.remEpoch(id_ko);
                end
                
                % 2) Insert data
                
                % Insert time
                time = this.time.getEpoch(1 : id_start-1);
                time.append(res.time);
                time.append(this.time.getEpoch(id_start : this.time.length));
                
                if isempty(this.prn)
                    code_old = [];
                else
                    code_old = Constellation_Collector.obsCode2num(this.obs_code, this.prn);
                    this.remEntry(code_old == 0);
                end
                if isempty(res.prn)
                    code_new = [];
                else
                    code_new = Constellation_Collector.obsCode2num(res.obs_code, res.prn);
                    res.remEntry(code_new == 0);
                end
                
                % new satellites to add
                [code_add, id_add] = setdiff(code_new, code_old);
                [code_common, id_new, id_old] = intersect( code_new, code_old);
                
                n_obs_new = size(res.value, 1);
                
                % resize data to add new epochs
                this.value = [this.value(1 : id_start-1, :); nan(n_obs_new, size(this.value, 2)); this.value(id_start : this.time.length, :)];
                
                % resize data to add columns for the new observables
                this.value = [this.value nan(size(this.value, 1), numel(code_add))];
                
                % add new data
                this.value(id_start + (0 : n_obs_new - 1) , [id_old; ((end - numel(code_add) + 1) : end)']) = res.value(:, [id_new; id_add]);
                
                
                nch_old = size(this.obs_code,2);
                nch_new = size(res.obs_code,2);
                if nch_old == nch_new % length of obs code might not be compatible
                    this.obs_code = [this.obs_code; res.obs_code(id_add, :)];
                elseif nch_old > nch_new
                    this.obs_code = [this.obs_code ; [res.obs_code(id_add, :) char(32*ones(size(res.obs_code,1),nch_old-nch_new,'uint8'))]];
                else
                    this.obs_code = [[this.obs_code char(32*ones(size(this.obs_code,1),nch_new-nch_old,'uint8'))]; res.obs_code(id_add, :)];
                end
                this.prn = [this.prn; res.prn(id_add)];
                
                this.time = time;
                
                this.type = res.type;       % save the last type
                this.rec_coo = res.rec_coo.getCopy; % save the last coo
            end
        end
        
        function remEpoch(this, lid_ko)
            % Remove an epoch from the residuals
            %
            % INPUT
            %   lid_ko  logical array of ko epochs
            %
            % SYNTAX
            %   this.remEpoch(lid_ko);
            
            if ~islogical(lid_ko)
                id_ko = lid_ko;
                lid_ko = false(1,this.time.length);
                lid_ko(id_ko) = true;
            end
            id_ok = ~lid_ko;
            this.time = this.time.getEpoch(id_ok);
            this.value(lid_ko, :) = [];
        end
        
        function remEntry(this, lid_ko)
            % Remove an enntry from the residuals
            %
            % INPUT
            %   lid_ko  logical array of ko entry
            %
            % SYNTAX
            %   this.remEntry(lid_ko);
            
            if ~islogical(lid_ko)
                id_ko = lid_ko;
                lid_ko = false(1, size(this.value,2));
                lid_ko(id_ko) = true;
            end
            this.value(:, lid_ko) = [];
            this.prn(lid_ko) = [];
            this.obs_code(lid_ko, :) = [];
        end
        
        function cutEpochs(this, new_lim)
            % Get the residual only in the time span given
            %
            % SYNTAX
            %   this.cutEpochs(new_limits)
            
            time_res = this.time.getNominalTime();
            sss_start = find(time_res.getMatlabTime >= round(new_lim.first.getMatlabTime * 86400 * time_res.getRate) / (86400 * time_res.getRate), 1, 'first');
            sss_stop = find(time_res.getMatlabTime > round(new_lim.last.getMatlabTime * 86400 * time_res.getRate) / (86400 * time_res.getRate), 1, 'first');
            
            lid = false(this.time.length(), 1);
            if ~isempty(sss_start)
                lid(1 : (sss_start - 1)) = true;
            end
            if ~isempty(sss_stop)
                lid(sss_stop : end) = true;
            end
            this.remEpoch(lid);
        end
    end
    
    % =========================================================================
    %%  GETTERS
    % =========================================================================
    methods
        function [res, obs_code, prn] = getU1(this, sys_c, freq_c)
            % Get residual matrix of the combined/single_freq processing
            %
            % SYNTAX
            %    [res, obs_code, prn] = this.getU1()
            if this.type < 3
                id_ok =  true(numel(this.prn), 1);
                if nargin > 1 && ~isempty(sys_c)
                    id_ok(this.obs_code(:,1) ~= sys_c) = false;
                end
                if nargin > 2 && ~isempty(freq_c)
                    id_ok(this.obs_code(:,3) ~= freq_c) = false;
                end
                
                prn = this.prn(id_ok);
                obs_code = this.obs_code(id_ok,:);
                res = this.value(:, id_ok);
            else
                prn = [];
                obs_code = '';
                res = [];
            end
        end
        
        function [res, obs_code, prn] = getPrU2(this, sys_c, freq_c)
            % Get residual matrix of the uncombined processing
            % Pseudo-codes residuals
            %
            % SYNTAX
            %    [res, obs_code, prn] =  = this.getPrU2()
            if this.type == 3
                id_ok =  this.obs_code(:,2) == 'C';
                if nargin > 1 && ~isempty(sys_c)
                    id_ok(this.obs_code(:,1) ~= sys_c) = false;
                end
                if nargin > 2 && ~isempty(freq_c)
                    id_ok(this.obs_code(:,3) ~= freq_c) = false;
                end
                
                prn = this.prn(id_ok);
                obs_code = this.obs_code(id_ok,:);
                res = this.value(:, id_ok);
            else
                prn = [];
                obs_code = '';
                res = [];
            end
        end
        
        function [is_ph] = isPhase(this)
            % get an index thta tell which resiual are phase
            %
            % SYNTAX:
            %    [is_ph] = this.isPhase()
            is_ph = this.obs_code(:,2) == 'L';
        end
        
        function [is_co] = isCombined(this)
            % get an index thta tell which resiual are phase
            %
            % SYNTAX:
            %    [is_ph] = this.isPhase()
            if size( this.obs_code,2) > 4
                is_co = this.obs_code(:,5) == ' ';
            else
                is_co = false(size( this.obs_code,1),1);
            end
        end
        
        function [res, obs_code, prn] = getPhU2(this, sys_c, freq_c)
            % Get residual matrix of the combined processing
            % Carrier-phase residuals
            %
            % SYNTAX
            %    [res, obs_code, prn] = this.getPhU2()
            if this.type == 3
                id_ok =  this.obs_code(:,2) == 'L';
                if nargin > 1 && ~isempty(sys_c)
                    id_ok(this.obs_code(:,1) ~= sys_c) = false;
                end
                if nargin > 2 && ~isempty(freq_c)
                    id_ok(this.obs_code(:,3) ~= freq_c) = false;
                end
                
                prn = this.prn(id_ok);
                obs_code = this.obs_code(id_ok,:);
                res = this.value(:, id_ok);
            else
                prn = [];
                obs_code = '';
                res = [];
            end
        end
        
        function [res, obs_code, prn, time] = getRangeResiduals(this, sys_c)
            % Get range residuals
            %
            % SYNTAX
            %    [res, obs_code, prn] = this.getU1()
            if this.type < 3
                if nargin == 1
                    [res, obs_code, prn] = this.getU1();
                else
                    [res, obs_code, prn] = this.getU1(sys_c);
                end
                time = this.time;
            else
                % To be done!!! in case of uncombined residuals
                prn = [];
                obs_code = '';
                res = [];
                time = GPS_Time;
            end
        end

        function [res, obs_code, prn, type] = get(this, sys_c, freq_c)
            % Get residual matrix stored in residuals
            %
            % INPUT
            %   sys_c   single character describing the constellation e.g. 'G', 'R', 'E', ...
            %   freq_c  single character describing the frequency number e.g. '1', '2', ....
            %
            % SYNTAX
            %    [res, obs_code, prn, type] = this.get(sys_c, freq_c)
            type = this.type;
            switch type
                case 0
                    prn = [];
                    obs_code = '';
                    res = [];
                case {1, 2} % Prepro or unconbined residuals (both uses U1)
                    if nargin == 1
                        [res, obs_code, prn] = this.getU1();
                    elseif nargin == 2
                        [res, obs_code, prn] = this.getU1(sys_c);
                    elseif nargin == 3
                        [res, obs_code, prn] = this.getU1(sys_c, freq_c);
                    end
                case 3 % if I have uncombined residuals, return just phases
                    if nargin == 1
                        [res, obs_code, prn] = this.getPhU2();
                        if isempty(res)
                            [res, obs_code, prn] = this.getPrU2();
                        end
                    elseif nargin == 2
                        [res, obs_code, prn] = this.getPhU2(sys_c);
                        if isempty(res)
                            [res, obs_code, prn] = this.getPrU2(sys_c);
                        end
                    elseif nargin == 3
                        [res, obs_code, prn] = this.getPhU2(sys_c, freq_c);
                        if isempty(res)
                            [res, obs_code, prn] = this.getPrU2(sys_c, freq_c);
                        end
                    end
            end
        end
        
        function res = getCopy(this)
            % Get a copy of the object
            %
            % SYNTAX
            %   res = this.getCopy();
            
            res = Residuals;
            res.importFromStruct(this.toStruct);
            res.time = this.time.getCopy;
            res.rec_coo = this.rec_coo.getCopy;
        end
        
        function sigma = getStd(this)
            % Get std of all the stored residuals
            % WARNING this is still a very rough estimation, 
            %         different frequencies have different noise
            %
            % SINTAX
            %   sigma = this.getStd()
            sigma = std(zero2nan(serialize(this.get())), 'omitnan');
        end
        
    end
    
    % =========================================================================
    %%  AUXILLIARY
    % =========================================================================
    methods
        function [az, el, sat_coo, sat_name, go_id] = getAzimuthElevation(this)
            % Get azimuth and elevation of each satellite stored in residuals
            %
            %
            %   [az, el, sat_coo, sat_name, go_id] = this.getAzimuthElevation();
            
            core = Core.getCurrentCore;
            sky = core.sky;
            if isempty(core.state.eph_name)
                fw = File_Wizard(Core.getCurrentSettings);
                fw.conjureNavFiles(this.time.first, this.time.last);
            end
            lim = this.time.first.getCopy;
            lim.append(this.time.last);
            core.initSkySession(lim);
            cc = core.getConstellationCollector;
            go_id = unique(cc.getIndex(this.obs_code(:,1), this.prn));
            sat_name = cc.getSatName(go_id);
            [az, el, sat_coo] = sky.getAzimuthElevation(this.rec_coo, this.time, go_id);
        end
    end
    
    % =========================================================================
    %%  MULTIPATH
    % =========================================================================
    methods
        function ant_mp = computeMultiPath(this, marker_name, l_max, flag_reg, is_ph, mode)
            % Get multi path maps in different modes
            %
            % z_map  Zernike
            % r_map  Zernike + (the methods specified on mode)
            % g_map  Simple Gridding of size [stk_grid_step]
            % c_map  Congruent cells gridding of size [stk_grid_step]
            % g1_map Simple Gridding of size [1x1]
            % c1_map Congruent cells gridding of size [1x1]
            %
            %
            % INPUT
            %   marker_name     name of the station
            %   l_max           maximum degree of the 3 steps of the zernike interpolation [l_max1, l_max2, l_max3]
            %                   to disable Zernike use l_max = 0
            %   flag_reg        add regularization points (pseudo obs at 0) in the empty areas of the sky
            %   is_ph           use phases instead of pseudo-ranged (default = true)
            %   mode            - 0 use Zernike + staking maps using congruent cells (variable azimuthal resolution)
            %                                   n_step_az = (360/max_n_step_az * cosd(el))
            %                   - [0, n, m] use Zernike + stacking maps using congruent cells with maximum size of [n x m] note that the output will always be 0.5 x 0.5 degrees 
            %                   - 1 use Zernike + stacking map of [5 x 1] 
            %                   - [1, n, m] use Zernike + stacking map of [n x m]
            % NOTE
            %   For mode 0 and 1 the output matrix will always be 0.5 x 0.5 degrees 
            %
            % SYNTAX
            %   this.computeMultiPath(marker_name, <l_max=[43,43,43]>, <flag_reg=true>, <is_ph=true>, <mode=[0 5 1]>)
            
            state = Core.getCurrentSettings;
            flag_discard_co = false;
            
            if nargin < 6 || isempty(mode)
                mode = 0; % Z + stacking
            end
            if numel(mode) == 3
                r_grid_step = mode(2,3);
                mode = mode(1);
            else
                r_grid_step = state.mp_zcongruent_up_nxm;
            end            
            
            n_min = state.mp_n_min; % minimum  number of points per cell
            n_min = 15; %  <== DEBUG
            
            % z_map  Zernike
            % r_map  Zernike + (the methods specified on mode)
            % g_map  Simple Gridding of size [stk_grid_step]
            % c_map  Congruent cells gridding of size [stk_grid_step]
            % g1_map Simple Gridding of size [1x1]
            % c1_map Congruent cells gridding of size [1x1]
            ltype_of_grids = [...
                sum(state.mp_l_max) > 0 ...
                state.mp_zcongruent_up_nxm(1) > 0 ...
                state.mp_regular_up_nxm(1) > 0 ...
                state.mp_congruent_up_nxm(1) > 0 ...
                state.mp_regular_nxm(1) > 0 ...
                state.mp_congruent_nxm(1) > 0];
                        
            log = Core.getLogger();
            ant_mp = struct();
            if this.isEmpty
                log.addWarning('Residuals have not been computed');
            else
                flag_debug = false;
                if flag_debug
                    % Enable all the grid types
                    ltype_of_grids = logical([1 1 1 1 1 1]); % All enabled
                end
                
                if nargin < 3 || isempty(l_max)
                    l_max = state.mp_l_max;
                end
                % Legacy support
                if numel(l_max) == 3
                    l_max = [l_max(1) l_max];
                end
                if numel(l_max) == 1
                    l_max = [l_max 0 0 0];
                end
                
                % Depending on the maximum zernike degree change the map resolution
                if max(l_max) > 180
                    grid_step = 0.1;
                elseif max(l_max) > 90                    
                    grid_step = 0.25;
                else
                    grid_step = 0.5;
                end            
                
                if nargin < 2
                    marker_name = 'UNKN';
                end
                
                if nargin < 4 || isempty(flag_reg)
                    flag_reg = true;
                end
                
                if nargin < 5 || isempty(is_ph)
                    % If there are phases use phases
                    is_ph = any((serialize(this.obs_code(:,2:3:end-1))) == 'L');
                end
                if is_ph
                    name = 'Carrier-phase residuals';
                    search_obs = 'L';
                else
                    name = 'Pseudo-ranges residuals';
                    search_obs = 'C';
                end
                
                deg2rad = pi/180;
                
                cc = Core.getConstellationCollector();
                [az, el, ~, ~, go_id] = this.getAzimuthElevation();
                sys_c_list = cc.getAvailableSys;
                
                log = Core.getLogger;
                log.addMarkedMessage(sprintf('Computing multipath mitigation coefficients for "%s"', marker_name));
                
                obs_code = this.obs_code;
                if Core.getCurrentSettings.FLAG_MP_IGNORE_TRK
                    for i = 1 : size(obs_code, 1)
                        obs_code(i, 4:3:end) = '_';
                    end
                end
                
                for sys_c = sys_c_list(:)'
                    ids = find(obs_code(:,1) == sys_c & any((obs_code(:,2:3:end-1)) == search_obs, 2));
                    if ~any(ids)
                        log.addWarning(sprintf('No %s found in %s for constellation %s', name, marker_name, cc.getSysName(sys_c)));
                    else                        
                        obs_id_num = cc.obsCode2num(obs_code(ids, :), zeros(size(ids, 1), 1)); % get all the data of the same frequency - all the satellites
                        uobs_id = unique(obs_id_num);
                        for  t = 1 : numel(uobs_id)
                            id = ids(obs_id_num == uobs_id(t)); % tracking for the specific obs_code
                            trk_code = obs_code(id(1),2:end);
                            
                            data_found = false;
                            
                            res = zero2nan(this.value(:, id));
                            res_go_id = cc.getIndex(obs_code(id, 1), this.prn(id));
                            
                            % Get all the data to interpolate
                            [~, id_sat] = ismember(res_go_id,go_id);

                            az_all = [];
                            el_all = [];
                            
                            % Propagate orbit nans
                            for s = 1 : numel(res_go_id)
                                id_ko = isnan(el(:,id_sat(s))) | isnan(az(:,id_sat(s)));
                                res(id_ko,s) = nan;
                            end
                            
                            res_all = res(~isnan(res(:)));
                            res_smt = Receiver_Commons.smoothMat(res, 'spline', 120/this.time.getRate);
                            res_smt = res_smt(~isnan(res(:)));
                            
                            go_id_list = [];
                            for s = 1 : numel(res_go_id)
                                id_ok = ~isnan(res(:,s));
                                if any(id_ok)
                                    data_found = true;
                                    % res_all = [res_all; serialize(res(id_ok, s))]; %#ok<AGROW>
                                    az_all = [az_all; az(id_ok, id_sat(s)) .* deg2rad]; %#ok<AGROW>
                                    el_all = [el_all; el(id_ok, id_sat(s)) .* deg2rad]; %#ok<AGROW>
                                    go_id_list = [go_id_list; res_go_id(s)]; %#ok<AGROW>
                                end
                            end
                            
                            if data_found
                                m_max = l_max;                                
                                % Remove outliers
                                id_ok = Core_Utils.polarCleaner(az_all, el_all, res_all, [360, 1]) & Core_Utils.polarCleaner(az_all, el_all, res_smt, [360, 1]);
                                log.addMessage(log.indent(sprintf('1. Outlier rejection (%.3f%%)', (sum(~id_ok) / numel(id_ok)) * 100), 9));
                                if flag_debug
                                    figure; plot(el_all/pi*180, res_all*1e3, '.'); hold on; plot(el_all(~id_ok)/pi*180, res_all(~id_ok)*1e3, 'o');
                                    legend('residuals', 'outliers');
                                    title((sprintf('Residuals of  %s %s%s [mm]', marker_name, sys_c, trk_code))); drawnow
                                    grid on;
                                end
                                clear res_smt;
                                az_all = az_all(id_ok);
                                el_all = el_all(id_ok);
                                res_all = res_all(id_ok);
                                n_obs = numel(res_all);
                                
                                % 3 sigma filter per latitude
                                if flag_reg
                                    log.addMessage(log.indent('2. Preparing regularization', 9));
                                    % Get regularization points based on empty sky areas
                                    [data_map, n_data_map, az_grid, el_grid] = Core_Utils.hemiGridder(az_all, el_all, res_all, [1 1]);
                                    [az_grid, el_grid] = meshgrid(az_grid, el_grid);
                                    az_reg = az_grid(n_data_map <= n_min);
                                    el_reg = el_grid(n_data_map <= n_min);
                                    
                                    % In the knots with few data add zero
                                    az_all = [az_all; az_reg];
                                    el_all = [el_all; el_reg];
                                    res_all = [res_all; zeros(size(el_reg))];
                                    
                                    if flag_discard_co
                                        % First approach
                                        % Do not consider observations under
                                        % cut-off => map the radius starting from cut-off
                                        id_ko = (el_all * 180/pi) < Core.getState.getCutOff;
                                        az_all(id_ko) = [];
                                        el_all(id_ko) = [];
                                        res_all(id_ko) = [];
                                    else
                                        % Second approach
                                        % Add additional points at the border close to radius 1
                                        % This regularization is needed if the mapping of the radius
                                        % have a cut-off
                                        for i = 0 : 0.5 : (Core.getState.getCutOff - 2.5)
                                            az_all = [az_all; (-pi : 0.05 : pi)'];
                                            el_all = [el_all; i/180*pi + (-pi : 0.05 : pi)'*0];
                                            res_all = [res_all; (-pi : 0.05 : pi)'*0];
                                        end
                                    end
                                end
                                
                                % Ignore data under cut-off
                                if flag_discard_co
                                    Zernike.setCutOff(max(0, (Core.getState.getCutOff - 0.5)));
                                end
                                res_work = res_all;
                                
                                if ~any(ltype_of_grids(1:2))
                                    % Zernike maps are not requested
                                    z_map = 0;
                                    r_map = 0;
                                else
                                    
                                    % Perform the first of 3 Zernike steps
                                    Zernike.setMode(0); % Set Zernike engine to recursive
                                    if l_max(1) > 0
                                        log.addMessage(log.indent(sprintf('%d. Zernike coef. estimation (l_max = %d) (1/4)', 2 + flag_reg*1, l_max(1)), 9));
                                        Zernike.setModeMF(0);
                                        el2radius = Zernike.getElFun;
                                        [z_par, l, m] = Zernike.analysisAllBlock(l_max(1), m_max(1), az_all, el2radius(el_all), res_work, 1e-5);
                                        [z_map1, az_grid, el_grid] = Zernike.synthesisGrid(l, m, z_par, grid_step);
                                        z_map1((el_grid * 180/pi) < Core.getState.getCutOff, :) = 0; % remove cutoff;
                                        res_work = res_work - Core_Utils.hgrid2scatter(az_all, el_all, z_map1, false, 'spline');
                                    else
                                        z_map1 = 0;
                                    end
                                    
                                    Zernike.setMode(0); % Set Zernike engine to recursive
                                    if l_max(2) > 0
                                        log.addMessage(log.indent(sprintf('%d. Zernike coef. estimation (l_max = %d) (1/4)', 2 + flag_reg*1, l_max(1)), 9));
                                        Zernike.setModeMF(1);
                                        el2radius = Zernike.getElFun;
                                        [z_par, l, m] = Zernike.analysisAllBlock(l_max(2), m_max(2), az_all, el2radius(el_all), res_work, 1e-5);
                                        [z_map2, az_grid, el_grid] = Zernike.synthesisGrid(l, m, z_par, grid_step);
                                        z_map2((el_grid * 180/pi) < Core.getState.getCutOff, :) = 0; % remove cutoff;
                                        res_work = res_work - Core_Utils.hgrid2scatter(az_all, el_all, z_map2, false, 'spline');
                                    else
                                        z_map2 = 0;
                                    end
                                    
                                    % Perform the second of 3 Zernike steps
                                    if l_max(3) > 0
                                        log.addMessage(log.indent(sprintf('%d. Zernike coef. estimation (l_max = %d) (2/4)', 2 + flag_reg*1, l_max(2)), 9));
                                        Zernike.setModeMF(2);
                                        Zernike.setCutOff(0); % This mapping function is more unstable at low elevations => use polar regularization
                                                                               
                                        el2radius = Zernike.getElFun;
                                        [z_par, l, m] = Zernike.analysisAllBlock(l_max(3), m_max(3), az_all, el2radius(el_all), res_work, 1e-5);
                                        [z_map3, az_grid, el_grid] = Zernike.synthesisGrid(l, m, z_par, grid_step);
                                        z_map3((el_grid * 180/pi) < Core.getState.getCutOff, :) = 0; % remove cutoff;
                                        res_work = res_work - Core_Utils.hgrid2scatter(az_all, el_all, z_map3, false, 'spline');
                                    else
                                        z_map3 = 0;
                                    end
                                    
                                    % Perform the third of 3 Zernike steps
                                    if l_max(4) > 0
                                        log.addMessage(log.indent(sprintf('%d. Zernike coef. estimation (l_max = %d) (3/4)', 2 + flag_reg*1, l_max(3)), 9));
                                        Zernike.setModeMF(3);
                                        Zernike.setCutOff(0); % This mapping function is more unstable at low elevations => use polar regularization
                                        [z_par, l, m] = Zernike.analysisAllBlock(l_max(4), m_max(4), az_all, el2radius(el_all), res_work, 1e-5);
                                        [z_map4, az_grid, el_grid] = Zernike.synthesisGrid(l, m, z_par, grid_step);
                                        z_map4((el_grid * 180/pi) < Core.getState.getCutOff, :) = 0; % remove cutoff;
                                        res_work = res_work - Core_Utils.hgrid2scatter(az_all, el_all, z_map4, false, 'spline');
                                    else
                                        z_map4 = 0;
                                    end
                                    
                                    % Generate maps
                                    log.addMessage(log.indent(sprintf('%d. Compute mitigation grids', 4 + flag_reg*1), 9));
                                    z_map = z_map1 + z_map2 + z_map3 + z_map4; % z_map  Zernike only
                                    
                                    if ~ltype_of_grids(2) % r_map  Zernike + (the methods specified on mode)
                                        r_map = 0;
                                    else
                                        % In this mode a gridding on the residuals is performed
                                        % to retrieve the high frequency multipath a double step procedure is performed, 
                                        % first low-res than high-res
                                        
                                         if mode == 1
                                            flag_congruent = false;
                                        else % if mode == 0
                                            flag_congruent = true;
                                        end
                                        
                                        
                                        % low-res step
                                        grid_step1 = [max(r_grid_step(1), 1.5) max(r_grid_step(end), 0.5)];
                                        out_step1 = [min([grid_step1(1), grid_step(1), 0.5]) min([grid_step1(end), grid_step(end) 0.5])];
                                        
                                        % high-res step
                                        grid_step2 = r_grid_step;
                                        out_step2 = [min([r_grid_step(1), grid_step(1)]) min([r_grid_step(end), grid_step(end)])];

                                        % Set the final size of the r_grid (it is generally [0.25 x 0.1] deg
                                        out_size = [90 360] ./ fliplr(out_step2);
                                        
                                        % Restore data with no regularization
                                        res_work = res_all(1 : n_obs);
                                        res_work = res_work - Core_Utils.hgrid2scatter(az_all(1 : n_obs), el_all(1 : n_obs), Core_Utils.resize2(z_map, out_size));
                                       
                                        % Compute low-res grid - STEP 1
                                        [r_map1, n_data_map, az_grid, el_grid] = Core_Utils.hemiGridder(az_all(1 : n_obs), el_all(1 : n_obs), res_work, grid_step1, out_step1, flag_congruent, n_min);
                                                                                
                                        % Resize it to the final grid size
                                        r_map = Core_Utils.resize2(r_map1, out_size);
   
                                        res_work = res_work - Core_Utils.hgrid2scatter(az_all(1 : n_obs), el_all(1 : n_obs), r_map);

                                        % Compute high-res grid - STEP 2 -- this is not CONGRUENT! 
                                        % note that high latitude cells are usually under n_min thr
                                        [r_map2, n_data_map, az_grid, el_grid] = Core_Utils.hemiGridder(az_all(1 : n_obs), el_all(1 : n_obs), res_work, grid_step2, grid_step2, false, 7);
                                        r_map = Core_Utils.resize2(z_map, out_size) + r_map + Core_Utils.resize2(r_map2, out_size);
                                    end
                                end
                                                                
                                % Compute normal and congruent maps as comparison (no regularization)
                                if ltype_of_grids(3) % g_map  Simple Gridding of size
                                    g_map = Core_Utils.hemiGridder(az_all, el_all, res_all, state.mp_regular_up_nxm , grid_step, false, n_min);
                                else
                                    g_map = 0;
                                end
                                if ltype_of_grids(4) % c_map  Congruent cells gridding of size [stk_grid_step]
                                    c_map = Core_Utils.hemiGridder(az_all, el_all, res_all, state.mp_congruent_up_nxm, grid_step, true, n_min);
                                else
                                    c_map = 0;
                                end
                                if ltype_of_grids(5) % g1_map Simple Gridding of size [1x1]
                                    g1_map = Core_Utils.hemiGridder(az_all(res_all ~= 0), el_all(res_all ~= 0), res_all(res_all ~= 0), state.mp_regular_nxm, [], false, n_min);
                                else
                                    g1_map = 0;
                                end
                                if ltype_of_grids(6) % c1_map Congruent cells gridding of size [1x1]
                                    c1_map = Core_Utils.hemiGridder(az_all(res_all ~= 0), el_all(res_all ~= 0), res_all(res_all ~= 0), state.mp_congruent_nxm, state.mp_congruent_nxm, true, n_min);
                                else
                                    c1_map = 0;
                                end
                                
                                if flag_debug
                                    clim = [-1 1] * max(-perc(1e3*(r_map(:)), 0.003),perc(1e3*(r_map(:)), 0.997));
                                    mp_map = z_map1;
                                    [az_grid, el_grid] = Core_Utils.getPolarGrid(360 / size(mp_map, 2), 90 / size(mp_map, 1));
                                    az_grid = Core_Utils.deg2rad(az_grid)';
                                    el_grid = Core_Utils.deg2rad(el_grid);
                                    
                                    %figure; imagesc(1e3*(z_map)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                    if numel(z_map1) > 1
                                        figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map1)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    end
                                    title((sprintf('Zernike expansion (1) of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow
                                    
                                    if numel(z_map2) > 1
                                        figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map2)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    end
                                    title((sprintf('Zernike expansion (2) of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow
                                    
                                    if numel(z_map3) > 1
                                        figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map3)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    end
                                    title((sprintf('Zernike expansion (3) of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow
                                    
                                    if numel(z_map4) > 1
                                        figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map4)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    end
                                    title((sprintf('Zernike expansion (4) of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow

                                    figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    title((sprintf('Zernike expansion of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow

                                    mp_map = r_map;
                                    [az_grid, el_grid] = Core_Utils.getPolarGrid(360 / size(mp_map, 2), 90 / size(mp_map, 1));
                                    az_grid = Core_Utils.deg2rad(az_grid)';
                                    el_grid = Core_Utils.deg2rad(el_grid);
                                                                        
                                    figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*r_map); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    title((sprintf('Final map of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow
                                    
                                    mp_map = c_map;
                                    [az_grid, el_grid] = Core_Utils.getPolarGrid(360 / size(mp_map, 2), 90 / size(mp_map, 1));
                                    az_grid = Core_Utils.deg2rad(az_grid)';
                                    el_grid = Core_Utils.deg2rad(el_grid);
                                    figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(mp_map)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    title((sprintf('Gridded map with congruent cells of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow

                                    mp_map = g_map;
                                    [az_grid, el_grid] = Core_Utils.getPolarGrid(360 / size(mp_map, 2), 90 / size(mp_map, 1));
                                    az_grid = Core_Utils.deg2rad(az_grid)';
                                    el_grid = Core_Utils.deg2rad(el_grid);
                                    figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(mp_map)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    title((sprintf('Gridded map of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow
                                    
                                    mp_map = c1_map;
                                    [az_grid, el_grid] = Core_Utils.getPolarGrid(360 / size(mp_map, 2), 90 / size(mp_map, 1));
                                    az_grid = Core_Utils.deg2rad(az_grid);
                                    el_grid = Core_Utils.deg2rad(el_grid);
                                    
                                    figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(mp_map)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    title((sprintf('Gridded map with congruent cells [1 x 1] of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow

                                    mp_map = g1_map;
                                    [az_grid, el_grid] = Core_Utils.getPolarGrid(360 / size(mp_map, 2), 90 / size(mp_map, 1));
                                    az_grid = Core_Utils.deg2rad(az_grid)';
                                    el_grid = Core_Utils.deg2rad(el_grid);
                                    figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(mp_map)); colormap((Cmap.get('PuOr', 2^11))); caxis(clim); colorbar;
                                    title((sprintf('Gridded map [1 x 1] of %s %s%s [mm]', marker_name, sys_c, trk_code)), 'interpreter', 'none'); drawnow
                                end
                                
                                if ~isfield(ant_mp, sys_c)
                                    ant_mp.(sys_c) = struct;
                                end
                                if ~isfield(ant_mp.(sys_c), trk_code)
                                    trk_code = strrep(trk_code, ' ', '_'); % structures do not support spaces
                                    ant_mp.(sys_c).(trk_code) = struct;
                                end
                                % Keep multiple solutions in the struct
                                % decide a-posteriori what it's better
                                
                                % Save grids of multi-path
                                ant_mp.(sys_c).(trk_code).z_map = single(z_map);             % Zernike map
                                ant_mp.(sys_c).(trk_code).r_map = single(r_map);             % Zernike map + gridded residuals
                                ant_mp.(sys_c).(trk_code).g_map = single(g_map);             % Simple Gridding of size [stk_grid_step]
                                ant_mp.(sys_c).(trk_code).c_map = single(c_map);             % Congruent cells gridding of size [stk_grid_step]
                                ant_mp.(sys_c).(trk_code).g1_map = single(g1_map);           % Simple Gridding of size [1x1]
                                ant_mp.(sys_c).(trk_code).c1_map = single(c1_map);           % c1_map Congruent cells gridding of size [1x1]
                                % Save Zernike coefficients
                                %ant_mp.(sys_c).(trk_code).z_par = [z_par1 z_par2];
                                %ant_mp.(sys_c).(trk_code).l = l;
                                %ant_mp.(sys_c).(trk_code).m = m;
                            else
                                if ~data_found
                                    log = Core.getLogger;
                                    log.addWarning(sprintf('No %s %s found in %s for constellation %s', name, trk_code, marker_name, cc.getSysName(sys_c)));
                                end
                                
                            end
                        end
                    end
                end
                % Get the time limit of the map solution
                ant_mp.time_lim = this.time.getEpoch([1 this.time.length]);
            end
        end
    end
    
    % =========================================================================
    %%  SHOW
    % =========================================================================
    methods
        function fh_list = showResSkyCartScatter(this, marker_name, sys_c_list, is_ph)
            % Plot residuals of the solution on cartesian axes
            %
            % SYNTAX
            %   this.showResSkyCartScatter(marker_name, sys_c_list, is_ph)
            log = Core.getLogger();
            fh_list = [];
            if this.isEmpty
                log.addWarning('Residuals have not been computed');
            else
                if nargin < 3 || isempty(sys_c_list)
                    sys_c_list = unique(this.obs_code(:,1));
                end
                if nargin < 4 || isempty(is_ph)
                    % If there are phases use phases
                    is_ph = any((serialize(this.obs_code(:,2:3:end-1))) == 'L');
                end
                if is_ph
                    name = 'Carrier-phase residuals';
                    search_obs = 'L';
                    scale = 1e3;
                else
                    name = 'Pseudo-ranges residuals';
                    search_obs = 'C';
                    scale = 1e2;
                end
                
                if nargin < 2
                    marker_name = '';
                end
                
                cc = Core.getConstellationCollector();
                [az, el, sat_coo, sat_name, go_id] = this.getAzimuthElevation();
                for sys_c = sys_c_list(:)'
                    ids = find(this.obs_code(:,1) == sys_c & any((this.obs_code(:,2:3:end-1)) == search_obs, 2));
                    if ~any(ids)
                        log.addWarning(sprintf('No %s found in %s for constellation %s', name, marker_name, cc.getSysName(sys_c)));
                    else
                        obs_id_num = cc.obsCode2num(this.obs_code(ids,:), zeros(size(ids, 1), 1));
                        uobs_id = unique(obs_id_num);
                        for  t = 1 : numel(uobs_id)
                            id = ids(obs_id_num == uobs_id(t)); % tracking for the specific obs_code
                            trk_code = this.obs_code(id(1),:);
                            
                            res = this.value(:, id) * scale;
                            res_go_id = cc.getIndex(this.obs_code(id, 1), this.prn(id));
                            
                            fh = figure('Visible', 'off'); fh.Name = sprintf('%03d: %s Res Cart %s', fh.Number, marker_name, strtrim(trk_code)); fh.NumberTitle = 'off';
                            Core_UI.beautifyFig(fh);
                            
                            fh_list = [fh_list; fh]; %#ok<AGROW>
                            fig_name = sprintf('Res_cart_%s_%s_%s_%s', marker_name, cc.getSysName(sys_c), this.time.first.toString('yyyymmdd_HHMM'));
                            fh.UserData = struct('fig_name', fig_name);
                            
                            [~, id_sat] = intersect(go_id, res_go_id);
                            
                            figure(fh); % get focus;
                            hold on;
                            go_id_list = [];
                            for s = 1 : numel(res_go_id)
                                id_ok = find(res(:,s) ~= 0);
                                if any(id_ok)
                                    [~, id_sort] = sort(abs(res(id_ok,s)));
                                    id_ok = id_ok(id_sort);
                                    line = scatter(az(id_ok, id_sat(s)), el(id_ok, id_sat(s)), 45, serialize(res(id_ok, s)), 'filled');
                                    line.UserData = res_go_id(s);
                                    go_id_list = [go_id_list; res_go_id(s)];
                                end
                            end
                            if isempty(go_id_list)
                                delete(fh);
                            else
                                ylim([0 90]); xlim([-180 180]);
                                caxis([-1 1] * max(2, min(6*std(noZero(res),'omitnan'), max(abs(noZero(res(:)))))));
                                colormap((Cmap.get('PuOr', 2^11)));
                                fh.Color = [.95 .95 .95]; cb = colorbar(); cbt = title(cb, iif(scale == 1e2, '[cm]', '[mm]')); cbt.Parent.UserData = cbt; ax = gca; ax.Color = 'none';
                                h = title(sprintf('Satellites residuals - receiver %s - %s\\fontsize{5} \n', strrep(marker_name,'_','\_'), cc.getSysExtName(sys_c)));  h.FontWeight = 'bold';
                                hl = xlabel('Azimuth [deg]'); hl.FontWeight = 'bold';
                                hl = ylabel('Elevation [deg]'); hl.FontWeight = 'bold';
                                
                                Core_UI.addSatMenu(fh, go_id_list);
                                Core_UI.beautifyFig(fh);
                                Core_UI.addExportMenu(fh);
                                Core_UI.addBeautifyMenu(fh);
                                fh.Visible = 'on'; drawnow;
                            end
                        end
                    end
                end
                if isempty(fh_list)
                    log.addWarning('Residuals have not been computed');
                end
            end
        end
        
        function fh_list = showResSkyPolarScatter(this, marker_name, sys_c_list, is_ph)
            % Plot residuals of the solution on polar axes
            %
            % SYNTAX
            %   this.showResSkyPolarScatter(marker_name, sys_c_list, is_ph)
            log = Core.getLogger();
            fh_list = [];
            if this.isEmpty
                log.addWarning('Residuals have not been computed');
            else
                if nargin < 3 || isempty(sys_c_list)
                    sys_c_list = unique(this.obs_code(:,1));
                end
                if nargin < 4 || isempty(is_ph)
                    % If there are phases use phases
                    is_ph = any((serialize(this.obs_code(:,2:3:end-1))) == 'L');
                end
                if is_ph
                    name = 'Carrier-phase residuals';
                    search_obs = 'L';
                    scale = 1e3;
                else
                    name = 'Pseudo-ranges residuals';
                    search_obs = 'C';
                    scale = 1e2;
                end
                
                if nargin < 2
                    marker_name = '';
                end
                
                cc = Core.getConstellationCollector();
                [az, el, ~, ~, go_id] = this.getAzimuthElevation();
                for sys_c = sys_c_list(:)'
                    ids = find(this.obs_code(:,1) == sys_c & any((this.obs_code(:,2:3:end-1)) == search_obs, 2));
                    if ~any(ids)
                        log.addWarning(sprintf('No %s found in %s for constellation %s', name, marker_name, cc.getSysName(sys_c)));
                    else
                        obs_id_num = cc.obsCode2num(this.obs_code(ids,:), zeros(size(ids, 1), 1));
                        uobs_id = unique(obs_id_num);
                        for  t = 1 : numel(uobs_id)
                            id = ids(obs_id_num == uobs_id(t)); % tracking for the specific obs_code
                            trk_code = this.obs_code(id(1),:);
                            
                            res = this.value(:, id) * scale;
                            res_go_id = cc.getIndex(this.obs_code(id, 1), this.prn(id));
                            
                            fh = figure('Visible', 'off'); fh.Name = sprintf('%03d: %s Res Polar %s', fh.Number, marker_name, strtrim(trk_code)); fh.NumberTitle = 'off';
                            Core_UI.beautifyFig(fh);
                            
                            fh_list = [fh_list; fh]; %#ok<AGROW>
                            fig_name = sprintf('Res_polar_%s_%s_%s_%s', marker_name, cc.getSysName(sys_c), this.time.first.toString('yyyymmdd_HHMM'));
                            fh.UserData = struct('fig_name', fig_name);
                            
                            [~, id_sat] = ismember(res_go_id,go_id);
                            
                            figure(fh); % get focus;
                            go_id_list = [];
                            for s = 1 : numel(res_go_id)
                                id_ok = find(res(:,s) ~= 0);
                                if any(id_ok)
                                    [~, id_sort] = sort(abs(res(id_ok,s)));
                                    id_ok = id_ok(id_sort);
                                    line = polarScatter(az(id_ok, id_sat(s))/180*pi, (90 -el(id_ok, id_sat(s)))/180*pi, 45, serialize(res(id_ok, s)), 'filled');
                                    line.UserData = res_go_id(s);
                                    hold on;
                                    go_id_list = [go_id_list; res_go_id(s)]; %#ok<AGROW>
                                end
                            end
                            if isempty(go_id_list)
                                delete(fh);
                            else
                                caxis([-1 1] * max(2, min(6*std(noZero(res),'omitnan'), max(abs(noZero(res(:)))))));
                                colormap((Cmap.get('PuOr', 2^11)));
                                fh.Color = [.95 .95 .95]; cb = colorbar(); cbt = title(cb, iif(scale == 1e2, '[cm]', '[mm]')); cbt.Parent.UserData = cbt; ax = gca; ax.Color = 'none';
                                h = title(sprintf('Satellites residuals\nreceiver %s - %s %s\\fontsize{5} \n', strrep(marker_name,'_','\_'), cc.getSysExtName(sys_c), strtrim(trk_code(2:end))));  h.FontWeight = 'bold';
                                
                                Core_UI.addSatMenu(fh, go_id_list);
                                Core_UI.beautifyFig(fh);
                                Core_UI.addExportMenu(fh);
                                Core_UI.addBeautifyMenu(fh);
                                fh.Visible = 'on'; drawnow;
                            end
                        end
                    end
                end
                if isempty(fh_list)
                    log.addWarning('Residuals have not been computed');
                end
            end
        end
        
        function fh_list = showRes(this, marker_name, sys_c_list, is_ph)
            % Plot residuals of the solution
            %
            % SYNTAX
            %   fh_list = this.showRes(marker_name, sys_c_list, is_ph)
            log = Core.getLogger();
            fh_list = [];
            if this.isEmpty
                log.addWarning('Residuals have not been computed');
            else
                if nargin < 3 || isempty(sys_c_list)
                    sys_c_list = unique(this.obs_code(:,1));
                end
                if nargin < 4 || isempty(is_ph)
                    % If there are phases use phases
                    is_ph = any((serialize(this.obs_code(:,2:3:end-1))) == 'L');
                end
                if is_ph
                    name = 'Carrier-phase residuals';
                    search_obs = 'L';
                    scale = 1e3;
                else
                    name = 'Pseudo-ranges residuals';
                    search_obs = 'C';
                    scale = 1e2;
                end
                
                if nargin < 2
                    marker_name = '';
                end
                
                cc = Core.getConstellationCollector();
                for sys_c = sys_c_list(:)'
                    ids = find(this.obs_code(:,1) == sys_c & any((this.obs_code(:,2:3:end-1)) == search_obs, 2));
                    if ~any(ids)
                        log.addWarning(sprintf('No %s found in %s for constellation %s', name, marker_name, cc.getSysName(sys_c)));
                    else
                        obs_id_num = cc.obsCode2num(this.obs_code(ids,:), zeros(size(ids, 1), 1));
                        uobs_id = unique(obs_id_num);
                        for  t = 1 : numel(uobs_id)
                            id = ids(obs_id_num == uobs_id(t)); % tracking for the specific obs_code
                            trk_code = this.obs_code(id(1),:);
                            res = zero2nan(this.value(:, id) * scale);
                            %res = Receiver_Commons.smoothMat(res, 'spline', 10);
                            res_go_id = cc.getIndex(this.obs_code(id, 1), this.prn(id));
                            
                            fh = figure('Visible', 'off'); fh.Name = sprintf('%03d: %s Res %s %s', fh.Number, marker_name, cc.getSysName(sys_c), strtrim(trk_code(2:end))); fh.NumberTitle = 'off';
                            Core_UI.beautifyFig(fh); Core_UI.beautifyFig(fh); drawnow;
                            
                            fig_name = sprintf('Res_polar_%s_%s_%s_%s', marker_name, cc.getSysName(sys_c), strtrim(trk_code(2:end)), this.time.first.toString('yyyymmdd_HHMM'));
                            fh.UserData = struct('fig_name', fig_name);
                            
                            time = this.time.getMatlabTime;
                            
                            go_id_list = [];
                            sat_name_list = {};
                            figure(fh); % get focus;
                            for s = 1 : numel(res_go_id)
                                id_ok = ~isnan(res(:,s));
                                if any(id_ok)
                                    line = Core_Utils.plotSep(time(id_ok), serialize(res(id_ok, s)), '.-', 'Color', Core_UI.getColor(s, numel(res_go_id)));
                                    line.UserData = res_go_id(s);
                                    hold on;
                                    go_id_list = [go_id_list; res_go_id(s)];
                                    sat_name_list = [sat_name_list {cc.getSatName(res_go_id(s))}];
                                end
                            end
                            if isempty(go_id_list)
                                delete(fh);
                            else
                                xlim(minMax(time));
                                ylim([-1 1] * max(abs(ylim)));
                                setTimeTicks();
                                h = title(sprintf('Receiver %s - %s %s\\fontsize{5} \n', strrep(marker_name,'_','\_'), cc.getSysExtName(sys_c), strtrim(trk_code(2:end))));  h.FontWeight = 'bold';
                                [~, icons] = legend(sat_name_list, 'Location', 'NorthEastOutside');
                                icons = icons(numel(sat_name_list) + 2 : 2 : end);
                                for i = 1 : numel(icons)
                                    icons(i).MarkerSize = 18;
                                    icons(i).LineWidth = 2;
                                end
                                ylabel(sprintf('Satellite Residuals %s', iif(scale == 1e2, '[cm]', '[mm]')));
                                Core_UI.addSatMenu(fh, go_id_list);
                                Core_UI.beautifyFig(fh);
                                Core_UI.addExportMenu(fh);
                                Core_UI.addBeautifyMenu(fh);
                                fh_list = [fh_list; fh]; %#ok<AGROW>
                                fh.Visible = 'on'; drawnow;
                            end
                        end
                    end
                end
                if isempty(fh_list)
                    log.addWarning('Residuals have not been computed');
                end
            end
        end
        
        function fh_list = showResPerSat(this, marker_name, sys_c_list, is_ph)
            % Plot the residuals of phase per tracking
            %
            % INPUT
            %   res     is the matrix of residuals satellite by satellite and can be passed from e.g. NET
            %
            % SYNTAX
            %   fh_list = this.showResPerSat(marker_name, sys_c_list, is_ph)
            
            log = Core.getLogger();
            fh_list = [];
            if this.isEmpty
                log.addWarning('Residuals have not been computed');
            else
                if nargin < 3 || isempty(sys_c_list)
                    sys_c_list = unique(this.obs_code(:,1));
                end
                if nargin < 4 || isempty(is_ph)
                    % If there are phases use phases
                    is_ph = any((serialize(this.obs_code(:,2:3:end-1))) == 'L');
                end
                if is_ph
                    name = 'Carrier-phase residuals';
                    search_obs = 'L';
                    scale = 1e3;
                else
                    name = 'Pseudo-ranges residuals';
                    search_obs = 'C';
                    scale = 1e2;
                end
                
                if nargin < 2
                    marker_name = '';
                end
                
                cc = Core.getConstellationCollector();
                for sys_c = sys_c_list(:)'
                    ids = find(this.obs_code(:,1) == sys_c & any((this.obs_code(:,2:3:end-1)) == search_obs, 2));
                    if ~any(ids)
                        log.addWarning(sprintf('No %s found in %s for constellation %s', name, marker_name, cc.getSysName(sys_c)));
                    else
                        obs_id_num = cc.obsCode2num(this.obs_code(ids,:), zeros(size(ids, 1), 1));
                        uobs_id = unique(obs_id_num);
                        for  t = 1 : numel(uobs_id)
                            id = ids(obs_id_num == uobs_id(t)); % tracking for the specific obs_code
                            trk_code = this.obs_code(id(1),:);
                            
                            res = this.value(:, id) * scale;
                            
                            fh = figure('Visible', 'off'); fh.Name = sprintf('%03d: %s Res %s', fh.Number, marker_name, trk_code); fh.NumberTitle = 'off';
                            Core_UI.beautifyFig(fh); drawnow;
                            
                            fh_list = [fh_list; fh]; %#ok<AGROW>
                            fig_name = sprintf('Res_Per_Sat_%s_%s_%s_%s_%s', trk_code(2:end), marker_name, cc.getSysName(sys_c), trk_code, this.time.first.toString('yyyymmdd_HHMM'));
                            fh.UserData = struct('fig_name', fig_name);
                            
                            ax2 = subplot(1, 24, 19:24);
                            ax1 = subplot(1, 24, 1:16);
                            
                            data_found = false;
                            figure(fh); % get focus;
                            for s = 1 : numel(id)
                                id_ok = find(~isnan(zero2nan(res(:,s))));
                                if any(id_ok)
                                    data_found = true;
                                    [~, id_sort] = sort(abs(res(id_ok, s)));
                                    scatter(ax1, id_ok(id_sort),  this.prn(id(s)) * ones(size(id_ok)), 80, (res(id_ok(id_sort), s)), 'filled');
                                    hold(ax1, 'on');
                                    err = std(zero2nan(res(:,s)), 'omitnan');
                                    if  verLessThan('matlab', '9.4')
                                        plot(ax2, mean(zero2nan(res(:,s)), 'omitnan') + [-err err], this.prn(id(s)) * [1 1], '.-', 'MarkerSize', 15, 'LineWidth', 3, 'Color', [0.6 0.6 0.6]);
                                        plot(ax2, mean(zero2nan(res(:,s)), 'omitnan'), this.prn(id(s)), '.', 'MarkerSize', 30, 'Color', [0.6 0.6 0.6]);
                                    else
                                        errorbar(ax2, mean(zero2nan(res(:,s)), 'omitnan'), this.prn(id(s)), err, '.', 'horizontal', 'MarkerSize', 30, 'LineWidth', 3, 'Color', [0.6 0.6 0.6]);
                                    end
                                    hold(ax2, 'on');
                                end
                            end
                            
                            if ~data_found
                                close(fh)
                                log = Core.getLogger;
                                log.addWarning(sprintf('No %s %s found in %s for constellation %s', name, trk_code, marker_name, cc.getSysName(sys_c)));
                            else
                                cax = caxis(ax1); caxis(ax1, [-1 1] * max(abs(cax)));
                                colormap(Cmap.get('PuOr', 2^11));
                                if min(abs(cax)) > 5
                                    setColorMap('PuOr', caxis(), 0.90, [-5 5])
                                end
                                cb = colorbar(ax1); cb.UserData = title(cb, iif(scale == 1e2, '[cm]', '[mm]')); ax1.Color = [0.9 0.9 0.9];
                                prn_ss = unique(cc.prn(cc.system == sys_c));
                                xlim(ax1, [1 size(res,1)]);
                                ylim(ax1, [min(prn_ss) - 1 max(prn_ss) + 1]);
                                h = ylabel(ax1, 'PRN'); h.FontWeight = 'bold';
                                ax1.YTick = prn_ss;
                                grid(ax1, 'on');
                                h = xlabel(ax1, 'epoch'); h.FontWeight = 'bold';
                                h = title(ax1, sprintf('%s %s %s\\fontsize{5} \n', cc.getSysName(sys_c), strrep(marker_name, '_', '\_'), trk_code(2:end)), 'interpreter', 'tex'); h.FontWeight = 'bold';
                                
                                ylim(ax2, [min(prn_ss) - 1 max(prn_ss) + 1]);
                                xlim(ax2, [-1 1] * (max(max(abs(mean(zero2nan(res(:,:)), 'omitnan'))), ...
                                    max(std(zero2nan(res(:,:)), 'omitnan'))) + 1));
                                ax2.YTick = prn_ss; ax2.Color = [1 1 1];
                                grid(ax2, 'on');
                                xlabel(ax2, sprintf('mean %s', iif(scale == 1e2, 'cm', 'mm')));
                                h = title(ax2, sprintf('mean\\fontsize{5} \n'), 'interpreter', 'tex'); h.FontWeight = 'bold';
                                linkaxes([ax1, ax2], 'y');
                                
                                Core_UI.beautifyFig(fh, 'dark');
                                Core_UI.addBeautifyMenu(fh);
                                fh.Visible = 'on'; drawnow;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % =========================================================================
    %%  PRIVATE
    % =========================================================================
    methods (Access = private)
        function init(this, type, time, value, prn, obs_code, rec_coo)
            % Init the residual object with new residuals (destroy the prevous content
            %
            % INPUT
            %   type        % 0,1,2,3 see RES_TYPE
            %   time        % time as GPS_Time                        GPS_Time [1 x 1] stores n_epoch
            %   pr          % matrix of pseudorange residuals
            %   ph          % matrix of carrier-phase residuals
            %
            %   value       % matrix of residuals
            %   obs_code    % type of tracking of the column (e.g. GL1C, GL1CL2WI, ...)
            %
            %   rec_coo     % <optional> Coordinates of the receiver
            %
            % SYNTAX
            %   this.import(type, time, pr, ph, prn, obs_code, rec_coo)
            
            this.type = type;
            this.time = time;
            this.value = value;
            this.prn = prn;
            this.obs_code = obs_code;
            this.rec_coo = rec_coo;
            
            % Remove entry with no obs_code
            code_ko = Constellation_Collector.obsCode2num(this.obs_code, this.prn);
            this.remEntry(code_ko == 0);
        end
    end
    
    % =========================================================================
    %%  STATIC
    % =========================================================================
    methods (Static)
        
    end
end
