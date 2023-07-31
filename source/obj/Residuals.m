classdef Residuals < Exportable
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 1.0
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
    %  Written by:        Andrea Gatti
    %  Contributors:      Andrea Gatti, ...
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
                    % uniform the size of obs_code
                    if size(res.obs_code,2) < size(this.obs_code,2)
                        n_el = size(this.obs_code,2);
                        res.obs_code = [res.obs_code, ' ' * char(ones(size(res.obs_code,1), n_el - size(res.obs_code,2), 'uint8'))];
                    end
                    if size(res.obs_code,2) > size(this.obs_code,2)
                        n_el = size(res.obs_code,2);
                        this.obs_code = [this.obs_code, ' ' * char(ones(size(this.obs_code,1), n_el - size(this.obs_code,2), 'uint8'))];
                    end
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
            
            if any(lid_ko)
                if ~islogical(lid_ko)
                    id_ko = lid_ko;
                    lid_ko = false(1,this.time.length);
                    lid_ko(id_ko) = true;
                end
                id_ok = ~lid_ko;
                this.time = this.time.getEpoch(id_ok);
                this.value(lid_ko, :) = [];
            end
        end
        
        function remEntry(this, lid_ko)
            % Remove an entry from the residuals
            %
            % INPUT
            %   lid_ko  logical array of ko entry
            %
            % SYNTAX
            %   this.remEntry(lid_ko);
            
            if any(lid_ko)
                if ~islogical(lid_ko)
                    id_ko = lid_ko;
                    lid_ko = false(1, size(this.value,2));
                    lid_ko(id_ko) = true;
                end
                this.value(:, lid_ko) = [];
                this.prn(lid_ko) = [];
                this.obs_code(lid_ko, :) = [];
            end
        end
        
        function cutEpochs(this, new_lim, end_lim)
            % Get the residual only in the time span given
            %
            % SYNTAX
            %   this.cutEpochs(new_limits)
            %   this.cutEpochs(lim_start, lim_stop)
            
            if nargin == 3
                new_lim = new_lim.getCopy;
                new_lim.append(end_lim);
            end
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
            % get an index that tell which residual are phase
            %
            % SYNTAX:
            %    [is_ph] = this.isPhase()
            if isempty(this.obs_code)
                is_ph = false(size(this.obs_code,1),1);
            else
                is_ph = this.obs_code(:,2) == 'L';
            end
        end
        
        function [is_co] = isCombined(this)
            % get an index that tell which residual are combined phases
            %
            % SYNTAX:
            %    [is_co] = this.isCombined()
            if size(this.obs_code,2) > 4
                is_co = this.obs_code(:,5) ~= ' ';
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
        function [az, el, sat_coo, sat_name, go_id] = getAzimuthElevation(this, id_ok)
            % Get azimuth and elevation of each satellite stored in residuals
            %
            %
            %   [az, el, sat_coo, sat_name, go_id] = this.getAzimuthElevation();
            
            core = Core.getCurrentCore;
            sky = core.sky;
            if isempty(core.state.eph_name)
                fw = File_Wizard();
                fw.conjureNavFiles(this.time.first, this.time.last);
            end
            if nargin == 2
                time = this.time.getEpoch(id_ok);
            else
                time = this.time;
            end
            lim = time.first.getCopy;
            lim.append(this.time.last);
            flag_no_clock = true;
            core.initSkySession(lim, flag_no_clock);
            cc = core.getConstellationCollector;
            go_id = unique(cc.getIndex(this.obs_code(:,1), this.prn));
            sat_name = cc.getSatName(go_id);
            [az, el, sat_coo] = sky.getAzimuthElevation(this.rec_coo, time, go_id);
        end
    end
    
    % =========================================================================
    %%  MULTIPATH
    % =========================================================================
    methods
        
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
                                    line = scatter(az(id_ok, id_sat(s)), el(id_ok, id_sat(s)), 25, serialize(res(id_ok, s)), 'filled');
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
                                    line = polarScatter(az(id_ok, id_sat(s))/180*pi, (90 -el(id_ok, id_sat(s)))/180*pi, 25, serialize(res(id_ok, s)), 'filled');
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
        
        function fh_list = showRes(res_list, marker_name, sys_c_list, is_ph)
            % Plot residuals of the solution
            %
            % SYNTAX
            %   fh_list = this.showRes(marker_name, sys_c_list, is_ph)
            log = Core.getLogger();
            fh_list = [];
            for this = res_list(:)'
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

                    if nargin < 2 || isempty(marker_name)
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
