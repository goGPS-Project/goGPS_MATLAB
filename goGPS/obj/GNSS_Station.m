%   CLASS GNSS_Station
% =========================================================================
%
% 
%   Class to store receiver data (observations, and characteristics)
%
% EXAMPLE
%   trg = GNSS_Station();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Receiver

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea, Giulio Tagliaferro ...
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
%--------------------------------------------------------------------------
classdef GNSS_Station < handle
    properties (SetAccess = public, GetAccess = public)
        marker_name    % marker name
        marker_type    % marker type
        number         % receiver number
        type           % receiver type
        version        % receiver version
        observer       % name of observer
        agency         % name of agency
        
        % ANTENNA ----------------------------------
        
        ant            % antenna number
        ant_type       % antenna type
        ant_delta_h    % antenna height from the ground [m]
        ant_delta_en   % antenna east/north offset from the ground [m]
        
        static         % static or dynamic receiver 1: static 0: dynamic
        
        
        work % handle to receiver Work Space
        out % handle to receiver outputs
        
        state
        log
        w_bar
    end
    
     % ==================================================================================================================================================
    %% PROPERTIES PLOTS
    % ==================================================================================================================================================
    
    properties
        slant_filter_win = 0;
    end
    % ==================================================================================================================================================
    %% METHODS INIT - CLEAN - RESET - REM -IMPORT
    % ==================================================================================================================================================
    methods
        function this = GNSS_Station()
            this.init();
        end
        
        function init(this)
            this.log = Logger.getInstance();
            this.state = Global_Configuration.getCurrentSettings();
            this.w_bar = Go_Wait_Bar.getInstance();
            this.reset();
        end
        
        
        function reset(this)
            this.marker_name  = 'unknown';  % marker name
            this.marker_type  = '';       % marker type
            this.number   = '000';
            this.type     = 'unknown';
            this.version  = '000';
        end
    end
     % ==================================================================================================================================================
    %% METHODS GETTER - TIME
    % ==================================================================================================================================================
    
    methods
        % standard utility
        function toString(this)
            % Display on screen information about the receiver
            % SYNTAX this.toString();
           this.work.toString();
           this.out.toString();
        end
        
        function req_rec = get(rec_list, marker_name)
            % Get the receivers with a certain Marker name
            % case unsensitive
            %
            % SYNTAX
            %   req_rec = rec_list.get(marker_name) 
            req_rec = [];
            for r = 1 : size(rec_list,2)
                rec = rec_list(~rec_list(:,r).isEmpty_mr ,r);
                if strcmpi(rec(1).getMarkerName, marker_name)
                    req_rec = [req_rec rec_list(:,r)]; %#ok<AGROW>
                end
            end                
        end
        
        function marker_name = getMarkerName(this)
            % Get the Marker name as specified in the RINEX file
            %
            % SYNTAX
            %   marker_name = getMarkerName(this)
            marker_name = this.marker_name;
            if isempty(marker_name)
                marker_name = this.getMarkerName4Ch();
            end
        end
        
        function marker_name = getMarkerName4Ch(this)
            % Get the Marker name as specified in the file name 
            % (first four characters)
            %
            % SYNTAX
            %   marker_name = getMarkerName(this)            
            marker_name = File_Name_Processor.getFileName(this.rinex_file_name);
            marker_name = marker_name(1 : min(4, length(marker_name)));
        end
        
        function out_prefix = getOutPrefix(this)
            % Get the name for exporting output (valid for dayly output)
            %   - marker name 4ch (from rinex file name)
            %   - 4 char year
            %   - 3 char doy
            %
            % SYNTAX
            [year, doy] = this.getCentralTime.getDOY();
            out_prefix = sprintf('%s_%04d_%03d_', this.getMarkerName4Ch, year, doy);
        end
        
        function is_empty = isEmpty_mr(this)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty();
            is_empty =  zeros(numel(this), 1);
            for r = 1 : numel(this)
                is_empty(r) =  this(r).out.isEmpty() && this(r).out.isEmpty();
            end
        end
        
        function is_empty = isEmptyOut_mr(this)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty();
            is_empty =  zeros(numel(this), 1);
            for r = 1 : numel(this)
                is_empty(r) =  this(r).out.isEmpty() && this(r).out.isEmpty();
            end
        end
        
        function n_epo = getNumEpochs(this)
            % Return the time span of the receiver
            % SYNTAX
            %   len = this.length();
            n_epo =  zeros(numel(this), 1);
            for r = 1 : numel(this)
                n_epo(r) =  this(r).work.time.length();
            end
            n_epo = sum(n_epo);
        end
        
        function n_sat = getMaxSat(this, sys_c)
            % get the number of satellites stored in the object
            %
            % SYNTAX 
            %   n_sat = getNumSat(<sys_c>)
            n_sat = zeros(size(this));
            
            for r = 1 : size(this, 2)
                rec = this(~this(:,r).isempty, r);
                
                if ~isempty(rec)
                    for s = 1 : size(rec, 1)
                        if nargin == 2
                            n_sat(s, r) = max(rec(s).work.go_id( (rec(s).work.system == sys_c)' & (rec(s).work.obs_code(:,1) == 'C' | rec(s).work.obs_code(:,1) == 'L') ));
                        else
                            n_sat(s, r) = max(rec(s).work.go_id(rec(s).work.obs_code(:,1) == 'C' | rec(s).work.obs_code(:,1) == 'L'));
                        end
                    end
                end
            end
            n_sat = max(n_sat, [], 1)';
        end
        
        function [time_lim_small, time_lim_large] = getWorkTimeSpan(this)
            % return a GPS_Time containing the first and last epoch stored in the Receiver
            %
            % OUTPUT
            %   time_lim_small     GPS_Time (first and last) epoch of the smaller interval
            %   time_lim_large     GPS_Time (first and last) epoch of the larger interval
            %
            % SYNTAX
            %   [time_lim_small, time_lim_large] = getTimeSpan(this);
            %
            time_lim_small = this(1).work.time.first;
            tmp_small = this(1).work.time.last;
            time_lim_large = time_lim_small.getCopy;
            tmp_large = tmp_small.getCopy;
            for r = 2 : numel(this)
                if time_lim_small < this(r).work.time.first
                    time_lim_small = this(r).work.time.first;
                end
                if time_lim_large > this(r).work.time.first
                    time_lim_large = this(r).work.time.first;
                end
                
                if tmp_small > this(r).work.time.last
                    tmp_small = this(r).work.time.last;
                end
                if tmp_large < this(r).work.time.last
                    tmp_large = this(r).work.time.last;
                end
            end
            time_lim_small.append(tmp_small);
            time_lim_large.append(tmp_large);
        end
        
        function [time_lim_small, time_lim_large] = getOutTimeSpan(this)
            % return a GPS_Time containing the first and last epoch stored in the Receiver
            %
            % OUTPUT
            %   time_lim_small     GPS_Time (first and last) epoch of the smaller interval
            %   time_lim_large     GPS_Time (first and last) epoch of the larger interval
            %
            % SYNTAX
            %   [time_lim_small, time_lim_large] = getTimeSpan(this);
            %
            time_lim_small = this(1).out.time.first;
            tmp_small = this(1).out.time.last;
            time_lim_large = time_lim_small.getCopy;
            tmp_large = tmp_small.getCopy;
            for r = 2 : numel(this)
                if time_lim_small < this(r).out.time.first
                    time_lim_small = this(r).out.time.first;
                end
                if time_lim_large > this(r).out.time.first
                    time_lim_large = this(r).out.time.first;
                end
                
                if tmp_small > this(r).out.time.last
                    tmp_small = this(r).out.time.last;
                end
                if tmp_large < this(r).out.time.last
                    tmp_large = this(r).out.time.last;
                end
            end
            time_lim_small.append(tmp_small);
            time_lim_large.append(tmp_large);
         end
        
        function enu = getPosENU_mr(this)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   enu     enu coordinates
            %
            % SYNTAX
            %   enu = this.getPosENU_mr()
            enu = this.getPosENU();
            n_rec = size(this, 2);
            n_sss = size(this, 1);
            enu = permute(reshape(enu', 3, n_sss, n_rec), [2 1 3]);
        end
    
        function xyz = getMedianPosXYZ_mr(this)
            % return the computed median position of the receiver
            %
            % OUTPUT
            %   xyz     geocentric coordinates
            %
            % SYNTAX
            %   xyz = this.getMedianPosXYZ_mr()
            
            xyz = [];
            for r = 1 : numel(this)
                if isempty(median(this(r).out.getPosXYZ(), 1))
                    xyz = [xyz; nan(1,3)]; %#ok<AGROW>
                else
                    xyz = [xyz; median(this(r).out.getPosXYZ(), 1)]; %#ok<AGROW>
                end
            end
        end
        
        function [lat, lon, h_ellips, h_ortho] = getMedianPosGeodetic_mr(this)
            % return the computed median position of the receiver
            % MultiRec: works on an array of receivers
            %
            % OUTPUT
            %   lat         latitude  [deg]
            %   lon         longitude [deg]
            %   h_ellips    ellipsoidical heigth [m]
            %   h_ortho     orthometric heigth [m]
            %
            % SYNTAX
            %   [lat, lon, h_ellips, h_ortho] = this.getMedianPosGeodetic();
            
            lat = nan(numel(this), 1);
            lon = nan(numel(this), 1);
            h_ellips = nan(numel(this), 1);
            h_ortho = nan(numel(this), 1);;
            for r = 1 : numel(this)
                xyz = this(r).out.xyz; %#ok<NASGU>
                try
                    xyz = median(this(r).xyz, 1);
                    [lat(r), lon(r), h_ellips(r)] = cart2geod(xyz);
                    if nargout == 4
                        gs = Global_Configuration.getInstance;
                        gs.initGeoid();
                        ondu = getOrthometricCorr(lat(r), lon(r), gs.getRefGeoid());
                        h_ortho(r) = h_ellips(r) - ondu; %#ok<AGROW>
                    end
                    lat(r) = lat(r) / pi * 180;
                    lon(r) = lon(r) / pi * 180;
                catch
                    % no position in the receiver
                end
            end
        end
             
        function getChalmersString(this)
            % get the string of the station to be used in http://holt.oso.chalmers.se/loading/
            % SYNTAX   this.getChalmersString();
            for r = 1 : size(this, 2)
                rec = this(~this(:,r).isempty, r);
                if ~isempty(rec)
                    xyz = rec.out.getMedianPosXYZ();
                    fprintf('%-24s %16.4f%16.4f%16.4f\n', rec(1).getMarkerName4Ch, xyz(1), xyz(2),xyz(3));
                end
            end
        end
        
        function [ztd, p_time, id_sync] = getZtd_mr(this)
            % MultiRec: works on an array of receivers
            % SYNTAX
            %  [ztd, p_time, id_sync] = this.getZtd_mr()
            [p_time, id_sync] = Receiver.getSyncTimeExpanded(this);
            
            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);
            
            n_rec = numel(this);
            ztd = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                ztd(~isnan(id_rec), r) = this(r).out.ztd(id_rec(~isnan(id_rec)));
            end
        end
        
        function [pwv, p_time, id_sync] = getPwv_mr(this)
            % MultiRec: works on an array of receivers
            % SYNTAX
            %  [pwv, p_time, id_sync] = this.getPwv_mr()
            [p_time, id_sync] = Receiver.getSyncTimeExpanded(this);
            
            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);
            
            n_rec = numel(this);
            pwv = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                pwv(~isnan(id_rec), r) = this(r).out.pwv(id_rec(~isnan(id_rec)));
            end
        end
        
        function [zwd, p_time] = getZwd_mr(this)
            % MultiRec: works on an array of receivers
            % SYNTAX
            %  [zwd, p_time, id_sync] = this.getZwd_mr()
            [p_time, id_sync] = Receiver.getSyncTimeTR(this);
            n_rec = numel(this);
            zwd = nan(size(id_sync{1}));
            for r = 1 : n_rec
                id_rec = id_sync{1}(:,r);
                zwd(~isnan(id_rec),r) = this(r).out.zwd(id_rec(~isnan(id_rec)));
            end
        end  
    end
     % ==================================================================================================================================================
    %% SETTER
    % ==================================================================================================================================================
    methods (Access = public)
        function setSlantFilterSize(this, win_size)
            % Setter multi_receiver to change filtering windows size for slant export
            % SYNTAX
            %   this.setSlantFilterSize(win_size)
            for r = 1:numel(this)
                this(r).slant_filter_win = win_size;
            end
    end
    end
     % ==================================================================================================================================================
    %% STATIC FUNCTIONS used as utilities
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function [p_time, id_sync] = getSyncTimeExpanded(rec, p_rate)
            % Get the common time among all the receivers
            %
            % SYNTAX
            %   [p_time, id_sync] = Receiver.getSyncTimeExpanded(rec, p_rate);
            %
            % EXAMPLE:
            %   [p_time, id_sync] = Receiver.getSyncTimeExpanded(rec, 30);
            
            if sum(~rec.isEmpty_mr) == 0
                % no valid receiver
                p_time = GPS_Time;
                id_sync = [];
            else
                if nargin == 1
                    p_rate = 1e-6;
                end
                
                % prepare reference time
                % processing time will start with the receiver with the last first epoch
                %          and it will stop  with the receiver with the first last epoch
                
                first_id_ok = find(~rec.isEmptyOut_mr, 1, 'first');
                if ~isempty(first_id_ok)
                    p_time_zero = round(rec(first_id_ok).out.time.first.getMatlabTime() * 24)/24; % get the reference time
                end
                
                % Get all the common epochs
                t = [];
                for r = 1 : numel(rec)
                    rec_rate = min(1, rec(r).out.time.getRate);
                    t = [t; round(rec(r).out.time.getRefTime(p_time_zero) * rec_rate) / rec_rate];
                    % p_rate = lcm(round(p_rate * 1e6), round(rec(r).out.time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
                end
                t = unique(t);
                
                % If p_rate is specified use it
                if nargin > 1
                    t = intersect(t, (0 : p_rate : t(end) + p_rate)');
                end
                
                % Create reference time
                p_time = GPS_Time(p_time_zero, t);
                id_sync = nan(p_time.length(), numel(rec));
                
                % Get intersected times
                for r = 1 : numel(rec)
                    rec_rate = min(1, rec(r).out.time.getRate);
                    [~, id1, id2] = intersect(t, round(rec(r).out.time.getRefTime(p_time_zero) * rec_rate) / rec_rate);
                    id_sync(id1 ,r) = id2;
                end
            end
        end
        
        function [pr, sigma_pr] = smoothCodeWithDoppler(pr, sigma_pr, pr_go_id, dp, sigma_dp, dp_go_id)
            % NOT WORKING due to iono
            % Smooth code with phase (aka estimate phase ambiguity and remove it)
            % Pay attention that the bias between phase and code is now eliminated
            % At the moment the sigma of the solution = sigma of phase
            %
            % SYNTAX
            %   [pr, sigma_ph] = smoothCodeWithDoppler(pr, sigma_pr, pr_go_id, ph, sigma_ph, ph_go_id, cs_mat)
            %
            % EXAMPLE
            %   [pr.obs, pr.sigma] = Receiver.smoothCodeWithDoppler(zero2nan(pr.obs), pr.sigma, pr.go_id, ...
            %                                                       zero2nan(dp.obs), dp.sigma, dp.go_id);
            %            

            for s = 1 : numel(dp_go_id)
                s_c = find(pr_go_id == dp_go_id(s));
                pr(isnan(dp(:,s)), s_c) = nan;
                
                lim = getOutliers(~isnan(pr(:,s)));
                for l = 1 : size(lim, 1)
                    id_arc = (lim(l,1) : lim(l,2))';
                    
                    len_a = length(id_arc);
                    A = [speye(len_a); [-speye(len_a - 1) sparse(len_a - 1, 1)] + [sparse(len_a - 1, 1) speye(len_a - 1)]];
                    Q = speye(2 * len_a - 1);
                    Q = spdiags([ones(len_a,1) * sigma_pr(s_c).^2; ones(len_a-1,1) * (2 * sigma_dp(s).^2)], 0, Q);
                    Tn = A'/Q;
                    data_tmp = [pr(id_arc, s_c); dp(id_arc, s)];
                    pr(id_arc, s) = (Tn*A)\Tn * data_tmp;
                end
            end
        end
        
        function [ph, sigma_ph] = smoothCodeWithPhase(pr, sigma_pr, pr_go_id, ph, sigma_ph, ph_go_id, cs_mat)
            % NOT WORKING when iono is present
            % Smooth code with phase (aka estimate phase ambiguity and remove it)
            % Pay attention that the bias between phase and code is now eliminated
            % At the moment the sigma of the solution = sigma of phase
            %
            % SYNTAX
            %   [ph, sigma_ph] = smoothCodeWithPhase(pr, sigma_pr, pr_go_id, ph, sigma_ph, ph_go_id, cs_mat)
            %
            % EXAMPLE
            %   [ph.obs, ph.sigma] = Receiver.smoothCodeWithPhase(zero2nan(pr.obs), pr.sigma, pr.go_id, ...
            %                                                   zero2nan(ph.obs), ph.sigma, ph.go_id, ph.cycle_slip);
            %            

            for s = 1 : numel(ph_go_id)
                s_c = find(pr_go_id == ph_go_id(s));
                pr(isnan(ph(:,s)), s_c) = nan;
                
                lim = getOutliers(~isnan(ph(:,s)), cs_mat(:,s));
                for l = 1 : size(lim, 1)
                    id_arc = (lim(l,1) : lim(l,2))';
                    
                    len_a = length(id_arc);
                    A = [speye(len_a); [-speye(len_a - 1) sparse(len_a - 1, 1)] + [sparse(len_a - 1, 1) speye(len_a - 1)]];
                    Q = speye(2 * len_a - 1);
                    Q = spdiags([ones(len_a,1) * sigma_pr(s_c).^2; ones(len_a-1,1) * (2 * sigma_ph(s).^2)], 0, Q);
                    Tn = A'/Q;
                    data_tmp = [pr(id_arc, s_c); diff(ph(id_arc, s))];
                    ph(id_arc, s) = (Tn*A)\Tn * data_tmp;
                end
            end
        end
        
        function obs_num = obsCode2Num(obs_code)
            % Convert a 3 char name into a numeric value (float)
            % SYNTAX
            %   obs_num = obsCode2Num(obs_code);            
            obs_num = Core_Utils.code3Char2Num(obs_code(:,1:3));
        end
        
        function obs_code = obsNum2Code(obs_num)
            % Convert a numeric value (float) of an obs_code into a 3 char marker
            % SYNTAX
            %   obs_code = obsNum2Code(obs_num)
            obs_code = Core_Utils.num2Code3Char(obs_num);            
        end
        
        function marker_num = markerName2Num(marker_name)
            % Convert a 4 char name into a numeric value (float)
            % SYNTAX
            %   marker_num = markerName2Num(marker_name);
            marker_num = Core_Utils.code4Char2Num(marker_name(:,1:4));
        end
        
        function marker_name = markerNum2Name(marker_num)
            % Convert a numeric value (float) of a station into a 4 char marker
            % SYNTAX
            %   marker_name = markerNum2Name(marker_num)
            marker_name = Core_Utils.num2Code4Char(marker_num);
        end
        
        function [y0, pc, wl, ref] = prepareY0(trg, mst, lambda, pivot)
            % prepare y0 and pivot_correction arrays (phase only)
            % SYNTAX [y0, pc] = prepareY0(trg, mst, lambda, pivot)
            % WARNING: y0 contains also the pivot observations and must be reduced by the pivot corrections
            %          use composeY0 to do it
            y0 = [];
            wl = [];
            pc = [];
            i = 0;
            for t = 1 : trg.n_epo
                for f = 1 : trg.n_freq
                    sat_pr = trg.p_range(t,:,f) & mst.p_range(t,:,f);
                    sat_ph = trg.phase(t,:,f) & mst.phase(t,:,f);
                    sat = sat_pr & sat_ph;
                    pc_epo = (trg.phase(t, pivot(t), f) - mst.phase(t, pivot(t), f));
                    y0_epo = ((trg.phase(t, sat, f) - mst.phase(t, sat, f)));
                    ref = median((trg.phase(t, sat, f) - mst.phase(t, sat, f)));
                    wl_epo = lambda(sat, 1);
                    
                    idx = i + (1 : numel(y0_epo))';
                    y0(idx) = y0_epo;
                    pc(idx) = pc_epo;
                    wl(idx) = wl_epo;
                    i = idx(end);
                end
            end
        end
        
        function y0 = composeY0(y0, pc, wl)
            % SYNTAX y0 = composeY0(y0, pc, wl)
            y0 = serialize((y0 - pc) .* wl);
            y0(y0 == 0) = []; % remove pivots
        end
        
        function [p_time, id_sync] = getSyncTimeTR(rec, obs_type, p_rate)
            % Get the common (shortest) time among all the used receivers and the target(s)
            % For each target (obs_type == 0) produce a different cella arrya with the sync of the other receiver
            % e.g.  Reference receivers @ 1Hz, trg1 @1s trg2 @30s
            %       OUTPUT 1 sync @1Hz + 1 sync@30s
            %
            % SYNTAX
            %   [p_time, id_sync] = Receiver.getSyncTimeTR(rec, obs_type, <p_rate>);
            %
            % SEE ALSO:
            %   this.getSyncTimeExpanded
            %
            if nargin < 3
                p_rate = 1e-6;
            end
            if nargin < 2
                obs_type = ones(1, numel(rec));
                obs_type(find(~rec.isEmpty_mr, 1, 'last')) = 0;
            end
            
            % Do the target(s) as last
            [~, id] = sort(obs_type, 'descend');
            
            % prepare reference time
            % processing time will start with the receiver with the last first epoch
            %          and it will stop  with the receiver with the first last epoch
            
            first_id_ok = find(~rec.isEmpty_mr, 1, 'first');
            p_time_zero = round(rec(first_id_ok).time.first.getMatlabTime() * 24)/24; % get the reference time
            p_time_start = rec(first_id_ok).time.first.getRefTime(p_time_zero);
            p_time_stop = rec(first_id_ok).time.last.getRefTime(p_time_zero);
            p_rate = lcm(round(p_rate * 1e6), round(rec(first_id_ok).time.getRate * 1e6)) * 1e-6;
            
            p_time = GPS_Time(); % empty initialization
            
            i = 0;
            for r = id
                ref_t{r} = rec(r).time.getRefTime(p_time_zero);
                if obs_type(r) > 0 % if it's not a target
                    if ~rec(r).isempty
                        p_time_start = max(p_time_start,  round(rec(r).time.first.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
                        p_time_stop = min(p_time_stop,  round(rec(r).time.last.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
                        p_rate = lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)) * 1e-6;
                    end
                else
                    % It's a target
                    
                    % recompute the parameters for the ref_time estimation
                    % not that in principle I can have up to num_trg_rec ref_time
                    % in case of multiple targets the reference times should be independent
                    % so here I keep the temporary rt0 rt1 r_rate var
                    % instead of ref_time_start, ref_time_stop, ref_rate
                    pt0 = max(p_time_start, round(rec(r).time.first.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
                    pt1 = min(p_time_stop, round(rec(r).time.last.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
                    pr = lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)) * 1e-6;
                    pt0 = ceil(pt0 / pr) * pr;
                    pt1 = floor(pt1 / pr) * pr;
                    
                    % return one p_time for each target
                    i = i + 1;
                    p_time(i) = GPS_Time(p_time_zero, (pt0 : pr : pt1)); %#ok<SAGROW>
                    p_time(i).toUnixTime();
                    
                    id_sync{i} = nan(p_time(i).length, numel(id));
                    for rs = id % for each rec to sync
                        if ~rec(rs).isempty && ~(obs_type(rs) == 0 && (rs ~= r)) % if it's not another different target
                            [~, id_ref, id_rec] = intersect(rec(rs).time.getRefTime(p_time_zero), (pt0 : pr : pt1));
                            id_sync{i}(id_rec, rs) = id_ref;
                        end
                    end
                end
            end
        end
                
        function sync(rec, rate)
            % keep epochs at a certain rate for a certain constellation
            %
            % SYNTAX
            %   this.keep(rate, sys_list)
            if nargin > 1 && ~isempty(rate)
                [~, id_sync] = Receiver.getSyncTimeExpanded(rec, rate);
            else
                [~, id_sync] = Receiver.getSyncTimeExpanded(rec);
            end
            % Keep the epochs in common
            % starting from the first when all the receivers are available
            % ending whit the last when all the receivers are available
            id_sync((find(sum(isnan(id_sync), 2) == 0, 1, 'first') : find(sum(isnan(id_sync), 2) == 0, 1, 'last')), :);
            
            % keep only synced epochs
            for r = 1 : numel(rec)
                rec(r).keepEpochs(id_sync(~isnan(id_sync(:, r)), r));
            end
            
        end
        
        function [res_ph1, mean_res, var_res] = legacyGetResidualsPh1(res_bin_file_name)
            %res_code1_fix  = [];                      % double differences code residuals (fixed solution)
            %res_code2_fix  = [];                      % double differences code residuals (fixed solution)
            %res_phase1_fix = [];                      % phase differences phase residuals (fixed solution)
            %res_phase2_fix = [];                      % phase differences phase residuals (fixed solution)
            %res_code1_float  = [];                    % double differences code residuals (float solution)
            %res_code2_float  = [];                    % double differences code residuals (float solution)
            res_phase1_float = [];                     % phase differences phase residuals (float solution)
            %res_phase2_float = [];                    % phase differences phase residuals (float solution)
            %outliers_code1 = [];                      % code double difference outlier? (fixed solution)
            %outliers_code2 = [];                      % code double difference outlier? (fixed solution)
            %outliers_phase1 = [];                     % phase double difference outlier? (fixed solution)
            %outliers_phase2 = [];                     % phase double difference outlier? (fixed solution)
            % observations reading
            log = Logger.getInstance();
            log.addMessage(['Reading: ' File_Name_Processor.getFileName(res_bin_file_name)]);
            d = dir(res_bin_file_name);
            fid_sat = fopen(res_bin_file_name,'r+');       % file opening
            num_sat = fread(fid_sat, 1, 'int8');                            % read number of satellites
            num_bytes = d.bytes-1;                                          % file size (number of bytes)
            num_words = num_bytes / 8;                                      % file size (number of words)
            num_packs = num_words / (2*num_sat*6);                          % file size (number of packets)
            buf_sat = fread(fid_sat,num_words,'double');                    % file reading
            fclose(fid_sat);                                                % file closing
            %res_code1_fix    = [res_code1_fix    zeros(num_sat,num_packs)]; % observations concatenation
            %res_code2_fix    = [res_code2_fix    zeros(num_sat,num_packs)];
            %res_phase1_fix   = [res_phase1_fix   zeros(num_sat,num_packs)];
            %res_phase2_fix   = [res_phase2_fix   zeros(num_sat,num_packs)];
            %res_code1_float  = [res_code1_float  zeros(num_sat,num_packs)];
            %res_code2_float  = [res_code2_float  zeros(num_sat,num_packs)];
            res_phase1_float = [res_phase1_float zeros(num_sat,num_packs)];
            %res_phase2_float = [res_phase2_float zeros(num_sat,num_packs)];
            %outliers_code1   = [outliers_code1   zeros(num_sat,num_packs)];
            %outliers_code2   = [outliers_code2   zeros(num_sat,num_packs)];
            %outliers_phase1  = [outliers_phase1  zeros(num_sat,num_packs)];
            %outliers_phase2  = [outliers_phase2  zeros(num_sat,num_packs)];
            i = 0;                                                           % epoch counter
            for j = 0 : (2*num_sat*6) : num_words-1
                i = i+1;                                                     % epoch counter increase
                %res_code1_fix(:,i)    = buf_sat(j + [1:num_sat]);            % observations logging
                %res_code2_fix(:,i)    = buf_sat(j + [1*num_sat+1:2*num_sat]);
                %res_phase1_fix(:,i)   = buf_sat(j + [2*num_sat+1:3*num_sat]);
                %res_phase2_fix(:,i)   = buf_sat(j + [3*num_sat+1:4*num_sat]);
                %res_code1_float(:,i)  = buf_sat(j + [4*num_sat+1:5*num_sat]);
                %res_code2_float(:,i)  = buf_sat(j + [5*num_sat+1:6*num_sat]);
                res_phase1_float(:,i) = buf_sat(j + (6*num_sat+1:7*num_sat));
                %res_phase2_float(:,i) = buf_sat(j + [7*num_sat+1:8*num_sat]);
                %outliers_code1(:,i)   = buf_sat(j + [8*num_sat+1:9*num_sat]);
                %outliers_code2(:,i)   = buf_sat(j + [9*num_sat+1:10*num_sat]);
                %outliers_phase1(:,i)  = buf_sat(j + [10*num_sat+1:11*num_sat]);
                %outliers_phase2(:,i)  = buf_sat(j + [11*num_sat+1:12*num_sat]);
            end
            res_ph1 = res_phase1_float';
            mean_res = mean(mean(zero2nan(res_ph1'), 'omitnan'), 'omitnan');
            var_res = max(var(zero2nan(res_ph1'), 'omitnan'));
        end
    end
    
    %% METHODS PLOTTING FUNCTIONS
    % ==================================================================================================================================================
    
    % Various debug images
    % name variant:
    %   c cartesian
    %   s scatter
    %   p polar
    %   m mixed
    methods (Access = public)
        function showAll(this)
            this.work.showAll@Receiver_Commons();
        end
        
        function showPositionENU(this, one_plot)
            % Plot East North Up coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionENU();
            if nargin == 1
                one_plot = true;
            end
            
            for r = 1 : size(this,2)
                rec = this(~this(:,r).isempty, r);
                if ~isempty(rec)
                    this.out.showPositionENU(one_plot);
                end
            end
        end
        
        function showPositionXYZ(this, one_plot)
            % Plot X Y Z coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionXYZ();
            if nargin == 1
                one_plot = true;
            end
            
            for r = 1 : size(this,2)
                rec = this(~this(:,r).isempty, r);
                if ~isempty(rec)
                    this.out.showPositionXYZ(one_plot);
                end
            end
        end
        
        function showMap(this, new_fig)
            if nargin < 2
                new_fig = true;
            end
            if new_fig
                f = figure;
            else
                f = gcf;
                hold on;
            end
            maximizeFig(f);
            [lat, lon] = cart2geod(this.getMedianPosXYZ_mr());
            
            plot(lon(:)./pi*180, lat(:)./pi*180,'.w','MarkerSize', 30);
            hold on;
            plot(lon(:)./pi*180, lat(:)./pi*180,'.k','MarkerSize', 10);
            plot(lon(:)./pi*180, lat(:)./pi*180,'ko','MarkerSize', 10, 'LineWidth', 2);
            
            if numel(this) == 1
                lon_lim = minMax(lon/pi*180);
                lat_lim = minMax(lat/pi*180);
                lon_lim(1) = lon_lim(1) - 0.05;
                lon_lim(2) = lon_lim(2) + 0.05;
                lat_lim(1) = lat_lim(1) - 0.05;
                lat_lim(2) = lat_lim(2) + 0.05;
            else
                lon_lim = xlim();
                lon_lim(1) = lon_lim(1) - diff(lon_lim)/3;
                lon_lim(2) = lon_lim(2) + diff(lon_lim)/3;
                lat_lim = ylim();
                lat_lim(1) = lat_lim(1) - diff(lat_lim)/3;
                lat_lim(2) = lat_lim(2) + diff(lat_lim)/3;
            end
            
            xlim(lon_lim);
            ylim(lat_lim);
            
            for r = 1 : numel(this)
                name = upper(this(r).getMarkerName());
                t = text(lon(r)./pi*180, lat(r)./pi*180, [' ' name ' '], ...
                    'FontWeight', 'bold', 'FontSize', 10, 'Color', [0 0 0], ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment','left');
                t.Units = 'pixels';
                t.Position(1) = t.Position(1) + 10 + 10 * double(numel(this) == 1);
                t.Units = 'data';
            end
            
            plot_google_map('alpha', 0.95, 'MapType', 'satellite');
            title('Receiver position');
            xlabel('Longitude [deg]');
            ylabel('Latitude [deg]');
        end % --> ask andrea
        
        function showDt(this)
            % Plot Clock error
            %
            % SYNTAX 
            %   this.plotDt
            
            for r = 1 : size(this, 2)
                rec = this(~this(:,r).isempty, r);
                if ~isempty(rec)
                    rec.work.showDt();
                end                
            end
        end
        
        function showZtdSlant(this, time_start, time_stop)
            %if isempty(this(1).ztd) || ~any(this(1).sat.slant_td(:))
            %    this(1).log.addWarning('ZTD and/or slants have not been computed');
            %else
            for r = 1 : size(this, 2)
                rec = this(~this(:,r).isempty, r);
                if isempty(rec)
                    this(1).log.addWarning('ZTD and/or slants have not been computed');
                else                    
                    this.out.showZtdSlant(time_start, time_stop);
                end
            end
            %end
        end
        
        function showTropoPar(this, par_name, new_fig)
            % one function to rule them all
            rec_ok = false(size(this,2), 1);
            for r = 1 : size(this, 2)
                switch lower(par_name)
                    case 'ztd'
                        rec_ok(r) = any(~isnan(this(:,r).out.getZtd));
                    case 'zwd'
                        rec_ok(r) = any(~isnan(this(:,r).out.getZwd));
                    case 'pwv'
                        rec_ok(r) = any(~isnan(this(:,r).out.getPwv));
                    case 'zhd'
                        rec_ok(r) = any(~isnan(this(:,r).out.getAprZhd));
                end               
            end
            rec_list = this(:, rec_ok);
            if numel(rec_list) == 0
                this(1).log.addError('No valid troposphere is present in the receiver list');
            else
                
                if nargin < 3
                    new_fig = true;
                end
                
                switch lower(par_name)
                    case 'ztd'
                        [tropo, t] = rec_list.out.getZtd();
                    case 'zwd'
                        [tropo, t] = rec_list.out.getZwd();
                    case 'pwv'
                        [tropo, t] = rec_list.out.getPwv();
                    case 'zhd'
                        [tropo, t] = rec_list.out.getAprZhd();
                end
                
                if ~iscell(tropo)
                    tropo = {tropo};
                    t = {t};
                end
                if isempty(tropo)
                    rec_list(1).out.log.addWarning([par_name ' and slants have not been computed']);
                else
                    if new_fig
                        f = figure; f.Name = sprintf('%03d: %s %s', f.Number, par_name, rec_list(1).out.cc.sys_c); f.NumberTitle = 'off';
                        old_legend = {};
                    else
                        l = legend;
                        old_legend = get(l,'String');
                    end
                    for r = 1 : size(rec_list, 2)
                        rec = rec_list(~rec_list(:,r).out.isempty, r);
                        if ~isempty(rec)
                            switch lower(par_name)
                                case 'ztd'
                                    [tropo, t] = rec.out.getZtd();
                                case 'zwd'
                                    [tropo, t] = rec.out.getZwd();
                                case 'pwv'
                                    [tropo, t] = rec.out.getPwv();
                                case 'zhd'
                                    [tropo, t] = rec.out.getAprZhd();
                            end
                            if new_fig
                                plot(t.getMatlabTime(), zero2nan(tropo'), '.', 'LineWidth', 4, 'Color', Core_UI.getColor(r, size(rec_list, 2))); hold on;
                            else
                                plot(t.getMatlabTime(), zero2nan(tropo'), '.', 'LineWidth', 4); hold on;
                            end
                            outm{r} = rec(1).getMarkerName();
                        end
                    end
                    
                    outm = [old_legend, outm];
                    [~, icons] = legend(outm, 'Location', 'NorthEastOutside', 'interpreter', 'none');
                    n_entry = numel(outm);
                    icons = icons(n_entry + 2 : 2 : end);
                    
                    for i = 1 : numel(icons)
                        icons(i).MarkerSize = 16;
                    end
                    
                    %ylim(yl);
                    %xlim(t(time_start) + [0 win_size-1] ./ 86400);
                    setTimeTicks(4,'dd/mm/yyyy HH:MMPM');
                    h = ylabel([par_name ' [m]']); h.FontWeight = 'bold';
                    grid on;
                    h = title(['Receiver ' par_name]); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                end
            end
        end  % --> ask andrea
        
        function showZhd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showTropoPar('ZHD', new_fig)
        end
        
        function showZwd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showTropoPar('ZWD', new_fig)
        end
        
        function showZtd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showTropoPar('ZTD', new_fig)
        end
        
        function showPwv(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showTropoPar('PWV', new_fig)
        end
        
        function showMedianTropoPar(this, par_name, new_fig)
            % one function to rule them all            
            rec_ok = false(size(this,2), 1);
            for r = 1 : size(this, 2)
                rec_ok(r) = any(~isnan(this(:,r).out.getZtd)); 
            end
            rec_list = this(:, rec_ok);
            
            if nargin < 3
                new_fig = true;
            end
            
            switch lower(par_name)
                case 'ztd'
                    [tropo] = rec_list(1).out.getZtd();
                case 'zwd'
                    [tropo] = rec_list(1).out.getZwd();
                case 'pwv'
                    [tropo] = rec_list(1).out.getPwv();
                case 'zhd'
                    [tropo] = rec_list(1).out.getAprZhd();
            end
                    
            if ~iscell(tropo)
                tropo = {tropo};
            end
            if isempty(tropo)
                rec_list(1).out.log.addWarning([par_name ' and slants have not been computed']);
            else
                if new_fig
                    f = figure; f.Name = sprintf('%03d: Median %s %s', f.Number, par_name, rec_list(1).out.cc.sys_c); f.NumberTitle = 'off';
                    old_legend = {};
                else
                    l = legend;
                    old_legend = get(l,'String');
                end
                for r = 1 : size(rec_list, 2)
                    rec = rec_list(~rec_list(:,r).isempty, r);
                    if ~isempty(rec)
                        switch lower(par_name)
                            case 'ztd'
                                [tropo] = rec.out.getZtd();
                            case 'zwd'
                                [tropo] = rec.out.getZwd();
                            case 'pwv'
                                [tropo] = rec.out.getPwv();
                            case 'zhd'
                                [tropo] = rec.out.getAprZhd();
                        end
                        [~, ~, ~, h_o] = rec(1).out.getPosGeodetic();
                        if new_fig
                            plot(h_o, median(tropo,'omitnan'), '.', 'MarkerSize', 25, 'LineWidth', 4, 'Color', Core_UI.getColor(r, size(rec_list, 2))); hold on;
                        else
                            plot(h_o, median(tropo,'omitnan'), '.', 'MarkerSize', 25, 'LineWidth', 4); hold on;
                        end
                        outm{r} = rec(1).getMarkerName();
                    end
                end
                
                outm = [old_legend, outm];
                [~, icons] = legend(outm, 'Location', 'NorthEastOutside', 'interpreter', 'none');
                n_entry = numel(outm);
                icons = icons(n_entry + 2 : 2 : end);
                
                for i = 1 : numel(icons)
                    icons(i).MarkerSize = 16;
                end
                
                %ylim(yl);
                %xlim(t(time_start) + [0 win_size-1] ./ 86400);
                h = ylabel([par_name ' [m]']); h.FontWeight = 'bold';
                h = xlabel('Elevation [m]'); h.FontWeight = 'bold';
                grid on;
                h = title(['Median Receiver ' par_name]); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
            end
        end

        function showMedianZhd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showMedianTropoPar('ZHD', new_fig)
        end
        
        function showMedianZwd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showMedianTropoPar('ZWD', new_fig)
        end
        
        function showMedianZtd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showMedianTropoPar('ZTD', new_fig)
        end
        
        function showMedianPwv(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showMedianTropoPar('PWV', new_fig)
        end
                
        function showZtdSlantRes_p(this, time_start, time_stop)
            for r = 1 : size(this, 2)
                if isempty(this(r).out.ztd) || ~any(this(r).out.sat.slant_td(:))
                    this.log.addWarning('ZTD and slants have not been computed');
                else
                    this.out.showZtdSlantRes_p(time_start, time_stop)
                end
            end
        end
    end 
end