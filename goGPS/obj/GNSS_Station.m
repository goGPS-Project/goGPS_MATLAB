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
        
        cc
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
        function this = GNSS_Station(cc, static)
            if nargin <= 1 || isempty(cc)
                cc = Constellation_Collector('G');
            end
            this.cc = cc;
            this.work = Receiver_Work_Space(cc, this);
            this.out = Receiver_Output(this.cc, this);
            if nargin >= 2 && ~isempty(static)
                this.static = logical(static);
            end
            this.init();
            this.reset();
        end
        
        function importRinexLegacy(this,rinex_file_name)
            if ~isempty(rinex_file_name) && (exist(rinex_file_name, 'file') == 2)
                this.work.rinex_file_name = rinex_file_name;                
            else
                this.work.rinex_file_name = '';
            end 
            this.work.load();
        end
        
        function importRinexes(this, rin_list, time_start, time_stop)
            % selct the files to be imported
            %
            % SYNTAX
            % this.importRinexes(rin_list, time_start, time_stop)
            rin_list.keepFiles(time_start, time_stop);
            this.work.importRinexFileList(rin_list, time_start, time_stop);

        end
        
        function init(this)
            this.log = Logger.getInstance();
            this.state = Global_Configuration.getCurrentSettings();
            this.w_bar = Go_Wait_Bar.getInstance();
            %this.work.resetWorkSpace();
            this.work = Receiver_Work_Space(this.cc, this);
        end
        
        function reset(this)
            this.marker_name  = 'unknown';  % marker name
            this.marker_type  = '';       % marker type
            this.number   = '000';
            this.type     = 'unknown';
            this.version  = '000';
        end
        
        function resetOut(this)
            this.out = Receiver_Output(this.cc, this);
        end
    end
     % ==================================================================================================================================================
    %% METHODS GETTER - TIME
    % ==================================================================================================================================================
    
    methods
        % standard utility
        function toString(sta_list)
            % Display on screen information about the receiver
            % SYNTAX this.toString();
            for i = 1:length(sta_list)
                if ~sta_list(i).work.isEmpty()
                    sta_list(i).work.toString();
                end
                if ~sta_list(i).out.isEmpty()
                    sta_list(i).out.toString();
                end
            end
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
            if ~isempty(this.work.rinex_file_name)
                marker_name = File_Name_Processor.getFileName(this.work.rinex_file_name);
            else
                marker_name = this.marker_name;
            end
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
                is_empty(r) =  this(r).work.isEmpty() && this(r).out.isEmpty();
            end
        end
        
        function is_empty = isEmpty(this)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty();
            
            is_empty =  this.work.isEmpty() && this.out.isEmpty();

        end
        
        function is_empty = isEmptyOut_mr(this)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty();
            is_empty =  zeros(numel(this), 1);
            for r = 1 : numel(this)
                is_empty(r) =  this(r).out.isEmpty();
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
        
        function n_sat = getMaxSat(sta_list, sys_c)
            % get the number of satellites stored in the object
            %
            % SYNTAX 
            %   n_sat = getNumSat(<sys_c>)
                  
            n_sat = zeros(numel(sta_list),1);
            for r = 1 : size(sta_list, 2)
                rec(r) = sta_list(r);
                if nargin == 2
                    n_sat(r) = rec.work.getMaxSat(sys_c);
                elseif nargin == 1
                    n_sat(r) = rec.work.getMaxSat();
                end
            end
        end
        
        function n_sat = getMaxNumSat(sta_list, sys_c)
            % get the number of maximum theoretical satellites stored in the object
            %
            % SYNTAX 
            %   n_sat = getMaxNumSat(<sys_c>)
                  
            n_sat = zeros(numel(sta_list),1);
            for r = 1 : size(sta_list, 2)
                rec(r) = sta_list(r);
                if nargin == 2
                    n_sat(r) = rec.work.getMaxNumSat(sys_c);
                elseif nargin == 1
                    n_sat(r) = rec.work.getMaxNumSat();
                end
            end
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
        
        function [tropo, time] = getTropoPar(sta_list, par_name)
            % get a tropo parameter among 'ztd', 'zwd', 'pwv', 'zhd'
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getAprZhd()
            
            tropo = {};
            time = {};
            for r = 1 : length(sta_list)
                time{r} = sta_list(r).out.getTime();
                switch lower(par_name)
                    case 'ztd'
                        [tropo{r}] = sta_list(r).out.getZtd();
                    case 'zwd'
                        [tropo{r}] = sta_list(r).out.getZwd();
                        if isempty(tropo{r}) || all(isnan(zero2nan(tropo{r})))
                            [tropo{r}] = sta_list(r).out.getAprZwd();
                        end
                    case 'pwv'
                        [tropo{r}] = sta_list(r).out.getPwv();
                    case 'zhd'
                        [tropo{r}] = sta_list(r).out.getAprZhd();
                end
            end
            
            if numel(tropo) == 1
                tropo = tropo{1};
                time = time{1};
            end
        end
        
        % wrappers to getTropoPar()
        function [tropo, time] = getZtd(sta_list)
            sta_list.getTropoPar('ztd');
        end
        
        function [tropo, time] = getZwd(sta_list)
            sta_list.getTropoPar('zwd');
        end
        
        function [tropo, time] = getPwv(sta_list)
            sta_list.getTropoPar('ztd');
        end
        
        function [tropo, time] = getAprZhd(sta_list)
            sta_list.getTropoPar('zhd');
        end
        
        function rec_works = getWork(sta_list, id)
            % return the working receiver for a gnss station array
            %
            % SYNTAX
            %  rec_works = sta_list.getWork(<id>)
            if nargin < 2
                rec_works= [sta_list.work];
            else
                id(id > numel(sta_list)) = [];
                rec_works= [sta_list(id).work];
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
                
                first_id_ok = find(~rec.isEmpty_mr, 1, 'first');
                if ~isempty(first_id_ok)
                    p_time_zero = round(rec(first_id_ok).time.first.getMatlabTime() * 24)/24; % get the reference time
                end
                
                % Get all the common epochs
                t = [];
                for r = 1 : numel(rec)
                    rec_rate = min(1, rec(r).time.getRate);
                    t = [t; round(rec(r).time.getRefTime(p_time_zero) * rec_rate) / rec_rate];
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
                    rec_rate = min(1, rec(r).time.getRate);
                    [~, id1, id2] = intersect(t, round(rec(r).time.getRefTime(p_time_zero) * rec_rate) / rec_rate);
                    id_sync(id1 ,r) = id2;
                end
            end
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
                    if ~rec(r).isEmpty
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
                        if ~rec(rs).isEmpty && ~(obs_type(rs) == 0 && (rs ~= r)) % if it's not another different target
                            [~, id_ref, id_rec] = intersect(rec(rs).time.getRefTime(p_time_zero), (pt0 : pr : pt1));
                            id_sync{i}(id_rec, rs) = id_ref;
                        end
                    end
                end
            end
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
        function showAll(sta_list)
            for i = 1:numel(sta_list)
                sta_list(i).out.showAll;
            end
        end
        
        function showPositionENU(this, one_plot)
            % Plot East North Up coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionENU();
            if nargin == 1
                one_plot = false;
            end
            
            for r = 1 : length(this)
                rec = this(r).out;
                if ~rec.isEmpty()
                    rec.showPositionENU(one_plot);
                end
            end
        end
        
        function showPositionXYZ(sta_list, one_plot)
            % Plot X Y Z coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionXYZ();
            if nargin == 1
                one_plot = false;
            end
            
            for r = 1 : size(sta_list,2)
                rec = sta_list(~sta_list(r).isEmpty, r);
                if ~isempty(rec)
                    sta_list.out.showPositionXYZ(one_plot);
                end
            end
        end
        
        function showMap(sta_list, new_fig)
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
            [lat, lon] = cart2geod(sta_list.getMedianPosXYZ_mr());
            
            plot(lon(:)./pi*180, lat(:)./pi*180,'.w','MarkerSize', 30);
            hold on;
            plot(lon(:)./pi*180, lat(:)./pi*180,'.k','MarkerSize', 10);
            plot(lon(:)./pi*180, lat(:)./pi*180,'ko','MarkerSize', 10, 'LineWidth', 2);
            
            if numel(sta_list) == 1
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
            
            for r = 1 : numel(sta_list)
                name = upper(sta_list(r).getMarkerName());
                t = text(lon(r)./pi*180, lat(r)./pi*180, [' ' name ' '], ...
                    'FontWeight', 'bold', 'FontSize', 10, 'Color', [0 0 0], ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment','left');
                t.Units = 'pixels';
                t.Position(1) = t.Position(1) + 10 + 10 * double(numel(sta_list) == 1);
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
                rec = this(r);
                if ~isempty(rec)
                    rec.out.showDt();
                end                
            end
        end
        
        function showZtdSlant(sta_list, time_start, time_stop)
            for r = 1 : size(sta_list, 2)
                rec = sta_list(~sta_list(r).isEmpty, r);
                if isempty(rec)
                    rec.log.addWarning('ZTD and/or slants have not been computed');
                else
                    if nargin < 3
                        rec.out.showZtdSlant();
                    else
                        rec.out.showZtdSlant(time_start, time_stop);
                    end
                end
            end
        end
        
        function showTropoPar(sta_list, par_name, new_fig)
            % one function to rule them all
            
            [tropo, t] = sta_list.getTropoPar(par_name);
            if ~iscell(tropo)
                tropo = {tropo};
                t = {t};
            end
            
            rec_ok = false(numel(sta_list), 1);
            for r = 1 : size(sta_list, 2)
                rec_ok(r) = ~isempty(tropo{r});
            end
            
            sta_list = sta_list(rec_ok);
            tropo = tropo(rec_ok);
            t = t(rec_ok);
            
            if numel(sta_list) == 0
                sta_list(1).log.addError('No valid troposphere is present in the receiver list');
            else
                if nargin < 3
                    new_fig = true;
                end                                
                
                if isempty(tropo)
                    sta_list(1).out.log.addWarning([par_name ' and slants have not been computed']);
                else
                    if new_fig
                        f = figure; f.Name = sprintf('%03d: %s %s', f.Number, par_name, sta_list(1).out.cc.sys_c); f.NumberTitle = 'off';
                        old_legend = {};
                    else
                        l = legend;
                        old_legend = get(l,'String');
                    end
                    for r = 1 : numel(sta_list)
                        rec = sta_list(r);
                        if new_fig
                            plot(t{r}.getMatlabTime(), zero2nan(tropo{r}'), '.', 'LineWidth', 4, 'Color', Core_UI.getColor(r, size(sta_list, 2))); hold on;
                        else
                            plot(t{r}.getMatlabTime(), zero2nan(tropo{r}'), '.', 'LineWidth', 4); hold on;
                        end
                        outm{r} = rec(1).getMarkerName();
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
        end
        
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
