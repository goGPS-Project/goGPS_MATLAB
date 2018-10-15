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
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 4 - nightly
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
        old_work % handle to the old Work Space
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
        
        function importRinexLegacy(this,rinex_file_name, rate)
            if ~isempty(rinex_file_name) && (exist(rinex_file_name, 'file') == 2)
                this.work.rinex_file_name = rinex_file_name;                
            else
                this.work.rinex_file_name = '';
            end 
            this.work.load(rate);
            this.work.out_start_time = this.work.time.first;
            this.work.out_stop_time = this.work.time.last;
        end
        
        function importRinexes(this, rin_list, time_start, time_stop, rate)
            % select the files to be imported
            %
            % SYNTAX
            % this.importRinexes(rin_list, time_start, time_stop, rate)
            this.work.importRinexFileList(rin_list, time_start, time_stop, rate);
        end
        
        function init(this)
            this.log = Core.getLogger();
            this.state = Core.getState();
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
        
        function id = getStationId(rec_list, marker_name)
            % Given a marker_name get the sequencial id of a station
            %
            % SYNTAX
            %   id = findStationId(this, marker_name)
            marker4ch_list = '';
            for r = 1 : numel(rec_list)
                try
                    marker4ch_list(r, :) = char(rec_list(r).getMarkerName4Ch);
                catch
                    % the name is shorter or missing => ignore
                end
            end
            id = find(Core_Utils.code4Char2Num(upper(marker4ch_list)) == Core_Utils.code4Char2Num(upper(marker_name)));
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
            if this.out.length == 0
                time = this.work.time.getCopy;
            else
                time = this.out.time.getCopy;
            end
            [year, doy] = time.getCentralTime.getDOY();
            out_prefix = sprintf('%s_%04d_%03d_', this.getMarkerName4Ch, year, doy);
        end
        
        function is_empty = isEmpty_mr(sta_list)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty_mr();
            is_empty =  zeros(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                is_empty(r) =  sta_list(r).work.isEmpty() && sta_list(r).out.isEmpty();
            end
        end
        
        function is_empty = isOutEmpty_mr(sta_list)
            % Return if the object (out) does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isOutEmpty_mr();
            is_empty =  zeros(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                is_empty(r) =  sta_list(r).out.isEmpty();
            end
        end        
        
        function is_empty = isEmpty(this)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty();
            
            is_empty =  this.work.isEmpty() && this.out.isEmpty();

        end
        
        function time = getTime(this)
            % return the time stored in the object
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getTime()
            time = this.out.getTime();
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
            %   [time_lim_small, time_lim_large] = getWorkTimeSpan(this);
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
            %   [time_lim_small, time_lim_large] = getOutTimeSpan(this);
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
         
        function [rate] = getRate(this)
            % SYNTAX
            %   rate = this.getRate();
            rate = this.out.getTime.getRate;
            if isnan(rate)            
                rate = this.work.getTime.getRate;
            end
        end
        
        function coo = getPos(sta_list)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   coo     Coordinate object array
            %
            % SYNTAX
            %   coo = sta_list.getPos()
            for r = 1 : numel(sta_list)
                coo(r) = sta_list(r).out.getPos();
            end
        end
        
        function xyz = getPosXYZ(sta_list)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   xyz     XYZ coordinates
            %
            % SYNTAX
            %   xyz = this.getPosENU()
            for r = 1 : numel(sta_list)
                xyz{r} = sta_list(r).out.getPosXYZ();
            end
        end
        
        function enu = getPosENU(sta_list)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   enu     enu coordinates
            %
            % SYNTAX
            %   enu = this.getPosENU()
            for r = 1 : numel(sta_list)
                enu{r} = sta_list(r).out.getPosENU();
            end
        end
        
        function [enu, p_time] = getPosENU_mr(sta_list)
            % return the positions computed for n receivers
            % multi_rec mode (synced)
            %
            % OUTPUT
            %   enu     enu synced coordinates
            %
            % SYNTAX
            %   enu = sta_list.getPosENU_mr()
            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list, [], true);
            
            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);
            
            n_rec = numel(sta_list);
            enu = nan(size(id_sync, 1), 3, n_rec);
            for r = 1 : n_rec
                enu_rec = sta_list(r).out.getPosENU();
                id_rec = id_sync(:,r);
                enu(~isnan(id_rec), :, r) = enu_rec(id_rec(~isnan(id_rec)), :);
            end
        end
    
        function xyz = getMedianPosXYZ(this)
            % return the computed median position of the receiver
            %
            % OUTPUT
            %   xyz     geocentric coordinates
            %
            % SYNTAX
            %   xyz = this.getMedianPosXYZ()
            
            xyz = [];
            for r = 1 : numel(this)
                if isempty(median(this(r).out.getPosXYZ(), 1))
                    xyz = [xyz; nan(1,3)]; %#ok<AGROW>
                else
                    xyz = [xyz; median(this(r).out.getPosXYZ(), 1)]; %#ok<AGROW>
                end
            end
        end
        
        function [lat, lon, h_ellips, h_ortho] = getMedianPosGeodetic(sta_list)
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
            %   [lat, lon, h_ellips, h_ortho] = sta_list.getMedianPosGeodetic();
            
            lat = nan(numel(sta_list), 1);
            lon = nan(numel(sta_list), 1);
            h_ellips = nan(numel(sta_list), 1);
            h_ortho = nan(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                if sta_list(1).static
                    [lat(r), lon(r), h_ellips(r), h_ortho(r)] = sta_list(r).out.getMedianPosGeodetic;
                else
                    [lat{r}, lon{r}, h_ellips{r}, h_ortho{r}] = sta_list(r).out.getMedianPosGeodetic;
                end
            end
        end
             
        function getChalmersString(this)
            % get the string of the station to be used in http://holt.oso.chalmers.se/loading/
            % SYNTAX   this.getChalmersString();
            this(1).log.addMarkedMessage('Chalmers ocean loading computation must be required manually:');            
            this(1).log.addMessage(this(1).log.indent('go to http://holt.oso.chalmers.se/loading/ and request a BLQ file'));
            this(1).log.addMessage(this(1).log.indent('using ocean tide model FES2004'));
            this(1).log.addMessage(this(1).log.indent('select also to compensate the values for the motion'));            
            this(1).log.addMessage(this(1).log.indent('Use the following string for the staion locations:'));
            this(1).log.addMessage([char(8) '//------------------------------------------------------------------------']);

            for r = 1 : size(this, 2)
                rec = this(~this(:,r).isEmpty, r);
                if ~isempty(rec)
                    xyz = rec.out.getMedianPosXYZ();
                    if isempty(xyz)
                        xyz = rec.work.getMedianPosXYZ();
                    end                        
                    this(1).log.addMessage([char(8) sprintf('%-24s %16.4f%16.4f%16.4f', rec(1).getMarkerName4Ch, xyz(1), xyz(2),xyz(3))]);
                end
            end
            
            this(1).log.addMessage([char(8) '//------------------------------------------------------------------------']);
        end
        
        function [ztd, p_time, id_sync] = getZtd_mr(sta_list)
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [ztd, p_time, id_sync] = this.getZtd_mr()
            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);
            
            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);
            
            n_rec = numel(sta_list);
            ztd = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec = id_rec(~isnan(id_rec) & id_rec <= length(sta_list(r).out.ztd));
                ztd(~isnan(id_rec), r) = sta_list(r).out.ztd(~isnan(id_rec));
            end
        end
                
        function [ztd_res, p_time, ztd_height] = getReducedZtd_mr(sta_list)
            % MultiRec: works on an array of receivers
            % Reduce the ZTD of all the stations removing the 
            % component dependent with the altitude 
            %
            % SYNTAX
            %  [ztd_res, p_time, ztd_height] = sta_list.getReducedZtd_mr()
            
            med_ztd = median(sta_list.getZtd_mr, 'omitnan')';
            degree = 5;
            [~, ~, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;

            ztd_height = Core_Utils.interp1LS(h_o, med_ztd, degree);
            [ztd, p_time] = sta_list.getZtd_mr();
            ztd_res = bsxfun(@minus, ztd', ztd_height)';
        end
        
        function [pwv, p_time, id_sync] = getPwv_mr(sta_list)
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [pwv, p_time, id_sync] = sta_list.getPwv_mr()
            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);
            
            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);
            
            n_rec = numel(sta_list);
            pwv = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                pwv(~isnan(id_rec), r) = sta_list(r).out.pwv(id_rec(~isnan(id_rec)));
            end
        end
        
        function [zwd, p_time] = getZwd_mr(sta_list)
            % MultiRec: works on an array of receivers
            % SYNTAX
            %  [zwd, p_time, id_sync] = sta_list.getZwd_mr()
            [p_time, id_sync] = GNSS_Station.getSyncTimeTR(sta_list);
            n_rec = numel(sta_list);
            zwd = nan(size(id_sync{1}));
            for r = 1 : n_rec
                id_rec = id_sync{1}(:,r);
                zwd(~isnan(id_rec),r) = sta_list(r).out.zwd(id_rec(~isnan(id_rec)));
            end
        end  
        
        function [tropo, time] = getTropoPar(sta_list, par_name)
            % get a tropo parameter among 'ztd', 'zwd', 'pwv', 'zhd'
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getAprZhd()
            
            tropo = {};
            time = {};
            for r = 1 : numel(sta_list)
                time{r} = sta_list(r).out.getTime();
                switch lower(par_name)
                    case 'ztd'
                        [tropo{r}] = sta_list(r).out.getZtd();
                    case 'zwd'
                        [tropo{r}] = sta_list(r).out.getZwd();
                        if isempty(tropo{r}) || all(isnan(zero2nan(tropo{r})))
                            [tropo{r}] = sta_list(r).out.getAprZwd();
                        end
                      case 'gn'
                        [tropo{r}] = sta_list(r).out.getGradient();
                    case 'ge'
                        [~,tropo{r}] = sta_list(r).out.getGradient();
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
        
        function [tropo, time] = getZtd(sta_list)
            % Get ZTD
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getZtd()
            [tropo, time] = sta_list.getTropoPar('ztd');
        end
        
        function [tropo, time] = getZwd(sta_list)
            % Get ZWD
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getZwd()
            [tropo, time] = sta_list.getTropoPar('zwd');
        end
        
        function [tropo, time] = getPwv(sta_list)
            % Get PWV
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getPwv()
            [tropo, time] = sta_list.getTropoPar('pwv');
        end
        
        function [tropo, time] = getAprZhd(sta_list)
            % Get ZHD
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getZhd()
            [tropo, time] = sta_list.getTropoPar('zhd');
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
            for r = 1 : numel(this)
                this(r).slant_filter_win = win_size;
            end
        end
    end
    
    % ==================================================================================================================================================
    %% STATIC FUNCTIONS used as utilities
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function [p_time, id_sync] = getSyncTimeExpanded(rec, p_rate, use_pos_time)
            % Get the common time among all the receivers
            %
            % SYNTAX
            %   [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(rec, <p_rate>, <use_pos_time>);
            %
            % EXAMPLE:
            %   [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(rec, 30);
            
            if nargin < 3 || isempty(use_pos_time)
                use_pos_time = false;
            end
            
            if sum(~rec.isEmpty_mr) == 0
                % no valid receiver
                p_time = GPS_Time;
                id_sync = [];
            else
                if nargin < 2 || isempty(p_rate)
                    p_rate = 1e-6;
                    
                    for r = 1 : numel(rec)
                        if use_pos_time
                            p_rate = lcm(round(p_rate * 1e6), round(rec(r).out.time_pos.getRate * 1e6)) * 1e-6; % enable this line to sync rates
                        else
                            p_rate = lcm(round(p_rate * 1e6), round(rec(r).out.time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
                        end
                    end
                end
                
                % prepare reference time
                % processing time will start with the receiver with the last first epoch
                %          and it will stop  with the receiver with the first last epoch
                
                first_id_ok = find(~rec.isEmpty_mr, 1, 'first');
                if ~isempty(first_id_ok)
                    if use_pos_time
                        p_time_zero = round(rec(first_id_ok).out.time_pos.first.getMatlabTime() * 24)/24; % get the reference time for positions
                    else
                        p_time_zero = round(rec(first_id_ok).out.time.first.getMatlabTime() * 24)/24; % get the reference time
                    end
                end
                
                % Get all the common epochs
                t = [];
                for r = 1 : numel(rec)
                    if use_pos_time
                        rec_rate = min(86400, iif(rec(r).out.time_pos.length == 1, 86400, rec(r).out.time_pos.getRate));
                        t = [t; round(rec(r).out.time_pos.getRefTime(p_time_zero) / rec_rate) * rec_rate];
                    else
                        rec_rate = min(1, rec(r).out.time.getRate);
                        t = [t; round(rec(r).out.time.getRefTime(p_time_zero) / rec_rate) * rec_rate];
                    end
                    % p_rate = lcm(round(p_rate * 1e6), round(rec(r).out.time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
                end
                t = unique(t);
                
                % If p_rate is specified use it
                if nargin > 1
                    t = intersect(t, (t(1) : p_rate : t(end) + p_rate)');
                end
                
                % Create reference time
                p_time = GPS_Time(p_time_zero, t);
                id_sync = nan(p_time.length(), numel(rec));
                
                % Get intersected times
                for r = 1 : numel(rec)
                    if use_pos_time
                        rec_rate = iif(rec(r).out.time_pos.length == 1, 86400, rec(r).out.time_pos.getRate);
                        [~, id1, id2] = intersect(t, round(rec(r).out.time_pos.getRefTime(p_time_zero) / rec_rate) * rec_rate);
                    else
                        rec_rate = min(1, rec(r).out.time.getRate);
                        [~, id1, id2] = intersect(t, round(rec(r).out.time.getRefTime(p_time_zero) / rec_rate) * rec_rate);
                    end
                    
                    id_sync(id1, r) = id2;
                end
            end
        end
        
        function [p_time, id_sync] = getSyncTimeTR(sta_list, obs_type, p_rate)
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
                % choose the longest as reference
                len = zeros(1, numel(sta_list));
                for r = 1 : numel(sta_list)
                    len(r) = sta_list(r).out.length;
                end
                obs_type = ones(1, numel(sta_list));
                obs_type(find(len == max(len), 1, 'first')) = 0;
            end
            
            % Do the target(s) as last
            [~, id] = sort(obs_type, 'descend');
            
            % prepare reference time
            % processing time will start with the receiver with the last first epoch
            %          and it will stop  with the receiver with the first last epoch
            
            first_id_ok = find(~sta_list.isOutEmpty_mr, 1, 'first');
            p_time_zero = round(sta_list(first_id_ok).out.time.first.getMatlabTime() * 24)/24; % get the reference time
            p_time_start = sta_list(first_id_ok).out.time.first.getRefTime(p_time_zero);
            p_time_stop = sta_list(first_id_ok).out.time.last.getRefTime(p_time_zero);
            p_rate = lcm(round(p_rate * 1e6), round(sta_list(first_id_ok).out.time.getRate * 1e6)) * 1e-6;
            
            p_time = GPS_Time(); % empty initialization
            
            i = 0;
            for r = id
                ref_t{r} = sta_list(r).out.time.getRefTime(p_time_zero);
                if obs_type(r) > 0 % if it's not a target
                    if ~sta_list(r).out.isEmpty
                        p_time_start = max(p_time_start,  round(sta_list(r).out.time.first.getRefTime(p_time_zero) * sta_list(r).out.time.getRate) / sta_list(r).out.time.getRate);
                        p_time_stop = min(p_time_stop,  round(sta_list(r).out.time.last.getRefTime(p_time_zero) * sta_list(r).out.time.getRate) / sta_list(r).out.time.getRate);
                        p_rate = lcm(round(p_rate * 1e6), round(sta_list(r).out.time.getRate * 1e6)) * 1e-6;
                    end
                else
                    % It's a target
                    
                    % recompute the parameters for the ref_time estimation
                    % not that in principle I can have up to num_trg_rec ref_time
                    % in case of multiple targets the reference times should be independent
                    % so here I keep the temporary rt0 rt1 r_rate var
                    % instead of ref_time_start, ref_time_stop, ref_rate
                    pt0 = max(p_time_start, round(sta_list(r).out.time.first.getRefTime(p_time_zero) * sta_list(r).out.time.getRate) / sta_list(r).out.time.getRate);
                    pt1 = min(p_time_stop, round(sta_list(r).out.time.last.getRefTime(p_time_zero) * sta_list(r).out.time.getRate) / sta_list(r).out.time.getRate);
                    pr = lcm(round(p_rate * 1e6), round(sta_list(r).out.time.getRate * 1e6)) * 1e-6;
                    pt0 = ceil(pt0 / pr) * pr;
                    pt1 = floor(pt1 / pr) * pr;
                    
                    % return one p_time for each target
                    i = i + 1;
                    p_time(i) = GPS_Time(p_time_zero, (pt0 : pr : pt1));
                    p_time(i).toUnixTime();
                    
                    id_sync{i} = nan(p_time(i).length, numel(id));
                    for rs = id % for each rec to sync
                        if ~sta_list(rs).out.isEmpty && ~(obs_type(rs) == 0 && (rs ~= r)) % if it's not another different target
                            [~, id_ref, id_rec] = intersect(round(sta_list(rs).out.time.getRefTime(p_time_zero) * 1e5)/1e5, (pt0 : pr : pt1));
                            id_sync{i}(id_rec, rs) = id_ref;
                        end
                    end
                end
            end
        end
        
        function bsl_ids = getBaselineId(n_rec)
            % Get id of all the combinations of the stations
            %
            % SYNTAX:
            %   bsl_id = GNSS_Station.getBaselineId(n_rec);
            [r1, r2] = meshgrid(1 : n_rec, 1 : n_rec);
            bsl_ids = [serialize(tril(r1, -1)) serialize(tril(r2, -1))];
            bsl_ids = bsl_ids(bsl_ids(:, 1) > 0 & bsl_ids(:, 2) > 0, :);
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
        
        function showObsStats(sta_list)
            % Show statistics about the observations stored in the object
            %
            % SYNTAX
            %   this.showObsStats()
            
            for s = 1 : numel(sta_list)
                sta_list(s).work.showObsStats();
            end
        end
         
        function showPositionENU(sta_list, one_plot)
            % Plot East North Up coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionENU();
            if nargin == 1
                one_plot = false;
            end
            
            for r = 1 : length(sta_list)
                rec = sta_list(r).out;
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
            
            for r = 1 : length(sta_list)
                rec = sta_list(r).out;
                if ~isempty(rec)
                    rec.showPositionXYZ(one_plot);
                end
            end
        end
        
        function showPositionSigmas(sta_list, one_plot)
            % Show Sigmas of the solutions
            %
            % SYNTAX
            %   this.showPositionSigmas();
            
            if nargin == 1
                one_plot = false;
            end
            
            for r = 1 : length(sta_list)
                rec = sta_list(r).out;
                if ~isempty(rec)
                    rec.showPositionSigmas(one_plot);
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
            [lat, lon] = cart2geod(sta_list.getMedianPosXYZ());
            
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
        end % ==> ask Andrea? Andrea say: "what?"
        
        function showScatteredMap(lat, lon, s_time, data) % todo
            
        end
            
        function showMapZwd(this, new_fig, epoch)
            if nargin < 2 || isempty(new_fig)
                new_fig = true;
            end
            if new_fig
                f = figure;
            else
                f = gcf;
                hold on;
            end
            maximizeFig(f);
            
            [res_tropo, s_time] = this.getReducedZtd_mr();
            %[res_tropo, s_time] = this.getZwd_mr();
            res_tropo = res_tropo * 1e2;
            
            if nargin < 3
                epoch = 1 : s_time.length();
            end
                        
            coo = Coordinates.fromXYZ(this.getMedianPosXYZ);
            [lat, lon] = coo.getGeodetic;
            
            sh = scatter(lon(:)./pi*180, lat(:)./pi*180, 250, res_tropo(epoch(1),:)', 'filled');
            hold on;
            % plot(lon(:)./pi*180, lat(:)./pi*180,'ko','MarkerSize', 15, 'LineWidth', 2);
            %caxis([-10 10]);
            colormap(gat);
            colorbar;
            
            lon_lim = xlim();
            lon_lim(1) = lon_lim(1) - 0.1;
            lon_lim(2) = lon_lim(2) + 0.1;
            lat_lim = ylim();
            lat_lim(1) = lat_lim(1) - 0.1;
            lat_lim(2) = lat_lim(2) + 0.1;
            
            xlim(lon_lim);
            ylim(lat_lim);
                        
            if new_fig
                if FTP_Downloader.checkNet
                    plot_google_map('alpha', 0.95, 'MapType', 'satellite');
                end
                xlabel('Longitude [deg]');
                ylabel('Latitude [deg]');
            end
            th = title(sprintf('Receiver %s', s_time.getEpoch(1).toString('yyyy-mm-dd HH:MM:SS')), 'FontSize', 30);
            
            drawnow
            if numel(epoch) > 1                
                for e = 2 : 10 : numel(epoch)
                    if any(res_tropo(epoch(e),:))
                        th.String = sprintf('Receiver %s', s_time.getEpoch(epoch(e)).toString('yyyy-mm-dd HH:MM:SS'));
                        sh.CData = res_tropo(epoch(e),:)';
                        drawnow
                    end
                end
            end
        end
        
        function showDt(this)
            % Plot Clock error
            %
            % SYNTAX 
            %   sta_list.plotDt
            
            for r = 1 : size(this, 2)
                rec = this(r);
                if ~isempty(rec)
                    rec.out.showDt();
                end                
            end
        end
        
        function f_handle = showQuality_p(sta_list, type, flag_smooth)
            % Plot Signal to Noise Ration in a skyplot
            % SYNTAX f_handles = this.plotSNR(sys_c)
            
            % SNRs
            if nargin < 2
                type = 'snr';
            end
            
            f_handle = [];
            
            for r = 1 : numel(sta_list)
                if ~sta_list(r).out.isEmpty
                    rec = sta_list(r).out;
                    [quality, az, el] = rec.getQuality();
                else
                    rec = sta_list(r).work;
                    [quality, az, el] = rec.getQuality(type);
                end
                
                if nargin > 2 && flag_smooth
                    quality = Receiver_Commons.smoothSatData([],[],zero2nan(quality), [], 'spline', 900 / this.getRate, 10); % smoothing Quality => to be improved
                end
                
                f = figure; f.Name = sprintf('%03d: %s', f.Number, upper(type)); f.NumberTitle = 'off';
                f_handle(r) = f;
                id_ok = (~isnan(quality));
                polarScatter(serialize(az(id_ok)) / 180 * pi, serialize(90 - el(id_ok)) / 180 * pi, 45, serialize(quality(id_ok)), 'filled');
                colormap(jet);  cax = caxis();
                switch type
                    case 'snr'
                        caxis([min(cax(1), 10), max(cax(2), 55)]);
                        setColorMap([10 55], 0.9);
                end
                colorbar();
                h = title(sprintf('%s - receiver %s', upper(type), sta_list(r).getMarkerName4Ch()), 'interpreter', 'none');
                h.FontWeight = 'bold'; h.Units = 'pixels';
                h.Position(2) = h.Position(2) + 20; h.Units = 'data';
                
            end
        end
        
        
        function showResiduals(sta_list)
            % Plot Satellite Residuals
            %
            % SYNTAX 
            %   sta_list.showResiduals
            
            for r = 1 : size(sta_list, 2)
                rec = sta_list(r);
                if ~isempty(rec)
                    if ~rec.out.isEmpty
                        rec.out.showResiduals();
                    else
                        rec.work.showResiduals();
                    end
                end                
            end
        end
                
        function showZtdSlant(sta_list, time_start, time_stop)
            for r = 1 : size(sta_list, 2)
                rec = sta_list(~sta_list(r).isEmpty, r);
                if isempty(rec)
                    log = Core.getLogger();
                    log.addWarning('ZTD and/or slants have not been computed');
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
                log = Core.getLogger();
                log.addError('No valid troposphere is present in the receiver list');
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
                            plot(t{r}.getMatlabTime(), zero2nan(tropo{r}'), '.-', 'LineWidth', 2, 'Color', Core_UI.getColor(r, size(sta_list, 2))); hold on;
                        else
                            plot(t{r}.getMatlabTime(), zero2nan(tropo{r}'), '.-', 'LineWidth', 2); hold on;
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
        
         
        
        function showGn(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showTropoPar('GN', new_fig)
        end
        
        
        function showGe(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            this.showTropoPar('GE', new_fig)
        end
        
        function showZtdVsHeight(sta_list)
            % Show Median ZTD of n_receivers vs Hortometric height
            % 
            % SYNTAX
            %   sta_list.showZtdVsHeight();
            figure; 
            med_ztd = median(sta_list.getZtd_mr * 1e2, 'omitnan')';
            subplot(2,1,1); 
            [~, ~, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;
            plot(h_o, med_ztd, '.', 'MarkerSize', 20); hold on;
            ylabel('Median ZTD [cm]');
            xlabel('Elevation [m]');
            title('ZTD vs Height')

            degree = 5;
            y_out = Core_Utils.interp1LS(h_o, med_ztd, degree, h_o);
            plot(sort(h_o), Core_Utils.interp1LS(h_o, med_ztd, degree, sort(h_o)), '-', 'Color', Core_UI.COLOR_ORDER(3,:), 'LineWidth', 2);
            subplot(2,1,2); 
            plot(h_o, med_ztd - y_out, '.', 'MarkerSize', 20);

            ylabel('residual [cm]');
            xlabel('Elevation [m]');
            title('reduced ZTD vs Height')            
            
            sta_strange = find(abs(med_ztd - y_out) > 8);
            if ~isempty(sta_strange)
                sta_list.log.addMessage('Strange station detected');
                for s = 1 : numel(sta_strange)
                    sta_list.log.addMessage(sprintf(' - %s', sta_list.station_name(sta_list.sta_ok(sta_strange(s)),:)));
                end
            end
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
                for r = 1 : numel(rec_list)
                    rec = rec_list(~rec_list(r).isEmpty, r);
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
 
        function showBaselineENU(sta_list, baseline_ids, plot_relative_variation, one_plot)
            % Function to plot baseline between 2 or more stations
            %
            % INPUT:
            %   sta_list                 list of GNSS_Station objects
            %   baseline_ids             n_baseline x 2 - couple of id in sta_list to be used
            %   plot_relative_variation  show full baseline dimension / variation wrt the median value 
            %   one_plot                 use subplots (E, N, U) or a single plot
            %
            % SYNTAX
            %   showBaselineENU(sta_list, <baseline_ids = []>, <plot_relative_variation = true>, <one_plot = false>)
            %
            
            if (nargin < 4) || isempty(one_plot)
                one_plot = false;
            end
            if (nargin < 3) || isempty(plot_relative_variation)
                plot_relative_variation = true;
            end
            
            % remove empty receivers
            sta_list = sta_list(~sta_list.isEmpty_mr);
            
            if nargin < 2 || isempty(baseline_ids)
                n_rec = numel(sta_list);
                baseline_ids = GNSS_Station.getBaselineId(n_rec);
            end
            
            for b = 1 : size(baseline_ids, 1)
                rec = sta_list(baseline_ids(b, :));
                if ~isempty(rec(1)) && ~isempty(rec(2))
                    [enu, time] = rec.getPosENU_mr();
                    if size(enu, 1) > 1
                        rec(1).log.addMessage('Plotting positions');
                        
                        % prepare data
                        baseline = diff(enu, 1, 3);
                        if plot_relative_variation
                            baseline = bsxfun(@minus, baseline, median(baseline, 'omitnan')) * 1e3;
                        end
                        t = time.getMatlabTime();
                        
                        f = figure; f.Name = sprintf('%03d: BSL ENU %s - %s', f.Number, rec(1).getMarkerName4Ch, rec(2).getMarkerName4Ch); f.NumberTitle = 'off';
                        color_order = handle(gca).ColorOrder;
                        
                        if ~one_plot, subplot(3,1,1); end
                        plot(t, baseline(:, 1), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                        ax(3) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); 
                        if plot_relative_variation
                            h = ylabel('East [mm]'); h.FontWeight = 'bold';
                        else
                            h = ylabel('East [m]'); h.FontWeight = 'bold';
                        end
                        grid on;
                        h = title(sprintf('Baseline %s - %s', rec(1).getMarkerName4Ch, rec(2).getMarkerName4Ch),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                        if ~one_plot, subplot(3,1,2); end
                        plot(t, baseline(:, 2), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                        ax(2) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); 
                        if plot_relative_variation
                            h = ylabel('North [mm]'); h.FontWeight = 'bold';
                        else
                            h = ylabel('North [m]'); h.FontWeight = 'bold';
                        end

                        grid on;
                        if ~one_plot, subplot(3,1,3); end
                        plot(t, baseline(:,3), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                        ax(1) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); 
                        if plot_relative_variation
                            h = ylabel('Up [mm]'); h.FontWeight = 'bold';
                        else
                            h = ylabel('Up [m]'); h.FontWeight = 'bold';
                        end
                        
                        grid on;
                        if one_plot
                            if plot_relative_variation
                                h = ylabel('ENU [mm]'); h.FontWeight = 'bold';
                            else
                                h = ylabel('ENU [m]'); h.FontWeight = 'bold';
                            end                            
                            legend({'East', 'North', 'Up'}, 'Location', 'NorthEastOutside', 'interpreter', 'none');
                            
                        else
                            linkaxes(ax, 'x');
                        end
                        grid on;
                        
                    else
                        rec(1).log.addMessage('Plotting a single point static position is not yet supported');
                    end
                end
                
            end
        end
    end 
end
