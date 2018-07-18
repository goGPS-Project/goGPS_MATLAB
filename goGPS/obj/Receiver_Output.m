%   CLASS Receiver_Output
% =========================================================================
%
% 
%   Class to store receiver outputs
%
% EXAMPLE
%   trg = Receiver_Output();
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
classdef Receiver_Output < Receiver_Commons
    % ==================================================================================================================================================
    %% PROPERTIES POSITION
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        time_pos % time of the positions
    end
    
    %% PROPERTIES CELESTIAL INFORMATIONS
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        sat = struct( ...
            'outlier_idx_ph',   [], ...    % logical index of outliers
            'cycle_slip_idx_ph',[], ...    % logical index of cycle slips
            'quality',          [], ...    % quality
            'az',               [], ...    % double  [n_epoch x n_sat] azimuth
            'el',               [], ...    % double  [n_epoch x n_sat] elevation
            'res',              [], ...    % residual per staellite
            'mfw',              [], ...    % mapping funvtion wet
            'mfh',              []  ...    % mapping funvtion hysdrostatic
            )
    end
    % ==================================================================================================================================================
    %% PROPERTIES TROPO
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        
        
        pressure      % pressure           double   [n_epoch x 1]
        temperature   % temperature           double   [n_epoch x 1]
        humidity      % humidity           double   [n_epoch x 1]
    end
    % ==================================================================================================================================================
    %% METHODS INIT - CLEAN - RESET - REM -IMPORT
    % ==================================================================================================================================================
    
    methods
        function this = Receiver_Output(cc, parent)
            if nargin < 2
                cc = Constellation_Collector('G');
            end
            this.cc = cc;
            this.parent = parent;
            this.init();
        end
        
        function init(this)
            this.log = Logger.getInstance();
            this.state = Global_Configuration.getCurrentSettings();
            this.rf = Core_Reference_Frame.getInstance();
            this.w_bar = Go_Wait_Bar.getInstance();
            
            this.reset();
        end
        
        function reset(this)
            this.reset@Receiver_Commons();
            
            this.sat = struct(  ...
                'outlier_idx_ph',   [], ...    % logical index of outliers
                'cycle_slip_idx_ph',[], ...    % logical index of cycle slips
                'quality',          [], ...    % quality
                'az',               [], ...    % double  [n_epoch x n_sat] azimuth
                'el',               [], ...    % double  [n_epoch x n_sat] elevation
                'res',              [], ...    % residual per staellite
                'mfw',              [], ...    % mapping funvtion wet
                'mfh',              []  ...    % mapping funvtion hysdrostatic
                );
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
            for r = 1 : numel(this)
                if ~this(r).isEmpty
                    fprintf('----------------------------------------------------------------------------------\n')
                    this(r).log.addMarkedMessage(sprintf('Receiver (%d) %s', r, this(r).parent.getMarkerName));
                    fprintf('----------------------------------------------------------------------------------\n')
                    this(r).log.addMessage(sprintf(' From     %s', this(r).time.first.toString()));
                    this(r).log.addMessage(sprintf(' to       %s', this(r).time.last.toString()));
                    this(r).log.newLine();
                    this(r).log.addMessage(sprintf(' Rate of the observations [s]:            %d', this(r).getRate()));
                    this(r).log.newLine();
                    
                    fprintf(' ----------------------------------------------------------\n')
                    if ~isempty(this(r).xyz)
                        enu = zero2nan(this(r).xyz); [enu(:, 1), enu(:, 2), enu(:, 3)] = cart2plan(zero2nan(this(r).xyz(:,1)), zero2nan(this(r).xyz(:,2)), zero2nan(this(r).xyz(:,3)));
                        xyz_m = median(zero2nan(this(r).xyz), 1, 'omitnan');
                        enu_m = median(enu, 1, 'omitnan');
                        this(r).log.newLine();
                        this(r).log.addMessage(' Receiver median position:');
                        this(r).log.addMessage(sprintf('     X = %+16.4f m        E = %+16.4f m\n     Y = %+16.4f m        N = %+16.4f m\n     Z = %+16.4f m        U = %+16.4f m', ...
                            xyz_m(1), enu_m(1), xyz_m(2), enu_m(2), xyz_m(3), enu_m(3)));
                    end
                    fprintf(' ----------------------------------------------------------\n')
                end
            end
        end
        
        function time = getTime(this)
            % return the time stored in the object
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getTime()
            
            time = this(1).time.getCopy();
        end
        
        function desync = getDesync(this)
            desync = this.desync;
        end        
        
        function dt_pp = getDtPrePro(this)
            dt_pp = this.dt_ip;
        end
        
        function time = getPositionTime(this)
            % return the time of the computed positions
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getPositionTime()
            time = this.time_pos.getCopy();
        end
        
        function [pwv, time] = getPwv(this)
            % SYNTAX
            %  [pwv, time] = this.getPwv()
            
            pwv = {};
            time = {};
            for r = 1 : size(this, 2)
                time{r} = this(1, r).time.getEpoch(this(1, r).getIdSync); %#ok<AGROW>
                pwv{r} = this(1, r).pwv(this(1, r).getIdSync); %#ok<AGROW>
                
                for s = 2 : size(this, 1)
                    pwv_tmp = this(s, r).pwv(this(s, r).getIdSync);
                    time_tmp = this(s, r).time.getEpoch(this(s, r).getIdSync);
                    pwv{r} = [pwv{r}; pwv_tmp];
                    time{r} = time{r}.append(time_tmp);
                end
            end
            
            if numel(pwv) == 1
                pwv = pwv{1};
                time = time{1};
            end
        end
        
        function [quality, az, el] = getQuality(this)
            % SYNTAX
            %  [quality, az, el] = this.getQuality()
            quality = this.sat.quality(this.getIdSync,:);
            az = this.sat.az(this.getIdSync,:);
            el = this.sat.el(this.getIdSync,:);
        end
        
        function missing_epochs = getMissingEpochs(this)
            % return a logical array of missing (code) epochs
            %
            % SYNTAX
            %   missing_epochs = this.getMissingEpochs()
            %
            missing_epochs = true(this.time.length,1);
        end
        
        function id_sync = getIdSync(this)
            id_sync = 1 : this.time.length;
        end
        
        function [mfh, mfw] = getSlantMF(this)
            mfh = this.sat.mfh;
            mfw = this.sat.mfw;
        end
        
        function [day_lim] = getDayLim(this)
            % get the start index and end index of the days presents in the object
            %
            % SYNTAX
            %  [day_lim] = this.getDayLim()
            mat_time = this.time.getMatlabTime();
            day = floor(mat_time);
            days = unique(day)';
            n_days = length(days);
            day_lim = zeros(n_days,2);
            for i = 1:n_days
                d = days(i);
                day_lim(i,1) = find(mat_time > d,1,'first');
                day_lim(i,2) = find(mat_time < (d+1),1,'last');
            end
        end
        
    end
    
    % ==================================================================================================================================================
    %% METHODS IMPORT / EXPORT
    % ==================================================================================================================================================
    
    methods
        function injectResult(this, rec_work)
            %inject the results of receiver work into receiver output
            %
            % SYNTAX
            %  this.injectResult(rec_work)
            
            % set the id_sync only to time in between out times
            basic_export = false;
            id_sync_old = rec_work.getIdSync();
            if isempty(id_sync_old)
                rec_work.id_sync = 1 : rec_work.time.length;
                basic_export = true;
            end
            is_last_session = rec_work.time.last >= this.state.sss_date_stop;
            rec_work.cropIdSync4out(true, ~this.state.isSmoothTropoOut() || is_last_session);
            
            work_time = rec_work.getTime();
            initial_len = this.time.length;
            is_this_empty = this.time.isempty || ~any(this.ztd);
            if is_this_empty
                idx1 = 1;
                idx2 = 0;
                this.time = work_time;
            else
                if this.state.isSmoothTropoOut()
                    smt_buf_rgt = this.time.getEpoch(numel(this.ztd));
                    smt_buf_lft = rec_work.time.first();
                    idx_smt1 = this.time.getEpoch(1 : numel(this.ztd)) >= smt_buf_lft;                    
                    idx_smt2 = rec_work.time.getEpoch(id_sync_old) <= smt_buf_rgt;
                    time_1 = this.time.getEpoch(idx_smt1);
                    time_2 = rec_work.time.getEpoch(id_sync_old(idx_smt2));
                end
                [this.time, idx1, idx2] = this.time.injectBatch(work_time);
            end
            %%% inject data
            if ~basic_export                
                % Inject times
                this.dt      = Core_Utils.injectData(this.dt, rec_work.getDt(), idx1, idx2);
                this.desync  = Core_Utils.injectData(this.desync, rec_work.getDesync(), idx1, idx2);
                this.dt_ip   = Core_Utils.injectData(this.dt_ip, rec_work.getDtIp(), idx1, idx2);
                this.apr_zhd = Core_Utils.injectData(this.apr_zhd, rec_work.getAprZhd(), idx1, idx2);
                this.apr_zwd = Core_Utils.injectData(this.apr_zwd, rec_work.getAprZwd(), idx1, idx2);
                
                % Inject used meteo parameters
                [p, t, h]         = rec_work.getPTH(true);
                this.pressure     = Core_Utils.injectData(this.pressure, p, idx1, idx2);
                this.temperature  = Core_Utils.injectData(this.temperature, t, idx1, idx2);
                this.humidity     = Core_Utils.injectData(this.humidity, h, idx1, idx2);
                
                % Inject mapping functions
                [mfh, mfw]            = rec_work.getSlantMF();
                this.sat.mfw          = Core_Utils.injectData(this.sat.mfw, mfw, idx1, idx2);
                this.sat.mfh          = Core_Utils.injectData(this.sat.mfh, mfh, idx1, idx2);
                                
                % Inject outliers and cs
                this.sat.outlier_idx_ph    = Core_Utils.injectData(this.sat.outlier_idx_ph, rec_work.getOOutPh(), idx1, idx2);
                this.sat.cycle_slip_idx_ph = Core_Utils.injectData(this.sat.cycle_slip_idx_ph, rec_work.getOCsPh(), idx1, idx2);
                
                if ~this.state.isSmoothTropoOut() || is_this_empty
                    % Inject tropo related parameters
                    this.ztd     = Core_Utils.injectData(this.ztd, rec_work.getZtd(), idx1, idx2);
                    this.zwd     = Core_Utils.injectData(this.zwd, rec_work.getZwd(), idx1, idx2);
                    this.pwv     = Core_Utils.injectData(this.pwv, rec_work.getPwv(), idx1, idx2);
                    [gn, ge]     = rec_work.getGradient();
                    this.tgn     = Core_Utils.injectData(this.tgn, gn, idx1, idx2);
                    this.tge     = Core_Utils.injectData(this.tge, ge, idx1, idx2);
                    this.sat.res = Core_Utils.injectData(this.sat.res, rec_work.getResidual(), idx1, idx2);
                else
                    % there is probblably smoothing
                    % save idx, they might be useful
                    bk_idx1 = idx1;
                    bk_idx2 = idx2;
                end
            end
            
            if (initial_len == size(this.sat.az,1)) && (idx2 == 0)
                [~, idx1, idx2] = this.time.injectBatch(work_time);                
            end
            [az, el] = rec_work.getAzEl;
            this.sat.az      = Core_Utils.injectData(this.sat.az, az, idx1, idx2);
            this.sat.el      = Core_Utils.injectData(this.sat.el, el, idx1, idx2);
            this.sat.quality = Core_Utils.injectData(this.sat.quality, rec_work.getQuality(), idx1, idx2);
            
            %%% single results
            if isempty(this.time_pos)
                idx1 = 1;
                idx2 = 0;
                this.time_pos = rec_work.getPositionTime();
                data_len  = rec_work.getPositionTime().length;
            else
                [this.time_pos, idx1, idx2] = this.time_pos.injectBatch(rec_work.getPositionTime());
                data_len  = rec_work.getPositionTime().length;
            end
            
            this.xyz      = Core_Utils.injectData(this.xyz, rec_work.getPosXYZ, idx1, idx2, [data_len, 3]);
            this.enu      = Core_Utils.injectData(this.enu, rec_work.getPosENU, idx1, idx2, [data_len, 3]);
            
            [lat, lon, h_ellips, h_ortho] = rec_work.getPosGeodetic();
            
            this.lat      = Core_Utils.injectData(this.lat, lat, idx1, idx2, [data_len, 1]);
            this.lon      = Core_Utils.injectData(this.lon, lon, idx1, idx2, [data_len, 1]);
            this.h_ellips = Core_Utils.injectData(this.h_ellips, h_ellips, idx1, idx2, [data_len, 1]);
            this.h_ortho  = Core_Utils.injectData(this.h_ortho, h_ortho, idx1, idx2, [data_len, 1]);
            
            this.s0_ip    = Core_Utils.injectData(this.s0_ip, rec_work.s0_ip, idx1, idx2, [data_len, 1]);
            this.s0       = Core_Utils.injectData(this.s0, rec_work.s0, idx1, idx2, [data_len, 1]);
            
            % reset the old  complete id_sync
            rec_work.id_sync = id_sync_old;
            % inject with smoothing
            if ~basic_export && ~is_this_empty && this.state.isSmoothTropoOut()
                rec_work.cropIdSync4out(false, is_last_session);
                idx_smt2 = idx_smt2(1 : numel(rec_work.getZtd));
                id_start     = find(time_1 >= rec_work.out_start_time, 1, 'first'); % The first id of the new session
                if ~isempty(id_start)
                    this.ztd     = Core_Utils.injectSmtData(this.ztd, rec_work.getZtd(), idx_smt1, idx_smt2, time_1, time_2, id_start);
                    this.zwd     = Core_Utils.injectSmtData(this.zwd, rec_work.getZwd(), idx_smt1, idx_smt2, time_1, time_2, id_start);
                    this.pwv     = Core_Utils.injectSmtData(this.pwv, rec_work.getPwv(), idx_smt1, idx_smt2, time_1, time_2, id_start);
                    [gn, ge]     = rec_work.getGradient();
                    this.tgn     = Core_Utils.injectSmtData(this.tgn, gn, idx_smt1, idx_smt2, time_1, time_2, id_start);
                    this.tge     = Core_Utils.injectSmtData(this.tge, ge, idx_smt1, idx_smt2, time_1, time_2, id_start);
                    res = nan(size(this.ztd,1),size(this.sat.res,2));
                    res_in = rec_work.getResidual();
                    for i = 1 : size(this.sat.res,2)
                        res(:,i)   = Core_Utils.injectSmtData(this.sat.res(:,i), res_in(:,i), idx_smt1, idx_smt2, time_1, time_2, id_start);
                    end
                    this.sat.res = res;
                else                    
                    % Inject tropo related parameters
                    tmp = rec_work.getZtd();
                    this.ztd     = [this.ztd; tmp(~idx_smt2)];
                    tmp = rec_work.getZwd();
                    this.zwd     = [this.zwd; tmp(~idx_smt2)];
                    tmp = rec_work.getPwv();
                    this.pwv     = [this.pwv; tmp(~idx_smt2)];
                    [gn, ge]     = rec_work.getGradient();
                    this.tgn     = [this.tgn; gn(~idx_smt2)];
                    this.tge     = [this.tge; ge(~idx_smt2)];
                    res_in = rec_work.getResidual();
                    this.sat.res = [this.sat.res; res_in(~idx_smt2,:)];
                end
                rec_work.id_sync = id_sync_old; % restore id_sync_old
            end
            
            %--- append additional coo
            if this.state.flag_coo_rate
                this.add_coo.coo =  [this.add_coo.coo ; rec_work.add_coo.coo];
                if isempty(this.add_coo.time)
                    this.add_coo.time = rec_work.add_coo.time.getCopy();
                else
                    this.add_coo.time.append(rec_work.add_coo.time);
                end
                this.add_coo.rate =  [this.add_coo.rate ; rec_work.add_coo.rate];
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
        
        function showAll(this)
            this.toString
            this.showAll@Receiver_Commons();
            this.showDt();
        end
        
        function showDt(this)
            % Plot Clock error
            %
            % SYNTAX
            %   this.plotDt
            
            rec = this;
            if ~isempty(rec)
                f = figure; f.Name = sprintf('%03d: Dt Err', f.Number); f.NumberTitle = 'off';
                t = rec.time.getMatlabTime();
                plot(t, rec.getDesync, '-k', 'LineWidth', 2);
                hold on;
                plot(t, (rec.getDtIP), '-', 'LineWidth', 2);
                if any(rec.getDt)
                    plot(t, rec.getDt, '-', 'LineWidth', 2);
                    plot(t, rec.getTotalDt, '-', 'LineWidth', 2);
                    legend('desync time', 'dt correction from LS on Code', 'residual dt from carrier phases', 'total dt', 'Location', 'NorthEastOutside');
                else
                    legend('desync time', 'dt correction from LS on Code', 'Location', 'NorthEastOutside');
                end
                xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('receiver clock error [s]'); h.FontWeight = 'bold';
                h = title(sprintf('dt - receiver %s', rec.parent.getMarkerName),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
            end
        end   
        
        function showOutliersAndCycleSlip(this, sys_c_list)
            % Plot the outliers found
            % SYNTAX this.showOutliers()
            
            if nargin == 1
                sys_c_list = this.cc.getAvailableSys;
            end
            f = figure; f.Name = sprintf('%03d: CS, Outlier', f.Number); f.NumberTitle = 'off';
            for sys_c = sys_c_list
                ss_id = find(this.cc.sys_c == sys_c);
                switch numel(sys_c_list)
                    case 2
                        subplot(1,2, ss_id);
                    case 3
                        subplot(2,2, ss_id);
                    case 4
                        subplot(2,2, ss_id);
                    case 5
                        subplot(2,3, ss_id);
                    case 6
                        subplot(2,3, ss_id);
                    case 7
                        subplot(2,4, ss_id);
                end
                
                ep = repmat((1: this.time.length)',1, size(this.sat.outlier_idx_ph, 2));
                
                for prn = this.cc.prn(this.cc.system == sys_c)'
                    s = this.cc.getIndex(sys_c, prn);
                    cs = ep(this.sat.cycle_slip_idx_ph(:, s) ~= 0);
                    plot(cs,  prn * ones(size(cs)), '.k', 'MarkerSize', 20);
                    hold on;
                    out = ep(this.sat.outlier_idx_ph(:, s) ~= 0);
                    plot(out,  prn * ones(size(out)), '.', 'MarkerSize', 20, 'Color', [1 0.4 0]);
                end
                prn_ss = this.cc.prn(this.cc.system == sys_c);
                xlim([1 this.time.length]);
                ylim([min(prn_ss) - 1 max(prn_ss) + 1]);
                h = ylabel('PRN'); h.FontWeight = 'bold';
                ax = gca(); ax.YTick = prn_ss;
                grid on;
                h = xlabel('epoch'); h.FontWeight = 'bold';
                h = title(sprintf('%s %s cycle-slip(b) & outlier(o)', this.cc.getSysName(sys_c), this.parent.marker_name), 'interpreter', 'none'); h.FontWeight = 'bold';
            end
        end
        
            function showOutliersAndCycleSlip_p(this, sys_c_list)
            % Plot Signal to Noise Ration in a skyplot
            % SYNTAX this.plotSNR(sys_c)
            
            % SNRs
            if nargin == 1
                sys_c_list = this.cc.getAvailableSys;
            end
            
            for sys_c = sys_c_list
                f = figure; f.Name = sprintf('%03d: CS, Out %s', f.Number, sys_c); f.NumberTitle = 'off';
                polarScatter([],[],1,[]);
                hold on;
                
                for s = this.cc.getGoIds(sys_c)
                    az = this.sat.az(:,s);
                    el = this.sat.el(:,s);
                    
                    cs = sum(this.sat.cycle_slip_idx_ph(:, s), 2) > 0;
                    out = sum(this.sat.outlier_idx_ph(:, s), 2) > 0;
                    
                    decl_n = (serialize(90 - el(cs)) / 180*pi) / (pi/2);
                    x = sin(az(cs)/180*pi) .* decl_n; x(az(cs) == 0) = [];
                    y = cos(az(cs)/180*pi) .* decl_n; y(az(cs) == 0) = [];
                    plot(x, y, '.k', 'MarkerSize', 20)
                    
                    decl_n = (serialize(90 - el(out)) / 180*pi) / (pi/2);
                    x = sin(az(out)/180*pi) .* decl_n; x(az(out) == 0) = [];
                    y = cos(az(out)/180*pi) .* decl_n; y(az(out) == 0) = [];
                    plot(x, y, '.', 'MarkerSize', 20, 'Color', [1 0.4 0]);
                end
                h = title(sprintf('%s %s cycle-slip(b) & outlier(o)', this.cc.getSysName(sys_c), this.parent.marker_name), 'interpreter', 'none'); h.FontWeight = 'bold';
            end
            
        end
        
        
    end
end
