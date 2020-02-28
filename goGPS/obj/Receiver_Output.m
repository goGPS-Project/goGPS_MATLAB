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
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti, Giulio Tagliaferro ...
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
        used_sys_c     % constellatioj used in the computataions
        selected_sys_c % active constellations in the computation
        sat = struct( ...
            'outliers',         [], ...    % logical index of outliers
            'cycle_slip',       [], ...    % logical index of cycle slips
            'quality',          [], ...    % quality
            'az',               [], ...    % double  [n_epoch x n_sat] azimuth
            'el',               [], ...    % double  [n_epoch x n_sat] elevation
            'mfw',              [], ...    % mapping function wet
            'mfh',              [], ...    % mapping function hysdrostatic
            'res',              [] ...    % processing residuals object
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
    %% METHODS MANIPULATION
    methods
        function applyRemAntennaOffset(sta_list, sgn)            
            % Add/Rem antennna offset stored in parent
            %
            % INPUT
            %   sgn:    +1 apply
            %           -1 remove
            % SYNTAX
            %    this.applyRemAntennaOffset(<sgn = 1>)
            
            if nargin == 1
                sgn = 1;
            end
            for r = 1 : numel(sta_list)
                coo_out = sta_list(r).getPos;
                coo_out.addOffset(-sgn .* [sta_list(r).parent.ant_delta_en sta_list(r).parent.ant_delta_h]);
                sta_list(r).xyz = coo_out.getXYZ;
            end
        end
    end
    
    %% METHODS INIT - CLEAN - RESET - REM -IMPORT
    % ==================================================================================================================================================
    
    methods
        function this = Receiver_Output(parent)
            this.parent = parent;
            this.initHandles();
            this.reset();
        end
        
        function reset(this)
            this.reset@Receiver_Commons();
            
            this.sat = struct(  ...
                'outliers',         [], ...    % logical index of outliers
                'cycle_slip',       [], ...    % logical index of cycle slips
                'quality',          [], ...    % quality
                'az',               [], ...    % double  [n_epoch x n_sat] azimuth
                'el',               [], ...    % double  [n_epoch x n_sat] elevation
                'res',              [], ...    % residual per satellite
                'res_pr_by_pr',     [], ...    % code residual per uncombined tracking
                'res_ph_by_ph',     [], ...    % phase residual per uncombined tracking
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
                    this(r).log.addMarkedMessage(sprintf('Receiver %s', this(r).parent.getMarkerName()));
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
        
        function time = getTimePositions(this)
            % return the time stored in the object
            % that correspond to position estimation epochs
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getTimePositions()
            
            time = this(1).time_pos.getCopy();
        end
        
        function [P,T,H] = getPTH(this)
            % get ztd
            %
            % SYNTAX
            %   ztd = this.getZtd()
            if max(this.getIdSync) > numel(this.ztd) || isempty(this.pressure) ||  isempty(this.temperature) ||  isempty(this.humidity)
                P = nan(size(this.getIdSync));
                T = nan(size(this.getIdSync));
                H = nan(size(this.getIdSync));
            else
                P = this.pressure(this.getIdSync);
                T = this.temperature(this.getIdSync);
                H = this.humidity(this.getIdSync);
            end
        end
        
        function desync = getDesync(this)
            desync = this.desync;
        end
        
        function dt_pp = getDtPrePro(this)
            dt_pp = this.dt_ip;
        end
        
        function n_sat = getMaxSat(sta_list, sys_c)
            % get the number of satellites stored in the object
            %
            % SYNTAX
            %   n_sat = getNumSat(<sys_c>)
            n_sat = zeros(size(sta_list));
            
            for r = 1 : size(sta_list, 2)
                rec = sta_list(~sta_list(r).isEmpty, r);
                if ~isempty(rec)
                    if nargin == 2
                        tmp = rec.getCC.getMaxNumSat(sys_c);
                        n_sat(r) = iif(isempty(tmp), 0, tmp);
                    else
                        tmp = rec.getCC.getMaxNumSat();
                        n_sat(r) = iif(isempty(tmp), 0, tmp);
                    end
                end
            end
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
                try
                    pwv{r} = this(1, r).pwv(this(1, r).getIdSync); %#ok<AGROW>
                catch 
                    pwv{r} = [];
                end
                
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
            try
                quality = this.sat.quality(this.getIdSync,:);
            catch
                quality = [];
            end
            
            try
                az = this.sat.az(this.getIdSync,:);
                el = this.sat.el(this.getIdSync,:);
            catch
                az = [];
                el = [];
            end
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
        
        function [mfh, mfw, cotan_term] = getSlantMF(this)
            mfh = this.sat.mfh;
            mfw = this.sat.mfw;
            if nargout == 3
                % at the moment this is needed for experimental Zernike estimation -> not yet implemented for out
                Core.getLogger.addWarning('at the moment Cotan_term in getSlantMF is needed for experimental Zernike estimation -> not yet implemented for out objects');
                cotan_term = nan; 
            end            
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
        
        function xyz = getIgsXYZ(this, mode)
            % Get the official IGS solution for either daily of weekly IGS
            % combinations
            %
            % SYNTAX:
            %    xyz = getIgsXYZ(this, mode)
            if nargin < 2
                mode = 'daily';
            end
            xyz = nan(size(this.xyz));
            if strcmpi(mode,'daily')
                data_dir = fullfile(Core.getInstallDir, '..' , 'data');
                fnp = File_Name_Processor();
                for e = 1: this.time_pos.length
                    c_time = this.time_pos.getEpoch(e);
                    filename = fnp.dateKeyRep(sprintf('%s/station/IGS_solutions/COO/${WWWW}/igs${YY}P${WWWWD}.ssc', data_dir), c_time);
                    if exist(filename, 'file') ~= 2
                        remote_file_name = fnp.dateKeyRep('gnss/products/${WWWW}/igs${YY}P${WWWWD}.ssc.Z',c_time);
                        ftp_dw = FTP_Downloader('cddis.nasa.gov',21);
                        [pathstr, name, ext] = fileparts(filename);
                        ftp_dw.downloadUncompress(remote_file_name, pathstr);
                    end
                    
                    % Multiplatform
                    if exist(filename, 'file') == 2
                        fid = fopen(filename, 'rb');
                        txt = fread(fid, '*char')';
                        fclose(fid);
                        cmdout = regexp(txt, sprintf('STAX   %s[^\n]*(?=\n)', upper(this.parent.getMarkerName4Ch)), 'match');
                        if ~isempty(cmdout)
                            x = sscanf(cmdout{2}(40:61),'%f');
                            cmdout = regexp(txt, sprintf('STAY   %s[^\n]*(?=\n)', upper(this.parent.getMarkerName4Ch)), 'match');
                            y = sscanf(cmdout{2}(40:61),'%f');
                            cmdout = regexp(txt, sprintf('STAZ   %s[^\n]*(?=\n)', upper(this.parent.getMarkerName4Ch)), 'match');
                            z = sscanf(cmdout{2}(40:61),'%f');
                            xyz(e,:) =[x y z];
                        end
                        % % UNIX only
                        % if exist(filename, 'file') == 2
                        %     [status,cmdout] = system(sprintf('grep ''STAX   %s'' %s',upper(this.parent.getMarkerName4Ch),filename));
                        %     if ~isempty(cmdout)
                        %         nl_id = find(cmdout==char(10));
                        %         x = sscanf(cmdout(nl_id(1)+(48:68)),'%f');
                        %         [status,cmdout] = system(sprintf('grep ''STAY   %s'' %s',upper(this.parent.getMarkerName4Ch),filename));
                        %         nl_id = find(cmdout==char(10));
                        %         y = sscanf(cmdout(nl_id(1)+(48:68)),'%f');
                        %         [status,cmdout] = system(sprintf('grep ''STAZ   %s'' %s',upper(this.parent.getMarkerName4Ch),filename));
                        %         nl_id = find(cmdout==char(10));
                        %         z = sscanf(cmdout(nl_id(1)+(48:68)),'%f');
                        %         xyz(e,:) =[x y z];
                        %     end
                    else
                        [pathstr, name, ext] = fileparts(remote_file_name);
                        
                        this.log.addWarning(sprintf(' File %s not found',[name, ext]));
                    end
                end
            end
        end
        
        function [ztd, gn ,ge] = getIgsTropo(this, mode)
            % get the official IGS solution fro tropospheric paramters
            %
            % SYNTAX:
            %    xyz = getIgsTropo(this, mode)
            if nargin < 2
                mode = 'interp_value';
            end
            ztd = nan(size(this.ztd));
            gn = nan(size(this.tgn));
            ge = nan(size(this.tge));
            tsc = Tropo_Sinex_Compare();
            sta_name = this.parent.getMarkerName4Ch;
            [missing_days] = tsc.addIGSOfficialStation(sta_name, this.time);
            lid_excl = false(this.time.length,1);
            for d = missing_days
                d1 = GPS_Time.fromMJD(d);
                d2 = d1.getCopy();
                d2.addSeconds(86400);
                lid_excl = lid_excl | (this.time> d1 & this.time < d2);
            end
            if ~isfield(tsc.results, 'r2') || ~isfield(tsc.results.r2, ['r' upper(sta_name)])
                this.log.addWarning(sprintf('No IGS solutions found for station %s', sta_name));
            else
                if strcmpi(mode,'interp_value')
                    ztd  = interp1(tsc.results.r2.(['r' upper(sta_name)]).time.getMatlabTime,tsc.results.r2.(['r' upper(sta_name)]).ztd,this.time.getMatlabTime,'linear');
                    gn  = interp1(tsc.results.r2.(['r' upper(sta_name)]).time.getMatlabTime,tsc.results.r2.(['r' upper(sta_name)]).tgn,this.time.getMatlabTime,'linear');
                    ge  = interp1(tsc.results.r2.(['r' upper(sta_name)]).time.getMatlabTime,tsc.results.r2.(['r' upper(sta_name)]).tge,this.time.getMatlabTime,'linear');
                    ztd(lid_excl) = nan;
                    gn(lid_excl) = nan;
                    ge(lid_excl) = nan;
                elseif strcmpi(mode,'difference')
                    ztd = timeSeriesComparison(tsc.results.r2.(['r' upper(sta_name)]).time.getMatlabTime, tsc.results.r2.(['r' upper(sta_name)]).ztd, this.time.getMatlabTime, this.ztd,'aggregate');
                    gn = timeSeriesComparison(tsc.results.r2.(['r' upper(sta_name)]).time.getMatlabTime, tsc.results.r2.(['r' upper(sta_name)]).tgn, this.time.getMatlabTime, this.tgn,'aggregate');
                    ge = timeSeriesComparison(tsc.results.r2.(['r' upper(sta_name)]).time.getMatlabTime, tsc.results.r2.(['r' upper(sta_name)]).tge, this.time.getMatlabTime, this.tge,'aggregate');
                    
                end
            end
        end
        
        
        function remEpoch(this, ep_idx)
            % remove epochs from the results
            %
            % SYNTAX
            % this.remEpoch(ep_idx)
            
            this.time.remEpoch(ep_idx);
            if ~isempty(this.zwd)
                this.zwd(ep_idx) = [];
            end
            if ~isempty(this.pressure)
                this.pressure(ep_idx) = [];
            end
            if ~isempty(this.temperature)
                this.temperature(ep_idx) = [];
            end
            if ~isempty(this.humidity)
                this.humidity(ep_idx) = [];
            end
            if ~isempty(this.desync)
                this.desync(ep_idx) = [];
            end
            if ~isempty(this.dt_ip)
                this.dt_ip(ep_idx) = [];
            end
            if ~isempty(this.dt)
                this.dt(ep_idx) = [];
            end
            if ~isempty(this.apr_zhd)
                this.apr_zhd(ep_idx) = [];
            end
            if ~isempty(this.ztd)
                this.ztd(ep_idx) = [];
            end
            if ~isempty(this.apr_zwd)
                this.apr_zwd(ep_idx) = [];
            end
            if ~isempty(this.pwv)
                this.pwv(ep_idx) = [];
            end
            if ~isempty(this.tgn)
                this.tgn(ep_idx) = [];
            end
            if ~isempty(this.tge)
                this.tge(ep_idx) = [];
            end
            if ~isempty(this.n_sat_ep)
                this.n_sat_ep(ep_idx) = [];
            end
            % remove from sat struct
            if ~isempty(this.sat.outliers)
                this.sat.outliers(ep_idx,:) = [];
            end
            if ~isempty(this.sat.cycle_slip)
                this.sat.cycle_slip(ep_idx,:) = [];
            end
            if ~isempty(this.sat.quality)
                this.sat.quality(ep_idx,:) = [];
            end
            if ~isempty(this.sat.az)
                this.sat.az(ep_idx,:) = [];
            end
            if ~isempty(this.sat.el)
                this.sat.el(ep_idx,:) = [];
            end
            if ~isempty(this.sat.mfw)
                this.sat.mfw(ep_idx,:) = [];
            end
            if ~isempty(this.sat.mfh)
                this.sat.mfh(ep_idx,:) = [];
            end
        end
        
    end
    
    % ==================================================================================================================================================
    %% METHODS IMPORT / EXPORT
    % ==================================================================================================================================================
    
    methods
        function injectResult(this, rec_work, rate)
            % inject the results of receiver work into receiver output
            %
            % SYNTAX
            %  this.injectResult(rec_work)
            
            state = Core.getCurrentSettings();
            if nargin < 3 || isempty(rate)
                if rec_work.time.getRate == state.getTropoOutRate
                    rate = 0;
                else
                    rate = state.getTropoOutRate();
                end
            end
                
            if ~(rec_work.isEmpty || rec_work.flag_currupted || not((rec_work.isPreProcessed && rec_work.quality_info.s0_ip < 2*1e2 && ~isempty(rec_work.quality_info.s0) && ~isnan(rec_work.quality_info.s0) && ~(rec_work.quality_info.s0 < 1e-5))))
                % set the id_sync only to time in between out times
                basic_export = false;
                id_sync_old = rec_work.getIdSync();
                
                this.used_sys_c = unique([this.used_sys_c rec_work.getActiveSys]);
                this.selected_sys_c = Core.getConstellationCollector.sys_c;
                if isempty(id_sync_old)
                    rec_work.id_sync = 1 : rec_work.time.length;
                    basic_export = true;
                end                
                is_last_session = rec_work.getTime.last >= this.state.sss_date_stop;
                % NOTE TROPO SMOOTHING.
                % in case of tropo smoothing we keep only the right buffer
                % because they will be used in the next session for
                % smoothing. The time of the left buffer is the one already
                % present in out. This means, that in the end the final
                % epochs of the output are the central ones of each
                % session and not the ones of the buffers.
                
                if any(rate)
                    id_sync_bk = rec_work.getIdSync;
                    id_sync = rec_work.getIdSync;
                    sync_time = round(rec_work.getTime.getNominalTime.getMatlabTime*86400,3);
                    id_ss = mod(sync_time, rate) == 0;
                    rec_work.id_sync = id_sync(id_ss);
                end
                
                rec_work.cropIdSync4out(true, ~this.state.isSmoothTropoOut() || is_last_session);
                work_time = rec_work.getTime();
                if ~work_time.isEmpty
                    initial_len = this.time.length;
                    is_this_empty = this.time.isEmpty;
                    if is_this_empty
                        idx1 = 1;
                        idx2 = 0;
                        this.time = work_time;
                    else
                        if this.state.isSmoothTropoOut()
                            time_old = this.time.getCopy();
                            re_time_bf = time_old.getNominalTime ;
                            smt_buf_rgt = re_time_bf.last;
                        end
                        [this.time, idx1, idx2] = this.time.injectBatch(work_time); % WARNING for tropo smoothing: the epoch before the end of the previous window will keep the their own epochs, the one after will keep the epoch of the work time
                        if this.state.isSmoothTropoOut()
                            smt_buf_lft = rec_work.time.getNominalTime().first();
                            idx_smt1 = re_time_bf >= smt_buf_lft;
                            idx_smt2 = rec_work.time.getEpoch(id_sync_old).getNominalTime <= smt_buf_rgt;
                            time_1 = time_old.getEpoch(idx_smt1);
                            time_2 = rec_work.time.getEpoch(id_sync_old(idx_smt2));
                        end
                    end
                    %%% inject data
                    if ~basic_export
                        % Inject times
                        if this.state.flag_out_dt
                            this.dt       = Core_Utils.injectData(this.dt, rec_work.getDt(), idx1, idx2);
                            this.desync   = Core_Utils.injectData(this.desync, rec_work.getDesync(), idx1, idx2);
                            this.dt_ip    = Core_Utils.injectData(this.dt_ip, rec_work.getDtIp(), idx1, idx2);
                            %this.n_sat_ep = Core_Utils.injectData(this.n_sat_ep, rec_work.getNSat(), idx1, idx2);
                        end
                        if this.state.flag_out_apr_tropo
                            this.apr_zhd = Core_Utils.injectData(this.apr_zhd, rec_work.getAprZhd(), idx1, idx2);
                            this.apr_zwd = Core_Utils.injectData(this.apr_zwd, rec_work.getAprZwd(), idx1, idx2);
                        end
                        
                        % Inject used meteo parameters
                        if this.state.flag_out_pth
                            [p, t, h]         = rec_work.getPTH(true);
                            this.pressure     = Core_Utils.injectData(this.pressure, p, idx1, idx2);
                            this.temperature  = Core_Utils.injectData(this.temperature, t, idx1, idx2);
                            this.humidity     = Core_Utils.injectData(this.humidity, h, idx1, idx2);
                        end
                        
                        % Inject mapping functions
                        if this.state.flag_out_mf
                            [mfh, mfw]            = rec_work.getSlantMF();
                            this.sat.mfw          = Core_Utils.injectData(this.sat.mfw, mfw, idx1, idx2);
                            this.sat.mfh          = Core_Utils.injectData(this.sat.mfh, mfh, idx1, idx2);
                        end
                        
                        % Inject outliers and cs
                        if this.state.flag_out_ocs
                            this.sat.outliers   = Core_Utils.injectData(this.sat.outliers, rec_work.getObsOutSat(), idx1, idx2);
                            this.sat.cycle_slip = Core_Utils.injectData(this.sat.cycle_slip, rec_work.getObsCsSat(), idx1, idx2);
                        end
                        
                        if this.state.isNSatOut()
                            cc = Core.getState.getConstellationCollector;
                            % all sats
                            if isempty(this.quality_info.n_spe)
                                this.quality_info.n_spe = struct('A', [], 'G', [], 'R', [], 'E', [], 'J', [], 'C', [], 'I', []);
                            end
                            this.quality_info.n_spe.A = Core_Utils.injectData(this.quality_info.n_spe.A, rec_work.quality_info.n_spe.A(rec_work.getIdSync), idx1, idx2);
                            % for each constellations
                            for sys_c = cc.getActiveSysChar
                                if ~isempty(rec_work.quality_info.n_spe.(sys_c))
                                    this.quality_info.n_spe.(sys_c) = Core_Utils.injectData(this.quality_info.n_spe.(sys_c), rec_work.quality_info.n_spe.(sys_c)(rec_work.getIdSync), idx1, idx2);
                                end
                            end
                        end
                        
                        if ~this.state.isSmoothTropoOut() || is_this_empty
                            % Inject tropo related parameters
                            if this.state.flag_out_ztd
                                this.ztd     = Core_Utils.injectData(this.ztd, rec_work.getZtd(), idx1, idx2);
                            end
                            if this.state.flag_out_zwd
                                this.zwd     = Core_Utils.injectData(this.zwd, rec_work.getZwd(), idx1, idx2);
                            end
                            if this.state.flag_out_pwv
                                this.pwv     = Core_Utils.injectData(this.pwv, rec_work.getPwv(), idx1, idx2);
                            end
                            if this.state.flag_out_tropo_g
                                [gn, ge]     = rec_work.getGradient();
                                this.tgn     = Core_Utils.injectData(this.tgn, gn, idx1, idx2);
                                this.tge     = Core_Utils.injectData(this.tge, ge, idx1, idx2);
                            end
                            if this.state.isResOut
                                if isempty(this.sat.res)
                                    this.sat.res = Residuals();
                                end
                                
                                % Get the residual only in the time span relative to the session (no buffers)
                                res = rec_work.sat.res.getCopy();
                                [~, lim] = this.state.getSessionLimits;                                
                                res.cutEpochs(lim);
                                this.sat.res.injest(res);
                            end
                        else
                            % there is probably smoothing
                            % save idx, they might be useful
                            bk_idx1 = idx1;
                            bk_idx2 = idx2;
                        end
                    end
                    
                    if (initial_len == size(this.sat.az,1)) && (idx2 == 0) && ~is_this_empty
                        [~, idx1, idx2] = this.time.injectBatch(work_time);
                    end
                    [az, el] = rec_work.getAzEl;
                    if this.state.flag_out_azel
                        this.sat.az      = Core_Utils.injectData(this.sat.az, az, idx1, idx2);
                        this.sat.el      = Core_Utils.injectData(this.sat.el, el, idx1, idx2);
                    end
                    if this.state.flag_out_quality
                        this.sat.quality = Core_Utils.injectData(this.sat.quality, rec_work.getQuality(), idx1, idx2);
                    end
                                        
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
                    
                    % Add antennna offset
                    coo_work = rec_work.getPos;
                    coo_work.addOffset(-[rec_work.parent.ant_delta_en rec_work.parent.ant_delta_h]);
                    xyz_work = coo_work.getXYZ;
                    
                    this.xyz      = Core_Utils.injectData(this.xyz, xyz_work, idx1, idx2, [data_len, 3]);
                    xyz_vcv = rec_work.getVCVXYZ;
                    if ~isempty(xyz_vcv)
                        this.xyz_vcv  = Core_Utils.injectData(this.xyz_vcv, xyz_vcv, idx1, idx2, [data_len, 6]);
                    end
                    this.enu      = Core_Utils.injectData(this.enu, rec_work.getPosENU, idx1, idx2, [data_len, 3]);
                    
                    this.quality_info.s0_ip     = Core_Utils.injectData(this.quality_info.s0_ip, rec_work.quality_info.s0_ip, idx1, idx2, [data_len, 1]);
                    this.quality_info.s0        = Core_Utils.injectData(this.quality_info.s0, rec_work.quality_info.s0, idx1, idx2, [data_len, 1]);
                    this.quality_info.n_epochs  = Core_Utils.injectData(this.quality_info.n_epochs, rec_work.quality_info.n_epochs, idx1, idx2, [data_len, 1]);
                    this.quality_info.n_obs     = Core_Utils.injectData(this.quality_info.n_obs, rec_work.quality_info.n_obs, idx1, idx2, [data_len, 1]);
                    this.quality_info.n_out     = Core_Utils.injectData(this.quality_info.n_out, rec_work.quality_info.n_out, idx1, idx2, [data_len, 1]);
                    this.quality_info.n_sat     = Core_Utils.injectData(this.quality_info.n_sat, rec_work.quality_info.n_sat, idx1, idx2, [data_len, 1]);
                                
                    this.quality_info.n_sat_max = Core_Utils.injectData(this.quality_info.n_sat_max, rec_work.quality_info.n_sat_max, idx1, idx2, [data_len, 1]);
                    this.quality_info.fixing_ratio = Core_Utils.injectData(this.quality_info.fixing_ratio, rec_work.quality_info.fixing_ratio, idx1, idx2, [data_len, 1]);
                    
                    
                    
                    % reset the old  complete id_sync
                    rec_work.id_sync = id_sync_old;
                    % inject with smoothing
                    if ~basic_export && ~is_this_empty && this.state.isSmoothTropoOut()
                        rec_work.cropIdSync4out(false, is_last_session); % if this is the last session cut the right part of the data
                        if is_last_session
                            idx_smt2 = idx_smt2(1 : numel(rec_work.getZtd));
                        end
                        % We have to find the first epoch in the receiver
                        % being pushed that is greater than the start of
                        % the session. This is done beacause might be that
                        % the first epoch are marked as outlier and thus
                        % not computed.
                        new_time = rec_work.getTime().getNominalTime;
                        first_new_time = new_time.getEpoch(find(new_time >= rec_work.out_start_time.getNominalTime, 1, 'first'));
                        clear new_time;
                        if ~isempty(time_1) && ~isempty(time_2)
                        id_stop     = find(time_1.getNominalTime >= first_new_time, 1, 'first'); % The first id of the new session
                        id_start    = find(time_2.getNominalTime >= first_new_time, 1, 'first'); % The first id of the new session
                        no_smooth = false;
                        else
                            id_stop = [];
                        end
                        if ~isempty(id_stop)
                            if this.state.flag_out_ztd
                                this.ztd     = Core_Utils.injectSmtData(zero2nan(this.ztd), zero2nan(rec_work.getZtd()), idx_smt1, idx_smt2, time_1, time_2, id_stop, id_start);
                            end
                            if this.state.flag_out_zwd
                                this.zwd     = Core_Utils.injectSmtData(zero2nan(this.zwd), zero2nan(rec_work.getZwd()), idx_smt1, idx_smt2, time_1, time_2, id_stop, id_start);
                            end
                            if this.state.flag_out_pwv
                                this.pwv     = Core_Utils.injectSmtData(zero2nan(this.pwv), zero2nan(rec_work.getPwv()), idx_smt1, idx_smt2, time_1, time_2, id_stop, id_start);
                            end
                            if this.state.flag_out_tropo_g
                                [gn, ge]     = rec_work.getGradient();
                                this.tgn     = Core_Utils.injectSmtData(zero2nan(this.tgn), zero2nan(gn), idx_smt1, idx_smt2, time_1, time_2, id_stop, id_start);
                                this.tge     = Core_Utils.injectSmtData(zero2nan(this.tge), zero2nan(ge), idx_smt1, idx_smt2, time_1, time_2, id_stop, id_start);
                            end
                            if this.state.isResOut
                                if isempty(this.sat.res)
                                    this.sat.res = Residuals();
                                end                                
                                
                                % Get the residual only in the time span relative to the session (no buffers)
                                res = rec_work.sat.res.getCopy();
                                [~, lim] = this.state.getSessionLimits;                                
                                res.cutEpochs(lim);
                                this.sat.res.injest(res);
                            end
                        else
                            % Inject tropo related parameters
                            if this.state.flag_out_ztd
                                tmp = rec_work.getZtd();
                                this.ztd     = [this.ztd; tmp(~idx_smt2)];
                            end
                            if this.state.flag_out_zwd
                                tmp = rec_work.getZwd();
                                this.zwd     = [this.zwd; tmp(~idx_smt2)];
                            end
                            if this.state.flag_out_pwv
                                tmp = rec_work.getPwv();
                                this.pwv     = [this.pwv; tmp(~idx_smt2)];
                            end
                            if this.state.flag_out_tropo_g
                                [gn, ge]     = rec_work.getGradient();
                                this.tgn     = [this.tgn; gn(~idx_smt2)];
                                this.tge     = [this.tge; ge(~idx_smt2)];
                            end
                            if this.state.isResCoOut
                                res_in = rec_work.getU1();
                                this.sat.res = [this.sat.res; res_in(~idx_smt2,:)];
                            end
                        end
                        rec_work.id_sync = id_sync_old; % restore id_sync_old
                    end
                    
                    %--- append additional coo
                    if this.state.flag_coo_rate
                        if isempty(this.add_coo)
                            is_empty_coo = true;
                            this.add_coo = struct('rate',[],'coo',[]);
                        else
                            is_empty_coo = false;
                        end
                        for i = 1 : length(rec_work.add_coo)
                            if is_empty_coo
                                this.add_coo(i) = struct('rate',[],'coo',[]);
                                this.add_coo(i).rate = rec_work.add_coo(i).rate;
                                this.add_coo(i).coo = rec_work.add_coo(i).coo.getCopy();
                            else
                                time_o = rec_work.add_coo(i).coo.time.getCopy();
                                coo_o = rec_work.add_coo(i).coo.getCopy();
                                discard_time = work_time.first;
                                discard_time.addSeconds(-this.add_coo(i).rate/2);
                                idx_rem = time_o < discard_time;
                                coo_o.rem(idx_rem);
                                [this.add_coo(i).coo.time, idx1, idx2] = this.add_coo(i).coo.time.injectBatch(coo_o.time);
                                this.add_coo(i).coo.xyz    = Core_Utils.injectData(this.add_coo(i).coo.xyz , coo_o.xyz , idx1, idx2);
                                if ~isempty(this.add_coo(i).coo.Cxx) && ~isempty(rec_work.add_coo(i).coo.Cxx)
                                    this.add_coo(i).coo.Cxx    = [this.add_coo(i).coo.Cxx(:,:,1 : idx1 - 1); coo_o.Cxx; this.add_coo(i).coo.Cxx(:,:,idx2 + 1 : end)];
                                end
                            end
                        end
                    end
                else
                    rec_work.id_sync = id_sync_old; % restore id_sync_old
                end
                
                if any(rate)
                    rec_work.id_sync = id_sync_bk;
                end
            end
            
            log = Core.getLogger();
            log.addMarkedMessage(sprintf('Computed results for receiver "%s" have been imported into out object', this.parent.getMarkerName4Ch()));
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
        
        function fh_list = showDt(this)
            % Plot Clock error
            %
            % SYNTAX
            %   this.plotDt
            
            fh_list = [];
            rec = this;
            if ~isempty(rec)
                t = rec.time.getMatlabTime();
                if isempty(t)
                    Core.getLogger.addError('No clock found in Receiver Output object\n');
                else
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: Dt Err', f.Number); f.NumberTitle = 'off';
                    
                    fh_list = f;
                    fig_name = sprintf('Dt_%s_%s', this.parent.getMarkerName4Ch, this.time.first.toString('yyyymmdd_HHMM'));
                    f.UserData = struct('fig_name', fig_name);
                    
                    l_list = {};
                    data = rec.getDesync;
                    if ~isempty(data)
                        l_list = [l_list {'De-Sync time (from RINEX)'}];
                        Core_Utils.plotSep(t, data * Core_Utils.V_LIGHT, '-k', 'LineWidth', 2);
                    end
                    hold on;
                    data = rec.getDtIP;
                    if ~isempty(data)
                        l_list = [l_list {'{\Delta}t from pre-processing'}];
                        Core_Utils.plotSep(t, data * Core_Utils.V_LIGHT, '-', 'LineWidth', 2);
                    end
                    data = rec.getDt();
                    if ~isempty(data)
                        l_list = [l_list {'{\Delta}t from the last step'}];
                        Core_Utils.plotSep(t, data * Core_Utils.V_LIGHT, '-', 'LineWidth', 2);
                    end
                    data = rec.getTotalDt();
                    if ~isempty(data)
                        l_list = [l_list {'{\Delta}t final '}];
                        Core_Utils.plotSep(t, data * Core_Utils.V_LIGHT, '-', 'LineWidth', 2);
                    end
                    if isempty(l_list)
                        Core.getLogger.addError('No clock found in Receiver Output object\n');
                        close(f);
                    else
                        legend(l_list, 'Location', 'NorthEastOutside');
                        
                        xlim([t(1) t(end)]); setTimeTicks(4); h = ylabel('receiver clock error [m]'); h.FontWeight = 'bold';
                        h = title(sprintf('dt - receiver %s', rec.parent.getMarkerName),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                        Core_UI.beautifyFig(f);
                        Core_UI.addExportMenu(f);
                        Core_UI.addBeautifyMenu(f);
                        f.Visible = 'on'; drawnow;
                    end
                end
            end
        end
        
        function fh_list = showProcessingQualityInfo(this)
            % Show quality info indexes for processing
            %   number of epochs
            %   number of observations
            %   max number of sat seen in one epoch
            %   sigma 0 of carrier phase processing
            %   sigma 0 of code pre-processing
            %
            % SYNTAX:
            %   this.showProcessingQualityInfo
            
            fh_list = [];
            if ~(this.isEmpty) && size(this.quality_info.n_obs, 1) > 1
                this.log.addMessage('Plotting processing quality info');
                
                win = figure('Visible', 'on', ...
                    'NumberTitle', 'off');
                win.Name = sprintf('%03d: %s, Quality Info', win.Number, this.parent.getMarkerName4Ch);
                
                fh_list = [fh_list; win]; 
                fig_name = sprintf('Quality_Info_%s_%s', this.parent.getMarkerName4Ch, this.time.first.toString('yyyymmdd_HHMM'));
                win.UserData = struct('fig_name', fig_name);
                
                % Single index
                n_plot = 6;
                
                win.Units = 'pixels';
                maximizeFig(win);

                %scroller = uix.ScrollingPanel('Parent', win);
                main_vb = uix.VBox('Parent', win, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                
                uicontrol('Parent', main_vb, ...
                    'Style', 'Text', ...
                    'String', sprintf('Processing quality info for rec %s\n', upper(this.parent.getMarkerName4Ch)), ...
                    'ForegroundColor', Core_UI.BLACK, ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', Core_UI.getFontSize(10), ...
                    'FontWeight', 'Bold', ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);

                container = uix.Grid('Parent', main_vb, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                
                %h = title(, 'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                                
                for i = 1 : n_plot
                    tmp_box(i) = uix.VBox('Parent', container, ...
                        'Padding', 5, ...
                        'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                end
                main_vb.Heights = [30 -1];
                container.Heights = [-1 -1 -1];
                container.Widths = [-1 -1];
                for i = 1 : n_plot
                    ax(i) = axes('Parent', tmp_box(i));
                end
                %scroller.Heights = sum(container.Heights);
                
                color_order = handle(gca).ColorOrder;
                t = this.getPositionTime().getMatlabTime();
                
                plot(ax(1), t, this.quality_info.n_epochs, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                h = ylabel(ax(1), sprintf('# epochs'), 'interpreter', 'none'); h.FontWeight = 'bold';
                h = title(ax(1), 'Number of valid epochs', 'interpreter', 'none'); h.FontWeight = 'bold';
                
                plot(ax(2), t, this.quality_info.n_obs, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:)); hold on;
                h = ylabel(ax(2), sprintf('# obs'), 'interpreter', 'none'); h.FontWeight = 'bold';
                h = title(ax(2), 'Total number of observations used', 'interpreter', 'none'); h.FontWeight = 'bold';
                                
                if isfield(this.quality_info, 'n_out')
                    plot(ax(3), t, this.quality_info.n_out, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(4,:)); hold on;
                    h = ylabel(ax(3), '# outliers'); h.FontWeight = 'bold';
                    h = title(ax(3), 'Number of observations removed as outliers', 'interpreter', 'none'); h.FontWeight = 'bold';
                end
                
                plot(ax(4), t, this.quality_info.n_sat_max, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:)); hold on;
                h = ylabel(ax(4), 'max # sat'); h.FontWeight = 'bold';
                h = title(ax(4), 'Maximum number of satellites seen in one epoch', 'interpreter', 'none'); h.FontWeight = 'bold';
                
                plot(ax(5), t, this.quality_info.s0_ip * 1e2, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(6,:)); hold on;
                h = ylabel(ax(5), 's0 pp [cm]'); h.FontWeight = 'bold';
                h = title(ax(5), 'Sigma0 as estimated from the Least Square solution (on pre-processing)', 'interpreter', 'none'); h.FontWeight = 'bold';
                
                plot(ax(6), t, this.quality_info.s0 * 1e2, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(5,:)); hold on;
                h = ylabel(ax(6), 's0 [cm]'); h.FontWeight = 'bold';
                h = title(ax(6), 'Final sigma0 as estimated from the Least Square solution', 'interpreter', 'none'); h.FontWeight = 'bold';
                
                for i = 1 : n_plot
                    if (t(end) > t(1))
                        xlim(ax(i), [t(1) t(end)]);
                    end
                    grid(ax(i), 'on');
                    setTimeTicks(ax(i), 4, 'auto');
                    ax(i).FontSize = Core_UI.getFontSize(9);
                end
                
                linkaxes(ax, 'x');
                
                fh = win; Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'dark');
                fh.Visible = 'on';
            else
                Core.getLogger.addMessage(sprintf('It seems that only one session is present - these plots are not supported for %s', this.parent.getMarkerName4Ch));
            end            
        end
                
        function fh_list = showOutliersAndCycleSlip(this, sys_c_list)
            % Plot the outliers found
            % SYNTAX this.showOutliers()
            
            fh_list = [];
            cc = Core.getState.getConstellationCollector;
            if nargin == 1
                sys_c_list = cc.getAvailableSys;
            end
            for sys_c = sys_c_list
                f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s CS, Out %s', f.Number, this.parent.getMarkerName4Ch, cc.getSysName(sys_c)); f.NumberTitle = 'off';
                
                fh_list = [fh_list; f]; %#ok<AGROW>
                fig_name = sprintf('OCS_%s_%s_%s', this.parent.getMarkerName4Ch, cc.getSysName(sys_c), this.time.first.toString('yyyymmdd_HHMM'));
                f.UserData = struct('fig_name', fig_name);
                
                ep = repmat((1: this.time.length)',1, size(this.sat.outliers, 2));
                if isempty(ep)
                    close(f);
                else
                    for prn = cc.prn(cc.system == sys_c)'
                        s = cc.getIndex(sys_c, prn);
                        cs = ep(this.sat.cycle_slip(:, s) ~= 0);
                        sat_on = ep(this.sat.az(:, s) ~= 0);
                        plot(sat_on,  prn * ones(size(sat_on)), '.', 'MarkerSize', 10, 'Color', [0.7 0.7 0.7]);
                        hold on;
                        plot(cs,  prn * ones(size(cs)), '.k', 'MarkerSize', 20);
                        out = ep(this.sat.outliers(:, s) ~= 0);
                        plot(out,  prn * ones(size(out)), '.', 'MarkerSize', 20, 'Color', [1 0.4 0]);
                    end
                    prn_ss = cc.prn(cc.system == sys_c);
                    xlim([1 this.time.length]);
                    ylim([min(prn_ss) - 1 max(prn_ss) + 1]);
                    h = ylabel('PRN'); h.FontWeight = 'bold';
                    ax = gca(); ax.YTick = prn_ss;
                    grid on;
                    h = xlabel('epoch'); h.FontWeight = 'bold';
                    h = title(sprintf('%s %s cycle-slip(b) & outlier(o)', cc.getSysName(sys_c), this.parent.marker_name), 'interpreter', 'none'); h.FontWeight = 'bold';
                    Core_UI.beautifyFig(f, 'dark');
                    Core_UI.addExportMenu(f);
                    Core_UI.addBeautifyMenu(f);
                    f.Visible = 'on'; drawnow;
                end
            end
        end
        
        function fh_list = showOutliersAndCycleSlip_p(this, sys_c_list)
            % Plot Signal to Noise Ration in a skyplot
            % SYNTAX this.plotSNR(sys_c)
            
            cc = Core.getState.getConstellationCollector;
            
            fh_list = [];
            % SNRs
            if nargin == 1
                sys_c_list = cc.getAvailableSys;
            end
            
            for sys_c = sys_c_list
                if isempty(this.sat.az) || isempty(this.sat.el)
                    Core.getLogger.addError('No azimuth and elevetion present in the receiver\n');
                else
                    is_empty = false;
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: CS, Out %s', f.Number, sys_c); f.NumberTitle = 'off';
                    polarScatter([],[],1,[]);
                    hold on;
                    
                    fh_list = [fh_list; f]; %#ok<AGROW>
                    fig_name = sprintf('OCS_polar_%s_%s_%s', this.parent.getMarkerName4Ch, cc.getSysName(sys_c), this.time.first.toString('yyyymmdd_HHMM'));
                    f.UserData = struct('fig_name', fig_name);
                    
                    for s = cc.getGoIds(sys_c)
                        az = this.sat.az(:,s);
                        el = this.sat.el(:,s);
                        is_empty = is_empty + ~(any(az));
                        if any(az)
                            cs = sum(this.sat.cycle_slip(:, s), 2) > 0;
                            out = sum(this.sat.outliers(:, s), 2) > 0;
                            sat_on = (this.sat.az(:, s) ~= 0);
                            
                            decl_n = (serialize(90 - el(sat_on)) / 180*pi) / (pi/2);
                            x = sin(az(sat_on)/180*pi) .* decl_n; x(az(sat_on) == 0) = [];
                            y = cos(az(sat_on)/180*pi) .* decl_n; y(az(sat_on) == 0) = [];
                            plot(x, y, '.', 'MarkerSize', 7, 'Color', [0.7 0.7 0.7]);
                            
                            decl_n = (serialize(90 - el(cs)) / 180*pi) / (pi/2);
                            x = sin(az(cs)/180*pi) .* decl_n; x(az(cs) == 0) = [];
                            y = cos(az(cs)/180*pi) .* decl_n; y(az(cs) == 0) = [];
                            plot(x, y, '.k', 'MarkerSize', 20);
                            
                            decl_n = (serialize(90 - el(out)) / 180*pi) / (pi/2);
                            x = sin(az(out)/180*pi) .* decl_n; x(az(out) == 0) = [];
                            y = cos(az(out)/180*pi) .* decl_n; y(az(out) == 0) = [];
                            plot(x, y, '.', 'MarkerSize', 20, 'Color', [1 0.4 0]);
                        end
                    end
                    if is_empty == numel(cc.getGoIds(sys_c))
                        close(f);
                        Core.getLogger.addError('No azimuth / elevation found in Receiver Output object\n');
                    else
                        h = title(sprintf('%s %s cycle-slip(b) & outlier(o)', cc.getSysName(sys_c), this.parent.marker_name), 'interpreter', 'none'); h.FontWeight = 'bold';
                        Core_UI.beautifyFig(f);
                        Core_UI.addExportMenu(f);
                        Core_UI.addBeautifyMenu(f);
                        f.Visible = 'on'; drawnow;
                    end
                end
            end
            
        end        
    end
end
