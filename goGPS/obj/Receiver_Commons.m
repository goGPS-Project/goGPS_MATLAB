%   CLASS Receiver_Commons
% =========================================================================
%
%
%   Class to store receiver common methods and abstract properties
%
% EXAMPLE
%   trg = Receiver_Commons();
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
classdef Receiver_Commons <  matlab.mixin.Copyable
    properties (SetAccess = public, GetAccess = public)
        parent         % habdle to parent object
        
        rid            % receiver interobservation biases
        flag_rid       % clock error for each obs code {num_obs_code}
    end
    
    % ==================================================================================================================================================
    %% PROPERTIES CELESTIAL INFORMATIONS
    % ==================================================================================================================================================
    
    properties (Abstract, SetAccess = public, GetAccess = public)
        sat
    end
    
    % ==================================================================================================================================================
    %% PROPERTIES POSITION
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        xyz            % position of the receiver (XYZ geocentric)
        xyz_vcv        % vcv matrix of the receiver (XYZ geocentric) var_x2 var_xy var_xz var_y2 var_yz var_z2
        enu            % position of the receiver (ENU local)
        
        lat            % ellipsoidal latitude
        lon            % ellipsoidal longitude
        h_ellips       % ellipsoidal height
        h_ortho        % orthometric height
        
        add_coo
        %         = struct( ...
        %             'coo',      [], ...    % additional estimated coo
        %             'time',         [], ...    % time of the coo
        %             'rate',          [] ...    % rate of the coo
        %             )
        
        
    end
    
    % ==================================================================================================================================================
    %% PROPERTIES TIME
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        time           % internal time ref of the stored epochs
        desync         % receiver clock desync (difference between nominal time and the time of the observation)
        dt_ip          % clock error correction estimated during init positioning
        dt             % reference clock error of the receiver [n_epochs x num_obs_code]
    end
    % ==================================================================================================================================================
    %% PROPERTIES TROPO
    % ==================================================================================================================================================
    
    properties (SetAccess = public, GetAccess = public)
        apr_zhd  % zenital hydrostatic delay           double   [n_epoch x 1]
        ztd      % total zenital tropospheric delay    double   [n_epoch x 1]
        zwd      % zenital wet delay                   double   [n_epoch x 1]
        apr_zwd  % apriori zenital wet delay           double   [n_epoch x 1]
        pwv      % precipitable water vapour           double   [n_epoch x 1]
        
        tgn      % tropospheric gradient north         double   [n_epoch x n_sat]
        tge      % tropospheric gradient east          double   [n_epoch x n_sat]
        
        tzer     % zernike modeling of the tropo
    end
    
    % ==================================================================================================================================================
    %% PROPERTIES QUALITY INDEXES
    % ==================================================================================================================================================
    
    properties
        quality_info = struct('s0', [], 's0_ip', [], 'n_epochs', [], 'n_obs', [], 'n_out', [], 'n_sat', [], 'n_sat_max', [], 'n_spe', [], 'fixing_ratio', [], 'C_pos_pos', []);
        a_fix
        s_rate
        n_sat_ep
    end
    
    % ==================================================================================================================================================
    %% PROPERTIES USEFUL HANDLES
    % ==================================================================================================================================================
    
    properties (SetAccess = protected, GetAccess = public)
        w_bar                                  % handle to waitbar
        state                                  % local handle of state;
        log                                    % handle to log
        %rf                                     % handle to reference farme
    end
    
    % ==================================================================================================================================================
    %% METHODS INIT - CLEAN - RESET - REM - IMPORT
    % ==================================================================================================================================================
    
    methods
        
        function clearHandles(this)
            this.log = [];
            this.state = [];
            this.w_bar = [];
        end
        
        function initHandles(this)
            this.log = Core.getLogger();
            this.state = Core.getState();
            this.w_bar = Go_Wait_Bar.getInstance();
        end
        
        function reset(this)
            this.time = GPS_Time();
            this.enu = [];
            this.lat = [];
            this.lon = [];
            
            this.h_ellips = [];
            this.h_ortho = [];
            
            this.quality_info = struct('s0', [], 's0_ip', [], 'n_epochs', [], 'n_obs', [], 'n_out', [], 'n_sat', [], 'n_sat_max', [], 'n_spe', [], 'fixing_ratio', [], 'C_pos_pos', []);
            
            this.a_fix = [];
            this.s_rate = [];
            
            this.xyz = [];
            
            this.apr_zhd  = [];
            this.zwd  = [];
            this.apr_zwd  = [];
            this.ztd  = [];
            this.pwv  = [];
            
            this.tgn = [];
            this.tge = [];
        end
    end
    % ==================================================================================================================================================
    %% METHODS GETTER - TIME
    % ==================================================================================================================================================
    
    methods
        function cc = getCC(this)
            % Get Constellation collector
            %
            % SYNTAX
            %   cc = this.getCC()
            cc = Core.getState.getConstellationCollector;
        end
        
        function toStringPos(this)
            % Display on screen information about the receiver position
            % SYNTAX this.toStringPos();
            for r = 1 : numel(this)
                if ~this(r).isEmpty && ~isempty(this(r).xyz)
                    [lat, lon, h_ellips, h_ortho] = this(r).getMedianPosGeodetic_mr();
                    this(r).log.addMarkedMessage(sprintf('Receiver %s   %11.7f  %11.7f    %12.7f m (ellipsoidal) - %12.7f (orthometric)', this(r).parent.getMarkerName, lat, lon, h_ellips, h_ortho));
                end
            end
        end
        
        function is_empty = isEmpty(this)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty();
            is_empty =  this.time.length() == 0;
        end
        
        function is_empty = isEmpty_mr(this)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty();
            is_empty =  zeros(numel(this), 1);
            for r = 1 : numel(this)
                is_empty(r) =  this(r).isEmpty();
            end
        end
        
        function len = length(this)
            % Return the time span of the receiver
            %
            % SYNTAX
            %   len = this.length();
            len = this.getTime.length();
        end
        
        function dt = getDt(this)
            dt =  this.dt;
        end
        
        function dt_ip = getDtIP(this)
            dt_ip = this.dt_ip;
        end
        
        % time
        function time = getCentralTime(this)
            % return the central epoch time stored in the a receiver
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getCentralTime()
            time = this(1).time.getCentralTime();
        end
        
        function [rate] = getRate(this)
            % SYNTAX
            %   rate = this.getRate();
            rate = this.time.getRate;
        end
        
        function dt = getTotalDt(this)
            dt = this.getDt + this.getDtPrePro;
        end
        
        function coo = getPos(this)
            % return the positions computed for the receiver
            % as Coordinates object
            %
            % OUTPUT
            %   coo     coordinates object
            %
            % SYNTAX
            %   coo = this.getPos()
            
            if ~isempty(this.xyz)
                coo = Coordinates.fromXYZ(this.xyz);
            elseif ~isempty(this.parent.work.xyz)
                coo = Coordinates.fromXYZ(this.parent.work.xyz);
            else
                coo = Coordinates.fromXYZ([0 0 0]);
            end
            
        end
        
        function xyz = getPosXYZ(this)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   xyz     geocentric coordinates
            %
            % SYNTAX
            %   xyz = this.getPosXYZ()
            xyz = [];
            for r = 1 : numel(this)
                if ~isempty(this(r).xyz)
                    xyz = [xyz; this(r).xyz]; %#ok<AGROW>
                else
                    xyz = [xyz; this(r).parent.work.xyz]; %#ok<AGROW>
                end
            end
        end
        
        function xyz_vcv = getVCVXYZ(this)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   xyz     geocentric coordinates
            %
            % SYNTAX
            %   xyz = this.getPosXYZ()
            xyz_vcv = [];
            for r = 1 : numel(this)
                if ~isempty(this(r).xyz_vcv)
                    xyz_vcv = [xyz_vcv; this(r).xyz_vcv]; %#ok<AGROW>
                else
                    xyz_vcv = [xyz_vcv; this(r).parent.work.xyz_vcv]; %#ok<AGROW>
                end
            end
        end
        
        
        
        % position
        function [lat, lon, h_ellips, h_ortho] = getPosGeodetic(this)
            % Return the positions computed for the receiver
            %
            % OUTPUT
            %   lat      = latitude                      [rad]
            %   lon      = longitude                     [rad]
            %   h_ellips = ellipsoidal height            [m]
            %   lat_geoc = geocentric spherical latitude [rad]
            %   h_ortho  = orthometric height            [m]
            %
            % SYNTAX
            %   [lat, lon, h_ellips, h_ortho] = this.getGeodetic()
            
            coo = this.getPos();
            
            if nargout > 3
                [lat, lon, h_ellips, h_ortho] = coo.getGeodetic();
            elseif nargout > 2
                [lat, lon, h_ellips] = coo.getGeodetic();
            else
                [lat, lon] = coo.getGeodetic();
            end
        end
        
        function enu = getPosENU(this)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   enu     enu coordinates
            %
            % SYNTAX
            %   enu = this.getPosENU()
            enu = this.getPos().getENU();
        end
        
        function [utm, utm_zone] = getPosUTM(this)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   utm          utm coordinates
            %   utm_zone     utm zone
            %
            % SYNTAX
            %   [utm, utm_zone] = this.getPosUTM()
            [utm, utm_zone] =  this.getPos().getENU;
        end
        
        function enu = getBaselineENU(this, rec)
            % return the baseline computed for the receiver wrt another
            %
            % OUTPUT
            %   enu     enu coordinates
            %
            % SYNTAX
            %   enu = this.getPosENU()
            enu = this.getPosENU() - rec.getPosENU();
        end
        
        function xyz = getMedianPosXYZ(this)
            % return the computed median position of the receiver
            %
            % OUTPUT
            %   xyz     geocentric coordinates
            %
            % SYNTAX
            %   xyz = this.getMedianPosXYZ()
            
            xyz = this.getPosXYZ();
            xyz = median(xyz, 1, 'omitnan');
        end
        
        function enu = getMedianPosENU(this)
            % return the computed median position of the receiver
            %
            % OUTPUT
            %   enu     geocentric coordinates
            %
            % SYNTAX
            %   enu = this.getMedianPosENU()
            
            enu = this.getPosENU();
            enu = median(enu, 1, 'omitnan');
        end
        
        function [lat_d, lon_d, h_ellips, h_ortho] = getMedianPosGeodetic(this)
            % return the computed median position of the receiver
            %
            % OUTPUT
            %   lat         latitude  [deg]
            %   lon         longitude [deg]
            %   h_ellips    ellipsoidical heigth [m]
            %   h_ortho     orthometric heigth [m]
            %
            % SYNTAX
            %   [lat, lon, h_ellips, h_ortho] = this.getMedianPosGeodetic();
            for r = 1 : numel(this)
                xyz = this(r).getPosXYZ();
                xyz = median(xyz, 1);
                if ~isempty(xyz) && ~isempty(this(r))
                    [lat_d(r), lon_d(r), h_ellips(r)] = cart2geod(xyz);
                    if nargout == 4
                        ondu = getOrthometricCorr(lat_d(r), lon_d(r), Core.getRefGeoid());
                        h_ortho(r) = h_ellips(r) - ondu;
                    end
                    lat_d(r) = lat_d(r) / pi * 180;
                    lon_d(r) = lon_d(r) / pi * 180;
                else
                    lat_d(r) = nan;
                    lon_d(r) = nan;
                    h_ellips(r) = nan;
                    h_ortho(r) = nan;
                end
            end
        end
        
        % tropo
        function [tropo, time] = getTropoPar(sta_list, par_name)
            % get a tropo parameter among 'ztd', 'zwd', 'pwv', 'zhd'
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getAprZhd()
            
            tropo = {};
            time = {};
            for r = 1 : numel(sta_list)
                time{r} = sta_list(r).getTime();
                switch lower(par_name)
                    case 'ztd'
                        [tropo{r}] = sta_list(r).getZtd();
                    case 'zwd'
                        [tropo{r}] = sta_list(r).getZwd();
                        if isempty(tropo{r}) || all(isnan(zero2nan(tropo{r})))
                            [tropo{r}] = sta_list(r).getAprZwd();
                        end
                        
                    case 'gn'
                        [tropo{r}] = sta_list(r).getGradient();
                    case 'ge'
                        [~,tropo{r}] = sta_list(r).getGradient();
                    case 'pwv'
                        [tropo{r}] = sta_list(r).getPwv();
                    case 'zhd'
                        [tropo{r}] = sta_list(r).getAprZhd();
                    case 'nsat'
                        [tropo{r}] = sta_list(r).getNSat();
                end
            end
            
            if numel(tropo) == 1
                tropo = tropo{1};
                time = time{1};
            end
        end
        
        function slant_td = getSlantTD(this)
            % Get the slant total delay
            %
            % SYNTAX
            %   slant_td = this.getSlantTD();
            
            [mfh, mfw] = this.getSlantMF();
            n_sat = size(mfh,2);
            zwd = this.getZwd();
            apr_zhd = this.getAprZhd();
            [az, el] = this.getAzEl();
            [tgn, tge] = this.getGradient();
            if  isempty(zwd) || isempty(apr_zhd)
                slant_td = [];
            else
                if isfield(this.sat, 'avail_index') && any(this.sat.avail_index(:))
                    mfh(~this.sat.avail_index) = nan;
                    mfw(~this.sat.avail_index) = nan;
                end
                mfw(el < Core.getState.getCutOff) = nan;
                mfh(el < Core.getState.getCutOff) = nan;
                slant_td = nan2zero(zero2nan(repmat(zwd,1,n_sat).*mfw) ...
                    + zero2nan(repmat(apr_zhd,1,n_sat).*mfh));
            end
            if ~any(tgn(:)) || ~any(tge(:)) || isempty(mfh(:)) || isempty(mfw(:))
            else
                if isfield(this.sat, 'avail_index') && any(this.sat.avail_index(:))
                    el(~this.sat.avail_index) = nan;
                end
                el(el < Core.getState.getCutOff) = nan;
                if this.state.mapping_function_gradient == 1
                    cotan_term = Atmosphere.chenHerringGrad(el/180*pi);
                elseif this.state.mapping_function_gradient == 2
                    cotan_term = Atmosphere.macmillanGrad(el/180*pi).*mfw;
                end
                cosaz = zero2nan(cosd(az));
                sinaz = zero2nan(sind(az));
                slant_td = slant_td ...
                    + nan2zero(repmat(tgn,1,n_sat) .* cotan_term .* cosaz ...
                    + repmat(tge,1,n_sat) .* cotan_term .* sinaz);
            end
            if Main_Settings.getNumZerTropoCoef > 0 & ~isempty(this.tzer)
                degree = ceil(-3/2 + sqrt(9/4 + 2*(Main_Settings.getNumZerTropoCoef  -1)));
                zer_val = nan(size(el,1),size(el,2),Main_Settings.getNumZerTropoCoef);
                idx = el >  Core.getState.getCutOff & this.sat.avail_index;
                 rho = (90 - el)/(90);
                 if state.mapping_function_gradient == 1
                     cotan_term = Atmosphere.chenHerringGrad(el);
                 elseif state.mapping_function_gradient == 2
                     cotan_term = Atmosphere.macmillanGrad(el).*mfw;
                 end
                for s = 1 : size(this.sat.el,2)
                    zer_val(idx(:,s),s,:) = Core_Utils.getAllZernike(degree, az(idx(:,s),s)/180*pi, rho);
                end
                for i = 2 : Main_Settings.getNumZerTropoCoef
                    slant_td = slant_td + zero2nan(repmat(this.tzer(:,i-1),1,n_sat).*cotan_term.*zer_val(:,:,i));
                end
            end
        end
        
        function ztd = getZtd(this)
            % get ztd
            %
            % SYNTAX
            %   ztd = this.getZtd()
            if max(this.getIdSync) > numel(this.ztd)
                ztd = nan(size(this.getIdSync));
            else
                ztd = this.ztd(this.getIdSync);
            end
        end
        
        function sztd = getSlantZTD(this, smooth_win_size, id_extract)
            % Get the "zenithalized" total delay
            % SYNTAX
            %   sztd = this.getSlantZTD(<flag_smooth_data = 0>)
            if nargin < 3
                id_extract = 1 : this.getTime.length;
            end
            
            if ~isempty(this(1).ztd)
                [mfh, mfw] = this.getSlantMF();
                if isempty(mfh)
                    this(1).log.addWarning('ZTD and slants have not been computed');
                    sztd = [];
                else
                    sztd = bsxfun(@plus, (zero2nan(this.getSlantTD) - bsxfun(@times, mfh, this.getAprZhd)) ./ mfw, this.getAprZhd);
                    sztd(sztd <= 0) = nan;
                    sztd = sztd(id_extract, :);
                    
                    if nargin >= 2 && smooth_win_size > 0
                        t = this.getTime.getEpoch(id_extract).getRefTime;
                        for s = 1 : size(sztd,2)
                            id_ok = ~isnan(sztd(:, s));
                            if sum(id_ok) > 3
                                lim = getOutliers(id_ok);
                                lim = limMerge(lim, 2 * smooth_win_size / this.getRate);
                                
                                %lim = [lim(1) lim(end)];
                                for l = 1 : size(lim, 1)
                                    if (lim(l, 2) - lim(l, 1) + 1) > 3
                                        id_ok = lim(l, 1) : lim(l, 2);
                                        ztd = this.getZtd();
                                        sztd(id_ok, s) = splinerMat(t(id_ok), sztd(id_ok, s) - zero2nan(ztd(id_ok)), smooth_win_size, 0.05) + zero2nan(ztd(id_ok));                                        
                                    end
                                end
                            end
                        end
                    end
                end
            else
                this(1).log.addWarning('ZTD and slants have not been computed');
                sztd = [];
            end
        end
        
        function apr_zhd = getAprZhd(this)
            % get a-priori ZHD
            %
            % SYNTAX
            %   zhd = this.getAprZhd()
            if max(this.getIdSync) > numel(this.apr_zhd)
                apr_zhd = nan(size(this.getIdSync));
            else
                apr_zhd = this.apr_zhd(this.getIdSync);
            end
        end
        
        function [n_sat, n_sat_ss] = getNSat(this)
            % get num sta per epoch
            %
            % OUTPUT
            %   n_sat       total number of sat in view
            %   n_sat_ss    struct(.G .E .R ...) number of sat per constellation
            %
            % SYNTAX
            %   [n_sat, n_sat_ss] = this.getNSat()
            
            n_sat = nan;
            n_sat_ss = struct();
            cc = this.getCC;
            if isfield(this.quality_info, 'n_spe') && ~isempty(this.quality_info.n_spe)  && ~isempty(this.quality_info.n_sat)
                n_sat = this.quality_info.n_spe.A(this.getIdSync);
                for sys_c = cc.getActiveSysChar()
                    if ~isempty(this.quality_info.n_spe.(sys_c))
                        n_sat_ss.(sys_c) = this.quality_info.n_spe.(sys_c)(this.getIdSync);
                    else
                        n_sat_ss.(sys_c) = [];
                    end
                end
            else                
                n_sat_ss.G = [];
                if isempty(this.n_sat_ep)
                    % retrieve the n_sat from residuals
                    if this.state.isResCoOut && ~isempty(this.sat.res)
                        n_sat = sum(~isnan(zero2nan(this.sat.res(this.getIdSync,:))), 2);
                        for sys_c = cc.sys_c
                            n_sat_ss.(sys_c) = sum(~isnan(zero2nan(this.sat.res(this.getIdSync, cc.system == sys_c))), 2);
                        end
                    else
                        n_sat = nan(size(this.getIdSync));
                        n_sat_ss.G = n_sat;
                    end
                else
                    if (max(this.getIdSync) > numel(this.n_sat_ep))
                        n_sat = nan(size(this.getIdSync), 'single');
                        n_sat_ss.G = n_sat;
                    else
                        n_sat = this.n_sat_ep(this.getIdSync);
                        if any(serialize(this.sat.res(this.getIdSync,:)))
                            % retrieve the n_sat from residuals
                            n_sat = sum(~isnan(zero2nan(this.sat.res(this.getIdSync,:))), 2);
                            for sys_c = cc.sys_c
                                n_sat_ss.(sys_c) = sum(~isnan(zero2nan(this.sat.res(this.getIdSync, cc.system == sys_c))), 2);
                            end
                        end
                    end
                end
                if isempty(n_sat_ss.G)
                    n_sat_ss.G = n_sat;
                end
            end
        end
        
        function zwd = getZwd(this)
            % get zwd
            %
            % SYNTAX
            %   zwd = this.getZwd()
            if max(this.getIdSync) > numel(this.zwd)
                zwd = nan(size(this.getIdSync), 'single');
            else
                zwd = this.zwd(this.getIdSync);
                if isempty(zwd) || all(isnan(zero2nan(zwd)))
                    zwd = this.getAprZwd();
                end
            end
        end
        
        function pwv = getPwv(this)
            % get pwv
            %
            % SYNTAX
            %   pwv = this.getPwv()
            if max(this.getIdSync) > numel(this.pwv)
                pwv = nan(size(this.getIdSync), 'single');
            else
                pwv = this.pwv(this.getIdSync);
            end
        end
        
        function [gn ,ge, time] = getGradient(this)
            % SYNTAX
            % [gn ,ge, time] = getGradient(this)
            if isempty(this.tgn)
                gn = nan(length(this.getIdSync),1, 'single');
            else
                gn = this.tgn(this.getIdSync);
            end
            if isempty(this.tgn)
                ge = nan(length(this.getIdSync),1, 'single');
            else
                ge = this.tge(this.getIdSync);
            end
            time = this.time.getSubSet(this.getIdSync);            
        end
        
        function [apr_zwd, time] = getAprZwd(this)
            % SYNTAX
            %  [apr_zwd, time] = this.getAprZwd()
            
            apr_zwd = this.apr_zwd(this.getIdSync);
            time = this.time.getEpoch(this.getIdSync);
        end
        
        function [az, el] = getAzEl(this)
            % Get the azimuth and elevation (on valid id_sync)
            %
            % SYNTAX
            %   [az, el] = this.getAzEl();
            az = this.getAz();
            el = this.getEl();
        end
        
        function [az] = getAz(this, go_id)
            % Get the azimuth (on valid id_sync)
            %
            % SYNTAX
            %   az = this.getAzEl();
            if isempty(this.sat.az)
                cc = this.getCC;
                this.sat.az = nan(this.time.length, cc.getNumSat, 'single');
            end
            if nargin < 2
                go_id = 1 : size(this.sat.az, 2);
            end
            
            az = this.sat.az(this.getIdSync, go_id);
        end
        
        function [el] = getEl(this, go_id)
            % Get the azimuth and elevation (on valid id_sync)
            %
            % SYNTAX
            %   el = this.getEl();
            if isempty(this.sat.el)
                cc = this.getCC;
                this.sat.el = nan(this.time.length, cc.getNumSat, 'single');
            end
            if nargin < 2
                go_id = 1 : size(this.sat.el, 2);
            end
            
            el = this.sat.el(this.getIdSync, go_id);
        end
        
        function res = getResidual(this)
            % get residual
            %
            % SYNTAX
            %   res = this.getResidual()
            try
                res = this.sat.res(this.getIdSync(),:);
            catch
                res = [];
            end
        end
        
        function out_prefix = getOutPrefix(this)
            % Get the name for exporting output (valid for dayly output)
            %   - marker name 4ch (from rinex file name)
            %   - 4 char year
            %   - 3 char doy
            %
            % SYNTAX
            time = this.time.getCopy;
            [year, doy] = time.getCentralTime.getDOY();
            out_prefix = sprintf('%s_%04d_%03d_', this.getMarkerName4Ch, year, doy);
        end
        
        function [sys_c, prn] = getSysPrn(this, go_id)
            % Return sys_c and prn for a given go_id
            %
            % SYNTAX
            %    [sys_c, prn] = this.getSysPrn(go_id)
            cc = this.getCC;
            [sys_c, prn] = cc.getSysPrn(go_id);
        end
        
    end
    
    % ==================================================================================================================================================
    %% METHODS UPDATERS
    % ==================================================================================================================================================
    methods
        function updateCoordinates(this)
            % upadte lat lon e ortometric height
            [this.lat, this.lon, this.h_ellips, this.h_ortho] = this.getMedianPosGeodetic();
        end
    end
    
    
    % ==================================================================================================================================================
    %% METHODS IMPORT / EXPORT
    % ==================================================================================================================================================
    
    methods
        function exportTropoSINEX(this, param_to_export)
            % exprot tropspheric product in a sinex file
            %
            % SYNTAX:
            %    exportTropoSinex(this, <param_to_export>)
            if nargin < 2
                param_to_export = [ 1 1 1 0 0 0 0 0];
            end
            for r = 1 : numel(this)
                if min(this(r).quality_info.s0) < 0.10 % If there is at least one good session export the data
                    try
                        rec = this(r);
                        if ~isempty(rec.getZtd)
                            time = this(r).getTime();
                            [year, doy] = this(r).getCentralTime.getDOY();
                            
                            % Detect session length
                            %flag_no_session = (time.last - time.first) > (3600 * 23); % if I have at least 23 hour is a daily session
                            flag_no_session = Core.getState.getSessionDuration() >= 86400;
                            
                            yy = num2str(year);
                            yy = yy(3:4);
                            if flag_no_session
                                sess_str = '0';
                            else
                                datevec_start = datevec(round(time.first.getMatlabTime()*24)/24); 
                                sess_str = 'abcdefghijklmnopqrstuvwx';
                                sess_str = sess_str(datevec_start(4) + 1);
                            end
                            fname = sprintf('%s',[rec.state.getOutDir() filesep rec.parent.getMarkerName4Ch sprintf('%03d', doy) sess_str '.' yy 'zpd']);
                            snx_wrt = SINEX_Writer(fname);
                            snx_wrt.writeTroSinexHeader( rec.time.first, rec.time.getSubSet(rec.time.length), rec.parent.getMarkerName4Ch)
                            snx_wrt.writeFileReference()
                            snx_wrt.writeAcknoledgments()
                            smpl_tropo = median(diff(rec.getIdSync)) * rec.time.getRate;
                            snx_wrt.writeTropoDescription(rec.state.cut_off, rec.time.getRate, smpl_tropo, snx_wrt.SINEX_MAPPING_FLAGS{this.state.mapping_function}, SINEX_Writer.SUPPORTED_PARAMETERS(param_to_export>0), false(length(param_to_export),1))
                            snx_wrt.writeSTACoo( rec.parent.marker_name, rec.xyz(1,1), rec.xyz(1,2), rec.xyz(1,3), 'UNDEF', 'GRD'); % The reference frame depends on the used orbit so it is generraly labled undefined a more intelligent strategy could be implemented
                            snx_wrt.writeTropoSolutionSt()
                            data = [];
                            if param_to_export(1)
                                data = [data rec.ztd(rec.getIdSync,:)*1e3 ];
                            end
                            if param_to_export(2)
                                data = [data rec.tgn(rec.getIdSync,:)*1e3 ];
                            end
                            if param_to_export(3)
                                data = [data rec.tge(rec.getIdSync,:)*1e3];
                            end
                            if param_to_export(4)
                                data = [data rec.getZwd*1e3];
                            end
                            if param_to_export(5)
                                data = [data rec.getPwv*1e3];
                            end
                            [P,T,H] = this.getPTH();
                            if param_to_export(6)
                                data = [data P];
                            end
                            if param_to_export(7)
                                data = [data T];
                            end
                            if param_to_export(8)
                                data = [data H];
                            end
                            snx_wrt.writeTropoSolutionStation(rec.parent.marker_name, rec.time.getEpoch(rec.getIdSync), data, [], param_to_export)
                            snx_wrt.writeTropoSolutionEnd()
                            snx_wrt.writeTroSinexEnd();
                            snx_wrt.close()
                            rec(1).log.addStatusOk(sprintf('Tropo saved into: %s', fname));
                        end
                    catch ex
                        rec(1).log.addError(sprintf('saving Tropo in sinex format failed: %s', ex.message));
                    end
                else
                    if isempty(max(this(r).quality_info.s0))
                        this(1).log.addWarning(sprintf('s02 no solution have been found, station skipped', max(this(r).quality_info.s0)));
                    else
                        this(1).log.addWarning(sprintf('s02 (%f m) too bad, station skipped', max(this(r).quality_info.s0)));
                    end
                end
            end
        end
        
        function exportTropoMat(this)
            % Export the troposphere into a MATLAB data format file
            % The data exported are:
            %  - lat
            %  - lon
            %  - h_ellips
            %  - h_ortho
            %  - ztd
            %  - zwd
            %  - time_utc in matlab format
            %
            % SYNTAX
            %   this.exportTropoMat
            
            for r = 1 : numel(this)
                if max(this(r).quality_info.s0) < 0.10
                    try
                        this(r).updateCoordinates;
                        time = this(r).getTime();
                        state = Core.getCurrentSettings;

                        t_start = state.sss_date_start;
                        t_end = state.sss_date_stop;
                        [year, doy] = t_start.getDOY();
                        t_start_str = t_start.toString('HHMM');
                        t_start.toUtc();
                        time.toUtc();
                        
                        lat = this(r).lat; %#ok<NASGU>
                        lon = this(r).lon; %#ok<NASGU>
                        h_ellips = this(r).h_ellips; %#ok<NASGU>
                        h_ortho = this(r).h_ortho; %#ok<NASGU>
                        ztd = this(r).getZtd(); %#ok<NASGU>
                        zwd = this(r).getZwd(); %#ok<NASGU>
                        utc_time = time.getMatlabTime; %#ok<NASGU>
                        
                        out_dir = fullfile(this(r).state.getOutDir(), sprintf('%4d', year), sprintf('%03d',doy));
                        if ~exist(out_dir, 'file')
                            mkdir(out_dir);
                        end
                        fname_old = sprintf('%s',[out_dir filesep this(r).parent.getMarkerName4Ch sprintf('%04d%03d_%4s_', year, doy, t_start_str) '*.m']);
                        old_file_list = dir(fname_old);

                        fname = sprintf('%s',[out_dir filesep this(r).parent.getMarkerName4Ch sprintf('%04d%03d_%4s_%d', year, doy, t_start_str,  round(t_end-t_start)) '.mat']);
                        save(fname, 'lat', 'lon', 'h_ellips', 'h_ortho', 'ztd', 'zwd', 'utc_time','-v6');
                        
                        this(1).log.addStatusOk(sprintf('Tropo saved into: "%s"', fname));
                        for f = 1 : numel(old_file_list)
                            % If I did not overwrite the old file, delete it
                            if ~strcmp([out_dir filesep old_file_list(f).name], fname)
                                this(1).log.addStatusOk(sprintf('Old tropo file deleted: "%s"', [out_dir filesep old_file_list(f).name]));
                                delete([out_dir filesep old_file_list(f).name]);
                            end
                        end
                    catch ex
                        this(1).log.addError(sprintf('saving Tropo in matlab format failed: "%s"', ex.message));
                    end
                else
                    if isempty(max(this(r).quality_info.s0))
                        this(1).log.addWarning(sprintf('s02 no solution have been found, station skipped'));
                    else
                        this(1).log.addWarning(sprintf('s02 (%f m) too bad, station skipped', max(this(r).quality_info.s0)));
                    end
                end
            end
        end        
        
        function exportTropoCSV(this)
            % Export the troposphere into a MATLAB data format file
            % The data exported are:
            %  - ztd
            %  - zwd
            %  - east gradient
            %  - north gradient
            %  - time_utc in matlab format
            %
            % SYNTAX
            %   this.exportTropoCSV
            
            for r = 1 : numel(this)
                if max(this(r).quality_info.s0) < 0.10
                    try
                        this(r).updateCoordinates;
                        state = Core.getCurrentSettings;
                        time = this(r).getTime();
                        t_start = state.sss_date_start;
                        t_end = state.sss_date_stop;
                        [year, doy] = t_start.getDOY();
                        t_start_str = t_start.toString('HHMM');
                        t_start.toUtc();
                        time.toUtc();
                        ztd = this(r).getZtd();
                        zwd = this(r).getZwd();
                        [gn,ge ] =  this(r).getGradient();
                        
                        out_dir = fullfile(this(r).state.getOutDir(), sprintf('%4d', year), sprintf('%03d',doy));
                        if ~exist(out_dir, 'file')
                            mkdir(out_dir);
                        end
                        
                        fname_old = sprintf('%s',[out_dir filesep this(r).parent.getMarkerName4Ch sprintf('%04d%03d_%4s_', year, doy, t_start_str) '*.csv']);
                        old_file_list = dir(fname_old);

                        fname = sprintf('%s',[out_dir filesep this(r).parent.getMarkerName4Ch sprintf('%04d%03d_%4s_%d', year, doy, t_start_str, round(t_end-t_start)) '.csv']);
                                                
                        fid = fopen(fname,'Wb');
                        n_data = time.length;
                        fprintf(fid,'Date               ,ZTD [m]     ,ZWD [m]     ,GE [m]      ,GN [m]      \n');
                        data = [time.toString('dd/mm/yyyy HH:MM:SS') char(44.*ones(n_data,1)) ...
                            reshape(sprintf('%12.6f',ztd),12,n_data)' char(44.*ones(n_data,1)) ...
                            reshape(sprintf('%12.6f',zwd),12,n_data)' char(44.*ones(n_data,1)) ...
                            reshape(sprintf('%12.6f',ge),12,n_data)' char(44.*ones(n_data,1)) ...
                            reshape(sprintf('%12.6f',gn),12,n_data)' char(10.*ones(n_data,1))];
                        fprintf(fid,data');
                        fclose(fid);                        
                        this(1).log.addStatusOk(sprintf('Tropo saved into: "%s"', fname));
                        for f = 1 : numel(old_file_list)
                            % If I did not overwrite the old file, delete it
                            if ~strcmp([out_dir filesep old_file_list(f).name], fname)
                                this(1).log.addStatusOk(sprintf('Old tropo file deleted: "%s"', [out_dir filesep old_file_list(f).name]));
                                delete([out_dir filesep old_file_list(f).name]);
                            end
                        end
                    catch ex
                        this(1).log.addError(sprintf('saving Tropo in csv format failed: "%s"', ex.message));
                    end
                else
                    if isempty(max(this(r).quality_info.s0))
                        this(1).log.addWarning(sprintf('s02 no solution have been found, station skipped'));
                    else
                        this(1).log.addWarning(sprintf('s02 (%f m) too bad, station skipped', max(this(r).quality_info.s0)));
                    end
                end
            end
        end
        
        function txt = exportWrfLittleR(this, save_on_disk)
            % export WRF-compatible file (LITTLE_R)
            if nargin == 1
                save_on_disk = true;
            end
            if save_on_disk
                [year, doy] = this.time.first.getDOY();
                yy = num2str(year);
                yy = yy(3:4);
                sess_str = '0'; % think how to get the right one from sss_id_list
                fname = sprintf([this.state.getOutDir() '/' this.parent.marker_name '%03d' sess_str '.' yy 'GPSZTD'], doy);
                fid = fopen(fname,'Wb');
            end
            this.updateCoordinates();
            meas_time = this.time.getSubSet(this.id_sync);
            meas_time.toUnixTime();
            txt = '';
            for i = 1 : length(this.id_sync)
                txt = sprintf(['%20.5f%20.5f%40s%40s%40s%40s%20.5f         0         0         0         0         0         F         F         F         0    ' ...
                    '     0%20s-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888' ...
                    '-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888-'...
                    '888888.00000-888888%13.5f      0-888888.00000-888888-888888.00000      0%13.5f      0-888888.00000      0-888888.00000-888888' ...
                    '-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888-888888.00000-888888', ...
                    '-777777.00000      0-777777.00000      0-888888.00000      0-888888.00000      0-888888.00000      0-888888.00000      0' ...
                    '-888888.00000      0-888888.00000      0-888888.00000      0-888888.00000      0\n'],...
                    this.lat, ...
                    this.lon, ...
                    this.parent.marker_name, ...
                    this.parent.marker_type, ...
                    'FM-114 GPSZTD', ...
                    'goGPS software', ...
                    this.h_ortho, ...
                    this.time.toString('yyyymmddHHMMSS'), ...
                    this.ztd(this.id_sync(i))*100, ...
                    this.h_ortho);
            end
            if save_on_disk
                fprintf(fid,'%s', tmp);
                fclose(fid);
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
            if size(this.xyz, 1) > 1
                this.showPositionENU();
                this.showPositionXYZ();
            end
            %this.showMap();
            this.showZtd();
            this.showZtdSlant();
            this.showZtdSlantRes_p();
            this.showResSky_p();
            this.showResSky_c();
            this.showOutliersAndCycleSlip();
            this.showOutliersAndCycleSlip_p();
            dockAllFigures();
        end
        
        function fh_list = showPositionENU(this, flag_one_plot, flag_add_coo)
            % Plot East North Up coordinates of the receiver
            %
            % SYNTAX 
            %   this.plotPositionENU(flag_one_plot, flag_add_coo);
            if nargin == 1 || isempty(flag_one_plot)
                flag_one_plot = false;
            end
            if ~(nargin >= 3 && ~isempty(flag_add_coo) && flag_add_coo > 0)
                flag_add_coo = 0;
            end
            
            fh_list = [];
            for r = 1 : numel(this)
                rec = this(r);
                if ~isempty(rec)
                    xyz = rec.getPosXYZ();
                    if size(xyz, 1) > 1 || flag_add_coo > 0
                        rec(1).log.addMessage('Plotting positions');
                        
                        f = figure('Visible', 'off'); f.Name = sprintf('%03d: PosENU', f.Number); f.NumberTitle = 'off';
                        fh_list = [fh_list; f]; %#ok<AGROW>
                        color_order = handle(gca).ColorOrder;
                        
                        if flag_add_coo == 0
                            xyz = rec.getPosXYZ();
                            
                            t = rec.getPositionTime().getMatlabTime();
                        else
                            xyz = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo.getXYZ;
                            t   = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).time.getMatlabTime();
                            if numel(t) < size(xyz,1)
                                fprintf('%s) There is a problem with add_coo, coordinates are incompatible with their times\n', rec.parent.getMarkerName4Ch);
                                % This is a problem: times and data should have the same dimension
                                xyz = xyz(1:numel(t),:);
                            end
                        end
                        xyz0 = rec.getMedianPosXYZ();
                        
                        fig_name = sprintf('ENU_at%gs_%s_%s', round(median(diff(t * 86400), 'omitnan') * 1e3) / 1e3, rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                        f.UserData = struct('fig_name', fig_name);
                        
                        [enu0(:,1), enu0(:,2), enu0(:,3)] = cart2plan(xyz0(:,1), xyz0(:,2), xyz0(:,3));
                        [enu(:,1), enu(:,2), enu(:,3)] = cart2plan(zero2nan(xyz(:,1)), zero2nan(xyz(:,2)), zero2nan(xyz(:,3)));
                        
                        if ~flag_one_plot, subplot(3,1,1); end
                        e = 1e3 * (enu(:,1) - enu0(1));
                        plot(t, e, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                        ax(3) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        yl = minMax(e);
                        ylim([min(-20, yl(1)) max(20, yl(2))]);
                        setTimeTicks(4); h = ylabel('East [mm]'); h.FontWeight = 'bold';
                        grid on;
                        h = title(sprintf('Position stability of the receiver %s @%gs\n std %.2f [mm]', rec(1).parent.marker_name, round(median(diff(t * 86400), 'omitnan')*1e3) / 1e3, sqrt(var(enu(:,1)*1e3))),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                        if ~flag_one_plot, subplot(3,1,2); end
                        n = 1e3 * (enu(:,2) - enu0(2));
                        plot(t, n, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                        ax(2) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        yl = minMax(n);
                        ylim([min(-20, yl(1)) max(20, yl(2))]);
                        setTimeTicks(4); h = ylabel('North [mm]'); h.FontWeight = 'bold';
                        h = title(sprintf('std %.2f [mm]',sqrt(var(enu(:,2)*1e3))),'interpreter', 'none'); h.FontWeight = 'bold';
                        grid on;
                        if ~flag_one_plot, subplot(3,1,3); end
                        up = 1e3 * (enu(:,3) - enu0(3));
                        plot(t, up, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                        ax(1) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        yl = minMax(up);
                        ylim([min(-20, yl(1)) max(20, yl(2))]);
                        setTimeTicks(4); h = ylabel('Up [mm]'); h.FontWeight = 'bold';
                        h = title(sprintf('std %.2f [mm]',sqrt(var(enu(:,3)*1e3))),'interpreter', 'none'); h.FontWeight = 'bold';
                        grid on;
                        if flag_one_plot
                            h = ylabel('ENU [cm]'); h.FontWeight = 'bold';
                        else
                            linkaxes(ax, 'x');
                        end                        
                        grid on;
                        Core_UI.beautifyFig(f);
                        Core_UI.addBeautifyMenu(f);
                        Core_UI.addExportMenu(f);
                        f.Visible = 'on'; drawnow;
                    else
                        rec(1).log.addWarning(sprintf('%s - Plotting a single point static position is not yet supported', rec.parent.getMarkerName4Ch));
                    end
                end
            end
        end

        function fh_list = showPositionPlanarUp(this, flag_add_coo)
            % Plot East North Up coordinates of the receiver
            %
            % SYNTAX 
            %   this.plotPositionENU(flag_one_plot);            
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo > 0)
                flag_add_coo = 0;
            end
            fh_list = [];
            for r = 1 : numel(this)
                rec = this(r);
                if ~isempty(rec)
                    xyz = rec.getPosXYZ();
                    if size(xyz, 1) > 1 || flag_add_coo > 0                        
                        rec(1).log.addMessage('Plotting positions');
                        
                        f = figure('Visible', 'off'); f.Name = sprintf('%03d: PosENU', f.Number); f.NumberTitle = 'off';
                        
                        fh_list = [fh_list; f]; %#ok<AGROW>
                        color_order = handle(gca).ColorOrder;
                        
                        if flag_add_coo == 0
                            xyz = rec.getPosXYZ();
                            
                            t = rec.getPositionTime().getMatlabTime();
                        else
                            xyz = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo.getXYZ;
                            t   = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).time.getMatlabTime();
                            if numel(t) < size(xyz,1)
                                fprintf('%s) There is a problem with add_coo, they are incompatible with their times\n', rec.parent.getMarkerName4Ch);
                                % This is a problem: times and data should have the same dimension
                                xyz = xyz(1:numel(t),:);
                            end
                        end                        
                        xyz0 = rec.getMedianPosXYZ();
                        
                        fig_name = sprintf('EN_U_at%gs_%s_%s', round(median(diff(t * 86400), 'omitnan') * 1e3) / 1e3, rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                        f.UserData = struct('fig_name', fig_name);

                        [enu0(:,1), enu0(:,2), enu0(:,3)] = cart2plan(xyz0(:,1), xyz0(:,2), xyz0(:,3));
                        [enu(:,1), enu(:,2), enu(:,3)] = cart2plan(zero2nan(xyz(:,1)), zero2nan(xyz(:,2)), zero2nan(xyz(:,3)));
                        
                        main_vb = uix.VBox('Parent', f, ...
                            'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                        
                        tmp_box1 = uix.VBox('Parent', main_vb, ...
                            'Padding', 5, ...
                            'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                        tmp_box2 = uix.VBox('Parent', main_vb, ...
                            'Padding', 5, ...
                            'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                        main_vb.Heights = [-2 -1];
                        Core_UI.beautifyFig(f);
                        f.Visible = 'on';
                        drawnow
                        f.Visible = 'off';
                        ax = axes('Parent', tmp_box1);

                        % Plot parallel
                        max_e = ceil(max(abs(1e3 * minMax(enu(:,1) - enu0(1))))/5) * 5;
                        max_n = ceil(max(abs(1e3 * minMax(enu(:,2) - enu0(2))))/5) * 5;
                        max_r = ceil(sqrt(max_e^2 + max_n^2) / 5) * 5;
                       
                        % Plot circles of precision
                        az_l = 0 : pi/200: 2*pi;
                        % dashed
                        id_dashed = serialize(bsxfun(@plus, repmat((0:20:395)',1,5), (1:5)));
                        az_l(id_dashed) = nan;
                        decl_s = ((10 : 10 : max_r));
                        for d = decl_s
                            x = cos(az_l).*d;
                            y = sin(az_l).*d;
                            plot(x,y,'color',[0.6 0.6 0.6], 'LineWidth', 2); hold on;
                            x = cos(az_l).*(d-5);
                            y = sin(az_l).*(d-5);
                            plot(x,y,'color',[0.75 0.75 0.75], 'LineWidth', 2); hold on;
                        end
                        
                        plot((enu(:,1) - enu0(1)) * 1e3, (enu(:,2) - enu0(2)) * 1e3, 'o', 'MarkerSize', 4, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                        axis equal;
                        h = ylabel('East [mm]'); h.FontWeight = 'bold';
                        h = xlabel('North [mm]'); h.FontWeight = 'bold';
                        ylim(max_r * [-1 1]);
                        xlim(max_r * [-1 1]);
                        grid on;
                        h = title(sprintf('Position Stability %s @%gs\nstd E %.2f mm - N %.2f mm\\fontsize{5} \n', strrep(rec.parent.getMarkerName4Ch, '_', '\_'), round(median(diff(t * 86400), 'omitnan')*1e3) / 1e3, std((enu(:,1) - enu0(1)) * 1e3, 'omitnan'), std((enu(:,2) - enu0(2)) * 1e3, 'omitnan')), 'FontName', 'Open Sans'); 
                        h.FontWeight = 'bold';
                        
                        ax = axes('Parent', tmp_box2);
                                   
                        up = 1e3 * (enu(:,3) - enu0(3));
                        plot(t, up, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                        ax(1) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        yl = minMax(up);
                        ylim([min(-20, yl(1)) max(20, yl(2))]);
                        setTimeTicks(4); h = ylabel('Up [cm]'); h.FontWeight = 'bold';
                        h = title(sprintf('Up std %.2f [mm]', std(enu(:,3)*1e3)), 'interpreter', 'none'); h.FontWeight = 'bold';
                        grid on;
                        Core_UI.beautifyFig(f);
                        Core_UI.addBeautifyMenu(f);
                        Core_UI.addExportMenu(f);
                        f.Visible = 'on'; drawnow;
                    else
                        rec(1).log.addMessage('Plotting a single point static position is not yet supported');
                    end
                end
            end
        end

        function fh_list = showPositionXYZ(this, flag_one_plot, flag_add_coo)
            % Plot X Y Z coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionXYZ();
            if nargin == 1 || isempty(flag_one_plot)
                flag_one_plot = false;
            end
            if ~(nargin >= 3 && ~isempty(flag_add_coo) && flag_add_coo > 0)
                flag_add_coo = 0;
            end
            
            fh_list = [];
            for r = 1 : numel(this)
                rec = this(r);
                if ~isempty(rec)
                    xyz = rec.getPosXYZ();
                    if size(xyz, 1) > 1 || flag_add_coo > 0
                        rec(1).log.addMessage('Plotting XYZ positions');
                        
                        f = figure('Visible', 'off'); f.Name = sprintf('%03d: PosXYZ', f.Number); f.NumberTitle = 'off';
                        fh_list = [fh_list; f]; %#ok<AGROW>
                        color_order = handle(gca).ColorOrder;
                        
                        if flag_add_coo == 0
                            xyz = rec.getPosXYZ();
                            
                            t = rec.getPositionTime().getMatlabTime();
                        else
                            xyz = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo.getXYZ;
                            t   = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).time.getMatlabTime();
                            if numel(t) < size(xyz,1)
                                fprintf('%s) There is a problem with add_coo, they are incompatible with their times\n', rec.parent.getMarkerName4Ch);
                                % This is a problem: times and data should have the same dimension
                                xyz = xyz(1:numel(t),:);
                            end
                        end
                        fig_name = sprintf('XYZ_at%gs_%s_%s', round(median(diff(t * 86400), 'omitnan') * 1e3) / 1e3, rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                        f.UserData = struct('fig_name', fig_name);
                        
                        xyz0 = rec.getMedianPosXYZ();
                        
                        x = 1e2 * bsxfun(@minus, zero2nan(xyz(:,1)), xyz0(1));
                        y = 1e2 * bsxfun(@minus, zero2nan(xyz(:,2)), xyz0(2));
                        z = 1e2 * bsxfun(@minus, zero2nan(xyz(:,3)), xyz0(3));
                        
                        if ~flag_one_plot, subplot(3,1,1); end
                        plot(t, x, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                        ax(3) = gca(); xlim([t(1) t(end)]); setTimeTicks(4); h = ylabel('X [cm]'); h.FontWeight = 'bold';
                        grid on;
                        h = title(sprintf('Position stability of the receiver %s @%gs\n std %.2f [cm]', rec(1).parent.marker_name, round(median(diff(t * 86400), 'omitnan') * 1e3) / 1e3, sqrt(var(x))),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                        if ~flag_one_plot, subplot(3,1,2); end
                        plot(t, y, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                        ax(2) = gca(); xlim([t(1) t(end)]); setTimeTicks(4); h = ylabel('Y [cm]'); h.FontWeight = 'bold';
                        grid on;
                        h = title(sprintf('std %.2f [cm]',sqrt(var(y))),'interpreter', 'none'); h.FontWeight = 'bold';
                        if ~flag_one_plot, subplot(3,1,3); end
                        plot(t, z, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                        ax(1) = gca(); xlim([t(1) t(end)]); setTimeTicks(4); h = ylabel('Z [cm]'); h.FontWeight = 'bold';
                        grid on;
                        if flag_one_plot
                            h = ylabel('XYZ [m]'); h.FontWeight = 'bold';
                        end
                        linkaxes(ax, 'x');                        
                        Core_UI.beautifyFig(f);
                        Core_UI.addBeautifyMenu(f);
                        f.Visible = 'on'; drawnow;
                        h = title(sprintf('std %.2f [cm]',sqrt(var(z))),'interpreter', 'none'); h.FontWeight = 'bold';

                    else
                        rec.log.addMessage('Plotting a single point static position is not yet supported');
                    end
                end
            end
        end
        
        function fh_list = showPositionSigmas(this, one_plot)
            % Show Sigmas of the solutions
            %
            % SYNTAX
            %   this.showPositionSigmas();
            
            if nargin == 1
                one_plot = false;
            end
            
            fh_list = [];
            rec = this;
            if ~isempty(rec)
                xyz = rec(1).getPosXYZ();
                if size(xyz, 1) > 1
                    rec(1).log.addMessage('Plotting ENU sigmas');
                    
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: sigma processing', f.Number); f.NumberTitle = 'off';
                    fh_list = [fh_list; f]; %#ok<AGROW>
                    fig_name = sprintf('ENU_s0_%s_%s', rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                    f.UserData = struct('fig_name', fig_name);
                    color_order = handle(gca).ColorOrder;
                    
                    s0 = rec.quality_info.s0;
                    s0_ip = rec.quality_info.s0_ip;
                    
                    t = rec.getPositionTime().getMatlabTime;
                    
                    if ~one_plot, subplot(2,1,2); end
                    plot(t, s0 * 1e2, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                    ax(2) = gca(); xlim([t(1) t(end)]); setTimeTicks(4); h = ylabel('s0 [cm]'); h.FontWeight = 'bold';
                    grid on;
                    if ~one_plot, subplot(2,1,1); end
                    plot(t, s0_ip * 1e2, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                    ax(1) = gca(); xlim([t(1) t(end)]); setTimeTicks(4); h = ylabel('s0 ip [cm]'); h.FontWeight = 'bold';
                    h = title(sprintf('Receiver %s', rec(1).parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                    grid on;
                    if one_plot
                        h = ylabel('Sigmas of the processing [cm]'); h.FontWeight = 'bold';
                    end
                    linkaxes(ax, 'x');
                    Core_UI.beautifyFig(f);
                    Core_UI.addBeautifyMenu(f);
                    f.Visible = 'on'; drawnow;
                else
                    rec.log.addMessage('Plotting a single point static position is not yet supported');
                end
            end
        end
                
        function fh_list = showMap(sta_list, new_fig)
            if nargin < 2
                new_fig = true;
            end
            if new_fig
                f = figure('Visible', 'off');
            else
                f = gcf;
                hold on;
            end
            fh_list = f;
            fig_name = sprintf('RecMap');
            f.UserData = struct('fig_name', fig_name);
            
            maximizeFig(f);
            [lat, lon] = sta_list.getMedianPosGeodetic();
            
            plot(lon(:), lat(:),'.k', 'MarkerSize', 5); hold on;
            % Label BG (in background w.r.t. the point)
            for r = 1 :  numel(sta_list)
                text(lon(r), lat(r), '               ', ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment','left');
            end
            
            for r = 1 : numel(sta_list)
                plot(lon(r), lat(r), '.', 'MarkerSize', 40, 'Color', Core_UI.getColor(r, numel(sta_list)));
            end
            plot(lon(:), lat(:), '.k', 'MarkerSize', 5);
            plot(lon(:), lat(:), 'ko', 'MarkerSize', 13, 'LineWidth', 2);
            
            if numel(sta_list) == 1
                lon_lim = minMax(lon);
                lat_lim = minMax(lat);
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
                name = upper(sta_list(r).parent.getMarkerName4Ch());
                text(lon(r), lat(r), ['   ' name], ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment','left');
            end
            
            plot_google_map('alpha', 0.95, 'MapType', 'satellite');
            title('Receiver position');
            xlabel('Longitude [deg]');
            ylabel('Latitude [deg]');
            Core_UI.beautifyFig(f);            
            Core_UI.addBeautifyMenu(f);
            f.Visible = 'on'; drawnow;            
        end
        
        
        function fh_list = showResSky_p(this, sys_c_list)
            % Plot residuals of the solution on polar scatter
            % SYNTAX this.plotResSkyPolar(sys_c)
            
            fh_list = [];
            cc = this.getCC;
            if isempty(this.sat.res)
                this.log.addWarning('Residuals have not been computed');
            else
                if nargin == 1 || isempty(sys_c_list)
                    sys_c_list = unique(cc.system);
                end
                
                for sys_c = sys_c_list
                    s = cc.getGoIds(sys_c);%this.go_id(this.system == sys_c);
                    res = abs(this.sat.res(:, s)) * 1e2;
                    
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s Res SkyC %s', f.Number, this.parent.getMarkerName4Ch, cc.getSysName(sys_c)); f.NumberTitle = 'off';

                    fh_list = [fh_list; f]; %#ok<AGROW>
                    fig_name = sprintf('SNR_polar_%s_%s_%s', this.parent.getMarkerName4Ch, cc.getSysName(sys_c), this.time.first.toString('yyyymmdd_HHMM'));
                    f.UserData = struct('fig_name', fig_name);
                    
                    id_ok = (res~=0);
                    az = this.sat.az(:, s);
                    el = this.sat.el(:, s);
                    polarScatter(serialize(az(id_ok))/180*pi,serialize(90-el(id_ok))/180*pi, 45, serialize(res(id_ok)), 'filled');
                    caxis([min(abs(this.sat.res(:)*1e2)) min(20, min(6*std(zero2nan(this.sat.res(:)*1e2),'omitnan'), max(abs(zero2nan(this.sat.res(:)*1e2)))))]);
                    colormap(flipud(hot)); f.Color = [.95 .95 .95]; cb = colorbar(); cbt = title(cb, '[cm]'); cbt.Parent.UserData = cbt;
                    h = title(sprintf('Satellites residuals - receiver %s - %s', this.parent.marker_name, cc.getSysExtName(sys_c)),'interpreter', 'none');  h.FontWeight = 'bold';
                    Core_UI.beautifyFig(f);
                    Core_UI.addBeautifyMenu(f);
                    f.Visible = 'on'; drawnow;
                end
            end
        end
        
        function fh_list = showResSky_c(this, sys_c_list)
            % Plot residuals of the solution on cartesian axes
            % SYNTAX this.plotResSkyCart()
            fh_list = [];
            cc = this.getCC;
            if isempty(this.sat.res)
                this.log.addWarning('Residuals have not been computed');
            else
                if nargin == 1 || isempty(sys_c_list)
                    sys_c_list = unique(cc.system);
                end
                
                for sys_c = sys_c_list
                    s  = cc.getGoIds(sys_c); %unique(this.go_id(this.system == sys_c));
                    res = abs(this.sat.res(:, s)) * 1e2;
                    
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s Res SkyC %s', f.Number, this.parent.getMarkerName4Ch, cc.getSysName(sys_c)); f.NumberTitle = 'off';
                    
                    fh_list = [fh_list; f]; %#ok<AGROW>
                    fig_name = sprintf('SNR_cartesian_%s_%s_%s', this.parent.getMarkerName4Ch, cc.getSysName(sys_c), this.time.first.toString('yyyymmdd_HHMM'));
                    f.UserData = struct('fig_name', fig_name);
                    
                    % this.updateAzimuthElevation()
                    id_ok = (res~=0);
                    az = this.sat.az(:, s);
                    el = this.sat.el(:, s);
                    scatter(serialize(az(id_ok)),serialize(el(id_ok)), 45, serialize(res(id_ok)), 'filled');
                    caxis([min(abs(this.sat.res(:)*1e2)) min(20, min(6*std(zero2nan(this.sat.res(:)*1e2),'omitnan'), max(abs(zero2nan(this.sat.res(:)*1e2)))))]);
                    colormap(flipud(hot)); f.Color = [.95 .95 .95]; cb = colorbar(); cbt = title(cb, '[cm]'); cbt.Parent.UserData = cbt; ax = gca; ax.Color = 'none';
                    h = title(sprintf('Satellites residuals - receiver %s - %s', this.parent.marker_name, cc.getSysExtName(sys_c)),'interpreter', 'none');  h.FontWeight = 'bold';
                    hl = xlabel('Azimuth [deg]'); hl.FontWeight = 'bold';
                    hl = ylabel('Elevation [deg]'); hl.FontWeight = 'bold';
                    Core_UI.beautifyFig(f);                    
                    Core_UI.addBeautifyMenu(f);
                    f.Visible = 'on'; drawnow;                    
                end
            end
        end
        
        function fh_list = showRes(sta_list)
            % In a future I could use multiple tabs for each constellation
            fh_list = [];
            for r = numel(sta_list)
                cc = sta_list(r).getCC;
                work = sta_list(r);
                if ~work.isEmpty
                    
                    win = figure('Visible', 'on', ...
                        'NumberTitle', 'off', ...
                        'units', 'normalized', ...
                        'outerposition', [0 0 1 1]); % create maximized figure
                    win.Name = sprintf('%03d - %s residuals', win.Number, work.parent.getMarkerName4Ch);
                    
                    fh_list = [fh_list; win]; %#ok<AGROW>
                    fig_name = sprintf('ENU_s0_%s_%s', work.parent.getMarkerName4Ch, work.time.first.toString('yyyymmdd_HHMM'));
                    win.UserData = struct('fig_name', fig_name);
                    
                    v_main = uix.VBoxFlex('Parent', win, ...
                        'Spacing', 5);
                    
                    % Axe with all the satellites
                    overview_box = uix.VBoxFlex('Parent', v_main, ...
                        'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                    
                    ax_all = axes('Parent', overview_box, 'Units', 'normalized');
                    
                    % Single sat axes
                    n_sat = work.getMaxSat;
                    
                    sat_box = uix.VBoxFlex('Parent', v_main, ...
                        'Padding', 5, ...
                        'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                    
                    v_main.Heights = [-2 -5];
                    
                    scroller = uix.ScrollingPanel('Parent', sat_box);
                    sat_grid = uix.Grid('Parent', scroller, ...
                        'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                    scroller.Heights = 120 * ceil(n_sat / 4);
                    for s = 1 : n_sat
                        single_sat(s) = uix.VBox('Parent', sat_grid, ...
                            'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                        uicontrol('Parent', single_sat(s), ...
                            'Style', 'Text', ...
                            'String', sprintf('Satellite %s', cc.getSatName(s)), ...
                            'ForegroundColor', Core_UI.BLACK, ...
                            'HorizontalAlignment', 'center', ...
                            'FontSize', Core_UI.getFontSize(7), ...
                            'FontWeight', 'Bold', ...
                            'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                        ax_sat(s) = axes('Parent', single_sat(s));
                    end
                    sat_grid.Heights = -ones(1, ceil(n_sat / 4));
                    for s = 1 : n_sat
                        single_sat(s).Heights = [18, -1];
                        drawnow
                    end
                    
                    %% fill the axes
                    win.Visible = 'on'; drawnow;
                    colors = Core_UI.getColor(1 : n_sat, n_sat);
                    ax_all.ColorOrder = colors; hold(ax_all, 'on');
                    plot(ax_all, work.time.getMatlabTime, zero2nan(work.sat.res)*1e3, '.-');
                    setTimeTicks(ax_all, 3, 'auto');
                    drawnow; ylabel('residuals [mm]'); grid on;
                    
                    id_ok = false(n_sat, 1);
                    for s = 1 : n_sat
                        res = zero2nan(work.sat.res(:,s))*1e3;
                        id_ok(s) = any(res);
                        if id_ok(s)
                            plot(ax_sat(s), work.time.getMatlabTime, res, '-', 'LineWidth', 2, 'Color', colors(s, :));
                            grid on; ax_sat(s).YMinorGrid = 'on'; ax_sat(s).XMinorGrid = 'on';
                            %setTimeTicks(ax_sat(s), 2,'HH:MM');
                        else
                            single_sat(s).Visible = 'off';
                        end
                    end
                    linkaxes([ ax_all,ax_sat(id_ok)]);
                    drawnow
                end
            end            
        end        
                                
        function ant_mp = computeMultiPath(this, type, l_max, flag_reg)
            % Get Zernike multi pth coefficients
            %
            % INPUT
            %   type    can be:
            %            'pr'   -> Uncombined pseudo-ranges residuals
            %            'ph'   -> Uncombined carrier-phase residuals
            %   l_max   maximum degree for of the Zernike polynomials
            %
            % SYNTAX
            %   this.computeMultiPath(type, l_max, flag_mask_reg)
            
            flag_debug = false;
            grid_step = 0.5;
            if nargin < 3 || isempty(l_max)
                l_max = [43 43 43];
            end
            if numel(l_max) == 1
                l_max = [l_max l_max l_max];
            end
            
            if nargin < 4 || isempty(flag_reg)
                flag_reg = true;
            end
            
            deg2rad = pi/180;
            
            cc = Core.getState.getConstellationCollector;
            sys_c_list = cc.getAvailableSys;
            if nargin < 2 || isempty(type)
                type = 'ph';
            end
            
            if type(2) == 'r'
                res = this.sat.res_pr_by_pr;
            else
                res = this.sat.res_ph_by_ph;
            end
            
            if type(2) == 'r'
                name = 'Uncombined pseudo-ranges residuals';
                id_obs = this.obs_code(:,1) == 'C';
            else
                name = 'Uncombined carrier-phase residuals';
                id_obs = this.obs_code(:,1) == 'L';
            end
            
            id_obs = find(id_obs);
            
            log = Core.getLogger;
            log.addMarkedMessage(sprintf('Computing multipath mitigation coefficients for "%s"', this.parent.getMarkerName));
            ant_mp = struct();
            if isempty(res)
                log.addError(sprintf('No %s residuals found in %s', name, this.parent.getMarkerName4Ch));
            else
                ss_ok = intersect(unique(this.system(id_obs)), sys_c_list);
                % for each satellite system
                for sys_c = ss_ok
                    id_obs_sys = id_obs(this.system(id_obs) == sys_c);
                    trk_ok = Core_Utils.unique2ch(this.obs_code(id_obs_sys,1:2));
                    % for each tracking
                    for t = 1 : size(trk_ok, 1)
                        cur_res_id = (this.system(id_obs)' == sys_c) & (this.obs_code(id_obs, 2) == trk_ok(t,2));
                        [~, prn] = cc.getSysPrn(this.go_id(id_obs(cur_res_id)));
                                                                                                
                        data_found = false;

                        cur_res = res(:, cur_res_id);
                        %res_tmp = Receiver_Commons.smoothMat(res_tmp, 'spline', 10);
                        
                        % Zernike filter
                        m_max = l_max;
                        
                        az_all = [];
                        el_all = [];
                        res_all = [];
                        for s = 1 : sum(cur_res_id)
                            id_ok = find(~isnan(zero2nan(cur_res(:, s))));
                            if any(id_ok)
                                data_found = true;
                                res_all = [res_all; cur_res(id_ok, s)];
                                
                                go_id = cc.getIndex(sys_c, prn(s));
                                az_all = [az_all; this.sat.az(id_ok, go_id) .* deg2rad];
                                el_all = [el_all; (this.sat.el(id_ok, go_id)) .* deg2rad];
                            end
                        end
                        if data_found                            
                            log.addMessage(log.indent(sprintf('1. preparing data for %s %s', sys_c, trk_ok(t,:)), 9));
                            az_all = [];
                            el_all = [];
                            res_all = cur_res(~isnan(zero2nan(cur_res(:))));
                            res_smt = Receiver_Commons.smoothMat(cur_res, 'spline', 300/this.getTime.getRate);
                            res_smt = res_smt(~isnan(zero2nan(cur_res(:))));                            
                            % az and el must be retrieved satellite by satellite
                            for s = 1 : sum(cur_res_id)
                                id_ok = find(~isnan(zero2nan(cur_res(:, s))));
                                if any(id_ok)
                                    data_found = true;
                                    
                                    go_id = cc.getIndex(sys_c, prn(s));
                                    az_all = [az_all; this.sat.az(id_ok, go_id) .* deg2rad];
                                    el_all = [el_all; (this.sat.el(id_ok, go_id)) .* deg2rad];
                                end
                            end
                            
                            % Remove outliers
                            id_ok = Core_Utils.polarCleaner(az_all, el_all, res_all, [360, 1]) & Core_Utils.polarCleaner(az_all, el_all, res_smt, [360, 1]);
                            log.addMessage(log.indent(sprintf('2. Outlier rejection (%.3f%%)', (sum(~id_ok) / numel(id_ok)) * 100), 9));
                            if flag_debug
                                figure; plot(el_all/pi*180, res_all*1e3, '.'); hold on; plot(el_all(~id_ok)/pi*180, res_all(~id_ok)*1e3, 'o');
                                legend('residuals', 'outliers');
                                title((sprintf('Residuals of  %s %s%s [mm]', this.parent.getMarkerName4Ch, sys_c, trk_ok(t,:)))); drawnow
                                grid on;
                            end
                            clear res_smt;
                            az_all = az_all(id_ok);
                            el_all = el_all(id_ok);
                            res_all = res_all(id_ok);
                            n_obs = numel(res_all);
                            
                            % 3 sigma filter per latitude
                            if flag_reg
                                log.addMessage(log.indent('3. Preparing regularization', 9));
                                % Get regularization points based on empty sky areas
                                [data_map, n_data_map, az_grid, el_grid] = Core_Utils.polarGridder(az_all, el_all, res_all, [1 1]);
                                [az_grid, el_grid] = meshgrid(az_grid, el_grid);
                                az_reg = az_grid(n_data_map <= 0);
                                el_reg = el_grid(n_data_map <= 0);
                                
                                [~, z_map, az_grid, el_grid] = Core_Utils.polarGridder(az_all, el_all, res_all, 0.5);
                                az_all = [az_all; az_reg];
                                el_all = [el_all; el_reg];
                                res_all = [res_all; zeros(size(el_reg))];
                                
                                % Add additional points at the board
                                for i = 0 : 0.5 : (Core.getState.getCutOff - 3)
                                    az_all = [az_all; (-pi : 0.05 : pi)'];
                                    el_all = [el_all; i/180*pi + (-pi : 0.05 : pi)'*0];
                                    res_all = [res_all; (-pi : 0.05 : pi)'*0];
                                end
                            else
                                az_grid = ((-180 + (grid_step(1) / 2)) : grid_step(1) : (180 - grid_step(1) / 2)) .* (pi/180);
                                el_grid = flipud(((grid_step(end) / 2) : grid_step(end) : 90 - (grid_step(end) / 2))' .* (pi/180));
                            end
                            
                            log.addMessage(log.indent(sprintf('%d. Zernike coef. estimation (l_max = %d) (1/3)', 3 + flag_reg*1, l_max(1)), 9));
                            %log.addMessage(log.indent('mapping r with m1 = pi/2 * (1-cos(el))', 12));
                            
                            % Use two
                            el2radius1 = @(el) cos(el).^2;
                            el2radius2 = @(el) sin(pi/2*cos(el).^2);
                            el2radius3 = @(el) sin(pi/2*cos(el));
                            %omf = @(el) 1 - sin(el);
                            %el2radius3 = @(el) omf(pi/2 * (1-cos(pi/2*(1-cos(el)))));
                            res_work = res_all;
                            
                            [res_filt, z_par1, l1, m1] = Core_Utils.zFilter(l_max(1), m_max(1), az_all, el2radius1(el_all), res_work, 1e-5);
                            res_work = res_work - res_filt;
                            
                            if l_max(2) > 0
                                log.addMessage(log.indent(sprintf('%d. Zernike coef. estimation (l_max = %d) (2/3)', 4 + flag_reg*1, l_max(2)), 9));
                                %log.addMessage(log.indent('mapping r with m2 = pi/2 * (1-cos(m1))', 12));
                                [res_filt, z_par2, l2, m2] = Core_Utils.zFilter(l_max(2), m_max(2), az_all, el2radius2(el_all), res_work, 1e-5);
                                res_work = res_work - res_filt;
                            end    
                            
                            if l_max(3) > 0
                                log.addMessage(log.indent(sprintf('%d. Zernike coef. estimation (l_max = %d) (3/3)', 4 + flag_reg*1, l_max(3)), 9));
                                [res_filt, z_par3, l3, m3] = Core_Utils.zFilter(l_max(3), m_max(3), az_all, el2radius3(el_all), res_work, 1e-5);
                                res_work = res_work - res_filt;
                            end

                            % Generate maps
                            log.addMessage(log.indent(sprintf('%d. Compute mitigation grids', 5 + flag_reg*1), 9));
                            [az_mgrid, el_mgrid] = meshgrid(az_grid, el_grid);
                            [z_map1] = Core_Utils.zSinthesys(l1, m1, az_mgrid, el2radius1(el_mgrid), z_par1);
                            if l_max(2) <= 0
                                z_map2 = 0;
                            else
                                [z_map2] = Core_Utils.zSinthesys(l2, m2, az_mgrid, el2radius2(el_mgrid), z_par2);
                            end
                            if l_max(3) <= 0
                                z_map3 = 0;
                            else
                                [z_map3] = Core_Utils.zSinthesys(l3, m3, az_mgrid, el2radius3(el_mgrid), z_par3);
                            end
                            z_map = z_map1 + z_map2 + z_map3;
                            res_work((n_obs + 1) : end) = 0; % Restore regularization to zero
                            [g_map, n_map, az_grid, el_grid] = Core_Utils.polarGridder(az_all, el_all, res_work, [4 1], grid_step);
                            
                            if flag_debug
                                %figure; imagesc(1e3*(z_map)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map1)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                title((sprintf('Zernike expansion (1) of %s %s%s [mm]', this.parent.getMarkerName4Ch, sys_c, trk_ok(t,:)))); drawnow

                                figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map2)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                title((sprintf('Zernike expansion (2) of %s %s%s [mm]', this.parent.getMarkerName4Ch, sys_c, trk_ok(t,:)))); drawnow

                                figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map3)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                title((sprintf('Zernike expansion (3) of %s %s%s [mm]', this.parent.getMarkerName4Ch, sys_c, trk_ok(t,:)))); drawnow

                                %zmap2scatter = griddedInterpolant(flipud([az_mgrid(:,end) - 2*pi, az_mgrid, az_mgrid(:,1) + 2*pi])', flipud([el_mgrid(:,end) el_mgrid el_mgrid(:,1)])', flipud([z_map(:,end) z_map z_map(:,1)])', 'linear');                                
                                %[g_map, n_map, az_grid_tmp, el_grid_tmp] = Core_Utils.polarGridder(az_all, el_all, zmap2scatter(az_all, el_all), [5 1]);
                                figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(g_map)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                title((sprintf('Gridder residuals of %s %s%s [mm]', this.parent.getMarkerName4Ch, sys_c, trk_ok(t,:)))); drawnow

                                figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                title((sprintf('Zernike expansion of %s %s%s [mm]', this.parent.getMarkerName4Ch, sys_c, trk_ok(t,:)))); drawnow

                                figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(z_map + g_map)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                title((sprintf('Final map of %s %s%s [mm]', this.parent.getMarkerName4Ch, sys_c, trk_ok(t,:)))); drawnow

                                [tmp] = Core_Utils.polarGridder(az_all, el_all, res_all, [4 1], grid_step);
                                figure; polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(tmp)); colormap((Cmap.get('PuOr', 2^11))); caxis([-5 5]); colorbar;
                                title((sprintf('Final map of %s %s%s [mm]', this.parent.getMarkerName4Ch, sys_c, trk_ok(t,:)))); drawnow
                            end
                            
                            if ~isfield(ant_mp, sys_c)
                                ant_mp.(sys_c) = struct;
                            end
                            if ~isfield(ant_mp.(sys_c), trk_ok(t,:))
                                ant_mp.(sys_c).(trk_ok(t,:)) = struct;
                            end
                            % Keep multiple solutions in the struct
                            % decide a-posteriori what it's better
                            
                            % Save grids of multi-path
                            ant_mp.(sys_c).(trk_ok(t,:)).az_grid = single(az_grid);
                            ant_mp.(sys_c).(trk_ok(t,:)).el_grid = single(el_grid);
                            ant_mp.(sys_c).(trk_ok(t,:)).z_map = single(z_map);             % Zernike math
                            ant_mp.(sys_c).(trk_ok(t,:)).g_map = single(z_map + g_map);     % Second gridding step
                            % Save Zernike coefficients
                            %ant_mp.(sys_c).(trk_ok(t,:)).z_par = [z_par1 z_par2];
                            %ant_mp.(sys_c).(trk_ok(t,:)).l = l;
                            %ant_mp.(sys_c).(trk_ok(t,:)).m = m;
                        end
                        if ~data_found
                            log = Core.getLogger;
                            log.addError(sprintf('No %s %s found in %s for constellation %s', name, trk_ok(t,:), this.parent.getMarkerName4Ch, cc.getSysName(sys_c)));
                        end
                    end
                end
                ant_mp.time_lim = this.time.getEpoch([1 this.time.length]);                
            end
        end

        function fh_list = showResScatter(this, sys_c_list, type, res)
            % Plot the residuals of phase per Satellite
            %
            % INPUT
            %   type    can be:
            %            'co'   -> Combined residuals (one set for each satellite) DEFAULT
            %            'pr'   -> Uncombined pseudo-ranges residuals
            %            'ph'   -> Uncombined carrier-phase residuals
            %   res     is the matrix of residuals satellite by satellite and can be passed from e.g. NET
            %
            % SYNTAX
            %   this.showResScatter(sys_c_list, type, res)
            
            fh_list = [];
            cc = Core.getState.getConstellationCollector;
            if nargin < 2 || isempty(sys_c_list)
                sys_c_list = cc.getAvailableSys;
            end
            if nargin < 3 || isempty(type)
                type = 'ph';
            end
            
            if nargin < 4 || isempty(res)                
                if type(2) == 'o'
                    res = this.sat.res;
                elseif type(2) == 'r'
                    res = this.sat.res_pr_by_pr;
                else
                    res = this.sat.res_ph_by_ph;
                end
            end
            if type(2) == 'o'
                name = 'Combined residuals sat. by sat.';
                scale = 1e3; % mm
                id_obs = 1 : size(res, 2);
            elseif type(2) == 'r'
                name = 'Uncombined pseudo-ranges residuals';
                id_obs = this.obs_code(:,1) == 'C';
                scale = 1e2; % cm
            else
                name = 'Uncombined carrier-phase residuals';
                id_obs = this.obs_code(:,1) == 'L';
                scale = 1e3; % mm
            end
            id_obs = find(id_obs);
            if isempty(res)
                log = Core.getLogger;
                log.addError(sprintf('No %s residuals found in %s', name, this.parent.getMarkerName4Ch));
            else
                if type(2) == 'o'
                    ss_ok = intersect(cc.sys_c, sys_c_list);
                else
                    ss_ok = intersect(unique(this.system(id_obs)), sys_c_list);                    
                end
                for sys_c = ss_ok
                    if type(2) == 'o'
                        trk_ok = 'COM'; % combined tracking -> only this is existing
                    else
                        id_obs_sys = id_obs(this.system(id_obs) == sys_c);
                        trk_ok = Core_Utils.unique3ch(this.obs_code(id_obs_sys,:));
                    end
                    for t = 1 : size(trk_ok, 1)
                        if type(2) == 'o'
                            cur_res_id = cc.system == sys_c;
                            prn = cc.prn(cur_res_id);
                        else
                            cur_res_id = (this.system(id_obs)' == sys_c) & (this.obs_code(id_obs, 2) == trk_ok(t,2)) & (this.obs_code(id_obs, 3) == trk_ok(t,3));
                            [~, prn] = cc.getSysPrn(this.go_id(id_obs(cur_res_id)));                            
                        end
                        
                        f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s ResMap %s', f.Number, this.parent.getMarkerName4Ch, cc.getSysName(sys_c)); f.NumberTitle = 'off';
                        
                        fh_list = [fh_list; f]; %#ok<AGROW>
                        fig_name = sprintf('Res_Map_%s_%s_%s_%s_%s', type, trk_ok(t,:), this.parent.getMarkerName4Ch, cc.getSysName(sys_c), this.time.first.toString('yyyymmdd_HHMM'));
                        f.UserData = struct('fig_name', fig_name);
                                                                        
                        data_found = false;

                        res_tmp = res(:, cur_res_id);
                        ax1 = axes;
                        deg2rad = pi/180;
                        % For each satellite
                        
                        %res_tmp = Receiver_Commons.smoothMat(res_tmp, 'spline', 10);
                        zfilter = false;
                        time = this.getTime();
                        [year, doy] = time.getDOY();                        
                        if zfilter
                            % Zernike filter
                            l_max = 37;
                            m_max = l_max;
                            az_all = [];
                            el_all = [];
                            res_all = [];
                            id_subset = []; offset = 0;
                            for s = 1 : sum(cur_res_id)
                                id_ok = find(~isnan(zero2nan(res_tmp(:, s))));
                                [~, tmp] = intersect(id_ok, find(mod(doy,2) == 0 | mod(doy,2) == 1));
                                id_subset = [id_subset; tmp + offset];
                                offset = offset + numel(id_ok);
                                if any(id_ok)
                                    data_found = true;
                                    res_all = [res_all; res_tmp(id_ok, s)];
                                    
                                    go_id = cc.getIndex(sys_c, prn(s));
                                    az_all = [az_all; this.sat.az(id_ok, go_id) .* deg2rad];
                                    el_all = [el_all; (this.sat.el(id_ok, go_id)) .* deg2rad];
                                end
                            end
                            if data_found
                                el2radius = @(el) 1 - el ./ (pi/2);
                                [z_par, l, m, A] = Core_Utils.zAnalisysAll(l_max, m_max, az_all(id_subset), el2radius(el_all(id_subset)), res_all(id_subset) * scale, 1e-8);
                                %figure; scatter(m, -l, 150, z_par, 'filled');                                
                                [res_allf] = Core_Utils.zSinthesys(l, m, az_all, el2radius(el_all), z_par);                                
                                %[res_allf, z_par, l, m,  A] = Core_Utils.zFilter(l_max, m_max, az_all, el_all, res_all * scale);
                                [std(res_all*scale) std(res_all*scale - res_allf)]
                                res_all = res_all*scale - res_allf;
                                [~, id_sort] = sort(abs(res_all));
                                polarScatter(az_all(id_sort), pi/2 - el_all(id_sort), 20, res_all(id_sort), 'filled');                                
                                %hold off; Core_Utils.showZernike3StylePCV(l, m, z_par); drawnow; colormap((Cmap.get('PuOr', 2^11)));
                                %subplot(ax1);                                
                            end
                        else                            
                            for s = 1 : sum(cur_res_id)
                                id_ok = find(~isnan(zero2nan(res_tmp(:, s))));
                                if any(id_ok)
                                    data_found = true;
                                    go_id = cc.getIndex(sys_c, prn(s));
                                    [~, id_sort] = sort(abs(res_tmp(id_ok, s)));
                                    polarScatter(this.sat.az(id_ok(id_sort), go_id).*deg2rad, (90 - this.sat.el(id_ok(id_sort), go_id)).*deg2rad, 20, scale * (res_tmp(id_ok(id_sort), s)), 'filled');
                                    %polarScatter(this.sat.az(id_ok, go_id).*deg2rad, (90 - this.sat.el(id_ok, go_id)).*deg2rad, 80, scale * (res_tmp(id_ok, s)), 'filled');
                                    hold on;
                                end
                            end
                            caxis(perc(abs(serialize(scale * (res_tmp(id_ok, :)))), 0.997) .* [-1.1 1.1]);
                        end
                        
                        if ~data_found
                            close(f)
                            log = Core.getLogger;
                            log.addError(sprintf('No %s %s found in %s for constellation %s', name, trk_ok(t,:), this.parent.getMarkerName4Ch, cc.getSysName(sys_c)));
                        else
                            colorbar;
                            cax = caxis(ax1);
                            caxis(ax1, [-1 1] * max(abs(cax)));
                            %colormap(Cmap.get('RdBu', 2^11));
                            colormap((Cmap.get('PuOr', 2^11)));
                            %if min(abs(cax)) > 5
                            %    setColorMap('RdBu', caxis(), 0.90, [-5 5])
                            %end
                            cb = colorbar(ax1); cb.UserData = title(cb, iif(scale == 1e2, '[cm]', '[mm]')); ax1.Color = [0.9 0.9 0.9];
                            h = ylabel(ax1, 'Elevation [deg]'); h.FontWeight = 'bold';
                            grid(ax1, 'on');
                            h = xlabel(ax1, 'Azimuth [deg]'); h.FontWeight = 'bold';
                            if type(2) == 'o'
                                h = title(ax1, sprintf('%s %s %s\\fontsize{5} \n', cc.getSysName(sys_c), strrep(this.parent.marker_name, '_', '\_'), name), 'interpreter', 'tex'); h.FontWeight = 'bold';
                            else
                                h = title(ax1, sprintf('%s %s %s %s\\fontsize{5} \n', cc.getSysName(sys_c), strrep(this.parent.marker_name, '_', '\_'), trk_ok(t,:), name), 'interpreter', 'tex'); h.FontWeight = 'bold';
                            end
                            
                            Core_UI.beautifyFig(f, 'dark');
                            Core_UI.addBeautifyMenu(f);
                            f.Visible = 'on'; drawnow;
                        end
                    end
                end
            end
        end
        
        function fh_list = showResPerSat(this, sys_c_list, type, res)
            % Plot the residuals of phase per Satellite
            %
            % INPUT
            %   type    can be:
            %            'co'   -> Combined residuals (one set for each satellite) DEFAULT
            %            'pr'   -> Uncombined pseudo-ranges residuals
            %            'ph'   -> Uncombined carrier-phase residuals
            %   res     is the matrix of residuals satellite by satellite and can be passed from e.g. NET
            %
            % SYNTAX
            %   this.showResPerSat(sys_c_list, type, res)
            
            fh_list = [];
            cc = Core.getState.getConstellationCollector;
            if nargin < 2 || isempty(sys_c_list)
                sys_c_list = cc.getAvailableSys;
            end
            if nargin < 3 || isempty(type)
                type = 'co';
            end
            
            if nargin < 4 || isempty(res)                
                if type(2) == 'o'
                    res = this.sat.res;
                elseif type(2) == 'r'
                    res = this.sat.res_pr_by_pr;
                else
                    res = this.sat.res_ph_by_ph;
                end
            end
            if type(2) == 'o'
                name = 'Combined residuals sat. by sat.';
                scale = 1e3; % mm
                id_obs = 1 : size(res, 2);
            elseif type(2) == 'r'
                name = 'Uncombined pseudo-ranges residuals';
                id_obs = this.obs_code(:,1) == 'C';
                scale = 1e2; % cm
            else
                name = 'Uncombined carrier-phase residuals';
                id_obs = this.obs_code(:,1) == 'L';
                scale = 1e3; % mm
            end
            id_obs = find(id_obs);
            if isempty(res)
                log = Core.getLogger;
                log.addError(sprintf('No %s residuals found in %s', name, this.parent.getMarkerName4Ch));
            else
                if type(2) == 'o'
                    ss_ok = intersect(cc.sys_c, sys_c_list);
                else
                    ss_ok = intersect(unique(this.system(id_obs)), sys_c_list);                    
                end
                for sys_c = ss_ok
                    if type(2) == 'o'
                        trk_ok = 'COM'; % combined tracking -> only this is existing
                    else
                        id_obs_sys = id_obs(this.system(id_obs) == sys_c);
                        trk_ok = Core_Utils.unique3ch(this.obs_code(id_obs_sys,:));
                    end
                    for t = 1 : size(trk_ok, 1)
                        if type(2) == 'o'
                            cur_res_id = cc.system == sys_c;
                            prn = cc.prn(cur_res_id);
                        else
                            cur_res_id = (this.system(id_obs)' == sys_c) & (this.obs_code(id_obs, 2) == trk_ok(t,2)) & (this.obs_code(id_obs, 3) == trk_ok(t,3));
                            [~, prn] = cc.getSysPrn(this.go_id(id_obs(cur_res_id)));                            
                        end
                        
                        f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s Res %s', f.Number, this.parent.getMarkerName4Ch, cc.getSysName(sys_c)); f.NumberTitle = 'off';
                        
                        fh_list = [fh_list; f]; %#ok<AGROW>
                        fig_name = sprintf('Res_Per_Sat_%s_%s_%s_%s_%s', type, trk_ok(t,:), this.parent.getMarkerName4Ch, cc.getSysName(sys_c), this.time.first.toString('yyyymmdd_HHMM'));
                        f.UserData = struct('fig_name', fig_name);
                        
                        ep = repmat((1: this.time.length)',1, size(this.sat.outliers, 2));
                        
                        fun = @(err) min(256,max(1, round(256 / max(zero2nan(std(this.sat.res(:,:), 'omitnan')).*1e3) * err)));
                        ax2 = subplot(1, 24, 19:24);
                        ax1 = subplot(1, 24, 1:16);
                        
                        data_found = false;

                        res_tmp = res(:, cur_res_id);
                        for s = 1 : sum(cur_res_id)
                            id_ok = find(~isnan(zero2nan(res_tmp(:, s))));
                            if any(id_ok)
                                data_found = true;
                                [~, id_sort] = sort(abs(res_tmp(id_ok, s)));
                                scatter(ax1, id_ok(id_sort),  prn(s) * ones(size(id_ok)), 80, scale * (res_tmp(id_ok(id_sort), s)), 'filled');
                                hold(ax1, 'on');
                                err = std(zero2nan(res_tmp(:,s)), 'omitnan') * scale;
                                if  verLessThan('matlab', '9.4') 
                                    plot(ax2, mean(zero2nan(res_tmp(:,s)), 'omitnan') + [-err err], prn(s) * [1 1], '.-', 'MarkerSize', 15, 'LineWidth', 3, 'Color', [0.6 0.6 0.6]);
                                    plot(ax2, mean(zero2nan(res_tmp(:,s)), 'omitnan'), prn(s), '.', 'MarkerSize', 30, 'Color', [0.6 0.6 0.6]);
                                else
                                    errorbar(ax2, mean(zero2nan(res_tmp(:,s)), 'omitnan') .* scale, prn(s), err, '.', 'horizontal', 'MarkerSize', 30, 'LineWidth', 3, 'Color', [0.6 0.6 0.6]);
                                end
                                hold(ax2, 'on');
                            end
                        end
                        
                        if ~data_found
                            close(f)
                            log = Core.getLogger;
                            log.addError(sprintf('No %s %s found in %s for constellation %s', name, trk_ok(t,:), this.parent.getMarkerName4Ch, cc.getSysName(sys_c)));
                        else
                            cax = caxis(ax1); caxis(ax1, [-1 1] * max(abs(cax)));
                            colormap(Cmap.get('RdBu', 2^11));
                            if min(abs(cax)) > 5
                                setColorMap('RdBu', caxis(), 0.90, [-5 5])
                            end
                            cb = colorbar(ax1); cb.UserData = title(cb, iif(scale == 1e2, '[cm]', '[mm]')); ax1.Color = [0.9 0.9 0.9];
                            prn_ss = unique(cc.prn(cc.system == sys_c));
                            xlim(ax1, [1 size(this.sat.res,1)]);
                            ylim(ax1, [min(prn_ss) - 1 max(prn_ss) + 1]);
                            h = ylabel(ax1, 'PRN'); h.FontWeight = 'bold';
                            ax1.YTick = prn_ss;
                            grid(ax1, 'on');
                            h = xlabel(ax1, 'epoch'); h.FontWeight = 'bold';
                            if type(2) == 'o'
                                h = title(ax1, sprintf('%s %s %s\\fontsize{5} \n', cc.getSysName(sys_c), strrep(this.parent.marker_name, '_', '\_'), name), 'interpreter', 'tex'); h.FontWeight = 'bold';
                            else
                                h = title(ax1, sprintf('%s %s %s %s\\fontsize{5} \n', cc.getSysName(sys_c), strrep(this.parent.marker_name, '_', '\_'), trk_ok(t,:), name), 'interpreter', 'tex'); h.FontWeight = 'bold';
                            end
                            
                            ylim(ax2, [min(prn_ss) - 1 max(prn_ss) + 1]);
                            xlim(ax2, [-1 1] * (max(max(abs(mean(zero2nan(res(:,:)), 'omitnan'))) * scale, ...
                                max(std(zero2nan(res(:,:)), 'omitnan')) * scale) + 1));
                            ax2.YTick = prn_ss; ax2.Color = [1 1 1];
                            grid(ax2, 'on');
                            xlabel(ax2, sprintf('mean %s', iif(scale == 1e2, 'cm', 'mm')));
                            h = title(ax2, sprintf('mean\\fontsize{5} \n'), 'interpreter', 'tex'); h.FontWeight = 'bold';
                            linkaxes([ax1, ax2], 'y');
                            
                            Core_UI.beautifyFig(f, 'dark');
                            Core_UI.addBeautifyMenu(f);
                            f.Visible = 'on'; drawnow;
                        end
                    end
                end
            end
        end
        
        function fh_list = showAniZtdSlant(this, time_start, time_stop, show_map, write_video)
            sztd = this.getSlantZTD(this.parent.slant_filter_win);
            
            fh_list = [];
            if isempty(this.ztd) || ~any(sztd(:))
                this.log.addWarning('ZTD and slants have not been computed');
            else
                f = figure; f.Name = sprintf('%03d: AniZtd', f.Number); f.NumberTitle = 'off';
                fh_list = [fh_list; f]; %#ok<AGROW>
                fig_name = sprintf('ZTD_Slant_ANI_%s_%s', rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                f.UserData = struct('fig_name', fig_name);
                
                if nargin >= 3
                    if isa(time_start, 'GPS_Time')
                        time_start = find(this.time.getMatlabTime >= time_start.first.getMatlabTime(), 1, 'first');
                        time_stop = find(this.time.getMatlabTime <= time_stop.last.getMatlabTime(), 1, 'last');
                    end
                    time_start = max(1, time_start);
                    time_stop = min(size(sztd,1), time_stop);
                else
                    time_start = 1;
                    time_stop = size(sztd,1);
                end
                
                id_sync = serialize(this(1).getIdSync);
                if isempty(id_sync)
                    id_sync(:, 1) = (1 : this.time.length())';
                end
                id_ok = id_sync(id_sync(:) > time_start & id_sync(:) < time_stop);
                t = this.time.getEpoch(id_ok).getMatlabTime;
                sztd = sztd(id_ok, :);
                
                if nargin < 4
                    show_map = true;
                end
                if nargin < 5
                    write_video = false;
                else
                    if write_video
                        vidObj = VideoWriter('./out.avi');
                        vidObj.FrameRate = 30;
                        vidObj.Quality = 100;
                        open(vidObj);
                    end
                end
                yl = (median(median(sztd, 'omitnan'), 'omitnan') + ([-6 6]) .* median(std(sztd, 'omitnan'), 'omitnan')) * 1e2;
                
                subplot(3,1,3);
                plot(t, sztd * 1e2,'.'); hold on;
                plot(t, this.ztd(id_ok) * 1e2,'k', 'LineWidth', 4);
                ylim(yl);
                hl = line('XData', t(1) * [1 1],'YData', yl, 'LineWidth', 2);
                xlim([t(1) t(end)]);
                setTimeTicks(4);
                h = ylabel('ZTD [cm]'); h.FontWeight = 'bold';
                grid on;
                
                % polar plot "true" Limits
                e_grid = [-1 : 0.2 : 1];
                n_grid = [-1 : 0.2 : 1];
                [ep, np] = meshgrid(e_grid, n_grid);
                fun = @(dist) exp(-((dist*1e5)/3e4).^2);
                
                ax_sky = subplot(3,1,1:2); i = time_start;
                az = (mod(this.sat.az(id_ok(i),:) + 180, 360) -180) ./ 180 * pi; az(isnan(az) | isnan(sztd(i,:))) = 1e10;
                el = (90 - this.sat.el(id_ok(i),:)) ./ 180 * pi; el(isnan(el) | isnan(sztd(i,:))) = 1e10;
                
                if show_map
                    td = nan(size(ep));
                    hm = imagesc(e_grid, n_grid, reshape(td(:), numel(n_grid), numel(e_grid))); hold on;
                    hm.AlphaData = 0.5;
                    ax_sky.YDir = 'normal';
                end
                hs = polarScatter(az, el, 250, sztd(i,:) * 1e2, 'filled');
                xlim([-1 1]); ylim([-1 1]);
                caxis(yl); colormap(jet(1024)); colorbar;
                
                subplot(3,1,3);
                for i = 2 : 2 : numel(id_ok)
                    % Move scattered points
                    az = (mod(this.sat.az(id_sync(i, 1),:) + 180, 360) -180) ./ 180 * pi; az(isnan(az) | isnan(sztd(i,:))) = 1e10;
                    el = (90 - this.sat.el(id_sync(i, 1),:)) ./ 180 * pi; el(isnan(el) | isnan(sztd(i,:))) = 1e10;
                    decl_n = el/(pi/2);
                    x = sin(az) .* decl_n;
                    y = cos(az) .* decl_n;
                    
                    id_ok = not(isnan(zero2nan(sztd(i,:))));
                    if show_map
                        if any(id_ok(:))
                            td = funInterp2(ep(:), np(:), x(1, id_ok)', y(1, id_ok)', sztd(i, id_ok)' * 1e2, fun);
                            hm.CData = reshape(td(:), numel(n_grid), numel(e_grid));
                        end
                    end
                    
                    hs.XData = x;
                    hs.YData = y;
                    hs.CData = sztd(i,:) * 1e2;
                    
                    % Move time line
                    hl.XData = t(i) * [1 1];
                    drawnow;
                    
                    if write_video
                        currFrame = export_fig(fig_h, '-nocrop', '-a1');
                        writeVideo(vidObj,currFrame);
                    end
                end
                if write_video
                    close(vidObj);
                end
            end
        end
        
        function fh_list = showAniZwdSlant(this, time_start, time_stop, show_map)
            fh_list = [];
            if isempty(this.zwd)
                this.log.addWarning('ZWD have not been computed');
            else
                f = figure; f.Name = sprintf('%03d: AniZwd', f.Number); f.NumberTitle = 'off';
                
                fh_list = [fh_list; f]; %#ok<AGROW>
                fig_name = sprintf('ZWD_Slant_ANI_%s_%s', rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                f.UserData = struct('fig_name', fig_name);
                
                szwd = this.getSlantZWD(this.parent.slant_filter_win);
                
                if nargin >= 3
                    if isa(time_start, 'GPS_Time')
                        time_start = find(this.time.getMatlabTime >= time_start.first.getMatlabTime(), 1, 'first');
                        time_stop = find(this.time.getMatlabTime <= time_stop.last.getMatlabTime(), 1, 'last');
                    end
                    time_start = max(1, time_start);
                    time_stop = min(size(szwd,1), time_stop);
                else
                    time_start = 1;
                    time_stop = size(szwd,1);
                end
                
                if isempty(this.id_sync(:))
                    this.id_sync = (1 : this.time.length())';
                end
                id_ok = this.id_sync(this.id_sync > time_start & this.id_sync < time_stop, 1);
                
                t = this.time.getEpoch(id_ok).getMatlabTime;
                szwd = szwd(id_ok, :);
                
                if nargin < 4
                    show_map = true;
                end
                yl = (median(median(szwd, 'omitnan'), 'omitnan') + ([-6 6]) .* median(std(szwd, 'omitnan'), 'omitnan')) * 1e2;
                
                subplot(3,1,3);
                plot(t, szwd * 1e2,'.'); hold on;
                plot(t, this.zwd(id_ok) * 1e2,'k', 'LineWidth', 4);
                ylim(yl);
                hl = line('XData', t(1) * [1 1],'YData', yl, 'LineWidth', 2);
                xlim([t(1) t(end)]);
                setTimeTicks(4);
                h = ylabel('ZWD [cm]'); h.FontWeight = 'bold';
                grid on;
                
                % polar plot "true" Limits
                e_grid = [-1 : 0.1 : 1];
                n_grid = [-1 : 0.1 : 1];
                [ep, np] = meshgrid(e_grid, n_grid);
                fun = @(dist) exp(-((dist*1e5)/3e4).^2);
                
                ax_sky = subplot(3,1,1:2); i = time_start;
                az = (mod(this.sat.az(id_ok(i),:) + 180, 360) -180) ./ 180 * pi; az(isnan(az) | isnan(szwd(i,:))) = 1e10;
                el = (90 - this.sat.el(id_ok(i),:)) ./ 180 * pi; el(isnan(el) | isnan(szwd(i,:))) = 1e10;
                
                if show_map
                    td = nan(size(ep));
                    hm = imagesc(e_grid, n_grid, reshape(td(:), numel(n_grid), numel(e_grid))); hold on;
                    hm.AlphaData = 0.5;
                    ax_sky.YDir = 'normal';
                end
                hs = polarScatter(az, el, 250, szwd(i,:) * 1e2, 'filled');
                caxis(yl); colormap(jet(1024)); colorbar;
                
                subplot(3,1,3);
                for i = 2 : numel(id_ok)
                    % Move scattered points
                    az = (mod(this.sat.az(this.id_sync(i, 1),:) + 180, 360) -180) ./ 180 * pi; az(isnan(az) | isnan(szwd(i,:))) = 1e10;
                    el = (90 - this.sat.el(this.id_sync(i, 1),:)) ./ 180 * pi; el(isnan(el) | isnan(szwd(i,:))) = 1e10;
                    decl_n = el/(pi/2);
                    x = sin(az) .* decl_n;
                    y = cos(az) .* decl_n;
                    
                    id_ok = not(isnan(zero2nan(szwd(i,:))));
                    if show_map
                        if any(id_ok(:))
                            td = funInterp2(ep(:), np(:), x(1, id_ok)', y(1, id_ok)', szwd(i, id_ok)' * 1e2, fun);
                            hm.CData = reshape(td(:), numel(n_grid), numel(e_grid));
                        end
                    end
                    
                    hs.XData = x;
                    hs.YData = y;
                    hs.CData = szwd(i,:) * 1e2;
                    
                    % Move time line
                    hl.XData = t(i) * [1 1];
                    drawnow;
                end
            end
        end
        
        function fh_list = showZtdSlant(this, time_start, time_stop)
            fh_list = [];
            rec = this;
            if isempty(rec)
                this(1).log.addWarning('ZTD and/or slants have not been computed');
            else
                cc = this.getCC;
                f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s Slant %s', f.Number, this.parent.getMarkerName4Ch, cc.sys_c); f.NumberTitle = 'off';
                fh_list = [fh_list; f];
                fig_name = sprintf('ZTD_Slant_%s_%s', rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                f.UserData = struct('fig_name', fig_name);
                
                t = rec(:).getTime.getMatlabTime;
                
                sztd = rec(:).getSlantZTD(rec(1).parent.slant_filter_win);
                if isempty(sztd)
                    close(f);
                else
                    if nargin >= 3
                        if isa(time_start, 'GPS_Time')
                            time_start = find(t >= time_start.first.getMatlabTime(), 1, 'first');
                            time_stop = find(t <= time_stop.last.getMatlabTime(), 1, 'last');
                        end
                        time_start = max(1, time_start);
                        time_stop = min(size(sztd,1), time_stop);
                    else
                        time_start = 1;
                        time_stop = size(sztd,1);
                    end
                    
                    if nargin < 4
                        win_size = (t(time_stop) - t(time_start)) * 86400;
                    end
                    
                    %yl = (median(median(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan') + ([-6 6]) .* median(std(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan'));
                    
                    plot(t, sztd,'.-'); hold on;
                    plot(t, zero2nan(rec(:).getZtd),'k', 'LineWidth', 4);
                    %ylim(yl);
                    %xlim(t(time_start) + [0 win_size-1] ./ 86400);
                    setTimeTicks(4);
                    h = ylabel('ZTD [m]'); h.FontWeight = 'bold';
                    grid on;
                    h = title(sprintf('Receiver %s ZTD', rec(1).parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                    drawnow;
                    Core_UI.beautifyFig(f, 'dark');
                    Core_UI.addBeautifyMenu(f);
                    f.Visible = 'on'; drawnow;
                end
            end
            
            
        end
        
        
        function fh_list = showTropoPar(sta_list, par_name, new_fig)
            % one function to rule them all
            
            fh_list = [];
            [tropo, t] = sta_list.getTropoPar(par_name);
            if ~iscell(tropo)
                tropo = {tropo};
                t = {t};
            end
            
            unit_convert = iif(ismember(lower(par_name), {'ztd', 'zhd', 'zwd', 'pwv', 'gn', 'ge'}), 1e2, 1);            
            
            rec_ok = false(numel(sta_list), 1);
            for r = 1 : size(sta_list, 2)
                rec_ok(r) = ~isempty(tropo{r});
            end
            
            sta_list = sta_list(rec_ok);
            tropo = tropo(rec_ok);
            t = t(rec_ok);
            
            if numel(sta_list) == 0
                log = Logger.getInstance();
                log.addError('No valid troposphere is present in the receiver list');
            else
                cc = Core.getState.getConstellationCollector;
                if nargin < 3
                    new_fig = true;
                end
                
                if isempty(tropo)
                    sta_list(1).log.addWarning([par_name ' and slants have not been computed']);
                else
                    if new_fig
                        f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s %s', f.Number, par_name, cc.sys_c); f.NumberTitle = 'off';
                        old_legend = {};
                    else
                        f = gcf;
                        l = legend;
                        old_legend = get(l,'String');
                    end
                    
                    fh_list = [fh_list; f];
                    if numel(sta_list) == 1 
                        % If I have only one receiver use as name the name of the receiver
                        fig_name = sprintf('%s_%s_%s', upper(par_name), sta_list.parent.getMarkerName4Ch, sta_list.time.first.toString('yyyymmdd_HHMM'));
                    else
                        % If I have more than one receiver use as name the name of the project
                        fig_name = sprintf('%s_%s', upper(par_name), strrep(Core.getState.getPrjName,' ', '_'));
                    end
                    f.UserData = struct('fig_name', fig_name);
                    
                    for r = 1 : numel(sta_list)
                        rec = sta_list(r);
                        if new_fig
                            Core_Utils.plotSep(t{r}.getMatlabTime(), zero2nan(tropo{r}') .* unit_convert, '.-', 'LineWidth', 2, 'Color', Core_UI.getColor(r, size(sta_list, 2))); hold on;
                        else
                            Core_Utils.plotSep(t{r}.getMatlabTime(), zero2nan(tropo{r}') .* unit_convert, '.-', 'LineWidth', 2); hold on;
                        end
                        outm{r} = rec(1).parent.getMarkerName();
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
                    setTimeTicks(4);
                    h = ylabel([par_name ' [cm]']); h.FontWeight = 'bold';
                    grid on;
                    
                     % Generate Name
                    switch lower(par_name)
                        case 'ztd'
                            ttl = 'Estimated Zenith Total Delay (ZTD)';
                        case 'zwd'
                            ttl = 'Estimated Zenith Wet Delay (ZWD)';
                        case 'gn'
                            ttl = 'Estimated Wet Delay North Gradient';
                        case 'ge'
                            ttl = 'Estimated Wet Delay East Gradient';
                        case 'pwv'
                            ttl = 'Estimated Precipitable Water Vapour (PWV)';
                        case 'zhd'
                            ttl = 'Estimated Zenith Hydrostatic Delay (ZHD)';
                        case 'nsat'
                            ttl = 'Number of used satellites';
                    end
                    
                    h = title(ttl); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                    Core_UI.beautifyFig(f);
                    Core_UI.addBeautifyMenu(f);
                    f.Visible = 'on'; drawnow;
                end
            end
        end
        
        function fh_list = showZhd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showTropoPar('ZHD', new_fig);
        end
        
        function fh_list = showZwd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showTropoPar('ZWD', new_fig);
        end
        
        function fh_list = showZtd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showTropoPar('ZTD', new_fig);
        end

        function fh_list = showZtdSlantRes_p(this, time_start, time_stop)
            
            fh_list = [];
            if isempty(this.sat.res)
                this.log.addWarning('ZTD and slants have not been computed');
            else                
                id_sync = this.getIdSync();
                
                t = this.time.getEpoch(id_sync).getMatlabTime;
                
                sztd = this.getSlantZTD(this.parent.slant_filter_win);
                if isempty(sztd)
                    this.log.addWarning('ZTD and slants have not been computed');
                else
                    sztd = bsxfun(@minus, sztd, this.ztd(id_sync)) .* 1e2;
                    if nargin >= 3
                        if isa(time_start, 'GPS_Time')
                            time_start = find(t >= time_start.first.getMatlabTime(), 1, 'first');
                            time_stop = find(t <= time_stop.last.getMatlabTime(), 1, 'last');
                        end
                        time_start = max(1, time_start);
                        time_stop = min(size(sztd,1), time_stop);
                    else
                        time_start = 1;
                        time_stop = size(sztd,1);
                    end
                    
                    %yl = (median(median(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan') + ([-6 6]) .* median(std(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan'));
                    
                    az = (mod(this.sat.az(id_sync,:) + 180, 360) -180) ./ 180 * pi; az(isnan(az) | isnan(sztd)) = 1e10;
                    el = (90 - this.sat.el(id_sync,:)) ./ 180 * pi; el(isnan(el) | isnan(sztd)) = 1e10;
                    
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s Slant Residuals', f.Number, this.parent.getMarkerName4Ch); f.NumberTitle = 'off';
                    
                    fh_list = f;
                    if numel(this) == 1 
                        % If I have only one receiver use as name the name of the receiver
                        fig_name = sprintf('%s_%s_%s', 'ZTDS', this.parent.getMarkerName4Ch, this.time.first.toString('yyyymmdd_HHMM'));
                    else
                        % If I have more than one receiver use as name the name of the project
                        fig_name = sprintf('%s_%s', 'ZTDS', strrep(Core.getState.getPrjName,' ', '_'));
                    end
                    f.UserData = struct('fig_name', fig_name);
                    
                    polarScatter(az(:), el(:), 25, abs(sztd(:)), 'filled'); hold on;
                    caxis(minMax(abs(sztd))); colormap(flipud(hot)); f.Color = [.95 .95 .95];
                    cb = colorbar(); cbt = title(cb, '[cm]'); cbt.Parent.UserData = cbt;
                    h = title(sprintf('Receiver %s ZTD - Slant difference', this.parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                    Core_UI.beautifyFig(f, 'dark');
                    Core_UI.addBeautifyMenu(f);
                    f.Visible = 'on'; drawnow;
                end
            end
        end

        function fh_list = showGn(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showTropoPar('GN', new_fig);
        end
         
        function fh_list = showGe(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showTropoPar('GE', new_fig);
        end
              
        function fh_list = showPwv(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showTropoPar('PWV', new_fig);
        end
        
        function fh_list = showNSat(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showTropoPar('nsat', new_fig);
        end
        
        function fh_list = showNSatSS(this, flag_smooth)
            % Show number of satellites in view per constellation
            if nargin == 1
                flag_smooth = true;
            end
            fh_list = [];
            cc = this.getCC;
            if ~this.isEmpty()
                [n_sat, n_sat_ss] = this.getNSat();
                f = figure('Visible', 'off'); f.Name = sprintf('%03d: nsat SS %s', f.Number, this.parent.getMarkerName4Ch); f.NumberTitle = 'off';
                
                fh_list = [fh_list; f];
                fig_name = sprintf('NSatSS_%s_%s', this.parent.getMarkerName4Ch, this.time.first.toString('yyyymmdd_HHMM'));
                f.UserData = struct('fig_name', fig_name);
                       
                %                 if isfield(n_sat_ss, 'S')
                %                     plot(this.getTime.getMatlabTime(), zero2nan(n_sat_ss.S), '.', 'MarkerSize', 40, 'LineWidth', 2); hold on;
                %                 end
                %                 if isfield(n_sat_ss, 'I')
                %                     plot(this.getTime.getMatlabTime(), zero2nan(n_sat_ss.I), '.', 'MarkerSize', 35, 'LineWidth', 2); hold on;
                %                 end
                %                 if isfield(n_sat_ss, 'C')
                %                     plot(this.getTime.getMatlabTime(), zero2nan(n_sat_ss.C), '.', 'MarkerSize', 30, 'LineWidth', 2); hold on;
                %                 end
                %                 if isfield(n_sat_ss, 'J')
                %                     plot(this.getTime.getMatlabTime(), zero2nan(n_sat_ss.J), '.', 'MarkerSize', 25, 'LineWidth', 2); hold on;
                %                 end
                %                 if isfield(n_sat_ss, 'E')
                %                     plot(this.getTime.getMatlabTime(), zero2nan(n_sat_ss.E), '.', 'MarkerSize', 20, 'LineWidth', 2); hold on;
                %                 end
                %                 if isfield(n_sat_ss, 'R')
                %                     plot(this.getTime.getMatlabTime(), zero2nan(n_sat_ss.R), '.', 'MarkerSize', 15, 'LineWidth', 2); hold on;
                %                 end
                %                 if isfield(n_sat_ss, 'G')
                %                     plot(this.getTime.getMatlabTime(), zero2nan(n_sat_ss.G), '.', 'MarkerSize', 10, 'LineWidth', 2); hold on;
                %                 end
                
                % If I'm plotting more than one day smooth the number of satellites
                if flag_smooth % && (this.getTime.last.getMatlabTime - this.getTime.first.getMatlabTime) > 1
                    for sys_c = cc.sys_c
                        Core_Utils.plotSep(this.getTime.getMatlabTime(), splinerMat(this.getTime.getRefTime, zero2nan(n_sat_ss.(sys_c)), 3600), '.-', 'MarkerSize', 10); hold on;
                    end
                    if numel(cc.sys_c) > 1
                        Core_Utils.plotSep(this.getTime.getMatlabTime(), splinerMat(this.getTime.getRefTime, zero2nan(n_sat), 3600), '.-k', 'MarkerSize', 10); hold on;
                    end
                else % If I'm plotting less than 24 hours of satellites number
                    Core_Utils.plotSep(this.getTime.getMatlabTime(), zero2nan(struct2array(n_sat_ss)), '.-', 'MarkerSize', 10); hold on;
                    if numel(cc.sys_c) > 1
                        Core_Utils.plotSep(this.getTime.getMatlabTime(), zero2nan(n_sat), '.-k', 'MarkerSize', 10);
                    end
                end
                setTimeTicks(4); h = ylabel('East [cm]'); h.FontWeight = 'bold';
                
                sys_list = {};
                for i = 1 : numel(cc.sys_c)
                    sys_list = [sys_list, {cc.sys_c(i)}];
                end
                if numel(cc.sys_c) > 1                    
                    legend([sys_list, {'All'}]);
                else
                    legend(sys_list);
                end
                %legend(sys_list);
                ylim([0 (max(serialize(n_sat)) + 1)]);
                grid minor;
                h = title(sprintf('N sat per constellation - %s', this.parent.getMarkerName4Ch),'interpreter', 'none'); h.FontWeight = 'bold';
                
                Core_UI.beautifyFig(f);
                Core_UI.addBeautifyMenu(f);
                f.Visible = 'on'; drawnow;
            end
        end
       
        
        function fh_list = showMedianTropoPar(this, par_name, new_fig);
            % one function to rule them all
            fh_list = [];
            rec_ok = false(size(this,2), 1);
            for r = 1 : size(this, 2)
                rec_ok(r) = any(~isnan(this(:,r).getZtd));
            end
            rec_list = this(:, rec_ok);
            
            if nargin < 3
                new_fig = true;
            end
            
            switch lower(par_name)
                case 'ztd'
                    [tropo] = rec_list.getZtd();
                case 'zwd'
                    [tropo] = rec_list.getZwd();
                case 'pwv'
                    [tropo] = rec_list.getPwv();
                case 'zhd'
                    [tropo] = rec_list.getAprZhd();
            end
            
            if ~iscell(tropo)
                tropo = {tropo};
            end
            if isempty(tropo)
                rec_list(1).log.addWarning([par_name ' and slants have not been computed']);
            else
                if new_fig
                    cc = this.getCC;
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: Median %s %s', f.Number, par_name, cc.sys_c); f.NumberTitle = 'off';
                    old_legend = {};
                else
                    f = gcf;
                    l = legend;
                    old_legend = get(l,'String');
                end
                fh_list = [fh_list; f]; %#ok<AGROW>
                % If I have more than one receiver use as name the name of the project
                fig_name = sprintf('%s_Median_%s', upper(par_name), strrep(Core.getState.getPrjName,' ', '_'));
                f.UserData = struct('fig_name', fig_name);
                
                for r = 1 : size(rec_list, 2)
                    rec = rec_list(~rec_list(:,r).isEmpty, r);
                    if ~isempty(rec)
                        switch lower(par_name)
                            case 'ztd'
                                [tropo] = rec.getZtd();
                            case 'zwd'
                                [tropo] = rec.getZwd();
                            case 'pwv'
                                [tropo] = rec.getPwv();
                            case 'zhd'
                                [tropo] = rec.getAprZhd();
                        end
                        [~, ~, ~, h_o] = rec(1).getPosGeodetic();
                        if new_fig
                            plot(h_o, median(tropo,'omitnan'), '.', 'MarkerSize', 25, 'LineWidth', 4, 'Color', Core_UI.getColor(r, size(rec_list, 2))); hold on;
                        else
                            plot(h_o, median(tropo,'omitnan'), '.', 'MarkerSize', 25, 'LineWidth', 4); hold on;
                        end
                        outm{r} = rec(1).parent.getMarkerName();
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
                Core_UI.beautifyFig(f);
                Core_UI.addBeautifyMenu(f);
                f.Visible = 'on'; drawnow;
            end
        end
        
        function fh_list = showMedianZhd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showMedianTropoPar('ZHD', new_fig);
        end
        
        function fh_list = showMedianZwd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showMedianTropoPar('ZWD', new_fig);
        end
        
        function fh_list = showMedianZtd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showMedianTropoPar('ZTD', new_fig);
        end
        
        function fh_list = showMedianPwv(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showMedianTropoPar('PWV', new_fig);
        end
                
        function plotResidual(this, flag_smooth)
            % Legacy plot residuals
            %
            % INPUT
            %   flag_smooth ift can be logical or a value (for spline_base_size) 
            %
            % SYNTAX
            %   this.plotResidual(flag_smooth)
            figure
            id_ok = any(this.sat.res);
            if nargin == 2 && flag_smooth
                if islogical(flag_smooth)
                    n_spline = 1800/30;
                else
                    n_spline = flag_smooth;
                end
                plot(this.smoothMat(zero2nan(this.sat.res(:,id_ok)), 'spline', n_spline),'.-');
            else
                plot(zero2nan(this.sat.res(:,id_ok)),'.');
            end
            cc = Core.getConstellationCollector();
            ant_ids = cc.getAntennaId();
            legend(ant_ids(id_ok));
            
        end
    end    
    %% METHODS UTILITIES FUNCTIONS
    % ==================================================================================================================================================
    
    methods (Access = public)
        function [map, map_fill, n_data_map, az_g, el_g] = getResMap(rec, step, size_conv, sys_c)
            % Export a map containing the residuals
            % This export works better when a sequence of session is passed to it
            % Note: it works one riceiver at a time
            %
            % OUTPUT
            %   snr_map         cartesian map of the mean observed SNR
            %   snr_map_fill    cartesian map of the mean observed SNR filled in polar view
            %   snr_mask        mask of all the position with snr < threshold
            %   n_data_map      number of data used for the mean (obs falling in the cell)
            %   out_map         map of the number of outliers flagged by goGPS per cell
            %
            % SYNTAX
            %   [map, map_fill, n_data_map, az_g, el_g] = this.getResMap(<step = 0.5>, <size_conv = 21>, <snr_thr = 45>, <sys_c_list>);
            
            use_work = false;
            cc = rec.getCC;
            
            if nargin < 2 || isempty(step)
                step = 0.5;
            end
            if nargin < 3 || isempty(size_conv)
                size_conv = 21;
            end
            
            log = Core.getLogger();
            if nargin < 4 || isempty(sys_c)
                sys_c = rec.getAvailableSys;
                sys_c = sys_c(1);
            end
            
            % Create an el/az grid
            [phi_g, az_g] = getGrid(step);
            az_g = sort(mod(az_g, 360));
            el_g = phi_g(phi_g > 0);
            
            % [num_occur] = hist(id_map, unique(id_map));
            map = zeros(numel(el_g), numel(az_g));
            n_data_map = zeros(numel(el_g), numel(az_g));
            
            rec.w_bar.createNewBar('Computing map');
            rec.w_bar.setBarLen(numel(rec));
            
            data = (rec.getResidual());
            if sys_c == 'A'
                [sys, prn] = cc.getSysPrn(1:size(data,2));
                id_keep = 1 : numel(sys);
            else
                [sys, prn] = cc.getSysPrn(1:size(data,2));
                id_keep = false(size(sys));
                for i = 1 : numel(sys)
                    id_keep(i) = instr(sys_c, sys(i));
                end
                data = data(:, id_keep);
            end
            az = rec.sat.az(:, id_keep);
            el = rec.sat.el(:, id_keep);
            
            if ~isempty(data)
                % Extract non NaN serialized data
                data = zero2nan(data);
                id_ok = (~isnan(data));
                
                % Eliminate empty epochs
                az = az(:, sum(id_ok) > 1);
                el = el(:, sum(id_ok) > 1);
                data = data(:, sum(id_ok) > 1);
                
                id_oo = (~isnan(data(:)));
                
                id_az = max(1, ceil(mod(az(id_oo), 360) / step)); % Get the index of the grid
                id_el = (numel(el_g)) - floor(max(0, min(90 - step/2, el(id_oo)) / step)); % Get the index of the grid
                
                id_map = (id_az - 1) * numel(el_g) + id_el;
                data = serialize(data(id_oo));
                for i = 1 : numel(id_map)
                    map(id_map(i)) = map(id_map(i)) + data(i);
                    n_data_map(id_map(i)) = n_data_map(id_map(i)) + 1;
                end
            end
            rec.w_bar.go();
            
            rec.w_bar.close();
            if isempty(data)
                log.addWarning(sprintf('No data found for %s', rec.parent.getMarkerName4Ch), 100);
                map_fill = nan(size(map));
            else
                map(n_data_map(:) > 0) = map(n_data_map(:) > 0) ./ n_data_map(n_data_map(:) > 0);
                map(n_data_map(:) == 0) = nan;
                
                % Convert map in polar coordinates to allow a better interpolation
                % In principle it is possible to compute directly this polar map
                % This can be implemented some lines above this point
                step_p = step;
                polar_map = ones(ceil(360 / step_p)) * 0;
                polar_map(2 : end - 1, 2 : end - 1) = nan;
                
                dc_g = (90 - el_g) / 90; % declination
                [az_mesh, dc_mesh] = meshgrid(az_g ./ 180 * pi, dc_g);
                x = round((sin(az_mesh) .* dc_mesh) / step_p * 180) + size(polar_map , 1) / 2;
                y = round((cos(az_mesh) .* dc_mesh) / step_p * 180) + size(polar_map , 1) / 2;
                id_p = size(polar_map,1) .* (x(:)) + (y(:)+1); % id of the polar projection corresponding to the cartesian one
                polar_map(id_p) = map(:);
                %polar_map_fill = max(min(snr_map(:)), min(max(snr_map(:)), simpleFill2D(polar_map, isnan(polar_map),  @(dist) exp(-((dist)).^2))));
                polar_map_fill = max(min(map(:)), min(max(map(:)), inpaint_nans(polar_map)));
                if numel(size_conv) > 0 && ~((numel(size_conv) == 1) && size_conv == 0)
                    polar_map_fill = circConv2(polar_map_fill, size_conv);
                end
                map_fill = nan(size(map));
                map_fill(:) = polar_map_fill(id_p);                
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
            %   [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(rec, p_rate);
            %
            % EXAMPLE:
            %   [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(rec, 30);
            
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
                    % p_rate = lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
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
        
        function [p_time, id_sync] = getSyncTimeTR(sta_list, obs_type, p_rate)
            % Get the common (shortest) time among all the used receivers and the target(s)
            % For each target (obs_type == 0) produces a different cell array with the sync of the other receivers
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
                    len(r) = sta_list(r).length;
                end
                obs_type = ones(1, numel(sta_list));
                obs_type(find(len == max(len), 1, 'first')) = 0;
            end
            
            % Do the target(s) as last
            [~, id] = sort(obs_type, 'descend');
            
            % prepare reference time
            % processing time will start with the receiver with the last first epoch
            %          and it will stop  with the receiver with the first last epoch
            
            first_id_ok = find(~sta_list.isEmpty_mr, 1, 'first');
            p_time_zero = round(sta_list(first_id_ok).time.first.getMatlabTime() * 24)/24; % get the reference time
            p_time_start = sta_list(first_id_ok).time.first.getRefTime(p_time_zero);
            p_time_stop = sta_list(first_id_ok).time.last.getRefTime(p_time_zero);
            p_rate = lcm(round(p_rate * 1e6), round(sta_list(first_id_ok).time.getRate * 1e6)) * 1e-6;
            
            p_time = GPS_Time(); % empty initialization
            
            i = 0;
            for r = id
                ref_t{r} = sta_list(r).time.getRefTime(p_time_zero);
                if obs_type(r) > 0 % if it's not a target
                    if ~sta_list(r).isEmpty
                        p_time_start = max(p_time_start,  round(sta_list(r).time.first.getRefTime(p_time_zero) * sta_list(r).time.getRate) / sta_list(r).time.getRate);
                        p_time_stop = min(p_time_stop,  round(sta_list(r).time.last.getRefTime(p_time_zero) * sta_list(r).time.getRate) / sta_list(r).time.getRate);
                        p_rate = lcm(round(p_rate * 1e6), round(sta_list(r).time.getRate * 1e6)) * 1e-6;
                    end
                else
                    % It's a target
                    
                    % recompute the parameters for the ref_time estimation
                    % not that in principle I can have up to num_trg_rec ref_time
                    % in case of multiple targets the reference times should be independent
                    % so here I keep the temporary rt0 rt1 r_rate var
                    % instead of ref_time_start, ref_time_stop, ref_rate
                    pt0 = max(p_time_start, round(sta_list(r).time.first.getRefTime(p_time_zero) * sta_list(r).time.getRate) / sta_list(r).time.getRate);
                    pt1 = min(p_time_stop, round(sta_list(r).time.last.getRefTime(p_time_zero) * sta_list(r).time.getRate) / sta_list(r).time.getRate);
                    pr = lcm(round(p_rate * 1e6), round(sta_list(r).time.getRate * 1e6)) * 1e-6;
                    pt0 = round(pt0*1e6)/1e6;
                    pt0 = ceil(pt0 / pr) * pr;
                    pt1 = floor(pt1 / pr) * pr;
                    
                    % return one p_time for each target
                    i = i + 1;
                    p_time(i) = GPS_Time(p_time_zero, (pt0 : pr : pt1));
                    p_time(i).toUnixTime();
                    
                    id_sync{i} = nan(p_time(i).length, numel(id));
                    for rs = id % for each rec to sync
                        if ~sta_list(rs).isEmpty && ~(obs_type(rs) == 0 && (rs ~= r)) % if it's not another different target
                            [~, id_ref, id_rec] = intersect(round(sta_list(rs).time.getRefTime(p_time_zero) * 1e1)/1e1, (pt0 : pr : pt1));
                            id_sync{i}(id_rec, rs) = id_ref;
                        end
                    end
                end
            end
        end
        
        function data_s = smoothMat(data_in, method, spline_base)
            % Smooth a matrix column by column using method
            %   "spline"
            %   "poly_quad"
            %
            % SYNTAX
            %   smoothMat(data_in, method, spline_base)
            if nargin < 3 || isempty(spline_base)
                spline_base = 10; % 5 min
            end
            if nargin < 2 || isempty(method)
                method = 'spline';
            end
            data_s = Receiver_Commons.smoothSatData([],[],data_in, [], method, spline_base);
        end
        
        function data_s = smoothSatData(data_az, data_el, data_in, cs_mat, method, spline_base, max_gap)
            % Smooth a matrix column by column using method
            %   "spline"
            %   "poly_quad"
            %
            % SYNTAX
            %  data_s = smoothMat(data_az, data_el, data_in, cs_mat, method, spline_base, max_gap)
            if nargin < 4
                cs_mat = [];
            end
            if nargin < 6 || isempty(spline_base)
                spline_base = 10; % 5 min
            end
            if nargin < 5 || isempty(method)
                method = 'spline';
            end
            if nargin < 7 || isempty(max_gap)
                max_gap = 0;
            end
            if strcmp(method,'spline')
                data_s = data_in;
                for s = 1 : size(data_s, 2)
                    if isempty(cs_mat)
                        lim = getOutliers(~isnan(data_s(:,s)));
                    else
                        lim = getOutliers(~isnan(data_s(:,s)), cs_mat(:,s));
                    end
                    if max_gap > 0
                        lim = limMerge(lim, max_gap);
                    end
                    
                    % remove small intervals
                    %lim((lim(:,2) - lim(:,1)) < spline_base, :) = [];
                    for l = 1 : size(lim, 1)
                        arc_size = lim(l,2) - lim(l,1) + 1;
                        id_arc = lim(l,1) : lim(l,2);
                        id_arc(isnan(data_s(id_arc, s))) = [];
                        data_tmp = data_s(id_arc, s);
                        if length(data_tmp) > 3
                            data_s(id_arc, s) = splinerMat([], data_tmp, min(arc_size, spline_base));
                        end
                    end
                end
            elseif strcmp(method,'poly_quad')
                data_s = data_in;
                for s = 1 : size(data_s, 2)
                    lim = getOutliers(~isnan(data_s(:,s)), cs_mat(:,s));
                    % remove small intervals
                    lim((lim(:,2) - lim(:,1)) < spline_base, :) = [];
                    for l = 1 : size(lim, 1)
                        n_ep =  lim(l,2) - lim(l,1) +1;
                        if n_ep > 4
                            data_tmp = data_s(lim(l,1) : lim(l,2), s);
                            el_tmp = data_el(lim(l,1) : lim(l,2), s);
                            az_tmp = data_az(lim(l,1) : lim(l,2), s);
                            t = [1:n_ep]';
                            A = [ones(n_ep,1) t t.^2];
                            %                             Qxx = cholinv(A'*A);
                            par = A\data_tmp;
                            %                             par = Qxx*A'*data_tmp;
                            quad_mod = A*par;
                            
                            %                             data_coll = data_tmp - quad_mod;
                            %                             s02 = mean(data_coll.^2) / (n_ep -3);
                            %                             Cy_haty_hat = s02*A'*Qxx;
                            %                             Cvv =  0;
                            data_s(lim(l,1) : lim(l,2), s) = splinerMat(te, data_coll , 1) + quad_mod;
                            data_s(lim(l,1) : lim(l,2), s) = quad_mod;
                        end
                    end
                end
            end
        end
    end
end
