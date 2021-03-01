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
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
%  Written by:        Andrea Gatti, Giulio Tagliaferro ...
%  Contributors:      Andrea Gatti, Giulio Tagliaferro ...
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
    % ==================================================================================================================================================
    %% CONSTANT PROTECTED
    % ==================================================================================================================================================
    properties (Constant, Access = protected)
        EMPTY_RES = struct('type', 0, 'time', GPS_Time, 'pr', [], 'ph', [], 'system', '', 'prn', [], 'obs_code', '')
    end
    
    % ==================================================================================================================================================
    %% OBJ POINTERS
    % ==================================================================================================================================================
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
        coo
        %         = struct( ...
        %             'coo',  [], ...    % additional estimated coo
        %             'time', [], ...    % time of the coo
        %             'rate', [] ...    % rate of the coo
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
        
        function is_fixed = isFixed(this)
            % Is the position of this station fixed?
            %
            % SYNTAX
            %   is_fixed = this.isFixed()
            rf = Core.getReferenceFrame;
            is_fixed = rf.isFixed(this.parent.getMarkerName4Ch);
        end
        
        function is_fixed = isFixedPrepro(this)
            % Is the position of this station fixed for pre-processing?
            %
            % SYNTAX
            %   is_fixed = this.isFixedPrepro()
            rf = Core.getReferenceFrame;
            is_fixed = rf.isFixedPrepro(this.parent.getMarkerName4Ch);
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
            if is_empty
                if not(isempty(this.coo))
                    is_empty =  this.coo.getTime.length() == 0;
                end
            end
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
        
        function coo = getPos(this, flag_add_coo)
            % return the positions computed for the receiver
            % as Coordinates object
            %
            % OUTPUT
            %   coo     coordinates object
            %
            % SYNTAX
            %   coo = this.getPos()
            
            if nargin == 1 || isempty(flag_add_coo)
                flag_add_coo = 0;
            end
            
            if flag_add_coo == -100
                % Get it from PrePro (or old PPP)
                if ~isempty(this.xyz)
                    coo = Coordinates.fromXYZ(this.xyz);
                elseif ~isempty(this.parent.work.xyz)
                    coo = Coordinates.fromXYZ(this.parent.work.xyz);
                else
                    coo = Coordinates.fromXYZ([0 0 0]);
                end
                if coo.time.isEmpty
                    tmp = this.getPositionTime();
                    if not(isempty(tmp))
                        coo.setTime(tmp);
                    end
                end
                coo.setRate(Core.getState.sss_duration);
                coo.info.n_epo = this.quality_info.n_epochs;
                coo.info.n_obs = this.quality_info.n_obs;
                coo.info.s0 = this.quality_info.s0;
                coo.info.s0_ip = this.quality_info.s0_ip;
                coo.info.flag = 0;
                coo.info.fixing_ratio = 0;
                if isempty(coo.info.s0_ip) || all(isnan(coo.info.s0_ip))
                    coo.info.coo_type = 'G';
                else
                    coo.info.coo_type = iif(this.isFixed || this.isFixedPrepro, 'F', 'G');
                end
                if coo.info.coo_type == 'F'
                    coo.Cxx = zeros(3); % Fixed covariance = 0
                else
                    coo.Cxx = nan(3); % Unknown coordinate covariance NaN
                end
                coo.info.rate = this.time.getRate;
                coo.info.master_name = categorical({this.parent.getMarkerName4Ch});
            elseif flag_add_coo == 0
                coo = this.coo;
                % Fallback (in case of empty coo)
                if isempty(coo)
                    coo = this.getPos(-100);
                end
            else
                if ~isempty(this.add_coo)
                    coo = this.add_coo(min(numel(this.add_coo), flag_add_coo)).coo;
                    if coo.time.isEmpty
                        coo.setTime(this.add_coo(min(numel(this.add_coo), flag_add_coo)).coo.time);
                    end
                else
                    log.addWarning(sprintf('No additional coordinates are present into %s', this.parent.getMarkerName4Ch));
                    coo = this.getPos();
                    if coo.time.isEmpty
                        tmp = this.getPositionTime();
                        if not(isempty(tmp))
                            coo.setTime(tmp);
                        end
                    end
                end
            end
            coo.setName(this.parent.getMarkerName4Ch, this.parent.getMarkerName);
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
            
            % this is completely wrong
            fprintf('Todo!!!');
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
                xyz = median(xyz, 1,'omitnan');
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
           
        function [mfh, mfw, grad_term] = getSlantMFGen(this, id_sync)
            % Get Mapping function for the satellite slant
            %
            % OUTPUT
            %   mfh: hydrostatic mapping function
            %   mfw: wet mapping function
            %
            % SYNTAX
            %   [mfh, mfw] = this.getSlantMF()            
            
            n_sat = this.parent.getMaxSat;
            if nargin == 1 || isempty(id_sync)
                try
                    id_sync = this.id_sync;
                catch
                    id_sync = [];
                end
            end
            if n_sat == 0
                mfh = [];
                mfw = [];
                grad_term = [];
            else
                if this.length > 0
                    atmo = Core.getAtmosphere();
                    [lat, lon, h_ellipse, h_ortho] = this.getMedianPosGeodetic();
                    lat = median(lat);
                    lon = median(lon);
                    h_ortho = median(h_ortho);
                    if ~isempty(this)
                        if this.state.mapping_function == Prj_Settings.MF_NIEL
                            [mfh, mfw] = atmo.niell(this.time, lat./180*pi, zero2nan(this.sat.el)./180*pi, h_ellipse);
                        elseif this.state.mapping_function == Prj_Settings.MF_VMF1
                            [mfh, mfw] = atmo.vmf_grd(this.time, lat./180*pi, lon./180*pi, (this.sat.el)./180*pi, h_ellipse,1);
                        elseif this.state.mapping_function == Prj_Settings.MF_VMF3_1
                            [mfh, mfw] = atmo.vmf_grd(this.time, lat./180*pi, lon./180*pi, (this.sat.el)./180*pi, h_ellipse,3);
                        elseif this.state.mapping_function == Prj_Settings.MF_VMF3_5
                            [mfh, mfw] = atmo.vmf_grd(this.time, lat./180*pi, lon./180*pi, (this.sat.el)./180*pi, h_ellipse,3);
                        elseif this.state.mapping_function == Prj_Settings.MF_GMF
                            [mfh, mfw] = atmo.gmf(this.time, lat./180*pi, lon./180*pi, h_ortho, zero2nan(this.sat.el)./180*pi);
                        end
                       
                        if ~isempty(id_sync)
                            mfh = mfh(id_sync, :);
                            mfw = mfw(id_sync, :);
                        end
                        
                        if nargout > 2
                            if this.state.mapping_function_gradient == 1
                                grad_term = Atmosphere.chenHerringGrad(zero2nan(this.sat.el)./180*pi);
                            elseif this.state.mapping_function_gradient == 2
                                grad_term = Atmosphere.macmillanGrad(zero2nan(this.sat.el)./180*pi);
                                % To be used with Hydrostatic mapping function and temperature fixed to 15 Celsius
                                % from Herring 1992 Modeling atmospheric delays in the analysis of space geodetic data (page 171, eq 8)
                                [mfh_grad] = atmo.herring(lat./180*pi, (this.sat.el)./180*pi, h_ellipse, 15);
                                grad_term = grad_term .* mfh_grad;
                            else
                                grad_term = [];
                            end
                            if ~isempty(id_sync) && ~isempty(grad_term)
                                grad_term = grad_term(id_sync, :);
                            end
                        end
                    end
                end
            end
        end
        
        function [slant_td, slant_wd, go_id] = getSlantTD(this, go_id)
            % Get the slant total delay
            %
            % SYNTAX
            %   [slant_td, go_id] = this.getSlantTD();
            
            id_sync  = this.getIdSync;
            
            slant_td = zeros(size(this.sat.el)); % Slant Total Delay
            slant_wd = zeros(size(this.sat.el)); % Slant Wet Delay
            
            zwd = nan2zero(this.getZwd);
            apr_zhd = single(this.getAprZhd);
            cosaz = zero2nan(cosd(this.sat.az(id_sync, :)));
            sinaz = zero2nan(sind(this.sat.az(id_sync, :)));
          
            [gn ,ge] = this.getGradient();
            if ~isempty(this.tzer)
                [mfh, mfw, grad_term] = this.getSlantMF();
            else
                if ~isempty(zwd)
                    if any(ge(:))
                        [mfh, mfw, grad_term] = this.getSlantMF();
                    else
                        [mfh, mfw] = this.getSlantMF();
                    end
                else
                    mfh = [];
                    mfw = [];
                end
            end

            if nargin < 2 || isempty(go_id) || strcmp(go_id, 'all')
                this.log.addMessage(this.log.indent('Updating tropospheric errors'))
                try
                    go_id = unique(this.go_id)';
                catch
                    % If it fails I'm in receiver output
                    go_id = 1 : size(mfh, 2);
                end
            else
                go_id = serialize(go_id)';
            end            
            
            % Computing delays
            if any(mfh(:)) && any(mfw(:))
                if any(ge(:))
                    for g = go_id
                        slant_wd(id_sync, g) = mfw(:,g) .* zwd + gn .* grad_term(:,g) .* cosaz(:,g) + ge .* grad_term(:,g) .* sinaz(:,g);
                        slant_td(id_sync, g) = mfh(:,g) .* apr_zhd + slant_wd(id_sync, g);
                    end
                else
                    for g = go_id
                        slant_wd(id_sync, g) = mfw(:,g).*zwd;
                        slant_td(id_sync, g) = mfh(:,g) .* apr_zhd + slant_wd(id_sync, g);
                    end
                end
                % If Zernike have been used rebuild delays in this way:
                if ~isempty(this.tzer)
                    degree = ceil(-3/2 + sqrt(9/4 + 2*(Prj_Settings.getNumZerTropoCoef  -1)));
                    el = this.sat.el(id_sync,:);
                    
                    rho = (90 - el)/(90);
                    zer_val = nan(size(this.sat.err_tropo,1),size(this.sat.err_tropo,2),Prj_Settings.getNumZerTropoCoef);
                    
                    az = this.sat.az(id_sync,:);
                    idx = el>0;
                    for  g = go_id
                        zer_val(idx(:,g),g,:) = Core_Utils.getAllZernike(degree, az(idx(:,g),g)/180*pi, rho(idx(:,g),g));
                    end
                    for i = 4 : Prj_Settings.getNumZerTropoCoef
                        slant_td(id_sync,:) = slant_td(id_sync,:) + zero2nan(repmat(this.tzer(:,i-3),1,size(el,2)).*grad_term.*zer_val(:,:,i));
                    end
                end
            end
        end
                
        function [slant_td] = getSlantTD_AzEl(this, time, az , el )
            % Get the slant total delay
            %
            % SYNTAX
            %   [slant_td, go_id] = this.getSlantTD();
            slant_td = zeros(size(az)); % Slant Total Delay
            
            zwd = interp1(this.getTime.getMatlabTime,this.getZwd,time.getMatlabTime);
            apr_zhd = interp1(this.getTime.getMatlabTime,this.getAprZhd,time.getMatlabTime);
            [gn ,ge] = this.getGradient();
            gn = interp1(this.getTime.getMatlabTime,gn,time.getMatlabTime);
            ge = interp1(this.getTime.getMatlabTime,ge,time.getMatlabTime);
            if any(ge(:))
                [mfh, mfw, grad_term] = this.getSlantMF();
            else
                [mfh, mfw] = this.getSlantMF();
            end
            
            cosaz = zero2nan(cosd(az));
            sinaz = zero2nan(sind(az));
            state = Core.getCurrentSettings();
            atmo = Core.getAtmosphere();
            this.updateCoordinates();
            
            if this.state.mapping_function == 3
                [mfh, mfw] = atmo.niell(time, this.lat./180*pi, el/180*pi,this.h_ellips);
            elseif this.state.mapping_function == 2
                [mfh, mfw] = atmo.vmf_grd(time, this.lat./180*pi, this.lon./180*pi, el./180*pi, this.h_ellips);
            elseif this.state.mapping_function == 1
                [mfh, mfw] = atmo.gmf(time, this.lat./180*pi, this.lon./180*pi, this.h_ortho, el./180*pi);
            end
            
            if any(ge(:))
                slant_td = mfw .* zwd + gn .* grad_term .* cosaz + ge .* grad_term .* sinaz;
            else
                slant_td = mfw .* zwd;
            end
            slant_td = mfh .* apr_zhd + slant_td;
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
        
        function [zwd, p_time, id_sync, tge, tgn] = getZwd_mr(sta_list)
            % Get synced data of zwd
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [zwd, p_time, id_sync] = this.getZwd_mr()
            %  [zwd, p_time, id_sync, tge, tgn] = this.getZwd_mr()
            
            [p_time, id_sync] = Receiver_Commons.getSyncTimeExpanded(sta_list);
            
            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);
            
            n_rec = numel(sta_list);
            zwd = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).zwd)) = nan;
                zwd(~isnan(id_rec), r) = sta_list(r).zwd(id_rec(~isnan(id_rec)));
            end
            
            if nargout == 5
                tge = nan(size(id_sync));
                tgn = nan(size(id_sync));
                for r = 1 : n_rec
                    id_rec = id_sync(:,r);
                    id_rec(id_rec > length(sta_list(r).zwd)) = nan;
                    tge(~isnan(id_rec), r) = sta_list(r).tge(id_rec(~isnan(id_rec)));
                    tgn(~isnan(id_rec), r) = sta_list(r).tgn(id_rec(~isnan(id_rec)));
                end
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
                                lim = getFlagsLimits(id_ok);
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
                if numel(this.quality_info.n_spe.A) < max(this.getIdSync)
                    n_sat = nan;
                    n_sat_ss.G = [];
                else
                    n_sat = this.quality_info.n_spe.A(this.getIdSync);
                    for sys_c = cc.getActiveSysChar()
                        if ~isempty(this.quality_info.n_spe.(sys_c))
                            n_sat_ss.(sys_c) = this.quality_info.n_spe.(sys_c)(this.getIdSync);
                        else
                            n_sat_ss.(sys_c) = [];
                        end
                    end
                end
            else
                n_sat = nan;
                n_sat_ss.G = [];
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
            
            if isempty(this.apr_zwd)
                apr_zwd = [];
                time = GPS_Time();
            else
                apr_zwd = this.apr_zwd(this.getIdSync);
                time = this.time.getEpoch(this.getIdSync);
            end
        end
        
        function [az, el] = getAzEl(this, go_id)
            % Get the azimuth and elevation (on valid id_sync)
            %
            % SYNTAX
            %   [az, el] = this.getAzEl();
            if nargin == 2
                az = this.getAz(go_id);
                el = this.getEl(go_id);
            else
                az = this.getAz();
                el = this.getEl();
            end
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
            
            id_sync = this.getIdSync;
            if isempty(id_sync)
                az = this.sat.az(:, go_id);
            else
                az = this.sat.az(this.getIdSync, go_id);
            end
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
            
            id_sync = this.getIdSync;
            if isempty(id_sync)
                el = this.sat.el(:, go_id);
            else
                el = this.sat.el(this.getIdSync, go_id);
            end
        end
        
        function [res, obs_code, prn] = getU1(this)
            % get residual
            %
            % SYNTAX
            %   [res, obs_code, prn] = this.getU1()
            res = this.sat.res.getU1;
            res = res(this.getIdSync(), :);
        end
        
        function [res, obs_code, prn] = getPhU2(this)
            % get phase residual
            %
            % SYNTAX
            %   res = this.getPhU2()
            [res, obs_code, prn] = this.sat.res.getPhU2;
            res = res(this.getIdSync(), :);
        end
        
        function [res, obs_code, prn] = getPrU2(this)
            % get code residual
            %
            % SYNTAX
            %   res = this.getPrU2()
            [res, obs_code, prn] = this.sat.res.getPrU2;
            res = res(this.getIdSync(), :);
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
    %% TESTER
    % ==================================================================================================================================================
    methods (Access = public, Static)
        function id_ko = checkTropo_mr(data, flag_reduce, thr_lev)
            % Check ZWD for possible outliers (works on the mr ZWD getter)
            %
            % INPUT
            %   tropo         matrix of delays [ n_obs x n_rec ]
            %   flag_reduce   reduce the entries using their correlated signal
            %   thr_lev       level of threshold (usually 5.5 - conservative)
            %
            % SYNTAX
            %   id_ko = Receiver_Commons.checkTropo_mr(tropo, flag_reduce)
            %
            % EXAMPLE
            %   % Prepare the data
            %   [zwd, p_time, id_sync] = sta_list.getZwd_mr;
            %   % Get outliers
            %   id_ko = Receiver_Commons.checkTropo_mr(zwd, true, thr_lev)
            
            % state = Core.getCurrentSettings; 
            if nargin < 2 || isempty(flag_reduce)
                flag_reduce = true;
            end
            if nargin < 3 || isempty(thr_lev)
                thr_lev = 5.5;
            end
            
            % Use moving median to remove outliers
            x = movmedian(Core_Utils.diffAndPred(data),5, 'omitnan');
            y = data;
            
            % I can reduce the input only if I  have enough data
            if size(y, 2) >= 2
                med = [0 cumsum(nan2zero(median(diff(y, 1, 2), 'omitnan')))];
                med(abs(med) > 0.5) = 0; % ZWD out of more than 0.5m is an outlier!
                y = bsxfun(@minus, y, med);
            end
            if size(y, 2) >= 2 && flag_reduce
                reduction = median(y, 2, 'omitnan');
                if size(y, 2) < 3
                    % Try to remove stronger spikes before it's too late
                    reduction = movmedian(reduction, 3);
                end
                id_ko = sum(not(isnan(y)), 2) <= 1;
                t = (1 : numel(reduction))';
                reduction(id_ko) = interp1q(t(not(id_ko)), reduction(not(id_ko)), t(id_ko));
                y = bsxfun(@minus, y, reduction);
            end
            id_ok = y < (6 * std(y(:), 'omitnan'));
            
            % Normalize y and its derivate x
            x = x(:) ./ std(x(id_ok),'omitnan');
            y = y(:) ./ std(y(id_ok),'omitnan');
            
            x = x ./ 1.5; % use a stronger treshold on the derivate
            
            % Recenter the outlier sensors
            x = x - median(x, 'omitnan');
            y = y - median(y, 'omitnan');
            
            y(y < 0) = y(y < 0) ./ 1.3; % use a stronger treshold on negative outliers
            
            % Test sensors against the treshold
            id_ko = reshape(hypot(nan2zero(x), y) > thr_lev, size(data, 1), size(data, 2));
            % expand flag to remove short arcs
            min_arc = 2; % max(1, round(state.getMinArc / 2))
            id_ko = ~isnan(data) & flagMerge(isnan(data) | id_ko, min_arc);
            
            % DEBUG: figure; hist(hypot(x,y), 1000);
            
            % DEBUG: figure; plot(x, y,'.'); axis equal
            % DEBUG: hold on; plot(0, 0, '*g');
            % DEBUG: hold on; plot(x(id_ko), y(id_ko), '.y');
            % DEBUG: ylim([-10 10]);
            % DEBUG: xlim([-10 10]);
            
            % DEBUG: figure; plotSep(data, '.-');
            % DEBUG: x = repmat((1:size(data, 1))', 1 ,size(data,2));
            % DEBUG: hold on; plot(x(id_ko), data(id_ko), '.', 'MarkerSize', 10, 'Color', 0.1*[0.9 0.9 0.9]);
            
            % DEBUG: data(id_ko) = nan;
            % DEBUG: figure; plotSep(data, '.-');
            % DEBUG: figure; imagesc(data)
        end
    end
    
    methods (Access = public)
        function err_status = cleanTropo(sta_list, thr_lev, flag_verbose) 
            % Check ZWD for possible outliers (works on the mr ZWD getter)
            % And remove all the outliers from the tropo parameters 
            %
            % INPUT
            %   thr_lev       level of threshold (usually 5.5 - conservative)
            %   flag_verbose  print cleaning stats
            %
            % SYNTAX
            %   id_ko = sta_list.cleanTropo(thr_lev, flag_reduce)
            
            if nargin < 2 || isempty(thr_lev)
                thr_lev = [];
            end
            
            if nargin < 3 || isempty(flag_verbose)
                flag_verbose = true;
            end
            
            err_status = 0; 
            if flag_verbose
                log = Core.getLogger;
            end
            
            % Prepare the data
            [tropo, p_time, id_sync] = sta_list.getZwd_mr;
            if any(tropo(:))
                % Get outliers
                id_ko = Receiver_Commons.checkTropo_mr(tropo, true, thr_lev);
                % Clean data Receiver by receiver
                if any(id_ko(:))
                    log.addMonoMessage(sprintf('\n--------------------------------------'));
                    log.addMonoMessage(sprintf(' Tropo cleaning (Based on ZWD values)'));
                    log.addMonoMessage(sprintf('\n--------------------------------------'));
                    log.addMonoMessage(sprintf('    GNSS    Valid   Outlier   Outlier'));
                    log.addMonoMessage(sprintf(' Station   epochs       num      perc'));
                    log.addMonoMessage(sprintf('--------------------------------------'));
                    for r = 1 : size(id_ko, 2)
                        log.addMonoMessage(sprintf('    %s %8d   %7d   %5.2f %%', sta_list(r).parent.getMarkerName4Ch, sum(not(isnan(tropo(:, r)))), sum(id_ko(:, r)), 100 * sum(id_ko(:, r)) / sum(not(isnan(tropo(:, r))))));  
                        % Get the id of the outlier for the current receiver
                        id_out = id_sync(id_ko(:,r),r);
                        if not(isempty(id_out))
                            sta_list(r).remEpoch(id_out);
                        end
                    end
                end
            else
                % No data available
                err_status = 1;
                if flag_verbose
                    log.addWarning('No valid ZWD found, cleaning tropo inactive')
                end
            end
            
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
                param_to_export = [ 1 1 1 1 1 0 0 0];
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
                            file_name = sprintf('%s',[rec.state.getOutDir() filesep rec.parent.getMarkerName4Ch sprintf('%03d', doy) sess_str '.' yy 'zpd']);
                            snx_wrt = SINEX_Writer(file_name);
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
                            rec(1).log.addStatusOk(sprintf('Tropo saved into: %s', file_name));
                        end
                    catch ex
                        Core_Utils.printEx(ex);
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
            
            log = Core.getLogger();
            for r = 1 : numel(this)
                if max(this(r).quality_info.s0) < 0.10
                    try
                        this(r).updateCoordinates;
                        time = this(r).getTime.getNominalTime();
                        state = Core.getCurrentSettings;

                        [~, t_start] = state.getSessionLimits();
                        t_start.toGps;
                        t_end = t_start.last;
                        t_start = t_start.first;
                        
                        [year, doy] = t_start.getDOY();
                        t_start_str = t_start.toString('HHMM');
                        t_start.toUtc();
                        time.toUtc();
                        
                        marker_name = this.parent.getMarkerName;
                        sys = this.parent.getActiveSys;
                        xyz = this.getPos.getXYZ;
                        [pressure, temperature, humidity] = this.getPTH();
                        
                        lat = this(r).lat; %#ok<NASGU>
                        lon = this(r).lon; %#ok<NASGU>
                        h_ellips = this(r).h_ellips; %#ok<NASGU>
                        h_ortho = this(r).h_ortho; %#ok<NASGU>
                        ztd = this(r).getZtd(); %#ok<NASGU>
                        zwd = this(r).getZwd(); %#ok<NASGU>
                        pwv = this(r).getPwv(); %#ok<NASGU>
                        utc_time = time.getMatlabTime; %#ok<NASGU>
                        
                        out_dir = fullfile(this(r).state.getOutDir(), sprintf('%4d', year), sprintf('%03d',doy));
                        if ~exist(out_dir, 'file')
                            mkdir(out_dir);
                        end
                        fname_old = sprintf('%s',[out_dir filesep this(r).parent.getMarkerName4Ch sprintf('%04d%03d_%4s_', year, doy, t_start_str) '*.m']);
                        old_file_list = dir(fname_old);

                        fname = sprintf('%s',[out_dir filesep this(r).parent.getMarkerName4Ch sprintf('%04d%03d_%4s_%d', year, doy, t_start_str,  round(t_end-t_start)+1) '.mat']);
                        save(fname, 'marker_name', 'sys', 'lat', 'lon', 'h_ellips', 'h_ortho', 'ztd', 'zwd', 'pwv', 'utc_time', 'xyz', 'pressure', 'temperature', 'humidity', '-v6');
                        
                        log.addStatusOk(sprintf('Tropo saved into: "%s"', fname));
                        for f = 1 : numel(old_file_list)
                            % If I did not overwrite the old file, delete it
                            if ~strcmp([out_dir filesep old_file_list(f).name], fname)
                                log.addStatusOk(sprintf('Old tropo file deleted: "%s"', [out_dir filesep old_file_list(f).name]));
                                delete([out_dir filesep old_file_list(f).name]);
                            end
                        end
                    catch ex
                        Core_Utils.printEx(ex);
                        log.addError(sprintf('saving Tropo in matlab format failed: "%s"', ex.message));
                    end
                else
                    if isempty(max(this(r).quality_info.s0))
                        log.addWarning(sprintf('s02 no solution have been found, station skipped'));
                    else
                        log.addWarning(sprintf('s02 (%f m) too bad, station skipped', max(this(r).quality_info.s0)));
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
            
            log = Core.getLogger();
            for r = 1 : numel(this)
                if max(this(r).quality_info.s0) < 0.10
                    try
                        this(r).updateCoordinates;
                        state = Core.getCurrentSettings;
                        time = this(r).getTime.getNominalTime;
                        [~, t_start] = state.getSessionLimits();
                        t_end = t_start.last;
                        t_start = t_start.first;
                        [year, doy] = t_start.getDOY();
                        t_start_str = t_start.toString('HHMM');
                        t_start.toUtc();
                        time.toUtc();
                        ztd = this(r).getZtd();
                        zwd = this(r).getZwd();
                        pwv = this(r).getPwv();
                        [gn, ge ] =  this(r).getGradient();
                        
                        out_dir = fullfile(this(r).state.getOutDir(), sprintf('%4d', year), sprintf('%03d',doy));
                        if ~exist(out_dir, 'file')
                            mkdir(out_dir);
                        end
                        
                        fname_old = sprintf('%s',[out_dir filesep this(r).parent.getMarkerName4Ch sprintf('%04d%03d_%4s_', year, doy, t_start_str) '*.csv']);
                        old_file_list = dir(fname_old);

                        fname = sprintf('%s',[out_dir filesep this(r).parent.getMarkerName4Ch sprintf('%04d%03d_%4s_%d', year, doy, t_start_str, round(t_end-t_start)+1) '.csv']);
                                                
                        n_data = time.length;
                        if n_data == 0
                            Core.getLogger.addError(sprintf('saving Tropo in csv format failed: "No exportable data found"'));
                        else
                            fid = fopen(fname,'Wb');
                            % Export any non empty parameter
                            fprintf(fid,'Date               %s%s%s%s%s\n', iif(isempty(zwd), '', ',ZTD [m]     '), iif(isempty(ztd), '', ',ZWD [m]     '), iif(isempty(pwv), '', ',PWV [mm]    '), iif(isempty(ge), '', ',GE [m]      '), iif(isempty(gn), '', ',GN [m]      '));
                            data = time.toString('dd/mm/yyyy HH:MM:SS');
                            if ~isempty(ztd)
                                data = [data reshape(sprintf(',%12.6f', ztd), 13, n_data)'];
                            end
                            if ~isempty(zwd)
                                data = [data reshape(sprintf(',%12.6f', zwd), 13, n_data)'];
                            end
                            if ~isempty(pwv)
                                data = [data reshape(sprintf(',%12.6f', pwv), 13, n_data)'];
                            end
                            if ~isempty(ge)
                                data = [data reshape(sprintf(',%12.6f', ge), 13, n_data)'];
                            end
                            if ~isempty(gn)
                                data = [data reshape(sprintf(',%12.6f', gn), 13, n_data)'];
                            end
                            data = [data char(10 .* ones(n_data,1))];
                            fprintf(fid,data');
                            fclose(fid);
                            log.addStatusOk(sprintf('Tropo saved into: "%s"', fname));
                            for f = 1 : numel(old_file_list)
                                % If I did not overwrite the old file, delete it
                                if ~strcmp([out_dir filesep old_file_list(f).name], fname)
                                    log.addStatusOk(sprintf('Old tropo file deleted: "%s"', [out_dir filesep old_file_list(f).name]));
                                    delete([out_dir filesep old_file_list(f).name]);
                                end
                            end
                        end
                    catch ex
                        Core_Utils.printEx(ex);
                        log.addError(sprintf('saving Tropo in csv format failed: "%s"', ex.message));
                    end
                else
                    if isempty(max(this(r).quality_info.s0))
                        log.addWarning(sprintf('s02 no solution have been found, station skipped'));
                    else
                        log.addWarning(sprintf('s02 (%f m) too bad, station skipped', max(this(r).quality_info.s0)));
                    end
                end
            end
        end
        
        function err_code = exportSlantDelay(this, time_rate, mode, format, flag_append, flag_wet)
            % Export Slant Total or Wet Delay
            % 
            % INPUT 
            %   time_rate   it could be a rate in seconds - default 21600 (6 hours)
            %               or an object of class GPS_Time
            %               (e.g. GPS_Time.fromString(2013/06/04 12:00:00))
            %
            %   mode        indicate if it is necessary to add the residuals
            %               0 - do not add residuals
            %               1 - add residuals 
            %               2 - add spline smoothed residuals (1 spline every 900 seconds)
            %
            %   format      'compact' - MJD SAT SWD EL AZ (default)
            %               'extended' - extended - Year DOY SOD SWD SAT SAT_X SAT_Y SAT_Z Receiver REC_X REC_Y REC_Z elAngle(rad) elAngle(deg) Azimuth(rad) Azimuth(deg)
            %
            %   flag_append if provided (true) go in append mode (default == true for compact mode, false for extended mode)
            %
            %   flag_wet    true if wet delay must be exported, otherwise export total delay
            %
            %
            % NOTE: at the moment adding residuals works only if the data have been porcessed 
            %       with the combined engine in iono-free mode
            %
            % SYNTAX
            %   this.exportSlantDelay(this, time_rate, mode, format, flag_append, flag_wet)
            
            err_code = 1;
            if nargin < 4 || isempty(format)
                format = 'compact';
            end
            if nargin < 3 || isempty(mode)
                mode = 0;
            end
            if nargin < 2 || isempty(time_rate)
                time_rate = 21600;
            end
            
            if nargin<5 || isempty(flag_append)
                if format(1) == 'c'
                    flag_append = true;
                else
                    flag_append = false;
                end
            end
            
            if nargin < 6 || isempty(flag_wet)
                flag_wet = true;
            end
            
            try
                log = Core.getLogger();
                cut_off = Core.getState.getCutOff;
                state = Core.getCurrentSettings;
                
                cc = Core.getConstellationCollector;
                go_id = [];
                if ~isempty(this.sat.res) && (isa(this.sat.res, 'Residuals') && ~this.sat.res.isEmpty)
                    go_id = cc.getIndex(this.sat.res.obs_code(:,1), this.sat.res.prn);
                end
                
                if isa(this, 'Receiver_Work_Space')
                    [~, t_start] = state.getSessionLimits();
                    t_end = t_start.last;
                    t_start = t_start.first;
                    
                    if isempty(go_id)
                        go_id = unique(this.go_id);
                    end
                else
                    t_start = state.sss_date_start;
                    t_end = state.sss_date_stop;
                    if isempty(go_id)
                        go_id = cc.index;
                    end
                end
                t_start_str = t_start.toString('HHMM');
                
                time_slant = this.getTime;
                if isempty(time_slant)
                    log.addWarning(sprintf('No solution found for %s', this.parent.getMarkerName4Ch));
                else
                    time0 = round(time_slant.first.getMatlabTime);
                    time_slant.changeRef(time0);  % Time referred to the beginning of the day;
                    time_slant = time_slant.getNominalTime;
                    if isa(time_rate, 'GPS_Time')
                        time_ref = time_rate;
                    else
                        time_ref = GPS_Time.fromRefTime(time0, unique(round(time_slant.getRefTime / time_rate)) * time_rate);
                    end
                    [slant_td, slant_wd, go_id] = this.getSlantTD(go_id);
                    if flag_wet
                        slant_d = slant_wd;
                    else
                        slant_d = slant_td;
                    end
                    slant_d = slant_d(:, go_id);
                    if ~isempty(time_slant) && ~isempty(slant_d)
                        sky = Core.getCoreSky;
                        coo = this.getPos;
                        time_pos = this.getPositionTime;
                        time_pos.changeRef(time0);
                        
                        core = Core.getCurrentCore;
                        sky = core.sky;
                        if isempty(state.eph_name)
                            fw = File_Wizard();
                            fw.conjureNavFiles(this.time.first, this.time.last);
                        end
                        lim = this.time.first.getCopy;
                        lim.append(this.time.last);
                        core.initSkySession(lim);
                        
                        [year, doy] = t_start.getDOY();
                        % Preparing the file_name
                        out_dir = fullfile(this.state.getOutDir(), sprintf('%4d', year), sprintf('%03d',doy));
                        if ~exist(out_dir, 'file')
                            mkdir(out_dir);
                        end
                        if format(1) == 'e'
                            file_path = sprintf('%s',[out_dir filesep this.parent.getMarkerName4Ch sprintf('_Slant%cD_%04d%03d_%4s_%d', iif(flag_wet, 'W', 'T'), year, doy, t_start_str, round(t_end-t_start)+1) '.SNX']);
                            fid = fopen(file_path, 'wb');
                        end
                        if format(1) == 'e' && fid <= 0
                            log.addError(sprintf('Writing "%s" seems impossible', file_path));
                        else
                            % Print header
                            if format(1) == 'e'
                                fprintf(fid, ['Year               Day           Seconds    Slant WV Delay(mm)     Satellite             SAT_X(Km)         SAT_Y(Km)         SAT_Z(Km)      ' ...
                                    'Receiver             REC_X(Km)         REC_Y(Km)         REC_Z(Km)  elAngle(rad)      elAngle(deg)      Azimuth(rad)      Azimuth(deg)']);
                            end
                            
                            % Get the data at the requested time
                            [az, el, sat_coo] = sky.getAzimuthElevation(coo.getMedianPos, time_ref, go_id);
                            
                            % Get residuals
                            if mode > 0
                                [res, obs_code_res, prn_res, time_res] = this.sat.res.getRangeResiduals();
                                go_id_res = cc.getIndex(obs_code_res(:,1), prn_res);
                                [go_id, idr, ids] = intersect(go_id_res, go_id);
                                
                                % Delete slants with no residuals (and res with no slant)
                                slant_d = slant_d(:, ids);
                                res = res(:, idr);
                                obs_code_res = obs_code_res(idr, :);
                                prn_res = prn_res(idr);
                                
                                if mode == 2
                                    res = Receiver_Commons.smoothMat(res, 'spline', 900 / this.getRate);
                                end
                            end
                            
                            sat_name = cc.getSatName(go_id);
                            err_code = time_ref.length;
                            for e = 1 : time_ref.length
                                
                                % Get id of the closest (in time) position
                                [time_dist, id_min] = min(abs(time_pos - time_ref.getEpoch(e)));
                                rec_xyz = coo.getElement(id_min).getXYZ;
                                
                                % Get id of the closest (in time) slant set
                                [time_dist, id_min] = min(abs(time_slant - time_ref.getEpoch(e)));
                                if time_dist > 1800 % (seconds)
                                    log.addWarning(sprintf('No solution found for %s at %s', this.parent.getMarkerName4Ch, time_ref.getEpoch(e).toString));
                                else
                                    slant_set = slant_d(id_min, :);
                                    slant_set(el(e, :) < cut_off) = nan;
                                    
                                    % Find residuals close to the required slant time
                                    if mode > 0
                                        [time_dist, id_min] = min(abs(time_res - time_ref.getEpoch(e)));
                                        if time_dist > 180 % (seconds)
                                            log.addWarning(sprintf('No residuals found for %s at %s', this.parent.getMarkerName4Ch, time_ref.getEpoch(e).toString));
                                        else
                                            slant_set = slant_set + double(zero2nan(res(id_min, :)));
                                        end
                                    end
                                    
                                    sat_ok = ~isnan(slant_set);
                                    
                                    [year, doy, sod] = time_ref.getEpoch(e).getDOY;
                                    
                                    str = '';
                                    if format(1) == 'c'
                                        for s = find(sat_ok)
                                            str = sprintf('%s%11.5f %4s %3s %8.4f %8.4f  %8.4f\n', ...
                                                str, time_ref.getEpoch(e).getMJD, this.parent.getMarkerName4Ch, sat_name(s,:), slant_set(s), el(e, s), mod(az(e, s), 360));
                                        end
                                    elseif format(1) == 'e'
                                        for s = find(sat_ok)
                                            sat_xyz = squeeze(sat_coo.xyz(e, s, :));
                                            str = sprintf(['%s%04d %03d %05d %.6f %3s %.6f %.6f %.6f ' ...
                                                '%4s %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n'], ...
                                                str, year, doy, sod, slant_set(s)*1e3, sat_name(s,:), sat_xyz(1)*1e-3, sat_xyz(2)*1e-3, sat_xyz(3)*1e-3, ...
                                                this.parent.getMarkerName4Ch, rec_xyz(1)*1e-3, rec_xyz(2)*1e-3, rec_xyz(3)*1e-3, el(e, s)/180*pi, el(e, s), mod(az(e, s), 360)/180*pi, mod(az(e, s), 360));
                                        end
                                    end
                                    if format(1) == 'c'
                                        file_path = sprintf('%s', fullfile(out_dir, sprintf('CE_%04d%03d%02d_hydr_gradients.S%cD', year, doy, round(sod/3600), iif(flag_wet, 'W', 'T'))));
                                        if flag_append
                                            fid = fopen(file_path, 'ab');
                                        else
                                            fid = fopen(file_path, 'wb');
                                        end
                                        if fid <= 0
                                            log.addError(sprintf('Writing "%s" seems impossible', file_path));
                                        else
                                            fprintf(fid, '%s', str);
                                            fclose(fid);
                                        end
                                        log.addStatusOk(sprintf('Slant %s delay for %s %s into: "%s"', iif(flag_wet, 'WV', 'Total'), this.parent.getMarkerName4Ch, iif(flag_append, 'appended', 'saved'), file_path));
                                    else
                                        fprintf(fid, '%s', str);
                                    end
                                end
                                err_code = err_code - 1;
                            end
                            if format(1) == 'e'
                                fclose(fid);
                                log.addStatusOk(sprintf('Slant %s delay saved into: "%s"', iif(flag_wet, 'WV', 'Total'), file_path));
                                err_code = 0;
                            end                            
                        end
                    end
                end
            catch ex
                Core_Utils.printEx(ex);
                log.addError(sprintf('saving slants in text format failed for %s: "%s"', this.parent.getMarkerName4Ch, ex.message));
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
            this.showResSkyPolarScatter();
            this.showResSkyCartScatter();
            this.showOutliersAndCycleSlip();
            this.showOutliersAndCycleSlip_p();
            dockAllFigures();
        end
        
        
        function  fh_list = showBaselineENU(recs, baseline_ids, flag_add_coo, n_obs)
            % Function to plot baseline between 2 or more stations
            %
            % INPUT:
            %   sta_list                 list of receiver commons objects
            %   baseline_ids/ref_id      n_baseline x 2 - couple of id in sta_list to be used
            %   flag_add_coo             use external / internal / high rate cordinates
            %   n_obs                    use only the last n_obs
            %   
            % SYNTAX
            %   showBaselineENU(sta_list, <baseline_ids = []>, flag_add_coo, n_obs))
            if nargin < 2 || isempty(baseline_ids)
                % remove empty receivers
                if flag_add_coo >= 0
                    recs = recs(~recs.isEmpty_mr);
                end
                
                n_rec = numel(recs);
                baseline_ids = GNSS_Station.getBaselineId(n_rec);
            end
            baseline_ids = baseline_ids(~any(isnan(baseline_ids)'),:);
            if numel(baseline_ids) == 1
                n_rec = numel(recs);
                ref_rec = setdiff((1 : n_rec)', baseline_ids);
                baseline_ids = [baseline_ids * ones(n_rec - 1, 1), ref_rec];
            end
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo ~= 0)
                flag_add_coo = 0;
            end
            
            if flag_add_coo >= 0
                recs = recs(~recs.isEmpty_mr);
            end
            log = Core.getLogger();
            
            fh_list =  [];
            for b = 1 : size(baseline_ids, 1)
                
                % Select input data
                flag_ready = false;
                if flag_add_coo == 0
                    % Use internal coordinates
                    coo_ref = recs(baseline_ids(b, 1)).getPos();
                    if coo_ref.time.isEmpty
                        coo_ref.setTime(recs(baseline_ids(b, 1)).getPositionTime());
                    end
                    coo_trg = recs(baseline_ids(b, 2)).getPos();
                    if coo_trg.time.isEmpty
                        coo_trg.setTime(recs(baseline_ids(b, 2)).getPositionTime());
                    end
                    flag_ready = true;
                elseif flag_add_coo > 0
                    % Use additional coordinates
                    if ~isempty(recs(baseline_ids(b, 1)).add_coo) &&  ~isempty(recs(baseline_ids(b, 2)).add_coo)
                        coo_ref = recs(baseline_ids(b, 1)).add_coo(min(numel(recs(baseline_ids(b, 1)).add_coo), flag_add_coo)).coo;
                        coo_trg = recs(baseline_ids(b, 2)).add_coo(min(numel(recs(baseline_ids(b, 2)).add_coo), flag_add_coo)).coo;
                        flag_ready = true;
                    else
                        log.addWarning(sprintf('No additional coordinates are present into %s or %s', recs(baseline_ids(b, 1)).parent.getMarkerName4Ch, recs(baseline_ids(b, 2)).parent.getMarkerName4Ch));
                    end
                    if isempty(coo_ref) || isempty(coo_trg)
                        log.addWarning(sprintf('No additional coordinates are present into %s or %s', recs(baseline_ids(b, 1)).parent.getMarkerName4Ch, recs(baseline_ids(b, 2)).parent.getMarkerName4Ch));
                        flag_ready = false;
                    end
                else
                    % Use coordinates from coo files
                    % If they exist
                    coo_path = recs(baseline_ids(b, 1)).getPos.getCooOutPath();
                    if exist(coo_path, 'file') == 2
                        coo_ref = Coordinates.fromCooFile(coo_path);
                        flag_ready = true;
                    else
                        log.addWarning(sprintf('Missing reference coordinate file: "%s"', coo_path));
                        coo_ref = Coordinates();
                    end
                    
                    coo_path = recs(baseline_ids(b, 2)).getPos.getCooOutPath();
                    if exist(coo_path, 'file') == 2
                        coo_trg = Coordinates.fromCooFile(coo_path);
                    else
                        log.addWarning(sprintf('Missing target coordinate file: "%s"', coo_path));
                        coo_trg = Coordinates();
                        flag_ready = false;
                    end
                end
                
                if flag_ready
                    flag_ready = coo_ref.length > 1 && coo_trg.length > 1;
                    if not(flag_ready)
                        log.addError('This plot is not available without a series of data (more than one point)')
                    end
                end
                
                % Effective plot of the baseline
                if flag_ready
                    fh = coo_trg.showCoordinatesENU(coo_ref, n_obs);
                    set(0, 'CurrentFigure', fh);
                    ax = fh.Children(end);
                    set(fh, 'CurrentAxes', ax);
                    bsl_str = [recs(baseline_ids(b, 2)).parent.getMarkerName4Ch ' - ' recs(baseline_ids(b, 1)).parent.getMarkerName4Ch];
                    fig_name = sprintf('BSL_EN_U_%s-%s_%s', recs(baseline_ids(b, 1)).parent.getMarkerName4Ch,recs(baseline_ids(b, 2)).parent.getMarkerName4Ch, recs(baseline_ids(b, 1)).getTime.first.toString('yyyymmdd_HHMM'));
                    fh.UserData = struct('fig_name', fig_name);
                    fh.Name = ['dENU ' bsl_str];
                    Core_UI.addExportMenu(fh);
                    drawnow
                    fh_list = [fh_list; fh];
                end
            end
        end

        function  fh_list = showBaselinePlanarUp(recs, baseline_ids, flag_add_coo, n_obs)
            % Function to plot baseline between 2 or more stations
            %
            % INPUT:
            %   sta_list                 list of receiver commons objects
            %   baseline_ids/ref_id      n_baseline x 2 - couple of id in sta_list to be used
            %   flag_add_coo             use external / internal / high rate cordinates
            %   n_obs                    use only the last n_obs
            %   
            % SYNTAX
            %   showBaselinePlanarUp(sta_list, <baseline_ids = []>, flag_add_coo, n_obs)
            if nargin < 2 || isempty(baseline_ids)
                % remove empty receivers
                recs = recs(~recs.isEmpty_mr);
                
                n_rec = numel(recs);
                baseline_ids = GNSS_Station.getBaselineId(n_rec);
            end

            if numel(baseline_ids) == 1
                n_rec = numel(recs);
                ref_rec = setdiff((1 : n_rec)', baseline_ids);
                baseline_ids = [baseline_ids * ones(n_rec - 1, 1), ref_rec];
            end
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo ~= 0)
                flag_add_coo = 0;
            end
            if nargin < 4 || isempty(n_obs)
                n_obs = 0;
            end
           
            recs = recs(~recs.isEmpty_mr);
            log = Core.getLogger();
            
            fh_list =  [];
            for b = 1 : size(baseline_ids, 1)
                
                % Select input data
                flag_ready = false;
                if flag_add_coo == 0
                    % Use internal coordinates
                    coo_ref = recs(baseline_ids(b, 1)).getPos();
                    if coo_ref.time.isEmpty
                        coo_ref.setTime(recs(baseline_ids(b, 1)).getPositionTime());
                    end
                    coo_trg = recs(baseline_ids(b, 2)).getPos();
                    if coo_trg.time.isEmpty
                        coo_trg.setTime(recs(baseline_ids(b, 2)).getPositionTime());
                    end
                    flag_ready = true;
                elseif flag_add_coo > 0
                    % Use additional coordinates
                    if ~isempty(recs(baseline_ids(b, 1)).add_coo) &&  ~isempty(recs(baseline_ids(b, 2)).add_coo)
                        coo_ref = recs(baseline_ids(b, 1)).add_coo(min(numel(recs(baseline_ids(b, 1)).add_coo), flag_add_coo)).coo;
                        coo_trg = recs(baseline_ids(b, 2)).add_coo(min(numel(recs(baseline_ids(b, 2)).add_coo), flag_add_coo)).coo;
                        flag_ready = true;
                    else
                        log.addWarning(sprintf('No additional coordinates are present into %s or %s', recs(baseline_ids(b, 1)).parent.getMarkerName4Ch, recs(baseline_ids(b, 2)).parent.getMarkerName4Ch));
                    end
                    if isempty(coo_ref) || isempty(coo_trg)
                        log.addWarning(sprintf('No additional coordinates are present into %s or %s', recs(baseline_ids(b, 1)).parent.getMarkerName4Ch, recs(baseline_ids(b, 2)).parent.getMarkerName4Ch));
                        flag_ready = false;
                    end
                else
                    % Use coordinates from coo files
                    % If they exist
                    coo_path = recs(baseline_ids(b, 1)).getCooOutPath();
                    if exist(coo_path, 'file') == 2
                        coo_ref = Coordinates.fromCooFile(coo_path);
                        flag_ready = true;
                    else
                        log.addWarning(sprintf('Missing reference coordinate file: "%s"', coo_path));
                        coo_ref = Coordinates();
                    end
                    
                    coo_path = recs(baseline_ids(b, 2)).getCooOutPath();
                    if exist(coo_path, 'file') == 2
                        coo_trg = Coordinates.fromCooFile(coo_path);
                    else
                        log.addWarning(sprintf('Missing target coordinate file: "%s"', coo_path));
                        coo_trg = Coordinates();
                        flag_ready = false;
                    end
                end
                
                if flag_ready
                    flag_ready = coo_ref.length > 1 && coo_trg.length > 1;
                    if not(flag_ready)
                        log.addError('This plot is not available without a series of data (more than one point)')
                    end
                end
                
                % Effective plot of the baseline
                if flag_ready
                    fh = coo_trg.showCoordinatesPlanarUp(coo_ref, n_obs);
                    set(0, 'CurrentFigure', fh);
                    ax = fh.Children(end).Children(end).Children;
                    set(fh, 'CurrentAxes', ax);
                    bsl_str = [recs(baseline_ids(b, 2)).parent.getMarkerName4Ch ' - ' recs(baseline_ids(b, 1)).parent.getMarkerName4Ch];
                    if numel(ax.Title.String) == 2
                        ax.Title.String{1} = [bsl_str  ax.Title.String{1}(9:end)];
                    else
                        ax.Title.String{1} = bsl_str;
                    end
                    ax.Title.String{1} = [bsl_str  ax.Title.String{1}(9:end)];
                    fig_name = sprintf('BSL_EN_U_%s-%s_%s', recs(baseline_ids(b, 1)).parent.getMarkerName4Ch,recs(baseline_ids(b, 2)).parent.getMarkerName4Ch, recs(baseline_ids(b, 1)).getTime.first.toString('yyyymmdd_HHMM'));
                    fh.UserData = struct('fig_name', fig_name);
                    fh.Name = ['dENU ' bsl_str];
                    Core_UI.addExportMenu(fh);
                    drawnow
                    fh_list = [fh_list; fh];
                end
            end
        end

        function fh_list = showPositionENU(this, flag_add_coo, n_obs)
            % Plot East North Up coordinates of the receiver
            %
            % INPUT
            %   flag_add_coo             use external / internal / high rate cordinates
            %   n_obs                    use only the last n_obs
            %
            % SYNTAX 
            %   this.plotPositionENU(flag_add_coo, flag_add_coo, n_obs)
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo ~= 0)
                flag_add_coo = 0;
            end
            
            if nargin < 3 || isempty(n_obs)
                n_obs = 0;
            end
            
            log = Core.getLogger();
            fh_list = [];
            for r = 1 : numel(this)
                rec = this(r);
                if ~isempty(rec)
                    xyz = rec.getPosXYZ();
                    if size(xyz, 1) > 1 || flag_add_coo ~= 0
                        log.addMessage('Plotting positions');
                        
                        if flag_add_coo == 0
                            coo = rec.getPos();
                            if coo.time.isEmpty
                                coo.setTime(rec.getPositionTime());
                            end
                        elseif flag_add_coo > 0
                            if ~isempty(rec.add_coo)
                                coo = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo;
                                if coo.time.isEmpty
                                    coo.setTime(rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo.time);
                                end
                            else
                                log.addWarning(sprintf('No additional coordinates are present into %s', rec.parent.getMarkerName4Ch));
                                coo = rec.getPos();
                                if coo.time.isEmpty
                                    coo.setTime(rec.getPositionTime());
                                end
                            end
                        else
                            % Use coordinates from coo files
                            % If they exist
                            coo_path = rec.getPos.getCooOutPath();
                            if exist(coo_path, 'file') == 2
                                coo = Coordinates.fromCooFile(coo_path);
                            else
                                log.addWarning(sprintf('Missing coordinate file: "%s"', coo_path));
                                coo = Coordinates();
                            end
                        end
                        
                        if not(coo.isEmpty)
                            fh = coo.showCoordinatesENU([], n_obs);
                            set(0, 'CurrentFigure', fh);
                            ax = subplot(3,1,1);
                            ax.Title.String{1} = [rec.parent.getMarkerName4Ch  ax.Title.String{1}(9:end)];
                            fig_name = sprintf('ENU_at%gs_%s_%s', coo.time.getRate, rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                            fh.UserData = struct('fig_name', fig_name);
                            Core_UI.addExportMenu(fh);
                            fh_list = [fh_list; fh]; %#ok<AGROW>
                        end
                    else
                        Core.getLogger.addWarning(sprintf('%s - Plotting a single point static position is not yet supported', rec.parent.getMarkerName4Ch));
                    end
                end
            end
        end

        function fh_list = showPositionPlanarUp(this, flag_add_coo, n_obs)
            % Plot East North Up coordinates of the receiver
            %
            % SYNTAX 
            %   this.plotPositionENU(flag_add_coo);            
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo ~= 0)
                flag_add_coo = 0;
            end
            
            if nargin < 3 || isempty(n_obs)
                n_obs = 0;
            end
            
            log = Core.getLogger();
            fh_list = [];
            for r = 1 : numel(this)
                rec = this(r);
                if ~isempty(rec)
                    xyz = rec.getPosXYZ();
                    if size(xyz, 1) > 1 || flag_add_coo ~= 0
                        log.addMessage('Plotting positions');
                        
                        if flag_add_coo == 0
                            coo = rec.getPos();
                            if coo.time.isEmpty
                                coo.setTime(rec.getPositionTime());
                            end
                        elseif flag_add_coo > 0
                            if ~isempty(rec.add_coo)
                                coo = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo;
                                if coo.time.isEmpty
                                    coo.setTime(rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo.time);
                                end
                            else
                                log.addWarning(sprintf('No additional coordinates are present into %s', rec.parent.getMarkerName4Ch));
                                coo = rec.getPos();
                                if coo.time.isEmpty
                                    coo.setTime(rec.getPositionTime());
                                end
                            end
                        else
                            % Use coordinates from coo files
                            % If they exist
                            coo_path = rec.getPos.getCooOutPath();
                            if exist(coo_path, 'file') == 2
                                coo = Coordinates.fromCooFile(coo_path);
                            else
                                log.addWarning(sprintf('Missing coordinate file: "%s"', coo_path));
                                coo = Coordinates();
                            end
                        end
                        if not(coo.isEmpty)
                            fh = coo.showCoordinatesPlanarUp([], n_obs);
                            set(0, 'CurrentFigure', fh);
                            ax = fh.Children(end).Children(end).Children(1);
                            ax.Title.String{1} = [rec.parent.getMarkerName4Ch  ax.Title.String{1}(9:end)];
                            fig_name = sprintf('PUP_at%gs_%s_%s', coo.time.getRate, rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                            fh.UserData = struct('fig_name', fig_name);
                            Core_UI.addExportMenu(fh);
                            fh_list = [fh_list; fh]; %#ok<AGROW>
                        end
                    else
                        Core.getLogger.addWarning(sprintf('%s - Plotting a single point static position is not yet supported', rec.parent.getMarkerName4Ch));
                    end
                end
            end
        end

        function fh_list = showPositionXYZ(this, flag_add_coo, n_obs)
            % Plot X Y Z coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionXYZ();            
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo ~= 0)
                flag_add_coo = 0;
            end
            
            if nargin < 3 || isempty(n_obs)
                n_obs = 0;
            end
            
            log = Core.getLogger();
            fh_list = [];
            for r = 1 : numel(this)
                rec = this(r);
                if ~isempty(rec)
                    xyz = rec.getPosXYZ();
                    if size(xyz, 1) > 1 || flag_add_coo ~= 0
                        log.addMessage('Plotting positions');
                        
                        if flag_add_coo == 0
                            coo = rec.getPos();
                            if coo.time.isEmpty
                                coo.setTime(rec.getPositionTime());
                            end
                        elseif flag_add_coo > 0
                            if ~isempty(rec.add_coo)
                                coo = rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo;
                                if coo.time.isEmpty
                                    coo.setTime(rec.add_coo(min(numel(rec.add_coo), flag_add_coo)).coo.time);
                                end
                            else
                                log.addWarning(sprintf('No additional coordinates are present into %s', rec.parent.getMarkerName4Ch));
                                coo = rec.getPos();
                                if coo.time.isEmpty
                                    coo.setTime(rec.getPositionTime());
                                end
                            end
                        else
                            % Use coordinates from coo files
                            % If they exist
                            coo_path = rec.getPos.getCooOutPath();
                            if exist(coo_path, 'file') == 2
                                coo = Coordinates.fromCooFile(coo_path);
                            else
                                log.addWarning(sprintf('Missing coordinate file: "%s"', coo_path));
                                coo = Coordinates();
                            end
                        end
                        if not(coo.isEmpty)
                            fh = coo.showCoordinatesXYZ([], n_obs);
                            set(0, 'CurrentFigure', fh);
                            ax = subplot(3,1,1);
                            ax.Title.String{1} = [rec.parent.getMarkerName4Ch  ax.Title.String{1}(9:end)];
                            fig_name = sprintf('XYZ_at%gs_%s_%s', coo.time.getRate, rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                            fh.UserData = struct('fig_name', fig_name);
                            Core_UI.addExportMenu(fh);
                            fh_list = [fh_list; fh]; %#ok<AGROW>
                        end
                    else
                        Core.getLogger.addWarning(sprintf('%s - Plotting a single point static position is not yet supported', rec.parent.getMarkerName4Ch));
                    end
                end
            end            
        end
        
        function fh_list = showPositionSigmas(this)
            % Show Sigmas of the solutions
            %
            % SYNTAX
            %   this.showPositionSigmas();
            
            fh_list = [];
            rec = this;
            if ~isempty(rec)
                xyz = rec(1).getPosXYZ();
                if size(xyz, 1) > 1
                    rec(1).log.addMessage('Plotting ENU sigmas');
                    
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: sigma processing', f.Number); f.NumberTitle = 'off';
                    fh_list = [fh_list; f];
                    fig_name = sprintf('ENU_s0_%s_%s', rec.parent.getMarkerName4Ch, rec.time.first.toString('yyyymmdd_HHMM'));
                    f.UserData = struct('fig_name', fig_name);
                    color_order = handle(gca).ColorOrder;
                    
                    s0 = rec.quality_info.s0;
                    s0_ip = rec.quality_info.s0_ip;
                    
                    t = rec.getPositionTime().getMatlabTime;
                    
                    subplot(2,1,2);
                    plot(t, s0 * 1e2, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                    ax(2) = gca(); xlim([t(1) t(end)]); setTimeTicks(4); h = ylabel('s0 [cm]'); h.FontWeight = 'bold';
                    grid on;
                    subplot(2,1,1);
                    plot(t, s0_ip * 1e2, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                    ax(1) = gca(); xlim([t(1) t(end)]); setTimeTicks(4); h = ylabel('s0 ip [cm]'); h.FontWeight = 'bold';
                    h = title(sprintf('Receiver %s', rec(1).parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                    grid on;
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
            
            Core_Utils.addGoogleMaps('alpha', 0.95, 'MapType', 'satellite');
            title('Receiver position');
            xlabel('Longitude [deg]');
            ylabel('Latitude [deg]');
            Core_UI.beautifyFig(f);            
            Core_UI.addBeautifyMenu(f);
            f.Visible = 'on'; drawnow;            
        end
        
        
        function fh_list = showResSkyPolarScatter(sta_list, sys_c_list, type)
            % Plot residuals of the solution on polar scatter
            %
            % INPUT
            %   sys_c_list  list of satellite system to show e.g. 'GRE'
            %   type        'pr' | 'ph'  pseudo-ranges or phases
            %
            % SYNTAX 
            %   this.plotResSkyPolar(sys_c)
            
            fh_list = [];
            for r = numel(sta_list)
                rec = sta_list(r);
                if ~isempty(rec.sat.res)
                    switch nargin
                        case 1
                            fh_list = [fh_list; rec.sat.res.showResSkyPolarScatter(rec.parent.getMarkerName4Ch)]; %#ok<AGROW>
                        case 2
                            fh_list = [fh_list; rec.sat.res.showResSkyPolarScatter(rec.parent.getMarkerName4Ch, sys_c_list)]; %#ok<AGROW>
                        case 3
                            fh_list = [fh_list; rec.sat.res.showResSkyPolarScatter(rec.parent.getMarkerName4Ch, sys_c_list, type(2) == 'h')]; %#ok<AGROW>
                    end
                end
            end
        end
        
        function fh_list = showResSkyCartScatter(sta_list, sys_c_list, type)
            % Plot residuals of the solution on cartesian axes
            %
            % INPUT
            %   sys_c_list  list of satellite system to show e.g. 'GRE'
            %   type        'pr' | 'ph'  pseudo-ranges or phases
            %
            % SYNTAX 
            %   sta_list.showResSkyCartScatter(sys_c_list, type)
            
            fh_list = [];
            for r = numel(sta_list)
                rec = sta_list(r);
                if ~rec.isEmpty && ~isempty(rec.sat.res)
                    switch nargin
                        case 1
                            fh_list = [fh_list; rec.sat.res.showResSkyCartScatter(rec.parent.getMarkerName4Ch)]; %#ok<AGROW>
                        case 2
                            fh_list = [fh_list; rec.sat.res.showResSkyCartScatter(rec.parent.getMarkerName4Ch, sys_c_list)]; %#ok<AGROW>
                        case 3
                            fh_list = [fh_list; rec.sat.res.showResSkyCartScatter(rec.parent.getMarkerName4Ch, sys_c_list, type(2) == 'h')]; %#ok<AGROW>
                    end
                end
            end
        end
        
        function fh_list = showRes(sta_list, sys_c_list, type)
            % Plot residuals 
            %
            % INPUT
            %   sys_c_list  list of satellite system to show e.g. 'GRE'
            %   type        'pr' | 'ph'  pseudo-ranges or phases
            %
            % SYNTAX
            %   fh_list = sta_list.showRes()            
            fh_list = [];
            for r = numel(sta_list)
                rec = sta_list(r);
                if ~rec.isEmpty && ~isempty(rec.sat.res)
                    switch nargin
                        case 1
                            fh_list = [fh_list; rec.sat.res.showRes(rec.parent.getMarkerName4Ch)]; %#ok<AGROW>
                        case 2
                            fh_list = [fh_list; rec.sat.res.showRes(rec.parent.getMarkerName4Ch, sys_c_list)]; %#ok<AGROW>
                        case 3
                            fh_list = [fh_list; rec.sat.res.showRes(rec.parent.getMarkerName4Ch, sys_c_list, type(2) == 'h')]; %#ok<AGROW>
                    end
                end
            end            
        end       
        
        function fh_list = showResPerSat(sta_list, sys_c_list, type)
            % Plot the residuals of phase per Satellite
            %
            % INPUT
            %   sys_c_list  list of satellite system to show e.g. 'GRE'
            %   type        'pr' | 'ph'  pseudo-ranges or phases
            %
            % SYNTAX
            %   sta_list.showResPerSat(sys_c_list, is_ph)
            fh_list = [];
            for r = numel(sta_list)
                rec = sta_list(r);
                if ~rec.isEmpty && ~isempty(rec.sat.res)                    
                    switch nargin
                        case 1
                            fh_list = [fh_list; rec.sat.res.showResPerSat(rec.parent.getMarkerName4Ch)]; %#ok<AGROW>
                        case 2
                            fh_list = [fh_list; rec.sat.res.showResPerSat(rec.parent.getMarkerName4Ch, sys_c_list)]; %#ok<AGROW>
                        case 3
                            fh_list = [fh_list; rec.sat.res.showResPerSat(rec.parent.getMarkerName4Ch, sys_c_list, type(2) == 'h')]; %#ok<AGROW>
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
                fh_list = [fh_list; f];
                fig_name = sprintf('ZTD_Slant_ANI_%s_%s', this.parent.getMarkerName4Ch, this.time.first.toString('yyyymmdd_HHMM'));
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
                Core.getLogger.addWarning('ZTD and/or slants have not been computed');
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
            if isempty(this.ztd)
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
                        time_start = 2;
                        time_stop = size(sztd,1);
                    end
                    
                    %yl = (median(median(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan') + ([-6 6]) .* median(std(sztd(time_start:time_stop, :), 'omitnan'), 'omitnan'));
                    
                    az = (mod(this.sat.az(id_sync,:) + 180, 360) -180) ./ 180 * pi;
                    el = (90 - this.sat.el(id_sync,:)) ./ 180 * pi;
                    id_ko = this.sat.el < Core.getState.getCutOff | abs(sztd) > 1;
                    sztd(id_ko) = nan;
                    
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
                    
                    polarScatter(az(~id_ko), el(~id_ko), 25, abs(sztd(~id_ko)), 'filled'); hold on;
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
            
            data = rec.getU1();
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
                %[map, n_data_map] = Core_Utils.hemiGridder(az(id_ok)./180*pi, el(id_ok)./180*pi, data(id_ok), [4 1], 0.5);
                
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
                        % p_time_start = max(p_time_start,  round(sta_list(r).time.first.getRefTime(p_time_zero) * sta_list(r).time.getRate) / sta_list(r).time.getRate);
                        % p_time_stop = min(p_time_stop,  round(sta_list(r).time.last.getRefTime(p_time_zero) * sta_list(r).time.getRate) / sta_list(r).time.getRate);
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
                        lim = getFlagsLimits(~isnan(data_s(:,s)));
                    else
                        lim = getFlagsLimits(~isnan(data_s(:,s)), cs_mat(:,s));
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
                    lim = getFlagsLimits(~isnan(data_s(:,s)), cs_mat(:,s));
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
