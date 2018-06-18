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
classdef Receiver_Commons < handle
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
        enu            % position of the receiver (ENU local)
        
        lat            % ellipsoidal latitude
        lon            % ellipsoidal longitude
        h_ellips       % ellipsoidal height
        h_ortho        % orthometric height
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
    end
    
    % ==================================================================================================================================================
    %% PROPERTIES QUALITY INDEXES
    % ==================================================================================================================================================
    
    properties
        s0_ip
        s0
        hdop
        vdop
        tdop
        a_fix
        s_rate
    end
    
    % ==================================================================================================================================================
    %% PROPERTIES USEFUL HANDLES
    % ==================================================================================================================================================
    
    properties (SetAccess = protected, GetAccess = public)
        cc = Constellation_Collector('G');      % local cc
        w_bar                                  % handle to waitbar
        state                                  % local handle of state;
        log                                    % handle to log
        rf                                     % handle to reference farme
    end
    
    
    
    
    % ==================================================================================================================================================
    %% METHODS INIT - CLEAN - RESET - REM -IMPORT
    % ==================================================================================================================================================
    
    methods
        
        function reset(this)
            this.time = GPS_Time();
            this.enu = [];
            this.lat = [];
            this.lon = [];
            
            this.h_ellips = [];
            this.h_ortho = [];
            
            this.s0_ip =  [];
            this.s0 =  [];
            this.hdop =  [];
            this.vdop =  [];
            this.tdop =  [];
            
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
            
            this.sat = struct( ...
                'avail_index',      [], ...    % boolean [n_epoch x n_sat] availability of satellites
                'outlier_idx_ph',   [], ...    % logical index of outliers
                'cycle_slip_idx_ph',[], ...    % logical index of cycle slips
                'err_tropo',        [], ...    % double  [n_epoch x n_sat] tropo error
                'err_iono',         [], ...    % double  [n_epoch x n_sat] iono error
                'solid_earth_corr', [], ...    % double  [n_epoch x n_sat] solid earth corrections
                'dtS',              [], ...    % double  [n_epoch x n_sat] staellite clok error at trasmission time
                'rel_clk_corr',     [], ...    % double  [n_epoch x n_sat] relativistic correction at trasmission time
                'tot',              [], ...    % double  [n_epoch x n_sat] time of travel
                'az',               [], ...    % double  [n_epoch x n_sat] azimuth
                'el',               [], ...    % double  [n_epoch x n_sat] elevation
                'cs',               [], ...    % Core_Sky
                'XS_tx',            [], ...    % compute Satellite postion a t transmission time
                'crx',              [], ...    % bad epochs based on crx file
                'res',              [], ...    % residual per staellite
                'slant_td',         []  ...    % slant total delay (except ionosphere delay)
                );
        end
    end
    % ==================================================================================================================================================
    %% METHODS GETTER - TIME
    % ==================================================================================================================================================
    
    methods
        function toStringPos(this)
            % Display on screen information about the receiver position
            % SYNTAX this.toStringPos();
            for r = 1 : numel(this)
                if ~this(r).isEmpty && ~isempty(this(r).xyz)
                    [lat, lon, h_ellips, h_ortho] = this(r).getMedianPosGeodetic_mr();
                    this(r).log.addMarkedMessage(sprintf('Receiver (%02d) %s   %11.7f  %11.7f    %12.7f m (ellipsoidal) - %12.7f (orthometric)', r, this(r).parent.marker_name, lat, lon, h_ellips, h_ortho));
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
            id = round(this(1).time.length()/2);
            time = this(1).time.getEpoch(id);
        end
        
        function [rate] = getRate(this)
            % SYNTAX
            %   rate = this.getRate();
            rate = this.getTime.getRate;
        end
        
        % position
        
        function dt = getTotalDt(this)
            dt = this.getDt + this.getDtIP;
        end
        
        function xyz = getPosXYZ_mr(this)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   xyz     geocentric coordinates
            %
            % SYNTAX
            %   xyz = this.getPosXYZ_mr()
            xyz = this.getPosXYZ();
            xyz = permute(reshape(xyz', 3, n_sss, n_rec), [2 1 3]);
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
                xyz = [xyz; this(r).xyz]; %#ok<AGROW>
            end
        end
        
        function [lat, lon, h_ellips, h_ortho] = getPosGeodetic(this)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   lat, lon, h_ellips, h_ortho     geodetic coordinates
            %
            % SYNTAX
            %   [lat, lon, h_ellips, h_ortho] = this.getPosGeodetic()
            [lat, lon, h_ellips] = cart2geod(this.getPosXYZ);
            if nargout == 4
                gs = Global_Configuration.getInstance;
                gs.initGeoid();
                ondu = getOrthometricCorr(lat, lon, gs.getRefGeoid());
                h_ortho = h_ellips - ondu;
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
            xyz = this.getPosXYZ();
            [enu(:,1), enu(:,2), enu(:,3)] = cart2plan(zero2nan(xyz(:,1)), zero2nan(xyz(:,2)), zero2nan(xyz(:,3)));
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
            for i = 1:numel(this)
                xyz = this(i).getMedianPosXYZ;
                [EAST, NORTH, h, utm_zone] = cart2plan(xyz(:,1), xyz(:,2), xyz(:,3));
                utm = [EAST, NORTH, h];
                utm_zone = utm_zone;
            end
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
        
        function [lat, lon, h_ellips, h_ortho] = getMedianPosGeodetic(this)
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
            xyz = this.getPosXYZ();
            xyz = median(xyz, 1);
            if ~isempty(this(1))
                [lat, lon, h_ellips] = cart2geod(xyz);
                if nargout == 4
                    gs = Global_Configuration.getInstance;
                    gs.initGeoid();
                    ondu = getOrthometricCorr(lat, lon, gs.getRefGeoid());
                    h_ortho = h_ellips - ondu;
                end
                lat = lat / pi * 180;
                lon = lon / pi * 180;
            else
                lat = [];
                lon = [];
                h_ellips = [];
                h_ortho = [];
            end
        end
        
        function ztd = getZtd(this)
            % get ztd
            %
            % SYNTAX
            %   ztd = this.getZtd()
            ztd = this.ztd(this.getIdSync);
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
            else
                this(1).log.addWarning('ZTD and slants have not been computed');
            end
        end
        
        function apr_zhd = getAprZhd(this)
            % get a-priori ZHD
            %
            % SYNTAX
            %   zhd = this.getAprZhd()
            apr_zhd = this.apr_zhd(this.getIdSync);
        end
        
        function zwd = getZwd(this)
            % get zwd
            %
            % SYNTAX
            %   zwd = this.getZwd()
            zwd = this.zwd(this.getIdSync);
            if isempty(zwd) || all(isnan(zero2nan(zwd)))
                zwd = this.apr_zwd(this.getIdSync);
            end
        end
        
        function pwv = getPwv(this)
            % get pwv
            %
            % SYNTAX
            %   pwv = this.getPwv()
            pwv = this.pwv(this.getIdSync);
        end
        
        function [gn ,ge, time] = getGradient(this)
            % SYNTAX
            % [gn ,ge, time] = getGradient(this)
            if isempty(this.tgn)
                gn = nan(length(this.getIdSync),1);
            else
                gn = this.tgn(this.getIdSync);
            end
            if isempty(this.tgn)
                ge = nan(length(this.getIdSync),1);
            else
                ge = this.tge(this.getIdSync);
            end
            time = this.time.getSubSet(this.getIdSync);
            
        end
        
        function [apr_zwd, time] = getAprZwd(this)
            % SYNTAX
            %  [apr_zwd, time] = this.getAprZwd()
            
            apr_zwd = this.apr_zwd(this.getIdSync);
            time = this.time.getSubSet(this.getIdSync);
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
            if nargin < 2
                go_id = 1 : this.cc.getMaxNumSat;
            end
            az = this.sat.az(this.getIdSync, go_id);
        end
        
        function [el] = getEl(this, go_id)
            % Get the azimuth and elevation (on valid id_sync)
            %
            % SYNTAX
            %   el = this.getEl();
            if nargin < 2
                go_id = 1 : this.cc.getMaxNumSat;
            end
            el = this.sat.el(this.getIdSync, go_id);
        end
        
        function res = getResidual(this)
            % get residual
            %
            % SYNTAX
            %   res = this.getResidual()
            res = this.sat.res(this.getIdSync,:);
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
        function exportTropoSINEX(this)
            for r = 1 : numel(this)
                try
                    rec = this(r);
                    if ~isempty(rec.getZtd)
                        [year, doy] = rec.time.first.getDOY();
                        yy = num2str(year);
                        yy = yy(3:4);
                        sess_str = '0'; %think how to get the right one from sss_id_list
                        fname = sprintf('%s',[rec.state.getOutDir() filesep rec.parent.marker_name sprintf('%03d', doy) sess_str '.' yy 'zpd']);
                        snx_wrt = SINEX_Writer(fname);
                        snx_wrt.writeTroSinexHeader( rec.time.first, rec.time.getSubSet(rec.time.length), rec.parent.marker_name)
                        snx_wrt.writeFileReference()
                        snx_wrt.writeAcknoledgments()
                        smpl_tropo = median(diff(rec.getIdSync)) * rec.time.getRate;
                        val_flags = {'TROTOT','TGNTOT','TGETOT'};
                        snx_wrt.writeTropoDescription(rec.state.cut_off, rec.time.getRate, smpl_tropo, snx_wrt.SINEX_MAPPING_FLAGS{this.state.mapping_function},val_flags, false(3,1))
                        snx_wrt.writeSTACoo( rec.parent.marker_name, rec.xyz(1,1), rec.xyz(1,2), rec.xyz(1,3), 'UNDEF', 'GRD'); % The reference frame depends on the used orbit so it is generraly labled undefined a more intelligent strategy could be implemented
                        snx_wrt.writeTropoSolutionSt()
                        snx_wrt.writeTropoSolutionStation(  rec.parent.marker_name, rec.time.getSubSet(rec.getIdSync), [rec.ztd(rec.getIdSync,:) rec.tgn(rec.getIdSync,:) rec.tge(rec.getIdSync,:)]*1000, [], {'TROTOT','TGNTOT','TGETOT'})
                        snx_wrt.writeTropoSolutionEnd()
                        snx_wrt.writeTroSinexEnd();
                        snx_wrt.close()
                        rec(1).log.addStatusOk(sprintf('Tropo saved into: %s', fname));
                    end
                catch ex
                    rec(1).log.addError(sprintf('saving Tropo in sinex format failed: %s', ex.message));
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
            %  - time_utc in matlab format
            %
            % SYNTAX
            %   this.exportTropoMat
            
            for t = 1 : numel(this)
                try
                    this(t).updateCoordinates;
                    time = this(t).getTime();
                    [year, doy] = this(t).getCentralTime.getDOY();
                    time.toUtc();
                    
                    lat = this(t).lat; %#ok<NASGU>
                    lon = this(t).lon; %#ok<NASGU>
                    h_ellips = this(t).h_ellips; %#ok<NASGU>
                    h_ortho = this(t).h_ortho; %#ok<NASGU>
                    ztd = this(t).getZtd(); %#ok<NASGU>
                    utc_time = time.getMatlabTime; %#ok<NASGU>
                    
                    fname = sprintf('%s',[this(t).state.getOutDir() filesep this(t).parent.marker_name sprintf('%04d%03d',year, doy) '.mat']);
                    save(fname, 'lat', 'lon', 'h_ellips', 'h_ortho', 'ztd', 'utc_time','-v6');
                    
                    this(1).log.addStatusOk(sprintf('Tropo saved into: %s', fname));
                catch ex
                    this(1).log.addError(sprintf('saving Tropo in matlab format failed: %s', ex.message));
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
                fid = fopen(fname,'w');
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
        
        function showPositionENU(this, one_plot)
            % Plot East North Up coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionENU();
            if nargin == 1
                one_plot = false;
            end
            
            for r = 1 : numel(this)
                rec = this(r);
                if ~isempty(rec)
                    xyz = rec.getPosXYZ();
                    if size(rec, 1) > 1 || size(xyz, 1) > 1
                        rec(1).log.addMessage('Plotting positions');
                        xyz0 = rec.getMedianPosXYZ();
                        enu0 = [];
                        [enu0(:,1), enu0(:,2), enu0(:,3)] = cart2plan(xyz0(:,1), xyz0(:,2), xyz0(:,3));
                        
                        f = figure; f.Name = sprintf('%03d: PosENU', f.Number); f.NumberTitle = 'off';
                        color_order = handle(gca).ColorOrder;
                        
                        xyz = rec.getPosXYZ();
                        xyz0 = rec.getMedianPosXYZ();
                        enu = [];
                        [enu(:,1), enu(:,2), enu(:,3)] = cart2plan(zero2nan(xyz(:,1)), zero2nan(xyz(:,2)), zero2nan(xyz(:,3)));
                        
                        t = rec.getPositionTime().getMatlabTime();
                        
                        [enu0(:,1), enu0(:,2), enu0(:,3)] = cart2plan(xyz0(:,1), xyz0(:,2), xyz0(:,3));
                        [enu(:,1), enu(:,2), enu(:,3)] = cart2plan(zero2nan(xyz(:,1)), zero2nan(xyz(:,2)), zero2nan(xyz(:,3)));
                        
                        if ~one_plot, subplot(3,1,1); end
                        plot(t, (1e3 * (enu(:,1) - enu0(1))), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                        ax(3) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('East [mm]'); h.FontWeight = 'bold';
                        grid on;
                        h = title(sprintf('Receiver %s', rec(1).parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                        if ~one_plot, subplot(3,1,2); end
                        plot(t, (1e3 * (enu(:,2) - enu0(2))), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                        ax(2) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('North [mm]'); h.FontWeight = 'bold';
                        grid on;
                        if ~one_plot, subplot(3,1,3); end
                        plot(t, (1e3 * (enu(:,3) - enu0(3))), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                        ax(1) = gca();
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('Up [mm]'); h.FontWeight = 'bold';
                        grid on;
                        if one_plot
                            h = ylabel('ENU [m]'); h.FontWeight = 'bold';
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
        
        function showPositionXYZ(this, one_plot)
            % Plot X Y Z coordinates of the receiver (as estimated by initDynamicPositioning
            % SYNTAX this.plotPositionXYZ();
            if nargin == 1
                one_plot = false;
            end
            
            rec = this;
            if ~isempty(rec)
                xyz = rec(:,1).getPosXYZ();
                if size(rec, 1) > 1 || size(xyz, 1) > 1
                    rec(1).log.addMessage('Plotting positions');
                    
                    f = figure; f.Name = sprintf('%03d: PosXYZ', f.Number); f.NumberTitle = 'off';
                    color_order = handle(gca).ColorOrder;
                    
                    xyz = rec(:).getPosXYZ();
                    xyz0 = rec(:).getMedianPosXYZ();
                    
                    t = rec.getPositionTime().getMatlabTime;                    
                    
                    x = 1e3 * bsxfun(@minus, zero2nan(xyz(:,1)), xyz0(1));
                    y = 1e3 * bsxfun(@minus, zero2nan(xyz(:,2)), xyz0(2));
                    z = 1e3 * bsxfun(@minus, zero2nan(xyz(:,3)), xyz0(3));
                    
                    if ~one_plot, subplot(3,1,1); end
                    plot(t, x, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                    ax(3) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('X [mm]'); h.FontWeight = 'bold';
                    grid on;
                    h = title(sprintf('Receiver %s', rec(1).parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                    if ~one_plot, subplot(3,1,2); end
                    plot(t, y, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                    ax(2) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('Y [mm]'); h.FontWeight = 'bold';
                    grid on;
                    if ~one_plot, subplot(3,1,3); end
                    plot(t, z, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                    ax(1) = gca(); xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('Z [mm]'); h.FontWeight = 'bold';
                    grid on;
                    if one_plot
                        h = ylabel('XYZ [m]'); h.FontWeight = 'bold';
                    end
                    linkaxes(ax, 'x');
                else
                    rec.log.addMessage('Plotting a single point static position is not yet supported');
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
            [lat, lon] = cart2geod(this.getMedianPosXYZ());
            
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
                name = upper(this(r).parent.getMarkerName());
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
        end
                       
        
        function showResSky_p(this, sys_c_list)
            % Plot residuals of the solution on polar scatter
            % SYNTAX this.plotResSkyPolar(sys_c)
            
            if isempty(this.sat.res)
                this.log.addWarning('Residuals have not been computed');
            else
                if nargin == 1
                    sys_c_list = unique(this.cc.system);
                end
                
                for sys_c = sys_c_list
                    s = this.cc.getGoIds(sys_c);%this.go_id(this.system == sys_c);
                    res = abs(this.sat.res(:, s));
                    
                    f = figure; f.Name = sprintf('%03d: Res P %s', f.Number, this.cc.getSysName(sys_c)); f.NumberTitle = 'off';
                    id_ok = (res~=0);
                    az = this.sat.az(:, s);
                    el = this.sat.el(:, s);
                    polarScatter(serialize(az(id_ok))/180*pi,serialize(90-el(id_ok))/180*pi, 45, serialize(res(id_ok)), 'filled');
                    caxis([min(abs(this.sat.res(:))) min(20, min(6*std(zero2nan(this.sat.res(:)),'omitnan'), max(abs(zero2nan(this.sat.res(:))))))]);
                    colormap(flipud(hot)); f.Color = [.95 .95 .95]; colorbar();
                    h = title(sprintf('Satellites residuals [m] - receiver %s - %s', this.parent.marker_name, this.cc.getSysExtName(sys_c)),'interpreter', 'none');  h.FontWeight = 'bold'; h.Units = 'pixels'; h.Position(2) = h.Position(2) + 20; h.Units = 'data';
                end
            end
        end
        
        function showResSky_c(this, sys_c_list)
            % Plot residuals of the solution on cartesian axes
            % SYNTAX this.plotResSkyCart()
            if isempty(this.sat.res)
                this.log.addWarning('Residuals have not been computed');
            else
                if nargin == 1
                    sys_c_list = unique(this.cc.system);
                end
                
                for sys_c = sys_c_list
                    s  = this.cc.getGoIds(sys_c);%unique(this.go_id(this.system == sys_c));
                    res = abs(this.sat.res(:, s));
                    
                    f = figure; f.Name = sprintf('%03d: Res C %s', f.Number, this.cc.getSysName(sys_c)); f.NumberTitle = 'off';
                    %this.updateAzimuthElevation()
                    id_ok = (res~=0);
                    az = this.sat.az(:, s);
                    el = this.sat.el(:, s);
                    scatter(serialize(az(id_ok)),serialize(el(id_ok)), 45, serialize(res(id_ok)), 'filled');
                    caxis([min(abs(this.sat.res(:))) min(20, min(6*std(zero2nan(this.sat.res(:)),'omitnan'), max(abs(zero2nan(this.sat.res(:))))))]);
                    colormap(flipud(hot)); f.Color = [.95 .95 .95]; colorbar(); ax = gca; ax.Color = 'none';
                    h = title(sprintf('Satellites residuals [m] - receiver %s - %s', this.parent.marker_name, this.cc.getSysExtName(sys_c)),'interpreter', 'none');  h.FontWeight = 'bold'; h.Units = 'pixels'; h.Position(2) = h.Position(2) + 20; h.Units = 'data';
                    hl = xlabel('Azimuth [deg]'); hl.FontWeight = 'bold';
                    hl = ylabel('Elevation [deg]'); hl.FontWeight = 'bold';
                end
            end
        end
        
        function showAniZtdSlant(this, time_start, time_stop, show_map, write_video)
            if isempty(this.ztd) || ~any(this.sat.slant_td(:))
                this.log.addWarning('ZTD and slants have not been computed');
            else
                f = figure; f.Name = sprintf('%03d: AniZtd', f.Number); f.NumberTitle = 'off';
                
                sztd = this.getSlantZTD(this.parent.slant_filter_win);
                
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
                
                if isempty(this.id_sync(:))
                    this.id_sync(:, 1) = 1 : this.time.length();
                end
                id_ok = this.id_sync(this.id_sync(:) > time_start & this.id_sync(:) < time_stop);
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
                setTimeTicks(4,'dd/mm/yy HH:MM');
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
                    az = (mod(this.sat.az(this.id_sync(i, 1),:) + 180, 360) -180) ./ 180 * pi; az(isnan(az) | isnan(sztd(i,:))) = 1e10;
                    el = (90 - this.sat.el(this.id_sync(i, 1),:)) ./ 180 * pi; el(isnan(el) | isnan(sztd(i,:))) = 1e10;
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
        
        function showAniZwdSlant(this, time_start, time_stop, show_map)
            if isempty(this.zwd) || ~any(this.sat.slant_td(:))
                this.log.addWarning('ZWD and slants have not been computed');
            else
                f = figure; f.Name = sprintf('%03d: AniZwd', f.Number); f.NumberTitle = 'off';
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
                setTimeTicks(4,'dd/mm/yyyy HH:MMPM');
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
        
        function showZtdSlant(this, time_start, time_stop)
            %if isempty(this(1).ztd) || ~any(this(1).sat.slant_td(:))
            %    this(1).log.addWarning('ZTD and/or slants have not been computed');
            %else
            rec = this;
            if isempty(rec)
                this(1).log.addWarning('ZTD and/or slants have not been computed');
            else
                f = figure; f.Name = sprintf('%03d: Ztd Slant %s', f.Number, rec(1).cc.sys_c); f.NumberTitle = 'off';
                t = rec(:).getTime.getMatlabTime;
                
                sztd = rec(:).getSlantZTD(rec(1).parent.slant_filter_win);
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
                setTimeTicks(4,'dd/mm/yyyy HH:MMPM');
                h = ylabel('ZTD [m]'); h.FontWeight = 'bold';
                grid on;
                h = title(sprintf('Receiver %s ZTD', rec(1).parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                drawnow;
            end
            
            
        end
        
        function showTropoPar(this, par_name, new_fig)
            % one function to rule them all
            rec_ok = false(size(this,2), 1);
            for r = 1 : size(this, 2)
                switch lower(par_name)
                    case 'ztd'
                        rec_ok(r) = any(~isnan(this(:,r).getZtd));
                    case 'zwd'
                        rec_ok(r) = any(~isnan(this(:,r).getZwd));
                    case 'pwv'
                        rec_ok(r) = any(~isnan(this(:,r).getPwv));
                    case 'zhd'
                        rec_ok(r) = any(~isnan(this(:,r).getAprZhd));
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
                        tropo = rec_list.getZtd();
                    case 'zwd'
                        tropo = rec_list.getZwd();
                    case 'pwv'
                        tropo = rec_list.getPwv();
                    case 'zhd'
                        tropo = rec_list.getAprZhd();
                end
                
                if ~iscell(tropo)
                    tropo = {tropo};
                end
                if isempty(tropo)
                    rec_list(1).log.addWarning([par_name ' and slants have not been computed']);
                else
                    if new_fig
                        f = figure; f.Name = sprintf('%03d: %s %s', f.Number, par_name, rec_list(1).cc.sys_c); f.NumberTitle = 'off';
                        old_legend = {};
                    else
                        l = legend;
                        old_legend = get(l,'String');
                    end
                    for r = 1 : size(rec_list, 2)
                        rec = rec_list(~rec_list(:,r).isEmpty, r);
                        if ~isempty(rec)
                            t = rec.getTime();
                            switch lower(par_name)
                                case 'ztd'
                                    tropo = rec.getZtd();
                                case 'zwd'
                                    tropo = rec.getZwd();
                                case 'pwv'
                                    tropo= rec.getPwv();
                                case 'zhd'
                                    tropo = rec.getAprZhd();
                            end
                            if new_fig
                                plot(t.getMatlabTime(), zero2nan(tropo'), '.', 'LineWidth', 4, 'Color', Core_UI.getColor(r, size(rec_list, 2))); hold on;
                            else
                                plot(t.getMatlabTime(), zero2nan(tropo'), '.', 'LineWidth', 4); hold on;
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
                    f = figure; f.Name = sprintf('%03d: Median %s %s', f.Number, par_name, rec_list(1).cc.sys_c); f.NumberTitle = 'off';
                    old_legend = {};
                else
                    l = legend;
                    old_legend = get(l,'String');
                end
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
            if isempty(this.ztd)
                this.log.addWarning('ZTD and slants have not been computed');
            else
                
                id_sync = this.getIdSync();
                
                t = this.time.getEpoch(id_sync).getMatlabTime;
                
                sztd = this.getSlantZTD(this.parent.slant_filter_win);
                sztd = bsxfun(@minus, sztd, this.ztd(id_sync));
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
                
                f = figure; f.Name = sprintf('%03d: Slant res', f.Number); f.NumberTitle = 'off';
                polarScatter(az(:), el(:), 25, abs(sztd(:)), 'filled'); hold on;
                caxis(minMax(abs(sztd))); colormap(flipud(hot)); f.Color = [.95 .95 .95]; colorbar();
                h = title(sprintf('Receiver %s ZTD - Slant difference', this.parent.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
            end
        end
        
        function plotResidual(this)
            figure
            plot(zero2nan(this.sat.res),'.');
        end
    end
end
