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
        function this = Receiver_Output(parent)
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
        
        function missing_epochs = getMissingEpochs(this)
            % return a logical array of missing (code) epochs
            %
            % SYNTAX
            %   missing_epochs = this.getMissingEpochs()
            %
            missing_epochs = true(this.time.length,1);
        end
        
        function id_sync = getIdSync(this)
            id_sync = true(this.time.length, 1);
        end
        
        function [mfh, mfw] = getSlantMF(this)
            mfh = this.sat.mfh;
            mfw = this.sat.mfw;
        end
        
         function slant_td = getSlantTD(this)
            % Get the slant total delay
            % SYNTAX
            %   slant_td = this.getSlantTD();

                
                
            [mfh, mfw] = this.getSlantMF();
            n_sat = size(mfh,2);
            zwd = this.getZwd();
            apr_zhd = this.getAprZhd();
            [az, el] = this.getAzEl();
            [tgn, tge] = this.getGradient();
            res = this.getResidual();
            
            cotel = zero2nan(cotd(el));
            cosaz = zero2nan(cosd(az));
            sinaz = zero2nan(sind(az));
            slant_td = nan2zero(zero2nan(res) ...
                     + zero2nan(repmat(zwd,1,n_sat).*mfw) ...
                     + zero2nan(repmat(apr_zhd,1,n_sat).*mfh) ...
                     + zero2nan(repmat(nan2zero(tgn),1,n_sat) .* mfw .* cotel .* cosaz) ...
                     + zero2nan(repmat(nan2zero(tge),1,n_sat) .* mfw .* cotel .* sinaz));
            
        end
        
    end
    
    % ==================================================================================================================================================
    %% METHODS IMPORT / EXPORT
    % ==================================================================================================================================================
    
    methods
        function injectResult(this, rec_work)
            work_time = rec_work.getTime();
            if isempty(this.time)
                idx1 = 1;
                idx2 = 0;
                this.time = work_time;
            else
                [this.time, idx1, idx2] = this.time.injectBatch(work_time);
            end
            %%% inject data
            this.dt      = Core_Utils.injectData(this.dt, rec_work.getDt(), idx1, idx2);
            this.dt_ip   = Core_Utils.injectData(this.dt_ip, rec_work.getDtIp(), idx1, idx2);
            this.apr_zhd = Core_Utils.injectData(this.apr_zhd, rec_work.getAprZhd(), idx1, idx2);
            this.apr_zwd = Core_Utils.injectData(this.apr_zwd, rec_work.getAprZwd(), idx1, idx2);
            this.ztd     = Core_Utils.injectData(this.ztd, rec_work.getZtd(), idx1, idx2);
            this.zwd     = Core_Utils.injectData(this.zwd, rec_work.getZwd(), idx1, idx2);
            this.pwv     = Core_Utils.injectData(this.pwv, rec_work.getPwv(), idx1, idx2);
            [gn, ge]     = rec_work.getGradient();
            this.tgn     = Core_Utils.injectData(this.tgn, gn, idx1, idx2);
            this.tge     = Core_Utils.injectData(this.tge, ge, idx1, idx2);
            [p, t, h]  = rec_work.getPTH(true);
            this.pressure     = Core_Utils.injectData(this.pressure, p, idx1, idx2);
            this.temperature     = Core_Utils.injectData(this.temperature, t, idx1, idx2);
            this.humidity     = Core_Utils.injectData(this.humidity, h, idx1, idx2);
            [az, el]   = rec_work.getAzEl();
            this.sat.az     = Core_Utils.injectData(this.sat.az, az, idx1, idx2);
            this.sat.el     = Core_Utils.injectData(this.sat.el, el, idx1, idx2);
            this.sat.res    = Core_Utils.injectData(this.sat.res, rec_work.getResidual(), idx1, idx2);
            [mfh, mfw]   = rec_work.getSlantMF();
            this.sat.mfw          = Core_Utils.injectData(this.sat.mfw, mfw, idx1, idx2);
            this.sat.mfh          = Core_Utils.injectData(this.sat.mfh, mfh, idx1, idx2);
            this.sat.outlier_idx_ph    = Core_Utils.injectData(this.sat.outlier_idx_ph, rec_work.getOOutPh(), idx1, idx2);
            this.sat.cycle_slip_idx_ph = Core_Utils.injectData(this.sat.cycle_slip_idx_ph, rec_work.getOCsPh(), idx1, idx2);
            this.sat.quality           = Core_Utils.injectData(this.sat.quality, rec_work.getQuality(), idx1, idx2);
            
            %%% single results
            if isempty(this.time_pos)
                idx1 = 1;
                idx2 = 0;
                this.time_pos = rec_work.getPositionTime();
                data_len = rec_work.getPositionTime().length;
            else
                [this.time_pos, idx1, idx2] = this.time_pos.injectBatch(rec_work.getPositionTime());
                data_len = rec_work.getPositionTime().length;
            end
            
            this.xyz      = Core_Utils.injectData(this.xyz, rec_work.getPosXYZ, idx1, idx2, [data_len, 3]);
            this.enu      = Core_Utils.injectData(this.enu, rec_work.getPosENU, idx1, idx2, [data_len, 3]);
            
            [lat, lon, h_ellips, h_ortho] = rec_work.getPosGeodetic();
            
            this.lat      = Core_Utils.injectData(this.lat, lat, idx1, idx2, [data_len, 1]);
            this.lon      = Core_Utils.injectData(this.lon, lon, idx1, idx2, [data_len, 1]);
            this.h_ellips = Core_Utils.injectData(this.h_ellips, h_ellips, idx1, idx2, [data_len, 1]);
            this.h_ortho  = Core_Utils.injectData(this.h_ortho, h_ortho, idx1, idx2, [data_len, 1]);
            
            this.s0_ip   = Core_Utils.injectData(this.s0_ip, rec_work.s0_ip, idx1, idx2, [data_len, 1]);
            this.s0      = Core_Utils.injectData(this.s0, rec_work.s0, idx1, idx2, [data_len, 1]);
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
                nans = zero2nan(double(~rec.getMissingEpochs()));
                plot(t, rec.getDesync .* nans, '-k', 'LineWidth', 2);
                hold on;
                plot(t, rec.getDtPr .* nans, ':', 'LineWidth', 2);
                plot(t, rec.getDtPh .* nans, ':', 'LineWidth', 2);
                plot(t, (rec.getDtIP - rec.getDtPr) .* nans, '-', 'LineWidth', 2);
                if any(rec.getDt)
                    plot(t, rec.getDt .* nans, '-', 'LineWidth', 2);
                    plot(t, rec.getTotalDt .* nans, '-', 'LineWidth', 2);
                    legend('desync time', 'dt pre-estimated from pseudo ranges', 'dt pre-estimated from phases', 'dt correction from LS on Code', 'residual dt from carrier phases', 'total dt', 'Location', 'NorthEastOutside');
                else
                    legend('desync time', 'dt pre-estimated from pseudo ranges', 'dt pre-estimated from phases', 'dt correction from LS on Code', 'Location', 'NorthEastOutside');
                end
                xlim([t(1) t(end)]); setTimeTicks(4,'dd/mm/yyyy HH:MMPM'); h = ylabel('receiver clock error [s]'); h.FontWeight = 'bold';
                h = title(sprintf('dt - receiver %s', rec.marker_name),'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
            end
        end        
    end
end
