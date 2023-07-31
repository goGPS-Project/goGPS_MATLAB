classdef Radiometer < handle
        
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 1.0
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
    %  Written by:       Andrea Gatti, Alice Bonfiglio
    %  Contributors:     Andrea Gatti, Alice Bonfiglio
    %  A list of all the historical goGPS contributors is in CREDITS.nfo
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
    
    properties
        time          % time as GPS_Time                                    GPS_Time [1 x 1] stores n_epoch
        
        zwd           % equivalent wet zenith delay [cm]                    double   [n_epoch x 1]
        sigma_zwd     % zwd sigma error [cm]                                double   [n_epoch x 1]
        zwd_21        % equivalent wet zenith delay - 21 GHz only [cm]      double   [n_epoch x 1]
        sigma_zwd_21  % zwd sigma error [cm]                                double   [n_epoch x 1]
        zwd_saast     % zenital wet delay from Saastamoinen model           double   [n_epoch x 1]
        air_mass      % airmass of the observation in the zenith           double   [n_epoch x 1]
        el            % elevation of the wvr observation [deg]              double   [n_epoch x 1]
        az            % azimuth of the wvr observation [deg]                double   [n_epoch x 1]
        
        bright_21     % observed brightness at 21.0 GHz [K]                 double   [n_epoch x 1]
        sigma_br_21   % sigma error of brightness at 21 GHz [K]             double   [n_epoch x 1]
        bright_314    % observed brightness at 31.4 GHz [K]                 double   [n_epoch x 1]
        sigma_br_314  % sigma error of brightness at 31.4 GHz [K]           double   [n_epoch x 1]
        zbright_21    % equivalent zenith brightness at 21.0 GHz [K]        double   [n_epoch x 1]
        sigma_zbr_21  % sigma error of eq zenith brightness at 21 GHz [K]   double   [n_epoch x 1]
        zbright_314   % equivalent zenith brightness at 31.4 GHz [K]        double   [n_epoch x 1]
        sigma_zbr_314 % sigma error of eq zenith brightness at 31.4 GHz [K] double   [n_epoch x 1]
        lwc           % liquid water content in zenith direction [mm]       double   [n_epoch x 1]
        sigma_lwc     % sigma error of liquid water content [mm]            double   [n_epoch x 1]
    end
    
    properties (Access = private)
        log     % logger
    end
    
    methods
        % Creator
        function this = Radiometer(file_name)
            % Core object creator
            this.log = Logger.getInstance();
            this.reset();
            if nargin == 1
                if iscell(file_name)
                    for f = 1 : length(file_name)
                        this.log.addMessage(sprintf('Importing %s', file_name{f}));
                        if exist(file_name{f}, 'file')
                            this.appendRadiometer(file_name{f});
                        else
                            this.log.addMessage(sprintf('Error loading the last file, the file does not exists'));
                        end
                    end
                else
                    this.log.addMessage(sprintf('Importing %s', file_name));
                    if exist(file_name, 'file')
                        this.appendRadiometer(file_name);
                    else
                        this.log.addMessage(sprintf('Error loading the file, it does not exists'));
                    end
                end
            end
        end
        
    end
    
    % =========================================================================
    %  METHODS
    % =========================================================================
    
    methods % Public Access
        function reset(this)
            this.time = GPS_Time();
            this.zwd  = [];
            this.sigma_zwd  = [];
            this.zwd_21  = [];
            this.sigma_zwd_21  = [];
            this.zwd_saast  = [];
            this.air_mass  = [];
            this.el  = [];
            this.az  = [];
            
            this.bright_21  = [];
            this.sigma_br_21  = [];
            this.bright_314  = [];
            this.sigma_br_314  = [];
            this.zbright_21  = [];
            this.sigma_zbr_21  = [];
            this.zbright_314  = [];
            this.sigma_zbr_314  = [];
            this.lwc  = [];
            this.sigma_lwc  = [];
        end
        
        function importRadiometer(this, file)
            % Import after reset a tropo file
            % SYNTAX: this.importTropo(file)
            this.reset();
            this.appendRadiometer(file)
        end
        
        function appendRadiometer (this, file)
            % import and append from a tropo file
            
            % Open tropo file as a string stream
            fid = fopen(file);
            txt = fread(fid,'*char')';
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  numel(txt)
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
            lim = [lim lim(:,2) - lim(:,1)];
            
            % importing header informations
            eoh = 3; % end of header (the header of radiometer files contains 3 lines)
            
            % corrupted lines
            ko_lines = find(lim(4:end, 3) ~= median(lim(4:end,3)));
            for l = numel(ko_lines) : -1 : 1
                txt(lim(ko_lines(l)+3, 1) : lim(ko_lines(l)+3, 2) + 1) = [];
            end
            
            % extract all the epoch lines
            data = sscanf(txt(lim(4,1):end)','%4d %2d %2d %2d %2d %2d %7f %6f %7f %6f %7f %6f %6f %6f %7f %6f %7f %6f %7f %6f %7f %6f %7f %6f\n');
            data = reshape(data, 24, numel(data)/24)';
            
            % import it as a GPS_Time obj
            if this.time.length() == 0
                this.time = GPS_Time(data(:,1:6));
            else
                this.time.append6ColDate(data(:,1:6));
            end
            
            zwd = data(:,7) * 1e-2;
            sigma_zwd = data(:,8) * 1e-2;
            zwd_21 = data(:,9) * 1e-2;
            sigma_zwd_21 = data(:,10) * 1e-2;
            zwd_saast = data(:,11) * 1e-2;
            air_mass = data(:,12) * 1e-2;
            el = data(:,13);
            az = data(:,14);
            az = az + 180 * ((90-el) < 0);
            for i = 1 : length(az)
                if az(i) > 360
                    az(i) = mod(az(i),360);
                end
            end
            el = 90 - abs(90-el);
            
            bright_21  = data(:,15);
            sigma_br_21 = data(:,16);
            bright_314 = data(:,17);
            sigma_br_314 = data(:,18);
            zbright_21 = data(:,19);
            sigma_zbr_21 = data(:,20);
            zbright_314 = data(:,21);
            sigma_zbr_314 = data(:,22);
            lwc = data(:,23);
            sigma_lwc = data(:,24);
            
            % Append in obj
            this.zwd = [this.zwd; zwd];
            this.sigma_zwd = [this.sigma_zwd; sigma_zwd];
            this.zwd_21 = [this.zwd_21; zwd_21];
            this.sigma_zwd_21 = [this.sigma_zwd_21; sigma_zwd_21];
            this.zwd_saast = [this.zwd_saast; zwd_saast];
            this.air_mass = [this.air_mass; air_mass];
            this.el = [this.el; el];
            this.az = [this.az; az];
            this.bright_21 = [this.bright_21; bright_21];
            this.sigma_br_21 = [this.sigma_br_21; sigma_br_21];
            this.bright_314 = [this.bright_314; bright_314];
            this.sigma_br_314 = [this.sigma_br_314; sigma_br_314];
            
            this.zbright_21 = [this.zbright_21; zbright_21];
            this.sigma_zbr_21 = [this.sigma_zbr_21; sigma_zbr_21];
            this.zbright_314 = [this.zbright_314; zbright_314];
            this.sigma_zbr_314 = [this.sigma_zbr_314; sigma_zbr_314];
            this.lwc = [this.lwc; lwc];
            this.sigma_lwc = [this.sigma_lwc; sigma_lwc];
            
            clear data;
        end
        
        function plotAniWD(this, spline_base)
            if nargin == 1
                spline_base = 2400;
            end
            
            % get unique combinations of azimuth and elevation
            temp = this.az*1e4 + this.el * 1e1;
            [a, id_a] = unique(temp);
            azel_unique = [this.az(id_a), this.el(id_a)];
            
            ref = cell(length(a),1);
            zwd = cell(length(a),1);
            wvr = cell(length(a),1);
            
            % elapsed time in seconds from the first observation (first obs in t = 0)
            t = round(this.time.getRefTime(this.time.first.getMatlabTime)*1e5)/1e5;
            
            % prediction time
            t_pred = [0:60:this.time.last-this.time.first];
            
            zwd_wvr = nan(length(t_pred), length(a));
            
            for i = 1 : length(a)
                ref{i,1} = find(temp == a(i));
                zwd{i,1} = this.zwd(ref{i,1});
                [~,~,~, zwd_wvr(:,i)] = splinerMat(t(ref{i,1}), zwd{i,1}, spline_base, 0.1, t_pred);
                t_obs = t(ref{i,1});
                id = find(diff(t_obs) > 1000);
                
                id_ko = false(size(t_pred));
                for j = 1 : length(id)
                    id_ko = id_ko | t_pred > t_obs(id(j))-spline_base/2 & t_pred <= t_obs(id(j) + 1)+spline_base/2;
                end
                id_ko(t_pred < t_obs(1) | t_pred > t_obs(end)) = 1;
                
                zwd_wvr(id_ko,i) = nan;
                %plot(t_pred, zwd_wvr(:,i),'.'); hold on;
                t_obs = [];
            end
            %hold on; plot(t, rad.zwd,'.k','LineWidth', 0.5);
            
            %try
            xaxis_grid = [-1 : 0.1 : 1];
            yaxis_grid = [-1 : 0.1 : 1];
            [xp, yp] = meshgrid(xaxis_grid, yaxis_grid);
            
            decl = deg2rad(azel_unique(:,2))/(pi/2);
            az = -deg2rad(azel_unique(:,1)) + pi/2;
            x = cos(az) .* decl;
            y = sin(az) .* decl;
            
            max_delta = max(abs(yaxis_grid(end)-yaxis_grid(1)),abs(xaxis_grid(end)-xaxis_grid(1)))/2;
            fun = @(dist) exp(-((dist/max_delta)*1e5/5e3).^2);
            
            for i = 1 : 5: size(zwd_wvr,1)
                id_ok = not(isnan(zero2nan(zwd_wvr(i,:))));
                delay = funInterp2(xp(:), yp(:), x(id_ok), y(id_ok), zwd_wvr(i,id_ok)', fun);
                if i == 1
                    figure;
                    subplot(3,1,3); plot(t_pred,zwd_wvr(:,:),'LineWidth', 1.5); hold on; plot(t, this.zwd,'.k','MarkerSize',4); yl = ylim;
                    hl = line('XData', t_pred(1) * [1 1],'YData', yl); ylim(yl); xlim([t_pred(1) t_pred(end)]);
                    subplot(3,1,1:2);
                    hm = imagesc(xaxis_grid, yaxis_grid, reshape(delay(:), numel(yaxis_grid), numel(xaxis_grid)));
                    hold on;
                    fh = gcf; fh.Children(2).YDir = 'normal'; ax = gca; ax.YDir = 'normal';
                    hm.AlphaData = 0.5;
                    h = polarScatter(deg2rad(azel_unique(:,1)), deg2rad(azel_unique(:,2)), 120, ones(size(azel_unique,1),1), 'filled'); hold on;
                    x_data = h.XData;
                    %focus of the polarscatter
                    caxis([min(this.zwd) max(this.zwd)]);
                    colormap(jet); colorbar();
                    xlim(xlim());
                    ylim(ylim());
                    h.CData = serialize(zwd_wvr(i,:));
                    h.XData = x_data + 1e10 .* isnan(serialize(zwd_wvr(i,:)))'; %don't show nan values
                    hold on;
                    colormap(jet); colorbar;
                    hl.XData = t_pred(i) * [1 1];
                    colormap(jet); colorbar();
                else
                    h.CData = serialize(zwd_wvr(i,:));
                    h.XData = x_data + 1e10 .* isnan(serialize(zwd_wvr(i,:)))'; %don't show nan values
                    hold on;
                    hl.XData = t_pred(i) * [1 1];
                    hm.CData = reshape(delay(:), numel(yaxis_grid), numel(xaxis_grid));
                end
                drawnow
            end
        end
    end
    
    methods (Static)
        function wvr = loadBatch(file_name, run_start, run_stop)
            % Load all the tropo files of a session
            %
            % SYNTAX:  
            %   tropo = Tropo.loadBatch(file_name, run_start, run_stop)
            %
            % INPUT:
            %   file_name     it should include the key ${RUN} that will be substituted with a 3 digits number containing the run, from run_start to run_stop
            %   run_start     number of the first run to load
            %   run_stop      number of the last run to load
            %
            % OUTPUT:
            %   tropo         troposphere object containing all the data loaded from the files
            %
            GPS_RUN = '${RUN}';
            file_name_list = {};
            r = 0;
            for run = run_start : run_stop
                cur_file_name = strrep(file_name, GPS_RUN, sprintf('%02d', run));
                if exist(cur_file_name, 'file') == 2
                    r = r + 1;
                    file_name_list{r} = cur_file_name;
                end
            end
            
            if numel(file_name_list) > 0
                wvr = Radiometer(file_name_list);
            else
                wvr = Radiometer();
            end
        end
    end
end
