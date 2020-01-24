classdef Tropo_Mat_Reader
%   CLASS Tropo_Mat_Reader
% =========================================================================
%
% DESCRIPTION
%   Class to load and visualize tropo_mat out file of goGPS
%
% EXAMPLE
%   tm = Tropo_Mat_Reader(marker_name);
%
% FOR A LIST OF CONSTANTS and METHODS use doc File_Rinex

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

    properties        
        base_dir
        dir_list
        marker
        data_set
    end
    
    methods (Static)
        % Creator
        function this = Tropo_Mat_Reader(file_list_base_dir, marker)
            % Create tropo mat object
            %
            % SINTAX
            %   this = Tropo_Mat_Reader(file_list_base_dir, marker);
            if nargin == 2
                % Tropo_Mat_Receiver object creator
                if ~iscell(file_list_base_dir)
                    base_dir = {file_list_base_dir};
                else
                    base_dir = file_list_base_dir;
                end
                file_list = {};
                for d = 1 : numel(base_dir)
                    this.base_dir = base_dir{d};
                    this.marker = marker;
                    
                    file_list = [file_list sort(File_Name_Processor.findFiles(base_dir{d}, marker))];
                end
                file_list = sort(file_list);
            else
                file_list = file_list_base_dir;
            end
            first_epoch = [];
            i = 0;
            for f = 1 : numel(file_list)
                [~, ~, ext] = fileparts(file_list{f});
                if ~strcmp(ext, '.mat')
                    fprintf('Skipping %5d/%5d "%s"\n', f, numel(file_list), file_list{f})
                else
                    i = i + 1;
                    fprintf('Loading %5d/%5d "%s"\n', f, numel(file_list), file_list{f})
                    tmp = load(file_list{f});
                    first_epoch(i) = tmp.utc_time(1);
                    tmp.coo = Coordinates.fromGeodetic(tmp.lat/180*pi, tmp.lon/180*pi, tmp.h_ellips);
                    tmp.time = GPS_Time(tmp.utc_time, [], false);
                    data_set(i) = tmp;
                end
            end
            % sort by time
            [~, id] = sort(first_epoch);
            if ~isempty(id)
                this.data_set = data_set(id);
            end
        end
    end

    methods
        function showBaselineENU(this, ref)
            % Show ENU baseline variations
            %
            % SINTAX
            %   this.showBaselineENU(tropo_mat_obj_ref);
            if ~isempty(this.data_set)
                f = figure('Visible', 'off'); f.Name = sprintf('%03d: BSL %s - %s', f.Number, ref.marker, this.marker); f.NumberTitle = 'off';
                        
                colormap(flipud(gat2));
                for i = 1 : numel(this.data_set)
                    enu(i,:) = this.data_set(i).coo.getENU();
                    tmp = this.data_set(i).time.getCentralTime.getCopy(); tmp.toUtc();
                    time(i) = round(tmp.getMatlabTime() * 24) / 24;                    
                end
                for i = 1 : numel(ref.data_set)
                    enu_ref(i,:) = ref.data_set(i).coo.getENU();
                    tmp = ref.data_set(i).time.getCentralTime.getCopy(); tmp.toUtc();
                    time_ref(i) = round(tmp.getMatlabTime() * 24) / 24;                   
                end
                [time, ida, idb] = intersect(time, time_ref);
                      
                enu = enu_ref(idb,:) - enu(ida,:);

                enu = enu - repmat(median(enu, 'omitnan'), size(enu,1), 1); % enu wrt first epoch

                for a = 1: 3
                    subplot(3,1,a);
                    plot(time, enu(:,a) * 1e3, '.', 'MarkerSize', 10, 'Color', Core_UI.getColor(a,3)); hold on;
                end
                ax = [];
                label = {'East [mm]', 'North [mm]', 'Up [mm]'};
                for a = 3:-1:1
                    ax(a) = subplot(3,1,a);
                    h = title(sprintf('std %.2f [cm]',sqrt(var(enu(:,a)*1e2))),'interpreter', 'none'); h.FontWeight = 'bold';
                    setTimeTicks(5,'mmm-dd HH:MM');
                    h = ylabel(label{a}); h.FontWeight = 'bold';
                    grid on;
                    xlim(minMax(time));
                    yl = ylim();
                    yl(1) = min(-20, yl(1));
                    yl(2) = max(20, yl(2));
                    ylim(yl);
                end
                linkaxes(ax, 'x');
                Core_UI.beautifyFig(f, 'light');
                Core_UI.addBeautifyMenu(f);
                f.Visible = 'on';
            end
        end
        
        function showPositionENU(this, color)
            % Show ENU position variations
            %
            % SINTAX
            %   this.showPositionENU(color);
            if ~isempty(this.data_set)
                if nargin < 2
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: ENU %s ', f.Number, this.marker); f.NumberTitle = 'off';
                end
                
                for i = 1 : numel(this.data_set)
                    enu(i,:) = this.data_set(i).coo.getENU();
                    tmp = this.data_set(i).time.getCentralTime.getCopy(); tmp.toUtc();
                    time(i) = tmp.getMatlabTime();                    
                end
                enu = enu - repmat(median(enu, 'omitnan'), size(enu,1), 1); % enu wrt first epoch                
                for a = 1: 3
                    subplot(3,1,a);
                    if nargin < 2
                        plot(time, enu(:,a) * 1e3, '.', 'MarkerSize', 10, 'Color', Core_UI.getColor(a,3)); hold on;
                    else
                        plot(time, enu(:,a) * 1e3, 'o', 'Color', color); hold on;
                    end
                end
                ax = [];
                label = {'East [mm]', 'North [mm]', 'Up [mm]'};
                for a = 3:-1:1
                    ax(a) = subplot(3,1,a);
                    cur_ax = get(ax(a));
                    if nargin < 2
                        h = title(sprintf('std %.2f [cm]',sqrt(var(enu(:,a)*1e2))),'interpreter', 'none'); h.FontWeight = 'bold';
                    else                        
                        h = title(sprintf('%s vs std %.2f [cm]', cur_ax.Title.String, sqrt(var(enu(:,a)*1e2))),'interpreter', 'none'); h.FontWeight = 'bold';
                    end
                    setTimeTicks(5,'mmm-dd HH:MM');
                    h = ylabel(label{a}); h.FontWeight = 'bold';
                    grid on;
                    xlim(minMax(time));
                end
                linkaxes(ax, 'x');
                Core_UI.beautifyFig(f, 'light');
                Core_UI.addBeautifyMenu(f);
                f.Visible = 'on';
            end
        end
        
        function [time_utc, tropo_opt, tropo_std] = showAllTropoPar(this, tropo_par, color)
            % Show Tropo parameter
            %
            % INPUT
            %   tropo_par   i.e. 'ztd', 'zwd', ...
            %   color       color of the lines / if not present color with temporal distance
            %            
            % SINTAX
            %   this.showAllTropoPar(tropo_par, color)
            if ~isempty(this.data_set)
                [time_utc, tropo_opt, tropo_std] = this.getBestTropo('ztd');
                if nargin < 3
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: %s %s ', f.Number, upper(tropo_par), this.marker); f.NumberTitle = 'off';
                else
                    f = gcf;
                end
                colormap(flipud(gat2));
                min_time = inf;
                max_time = -inf;
                for i = 1 : numel(this.data_set)
                    min_time = min(min_time, min(this.data_set(i).utc_time));
                    max_time = max(max_time, max(this.data_set(i).utc_time));
                    if nargin < 3 || isempty(color)
                        scatter(this.data_set(i).utc_time, this.data_set(i).(lower(tropo_par)) .* 1e2, 10, (this.data_set(i).utc_time(end) - this.data_set(i).utc_time) .* 24, 'filled'); hold on
                    else
                        plot(this.data_set(i).utc_time, this.data_set(i).(lower(tropo_par)) .* 1e2, 'color', color); hold on
                    end
                end
                plot(time_utc(tropo_std .* 1e2 <= 0.6), tropo_opt(tropo_std .* 1e2 <= 0.6) .* 1e2, '-k', 'LineWidth', 2); % plot only "good" prediction values
                if nargin < 3 || isempty(color)
                    caxis = [0 12];
                    colorbar('Location', 'south');
                end
                
                axis tight
                grid minor;
                xlim([min_time, max_time]);
                setTimeTicks();
                ylim(ylim + [-1 1]);
                title(this.marker);
                Core_UI.beautifyFig(f, 'light');
                Core_UI.addBeautifyMenu(f);
                f.Visible = 'on';
            end
        end
        
        function [time_utc, tropo, tropo_std, tropo_dt_std] = getBestTropo(this, tropo_par, show_fig)            
            % Get optimal tropospheric parameter
            %
            % INPUT
            %   tropo_par   i.e. 'ztd', 'zwd', ...
            %   tropo_par   flag for debugging
            %            
            %
            % OUTPUT
            %   time_utc        time of tropo (UTC time)
            %   tropo           optima tropo requested parameter from set
            %   tropo_std       std of the used data, use it to filter tropo 
            %                   e.g. ZTD tropo(tropo_std > 0.006) = nan; % remove under 1mm of PWV error
            %   tropo_dt_std    std of the data samples vs Dt (epoch since the last estimation of the sample)
            %                   err = abs(tropo_sample - tropo) 
            %                   data is considered when err < 6cm
            %                   Each epoch = <rate> seconds
            %
            % SINTAX
            %   [time_utc, tropo, tropo_std] = this.getBestTropo(tropo_par, show_fig)
            %
            % EXAMPLE
            %   tic; [time_utc, tropo, tropo_std, tropo_dt_std] = tm1(2).getBestTropo('ztd'); toc; 
            %   fh = figure; plot((1:numel(tropo_dt_std)) * 30 / 3600, tropo_dt_std * 1e2, 'LineWidth', 2); title('STD vs \Deltat');
            %   xlabel('\Deltat  [hours]'); ylabel('std [cm]');
            %   Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'true');
            if nargin < 3 || isempty(show_fig)
                show_fig = false;
            end
            
            win_size = 11;
            min_time = inf;
            max_time = -inf;
            rate = [];
            neps = [];
            time_utc = [];
            tropo = [];
            tropo_std = [];
            n_samples = numel(this.data_set);
            xyz = nan(n_samples, 3);
            enu = nan(n_samples, 3);
            if ~isempty(this.data_set)
                % Loop on samples to get size informations
                for i = 1 : n_samples
                    min_time = min(min_time, round(this.data_set(i).utc_time(1) * 86400));
                    max_time = max(max_time, round(this.data_set(i).utc_time(end) * 86400));
                    rate = [rate; diff(round(this.data_set(i).utc_time(:) * 86400))];                    
                    neps = [neps; numel(this.data_set(i).utc_time)]; % number of epochs per sample
                    
                    xyz(i, :) = this.data_set(i).coo.getXYZ();
                    enu(i, :) = this.data_set(i).coo.getENU();
                end
                
                % loop on samples to generate a synced matrix of all the results
                rate = median(rate);
                n_epochs = (max_time - min_time + rate) / rate;
                time_utc = (min_time : rate : max_time)' ./ 86400;
                if show_fig
                    fh = figure('Visible', 'off');
                end
                tropo_win = nan(n_epochs, win_size);
                tropo_mean = nan(n_samples, 1);
                max_ep = 0;
                for i = 1 : n_samples
                    s = mod(i - 1, win_size) + 1;
                    idt = (round(this.data_set(i).utc_time * 86400) - min_time) / rate + 1;
                    max_ep = max(max_ep, idt(end) - idt(1) + 1);
                    tropo_win(idt, s) =  this.data_set(i).(lower(tropo_par));
                    if show_fig
                        plot(idt, tropo_win(idt, s), ':k'); hold on;
                    end

                    idt = idt([1 : min(3600 / rate, numel(idt))  min(3600/rate + 1, numel(idt)) : (length(idt) - 3600 / rate)]); % discard the last hour but keep at least one hour                    
                    tropo_mean(i) = mean(tropo_win(idt, s), 'omitnan');
                    
                    if show_fig
                        plot(idt, tropo_win(idt, s), 'LineWidth', 2); hold on;
                    end
                end
                
                if show_fig
                    drawnow
                    title('All the estimated Tropo');
                    %figure; plot(tropo_mean.*1e2); hold on; yyaxis right; plot((enu(:,3)- median(enu(:,3))).*1e2); title('Tropo vs UP');
                end
                
                tropo = nan(n_epochs, 1);
                tropo_std = nan(n_epochs, 1);
                tropo_win = tropo_win';
                fprintf('Compute optimal tropo - epoch: %6d / %6d\n', 0, n_epochs);
                for i = 1 : n_epochs
                    fprintf('%s%6d / %6d', 8 .* ones(1,15, 'uint8'), i, n_epochs);
                    tmp = noNaN(zero2nan(tropo_win(:,i)));
                    [res_tmp, id_res] = sort(abs(tmp - median(tmp))); % sorted residuals w.r.t. median
                    min_n_data = 3;
                    n_out = 3;
                    
                    tmp = res_tmp([1 : min(numel(res_tmp), min_n_data) min(numel(res_tmp), min_n_data + 1) : numel(res_tmp) - (n_out + 1)]) + median(tmp);
                    tropo(i) = mean(tmp);
                    tropo_std(i) = std(tmp);
                end
                fprintf('\n');
                tropo_std(tropo_std == 0) = max(tropo_std);
                
                if nargout > 3
                    tropo_dt_std = zeros(max_ep, 1);
                    tropo_dt_std_n = zeros(max_ep, 1);
                    for i = 1 : n_samples
                        s = mod(i - 1, win_size) + 1;
                        idt = (round(this.data_set(i).utc_time * 86400) - min_time) / rate + 1;
                        dt_var = nan2zero(tropo_win(s, idt)' - tropo(idt));
                        % If the error is too high this is clearly an outlier (exclude  over 6 cm)
                        if any(abs(dt_var) > 6 ./ 1e2)
                            dt_var(abs(dt_var) > 6 ./ 1e2) = nan;
                        end
                        tropo_dt_std(idt(end) - idt + 1) = tropo_dt_std(idt(end) - idt + 1) + nan2zero(dt_var).^2;
                        tropo_dt_std_n(idt(end) - idt + 1) = tropo_dt_std_n(idt(end) - idt + 1) + ~isnan(dt_var);
                    end
                    tropo_dt_std = sqrt(tropo_dt_std ./ (tropo_dt_std_n + eps));
                end
                
                if show_fig
                    Core_UI.beautifyFig(fh, 'true');
                    fh.Visible = 'on';
                    
                    figure(fh); 
                    plot(tropo, 'k', 'LineWidth', 2);
                    
                    fh = figure('Visible', 'off');
                    nnan = ~isnan(tropo);
                    x = [time_utc(nnan); flipud(time_utc(nnan))];
                    patch = fill(x, [tropo(nnan) + tropo_std(nnan); flipud(tropo(nnan) - tropo_std(nnan))], 'r');
                    set(patch, 'edgecolor', 'none');
                    set(patch, 'FaceAlpha', 0.3);
                    hold on;
                    plot(time_utc, tropo, 'LineWidth', 2);
                    hold off;
                    title('Tropo with solution std')
                    pause(0.1);
                    Core_UI.beautifyFig(fh, 'true');
                    fh.Visible = 'on';
                end
            end
        end
    end
    
    methods (Static)
        function tm = test
            tm = Tropo_Mat_Reader('/Volumes/Data/goGPS_data/project/Lampo/OUT_NRT', 'LP07');
            tm.showAllTropoPar('ztd'); tm.showENU();
        end
    end
end
