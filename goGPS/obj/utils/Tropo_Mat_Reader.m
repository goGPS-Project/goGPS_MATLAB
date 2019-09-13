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
%    |___/                    v 1.0 beta 3jp
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
            % Show ENU position variations
            if ~isempty(this.data_set)
                f = figure; f.Name = sprintf('%03d: BSL %s - %s', f.Number, ref.marker, this.marker); f.NumberTitle = 'off';
                        
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
                    setTimeTicks(10,'mmm-dd HH:MM');
                    h = ylabel(label{a}); h.FontWeight = 'bold';
                    grid on;
                    xlim(minMax(time));
                    yl = ylim();
                    yl(1) = min(-20, yl(1));
                    yl(2) = max(20, yl(2));
                    ylim(yl);
                end
                linkaxes(ax, 'x');
            end
        end
        
        function showENU(this, color)
            % Show ENU position variations
            if ~isempty(this.data_set)
                if nargin < 2
                    f = figure; f.Name = sprintf('%03d: ENU %s ', f.Number, this.marker); f.NumberTitle = 'off';
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
                    setTimeTicks(10,'mmm-dd HH:MM');
                    h = ylabel(label{a}); h.FontWeight = 'bold';
                    grid on;
                    xlim(minMax(time));
                end
                linkaxes(ax, 'x');
            end
        end
        
        function showAllTropoPar(this, tropo_par, color)
            if ~isempty(this.data_set)
                if nargin < 3
                    f = figure; f.Name = sprintf('%03d: %s %s ', f.Number, upper(tropo_par), this.marker); f.NumberTitle = 'off';
                end
                colormap(flipud(gat2));
                min_time = inf;
                max_time = -inf;
                for i = 1 : numel(this.data_set)
                    min_time = min(min_time, min(this.data_set(i).utc_time));
                    max_time = max(max_time, max(this.data_set(i).utc_time));
                    if nargin < 3
                        scatter(this.data_set(i).utc_time, this.data_set(i).(lower(tropo_par)) .* 1e2, 10, this.data_set(i).utc_time(end) - this.data_set(i).utc_time, 'filled'); hold on
                    else
                        plot(this.data_set(i).utc_time, this.data_set(i).(lower(tropo_par)) .* 1e2, 'color', color); hold on
                    end
                end
                xlim([min_time, max_time]);
                grid minor;
                setTimeTicks();
                title(this.marker);
            end
        end
        
        function [time_utc, tropo] = getBestTropo(this, tropo_par)
            
            time_utc = [];
            tropo = [];
            tropo_w = [];
            if ~isempty(this.data_set)
                for i = 1 : numel(this.data_set)
                    cur_time = round(this.data_set(i).utc_time * 86400);                    
                    cur_tropo = this.data_set(i).(lower(tropo_par));
                    
                    time_utc = [time_utc cur_time];
                    tropo_w = [tropo_w (cur_time(end) - cur_time)];
                    tropo = [tropo cur_tropo];                    
                end                
            end
            time_utc = time_utc / 86400;            
        end
    end
    
    methods (Static)
        function tm = test
            tm = Tropo_Mat_Reader('/Volumes/Data/goGPS_data/project/Lampo/OUT_NRT', 'LP07');
            tm.showAllTropoPar('ztd'); tm.showENU();
        end
    end
end
