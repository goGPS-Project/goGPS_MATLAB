%   CLASS Baseline
% =========================================================================
%
%
%   Class to manage baselines
%
% EXAMPLE
%   trg = Baseline();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Receiver

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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
classdef Baseline < handle
%% DRIVERS
    properties (SetAccess = public, GetAccess = public)
        name     % name
        time     % epoch_time [GPS_Time]
        enu      % ENU coordinates
        flag_q   % quality flag Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp
        nspe     % n_sat per epoch
        std_enu  % std of the ENU coordinates
        ratio    % ratio ???
    end
    
    methods (Static)
        function this = Baseline()
            % Creator method
            %
            % SYNTAX
            %   this = Baseline()
            this.init();
        end       
        
        function this = fromPosFile(pos_file, flag_mean)
            % EXAMPLE
            %   BaseLine.fromPosFile('/Volumes/Data/goGPS_data/project/GIMS_Test/out/rtklib/001N_002N_0001.pos');
            narginchk(1,2);
            if ~iscell(pos_file)
                [cur_dir, cur_name, cur_ext] = fileparts(pos_file);
                existing_files = dir(pos_file);
                pos_file = {};
                for i = 1 : numel(existing_files)
                    pos_file{i} = [fullfile(cur_dir, existing_files(i).name)];
                end
            end
            
            this = Baseline();
            log = Logger.getInstance();
            for i = 1 : numel(pos_file)
                cur_file = pos_file{i};
                log.addMarkedMessage(sprintf('Parsing %s', cur_file));
                [cur_dir, cur_name, cur_ext] = fileparts(cur_file);
                this.name = cur_name(1:9); % 'It is generated as <marker_name_1>_<marker_name_2>'
                
                fid = fopen(cur_file, 'r');
                if fid <= 0
                   log.addError(sprintf('"%s" not found!!!', cur_file)); 
                else
                    try
                        txt = fread(fid, '*char')';
                        fclose(fid);
                        
                        if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                            has_cr = true;  % The file has carriage return - I hate you Bill!
                        else
                            has_cr = false;  % The file is UNIX standard
                        end
                        
                        % get new line separators
                        nl = regexp(txt, '\n')';
                        if nl(end) <  (numel(txt) - double(has_cr))
                            nl = [nl; numel(txt)];
                        end
                        lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                        lim = [lim lim(:,2) - lim(:,1)];
                        while lim(end,3) < 3
                            lim(end,:) = [];
                        end
                        
                        % removing empty lines at end of file
                        while (lim(end,1) - lim(end-1,1))  < 2
                            lim(end,:) = [];
                        end
                        
                        % Ignore comment lines
                        lko = txt(lim(:,1)) == '%';
                        lim(lko, :) = [];
                        
                        if ~isempty(lim)
                            % find all the observation lines
                            n_epo = size(lim, 1);
                            % extract all the epoch lines
                            string_time = txt(repmat(lim(:,1),1,23) + repmat(0:22, n_epo, 1))';
                            % convert the times into a 6 col time
                            date = cell2mat(textscan(string_time,'%4f/%2f/%2f %2f:%2f:%6.3f'));
                            % import it as a GPS_Time obj
                            time = GPS_Time(date, [], true);
                            
                            string_data = txt(repmat(lim(:,1), 1, 122) + repmat(23:(lim(1,3)+1), n_epo, 1))';
                            data = reshape(sscanf(string_data, '%f'), 13, size(lim, 1))';
                            
                            % Outlier rejection
                            if size(data, 1) > 1
                                id_ok = (data(:, 4) == 1) & ... % Keep only fixed solutions
                                    max(data(:, 6:8)')' < 0.015 & ... % Keep formal error less than 1.5cm
                                    all(data(:, 6:8) > 0, 2) & ... % Keep formal error greater than 0
                                    all(abs(data(:, 1:3) - movmedian(data(:, 1:3), 13, 'omitnan')) < 0.01, 2); % Keep derivative < 1 cm
                                data = data(id_ok, :);
                                time = time.getEpoch(id_ok);
                            end
                            if nargin >= 2 && flag_mean
                                this.time.append(time.getCentralTime);
                                
                                %this.enu = [this.enu; sum(data(:, 1:3) ./ data(:, 6:8).^2) ./ sum(1 ./ data(:, 6:8).^2)];
                                this.enu = [this.enu; median(data(:, 1:3))];
                                this.flag_q = [this.flag_q; mean(data(:, 4))];
                                this.nspe = [this.nspe; min(data(:, 5))];
                                this.std_enu = [this.std_enu; max(data(:, 6:8))];
                                this.ratio = [this.ratio; max(data(:, 13))];
                            else
                                this.time.append(time);
                                this.enu = [this.enu; data(:, 1:3)];
                                this.flag_q = [this.flag_q; data(:, 4)];
                                this.nspe = [this.nspe; data(:, 5)];
                                this.std_enu = [this.std_enu; data(:, 6:8)];
                                this.ratio = [this.ratio; data(:, 13)];
                            end
                        end
                        
                    catch ex
                        log.addError(sprintf('Failed to import %s', cur_file));
                    end
                end
            end
            
        end
    end
    
    methods
        function init(this)
            this.name = '';    % name
            this.time = GPS_Time(); % time
            this.enu = [];     % ENU coordinates
            this.flag_q  = []; % quality flag 1 = best fixed, 2 = float, ...
            this.nspe = [];    % n_sat per epoch
            this.std_enu = []; % std of the ENU coordinates
            this.ratio = [];   % ratio ???
        end
        
        function [enu, time] = getENU(this)
            enu = this.enu;
            time = this.time.getCopy;
        end
    end
    
    methods
        function showBaselineENU(this, plot_relative_variation, one_plot)
            if (nargin < 3) || isempty(one_plot)
                one_plot = false;
            end
            if (nargin < 2) || isempty(plot_relative_variation)
                plot_relative_variation = true;
            end

            log = Logger.getInstance;
            [enu, time] = this.getENU();
            if size(enu, 1) > 1
                log.addMessage('Plotting positions');
                
                % prepare data
                if plot_relative_variation
                    enu = bsxfun(@minus, enu, median(enu, 'omitnan')) * 1e3;
                end
                t = time.getMatlabTime();
                
                f = figure; f.Name = sprintf('%03d: BSL ENU %s', f.Number, this.name); f.NumberTitle = 'off';
                color_order = handle(gca).ColorOrder;
                
                if ~one_plot, subplot(3,1,1); end
                plot(t, enu(:, 1), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                ax(3) = gca();
                if (t(end) > t(1))
                    xlim([t(1) t(end)]);
                end
                setTimeTicks(4,'dd/mm/yyyy HH:MM');
                if plot_relative_variation
                    h = ylabel('East [mm]'); h.FontWeight = 'bold';
                else
                    h = ylabel('East [m]'); h.FontWeight = 'bold';
                end
                grid minor;
                h = title(sprintf('Baseline %s \t\tstd E %.2f - N %.2f - U%.2f -', this.name, std(enu, 'omitnan')), 'interpreter', 'none'); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                
                if ~one_plot, subplot(3,1,2); end
                plot(t, enu(:, 2), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:));
                ax(2) = gca();
                if (t(end) > t(1))
                    xlim([t(1) t(end)]);
                end
                setTimeTicks(4,'dd/mm/yyyy HH:MM');
                if plot_relative_variation
                    h = ylabel('North [mm]'); h.FontWeight = 'bold';
                else
                    h = ylabel('North [m]'); h.FontWeight = 'bold';
                end
                
                grid minor;
                if ~one_plot, subplot(3,1,3); end
                plot(t, enu(:,3), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:));
                ax(1) = gca();
                if (t(end) > t(1))
                    xlim([t(1) t(end)]);
                end
                setTimeTicks(4,'dd/mm/yyyy HH:MM');
                if plot_relative_variation
                    h = ylabel('Up [mm]'); h.FontWeight = 'bold';
                else
                    h = ylabel('Up [m]'); h.FontWeight = 'bold';
                end
                
                grid minor;
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
                this.log.addMessage('Plotting a single point static position is not yet supported');
            end
        end
    end
end
