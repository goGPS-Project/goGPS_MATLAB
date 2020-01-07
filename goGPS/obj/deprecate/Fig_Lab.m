%   CLASS Fig_Lab
% =========================================================================
%
% DESCRIPTION
%   Set of useful functions for satellite related computations
%
% EXAMPLE
%
% FOR A LIST OF CONSTANTs and METHODS use doc Fig_Lab


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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

classdef Fig_Lab < handle

    methods (Static)
        function this = Fig_Lab()
            % empty creator
        end

    end
    % =========================================================================
    %  CONSTELLATION MANAGEMENT
    % =========================================================================

    methods (Static) % Public Access Generic utilities

        function plotENU(time, enu, spline_base, color_type)
            narginchk(2,4);

            if nargin == 2
                spline_base = 0;
            end
            if nargin < 4
                color_type = 0;
            end
            
            if ~isempty(time)
                if time.getRate > 86000
                    date_style = 'dd mmm yyyy';
                else
                    date_style = 'dd mmm yyyy HH:MM';
                end
            end

            if size(enu, 3) == 2
                enu_var = enu(:,:,2);
                enu = enu(:,:,1);
            else
                enu_var = [];
            end
            
            spline_base = spline_base * (size(enu,1) > spline_base);

            persistent m_enu
            
            if size(color_type, 2) == 3
                if size(color_type, 1) < 3
                    color_order = repmat(color_type, 3, 1);
                else
                    color_order = color_type;
                end
                color_type = true;
            else
                if (color_type == 0)
                    fh = figure();
                    m_enu = [];
                end
                %maximizeFig(gcf);
                color_order = handle(gca).ColorOrder;
                if color_type == 1
                    color_order = (min(1, max(0, 0.0 + (color_order - repmat(mean(color_order, 2),1,3)) * 1 + repmat(mean(color_order, 2),1,3))));
                end
                if color_type == -1
                    color_order = color_order * 0;
                end
            end
            
            % prepare data
            if isempty(m_enu)
                m_enu = [mean(enu(~isnan(enu(:,1)),1)) mean(enu(~isnan(enu(:,2)),2)) mean(enu(~isnan(enu(:,3)),3))];
            end
            data_e = (enu(:,1) - m_enu(:,1))*1e2;
            data_n = (enu(:,2) - m_enu(:,2))*1e2;
            data_u = (enu(:,3) - m_enu(:,3))*1e2;
            %data_e = enu(:,1) * 1e2;
            %data_n = enu(:,2) * 1e2;
            %data_u = enu(:,3) * 1e2;

            if spline_base > 0
                if isempty(time) || (isa(time, 'GPS_Time') && time.isempty())
                    data_e_s = splinerMat(1:numel(data_e), data_e(:), spline_base,0);
                    data_n_s = splinerMat(1:numel(data_e), data_n(:), spline_base,0);
                    data_u_s = splinerMat(1:numel(data_e), data_u(:), spline_base,0);

                else
                    [~, ~, ~, data_e_s] = splinerMat(time.getEpoch(~isnan(data_e)).getGpsTime(), data_e(~isnan(data_e)), spline_base * time.getRate(),0, time.getGpsTime());
                    [~, ~, ~, data_n_s] = splinerMat(time.getEpoch(~isnan(data_n)).getGpsTime(), data_n(~isnan(data_e)), spline_base * time.getRate(),0, time.getGpsTime());
                    [~, ~, ~, data_u_s] = splinerMat(time.getEpoch(~isnan(data_u)).getGpsTime(), data_u(~isnan(data_e)), spline_base * time.getRate(),0, time.getGpsTime());
                end
            end

            if isempty(time) || (isa(time, 'GPS_Time') && time.isempty())
                subplot(3,1,1); if color_type; hold on; end
                plot(data_e, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                if spline_base > 0
                    plot(data_e_s, '--', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
                end
                if ~isempty(enu_var)
                    Fig_Lab.plotError(data_e, enu_var(:,1), color_order(1,:));
                end

                subplot(3,1,2);
                plot(data_n, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(2,:));  hold on;
                if spline_base > 0
                    plot(data_n_s, '--', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
                end
                if ~isempty(enu_var)
                    Fig_Lab.plotError(data_n, enu_var(:,2), color_order(2,:));
                end

                subplot(3,1,3);
                plot(data_u, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(3,:));  hold on;
                if spline_base > 0
                    plot(data_u_s, '--', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
                end
                if ~isempty(enu_var)
                    Fig_Lab.plotError(data_u, enu_var(:,3), color_order(3,:));
                end
            else
                t = time.getMatlabTime;
                t = t(1:numel(data_e));
                subplot(3,1,1);
                plot(t, data_e, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                if spline_base > 0
                    plot(t, data_e_s, '--', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
                end
                if ~isempty(enu_var)
                    Fig_Lab.plotError(t, data_e, enu_var(:,1), color_order(1,:));
                end
                setTimeTicks(4,date_style);

                subplot(3,1,2);
                plot(t, data_n, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(2,:));  hold on;
                if spline_base > 0
                    plot(t, data_n_s, '--', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
                end
                if ~isempty(enu_var)
                    Fig_Lab.plotError(t, data_n, enu_var(:,2), color_order(2,:));
                end
                setTimeTicks(4,date_style);

                subplot(3,1,3);
                plot(t, data_u, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(3,:));  hold on;
                if spline_base > 0
                    plot(t, data_u_s, '--', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
                end
                if ~isempty(enu_var)
                    Fig_Lab.plotError(t, data_u, enu_var(:,3), color_order(3,:));
                end
                setTimeTicks(4,date_style);
            end
            subplot(3,1,1); ax(1) = gca;
            grid on;
            if color_type
                title(sprintf('%s vs std %.2f [cm]', handle(gca).Title.String, std(data_e(~isnan(data_e)))));
            else
                title(sprintf('std %.2f [cm]', std(data_e(~isnan(data_e)))));
            end
            ylabel('\Delta EAST [cm]', 'FontWeight', 'Bold')

            subplot(3,1,2); ax(2) = gca;
            grid on;
            if color_type
                title(sprintf('%s vs std %.2f [cm]', handle(gca).Title.String, std(data_n(~isnan(data_n)))));
            else
                title(sprintf('std %.2f [cm]', std(data_n(~isnan(data_n)))));
            end
            ylabel('\Delta NORTH [cm]', 'FontWeight', 'Bold')

            subplot(3,1,3); ax(3) = gca;
            grid on;
            if color_type
                title(sprintf('%s vs std %.2f [cm]', handle(gca).Title.String, std(data_u(~isnan(data_u)))));
            else
                title(sprintf('std %.2f [cm]', std(data_u(~isnan(data_u)))));
            end
            linkaxes(ax,'x');
            ylabel('\Delta UP [cm]', 'FontWeight', 'Bold')
        end
    end
    
    methods (Static) % Public Access
        function plotError(x, y, e, color)
            % SYNTAX:
            %   plotError(x, y, e, color)
            %   plotError(y, e, color)
            %   plotError(y, e)
            
            switch nargin
                case 4
                    patchColor = min(color + 0.3, 1);
                case 3
                    patchColor = min(e + 0.3, 1);
                    e = y;
                    y = x;
                    x = 1:length(y);
                case 2
                    patchColor = min([0.5 0.5 0.5] + 0.3, 1);
                    e = y;
                    y = x;
                    x = 1:length(y);
            end
            % plotError(x,y, hline(1).Color)
            plot(x, y + e, x, y - e, 'color', patchColor);
            patch([x(~isnan(e)) fliplr(x(~isnan(e)))], [y(~isnan(e)) + e(~isnan(e)); flipud(y(~isnan(e)) -e(~isnan(e)))], 1, ...
                'facecolor',patchColor, ...
                'edgecolor','none', ...
                'facealpha', 0.1);
            
        end
        
        function [time, enu, xyz] = plotExtractionPos(file_name_extraction, hold_on, plot_list)
            % SYNTAX:
            %   plotExtractionPos(<file_name_extraction>, plot_list)
            %
            % INPUT:
            %   file_name_extraction = name of the extraction file
            %   plot_list = list of plots to show: valid values [1 2]
            %      1:    plot ENU
            %      2:    plot modulus (3D) error
            %
            % EXAMPLE:
            %    state = Core.getCurrentSettings();
            %
            %    file_name_base = fnp.dateKeyRep(fnp.checkPath(fullfile(state.getOutDir(), sprintf('%s_%s${YYYY}${DOY}', marker_trg, marker_mst))), sss_date_start);
            %    file_name_base = fnp.dateKeyRep(sprintf('%s_${YYYY}${DOY}',file_name_base), sss_date_stop);
            %    file_name = sprintf('%s_extraction.txt', file_name_base);
            %    Fig_Lab.plotExtractionPos(file_name);
            %
            % DESCRIPTION:
            %   Plot the results contained into the extraction file of a batch execution

            narginchk(1,3);
            if nargin < 2
                hold_on = 0;
            end
            if nargin < 3
                plot_list = 1;
            end

            log = Core.getLogger();

            fid = fopen(file_name_extraction,'r');

            if (fid < 0)
                log.addError(['Failed to open ', file_name_extraction]);
            else
                txt = fread(fid,'*char')';
                log.addMessage(['Reading ', file_name_extraction]);
                fclose(fid);

                data = sscanf(txt','%4d-%3d  %2d/%2d/%2d    %2d:%2d:%6f   %14f   %14f   %14f   %14f   %14f   %14f\n');
                if numel(data) < 14
                    data = sscanf(txt','%4d-%3d  %4d/%2d/%2d    %2d:%2d:%6f   %14f   %14f   %14f   %14f   %14f   %14f\n');
                end
                data = reshape(data, 14, numel(data)/14)';
                time = GPS_Time(datenum([data(:,1), data(:,4:8)]));
                xyz = data(:,9:11);
                enu = data(:,12:14);
                clear data;

                m_xyz = mean(xyz);

                if ~isempty(intersect(plot_list, 1))
                    if time.getRate < 3600
                        %spline_base = double(((time.last.getUnixTime - time.first.getUnixTime) / 86400)) * round(86400 / time.getRate() / 24 * 1);
                        spline_base = 0;
                    else
                        spline_base = 0;
                    end
                    Fig_Lab.plotENU(time, enu, spline_base, hold_on);
                end

                if ~isempty(intersect(plot_list, 2))
                    color_order = handle(gca).ColorOrder;
                    fh2 = figure(); maximizeFig(fh2);
                    data = sqrt((xyz(:,1) - m_xyz(:,1)).^2 + (xyz(:,2) - m_xyz(:,2)).^2 + (xyz(:,3) - m_xyz(:,3)).^2) * 1e2;
                    data_s = iif(numel(data)>14, splinerMat(time.getGpsTime(), data(:), min(numel(data),14) * time.getRate(),0), data);
                    plot(time.getMatlabTime, data, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(4,:));  hold on;
                    plot(time.getMatlabTime, data_s, 'k--');
                    grid on; setTimeTicks(4,'dd mmm yyyy');
                    title(sprintf('3D - mean %.2f [cm]', mean(data)));
                    ylabel('displacement [cm]', 'FontWeight', 'Bold')
                end
            end
        end

        function [time, enu, err] = plotBatchPos(file_name_pos, hold_on)
            % SYNTAX:
            %   plotBatchPos(<file_name_pos>, plot_list)
            %
            % INPUT:
            %   file_name_pos = name of the batch position file
            %
            % EXAMPLE:
            %    state = Core.getCurrentSettings();
            %
            %    file_name_base = fnp.dateKeyRep(fnp.checkPath(fullfile(state.getOutDir(), sprintf('%s_%s${YYYY}${DOY}', marker_trg, marker_mst))), sss_date_start);
            %    file_name_base = fnp.dateKeyRep(sprintf('%s_${YYYY}${DOY}',file_name_base), sss_date_stop);
            %    file_name = sprintf('%s_position.txt', file_name_base);
            %    Fig_Lab.plotBatchPos(file_name);
            %
            % DESCRIPTION:
            %   Plot the results contained into the position file of a batch execution

            if nargin < 2
                hold_on = 0;
            end

            log = Core.getLogger();

            fid = fopen(file_name_pos,'r');

            if (fid < 0)
                log.addError(['Failed to open ', file_name_pos]);
            else
                txt = fread(fid,'*char')';
                log.addMessage(['Reading ', file_name_pos]);
                fclose(fid);

                data = sscanf(txt(207:end)',' %4d-%3d  %2d/%2d/%2d    %2d:%2d:%6f   %14f   %14f   %14f  %14f\n');
                data = reshape(data, 12, numel(data)/12)';
                time = GPS_Time(datenum([data(:,1), data(:,4:8)]));
                enu = data(:,9:11);
                err = data(:,12);
                clear data;

                Fig_Lab.plotENU(time, enu, 0*round(1*86400/time.getRate()), hold_on);
            end
        end

        function test()
            % test plots
            Fig_Lab.plotExtractionPos('../data/project/default_DD_batch/out/CAC2_CAC3_2016333_2016363_extraction.txt');
            Fig_Lab.plotBatchPos('../data/project/default_DD_batch/out/CAC2_CAC3_2016333_2016363_position.txt');
        end
    end
end
