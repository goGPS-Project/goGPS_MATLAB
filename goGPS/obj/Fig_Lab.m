%   CLASS Fig_Lab
% =========================================================================
%
% DESCRIPTION
%   Set of useful functions for satellite related computations
%
% EXAMPLE
%
% FOR A LIST OF CONSTANTs and METHODS use doc goGNSS


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 3
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
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

    methods (Static) % Public Access

        function [time, enu, xyz] = plotExtractionPos(file_name_extraction, plot_list)
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
            %    state = Go_State.getCurrentSettings();
            %
            %    file_name_base = fnp.dateKeyRep(fnp.checkPath(fullfile(state.getOutDir(), sprintf('%s_%s${YYYY}${DOY}', marker_trg, marker_mst))), sss_date_start);
            %    file_name_base = fnp.dateKeyRep(sprintf('%s_${YYYY}${DOY}',file_name_base), sss_date_stop);
            %    file_name = sprintf('%s_extraction.txt', file_name_base);
            %    Fig_Lab.plotExtractionPos(file_name);
            %
            % DESCRIPTION:
            %   Plot the results contained into the extraction file of a batch execution

            narginchk(1,2);
            if nargin == 1
                plot_list = [1 2];
            end

            logger = Logger.getInstance();

            fid = fopen(file_name_extraction,'r');

            if (fid < 0)
                logger.addError(['Failed to open ', file_name_extraction]);
            else
                txt = fread(fid,'*char')';
                logger.addMessage(['Reading ', file_name_extraction]);
                fclose(fid);

                data = sscanf(txt','%4d-%3d  %2d/%2d/%2d    %2d:%2d:%6f   %14f   %14f   %14f   %14f   %14f   %14f\n');
                data = reshape(data, 14, numel(data)/14)';
                time = GPS_Time(datenum([data(:,1), data(:,4:8)]));
                xyz = data(:,9:11);
                enu = data(:,[12 13 14]);
                clear data;

                m_xyz = mean(xyz);
                m_enu = mean(enu);

                if ~isempty(intersect(plot_list, 1))
                    fh = figure(); maximizeFig(fh);
                    color_order = handle(gca).ColorOrder;

                    subplot(3,1,1); ax(1) = gca;
                    data = (enu(:,1) - m_enu(:,1))*1e3;
                    data_s = iif(numel(data)>4, splinerMat(time.getGpsTime(), data(:), min(numel(data),14) * time.getRate(),0), data);
                    plot(time.getMatlabTime, data, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(1,:));  hold on;
                    plot(time.getMatlabTime, data_s, 'k--');
                    grid on; setTimeTicks(4,'dd mmm yyyy');
                    title(sprintf('East - std %.2f [mm]', std(data)));
                    ylabel('displacement [mm]', 'FontWeight', 'Bold')

                    subplot(3,1,2); ax(2) = gca;
                    data = (enu(:,2) - m_enu(:,2))*1e3;
                    data_s = iif(numel(data)>14, splinerMat(time.getGpsTime(), data(:), min(numel(data),14) * time.getRate(),0), data);
                    plot(time.getMatlabTime, data, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(2,:));  hold on;
                    plot(time.getMatlabTime, data_s, 'k--');
                    grid on; setTimeTicks(4,'dd mmm yyyy');
                    title(sprintf('North - std %.2f [mm]', std(data)));
                    ylabel('displacement [mm]', 'FontWeight', 'Bold')

                    subplot(3,1,3); ax(3) = gca;
                    data = (enu(:,3) - m_enu(:,3))*1e3;
                    data_s = iif(numel(data)>4, splinerMat(time.getGpsTime(), data(:), min(numel(data),14) * time.getRate(),0), data);
                    plot(time.getMatlabTime, data, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(3,:));  hold on;
                    plot(time.getMatlabTime, data_s, 'k--');
                    grid on; setTimeTicks(4,'dd mmm yyyy');
                    title(sprintf('Up - std %.2f [mm]', std(data)));
                    linkaxes(ax,'x');
                    ylabel('displacement [mm]', 'FontWeight', 'Bold')
                end

                if ~isempty(intersect(plot_list, 2))
                    fh2 = figure(); maximizeFig(fh2);
                    data = sqrt((xyz(:,1) - m_xyz(:,1)).^2 + (xyz(:,2) - m_xyz(:,2)).^2 + (xyz(:,3) - m_xyz(:,3)).^2) * 1e3;
                    data_s = iif(numel(data)>14, splinerMat(time.getGpsTime(), data(:), min(numel(data),14) * time.getRate(),0), data);
                    plot(time.getMatlabTime, data, '.-', 'MarkerSize', 20, 'LineWidth', 2, 'Color', color_order(4,:));  hold on;
                    plot(time.getMatlabTime, data_s, 'k--');
                    grid on; setTimeTicks(4,'dd mmm yyyy');
                    title(sprintf('3D - mean %.2f [mm]', mean(data)));
                    ylabel('displacement [mm]', 'FontWeight', 'Bold')
                end
            end
        end

        function test()
            % test plots
            Fig_Lab.plotExtractionPos('../data/project/default_DD_batch/out/CAC2_CAC3_2016333_2016363_extraction.txt');
        end
    end
end
