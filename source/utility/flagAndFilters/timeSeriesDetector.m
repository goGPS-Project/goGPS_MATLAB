% timeSeriesDetector - Detect and filter outliers in a time series data
%
% SYNTAX:
%   [smoothed_data, smoothed_spline, data_filtered] = timeSeriesDetector(time, data, big_win, thr_big, small_win, thr_small)
%
% INPUT:
%   time        : Time vector (can be in GPS_Time format or MATLAB time format)
%   data        : Time series data as a column vector second column can be variances
%   big_win     : Window size for detecting large outliers (in days)
%   thr_big     : Threshold multiplier for detecting large outliers
%   small_win   : (Optional) Window size for detecting small outliers (in days)
%   thr_small   : (Optional) Threshold multiplier for detecting small outliers
%
% OUTPUTS:
%   smoothed_data   : Data after removing detected outliers
%   smoothed_spline : Smoothed spline representation of the data
%   data_filtered   : Data with detected outliers set to NaN
%
% DESCRIPTION:
%   This function detects and filters outliers in a time series data using
%   a two-stage robust filtering approach. In the first stage, large outliers
%   are detected and removed using a large window size. In the optional second
%   stage, smaller outliers are detected using a smaller window size. The function
%   also computes a smoothed spline representation of the data.
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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

function [smoothed_data, id_ko, smoothed_spline, data_filtered] = timeSeriesDetector(time, data, big_win, thr_big, small_win, thr_small, varargin)

    % Convert time if in GPS_Time format
    if isa(time,'GPS_Time')        
        time = time.getMatlabTime;
    end
    
    % Initialize
    data_filtered = data(:,1);
    
    % Large outliers detection and filtering
    smoothed_data = robFilt(time, [data_filtered data(:,2)], big_win);
    residuals = data_filtered - smoothed_data;
    id_ko = (abs(residuals) / robStd(residuals)) > thr_big;
    data_filtered(id_ko) = nan;
    
    flag_2step = nargin > 5 && ~isempty(small_win) && ~isempty(thr_small);
    flag_plot = (nargin == 7  && strcmp(varargin{1}, '-plot')) || (nargin == 5  && strcmp(small_win, '-plot'));

    % Small outliers detection (if parameters provided)
    if flag_2step
        smoothed_data = robFilt(time, [data_filtered data(:,2)], small_win);
        residuals = data(:,1) - smoothed_data;
        id_ko = (abs(residuals) / robStd(data_filtered - smoothed_data)) > thr_small;
    end
        
    % Compute smoothed_spline (if output requested)
    if nargout > 2 || flag_plot
        if flag_2step
            smoothed_spline = splinerMat(time, smoothed_data, small_win, 1e-5);
        else
            smoothed_spline = splinerMat(time, smoothed_data, big_win, 1e-5);
        end
    end
    
    % Prepare data_filtered for output (if output requested)
    if nargout > 3 || flag_plot
        data_filtered = data(:,1);
        data_filtered(id_ko) = nan;
    end

    % Plotting
    if flag_plot
        fh = figure;
        plot(time, data(:,1), '.', 'color', [0.5 0.5 0.5]);
        hold on;
        plotSep(time, data_filtered, '.-', 'color', Core_UI.getColor(1));
        plotSep(time, smoothed_data,'.-');
        plot(time, smoothed_spline,'.-');
        setAllLinesWidth(2);
        ylim(minMax(data_filtered(:)) + 0.1*[-1 1] * diff(minMax(data_filtered(:))));
        legend('Original Data', 'FilteredData', 'robFilt', 'robSpline');
        Core_UI.addLineMenu(fh);
    end
end
