function [detrended_data, trend] = strongDeTrend(data, thr_perc, thr_perc_global, n_sigma)
% Returns the detrend per column with a certaint percentile of values in data
% The code uses only the data filtered by the requested percentile to estimate 
% an detred column by column
% With this estimation only the data within n_sigma range are used for the 
% robust detrend estimation
%
% SYNTAX:
%   smean = strongDeTrend(data, thr_perc, thr_perc_global, n_sigma)
%
% INPUT:
%   data                matrix of values
%   thr_perc            percentile requested [0 1] per column
%   thr_perc_global     percentile requested [0 1] global
%   n_sigma             number of sigmas
%
% OUTPUT:
%   detrended_data   strong detrend per column

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 (GReD srl) Andrea Gatti, Giulio Tagliaferro, Eugenio Realini
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti ...
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

    if nargin < 2 || isempty(thr_perc) 
        thr_perc = 0.99;
    end
    
    if nargin < 3 || isempty(thr_perc_global) 
        thr_perc_global = 0.97;
    end
    
    if nargin < 4 || isempty(n_sigma) 
        n_sigma = 5;
    end
    
    % filter data
    detrended_data = data;
    data(abs(data) >= perc(abs(data(:)), thr_perc_global)) = nan;
    x_in = (1 : size(data,1))';
    for c = 1 : size(data, 2)
        data(:,c) = data(:,c) - movmean(data(:,c), 5, 'omitnan'); % reduce the signal with moving mean(to remove the trend)
        data((abs(data(:, c)) >= perc(abs(data(:, c)), thr_perc)), c) = nan; % keep thr_perc of data
        thr = n_sigma * std(data(:, c), 'omitnan');
        lid_ok =  abs(data(:,c)) <= thr;
        try
            warning off
            trend = Core_Utils.interp1LS(x_in(lid_ok), detrended_data(lid_ok, c), 1, x_in); % use the current filtered data to compute the trend
            warning on;
        catch
            trend = 0;
        end
        data(:,c) = detrended_data(:,c) - trend;
        thr = n_sigma * std(data(lid_ok, c), 'omitnan'); % repeat the filtering with the new detrended data
        lid_ok =  abs(data(:,c)) <= thr;
        try
            warning off
            trend = Core_Utils.interp1LS(x_in(lid_ok), detrended_data(lid_ok, c), 1, x_in); % recompute trend with the new
            warning on;
        catch
            trend = 0;
        end
        detrended_data(:,c) = detrended_data(:,c) - trend;
    end    
end
