% ROBFILT applies a moving window filter to a time series, ignoring NaNs.
%
% INPUT:
%   ts: An n x 1 array representing the time series. May contain NaNs.
%   filter_size: An integer representing the size of the moving window.
%
% OUTPUT:
%   filtered: An n x 1 array representing the filtered time series.
%
% SYNTAX: filtered = robfilt(ts, filter_size)
%
% TEST:
% n_obs = 10000;
% n_bad = n_obs/4;
%
% % Original noise
% a = 5*randn(n_obs,1);
%
% % Add a random walk signal
% a = a + cumsum(randn(n_obs,1))/10 + 50*sin((1:n_obs)'*2*pi/n_obs*10) - 100*sin((1:n_obs)'*2*pi/n_obs*3) + - 200*sin((1:n_obs)'*2*pi/n_obs);
%
% % Keep original signal
% a0 = a;
%
% % Introduce some large outliers
% bad_id = randi(n_obs,n_bad,1);
% a(bad_id) = a(bad_id) + randn(n_bad,1)*randi(100,1,1);
% a(randi(n_obs,n_bad/4,1)) = randn(n_bad/4,1)*80;
%
% % Plot the data
% figure; tiledlayout(2,1);
% ax = nexttile();
% plot(a);
% hold on;
% plot(a0);
%
% % Apply the robust filter and time it
% tic;
% tmp = robFilt(a,31);
% toc;
%
% % Plot the filtered data
% plot(tmp, 'linewidth', 2);
% yl = ylim();
% legend('data + noise', 'original data', 'filtered data');
% ax(2) = nexttile();
% plot(a0-tmp, 'linewidth', 2);
% ylim(yl);
% legend('original - filtered data');
% linkaxes(ax, 'x');

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

function filtered = robFilt(ts, filter_size)
    % Get the size of the time series
    n = length(ts);
    filter_size = min(filter_size,n);

    % Determine the padding size and initial values for padding
    pad_size = floor(filter_size / 2);
    initial_values = initPad(ts(1:pad_size), pad_size);
    final_values = initPad(flipud(ts(end-pad_size+1:end)), pad_size);

    % Pad the data
    ts = [2*initial_values-flipud(ts(2:pad_size+1)); ts; 2*final_values-flipud(ts(end-pad_size:end-1))];

    % Initialize the output
    filtered = zeros(n, 1);

    % Set the maximum number of iterations for the robust adjustment
    max_iter = 50;
    
    % Loop over the time series
    for i = 1:n
        % Get the data in the window
        window = ts(i:i+filter_size-1);

        % Ignore NaNs
        window = window(~isnan(window));
        
        % If window is empty, continue to the next iteration
        if isempty(window)
            continue;
        end

        % Set the threshold for the robust adjustment
        thrs = 1.4826 * mad(window, 1);
        
        % Perform the robust adjustment
        j = 0;
        dt = 1e9;
        dt_prev = -1e9;
        while (j < max_iter && abs(dt - dt_prev) > 0.005)
            dt_prev = dt;
            w = ones(size(window));
            ares_n = abs(window - dt_prev) / thrs;
            idx_rw = ares_n > 1;
            if any(idx_rw)
                w(idx_rw) = 1 ./ ares_n(idx_rw).^2;
            end
            tmp = window .* w;
            dt = sum(tmp) / sum(w);
            j = j + 1;
        end
        
        % Save the result
        filtered(i) = dt;
    end
end

function initial_value = initPad(data, pad_size)
    % Determine the number of data points to use for the initial value    
    % Ignore NaNs
    data = data(~isnan(data));

    % If data is empty, return NaN
    if isempty(data)
        initial_value = NaN;
        return;
    end

    num_points = min(pad_size, 3);
    dt_prev = data(1);
    dt = robAdj(data(1:num_points)');
    num_points = num_points + 1;
    norm_factor = std(data);
    thr = min(pad_size,32);
    while abs(dt - dt_prev)/norm_factor < 1 && num_points < thr
        dt_prev = dt;
        dt = robAdj(data(1:num_points)');
        num_points = num_points + 1;
    end
    num_points = max(1,num_points - 1);
    dt = robAdj(data(1:num_points)');    
    initial_value = dt;
end

function m = mad(x, flag)
    % mad calculates the median absolute deviation of an array, ignoring NaNs.
    % It's implemented to avoid dependence on the Statistics and Machine Learning Toolbox.
    x = x(~isnan(x));
    m = median(abs(x - median(x)), flag);
end
