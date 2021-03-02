function [lid_ko, clean_data] = getOutliers(data, sigma, win_size, n_sigma_win, n_sigma_jmp, range)
% Return outliers in a relatively constant time series
%
% INPUT
%   x             flag array (as logical)
%   sigma         points with forced splits (as logical) (default auto)
%   win_size      size of the moving window (default 7)
%   n_sigma_win   n_sigma thr to consider the moving window as affected by outliers (default 2)
%   n_sigma_jmp   n_sigma to detect a jmp (outlier) (default 4)
%   range         min max acceptable values (default auto)
%
% SYNTAX
%   [id_ko, clean_data] = getOutliers(data, sigma, win_size, n_sigma_win, n_sigma_jmp, range)
%


%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti
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

    % Compute thresholds - empirical but it works (outliers should be less than 75% percent)
    if nargin < 2 || isempty(sigma)
        thr_perc = 0.75;
        tmp = diff(x_orig);
        tmp(abs(tmp) > 3 * perc(noNaN(abs(tmp)), thr_perc)) = nan;
        sigma = std(tmp, 'omitnan');
        sigma(isnan(sigma)) = median(sigma, 'omitnan');
        clear tmp;
    end
    if nargin < 3 || isempty(win_size)
        win_size = 7;                % window size for outlier detection
    end
    if nargin < 4 || isempty(n_sigma_win)
        n_sigma_win = 3;
    end
    if nargin < 5 || isempty(n_sigma_jmp)
        n_sigma_jmp = 4;
    end
    mu_win = 3 * win_size;
    
    clean_data = data(:,1);
    
    % Initialization
    
    % Determine maximum acceptable value range
    if nargin < 6 || isempty(range)
        mov_range = perc(noNaN(diff(clean_data)), [0.2 0.5 0.8]);
        mov_range = mov_range([1 3]) + [-40 40]' .* diff(mov_range);
        
        tmp = nan2zero(medfilt_mat(clean_data, round(win_size)));
        lid_ko = clean_data < (tmp + mov_range(1)) | (clean_data > (tmp + mov_range(2)));
        
        range = perc(clean_data, [0.2 0.5 0.8]);
        range = range([1 3]) + [-6 6]' .* diff(range);
    else
        lid_ko = true(size(clean_data));
    end
    lid_ko = (clean_data < -1e4 | clean_data > 1e+4) | (lid_ko & (clean_data < range(1) | clean_data > range(2))); % if it is around the normal range 80% of data it's ok
    clean_data(lid_ko) = tmp(lid_ko);

    % Detection of outliers from formal sigma (usually underestimated)
    lid_ko = lid_ko | isnan(clean_data(:,1));
    if size(data,2) == 2
        % there is also std
        lid_ko(data(:,2) > 6 * sigma) = true;
        clean_data(lid_ko) = nan;
    end
    
    % Outlier parameters
    jmp_ok = n_sigma_jmp * sigma;     % A jump smaller than this threshold is ok
    sigma_ok = n_sigma_win * sigma;   % A sigma under this level in a moving window of min_win size is acceptable

    % compute running std
    tmp_win_size = min(numel(clean_data), win_size);
    std_ko = movstd(clean_data, tmp_win_size, 'Endpoints', 'discard');
    std_ko(isnan(std_ko)) = 1e10;
    std_ko = std_ko > sigma_ok;
    std_ko = [std_ko; repmat(std_ko(end), tmp_win_size - 1, 1)];
    % Start outlier detection
    for k = 1:2
        for i = 1:length(clean_data)
            if ~lid_ko(i)
                % Check for outliers ------------------------------------------
                if std_ko(i)
                    win = clean_data(i : min(length(clean_data), i -1 + win_size)); % window in the future
                    while (std(win, 'omitnan') > sigma_ok) % if future obs are unstable shrink the window
                        win(end) = [];
                    end
                    
                    % If the window have been shortened there is probably on outlier in the window
                    if sum(not(isnan(win))) < 5
                        % check the current epoch with a window in the past
                        win = clean_data(max(1, i - win_size + 1) : i);
                        if any(isnan(win))
                            mu = median(clean_data(max(1, i - mu_win + 1):(i-1)), 'omitnan');
                            if not(any(mu))
                                mu = median(clean_data(max(1, i - 5*mu_win + 1):(i-1)), 'omitnan');
                            end
                            win(isnan(win)) = mu; % set nan values t the actual running mean
                        end
                        % if the sigma is above thr or if the last epoch is jumping too much
                        if length(win) == 1 || std(win, 'omitnan') > sigma_ok || ((length(win) > 1) && abs(diff(win(end-1:end))) > jmp_ok)
                            % the current epoch is an outlier
                            %fprintf('Outlier detected %d\n', i);
                            clean_data(i) = nan;
                            lid_ko(i) = true;
                        end
                    end
                end
                % -------------------------------------------------------------
            end
        end
    end
    clean_data(lid_ko) = nan;
    %fprintf('%d outliers found\n', sum(id_ko))
    %figure; plot(clean_data, 'o');
end