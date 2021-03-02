function [jump_list, lid_ko, data,  running_mean] = getJumps(data, rate_scale)
% Return outliers in a relatively constant time series
%
% INPUT
%   data          data (column, eventually it can have a second column containing std)
%   rate_scale    24/rate_scale => used for scaling window size (1 => daily, 24 => hourly ...)
%
% OUTPUT
%   
% SYNTAX
%   [jump_list, id_ko, data, running_mean] = getJumps(data, rate_scale)
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
    flag_smooth = true; % Apply a spline smoothing to rate (duration 3 days)
    spline_base = 7; % in days
    
    if size(data,2) == 2
        data_std = data(:,2);
        data = data(:,1);
    else
        data_std = 0;
    end
    med_data = median(data, 'omitnan');
    data = data - med_data;

    % Compute thresholds - empirical but it works (outliers should be less than 75% percent)
    thr_perc = 0.75;
    tmp = diff(data);
    tmp(abs(tmp) > 3 * perc(noNaN(abs(tmp)), thr_perc)) = nan;
    thr = std(tmp, 'omitnan');
    clear tmp;

    % Outlier detection
    if any(data_std)
        [lid_ko, data] = getOutliers([data data_std], thr);
    else
        [lid_ko, data] = getOutliers(data, thr);
    end

    flag_running = true;
    % Initialization
    
    % Jump parameters
    run_win = 5 * rate_scale;   % window size for running mean
    cumsum_thr = 2;             % above this threshold, cumulate residuals
    cum_thr = 6;                % above this threshold the cumulate TRIGGERS A JUMP

    % Recompute a running_threshold
    % when the noise is higher the discrimination of jumps is lower
    thr = movstd(diff(data), run_win, 'omitnan');
    if not(any(thr))
        thr(:) = nan2zero(std(data, 'omitnan'));
    end
        
    q_win = floor(run_win / 2);
    thr((q_win + 1) : end - q_win) = min(thr((2*q_win +1):end), thr(1:end - (2*q_win)));
    
    mthr = median(thr, 'omitnan');
    
    thr(isnan(zero2nan(thr))) = mthr;
    run_win_variable = round(max(1,sqrt((thr - mthr)/mthr)) * run_win);
    run_win_variable = [run_win_variable; run_win_variable(end)];
    %thr = mthr + 0.1 * (thr - mthr);
    
    mu = data(1); % initial mean
    i = 1 + find(not(isnan(data(2:end))), 1, 'first');  % running index on observations
    [s_plus, s_minus] = deal(zeros(size(data,1),1)); % sensors
    [n_plus, n_minus] = deal(uint16(zeros(size(data,1),1))); % epochs since sensor cumulation started

    running_mean = data; % running mean

    j = 1;           % running index on jumps
    jump_list = 0;   % list of jumps
    while i <= length(data)
        % if jumps are not detected
        if isnan(data(i))
            s_plus(i) = s_plus(i-1);        % keep the old sensor
            s_minus(i) = s_minus(i-1);      % keep the old sensor
            n_plus(i)  = iif(s_plus(i)  > 0, n_plus(i-1)  + 1, 0);
            n_minus(i) = iif(s_minus(i) < 0, n_minus(i-1) + 1, 0);
            running_mean(i) = mu;
            i = i + 1;
        else
            % while no jumps are detected and i is running on observations
            while i <= length(data) && s_plus(i-1) <= (cum_thr * thr(i-1)) && s_minus(i-1) >= -(cum_thr * thr(i-1))
                val = iif(abs(data(i) - mu) > 0.9*(cum_thr * thr(i-1)), mu, data(i));
                s_plus(i)  = max(0, s_plus(i-1)  + (data(i) - mu) - (cumsum_thr * thr(i-1)));
                s_minus(i) = min(0, s_minus(i-1) + (data(i) - mu) + (cumsum_thr * thr(i-1)));

                n_plus(i)  = iif(s_plus(i)  > 0 , n_plus(i-1)  + 1, 0);
                n_minus(i) = iif(s_minus(i) < 0 , n_minus(i-1) + 1, 0);

                if flag_running
                    mu = nan2zero(mean(data(max(jump_list(j)+1, i + 1 - run_win_variable(i)):i), 'omitnan'));
                end
                running_mean(i) = mu;
                i = i + 1;
            end
            if i <= length(data)
                n_d = i - 1;
                j = j + 1;
                jump_list(j) = iif(s_plus(n_d) > (cum_thr * thr(i-1)), n_d - n_plus(n_d), n_d - n_minus(i-1));
                mu = nan2zero(mean(data(max(jump_list(j)+1, n_d - run_win_variable(n_d)):(n_d)), 'omitnan'));

                jmp_magnitude = running_mean(jump_list(j)) - mu;
                Core.getLogger.addMessage(sprintf('Jump detected at epoch %d, it happened in epoch %d, length: %d, old mean %f.2, new mean %f.2, for a jump of %.2f\n', n_d, jump_list(j), n_d - jump_list(j) + 1, running_mean(jump_list(j)), mu, jmp_magnitude),100);
                if (i <= length(data))
                    [s_plus(i-1), s_minus(i-1)] = deal(0);
                    [n_plus(i-1), n_minus(i-1)] = deal(0);
                end
            end
        end
        mu = nan2zero(mean(data(max(jump_list(j)+1, i - run_win_variable(i-1)):(i-1)), 'omitnan'));
        running_mean(i-1) = mu;
    end

    % compute running mean backwords
    if flag_running
        running_smooth = running_mean;
        robustness_perc = 0.90;
        n_sigma = 6;
        
        % flag small jumps (1-2 epochs are too small to be considered "real jumps"
        for j = fliplr(find(diff(jump_list) < 3))
            id_ko = jump_list(j)+1:jump_list(j+1);
            data(id_ko) = nan;
            lid_ko(id_ko) = true;
            jump_list(j) = jump_list(j+1);
            jump_list(j+1) = [];
        end
        
        last_win_id = jump_list(j)+1 : size(data,1); % last window
        % Compute a trend "robust" using the robustness_perc of data
        %step_median(last_win_id) = median(data(last_win_id), 'omitnan');
        [tmp, trend] = strongDeTrend(data(last_win_id), robustness_perc, 1-((1-robustness_perc)/2), n_sigma);
        if flag_smooth && sum(not(isnan((data(last_win_id))))) > 5
            if any(trend); running_mean(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
            try
                [~, ~, ~, running_smooth(last_win_id)] = splinerMat(1:numel(last_win_id), data(last_win_id), ceil(spline_base * rate_scale), 1e-5, 1:numel(last_win_id)); % medium splines
            catch
                if any(trend); running_smooth(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
            end
        else
            if any(trend); running_smooth(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
        end
        for j = 2 : numel(jump_list)
            last_win_id = jump_list(j-1)+1 : (jump_list(j));
            %step_median(last_win_id) = median(data(last_win_id), 'omitnan');
            [tmp, trend] = strongDeTrend(data(last_win_id), robustness_perc, 1-((1-robustness_perc)/2), n_sigma);
            if flag_smooth && sum(not(isnan((data(last_win_id))))) > 5
                if any(trend); running_mean(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
                try
                    [~, ~, ~, running_smooth(last_win_id)] = splinerMat(1:numel(last_win_id), data(last_win_id), ceil(spline_base * rate_scale), 1e-5, 1:numel(last_win_id)); % medium splines
                catch
                    if any(trend); running_smooth(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
                end
            else
                if any(trend); running_smooth(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
            end
        end
        
        %  remove very small jumps - recursive
        flag_rem_jmp = true;
        while flag_rem_jmp
            % Check the magnitude of the jump, with both moving mean and splines, if it is too small is a false alarm
            jmp_magnitude = min(diff(running_mean(jump_list(2:end) + repmat([0 1]', 1, numel(jump_list)-1))), ...
                                diff(running_smooth(jump_list(2:end) + repmat([0 1]', 1, numel(jump_list)-1))));
            % Jump smaller than half the noise in that window of data will be ignored
            id_small_jmp = abs(jmp_magnitude)' < max(0.4,(thr(jump_list(2:end))));
            if any(id_small_jmp)
                jump_list(find(id_small_jmp) + 1) = [];  % remove jmp
                
                % recompute running mean
                last_win_id = jump_list(end)+1 : size(data,1);
                %step_median(last_win_id) = median(data(last_win_id), 'omitnan');
                % remove the trend before computing the running mean (this help in case of very steep intervals of data)
                [tmp, trend] = strongDeTrend(data(last_win_id), robustness_perc, 1-((1-robustness_perc)/2), n_sigma);
                if flag_smooth && sum(not(isnan((data(last_win_id))))) > 5
                    if any(trend); running_mean(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
                    try
                        [~, ~, ~, running_smooth(last_win_id)] = splinerMat(1:numel(last_win_id), data(last_win_id), ceil(spline_base * rate_scale), 1e-5, 1:numel(last_win_id)); % medium splines
                    catch
                        running_smooth(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend;
                    end
                else
                    if any(trend); running_smooth(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
                end
                for j = 2 : numel(jump_list)
                    last_win_id = jump_list(j-1)+1 : (jump_list(j));
                    [tmp, trend] = strongDeTrend(data(last_win_id), robustness_perc, 1-((1-robustness_perc)/2), n_sigma);
                    if flag_smooth && sum(not(isnan((data(last_win_id))))) > 5
                        if any(trend); running_mean(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
                        try
                            [~, ~, ~, running_smooth(last_win_id)] = splinerMat(1:numel(last_win_id), data(last_win_id), ceil(spline_base * rate_scale), 1e-5, 1:numel(last_win_id)); % medium splines
                        catch
                            if any(trend); running_smooth(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
                        end
                    else
                        if any(trend); running_smooth(last_win_id) = movmean(tmp, run_win, 'omitnan') + trend; end
                    end
                end
            else
                flag_rem_jmp = false;
            end
        end
        running_mean = running_smooth;
    else
        last_win_id = jump_list(j)+1 : size(data,1);
        running_mean(last_win_id) = mean(data(last_win_id), 'omitnan');
    end
    
    running_mean = running_mean + med_data;
    data = data + med_data;
    % figure; plot(-s_minus); hold on; plot(s_plus); hold on; plot(cum_thr * thr)
end
