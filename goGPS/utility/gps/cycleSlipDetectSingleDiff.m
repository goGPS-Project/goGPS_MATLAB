function [ph_R, ph_M] = cycleSlipDetectSingleDiff(ph_R, ph_M)
% SYNTAX:
%    [ph_R, ph_M] = cycleSlipDetectSingleDiff(ph_R, ph_M);
%
% DESCRIPTION:
%    clean phases observations using single differences R - M
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti
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

    ph_R(ph_R == 0) = NaN;
    ph_M(ph_M == 0) = NaN;

    [~, flag_array] = cleanPhaseObsSingleDiff_v1(ph_R - ph_M);

    ph_R(flag_array) = NaN;
    ph_M(flag_array) = NaN;
    ph_R(isnan(ph_R)) = 0;
    ph_M(isnan(ph_M)) = 0;
end

function [data, flag_array] = cleanPhaseObsSingleDiff_v1(data, thr)
    data = data';
    state = Core.getCurrentSettings();

    % remove big outliers
    if nargin == 1
        thr = [5 2 7];
    end

    % get dimension of the set
    n_obs = size(data, 1);
    n_set = size(data, 2);

    %% remove short arcs --------------------------
    flag_array = ~isnan(data);
    data_out = flag_array & ~remove_short_arcs(~isnan(data)', state.getMinArc())';
    data(data_out) = nan;
  
    % working on third derivative ---------------
    d3data = diff(data, 3);

    % data mean / std
    d3data_m = nan(n_set, 1);
    d3data_s = nan(n_set, 1);

    % get the biggest outliers
    for s = 1 : n_set
        d2data_ok = d3data(~isnan(d3data(:, s)),s);
        if ~isempty(d2data_ok)
            d3data_m(s) = mean(d2data_ok);
            d3data_s(s) = std(d2data_ok);
        end
    end
    d3data_ref = repmat(median(d3data, 2, 'omitnan'), 1, n_set);

    data_out(4:end, :) = data_out(4:end, :) | abs(d3data - d3data_ref) > thr(2) * median(d3data_s(~isnan(d3data_m)));
    data(data_out) = nan;

    data = data';
    data(~remove_short_arcs(~isnan(data), state.getMinArc())) = nan;
    flag_array = isnan(data);
end

function [data, flag_array] = cleanPhaseObsSingleDiff(data, thr)
    data = data';
    state = Core.getCurrentSettings();

    % get dimension of the set
    n_obs = size(data, 1);
    n_set = size(data, 2);

    % remove short arcs --------------------------
    flag_array = ~isnan(data);
    data_out = flag_array & ~remove_short_arcs(~isnan(data)', state.getMinArc())';
    data(data_out) = nan;

    % remove big outliers
    if nargin == 1
        thr = [2 6 7];
    end

    spline_base = 210;

    % working on first derivative ----------------
    ddata = diff(data);

    % data mean / std
    ddata_m = nan(n_set, 1);
    ddata_s = nan(n_set, 1);

    % get the biggest outliers
    for s = 1 : n_set
        ddata_ok = ddata(~isnan(ddata(:, s)),s);
        if ~isempty(ddata_ok)
            ddata_m(s) = mean(ddata_ok);
            ddata_s(s) = std(ddata_ok);
        end
    end
    data_out(2:end, :) = data_out(2:end, :) | abs(ddata) - median(ddata_m(~isnan(ddata_m))) > thr(1) * median(ddata_s(~isnan(ddata_m)));

    ddata(data_out(2:end,:)) = NaN;

    % working on second derivative ---------------

    d2data = diff(ddata);

    % data mean / std
    d2data_m = nan(n_set, 1);
    d2data_s = nan(n_set, 1);

    % get the biggest outliers
    for s = 1 : n_set
        d2data_ok = d2data(~isnan(d2data(:, s)),s);
        if ~isempty(d2data_ok)
            d2data_m(s) = mean(d2data_ok);
            d2data_s(s) = std(d2data_ok);
        end
    end
    data_out(3:end, :) = data_out(3:end, :) | abs(d2data) - mean(d2data_m(~isnan(d2data_m))) > thr(2) * median(d2data_s(~isnan(d2data_m)));

    % cleaning data ------------------------------

    d2data(data_out(3:end,:)) = NaN;
    ddata([false(1, n_set); isnan(d2data)]) = NaN;
    data([false(1, n_set); isnan(ddata)]) = NaN;

    data_out = ~isnan(data) & ~remove_short_arcs(~isnan(data)', state.getMinArc())';
    d2data(data_out(3:end, :)) = nan;
    ddata(data_out(2:end, :)) = nan;
    data(data_out) = nan;

    %figure; plot(d2data, '.-'); set(gcf, 'Name', 'data diff 2');
    %figure; plot(ddata, '.-'); set(gcf, 'Name', 'data diff');
    %figure; plot(data, '.-'); set(gcf, 'Name', 'data');
    %hold on; plot(tmp, 'ok', 'MarkerSize', 3);

    % removing splines from data derivatives -----

    % size of the win in seconds
    spline_win_size = 1800;
    data_smooth = d2data;
    t = 1 : n_obs;
    for s = 1 : n_set
        valid_obs = ~isnan(data_smooth(:, s));
        if ~isnan(d2data_m(s))
            data_smooth(valid_obs, s) = splinerMat(t(valid_obs)', data_smooth(valid_obs,s), ceil(spline_win_size/state.getProcessingRate()), 1e-9);
        end
    end

    %figure; plot(d2data, 'LineWidth', 2);
    %hold on; plot(data_smooth, 'k', 'LineWidth', 2);
    %figure; plot(d2data - data_smooth, '.-');

    % estimate the mediam behaviour of the differencial observations

    data_ref = (median((data_smooth), 2, 'omitnan'));
    data_ref(abs((diff(data_ref) - mean(diff(data_ref(~isnan(data_ref)))))) > 3*std(diff(data_ref(~isnan(data_ref))))) = NaN;
    data_ref = remove_short_arcs(data_ref', state.getMinArc())';

    % apply smoothing
    ttmp = t(1:end-2);
    ttmp(isnan(data_ref)) =  [];
    [~, ~, ~, data_ref] = splinerMat(ttmp', data_ref(~isnan(data_ref)), ceil(spline_base / state.getProcessingRate()), 1e-9, t(1:end-2)');

    %hold on; plot(data_ref, 'k', 'LineWidth', 2);

    % remove the median behaviour from the observations
    clean_ph = d2data - repmat(data_ref, 1, n_set);

    %figure; plot(clean_ph, 'LineWidth', 2);

    % last thr on the residuals
    data_out(3 : end, :) = data_out(3 : end, :) | abs(clean_ph) > thr(3) * std(clean_ph(~isnan(clean_ph)));
    data(data_out) = nan;
    data_out = ~isnan(data) & ~remove_short_arcs(~isnan(data)', state.getMinArc())';
    data(data_out) = nan;
    %figure; plot(data, '.-'); set(gcf, 'Name', 'data clean');

    %flag_array = (~flag_array & data_out)';
    data = data';
    flag_array = isnan(data);
end
