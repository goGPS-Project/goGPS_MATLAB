function [data, flag_array] = cleanPhaseObs(data_in, thr)
% SYNTAX:
%    [flagIntervals] = getOutliers(flags)
%
% DESCRIPTION:
%    returns start and end of flagged intervals
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

data = zero2nan(data_in');
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
    thr = [2 5 3];
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

data_ref = (median((d2data - data_smooth), 2, 'omitnan'));

%hold on; plot(data_ref, 'k', 'LineWidth', 2);

% remove the median behaviour from the observations
clean_ph = d2data - repmat(data_ref, 1, n_set);

% size of the win in seconds
med_win_size = 150;
spline_win_size = 180;
data_smooth = clean_ph;
for s = 1 : n_set
    valid_obs = ~isnan(data_smooth(:, s));
    if ~isnan(d2data_m(s))
        data_smooth(valid_obs, s) = medfilt_mat(data_smooth(valid_obs,s), ceil(med_win_size/state.getProcessingRate()));
        %data_smooth(valid_obs, s) = splinerMat(t(valid_obs)', data_smooth(valid_obs,s), ceil(spline_win_size/state.getProcessingRate()));
    end
end

%hold on; plot(data_smooth, 'k', 'LineWidth', 2);

clean_ph = clean_ph - data_smooth;
data_out(3 : end, :) = data_out(3 : end, :) | abs(clean_ph) > thr(3) * std(clean_ph(~isnan(clean_ph)));
data(data_out) = nan;
data_out = ~isnan(data) & ~remove_short_arcs(~isnan(data)', state.getMinArc())';
data(data_out) = nan;
%figure; plot(data, '.-'); set(gcf, 'Name', 'data clean');

data = nan2zero(data');
flag_array = (data ~= data_in);
