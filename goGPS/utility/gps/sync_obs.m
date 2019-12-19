function [time_zero, time_GPS, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, codeC1, max_int] = ...
          sync_obs(time_i, date_i, pr1_i, ph1_i, pr2_i, ph2_i, dop1_i, dop2_i, snr1_i, snr2_i, codeC1_i, interval, processing_interval)
% SYNTAX:
%   [time_GPS, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, codeC1, max_int] = ...
%   sync_obs(time_i, date_i, pr1_i, ph1_i, pr2_i, ph2_i, dop1_i, dop2_i, snr1_i, snr2_i, codeC1_i, interval, processing_interval);
%
% INPUT:
%   time_i = receiver GPS_Time
%   date_i = date (year,month,day,hour,minute,second)
%   pr1_i = code observation (L1 carrier)
%   ph1_i = phase observation (L1 carrier)
%   pr2_i = code observation (L2 carrier)
%   ph2_i = phase observation (L2 carrier)
%   dop1_i = Doppler observation (L1 carrier)
%   dop2_i = Doppler observation (L2 carrier)
%   snr1_i = signal-to-noise ratio (L1 carrier)
%   snr2_i = signal-to-noise ratio (L2 carrier)
%   interval = observation time interval [s]
%   processing_interval = processing time interval [s]
%
% OUTPUT:
%   time_zero = reference epoch as seconds from 6 Jan 1980
%   time_GPS = nominal seconds from time_zero
%   time = receiver seconds from time_zero
%   week = GPS week
%   date = date (year,month,day,hour,minute,second)
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   dop1 = Doppler observation (L1 carrier)
%   dop2 = Doppler observation (L2 carrier)
%   snr1 = signal-to-noise ratio (L1 carrier)
%   snr2 = signal-to-noise ratio (L2 carrier)
%
% DESCRIPTION:
%   Synchronize different sets of observations. Zeros where not available.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
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

% number of satellite slots
nSatTot = size(pr1_i,1);

% number of observation datasets (e.g. number of read RINEX files)
nObsSet = size(pr1_i,3);

% logical indexes of all the receivers [n_epochs x 1 x n_rec] where pr is valid
id_ok = (permute(sum(pr1_i,1),[2 1 3]) > 0);

min_time = -inf;
max_time = inf;
for s = 1 : nObsSet
    min_time = max(min_time, time_i(s).first(squeeze(id_ok(:,1,s))).getGpsTime());
    max_time = min(max_time, time_i(s).last(squeeze(id_ok(:,1,s))).getGpsTime());
end
clear id_ok;

% find the largest interval
max_int = max(interval(:));
if (processing_interval > max_int)
    max_int = processing_interval;
end

% define the reference time
time_zero = ceil(min_time/max_int)*max_int;
max_time_diff = round((floor(max_time./max_int).*max_int - time_zero) / max_int);
time_GPS = (0 : 1 : max_time_diff)' * max_int;

% -1 needed for preprocessing, there is a bug somewere in there
time_zero = time_zero - 1;
time_GPS = time_GPS + 1;

% number of reference epochs
ref_len = length(time_GPS);

% create containers
time = zeros(ref_len, 1, nObsSet);
week = zeros(ref_len, 1, nObsSet);
date = zeros(ref_len, 6, nObsSet);
pr1  = zeros(nSatTot, ref_len, nObsSet);
ph1  = zeros(nSatTot, ref_len, nObsSet);
pr2  = zeros(nSatTot, ref_len, nObsSet);
ph2  = zeros(nSatTot, ref_len, nObsSet);
dop1 = zeros(nSatTot, ref_len, nObsSet);
dop2 = zeros(nSatTot, ref_len, nObsSet);
snr1 = zeros(nSatTot, ref_len, nObsSet);
snr2 = zeros(nSatTot, ref_len, nObsSet);
codeC1 = zeros(nSatTot, ref_len, nObsSet);

for s = 1 : nObsSet
    time_prog = time_i(s).getGpsTime(time_zero); % substract the first element to reduce the magnitude of all the values
    [~, idx_t, idx_z] = intersect(roundmod(time_GPS - 1, max_int), roundmod(time_prog - 1, interval(s)));
    tmp_time = time_i(s).getEpoch(idx_z);
    time(idx_t, s) = tmp_time.getGpsTime(time_zero);
    week(idx_t, s) = tmp_time.getGpsWeek();
    date(idx_t, :, s) = date_i(idx_z, :, s);
    pr1(:,  idx_t, s) = pr1_i(:, idx_z, s);
    ph1(:,  idx_t, s) = ph1_i(:, idx_z, s);
    pr2(:,  idx_t, s) = pr2_i(:, idx_z, s);
    ph2(:,  idx_t, s) = ph2_i(:, idx_z, s);
    dop1(:, idx_t, s) = dop1_i(:, idx_z, s);
    dop2(:, idx_t, s) = dop2_i(:, idx_z, s);
    snr1(:, idx_t, s) = snr1_i(:, idx_z, s);
    snr2(:, idx_t, s) = snr2_i(:, idx_z, s);
    codeC1(:, idx_t, s) = codeC1_i(:, idx_z, s);
end
