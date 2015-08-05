function [time_ref, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, max_int] = ...
          sync_obs(time_i, week_i, date_i, pr1_i, ph1_i, pr2_i, ph2_i, dop1_i, dop2_i, snr1_i, snr2_i, interval)

% SYNTAX:
%   [time_ref, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2] = ...
%   sync_obs(time_i, week_i, date_i, pr1_i, ph1_i, pr2_i, ph2_i, dop1_i, dop2_i, snr1_i, snr2_i, interval);
%
% INPUT:
%   time_i = receiver seconds-of-week
%   week_i = GPS week
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
%
% OUTPUT:
%   time_ref = reference seconds-of-week
%   time = receiver seconds-of-week
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

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
%----------------------------------------------------------------------------------------------
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
%----------------------------------------------------------------------------------------------

%number of satellite slots
nSatTot = size(pr1_i,1);

%number of observation datasets (e.g. number of read RINEX files)
nObsSet = size(pr1_i,3);

%find min and max time tags (in common among all observation datasets)
time_i_nan = time_i;
time_i_nan(permute(sum(pr1_i,1),[2 1 3])==0) = NaN; %set NaN to epochs which don't have any pseudorange
min_time = max(min(time_i_nan,[],1));
min_time_prog = min(min(time_i_nan,[],1));
max_time = min(max(time_i_nan,[],1));

%find the largest interval
max_int = max(interval(:));
%max_int = 30;

%define the reference time
time_ref = (roundmod(min_time,max_int) : max_int : roundmod(max_time,max_int))';
tow_ref = mod(time_ref,60*60*24*7);
tow_ref=roundmod(tow_ref,max_int);

%number of reference epochs
ref_len = length(tow_ref);

%create containers
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

time_prog = time_i - min_time_prog; % substract the first element to reduce the magnitude of all the values
time_ref_prog = time_ref - min_time_prog;

for s = 1 : nObsSet
%    [~, idx_t, idx_z] = intersect(roundmod(time_ref_prog, max_int), roundmod(time_prog(:,1,s), max_int));
    [~, idx_t, idx_z] = intersect(roundmod(time_ref_prog, max_int), roundmod(time_prog(:,1,s), interval(s)));    
    time(idx_t, s) = time_i(idx_z, 1, s);
    week(idx_t, s) = week_i(idx_z, 1, s);
    date(idx_t, :, s) = date_i(idx_z, :, s);
    pr1(:,  idx_t, s) = pr1_i(:, idx_z, s);
    ph1(:,  idx_t, s) = ph1_i(:, idx_z, s);
    pr2(:,  idx_t, s) = pr2_i(:, idx_z, s);
    ph2(:,  idx_t, s) = ph2_i(:, idx_z, s);
    dop1(:, idx_t, s) = dop1_i(:, idx_z, s);
    dop2(:, idx_t, s) = dop2_i(:, idx_z, s);
    snr1(:, idx_t, s) = snr1_i(:, idx_z, s);
    snr2(:, idx_t, s) = snr2_i(:, idx_z, s);
end
