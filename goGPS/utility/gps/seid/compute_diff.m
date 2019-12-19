function [diff_L4, diff_P4, commontime, stations_idx, L1_series, L2_series, L4_series, P4_series] = compute_diff(L1_series, L2_series, P1_series, P2_series, name_series, time_series)

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

cc = Core.getCurrentSettings().getConstellationCollector();
lambda1 = cc.gps.L_VEC(1);
lambda2 = cc.gps.L_VEC(2);

%compute dL4
n_sta = length(name_series);

%choose reference time for interpolation
commontime = time_series{1};
for i = 2 : n_sta
    commontime = intersect(commontime,time_series{i});
end

%matrices same length as L4
stations_idx = zeros(n_sta, length(commontime));
L4_series = zeros(32, length(commontime), n_sta);
P4_series = zeros(32, length(commontime), n_sta);

%matrices same size as diff_L4
diff_L4 = zeros(32, length(commontime) - 1, n_sta);
diff_P4 = zeros(32, length(commontime) - 1, n_sta);

%compute L4 = L1 - L2, dL4, P4 = P1 - P2, dP4
for sta = 1 : n_sta
    [~, ~, stations_idx(sta,:)] = intersect(commontime, time_series{sta});
    L4_series(:,:,sta) = L1_series{sta} * lambda1 - L2_series{sta} * lambda2;
    P4_series(:,:,sta) = P2_series{sta}           - P1_series{sta};
    diff_L4(:,:,sta) = diff(L4_series(:,:,sta)')';
    diff_P4(:,:,sta) = diff(P4_series(:,:,sta)')';
end
