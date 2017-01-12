function [obs_out] = remove_short_arcs(obs_in, threshold_gap)

% SYNTAX:
%   [obs_out] = remove_short_arcs(obs_in, threshold_gap);
%
% INPUT:
%   obs_in = input observation matrix (num_sat x epochs)
%   threshold_gap = threshold on the number of epochs for an arc to be
%                   considered "short"
%
% OUTPUT:
%   obs_out = input observation matrix (num_sat x epochs)
%
% DESCRIPTION:
%   Removal of observation arcs shorter than given threshold.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

if (nargin < 2)
    threshold_gap = 10;
end

sat_conf_map = diff(obs_in ~= 0, 1, 2);

obs_out = obs_in;
n_epochs = size(obs_out,2);

%remove short arcs from beginning of observation timespan
[row_set, col_set] = find(sat_conf_map == -1);
idx = find(col_set <= threshold_gap);
epoch = col_set(idx);
sat = row_set(idx);
for s = 1 : length(sat)
    obs_out(sat(s), 1:epoch(s)) = 0;
end

%remove short arcs from end of observation timespan
[row_rise, col_rise] = find(sat_conf_map == 1);
idx = find(col_rise >= (n_epochs - threshold_gap));
epoch = col_rise(idx);
sat = row_rise(idx);
for s = 1 : length(sat)
    obs_out(sat(s), epoch(s):n_epochs) = 0;
end

%remove the remaining short arcs
sat = union(row_set, row_rise);
for s = 1 : length(sat)
    col = sort([col_set(row_set == sat(s));col_rise(row_rise == sat(s))]);
    dd = diff(col);
    idx = find(diff(col) <= threshold_gap);
    for i = 1 : length(idx)
        obs_out(sat(s),col(idx(i))+1:col(idx(i))+dd(idx(i))) = 0;
    end
end
