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

if (nargin < 2)
    threshold_gap = 10;
end

obs_out = obs_in;

for s = 1 : size(obs_in,1)

    % find intervals of zeros
    fI = getOutliers(obs_in(s,:) == 0 | isnan(obs_in(s, :)));

    % find the intervals of good obs_in
    vI = [[1; fI(:,2)+1] [fI(:,1)-1; size(obs_in,2)]];
    vI = [vI (vI(:,2) - vI(:,1) +1)];

    vI = vI(vI(:,3) <= threshold_gap, :);

    for i = 1 : size(vI,1)
        obs_out(s,vI(i,1) : vI(i,2)) = 0;
    end
end
