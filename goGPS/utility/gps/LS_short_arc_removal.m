function [A, y0, b, Q, obs_track, amb_num, amb_prn_track, rem_amb] = LS_short_arc_removal(A, y0, b, Q, obs_track, amb_num, amb_prn_track, min_arc)
% Remove ambiguity unkowns with arcs shorter than given threshold
% SYNTAX: [A, y0, b, Q, obs_track, amb_num, amb_prn_track, rem_amb] = LS_short_arc_removal(A, y0, b, Q, obs_track, [], amb_prn_track, min_arc)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Eugenio Realini
%  Contributors:     Eugenio Realini, Andrea Gatti, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can REDistribute it and/or modify
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

% %remove observations without pivot OR slave satellite
% rem_obs = find(sum(A(:,4:end),2));
% rem_obs = setdiff(rem_obs, ok_obs);
% if (~isempty(rem_obs))
%     A(rem_obs,:) = [];
%     y0(rem_obs) = [];
%     b(rem_obs) = [];
%     Q(rem_obs,:) = []; Q(:,rem_obs) = [];
%     obs_track(rem_obs,:) = [];
%     ok_obs = find(sum(A(:,4:end),2));
% end

amb_num = numel(amb_prn_track);
pos_num = size(A,2) - amb_num;
rem_amb = setdiff(find(sum(A~=0,1) < min_arc), 1 : pos_num);

if (~isempty(rem_amb))
    rem_obs = [];
    for r = 1 : length(rem_amb)
        rem_obs = [rem_obs; find(A(:,rem_amb(r))~=0)]; %#ok<AGROW>
    end
    A(rem_obs,:) = [];
    y0(rem_obs) = [];
    b(rem_obs) = [];
    Q(rem_obs,:) = []; Q(:,rem_obs) = [];
    obs_track(rem_obs,:) = [];

    A(:,rem_amb) = [];
    amb_num = amb_num - length(rem_amb);
    amb_prn_track(rem_amb - pos_num) = [];

%     %check again and remove observations without pivot OR slave satellite
%     rem_obs = find(sum(A(:,4:end),2));
%     A(rem_obs,:) = [];
%     y0(rem_obs) = [];
%     b(rem_obs) = [];
%     Q(rem_obs,:) = []; Q(:,rem_obs) = [];
%     obs_track(rem_obs,:) = [];
%
%     rem_amb = find(sum(A~=0,1) < min_arc);
end
