function [pr_out, ph_out] = correct_time_desync(time_ref, time_in, pr_in, ph_in, lambda)

% SYNTAX:
%   [pr_out, ph_out] = correct_time_desync(time_ref, time_in, pr_in, ph_in, lambda);
%
% INPUT:
%   obs_in =
%
% OUTPUT:
%   obs_out =
%
% DESCRIPTION:
%   Correction of nominal time (e.g. RINEX time tag) desynchronization.

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

pr_out = pr_in;
ph_out = ph_in;

n_sat = size(pr_in,1);
n_epochs = size(pr_in,2);

time_desync = time_ref - time_in;

time_desync = time_desync(:,ones(n_sat,1))';
lambda = lambda(:,ones(n_epochs,1));

idx_pr = pr_in ~= 0;
idx_ph = ph_in ~= 0;
pr_out(idx_pr) = pr_in(idx_pr) + time_desync(idx_pr)*Core_Utils.V_LIGHT;
ph_out(idx_ph) = ph_in(idx_ph) + time_desync(idx_ph)*Core_Utils.V_LIGHT./lambda(idx_ph);
