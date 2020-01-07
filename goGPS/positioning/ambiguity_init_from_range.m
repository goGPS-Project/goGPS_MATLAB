function [N_stim] = ambiguity_init_from_range(distR, distM, pivot_index, ph_R, ph_M, lambda)

% SYNTAX:
%   [N_stim] = ambiguity_init_from_range(distR, distM, pivot_index, ph_R, ph_R)
%
% INPUT:
%   distR = ROVER-SATELLITE range
%   distM = MASTER-SATELLITE range
%   pivot_index = pivot satellite index
%   XR_approx = receiver approximate position (X,Y,Z)
%   ph_R = ROVER-SATELLITE phase measurement
%   ph_M = MASTER-SATELLITE phase measurement
%   lambda = vector containing GNSS wavelengths for available satellites
%
%
% OUTPUT:
%   N_stim = phase ambiguity estimation
%
% DESCRIPTION:
%   This function estimates float phase ambiguities from geometrical range
%   and phase observations

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     Stefano Caldera, ...
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
% range
comb_pr = (distR - distM) - (distR(pivot_index) - distM(pivot_index));

%observed phase double differences
comb_ph = (ph_R - ph_M) - (ph_R(pivot_index) - ph_M(pivot_index));

%linear combination of initial ambiguity estimate
N_stim = comb_pr ./ lambda - comb_ph;
%sigmaq_N_stim = sum(sigmaq_pos_R) ./ lambda.^2;


