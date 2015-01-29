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

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
% Portion of code by Stefano Caldera
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


N_stim = [];

% range
comb_pr = (distR - distM) - (distR(pivot_index) - distM(pivot_index));

%observed phase double differences
comb_ph = (ph_R - ph_M) - (ph_R(pivot_index) - ph_M(pivot_index));

%linear combination of initial ambiguity estimate
N_stim = comb_pr ./ lambda - comb_ph;
%sigmaq_N_stim = sum(sigmaq_pos_R) ./ lambda.^2;


