function [N_stim, sigmaq_N_stim] = amb_estimate_observ_SA(pr_Rsat, ph_Rsat, lambda)

% SYNTAX:
%   [N_stim, sigmaq_N_stim] = amb_estimate_observ_SA(pr_Rsat, ph_Rsat, lambda);
%
% INPUT:
%   pr_Rsat = ROVER-SATELLITE code-pseudorange
%   ph_Rsat = ROVER-SATELLITE phase-pseudorange
%   lambda = vector containing GNSS wavelengths for available satellites
%
% OUTPUT:
%   N_stim = linear combination ambiguity estimate
%   sigmaq_N_stim = assessed variances of combined ambiguity
%
% DESCRIPTION:
%   Estimation of phase ambiguities (and of their error variance) by
%   using both phase and code observations (satellite-receiver distance) in
%   undifferenced mode.

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

%variables initialization
global sigmaq_cod1

N_stim = ((pr_Rsat - ph_Rsat .* lambda)) ./ lambda;
sigmaq_N_stim = sigmaq_cod1 ./ lambda.^2;
