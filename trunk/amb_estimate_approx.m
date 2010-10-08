function [N_stim, sigmaq_N_stim] = amb_estimate_approx(pos_R, pos_M, sigmaq_pos_R, ...
         pr_Rsat, pr_Msat, ph_Rsat, ph_Msat, Eph, time, pivot, sat, phase) %#ok<INUSL>

% SYNTAX:
%   [N_stim, sigmaq_N_stim] = amb_estimate_approx(pos_R, pos_M, sigmaq_pos_R, ...
%   pr_Rsat, pr_Msat, ph_Rsat, ph_Msat, Eph, time, pivot, sat, phase);
%
% INPUT:
%   pos_R = ROVER assessed position (X,Y,Z)
%   pos_M = MASTER assessed position (X,Y,Z)
%   sigmaq_pos_R = rounded ROVER position variance
%   pr_Rsat = ROVER-SATELLITE code-pseudorange
%   pr_Msat = MASTER-SATELLITE code-pseudorange (not used)
%   ph_Rsat = ROVER-SATELLITE phase-pseudorange
%   ph_Msat = MASTER-SATELLITE phase-pseudorange
%   Eph = ephemerides matrix
%   time = GPS time
%   pivot = PIVOT satellite
%   sat = configuration of satellites in view
%   phase = carrier L1 (phase=1), carrier L2 (phase=2)
%
% OUTPUT:
%   N_stim = linear combination of ambiguity estimate
%   sigmaq_N_stim = assessed variances of combined ambiguity
%
% DESCRIPTION:
%   Estimation of combined (double difference) phase ambiguities
%   (and of their error variance) by using both phase observations
%   and satellite-receiver distance, based on the ROVER
%   approximate position.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
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
global lambda1
global lambda2

%PIVOT position research
i = find(pivot == sat);

%PIVOT position correction (clock and rotation)
Rot_Xpivot = sat_corr(Eph, sat(i), time, pr_Rsat(i), pos_R);

%ROVER,MASTER-PIVOT approximate pseudorange estimate
pr_stim_RP = sqrt(sum((pos_R - Rot_Xpivot).^2));
pr_stim_MP = sqrt(sum((pos_M - Rot_Xpivot).^2));

%phase observations
ph_RP = ph_Rsat(i);
ph_MP = ph_Msat(i);

%loop on all used satellites
for m = 1 : size(sat,1)
    
    %new satellites position correction (clock and Earth rotation)
    Rot_X = sat_corr(Eph, sat(m), time, pr_Rsat(m), pos_R);
    
    %ROVER,MASTER-SATELLITES pseudorange estimate
    pr_stim_Rsat(m,1) = sqrt(sum((pos_R - Rot_X).^2));
    pr_stim_Msat(m,1) = sqrt(sum((pos_M - Rot_X).^2));
end

%estimated code pseudorange double differences
comb_pr = (pr_stim_Rsat - pr_stim_Msat) - (pr_stim_RP - pr_stim_MP);

%observed phase pseudorange double differences
comb_ph = (ph_Rsat - ph_Msat) - (ph_RP - ph_MP);

%linear combination of initial ambiguity estimate
if (phase == 1)
    N_stim = ((comb_pr - comb_ph * lambda1)) / lambda1;
    sigmaq_N_stim = sum(sigmaq_pos_R) / lambda1^2;
else
    N_stim = ((comb_pr - comb_ph * lambda2)) / lambda2;
    sigmaq_N_stim = sum(sigmaq_pos_R) / lambda2^2;
end