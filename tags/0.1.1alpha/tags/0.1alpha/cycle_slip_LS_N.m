function [slip, comb_N_slip, sat_slip, sigmaq_comb_N_slip] = cycle_slip_LS_N(posM, posR, comb_N_kalman, ...
         ph_Rsat, ph_Msat, pr_Rsat, pr_Msat, snr_R, snr_M, Eph, time, pivot, sat, iono, alpha, phase)

% SYNTAX:
%   [slip, comb_N_slip, sat_slip, sigmaq_comb_N_slip] = cycle_slip_LS_N(posM, posR, comb_N_kalman, ...
%   ph_Rsat, ph_Msat, pr_Rsat, pr_Msat, snr_R, snr_M, Eph, time, pivot, sat, iono, alpha, phase);
%
% INPUT:
%   posM = MASTER position (X,Y,Z)
%   posR = ROVER position (X,Y,Z) estimated by the Kalman filter
%   comb_N_kalman = phase ambiguities (double difference) estimated by the Kalman filter
%   ph_Rsat = ROVER-SATELLITE phase observation
%   ph_Msat = MASTER-SATELLITE phase observation
%   pr_Rsat = ROVER-SATELLITE code pseudorange
%   pr_Msat = MASTER-SATELLITE code pseudorange
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   Eph = ephemerides matrix
%   time = GPS time
%   pivot = PIVOT satellite
%   sat = visible satellites configuration
%   iono = ionosphere parameters
%   alpha = cycle-slip detection threshold
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% OUTPUT:
%   slip = boolean variable (slip=1 if there is a cycle-slip)
%   comb_N_slip = phase ambiguity estimation after the cycle-slip
%   sat_slip = cycle-slip-affected satellite
%   sigmaq_comb_N_slip = variance of the new ambiguity estimates
%
% DESCRIPTION:
%   Detects the presence of cycle-slips and in case estimates
%   the new phase ambiguities by least squares method.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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

%variable initialization
% global lambda1 lambda2
global sigmaq0_N

%number of visible satellites
nsat = size(sat,1);

%N combination estimation (least squares)
if (phase == 1)
    [null_pos_R, null_cov_pos_R, comb_N_stim, cov_comb_N_stim] = code_phase_double_diff(posR, pr_Rsat, ph_Rsat, snr_R, posM, pr_Msat, ph_Msat, snr_M, time, sat, pivot, Eph, 1, iono); %#ok<ASGLU>
else
    [null_pos_R, null_cov_pos_R, comb_N_stim, cov_comb_N_stim] = code_phase_double_diff(posR, pr_Rsat, ph_Rsat, snr_R, posM, pr_Msat, ph_Msat, snr_M, time, sat, pivot, Eph, 2, iono); %#ok<ASGLU>
end

if isempty(cov_comb_N_stim) %if it was not possible to compute the covariance matrix
    cov_comb_N_stim = sigmaq0_N * eye(nsat);
end
sigmaq_comb_N = diag(cov_comb_N_stim);

%initialization
comb_N_slip = [];
sat_slip = [];
sigmaq_comb_N_slip = [];
slip = 0;

%cycle-slip detection
for i = 1 : nsat
    %test on the estimated value of the phase ambiguities
    if (sat(i) ~= pivot)

        %alpha = max(3*Cee(o3+sat(i),o3+sat(i)),3*sigmaq_comb_N(i));

        if (abs(comb_N_kalman(sat(i)) - comb_N_stim(i)) > alpha)

            %save the new phase ambiguity estimation
            comb_N_slip = [comb_N_slip; comb_N_stim(i)];

            %save the slipped satellite
            sat_slip = [sat_slip; sat(i)];

            %save the variance of the new phase ambiguity estimation
            sigmaq_comb_N_slip = [sigmaq_comb_N_slip; sigmaq_comb_N(i)];

            %flag identifying a cycle-slip
            slip = 1;
        end
    end
end
