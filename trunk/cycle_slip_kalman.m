function [slip, comb_N_slip, sat_slip, sigmaq_comb_N_slip] = cycle_slip_kalman(posM, posR, comb_N_kalman, ...
         ph_Rsat, ph_Msat, pr_Rsat, pr_Msat, snr_R, snr_M, Eph, time, pivot, sat, iono, alpha, phase)

% SYNTAX:
%   [slip, comb_N_slip, sat_slip, sigmaq_comb_N_slip] = cycle_slip_kalman(posM, posR, comb_N_kalman, ...
%   ph_Rsat, ph_Msat, pr_Rsat, pr_Msat, Eph, time, pivot, sat, phase);
%
% INPUT:
%   posM = MASTER position (X,Y,Z)
%   posR = ROVER position (X,Y,Z) estimated by the Kalman filter
%   comb_N_kalman = phase ambiguities (double difference) estimated by the Kalman filter
%   ph_Rsat = ROVER-SATELLITE phase observation
%   ph_Msat = MASTER-SATELLITE phase observation
%   pr_Rsat = ROVER-SATELLITE code pseudorange
%   pr_Msat = MASTER-SATELLITE code pseudorange
%   Eph = ephemerides matrix
%   time = GPS time
%   pivot = PIVOT satellite
%   sat = visible satellites configuration
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
%   the new phase ambiguities. Estimation of the satellite-receiver
%   range on the basis of the Kalman filter.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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
global o3 Cee

% %number of unknown phase ambiguities 
% if (length(phase) == 1)
%     nN = 32;
% else
%     nN = 64;
% end

%number of visible satellites
nsat = size(sat,1);

% %PIVOT search
% i = find(pivot == sat);

% %PIVOT position (with clock error and Earth rotation corrections)
% posP = sat_corr(Eph, sat(i), time, pr_Rsat(i), posR);

% %estimation of ROVER-PIVOT and MASTER-PIVOT pseudoranges
% pr_RP_app = sqrt(sum((posR - posP).^2));
% pr_MP_app = sqrt(sum((posM - posP).^2));
% 
% %phase observations
% ph_RP = ph_Rsat(i);
% ph_MP = ph_Msat(i);

% %initialization
% pr_app = [];
% 
% %computation for all the satellites, PIVOT included
% for i = 1 : nsat
% 
%     %satellite position (with clock error and Earth rotation corrections)
%     posS = sat_corr(Eph, sat(i), time, pr_Rsat(i), posR);
% 
%     %estimation of ROVER-SATELLITE and MASTER-SATELLITE pseudoranges
%     pr_RS_app = sqrt(sum((posR - posS).^2));
%     pr_MS_app = sqrt(sum((posM - posS).^2));
% 
%     %computation of crossed pseudoranges
%     pr_app = [pr_app; (pr_RS_app - pr_MS_app) - (pr_RP_app - pr_MP_app)];
% end
% 
% %observed phase double difference
% ddp = (ph_Rsat - ph_Msat) - (ph_RP - ph_MP);
% 
% %phase ambiguities estimation
% if (phase == 1)
%     comb_N_app = pr_app / lambda1 - ddp;
% else
%     comb_N_app = pr_app / lambda2 - ddp;
% end

%N combination estimation (least squares)
if (phase == 1)
    [null_pos_R, null_cov_pos_R, comb_N_stim, cov_comb_N_stim] = code_phase_double_diff(posR, pr_Rsat, ph_Rsat, snr_R, posM, pr_Msat, ph_Msat, snr_M, time, sat, pivot, Eph, iono, 1);
else
    [null_pos_R, null_cov_pos_R, comb_N_stim, cov_comb_N_stim] = code_phase_double_diff(posR, pr_Rsat, ph_Rsat, snr_R, posM, pr_Msat, ph_Msat, snr_M, time, sat, pivot, Eph, iono, 2);
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

%-------------------------------------------------------------------------------