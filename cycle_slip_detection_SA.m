function [slip, N_slip, sat_slip] = cycle_slip_detection_SA(N_kalman, ...
         pr_Rsat, ph_Rsat, err_iono_RS, doppler_pred_range, sat, sat_born, alfa, phase)

% SYNTAX:
%   [slip, N_slip, sat_slip] = cycle_slip_detection_SA(N_kalman, ...
%    pr_Rsat, ph_Rsat, err_iono_RS, doppler_pred_range, sat, sat_born, alfa, phase);
%
% INPUT:
%   N_kalman = phase ambiguities (double difference) estimated by the Kalman filter
%   pr_Rsat = ROVER-SATELLITE code observation
%   ph_Rsat = ROVER-SATELLITE phase observation
%   err_iono_RS = ionospheric error
%   doppler_pred_range = predicted range based on phase and Doppler
%   observations from previous epoch
%   sat = visible satellites configuration
%   sat_born = new satellites (added in this epoch)
%   alfa = cycle-slip detection threshold
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% OUTPUT:
%   slip = boolean variable (slip=1 if there is a cycle-slip)
%   N_slip = phase ambiguity estimation after the cycle-slip
%   sat_slip = cycle-slip-affected satellite
%
% DESCRIPTION:
%   Detects the presence of cycle-slips and in case estimates
%   the new phase ambiguities. Estimation of the satellite-receiver
%   range on the basis of the Kalman filter.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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
global lambda1
global lambda2

%number of visible satellites
nsat = size(sat,1);

%phase ambiguities estimation
if (phase == 1)
    N_stim = (pr_Rsat - lambda1 * ph_Rsat - 2 * err_iono_RS) / lambda1;
else
    N_stim = (pr_Rsat - lambda2 * ph_Rsat - 2 * err_iono_RS) / lambda2;
end

%initialization
N_slip = [];
sat_slip = [];
slip = 0;

%cycle-slip detection
for i = 1 : nsat

    %test on the estimated value of the phase ambiguities
    if (~ismember(sat(i),sat_born))

        %test on:
        % - Kalman-estimated phase ambiguities compared with code-phase range (double differences)
        % - Doppler-predicted phase range compared to observed phase range
        if ((abs(N_kalman(sat(i)) - N_stim(i)) > 50) || (doppler_pred_range(i)) & (abs(doppler_pred_range(i) - ph_Rsat(i)) > alfa))

            %save of the new phase ambiguity estimation
            N_slip = [N_slip; N_stim(i)];

            %save of the slipped satellite
            sat_slip = [sat_slip; sat(i)];

            %flag identifying a cycle-slip
            slip = 1;
%             fprintf('Cycle-slip on satellite %d: range difference = %.3f\n', sat(i), doppler_pred_range(i) - ph_Rsat(i));        
%             pause
        end
    end
end
