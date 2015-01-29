function [slip, N_slip, sat_slip] = cycle_slip_detection_SA(N_kalman, ...
         pr, ph, err_iono, doppler_pred_range, sat, sat_born, alpha, lambda)

% SYNTAX:
%   [slip, N_slip, sat_slip] = cycle_slip_detection_SA(N_kalman, ...
%    pr, ph, err_iono, doppler_pred_range, sat, sat_born, alpha, lambda);
%
% INPUT:
%   N_kalman = phase ambiguities (double difference) estimated by the Kalman filter
%   pr = ROVER-SATELLITE code observation
%   ph = ROVER-SATELLITE phase observation
%   err_iono = ionospheric error
%   doppler_pred_range = predicted range based on phase and Doppler observations from previous epoch
%   sat = visible satellites configuration
%   sat_born = new satellites (added in this epoch)
%   alpha = cycle-slip detection threshold
%   lambda = vector containing GNSS wavelengths for available satellites
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

%variable initialization
global flag_doppler_cs

%number of visible satellites
nsat = size(sat,1);

%phase ambiguities estimation
N_stim = (pr - lambda .* ph - 2 * err_iono) ./ lambda;

%initialization
N_slip = [];
sat_slip = [];
slip = 0;

%disable Doppler-based cycle-slip detection if Doppler observations are not available either for the rover or the master
if (flag_doppler_cs & doppler_pred_range == 0)
    fprintf('Rover Doppler observations are not available: disabling Doppler-based cycle-slip detection.\n');
    flag_doppler_cs = 0;
end

%cycle-slip detection
for i = 1 : nsat
    
    cs = 0;

    %test on the estimated value of the phase ambiguities
    if (~ismember(sat(i),sat_born))

        %Kalman-estimated phase ambiguities compared with ambiguities estimated by using the observed pseudorange
        if (~flag_doppler_cs & (abs(N_kalman(sat(i)) - N_stim(i)) > alpha))
            cs = 1;
        end

        %Doppler-predicted phase range compared to observed phase range
        if (flag_doppler_cs & (abs(doppler_pred_range(i) - ph(i)) > alpha))
           cs = 1;
        end

        if (cs)
            %save of the new phase ambiguity estimation
            N_slip = [N_slip; N_stim(i)];

            %save of the slipped satellite
            sat_slip = [sat_slip; sat(i)];

            %flag identifying a cycle-slip
            slip = 1;
        end
    end
end
