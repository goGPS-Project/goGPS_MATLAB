function [slip, N_slip, sat_slip] = cycle_slip_detection_SA(N_kalman, ph, distR, dtS, dtR, err_tropo, err_iono, phwindup, doppler_pred_range, sat_pr, sat_ph, sat_born, alpha, lambda)

% SYNTAX:
%   [slip, N_slip, sat_slip] = cycle_slip_detection_SA(N_kalman, ...
%    pr, ph, err_iono, doppler_pred_range, sat, sat_born, alpha, lambda);
%
% INPUT:
%   N_kalman = phase ambiguities (double difference) estimated by the Kalman filter
%   ph = ROVER-SATELLITE phase observation
%   distR = ROVER-SATELLITE geometric range
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

%variable initialization
global flag_doppler_cs

v_light = Core_Utils.V_LIGHT;

%number of visible satellites
nsat = size(sat_ph,1);

%sat_ph is a subset of sat_pr
[~, index] = intersect(sat_pr,sat_ph);

%phase ambiguities estimation
rho = distR + sum(dtR) - v_light*dtS + err_tropo - err_iono + lambda.*phwindup;
N_stim = (rho(index) - lambda(index,1).*ph(index))./lambda(index,1);

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
    if (~ismember(sat_ph(i),sat_born))

        %Kalman-estimated phase ambiguities compared with ambiguities estimated by using the observed pseudorange
        if (~flag_doppler_cs & (abs(N_kalman(sat_ph(i)) - N_stim(i)) > alpha))
            cs = 1;
        end

        %Doppler-predicted phase range compared to observed phase range
        if (flag_doppler_cs & (abs(doppler_pred_range(i) - ph(i)) > alpha))
           cs = 1;
        end

        if (cs)
            %save the new phase ambiguity estimation
            N_slip = [N_slip; N_stim(i)];

            %save the slipped satellite
            sat_slip = [sat_slip; sat_ph(i)];

            %flag identifying a cycle-slip
            slip = 1;
        end
    end
end
