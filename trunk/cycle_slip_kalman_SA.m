function [slip, N_slip, sat_slip] = cycle_slip_kalman_SA(posR, N_kalman, ...
         ph_Rsat, pr_Rsat, Eph, time, sat, alfa, phase)

% SYNTAX:
%   [slip, N_slip, sat_slip] = cycle_slip_kalman_SA(posR, N_kalman, ...
%   ph_Rsat, pr_Rsat, Eph, time, sat, phase);
%
% INPUT:
%   posR = ROVER position (X,Y,Z) estimated by the Kalman filter
%   N_kalman = phase ambiguities (double difference) estimated by the Kalman filter
%   ph_Rsat = ROVER-SATELLITE phase observation
%   pr_Rsat = ROVER-SATELLITE code pseudorange
%   Eph = ephemerides matrix
%   time = GPS time
%   sat = visible satellites configuration
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
%                           goGPS v0.1 beta
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
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
global lambda1;
global lambda2;

%number of visible satellites
nsat = size(sat,1);

%initialization
pr_stim = [];

%computation for all the satellites, PIVOT included
for i = 1 : nsat

    %satellite position (with clock error and Earth rotation corrections)
    posS = sat_corr(Eph, sat(i), time, pr_Rsat(i), posR);

    %estimation of ROVER-PIVOT and MASTER-PIVOT pseudoranges
    pr_RS_stim = sqrt(sum((posR - posS).^2));

    %store pseudoranges
    pr_stim = [pr_stim; pr_RS_stim];
end

%phase ambiguities estimation
if (phase == 1)
    N_stim = pr_stim / lambda1 - ph_Rsat;
else
    N_stim = pr_stim / lambda2 - ph_Rsat;
end

%initialization
N_slip = [];
sat_slip = [];
slip = 0;

%cycle-slip detection
for i = 1 : nsat

    %test on the estimated value of the phase ambiguities
    if (abs(N_kalman(sat(i)) - N_stim(i)) > alfa)

        %save of the new phase ambiguity estimation
        N_slip = [N_slip; N_stim(i)];

        %save of the slipped satellite
        sat_slip = [sat_slip; sat(i)];

        %flag identifying a cycle-slip
        slip = 1;
    end
end
