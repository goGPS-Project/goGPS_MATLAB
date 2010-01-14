function [slip, N_slip, sat_slip] = cycle_slip_kalman(posM, posR, N_kalman, ...
         ph_Rsat, ph_Msat, pr_Rsat, pr_Msat, Eph, time, pivot, sat, alfa, phase)

% SYNTAX:
%   [slip, N_slip, sat_slip] = cycle_slip_kalman(posM, posR, N_kalman, ...
%   ph_Rsat, ph_Msat, pr_Rsat, pr_Msat, Eph, time, pivot, sat, phase);
%
% INPUT:
%   posM = MASTER position (X,Y,Z)
%   posR = ROVER position (X,Y,Z) estimated by the Kalman filter
%   N_kalman = phase ambiguities (double difference) estimated by the Kalman filter
%   ph_Rsat = ROVER-SATELLITE phase observation
%   ph_Msat = MASTER-SATELLITE phase observation
%   pr_Rsat = ROVER-SATELLITE code pseudorange
%   pr_Msat = MASTER-SATELLITE code pseudorange
%   Eph = ephemerides matrix
%   time = GPS time
%   pivot = PIVOT satellite
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
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
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

%PIVOT search
i = find(pivot == sat);

%PIVOT position (with clock error and Earth rotation corrections)
posP = sat_corr(Eph, sat(i), time, pr_Rsat(i), posR);

%estimation of ROVER-PIVOT and MASTER-PIVOT pseudoranges
pr_RP_stim = sqrt(sum((posR - posP).^2));
pr_MP_stim = sqrt(sum((posM - posP).^2));

%phase observations
ph_RP = ph_Rsat(i);
ph_MP = ph_Msat(i);

%initialization
pr_stim = [];

%computation for all the satellites, PIVOT included
for i = 1 : nsat

    %satellite position (with clock error and Earth rotation corrections)
    posS = sat_corr(Eph, sat(i), time, pr_Rsat(i), posR);

    %estimation of ROVER-PIVOT and MASTER-PIVOT pseudoranges
    pr_RS_stim = sqrt(sum((posR - posS).^2));
    pr_MS_stim = sqrt(sum((posM - posS).^2));

    %computation of crossed pseudoranges
    pr_stim = [pr_stim; (pr_RS_stim - pr_MS_stim) - (pr_RP_stim - pr_MP_stim)];
end

%observed phase double difference
ddp = (ph_Rsat - ph_Msat) - (ph_RP - ph_MP);

%phase ambiguities estimation
if (phase == 1)
    N_stim = pr_stim / lambda1 - ddp;
else
    N_stim = pr_stim / lambda2 - ddp;
end

%initialization
N_slip = [];
sat_slip = [];
slip = 0;

%cycle-slip detection
for i = 1 : nsat

    %test on the estimated value of the phase ambiguities
    if (sat(i) ~= pivot) & (abs(N_kalman(sat(i)) - N_stim(i)) > alfa)

        %cycle-slip detected
        %['Detected CYCLE-SLIP for the satellite ' num2str(sat(i)) ' at epoch ' num2str(time)];

        %save of the new phase ambiguity estimation
        N_slip = [N_slip; N_stim(i)];

        %save of the slipped satellite
        sat_slip = [sat_slip; sat(i)];

        %flag identifying a cycle-slip
        slip = 1;
    end
end

%-------------------------------------------------------------------------------