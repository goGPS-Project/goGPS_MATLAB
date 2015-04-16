function [slip, N_slip, sat_slip] = cycle_slip_detection(N_kalman, ...
         ph_Rsat, ph_Msat, pr_RSapp, pr_MSapp, doppler_pred_range_R, doppler_pred_range_M, pivot, sat, sat_born, alpha, lambda)

% SYNTAX:
%   [slip, N_slip, sat_slip] = cycle_slip_detection(N_kalman, ...
%   ph_Rsat, ph_Msat, prRS_app, prMS_app, doppler_pred_range_R, doppler_pred_range_M, pivot, sat, sat_born, alpha, lambda);
%
% INPUT:
%   N_kalman = phase ambiguities (double difference) estimated by the Kalman filter
%   ph_Rsat = ROVER-SATELLITE phase observation
%   ph_Msat = MASTER-SATELLITE phase observation
%   pr_Rsat = ROVER-SATELLITE code pseudorange observation
%   pr_Msat = MASTER-SATELLITE code pseudorange observation
%   prRS_app = ROVER-SATELLITE approximate pseudorange
%   prMS_app = MASTER-SATELLITE approximate pseudorange
%   doppler_pred_range_R = ROVER-SATELLITE Doppler-based predicted range
%   doppler_pred_range_M = MASTER-SATELLITE Doppler-based predicted range
%   pivot = PIVOT satellite
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
%   the new phase ambiguities.

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
global flag_doppler_cs

%number of visible satellites
nsat = size(sat,1);

%PIVOT search
i = find(pivot == sat);

%ROVER-PIVOT and MASTER-PIVOT approximate pseudoranges
pr_RPapp = pr_RSapp(i);
pr_MPapp = pr_MSapp(i);

%phase observations
ph_RP = ph_Rsat(i);
ph_MP = ph_Msat(i);

%approximate code double differences
comb_pr = (pr_RSapp - pr_MSapp) - (pr_RPapp - pr_MPapp);

%observed phase double differences
comb_ph = (ph_Rsat - ph_Msat) - (ph_RP - ph_MP);

%phase ambiguities estimation
N_stim = comb_pr ./ lambda - comb_ph;

%initialization
N_slip = [];
sat_slip = [];
slip = 0;

%disable Doppler-based cycle-slip detection if Doppler observations are not available either for the rover or the master
if (flag_doppler_cs & doppler_pred_range_R == 0)
    fprintf('Rover Doppler observations are not available: disabling Doppler-based cycle-slip detection.\n');
    flag_doppler_cs = 0;
end
if (flag_doppler_cs & doppler_pred_range_M == 0)
    fprintf('Master Doppler observations are not available: disabling Doppler-based cycle-slip detection.\n');
    flag_doppler_cs = 0;
end

%cycle-slip detection
for i = 1 : nsat
    
    cs = 0;

    if (~ismember(sat(i),sat_born))

        %Kalman-estimated phase ambiguities compared with ambiguities estimated by using approximate pseudorange (double differences)
        if (~flag_doppler_cs & (abs(N_kalman(sat(i)) - N_stim(i)) > alpha) & (sat(i) ~= pivot))
            cs = 1;
        end

        %Doppler-predicted phase range compared to observed phase range (ROVER and MASTER)
        if (flag_doppler_cs & (abs(doppler_pred_range_R(i) - ph_Rsat(i)) > alpha | ...
                               abs(doppler_pred_range_M(i) - ph_Msat(i)) > alpha))
           cs = 1;

           %if the pivot slips, all ambiguity combinations must be re-estimated
           if (sat(i) == pivot)
               N_slip = N_stim;
               sat_slip = sat;
               sat_slip(sat_slip == pivot) = [];
               sat_slip = setdiff(sat_slip,sat_born);
               slip = 1;
               return
           end
        end

        if (cs)
            %new phase ambiguity
            N_slip = [N_slip; N_stim(i)];

            %slipped satellite
            sat_slip = [sat_slip; sat(i)];

            %flag identifying one or more cycle-slips
            slip = 1;
        end
    end
end

%-------------------------------------------------------------------------------