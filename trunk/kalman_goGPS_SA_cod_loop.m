function [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_SA_cod_loop(time, Eph, iono, pr1_Rsat, pr2_Rsat, snr_R, phase)

% SYNTAX:
%   kalman_goGPS_SA_cod_loop(time, Eph, iono, pr1_Rsat, pr2_Rsat, snr_R, phase);
%
% INPUT:
%   time = GPS time
%   Eph = satellite ephemerides
%   iono = ionospheric parameters
%   pr1_Rsat = ROVER-SATELLITE code pseudorange (L1 carrier)
%   pr2_Rsat = ROVER-SATELLITE code pseudorange (L2 carrier)
%   snr_R = signal-to-noise ratio for ROVER observations
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% OUTPUT:
%   check_on = boolean variable for satellite addition
%   check_off = boolean variable for satellite loss
%   check_pivot = boolean variable for pivot change
%   check_cs = boolean variable for cycle-slip
%
% DESCRIPTION:
%   Kalman filter for the ROVER trajectory computation.
%   Standalone positioning using code.

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

global sigmaq_velx sigmaq_vely sigmaq_velz
global min_nsat cutoff snr_threshold o1 o2 o3

global Xhat_t_t X_t1_t T I Cee nsat conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM

%----------------------------------------------------------------------------------------
% INITIALIZATION
%----------------------------------------------------------------------------------------

%output variables to point out events (satellite addition, losses, etc)
check_on = 0;
check_off = 0;
check_pivot = 0;
check_cs = 0;

%azimuth, elevation and ROVER-satellite and MASTER-satellite distances
azR = zeros(32,1);
azM = zeros(32,1);
elR = zeros(32,1);
elM = zeros(32,1);
distR = zeros(32,1);
distM = zeros(32,1);

%----------------------------------------------------------------------------------------
% MODEL ERROR COVARIANCE MATRIX
%----------------------------------------------------------------------------------------

%re-initialization of Cvv matrix of the model error
Cvv = zeros(o3);
Cvv(o1,o1) = sigmaq_velx;
Cvv(o2,o2) = sigmaq_vely;
Cvv(o3,o3) = sigmaq_velz;

%------------------------------------------------------------------------------------
% SATELLITE SELECTION
%------------------------------------------------------------------------------------

%visible satellites
if (length(phase) == 2)
    sat = find( (pr1_Rsat ~= 0) & (pr2_Rsat ~= 0) );
else
    if (phase == 1)
        sat = find( pr1_Rsat ~= 0 );
    else
        sat = find( pr2_Rsat ~= 0 );
    end
end

%------------------------------------------------------------------------------------
% SATELLITE ELEVATION, PIVOT AND CUT-OFF
%------------------------------------------------------------------------------------

j = 1;
bad_sat = [];

for i = 1:size(sat)

    %satellite position correction (clock and Earth rotation)
    Rot_X = sat_corr(Eph, sat(i), time, pr1_Rsat(i), X_t1_t([1,o1+1,o2+1])');

    %azimuth, elevation, ROVER-SATELLITE distance computation
    [azR(sat(i)), elR(sat(i)), distR(sat(i))] = topocent(X_t1_t([1,o1+1,o2+1]), Rot_X');

    %test on elevation and on signal-to-noise ratio
    if (elR(sat(i)) < cutoff) | (snr_R(sat(i)) < snr_threshold)
        bad_sat(j,1) = sat(i);
        j = j + 1;
    end
end

%removal of satellites with elevation or SNR lower than the respective threshold
sat(ismember(sat,bad_sat) == 1) = [];

%previous pivot
if (pivot ~= 0)
    pivot_old = pivot;
end

%current pivot
[null_max_elR, i] = max(elR(sat)); %#ok<ASGLU>
pivot = sat(i);

%----------------------------------------------------------------------------------------
% SATELLITE CONFIGURATION
%----------------------------------------------------------------------------------------

%previous satellite configuration (with phase measurements)
sat_old = find(conf_sat == 1);

%satellite configuration
conf_sat = zeros(32,1);
conf_sat(sat) = +1;

%no cycle-slips working with code only
conf_cs = zeros(32,1);

%number of visible satellites
nsat = size(sat,1);

%------------------------------------------------------------------------------------
% OBSERVATION EQUATIONS
%------------------------------------------------------------------------------------

%if the number of visible satellites is equal or greater than min_nsat
if (nsat(end) >= min_nsat)

    %standalone code ROVER positioning
    if (phase(1) == 1)
        [pos_R, cov_pos_R] = code_SA(X_t1_t([1,o1+1,o2+1]), pr1_Rsat, time, Eph, iono);
    else
        [pos_R, cov_pos_R] = code_SA(X_t1_t([1,o1+1,o2+1]), pr2_Rsat, time, Eph, iono);
    end

    %zeroes vector useful in matrix definitions
    Z_1_om = zeros(1,o1-1);

    %computation of the H matrix for code observations
    H = [1 Z_1_om 0 Z_1_om 0 Z_1_om;
         0 Z_1_om 1 Z_1_om 0 Z_1_om;
         0 Z_1_om 0 Z_1_om 1 Z_1_om];

    %Y0 vector for code observations
    y0 = pos_R;

    %covariance matrix of observations
    Cnn = cov_pos_R(1:3,1:3);
end

%------------------------------------------------------------------------------------
% SATELLITE ADDITION/LOSS
%------------------------------------------------------------------------------------

%search for a lost satellite
if (length(sat) < length(sat_old))

    check_off = 1;
end

%search for a new satellite
if (length(sat) > length(sat_old))

    check_on = 1;
end

%------------------------------------------------------------------------------------
% PIVOT CHANGE
%------------------------------------------------------------------------------------

%search for a possible PIVOT change
if (pivot ~= pivot_old)

    check_pivot = 1;
end

%----------------------------------------------------------------------------------------
% KALMAN FILTER
%----------------------------------------------------------------------------------------

%Kalman filter equations
if (nsat >= min_nsat)

    K = T*Cee*T' + Cvv;

    G = K*H' * (H*K*H' + Cnn)^(-1);

    Xhat_t_t = (I-G*H)*X_t1_t + G*y0;

    X_t1_t = T*Xhat_t_t;

    Cee = (I-G*H)*K;

else
    %positioning done only by the system dynamics

    Xhat_t_t = X_t1_t;

    X_t1_t = T*Xhat_t_t;

    Cee = T*Cee*T';

end

%   vvX = Xhat_t_t(2,end);
%   vvY = Xhat_t_t(o1+2,end);
%   vvZ = Xhat_t_t(o2+2,end);
%   vvv = sqrt(vvX(end)^2 + vvY(end)^2 + vvZ(end)^2);

%positioning error
%sigma_rho = sqrt(Cee(1,1,end) + Cee(o1+1,o1+1,end) + Cee(o2+1,o2+1,end));
