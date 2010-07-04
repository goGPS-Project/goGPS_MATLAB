function kalman_goGPS_cod_init (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ...
         pr2_Rsat, pr2_Msat, snr_R, snr_M, phase)

% SYNTAX:
%   kalman_goGPS_cod_init (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ...
%   pr2_Rsat, pr2_Msat, snr_R, snr_M, phase);
%
% INPUT:
%   pos_M = MASTER position (X,Y,Z)
%   time = GPS time
%   Eph = satellite ephemerides
%   iono = ionosphere parameters
%   pr1_Rsat = ROVER-SATELLITE code pseudorange (L1 carrier)
%   pr1_Msat = MASTER-SATELLITE code pseudorange (L1 carrier)
%   pr2_Rsat = ROVER-SATELLITE code pseudorange (L2 carrier)
%   pr2_Msat = MASTER-SATELLITE code pseudorange (L2 carrier)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   phase = L1 carrier (phase=1) L2 carrier (phase=2)
%
% DESCRIPTION:
%   Code-only Kalman filter initialization.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
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

global sigmaq0 sigmaq_velx sigmaq_vely sigmaq_velz
global cutoff o1 o2 o3

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP

%--------------------------------------------------------------------------------------------
% KALMAN FILTER DYNAMIC MODEL
%--------------------------------------------------------------------------------------------

%zero vector useful in matrix definitions
Z_o1_o1 = zeros(o1);

%T matrix construction - system dynamics
T0 = eye(o1) + diag(ones(o1-1,1),1);

%second degree polynomial
% T0 = [1 1; 0 1];
%third degree polynomial
% T0 = [1 1 0; 0 1 1; 0 0 1]

%system dynamics
%X(t+1)  = X(t) + Vx(t)
%Vx(t+1) = Vx(t)
%... <-- same for Y and Z

T = [T0      Z_o1_o1 Z_o1_o1;
     Z_o1_o1 T0      Z_o1_o1;
     Z_o1_o1 Z_o1_o1 T0];

%identity matrix for following computations
I = eye(o3);

%model error covariance matrix
Cvv = zeros(o3);
Cvv(o1,o1) = sigmaq_velx;
Cvv(o2,o2) = sigmaq_vely;
Cvv(o3,o3) = sigmaq_velz; %#ok<NASGU>

%--------------------------------------------------------------------------------------------
% SATELLITE SELECTION
%--------------------------------------------------------------------------------------------

if (length(phase) == 2)
    sat = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) & ...
                (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) );
else
    if (phase == 1)
        sat = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) );
    else
        sat = find( (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) );
    end
end

%--------------------------------------------------------------------------------------------
% ESTIMATION OF INITIAL POSITION BY BANCROFT ALGORITHM
%--------------------------------------------------------------------------------------------

if (length(sat) >= 4)
    [pos_R, pos_SAT] = input_bancroft(pr1_Rsat(sat), sat, time(1), Eph);
else
    error('%d satellites are not enough to apply Bancroft algorithm\n', length(sat));
end

pos_R = pos_R(1:3);
pos_SAT = pos_SAT(:,1:3);

%------------------------------------------------------------------------------------
% CHECK SATELLITE ELEVATION, PIVOT AND CUT-OFF
%-----------------------------------------------------------------------------------

%initialization
azR = zeros(32,1);
elR = zeros(32,1);
distR = zeros(32,1);
azM = zeros(32,1);
elM = zeros(32,1);
distM = zeros(32,1);

%satellite azimuth, elevation, ROVER-SATELLITE distance
[azR(sat), elR(sat), distR(sat)] = topocent(pos_R, pos_SAT);

%elevation cut-off
sat_cutoff = find(elR > cutoff);
sat = intersect(sat,sat_cutoff);

%previous pivot
pivot_old = 0;

%actual pivot
[null_max_elR, i] = max(elR(sat)); %#ok<ASGLU>
pivot = sat(i);

%--------------------------------------------------------------------------------------------
% SATELLITE CONFIGURATION
%--------------------------------------------------------------------------------------------

%satellite configuration
conf_sat = zeros(32,1);
conf_sat(sat,1) = +1;

%no cycle-slips when working with code only
conf_cs = zeros(32,1);

%--------------------------------------------------------------------------------------------
% KALMAN FILTER INITIAL STATE
%--------------------------------------------------------------------------------------------

%zero vector useful in matrix definitions
Z_om_1 = zeros(o1-1,1);

%ROVER positioning by code double differences
if (phase(1) == 1)
    [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat), snr_R(sat), pos_M, pr1_Msat(sat), snr_M(sat), time, sat, pivot, Eph, iono); %#ok<NASGU>
else
    [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat), snr_R(sat), pos_M, pr2_Msat(sat), snr_M(sat), time, sat, pivot, Eph, iono); %#ok<NASGU>
end

%second iteration to improve the accuracy
%obtained in the previous step (from some meters to some centimeters)
if (phase(1) == 1)
    [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_double_diff(pos_R, pr1_Rsat(sat), snr_R(sat), pos_M, pr1_Msat(sat), snr_M(sat), time, sat, pivot, Eph, iono);
else
    [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_double_diff(pos_R, pr2_Rsat(sat), snr_R(sat), pos_M, pr2_Msat(sat), snr_M(sat), time, sat, pivot, Eph, iono);
end

if isempty(cov_pos_R) %if it was not possible to compute the covariance matrix
    cov_pos_R = sigmaq0 * eye(3);
end
sigmaq_pos_R = diag(cov_pos_R);

%initial state (position and velocity)
Xhat_t_t = [pos_R(1); Z_om_1; pos_R(2); Z_om_1; pos_R(3); Z_om_1];

%estimation at time t+1
X_t1_t = T*Xhat_t_t;

%--------------------------------------------------------------------------------------------
% INITIAL STATE COVARIANCE MATRIX
%--------------------------------------------------------------------------------------------

Cee(:,:) = zeros(o3);
Cee(1,1) = sigmaq_pos_R(1);
Cee(o1+1,o1+1) = sigmaq_pos_R(2);
Cee(o2+1,o2+1) = sigmaq_pos_R(3);
Cee(2:o1,2:o1) = sigmaq0 * eye(o1-1);
Cee(o1+2:o2,o1+2:o2) = sigmaq0 * eye(o1-1);
Cee(o2+2:o3,o2+2:o3) = sigmaq0 * eye(o1-1);

%--------------------------------------------------------------------------------------------
% INITIAL KALMAN FILTER DOP
%--------------------------------------------------------------------------------------------

%covariance propagation
Cee_XYZ = Cee([1 o1+1 o2+1],[1 o1+1 o2+1]);
Cee_ENU = global2localCov(Cee_XYZ, Xhat_t_t([1 o1+1 o2+1]));

%KF DOP computation
KPDOP = sqrt(Cee_XYZ(1,1) + Cee_XYZ(2,2) + Cee_XYZ(3,3));
KHDOP = sqrt(Cee_ENU(1,1) + Cee_ENU(2,2));
KVDOP = sqrt(Cee_ENU(3,3));
