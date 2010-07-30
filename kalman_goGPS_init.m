function kalman_goGPS_init (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ...
         ph1_Rsat, ph1_Msat, pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, snr_R, snr_M, phase)

% SYNTAX:
%   kalman_goGPS_init (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ...
%   ph1_Rsat, ph1_Msat, pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, snr_R,
%   snr_M, phase);
%
% INPUT:
%   pos_M = master known position (X,Y,Z)
%   time = GPS time
%   Eph = satellites ephemerides
%   iono =  ionospheric parameters (vector of zeroes if not available)
%   pr1_Rsat = ROVER-SATELLITE code pseudorange (carrier L1)
%   pr1_Msat = MASTER-SATELLITE code pseudorange (carrier L1)
%   ph1_Rsat = ROVER-SATELLITE phase observation (carrier L1)
%   ph1_Msat = MASTER-SATELLITE phase observation (carrier L1)
%   pr2_Rsat = ROVER-SATELLITE code pseudorange (carrier L2)
%   pr2_Msat = MASTER-SATELLITE code pseudorange (carrier L2)
%   ph2_Rsat = ROVER-SATELLITE phase observation (carrier L2)
%   ph2_Msat = MASTER-SATELLITE phase observation (carrier L2)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   phase = carrier L1 (phase=1) carrier L2 (phase=2)
%
% DESCRIPTION:
%   Kalman filter initialization with the computation of the ROVER
%   initial position (X,Y,Z).

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

%--------------------------------------------------------------------------------------------
% KALMAN FILTER PARAMETERS
%--------------------------------------------------------------------------------------------

global sigmaq0 sigmaq_velx sigmaq_vely sigmaq_velz sigmaq0_N
global cutoff snr_threshold o1 o2 o3 nN

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global flag_LS_N_estim

%--------------------------------------------------------------------------------------------
% SELECTION SINGLE / DOUBLE FREQUENCY
%--------------------------------------------------------------------------------------------

%number of unknown phase ambiguities
if (length(phase) == 1)
    nN = 32;
else
    nN = 64;
end

%--------------------------------------------------------------------------------------------
% DYNAMIC MODEL OF THE KALMAN FILTER
%--------------------------------------------------------------------------------------------

%zeroes vectors useful in the matrices definition
Z_nN_o1 = zeros(nN,o1);
Z_o1_nN = zeros(o1,nN);
Z_o1_o1 = zeros(o1);

%T matrix construction - system dynamics
%position and velocity equations
T0 = eye(o1) + diag(ones(o1-1,1),1);

%second degree polynomial
% T0 = [1 1; 0 1];
%third degree polynomial
% T0 = [1 1 0; 0 1 1; 0 0 1]

%matrix structure of initial comb_N
N0 = eye(nN);

%system dynamics
%X(t+1)  = X(t) + Vx(t)
%Vx(t+1) = Vx(t)
%... <-- for the other two variables Y e Z
%comb_N(t+1) = comb_N(t)

T = [T0      Z_o1_o1 Z_o1_o1 Z_o1_nN;
     Z_o1_o1 T0      Z_o1_o1 Z_o1_nN;
     Z_o1_o1 Z_o1_o1 T0      Z_o1_nN;
     Z_nN_o1 Z_nN_o1 Z_nN_o1 N0];

%construction of an identity matrix of 38 variables (6 for position and
%velocity + 32 or 64 for the satellites number) for the further computations
I = eye(o3+nN);

%model error covariance matrix
Cvv = zeros(o3+nN);
Cvv(o1,o1) = sigmaq_velx;
Cvv(o2,o2) = sigmaq_vely;
Cvv(o3,o3) = sigmaq_velz; %#ok<NASGU>

%--------------------------------------------------------------------------------------------
% SELECTION OF THE SATELLITES
%--------------------------------------------------------------------------------------------

if (length(phase) == 2)
    sat_pr = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) & (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) );
    sat = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) & (ph1_Rsat ~= 0) & (ph1_Msat ~= 0) & ...
                (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) & (ph2_Rsat ~= 0) & (ph2_Msat ~= 0) );
else
    if (phase == 1)
        sat_pr = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) );
        sat = find( (pr1_Rsat ~= 0) & (pr1_Msat ~= 0) & ...
                    (ph1_Rsat ~= 0) & (ph1_Msat ~= 0) );
    else
        sat_pr = find( (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) );
        sat = find( (pr2_Rsat ~= 0) & (pr2_Msat ~= 0) & ...
                    (ph2_Rsat ~= 0) & (ph2_Msat ~= 0) );
    end
end

%only satellites with code and phase
%sat_pr = sat;

%--------------------------------------------------------------------------------------------
% ESTIMATION OF THE INITIAL POSITION BY BANCROFT ALGORITHM
%--------------------------------------------------------------------------------------------

if (length(sat_pr) >= 4)
    [pos_R, pos_SAT] = input_bancroft(pr1_Rsat(sat_pr), sat_pr, time(1), Eph);
else
    error('%d satellites are not enough to apply Bancroft algorithm\n', length(sat_pr));
end

pos_R = pos_R(1:3);
pos_SAT = pos_SAT(:,1:3);

%------------------------------------------------------------------------------------
% SATELLITE ELEVATION CONTROL, PIVOT AND CUT-OFF
%-----------------------------------------------------------------------------------

%initialization
azR = zeros(32,1);
azM = zeros(32,1);
elR = zeros(32,1);
elM = zeros(32,1);
distR = zeros(32,1);
distM = zeros(32,1);

%azimuth, elevation, ROVER-SATELLITE distance computation
[azR(sat_pr), elR(sat_pr), distR(sat_pr)] = topocent(pos_R, pos_SAT);

%azimuth, elevation, MASTER-SATELLITE distance computation
[azM(sat_pr), elM(sat_pr), distM(sat_pr)] = topocent(pos_M, pos_SAT);

%elevation cut-off and signal-to-noise ratio threshold
sat_cutoff = find(elR > cutoff | (snr_R ~= 0 & snr_R < snr_threshold));
sat_pr = intersect(sat_pr,sat_cutoff);
sat = intersect(sat,sat_cutoff);

%previous pivot
pivot_old = 0;

%actual pivot
if ~isempty(sat)
    [max_elR, i] = max(elR(sat)); %#ok<ASGLU>
    pivot = sat(i);
else
    [max_elR, i] = max(elR(sat_pr)); %#ok<ASGLU>
    pivot = sat_pr(i);
end

%--------------------------------------------------------------------------------------------
% SATELLITES CONFIGURATION
%--------------------------------------------------------------------------------------------

%satellites configuration: code only (-1), both code and phase (+1);
conf_sat = zeros(32,1);
conf_sat(sat_pr) = -1;
conf_sat(sat) = +1;

%cycle-slip configuration (no cycle-slip)
conf_cs = zeros(32,1);

%--------------------------------------------------------------------------------------------
% KALMAN FILTER INITIAL STATE
%--------------------------------------------------------------------------------------------

%zeroes vectors useful in the matrices definition
Z_om_1 = zeros(o1-1,1);
sigmaq_comb_N = zeros(nN,1);

%ROVER positioning with code double differences
if (phase(1) == 1)
    if (sum(abs(iono)) == 0) %if ionospheric parameters are not available they are set equal to 0
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr1_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph);
    else
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr1_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono);
    end
else
    if (sum(abs(iono)) == 0) %if ionospheric parameters are not available they are set equal to 0
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr2_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph);
    else
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr2_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono);
    end
end

%do not use least squares ambiguity estimation
% NOTE: LS amb. estimation is automatically switched off if the number of
% satellites with phase available is not sufficient
if (~flag_LS_N_estim) | (size(sat) < 4)
%if (size(sat) < 4)
    
    %ROVER positioning with code double differences
    if (phase(1) == 1)
        if (sum(abs(iono)) == 0) %if ionospheric parameters are not available they are set equal to 0
            [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_double_diff(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr1_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph);
        else
            [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_double_diff(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr1_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono);
        end
    else
        if (sum(abs(iono)) == 0) %if ionospheric parameters are not available they are set equal to 0
            [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_double_diff(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr2_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph);
        else
            [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_double_diff(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr2_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono);
        end
    end
    if isempty(cov_pos_R) %if it was not possible to compute the covariance matrix
        cov_pos_R = sigmaq0 * eye(3);
    end
    sigmaq_pos_R = diag(cov_pos_R);
    
    %satellite combinations initialization: initialized value
    %if the satellite is visible, 0 if the satellite is not visible
    comb_N1_stim = zeros(32,1);
    comb_N2_stim = zeros(32,1);
    sigmaq_comb_N1 = zeros(32,1);
    sigmaq_comb_N2 = zeros(32,1);
    
    %computation of the phase double differences in order to estimate N
    if ~isempty(sat)
        [comb_N1_stim(sat), sigmaq_comb_N1(sat)] = amb_estimate_observ(pr1_Rsat(sat), pr1_Msat(sat), ph1_Rsat(sat), ph1_Msat(sat), pivot, sat, 1);
        [comb_N2_stim(sat), sigmaq_comb_N2(sat)] = amb_estimate_observ(pr2_Rsat(sat), pr2_Msat(sat), ph2_Rsat(sat), ph2_Msat(sat), pivot, sat, 2);
    end
    
    if (length(phase) == 2)
        comb_N_stim = [comb_N1_stim; comb_N2_stim];
        sigmaq_comb_N = [sigmaq_comb_N1; sigmaq_comb_N2];
    else
        if (phase == 1)
            comb_N_stim = comb_N1_stim;
            sigmaq_comb_N = sigmaq_comb_N1;
        else
            comb_N_stim = comb_N2_stim;
            sigmaq_comb_N = sigmaq_comb_N2;
        end
    end

%use least squares ambiguity estimation
else
    
    %satellite combinations initialization: initialized value
    %if the satellite is visible, 0 if the satellite is not visible
    comb_N1_stim = zeros(32,1);
    comb_N2_stim = zeros(32,1);

    %ROVER positioning improvement with code and phase double differences
    if ~isempty(sat)
        [     pos_R,      cov_pos_R, comb_N1_stim(sat), cov_comb_N1_stim, PDOP, HDOP, VDOP] = code_phase_double_diff(pos_R, pr1_Rsat(sat), ph1_Rsat(sat), snr_R(sat), pos_M, pr1_Msat(sat), ph1_Msat(sat), snr_M(sat), time, sat, pivot, Eph, 1, iono);
        [null_pos_R, null_cov_pos_R, comb_N2_stim(sat), cov_comb_N2_stim]    = code_phase_double_diff(pos_R, pr2_Rsat(sat), ph2_Rsat(sat), snr_R(sat), pos_M, pr2_Msat(sat), ph2_Msat(sat), snr_M(sat), time, sat, pivot, Eph, 2, iono); %#ok<ASGLU>
    end
    
    if isempty(cov_pos_R) %if it was not possible to compute the covariance matrix
        cov_pos_R = sigmaq0 * eye(3);
    end
    sigmaq_pos_R = diag(cov_pos_R);
    
    if isempty(cov_comb_N1_stim) %if it was not possible to compute the covariance matrix
        cov_comb_N1_stim = sigmaq0_N * eye(length(sat));
    end
    
    if isempty(cov_comb_N2_stim) %if it was not possible to compute the covariance matrix
        cov_comb_N2_stim = sigmaq0_N * eye(length(sat));
    end
    
    if (length(phase) == 2)
        comb_N_stim = [comb_N1_stim; comb_N2_stim];
        sigmaq_comb_N(sat) = diag(cov_comb_N1_stim);
        sigmaq_comb_N(sat+nN) = diag(cov_comb_N2_stim);
    else
        if (phase == 1)
            comb_N_stim = comb_N1_stim;
            sigmaq_comb_N(sat) = diag(cov_comb_N1_stim);
        else
            comb_N_stim = comb_N2_stim;
            sigmaq_comb_N(sat) = diag(cov_comb_N2_stim);
        end
    end
end

%initialization of the initial point with 6(positions and velocities) +
%32 or 64 (N combinations) variables
Xhat_t_t = [pos_R(1); Z_om_1; pos_R(2); Z_om_1; pos_R(3); Z_om_1; comb_N_stim];

%point estimation at step t+1 X Vx Y Vy Z Vz comb_N
%estimation at step t, because the initial velocity is equal to 0
X_t1_t = T*Xhat_t_t;

%--------------------------------------------------------------------------------------------
% INITIAL STATE COVARIANCE MATRIX
%--------------------------------------------------------------------------------------------

%initial state covariance matrix
Cee(:,:) = zeros(o3+nN);
Cee(1,1) = sigmaq_pos_R(1);
Cee(o1+1,o1+1) = sigmaq_pos_R(2);
Cee(o2+1,o2+1) = sigmaq_pos_R(3);
Cee(2:o1,2:o1) = sigmaq0 * eye(o1-1);
Cee(o1+2:o2,o1+2:o2) = sigmaq0 * eye(o1-1);
Cee(o2+2:o3,o2+2:o3) = sigmaq0 * eye(o1-1);
Cee(o3+1:o3+nN,o3+1:o3+nN) = diag(sigmaq_comb_N);

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