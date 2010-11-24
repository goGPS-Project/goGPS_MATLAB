function kalman_goGPS_SA_init (time, Eph, iono, pr1_Rsat, ph1_Rsat, pr2_Rsat, ph2_Rsat, snr_R, phase)

% SYNTAX:
%   kalman_goGPS_SA_init (time, Eph, iono, pr1_Rsat, ph1_Rsat, pr2_Rsat, ph2_Rsat, snr_R, phase);
%
% INPUT:
%   time = GPS time
%   Eph = satellite ephemerides
%   iono = ionosphere parameters
%   pr1_Rsat = ROVER-SATELLITE code pseudorange (L1 carrier)
%   ph1_Rsat = ROVER-SATELLITE phase observation (carrier L1)
%   pr2_Rsat = ROVER-SATELLITE code pseudorange (L2 carrier)
%   ph2_Rsat = ROVER-SATELLITE phase observation (carrier L2)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   phase = L1 carrier (phase=1) L2 carrier (phase=2)
%
% DESCRIPTION:
%   Standalone phase and code Kalman filter initialization.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni, Eugenio Realini
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

global sigmaq0 sigmaq0_N
global cutoff o1 o2 o3 nN

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP

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
% KALMAN FILTER DYNAMIC MODEL
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

%--------------------------------------------------------------------------------------------
% SATELLITE SELECTION
%--------------------------------------------------------------------------------------------

if (length(phase) == 2)
    sat_pr = find( (pr1_Rsat ~= 0) & (pr2_Rsat ~= 0) );
    sat = find( (pr1_Rsat ~= 0) & (ph1_Rsat ~= 0) & ...
                (pr2_Rsat ~= 0) & (ph2_Rsat ~= 0) );
else
    if (phase == 1)
        sat_pr = find( (pr1_Rsat ~= 0) );
        sat = find( (pr1_Rsat ~= 0) & ...
                    (ph1_Rsat ~= 0) );
    else
        sat_pr = find( (pr2_Rsat ~= 0) );
        sat = find( (pr2_Rsat ~= 0) & ...
                    (ph2_Rsat ~= 0) );
    end
end

%only satellites with code and phase
%sat_pr = sat;

%--------------------------------------------------------------------------------------------
% ESTIMATION OF INITIAL POSITION BY BANCROFT ALGORITHM
%--------------------------------------------------------------------------------------------

if (length(sat_pr) >= 4)
    [pos_R, pos_SAT] = input_bancroft(pr1_Rsat(sat_pr), sat_pr, time(1), Eph);
else
    error('%d satellites are not enough to apply Bancroft algorithm\n', length(sat_pr));
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
[azR(sat_pr), elR(sat_pr), distR(sat_pr)] = topocent(pos_R, pos_SAT);

%elevation cut-off 
sat_cutoff = find(elR > cutoff);
sat_pr = intersect(sat_pr,sat_cutoff);
sat = intersect(sat,sat_cutoff);

pivot_old = 0;

%current pivot
[null_max_elR, i] = max(elR(sat_pr)); %#ok<ASGLU>
pivot = sat_pr(i);

%--------------------------------------------------------------------------------------------
% SATELLITE CONFIGURATION
%--------------------------------------------------------------------------------------------

%satellites configuration: code only (-1), both code and phase (+1);
conf_sat = zeros(32,1);
conf_sat(sat_pr) = -1;
conf_sat(sat) = +1;

%cycle-slip configuration
conf_cs = zeros(32,1);

%--------------------------------------------------------------------------------------------
% KALMAN FILTER INITIAL STATE
%--------------------------------------------------------------------------------------------

%zero vector useful in matrix definitions
Z_om_1 = zeros(o1-1,1);
sigmaq_N = zeros(nN,1);

if (length(sat_pr) >= 4)
    
    %stand-alone ROVER positioning with code
    if (phase(1) == 1)
        [pos_R, cov_pos_R] = code_SA(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph, iono);
    else
        [pos_R, cov_pos_R] = code_SA(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph, iono);
    end
    
    %second iteration to improve the accuracy
    %obtained in the previous step (from some meters to some centimeters)
    if (phase(1) == 1)
        [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_SA(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph, iono);
    else
        [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_SA(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph, iono);
    end
    
    if isempty(cov_pos_R) %if it was not possible to compute the covariance matrix
        cov_pos_R = sigmaq0 * eye(3);
    end
    sigmaq_pos_R = diag(cov_pos_R);
    
else
    error('%d satellites are not enough to make ROVER positioning in stand-alone\n', length(sat_pr));
end

%do not use least squares ambiguity estimation
% NOTE: LS amb. estimation is automatically switched off if the number of
% satellites with phase available is not sufficient
if (length(sat) < 4)
    
    %ROVER positioning with code double differences
    if (phase(1) == 1)
         [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_SA(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph, iono);
    else
         [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_SA(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), sat_pr, time, Eph, iono);
    end
    if isempty(cov_pos_R) %if it was not possible to compute the covariance matrix
        cov_pos_R = sigmaq0 * eye(3);
    end
    sigmaq_pos_R = diag(cov_pos_R);
    
    %satellite combinations initialization: initialized value
    %if the satellite is visible, 0 if the satellite is not visible
    N1_stim = zeros(32,1);
    N2_stim = zeros(32,1);
    sigmaq_N1 = zeros(32,1);
    sigmaq_N2 = zeros(32,1);
    
    %computation of the phase double differences in order to estimate N
    if ~isempty(sat)
        [N1_stim(sat), sigmaq_N1(sat)] = amb_estimate_observ_SA(pr1_Rsat(sat), ph1_Rsat(sat), 1);
        [N2_stim(sat), sigmaq_N2(sat)] = amb_estimate_observ_SA(pr2_Rsat(sat), ph2_Rsat(sat), 2);
    end

    if (length(phase) == 2)
        N_stim = [N1_stim; N2_stim];
        sigmaq_N = [sigmaq_N1; sigmaq_N2];
    else
        if (phase == 1)
            N_stim = N1_stim;
            sigmaq_N = sigmaq_N1;
        else
            N_stim = N2_stim;
            sigmaq_N = sigmaq_N2;
        end
    end

%use least squares ambiguity estimation
else
    
    %satellite combinations initialization: initialized value
    %if the satellite is visible, 0 if the satellite is not visible
    N1_stim = zeros(32,1);
    N2_stim = zeros(32,1);

    %ROVER positioning improvement with code and phase double differences
    if ~isempty(sat)
        [     pos_R,      cov_pos_R, N1_stim(sat), cov_N1_stim, PDOP, HDOP, VDOP] = code_phase_SA(pos_R, pr1_Rsat(sat), ph1_Rsat(sat), snr_R(sat), sat, time, Eph, 1, iono);
        [null_pos_R, null_cov_pos_R, N2_stim(sat), cov_N2_stim] = code_phase_SA(pos_R, pr2_Rsat(sat), ph2_Rsat(sat), snr_R(sat), sat, time, Eph, 2, iono); %#ok<ASGLU>
    end
    
    if isempty(cov_pos_R) %if it was not possible to compute the covariance matrix
        cov_pos_R = sigmaq0 * eye(3);
    end
    sigmaq_pos_R = diag(cov_pos_R);
    
    if isempty(cov_N1_stim) %if it was not possible to compute the covariance matrix
        cov_N1_stim = sigmaq0_N * eye(length(sat));
    end
    
    if isempty(cov_N2_stim) %if it was not possible to compute the covariance matrix
        cov_N2_stim = sigmaq0_N * eye(length(sat));
    end
    
    if (length(phase) == 2)
        N_stim = [N1_stim; N2_stim];
        sigmaq_N(sat) = diag(cov_N1_stim);
        %sigmaq_N(sat) = (sigmaq_cod1 / lambda1^2) * ones(length(sat),1);
        sigmaq_N(sat+nN) = diag(cov_N2_stim);
        %sigmaq_N(sat+nN) = (sigmaq_cod2 / lambda2^2) * ones(length(sat),1);
    else
        if (phase == 1)
            N_stim = N1_stim;
            sigmaq_N(sat) = diag(cov_N1_stim);
            %sigmaq_N(sat) = (sigmaq_cod1 / lambda1^2) * ones(length(sat),1);
        else
            N_stim = N2_stim;
            sigmaq_N(sat) = diag(cov_N2_stim);
            %sigmaq_N(sat) = (sigmaq_cod2 / lambda2^2) * ones(length(sat),1);
        end
    end
end

%initialization of the initial point with 6(positions and velocities) +
%32 or 64 (N combinations) variables
Xhat_t_t = [pos_R(1); Z_om_1; pos_R(2); Z_om_1; pos_R(3); Z_om_1; N_stim];

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
Cee(o3+1:o3+nN,o3+1:o3+nN) = diag(sigmaq_N);

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
