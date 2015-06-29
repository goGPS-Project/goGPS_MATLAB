function [kalman_initialized] = goGPS_KF_SA_code_init(XR0, time_rx, pr1, pr2, snr, Eph, SP3, iono, sbas, lambda, phase)

% SYNTAX:
%   [kalman_initialized] = goGPS_KF_SA_code_init(XR0, time_rx, pr1, pr2, snr, Eph, SP3, iono, sbas, lambda, phase);
%
% INPUT:
%   XR0 = rover approximate position (X,Y,Z)
%   time_rx = GPS reception time
%   pr1 = ROVER-SATELLITE code pseudorange (L1 carrier)
%   pr2 = ROVER-SATELLITE code pseudorange (L2 carrier)
%   snr = ROVER-SATELLITE signal-to-noise ratio
%   Eph = satellite ephemerides
%   SP3 = structure containing precise ephemeris data
%   iono = ionosphere parameters
%   sbas = SBAS corrections
%   lambda = wavelength matrix (depending on the enabled constellations)
%   phase = L1 carrier (phase=1) L2 carrier (phase=2)
%
% OUTPUT:
%   kalman_initialized = flag to point out whether Kalman has been successfully initialized
%
% DESCRIPTION:
%   Standalone code-only Kalman filter initialization.

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

global sigmaq0
global cutoff snr_threshold cond_num_threshold o1 o2 o3

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old interval
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global n_sys

kalman_initialized = 0;

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1,1);

%topocentric coordinates initialization
azR = zeros(nSatTot,1);
elR = zeros(nSatTot,1);
distR = zeros(nSatTot,1);
azM = zeros(nSatTot,1);
elM = zeros(nSatTot,1);
distM = zeros(nSatTot,1);

%--------------------------------------------------------------------------------------------
% KALMAN FILTER DYNAMIC MODEL
%--------------------------------------------------------------------------------------------

%zero vector useful in matrix definitions
Z_o1_o1 = zeros(o1);

%T matrix construction - system dynamics
T0 = eye(o1) + diag(ones(o1-1,1),1)*interval;

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

%--------------------------------------------------------------------------------------------
% SATELLITE SELECTION
%--------------------------------------------------------------------------------------------

if (length(phase) == 2)
    sat = find( (pr1 ~= 0) & (pr2 ~= 0) );
else
    if (phase == 1)
        sat = find( pr1 ~= 0 );
    else
        sat = find( pr2 ~= 0 );
    end
end
sat = sat(ismember(sat, Eph(30,:)));

%--------------------------------------------------------------------------------------------
% SBAS FAST CORRECTIONS
%--------------------------------------------------------------------------------------------

if (~isempty(sbas))
    %apply SBAS fast (pseudorange) corrections
    pr1(sat) = pr1(sat) + sbas.prc(sat)';
end

%--------------------------------------------------------------------------------------------
% APPROXIMATE POSITION
%--------------------------------------------------------------------------------------------

if ((sum(abs(XR0)) == 0) || isempty(XR0))
    %approximate position not available
    flag_XR = 0;
else
    %approximate position available
    flag_XR = 1;
end

%--------------------------------------------------------------------------------------------
% KALMAN FILTER INITIAL STATE
%--------------------------------------------------------------------------------------------

%zero vector useful in matrix definitions
Z_om_1 = zeros(o1-1,1);

min_nsat_LS = 3 + n_sys;

if (length(sat) >= min_nsat_LS)
    
    if (phase == 1)
        [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, elR(sat), azR(sat), distR(sat), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1(sat), snr(sat), Eph, SP3, iono, sbas, XR0, [], [], sat, [], lambda(sat,:), cutoff, snr_threshold, phase, flag_XR, 0); %#ok<ASGLU>
    else
        [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, elR(sat), azR(sat), distR(sat), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2(sat), snr(sat), Eph, SP3, iono, sbas, XR0, [], [], sat, [], lambda(sat,:), cutoff, snr_threshold, phase, flag_XR, 0); %#ok<ASGLU>
    end
    
    %--------------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING
    %--------------------------------------------------------------------------------------------
    
    conf_sat(:,1) = zeros(nSatTot,1);
    conf_sat(sat,1) = +1;
    
    %no cycle-slips when working with code only
    conf_cs = zeros(nSatTot,1);
    
    pivot_old = 0;
    
    %current pivot
    [null_max_elR, i] = max(elR(sat)); %#ok<ASGLU>
    pivot = sat(i);
    
    %if at least min_nsat_LS satellites are available after the cutoffs, and if the 
    % condition number in the least squares does not exceed the threshold
    if (size(sat,1) >= min_nsat_LS && cond_num < cond_num_threshold)
        
        if isempty(cov_XR) %if it was not possible to compute the covariance matrix
            cov_XR = sigmaq0 * eye(3);
        end
        sigma2_XR = diag(cov_XR);
    else
        return
    end
else
    return
end

%initial state (position and velocity)
Xhat_t_t = [XR(1); Z_om_1; XR(2); Z_om_1; XR(3); Z_om_1];

%estimation at time t+1
X_t1_t = T*Xhat_t_t;

%--------------------------------------------------------------------------------------------
% INITIAL STATE COVARIANCE MATRIX
%--------------------------------------------------------------------------------------------

Cee(:,:) = zeros(o3);
Cee(1,1) = sigma2_XR(1);
Cee(o1+1,o1+1) = sigma2_XR(2);
Cee(o2+1,o2+1) = sigma2_XR(3);
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

kalman_initialized = 1;
