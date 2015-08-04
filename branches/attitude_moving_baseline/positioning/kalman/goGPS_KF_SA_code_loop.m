function [check_on, check_off, check_pivot, check_cs] = goGPS_KF_SA_code_loop(time_rx, pr1, pr2, snr, Eph, SP3, iono, sbas, lambda, phase)

% SYNTAX:
%   [check_on, check_off, check_pivot, check_cs] = goGPS_KF_SA_code_loop(time_rx, pr1, pr2, snr, Eph, SP3, iono, sbas, lambda, phase);
%
% INPUT:
%   time_rx = GPS reception time
%   pr1 = ROVER-SATELLITE code pseudorange (L1 carrier)
%   pr2 = ROVER-SATELLITE code pseudorange (L2 carrier)
%   snr = signal-to-noise ratio for ROVER observations
%   Eph = satellite ephemerides
%   SP3 = structure containing precise ephemeris data
%   iono = ionospheric parameters
%   sbas = SBAS corrections
%   lambda = wavelength matrix (depending on the enabled constellations)
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

global sigmaq0 sigmaq_vE sigmaq_vN sigmaq_vU
global cutoff snr_threshold cond_num_threshold o1 o2 o3

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global n_sys

%----------------------------------------------------------------------------------------
% INITIALIZATION
%----------------------------------------------------------------------------------------

%output variables to point out events (satellite addition, losses, etc)
check_on = 0;
check_off = 0;
check_pivot = 0;
check_cs = 0;

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1,1);

%azimuth, elevation and ROVER-satellite and MASTER-satellite distances
azR = zeros(nSatTot,1);
azM = zeros(nSatTot,1);
elR = zeros(nSatTot,1);
elM = zeros(nSatTot,1);
distR = zeros(nSatTot,1);
distM = zeros(nSatTot,1);

%----------------------------------------------------------------------------------------
% MODEL ERROR COVARIANCE MATRIX
%----------------------------------------------------------------------------------------

%re-initialization of Cvv matrix of the model error
% (if a static model is used, no noise is added)
Cvv = zeros(o3);
if (o1 > 1)
    Cvv(o1,o1) = sigmaq_vE;
    Cvv(o2,o2) = sigmaq_vN;
    Cvv(o3,o3) = sigmaq_vU;
    
    %propagate diagonal local cov matrix to global cov matrix
    Cvv([o1 o2 o3],[o1 o2 o3]) = local2globalCov(Cvv([o1 o2 o3],[o1 o2 o3]), X_t1_t([1 o1+1 o2+1]));
end

%------------------------------------------------------------------------------------
% SATELLITE SELECTION
%------------------------------------------------------------------------------------

%visible satellites
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

%previous satellite configuration
sat_old = find(conf_sat == 1);

%--------------------------------------------------------------------------------------------
% SBAS FAST CORRECTIONS
%--------------------------------------------------------------------------------------------

if (~isempty(sbas))
    %apply SBAS fast (pseudorange) corrections
    pr1(sat) = pr1(sat) + sbas.prc(sat)';
end

%------------------------------------------------------------------------------------
% OBSERVATION EQUATIONS
%------------------------------------------------------------------------------------

%if the number of visible satellites is sufficient
min_nsat_LS = 3 + n_sys;
if (size(sat,1) >= min_nsat_LS)
    
    %approximate position
    XR0 = X_t1_t([1,o1+1,o2+1]);
    flag_XR = 1;
    
    if (phase == 1)
        [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, elR(sat), azR(sat), distR(sat), is_GLO, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1(sat), snr(sat), Eph, SP3, iono, sbas, XR0, [], [], sat, [], lambda(sat,:), cutoff, snr_threshold, phase, flag_XR, 0); %#ok<ASGLU>
    else
        [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, elR(sat), azR(sat), distR(sat), is_GLO, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2(sat), snr(sat), Eph, SP3, iono, sbas, XR0, [], [], sat, [], lambda(sat,:), cutoff, snr_threshold, phase, flag_XR, 0); %#ok<ASGLU>
    end
    
    %----------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION (NOT USED)
    %----------------------------------------------------------------------------------------

    %satellite configuration
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat) = +1;
    
    %no cycle-slips working with code only
    conf_cs = zeros(nSatTot,1);
    
    %previous pivot
    if (pivot ~= 0)
        pivot_old = pivot;
    end
    
    %current pivot
    [null_max_elR, i] = max(elR(sat)); %#ok<ASGLU>
    pivot = sat(i);
    
    %if at least min_nsat_LS satellites are available after the cutoffs, and if the 
    % condition number in the least squares does not exceed the threshold
    if (size(sat,1) >= min_nsat_LS && cond_num < cond_num_threshold)

        if isempty(cov_XR) %if it was not possible to compute the covariance matrix
            cov_XR = sigmaq0 * eye(3);
        end
        
        %zeroes vector useful in matrix definitions
        Z_1_om = zeros(1,o1-1);
        
        %computation of the H matrix for code observations
        H = [1 Z_1_om 0 Z_1_om 0 Z_1_om;
            0 Z_1_om 1 Z_1_om 0 Z_1_om;
            0 Z_1_om 0 Z_1_om 1 Z_1_om];
        
        %Y0 vector for code observations
        y0 = XR;
        
        %covariance matrix of observations
        Cnn = cov_XR(1:3,1:3);
    else
        %to point out that notwithstanding the satellite configuration,
        %data were not analysed (motion by dynamics only).
        pivot = 0;
    end
else
    %to point out that notwithstanding the satellite configuration,
    %data were not analysed (motion by dynamics only).
    pivot = 0;
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
if (pivot ~= pivot_old && pivot_old ~= 0)

    check_pivot = 1;
end

%----------------------------------------------------------------------------------------
% KALMAN FILTER
%----------------------------------------------------------------------------------------

%Kalman filter equations
if (size(sat,1) >= min_nsat_LS && cond_num < cond_num_threshold)

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

%--------------------------------------------------------------------------------------------
% KALMAN FILTER DOP
%--------------------------------------------------------------------------------------------

%covariance propagation
Cee_XYZ = Cee([1 o1+1 o2+1],[1 o1+1 o2+1]);
Cee_ENU = global2localCov(Cee_XYZ, Xhat_t_t([1 o1+1 o2+1]));

%modified DOP computation
KPDOP = sqrt(Cee_XYZ(1,1) + Cee_XYZ(2,2) + Cee_XYZ(3,3));
KHDOP = sqrt(Cee_ENU(1,1) + Cee_ENU(2,2));
KVDOP = sqrt(Cee_ENU(3,3));

%   vvX = Xhat_t_t(2,end);
%   vvY = Xhat_t_t(o1+2,end);
%   vvZ = Xhat_t_t(o2+2,end);
%   vvv = sqrt(vvX(end)^2 + vvY(end)^2 + vvZ(end)^2);

%positioning error
%sigma_rho = sqrt(Cee(1,1,end) + Cee(o1+1,o1+1,end) + Cee(o2+1,o2+1,end));
