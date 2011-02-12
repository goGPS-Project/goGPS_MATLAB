function [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop ...
         (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ph1_Rsat, ph1_Msat, dop1_Rsat, ...
         pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, dop2_Rsat, snr_R, snr_M, phase, ref)

% SYNTAX:
%   [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop ...
%   (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ph1_Rsat, ph1_Msat, dop1_Rsat, ...
%   pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, dop2_Rsat, snr_R, snr_M, phase, ref);
%
% INPUT:
%   pos_M = Master given cooridnates (X,Y,Z)
%   time = GPS time
%   Eph = satellites ephemerides
%   iono = ionosphere parameters  (not used)
%   pr1_Rsat  = ROVER-SATELLITE code-pseudorange (carrier L1)
%   pr1_Msat  = MASTER-SATELLITE code-pseudorange (carrier L1)
%   ph1_Rsat  = ROVER-SATELLITE phase observations (carrier L1)
%   ph1_Msat  = MASTER-SATELLITE phase observations (carrier L1)
%   dop1_Rsat = ROVER_SATELLITE Doppler observation (carrier L1)
%   pr2_Rsat  = ROVER-SATELLITE code-pseudorange (carrier L2)
%   pr2_Msat  = MASTER-SATELLITE code-pseudorange (carrier L2)
%   ph2_Rsat  = ROVER-SATELLITE phase observations (carrier L2)
%   ph2_Msat  = MASTER-SATELLITE phase observations (carrier L2)
%   dop2_Rsat = ROVER_SATELLITE Doppler observation (carrier L2)
%   snr_R = signal-to-noise ratio for ROVER observations
%   snr_M = signal-to-noise ratio for MASTER observations
%   phase = carrier L1 (phase=1), carrier L2 (phase=2)
%   ref = reference line
%
% OUTPUT:
%   check_on = boolean variable for satellite births
%   check_off = boolean variable for satellite deaths
%   check_pivot = boolean variable for changes of pivot
%   check_cs = boolean variable for cycle-slips
%
% DESCRIPTION:
%   Kalman filter for estimating the ROVER path.
%   Birth and death of satellites are considered, besides cycle slips and changes of pivot.
%   Constrained path.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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

global lambda1 lambda2

global sigmaq_vel sigmaq0_N
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph
global min_nsat cutoff snr_threshold cs_threshold o1 nN
global s0 ax ay az

global Xhat_t_t X_t1_t Yhat_t_t Y_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP
global doppler_pred_range1_R doppler_pred_range2_R

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

%re-initialization of the model error Cvv matrix
Cvv = zeros(o1+nN);
Cvv(o1,o1) = sigmaq_vel;

%------------------------------------------------------------------------------------
% RECEIVER PREDICTED POSITION
%------------------------------------------------------------------------------------
posR_app = Y_t1_t';

%------------------------------------------------------------------------------------
% SATELLITE SELECTION
%------------------------------------------------------------------------------------

%visible satellites
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

%----------------------------------------------------------------------------------------
% CONVERSION FROM CARTESIAN TO GEODETIC COORDINATES
%----------------------------------------------------------------------------------------

%approximated coordinates X Y Z
X_app = posR_app(1);
Y_app = posR_app(2);
Z_app = posR_app(3);

[phiR_app, lamR_app, hR_app] = cart2geod(X_app, Y_app, Z_app);
[phiM, lamM, hM] = cart2geod(pos_M(1), pos_M(2), pos_M(3));

%------------------------------------------------------------------------------------
% SATELLITE ELEVATION, PIVOT AND CUT-OFF
%------------------------------------------------------------------------------------

j = 1;
bad_sat = [];

posS = zeros(3,32);
dtS = zeros(32,1);
prRS_app = zeros(32,1);
prMS_app = zeros(32,1);
err_tropo_RS = zeros(32,1);
err_tropo_MS = zeros(32,1);
err_iono_RS = zeros(32,1);
err_iono_MS = zeros(32,1);

for i = 1:size(sat_pr)

    i_sat = sat_pr(i);

    %satellite position (with clock error and Earth rotation corrections)
    [posS_tmp, dtS_tmp] = sat_corr(Eph, sat_pr(i), time, pr1_Rsat(sat_pr(i)));

    if (~isempty(posS_tmp))
        
        posS(:,i_sat) = posS_tmp;
        dtS(i_sat,1) = dtS_tmp;
        
        %computation of the satellite azimuth and elevation
        [azR(i_sat), elR(i_sat), distR(i_sat)] = topocent(posR_app, posS(:,i_sat)');
        [azM(i_sat), elM(i_sat), distM(i_sat)] = topocent(pos_M, posS(:,i_sat)');

        %computation of ROVER-SATELLITE approximated pseudorange
        prRS_app(i_sat) = sqrt(sum((posR_app - posS(:,i_sat)).^2));
        prMS_app(i_sat) = sqrt(sum((pos_M - posS(:,i_sat)).^2));
        
        %computation of tropospheric errors
        err_tropo_RS(i_sat) = err_tropo(elR(i_sat), hR_app);
        err_tropo_MS(i_sat) = err_tropo(elM(i_sat), hM);
        
        %computation of ionospheric errors
        err_iono_RS(i_sat) = err_iono(iono, phiR_app*180/pi, lamR_app*180/pi, azR(i_sat), elR(i_sat), time);
        err_iono_MS(i_sat) = err_iono(iono, phiM*180/pi, lamM*180/pi, azM(i_sat), elM(i_sat), time);
    end

    %test ephemerides availability, elevation and signal-to-noise ratio
    if (isempty(posS_tmp) | elR(i_sat) < cutoff | snr_R(i_sat) < snr_threshold)
        bad_sat(j,1) = i_sat;
        j = j + 1;
    end
end

%removal of satellites without ephemerides or with elevation or SNR lower than the respective threshold
sat_pr(ismember(sat_pr,bad_sat) == 1) = [];
sat(ismember(sat,bad_sat) == 1) = [];

%previous pivot
if (pivot ~= 0)
    pivot_old = pivot;
end

%current pivot
if ~isempty(sat)
    [max_elR, i] = max(elR(sat)); %#ok<ASGLU>
    pivot = sat(i);
else
    [max_elR, i] = max(elR(sat_pr)); %#ok<ASGLU>
    pivot = sat_pr(i);
end
%pivot = find(elR == max(elR));

%----------------------------------------------------------------------------------------
% SATELLITE CONFIGURATION
%----------------------------------------------------------------------------------------

%previous satellite configuration (with phase measurements)
sat_old = find(conf_sat == 1);

%satellite configuration: code only (-1), both code and phase (+1);
conf_sat = zeros(32,1);
conf_sat(sat_pr) = -1;
conf_sat(sat) = +1;

%cycle-slip configuration
conf_cs = zeros(32,1);

%total number of visible satellites
nsat = size(sat_pr,1);
n = nsat - 1;

%if the number of visible satellites is equal or greater than min_nsat
if (nsat >= min_nsat)
    
    %------------------------------------------------------------------------------------
    % SATELLITE ADDITION/LOSS
    %------------------------------------------------------------------------------------

    sat_dead = []; %#ok<NASGU>
    sat_born = [];

    %search for a lost satellite
    if (length(sat) < length(sat_old))

        check_off = 1;

        %save lost satellites
        sat_dead = setdiff(sat_old,sat);

        %for lost satellites it is fundamental to set their N-PIVOT
        % combinations to 0. Furthermore it could be convenient to raise
        %their uncertainty (not necessary - done when a new satellite is
        %added)
        comb_N1 = 0;
        comb_N2 = 0;

        if (length(phase) == 2)
            X_t1_t(o1+sat_dead,1) = comb_N1;
            X_t1_t(o1+32+sat_dead,1) = comb_N2;
        else
            if (phase == 1)
                X_t1_t(o1+sat_dead,1) = comb_N1;
            else
                X_t1_t(o1+sat_dead,1) = comb_N2;
            end
        end
    end

    %search for a new satellite
    if (length(sat) > length(sat_old))
        
        check_on = 1;
        
        %new satellites
        sat_born = setdiff(sat,sat_old);
        sat_born = setdiff(sat_born,pivot);
    end

    %------------------------------------------------------------------------------------
    % PIVOT CHANGE
    %------------------------------------------------------------------------------------

    %search for a possible PIVOT change
    if (pivot ~= pivot_old)

        check_pivot = 1;

        %matrix construction to update the PIVOT change
        %sat: vector with the current visible satellites
        %nsat: current satellites vector dimension
        R = zeros(32);
        R(sat,sat) = eye(length(sat));
        R(sat,pivot) = -1;
        R(pivot_old,pivot_old) = 0;
        R(pivot,pivot) = 0;

        I0 = eye(o1);
        Z_32_o1 = zeros(32,o1);
        Z_o1_32 = zeros(o1,32);
        Z_32_32 = zeros(32,32);

        %total matrix construction
        %sat_old, sat
        if (length(phase) == 2)
            A = [I0 Z_o1_32 Z_o1_32; Z_32_o1 R Z_32_32; Z_32_o1 Z_32_32 R];
        else
            A = [I0 Z_o1_32; Z_32_o1 R];
        end

        %new state estimate
        X_t1_t = A*X_t1_t;

        %re-computation of the Cee covariance matrix at the previous epoch
        Cee = A*Cee*A';
    end

    %------------------------------------------------------------------------------------
    % CYCLE-SLIP
    %------------------------------------------------------------------------------------

    if ~isempty(sat)

        %Test presence/absence of a cycle-slip at the current epoch.
        %The state of the system is not changed yet
        if (length(phase) == 2)
            [check_cs1, N_slip1, sat_slip1] = cycle_slip_detection(X_t1_t(o1+1:o1+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), doppler_pred_range1_R(sat), pivot, sat, sat_born, cs_threshold, 1); %#ok<ASGLU>
            [check_cs2, N_slip2, sat_slip2] = cycle_slip_detection(X_t1_t(o1+33:o1+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), doppler_pred_range2_R(sat), pivot, sat, sat_born, cs_threshold, 2); %#ok<ASGLU>

            if (check_cs1 | check_cs2)
                check_cs = 1;
            end
        else
            if (phase == 1)
                [check_cs, N_slip, sat_slip] = cycle_slip_detection(X_t1_t(o1+1:o1+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), doppler_pred_range1_R(sat), pivot, sat, sat_born, cs_threshold, 1); %#ok<ASGLU>
            else
                [check_cs, N_slip, sat_slip] = cycle_slip_detection(X_t1_t(o1+1:o1+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), doppler_pred_range2_R(sat), pivot, sat, sat_born, cs_threshold, 2); %#ok<ASGLU>
            end
        end
    else
        sat_slip1 = [];
        sat_slip2 = [];
        sat_slip = [];
        check_cs1 = 0;
        check_cs2 = 0;
        check_cs = 0;
    end

    %------------------------------------------------------------------------------------
    % PHASE AMBIGUITY ESTIMATION
    %------------------------------------------------------------------------------------

    if (check_on | check_cs)
        if (length(phase) == 2)
            %[N1_slip, N1_born] = amb_estimate_LS(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat_slip1, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o1+sat_pr), Cee(o1+sat_pr, o1+sat_pr));
            %[N2_slip, N2_born] = amb_estimate_LS(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat_slip2, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot, phase, X_t1_t(o1+sat_pr), Cee(o1+sat_pr, o1+sat_pr));
            [N1_slip, N1_born] = amb_estimate_LS(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat_slip1, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o1+sat_pr));
            [N2_slip, N2_born] = amb_estimate_LS(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat_slip2, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot, phase, X_t1_t(o1+sat_pr));
            
            if (check_on)
                X_t1_t(o1+sat_born,1) = N1_born;
                X_t1_t(o1+32+sat_born,1) = N2_born;
                %Cvv(o1+sat_born,o1+sat_born) = sigmaq_N1_born * eye(size(sat_born,1));
                %Cvv(o1+32+sat_born,o1+32+sat_born) = sigmaq_N2_born * eye(size(sat_born,1));
                Cvv(o1+sat_born,o1+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                Cvv(o1+32+sat_born,o1+32+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            end
            
            if (check_cs1)
                conf_cs(sat_slip1) = 1;
                X_t1_t(o1+sat_slip1) = N1_slip;
                Cvv(o1+sat_slip1,o1+sat_slip1) = sigmaq0_N * eye(size(sat_slip1,1));
            end
            
            if (check_cs2)
                conf_cs(sat_slip2) = 1;
                X_t1_t(o1+32+sat_slip2) = N2_slip;
                Cvv(o1+32+sat_slip2,o1+32+sat_slip2) = sigmaq0_N * eye(size(sat_slip2,1));
            end
        else
            if (phase == 1)
                [N_slip, N_born] = amb_estimate_LS(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o1+sat_pr), Cee(o1+sat_pr, o1+sat_pr));
                %[N_slip, N_born] = amb_estimate_LS(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o1+sat_pr));
            else
                [N_slip, N_born] = amb_estimate_LS(posR_app, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), (lambda2/lambda1)^2 * err_iono_RS(sat_pr), (lambda2/lambda1)^2 * err_iono_MS(sat_pr), pivot, phase, X_t1_t(o1+sat_pr), Cee(o1+sat_pr, o1+sat_pr));
                %[N_slip, N_born] = amb_estimate_LS(posR_app, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_RS(sat_pr), err_iono_MS(sat_pr), pivot, phase, X_t1_t(o1+sat_pr));
            end
            
            if (check_on)
                X_t1_t(o1+sat_born,1) = N_born;
                %Cvv(o1+sat_born,o1+sat_born) = sigmaq_N_born * eye(size(sat_born,1));
                Cvv(o1+sat_born,o1+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            end
            
            if (check_cs)
                conf_cs(sat_slip) = 1;
                X_t1_t(o1+sat_slip) = N_slip;
                Cvv(o1+sat_slip,o1+sat_slip) = sigmaq0_N * eye(size(sat_slip,1));
            end
        end
    end
    
    %------------------------------------------------------------------------------------
    % OBSERVATION EQUATIONS
    %------------------------------------------------------------------------------------

    %rows in which the phase observation is available
    p = find(ismember(setdiff(sat_pr,pivot),setdiff(sat,pivot))==1);

    %function that allows to estimate the Kalman filter parameters                                                            
    [alfa1, prstim_pr1, prstim_ph1, ddc1, ddp1, alfa0] = input_kalman_vinc(posR_app, posS(:,sat_pr), prRS_app(sat_pr), prMS_app(sat_pr), pr1_Rsat(sat_pr), ph1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Msat(sat_pr), err_tropo_RS(sat_pr), err_iono_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_MS(sat_pr), sat_pr, pivot, 1);
    [alfa2, prstim_pr2, prstim_ph2, ddc2, ddp2       ] = input_kalman_vinc(posR_app, posS(:,sat_pr), prRS_app(sat_pr), prMS_app(sat_pr), pr2_Rsat(sat_pr), ph2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Msat(sat_pr), err_tropo_RS(sat_pr), err_iono_RS(sat_pr), err_tropo_MS(sat_pr), err_iono_MS(sat_pr), sat_pr, pivot, 2);

    %vectors of zeros useful for the matrix declaration
    Z_n_nN = zeros(n,nN);
    Z_n_om = zeros(n,o1-1);

    %H matrix estimation for the code
    H_cod1 = [alfa1 Z_n_om Z_n_nN];
    H_cod2 = [alfa2 Z_n_om Z_n_nN];
    if (length(phase) == 2)
        H_cod = [H_cod1; H_cod2];
    else
        if (phase == 1)
            H_cod = H_cod1;
        else
            H_cod = H_cod2;
        end
    end

    %lambda positions estimate
    L_fas1 = zeros(n,32);
    L_fas2 = zeros(n,32);
    v = 1;
    for u = 1 : n+1 % pivot included
        if (sat_pr(u) ~= pivot)
            %minus sign becaude of the particular data processing
            L_fas1(v,sat_pr(u)) = -(lambda1);
            L_fas2(v,sat_pr(u)) = -(lambda2);
            v = v+1;
        end
    end

    %H matrix H estimation for the phase
    if ~isempty(p)
        H_fas1 = [alfa1(p) Z_n_om(p,:) Z_n_nN(p,:)];
        H_fas2 = [alfa2(p) Z_n_om(p,:) Z_n_nN(p,:)];
        if (length(phase) == 2)
            H_fas1(:,o1+1:o1+32) = L_fas1(p,:);
            H_fas2(:,o1+33:o1+64) = L_fas2(p,:);
            H_fas = [H_fas1; H_fas2];
        else
            if (phase == 1)
                H_fas1(:,o1+1:o1+32) = L_fas1(p,:);
                H_fas = H_fas1;
            else
                H_fas2(:,o1+1:o1+32) = L_fas2(p,:);
                H_fas = H_fas2;
            end
        end
    else
        H_fas = [];
    end

    %complete H matrix construction
    H = [H_cod; H_fas];

    %Y0 vector estimate for the code
    y0_cod1 = ddc1 - prstim_pr1 + alfa1*X_t1_t(1);
    y0_cod2 = ddc2 - prstim_pr2 + alfa2*X_t1_t(1);

    %Y0 vector estimate for the phase
    if ~isempty(p)
        y0_fas1 = ddp1(p) - prstim_ph1(p) + alfa1(p)*X_t1_t(1);
        y0_fas2 = ddp2(p) - prstim_ph2(p) + alfa2(p)*X_t1_t(1);
    else
        y0_fas1 = [];
        y0_fas2 = [];
    end

    %complete Y0 vector construction
    if (length(phase) == 2)
        y0_cod = [y0_cod1; y0_cod2];
        y0_fas = [y0_fas1; y0_fas2];
    else
        if (phase == 1)
            y0_cod = y0_cod1;
            y0_fas = y0_fas1;
        else
            y0_cod = y0_cod2;
            y0_fas = y0_fas2;
        end
    end
    y0 = [y0_cod; y0_fas];

    %------------------------------------------------------------------------------------
    % OBSERVATION COVARIANCE MATRIX
    %------------------------------------------------------------------------------------

    %construction of the cofactor matrix
    Q = cofactor_matrix(elR(sat_pr), elM(sat_pr), snr_R(sat_pr), snr_M(sat_pr), sat_pr, pivot);

    %zeroes vector useful in matrix definitions
    Z_n_n = zeros(n,n);

    %matrix construction by conveniently multiplying by
    %the code variance and the phase variance
    if ~isempty(p)
        if (length(phase) == 2)
            Cnn = [sigmaq_cod1*Q(:,:) Z_n_n(:,:) Z_n_n(:,p) Z_n_n(:,p); Z_n_n(:,:) sigmaq_cod2*Q(:,:) Z_n_n(:,p) Z_n_n(:,p);
                   Z_n_n(p,:) Z_n_n(p,:) sigmaq_ph*Q(p,p) Z_n_n(p,p); Z_n_n(p,:) Z_n_n(p,:) Z_n_n(p,p) sigmaq_ph*Q(p,p)];
        else
            if (phase == 1)
                Cnn = [sigmaq_cod1*Q(:,:) Z_n_n(:,p); Z_n_n(p,:) sigmaq_ph*Q(p,p)];
            else
                Cnn = [sigmaq_cod2*Q(:,:) Z_n_n(:,p); Z_n_n(p,:) sigmaq_ph*Q(p,p)];
            end
        end
    else
        if (length(phase) == 2)
            Cnn = [sigmaq_cod1*Q Z_n_n; Z_n_n sigmaq_cod2*Q];
        else
            if (phase == 1)
                Cnn = sigmaq_cod1*Q;
            else
                Cnn = sigmaq_cod2*Q;
            end
        end
    end
    
    %------------------------------------------------------------------------------------
    % DILUTION OF PRECISION
    %------------------------------------------------------------------------------------
    
    cov_XYZ = (alfa0'*alfa0)^-1;
    cov_ENU = global2localCov(cov_XYZ, posR_app);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
    
    %--------------------------------------------------------------------------------------------
    % DOPPLER-BASED PREDICTION OF PHASE RANGES
    %--------------------------------------------------------------------------------------------
    doppler_pred_range1_R = zeros(32,1);
    doppler_pred_range2_R = zeros(32,1);
    if (dop1_Rsat(sat))
        doppler_pred_range1_R(sat,1) = ph1_Rsat(sat) - dop1_Rsat(sat);
    end
    if (dop2_Rsat(sat))
        doppler_pred_range2_R(sat,1) = ph2_Rsat(sat) - dop2_Rsat(sat);
    end

else
    %to point out that notwithstanding the satellite configuration,
    %data were not analysed (motion by dynamics only).
    pivot = 0;
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
    %positioning made by means of the kalman filter dynamic only
    %at constant velocity

    Xhat_t_t = X_t1_t;

    X_t1_t = T*Xhat_t_t;

    Cee = T*Cee*T';
end

%----------------------------------------------------------------------------------------
% CARTESIAN COORDINATES ESTIMATE
%----------------------------------------------------------------------------------------

%curvilinear coordinate localization
i = find((Xhat_t_t(1) >= s0(1:end-1)) & (Xhat_t_t(1) < s0(2:end)));

%cartesian coordinates estimate
Yhat_t_t(1,1) = ref(i,1) + ax(i) * (Xhat_t_t(1) - s0(i));
Yhat_t_t(2,1) = ref(i,2) + ay(i) * (Xhat_t_t(1) - s0(i));
Yhat_t_t(3,1) = ref(i,3) + az(i) * (Xhat_t_t(1) - s0(i));

%cartesian coordinates position
if ((Yhat_t_t(1) < min(ref(i,1),ref(i+1,1))) | (Yhat_t_t(1) > max(ref(i,1),ref(i+1,1))) | ...
    (Yhat_t_t(2) < min(ref(i,2),ref(i+1,2))) | (Yhat_t_t(2) > max(ref(i,2),ref(i+1,2))) | ...
    (Yhat_t_t(3) < min(ref(i,3),ref(i+1,3))) | (Yhat_t_t(3) > max(ref(i,3),ref(i+1,3))))

    d1 = sqrt(sum((ref(i,:) - Yhat_t_t).^2));
    d2 = sqrt(sum((ref(i+1,:) - Yhat_t_t).^2));

    if (d1 < d2)
        Xhat_t_t(1) = s0(i);
        Yhat_t_t(:) = ref(i,:)';
    else
        Xhat_t_t(1) = s0(i+1);
        Yhat_t_t(:) = ref(i+1,:)';
    end
end

%curvilinear coordinate localization
i = find((X_t1_t(1) >= s0(1:end-1)) & (X_t1_t(1) < s0(2:end)));

%cartesian coordinates estimate
Y_t1_t(1,1) = ref(i,1) + ax(i) * (X_t1_t(1) - s0(i));
Y_t1_t(1,2) = ref(i,2) + ay(i) * (X_t1_t(1) - s0(i));
Y_t1_t(1,3) = ref(i,3) + az(i) * (X_t1_t(1) - s0(i));

%cartesian coordinates position
if ((Y_t1_t(1) < min(ref(i,1),ref(i+1,1))) | (Y_t1_t(1) > max(ref(i,1),ref(i+1,1))) | ...
    (Y_t1_t(2) < min(ref(i,2),ref(i+1,2))) | (Y_t1_t(2) > max(ref(i,2),ref(i+1,2))) | ...
    (Y_t1_t(3) < min(ref(i,3),ref(i+1,3))) | (Y_t1_t(3) > max(ref(i,3),ref(i+1,3))))

    d1 = sqrt(sum((ref(i,:) - Y_t1_t).^2));
    d2 = sqrt(sum((ref(i+1,:) - Y_t1_t).^2));

    if (d1 < d2)
        X_t1_t(1) = s0(i);
        Y_t1_t(:) = ref(i,:);
    else
        X_t1_t(1) = s0(i+1);
        Y_t1_t(:) = ref(i+1,:);
    end
end