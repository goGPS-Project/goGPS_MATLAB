function [kalman_initialized] = goGPS_KF_SA_code_phase_init(XR0, time_rx, pr1, ph1, dop1, pr2, ph2, dop2, snr, Eph, SP3, iono, sbas, lambda, frequencies, obs_comb, flag_IAR, flag_XR, flag_tropo)

% SYNTAX:
%   [kalman_initialized] = goGPS_KF_SA_code_phase_init(XR0, time_rx, pr1, ph1, dop1, pr2, ph2, dop2, snr, Eph, SP3, iono, sbas, lambda, frequencies, obs_comb, flag_IAR, flag_XR, flag_tropo);
%
% INPUT:
%   pos_R = rover approximate coordinates (X, Y, Z)
%   time_rx = GPS time
%   pr1  = ROVER-SATELLITE code pseudorange (L1 carrier)
%   ph1  = ROVER-SATELLITE phase observation (carrier L1)
%   dop1 = ROVER_SATELLITE Doppler observation (carrier L1)
%   pr2  = ROVER-SATELLITE code pseudorange (L2 carrier)
%   ph2  = ROVER-SATELLITE phase observation (carrier L2)
%   dop2 = ROVER_SATELLITE Doppler observation (carrier L2)
%   snr  = ROVER-SATELLITE signal-to-noise ratio
%   Eph  = satellite ephemerides
%   SP3 = structure containing precise ephemeris data
%   iono =  ionospheric parameters (vector of zeroes if not available)
%   sbas = SBAS corrections
%   lambda = wavelength matrix (depending on the enabled constellations)
%   frequencies = L1 carrier (phase=1) L2 carrier (phase=2)
%   obs_comb = observations combination (e.g. iono-free: obs_comb = 'IONO_FREE')
%   flag_IAR = boolean variable to enable/disable integer ambiguity resolution
%   flag_XR  = 0: unknown
%              1: approximated
%              2: fixed
%   flag_tropo = boolean variable to enable/disable tropospheric delay estimation
%
% OUTPUT:
%   kalman_initialized = flag to point out whether Kalman has been successfully initialized
%
% DESCRIPTION:
%   Standalone phase and code Kalman filter initialization.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Andrea Nardo
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

global sigmaq0 sigmaq0_N sigmaq0_tropo zero_time
global cutoff snr_threshold cond_num_threshold o1 o2 o3 nN nT nC

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old interval
global azR elR distR azM elM distM phwindup
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global doppler_pred_range1_R doppler_pred_range2_R
global ratiotest mutest succ_rate fixed_solution
global n_sys

kalman_initialized = 0;

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1,1);

%compute inter-frequency factors (for the ionospheric delay)
ionoFactor = goGNSS.getInterFreqIonoFactor(lambda);

%iono-free coefficients
alpha1 = (goGNSS.F1^2/(goGNSS.F1^2 - goGNSS.F2^2));
alpha2 = (goGNSS.F2^2/(goGNSS.F1^2 - goGNSS.F2^2));
alphat = 77;
alphan = 60;
lambdaIF = alphat/(alphat^2-alphan^2)*lambda(:,1);

%topocentric coordinates initialization
azR = zeros(nSatTot,1);
elR = zeros(nSatTot,1);
distR = zeros(nSatTot,1);
azM = zeros(nSatTot,1);
elM = zeros(nSatTot,1);
distM = zeros(nSatTot,1);

%phase wind-up matrix initialization
phwindup = zeros(nSatTot,1);

%--------------------------------------------------------------------------------------------
% SELECTION SINGLE / DUAL FREQUENCY
%--------------------------------------------------------------------------------------------

%number of unknown phase ambiguities
if (length(frequencies) == 2 && strcmp(obs_comb,'NONE'))
    nN = nSatTot*2;
else
    nN = nSatTot;
end

%--------------------------------------------------------------------------------------------
% NUMBER OF TROPOSPHERIC PARAMETERS
%--------------------------------------------------------------------------------------------

nT = 1;

%--------------------------------------------------------------------------------------------
% NUMBER OF CLOCK PARAMETERS
%--------------------------------------------------------------------------------------------

nC = 1;

%--------------------------------------------------------------------------------------------
% KALMAN FILTER DYNAMIC MODEL
%--------------------------------------------------------------------------------------------

%zeroes vectors useful in the matrices definition
Z_nN_o1 = zeros(nN,o1);
Z_o1_nN = zeros(o1,nN);
Z_nT_nN = zeros(nT,nN);
Z_nN_nT = zeros(nN,nT);
Z_nC_nN = zeros(nC,nN);
Z_nN_nC = zeros(nN,nC);
Z_nT_o1 = zeros(nT,o1);
Z_o1_nT = zeros(o1,nT);
Z_nC_o1 = zeros(nC,o1);
Z_o1_nC = zeros(o1,nC);
Z_nT_nC = zeros(nT,nC);
Z_nC_nT = zeros(nC,nT);
Z_o1_o1 = zeros(o1);

%T matrix construction - system dynamics
%position and velocity equations
D0 = eye(o1) + diag(ones(o1-1,1),1)*interval;

%second degree polynomial
% T0 = [1 1; 0 1];
%third degree polynomial
% T0 = [1 1 0; 0 1 1; 0 0 1]

%matrix structure of initial comb_N
N0 = eye(nN);

%matrix structure of tropospheric parameters
T0 = eye(nT);

%matrix structure of clock parameters
C0 = eye(nC);

%system dynamics
%X(t+1)  = X(t) + Vx(t)
%Vx(t+1) = Vx(t)
%... <-- for the other two variables Y e Z
%comb_N(t+1) = comb_N(t)

T = [D0      Z_o1_o1 Z_o1_o1 Z_o1_nN Z_o1_nT Z_o1_nC;
     Z_o1_o1 D0      Z_o1_o1 Z_o1_nN Z_o1_nT Z_o1_nC;
     Z_o1_o1 Z_o1_o1 D0      Z_o1_nN Z_o1_nT Z_o1_nC;
     Z_nN_o1 Z_nN_o1 Z_nN_o1 N0      Z_nN_nT Z_nN_nC;
     Z_nT_o1 Z_nT_o1 Z_nT_o1 Z_nT_nN T0      Z_nT_nC
     Z_nC_o1 Z_nC_o1 Z_nC_o1 Z_nC_nN Z_nC_nT C0];

%construction of an identity matrix
I = eye(o3+nN+nT+nC);

%--------------------------------------------------------------------------------------------
% SATELLITE SELECTION
%--------------------------------------------------------------------------------------------

if (length(frequencies) == 2)
    sat_pr = find( (pr1 ~= 0) & (pr2 ~= 0) );
    sat = find( (pr1 ~= 0) & (ph1 ~= 0) & (pr2 ~= 0) & (ph2 ~= 0) );
else
    if (frequencies == 1)
        sat_pr = find( (pr1 ~= 0) );
        sat = find( (pr1 ~= 0) & (ph1 ~= 0) );
    else
        sat_pr = find( (pr2 ~= 0) );
        sat = find( (pr2 ~= 0) & (ph2 ~= 0) );
    end
end
if (isempty(SP3))
    eph_avail = Eph(30,:);
else
    eph_avail = SP3.avail;
end
sat_pr = sat_pr(ismember(sat_pr, eph_avail));
sat = sat(ismember(sat, eph_avail));

%only satellites with code and phase
%sat_pr = sat;

%--------------------------------------------------------------------------------------------
% SBAS FAST CORRECTIONS
%--------------------------------------------------------------------------------------------

if (~isempty(sbas))
    %apply SBAS fast (pseudorange) corrections
    pr1(sat_pr) = pr1(sat_pr) + sbas.prc(sat_pr)';
end

%------------------------------------------------------------------------------------
% APPROXIMATE POSITION
%-----------------------------------------------------------------------------------

if ((sum(abs(XR0)) == 0) || isempty(XR0))
    %approximate position not available
    flag_XR = 0;
% else
%     %approximate position available
%     flag_XR = 1;
end

%--------------------------------------------------------------------------------------------
% KALMAN FILTER INITIAL STATE
%--------------------------------------------------------------------------------------------

%zero vector useful in matrix definitions
Z_om_1 = zeros(o1-1,1);
sigma2_N = zeros(nN,1);

min_nsat_LS = 3 + n_sys;

if (length(sat_pr) >= min_nsat_LS)

    sat_pr_old = sat_pr;

    if (frequencies(1) == 1)
        if (length(frequencies) < 2 || ~strcmp(obs_comb,'IONO_FREE'))
            [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono1, sat_pr, elR(sat_pr), azR(sat_pr), distR(sat_pr), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1(sat_pr), snr(sat_pr), Eph, SP3, iono, sbas, XR0, [], [], sat_pr, [], lambda(sat_pr,:), cutoff, snr_threshold, frequencies, flag_XR, 0); %#ok<ASGLU>
        else
            [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono1, sat_pr, elR(sat_pr), azR(sat_pr), distR(sat_pr), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, alpha1*pr1(sat_pr) - alpha2*pr2(sat_pr), snr(sat_pr), Eph, SP3, zeros(8,1), sbas, XR0, [], [], sat_pr, [], zeros(length(sat_pr),2), cutoff, snr_threshold, frequencies, flag_XR, 0); %#ok<ASGLU>
        end
        
        err_iono2 = err_iono1 .* ionoFactor(sat_pr,2);
    else
        [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono2, sat_pr, elR(sat_pr), azR(sat_pr), distR(sat_pr), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2(sat_pr), snr(sat_pr), Eph, SP3, iono, sbas, XR0, [], [], sat_pr, [], lambda(sat_pr,:), cutoff, snr_threshold, frequencies, flag_XR, 0); %#ok<ASGLU>
        
        err_iono1 = err_iono2 ./ ionoFactor(sat_pr,2);
    end

    %apply cutoffs also to phase satellites
    sat_removed = setdiff(sat_pr_old, sat_pr);
    sat(ismember(sat,sat_removed)) = [];
    
    %if multi-system observations, then an additional parameter to estimate the inter-system bias
    %for each additional system is needed
    uni_sys = unique(sys(sys ~= 0));
    num_sys = length(uni_sys);
    min_nsat = 3 + num_sys;

    %--------------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING
    %--------------------------------------------------------------------------------------------
    
    %satellites configuration: code only (-1), both code and phase (+1);
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat_pr) = -1;
    conf_sat(sat) = +1;
    
    %cycle-slip configuration
    conf_cs = zeros(nSatTot,1);
    
    pivot_old = 0;
    
    %current pivot
    [null_max_elR, i] = max(elR(sat_pr)); %#ok<ASGLU>
    pivot = sat_pr(i);
    
    %if the number of satellites is not sufficient after the cutoffs, and
    %if the condition number in the least squares exceeds the threshold
    if (size(sat_pr,1) >= min_nsat && cond_num < cond_num_threshold)
        
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

%compute phase wind-up correction
phwindup(sat,1) = phase_windup_correction(time_rx, XR, XS, SP3, phwindup(sat,1));

%do not use least squares ambiguity estimation
% NOTE: LS amb. estimation is automatically switched off if the number of
% satellites with phase available is not sufficient
if (length(sat) < min_nsat)
    
    %ambiguity initialization: initialized value
    %if the satellite is visible, 0 if the satellite is not visible
    N1 = zeros(nSatTot,1);
    N2 = zeros(nSatTot,1);
    sigma2_N1 = zeros(nSatTot,1);
    sigma2_N2 = zeros(nSatTot,1);
    
    %computation of the phase double differences in order to estimate N
    if ~isempty(sat)
        [N1(sat), sigma2_N1(sat)] = amb_estimate_observ_SA(pr1(sat), ph1(sat), lambda(sat,1));
        [N2(sat), sigma2_N2(sat)] = amb_estimate_observ_SA(pr2(sat), ph2(sat), lambda(sat,2));
    end

    if (length(frequencies) == 2)
        N = [N1; N2];
        sigma2_N = [sigma2_N1; sigma2_N2];
    else
        if (frequencies == 1)
            N = N1;
            sigma2_N = sigma2_N1;
        else
            N = N2;
            sigma2_N = sigma2_N2;
        end
    end

%use least squares ambiguity estimation
else
    
    %ambiguity initialization: initialized value
    %if the satellite is visible, 0 if the satellite is not visible
    N1 = zeros(nSatTot,1);
    N2 = zeros(nSatTot,1);
    N_IF = zeros(nSatTot,1);
    cov_N1 = [];
    cov_N2 = [];
    cov_N_IF = [];

    %ROVER positioning improvement with code and phase double differences
    if ~isempty(sat)
        if (~strcmp(obs_comb,'IONO_FREE'))
            [XR, dtR, N1(sat), cov_XR, var_dtR, cov_N1, PDOP, HDOP, VDOP] = LS_SA_code_phase(XR, XS, pr1(sat_pr), ph1(sat_pr), snr(sat_pr), elR(sat_pr), distR(sat_pr), sat_pr, sat, dtS, err_tropo, err_iono1, phwindup(sat_pr), sys, lambda(sat_pr,1));
            [ ~,   ~, N2(sat),      ~,       ~, cov_N2]                   = LS_SA_code_phase(XR, XS, pr2(sat_pr), ph2(sat_pr), snr(sat_pr), elR(sat_pr), distR(sat_pr), sat_pr, sat, dtS, err_tropo, err_iono2, phwindup(sat_pr), sys, lambda(sat_pr,2));
        else
            [XR, dtR, N_IF(sat), cov_XR, var_dtR, cov_N_IF, PDOP, HDOP, VDOP] = LS_SA_code_phase(XR, XS, alpha1*pr1(sat_pr) - alpha2*pr2(sat_pr), alphat*ph1(sat_pr) - alphan*ph2(sat_pr), snr(sat_pr), elR(sat_pr), distR(sat_pr), sat_pr, sat, dtS, err_tropo, zeros(size(sat_pr)), phwindup(sat_pr), sys, lambdaIF(sat_pr,1));
        end
    end
    
    if isempty(cov_XR) %if it was not possible to compute the covariance matrix
        cov_XR = sigmaq0 * eye(3);
    end
    sigma2_XR = diag(cov_XR);
    
    if isempty(cov_N1) %if it was not possible to compute the covariance matrix
        cov_N1 = sigmaq0_N * eye(length(sat));
    end
    
    if isempty(cov_N2) %if it was not possible to compute the covariance matrix
        cov_N2 = sigmaq0_N * eye(length(sat));
    end
    
    if isempty(cov_N_IF) %if it was not possible to compute the covariance matrix
        cov_N_IF = sigmaq0_N * eye(length(sat));
    end
    
    if (length(frequencies) == 2)
        if (strcmp(obs_comb,'NONE'))
            N = [N1; N2];
            sigma2_N(sat) = diag(cov_N1);
            %sigma2_N(sat) = (sigmaq_cod1 / lambda1^2) * ones(length(sat),1);
            sigma2_N(sat+nSatTot) = diag(cov_N2);
            %sigma2_N(sat+nSatTot) = (sigmaq_cod2 / lambda2^2) * ones(length(sat),1);
        elseif (strcmp(obs_comb,'IONO_FREE'))
            N = N_IF;
            sigma2_N(sat) = diag(cov_N_IF);
        end
    else
        if (frequencies == 1)
            N = N1;
            sigma2_N(sat) = diag(cov_N1);
            %sigma2_N(sat) = (sigmaq_cod1 / lambda1^2) * ones(length(sat),1);
        else
            N = N2;
            sigma2_N(sat) = diag(cov_N2);
            %sigma2_N(sat) = (sigmaq_cod2 / lambda2^2) * ones(length(sat),1);
        end
    end
end

%a-priori tropospheric delay
if (flag_tropo)
    [week, sow] = time2weektow(time_rx + zero_time);
    date = gps2date(week, sow);
    [~, mjd] = date2jd(date);
    
    [phi_R, lam_R, h_R] = cart2geod(XR(1), XR(2), XR(3));
    [pres_R, temp_R, undu_R] = gpt(mjd, phi_R, lam_R, h_R);
    ZHD_R = saast_dry(pres_R, h_R - undu_R, phi_R*180/pi);
    ZWD_R = saast_wet(temp_R, goGNSS.STD_HUMI, h_R - undu_R);
else
    ZWD_R = 0;
end

%initialization of the state vector
Xhat_t_t = [XR(1); Z_om_1; XR(2); Z_om_1; XR(3); Z_om_1; N; ZWD_R; goGNSS.V_LIGHT*dtR];

%state update at step t+1 X Vx Y Vy Z Vz comb_N
%estimation at step t, because the initial velocity is equal to 0
X_t1_t = T*Xhat_t_t;

%--------------------------------------------------------------------------------------------
% RECONSTRUCTION OF FULL ZTD
%--------------------------------------------------------------------------------------------
if (flag_tropo)
    Xhat_t_t(o3+nN+(1:nT)) = ZHD_R + Xhat_t_t(o3+nN+(1:nT));
end

%--------------------------------------------------------------------------------------------
% INITIAL STATE COVARIANCE MATRIX
%--------------------------------------------------------------------------------------------

%initial state covariance matrix
Cee(:,:) = zeros(o3+nN+nT+nC);
Cee(1,1) = sigma2_XR(1);
Cee(o1+1,o1+1) = sigma2_XR(2);
Cee(o2+1,o2+1) = sigma2_XR(3);
Cee(2:o1,2:o1) = sigmaq0 * eye(o1-1);
Cee(o1+2:o2,o1+2:o2) = sigmaq0 * eye(o1-1);
Cee(o2+2:o3,o2+2:o3) = sigmaq0 * eye(o1-1);
Cee(o3+1:o3+nN,o3+1:o3+nN) = diag(sigma2_N);
Cee(o3+nN+1:o3+nN+nT,o3+nN+1:o3+nN+nT) = sigmaq0_tropo * eye(nT);
Cee(o3+nN+nT+1:o3+nN+nT+nC,o3+nN+nT+1:o3+nN+nT+nC) = goGNSS.V_LIGHT^2*var_dtR * eye(nC);

%--------------------------------------------------------------------------------------------
% INTEGER AMBIGUITY SOLVING BY LAMBDA METHOD
%--------------------------------------------------------------------------------------------

% if (flag_IAR && ~isempty(sat))
%     %try to solve integer ambiguities
%     [Xhat_t_t([1 o1+1 o2+1]), Xhat_t_t(o3+sat)] = lambdafix(Xhat_t_t([1 o1+1 o2+1]), Xhat_t_t(o3+sat), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), Cee(o3+sat,o3+sat), Cee([1 o1+1 o2+1],o3+sat));
% else
    ratiotest = [ratiotest NaN];
    mutest    = [mutest NaN];
    succ_rate = [succ_rate NaN];
    fixed_solution = [fixed_solution 0];
% end

%--------------------------------------------------------------------------------------------
% DOPPLER-BASED PREDICTION OF PHASE RANGES
%--------------------------------------------------------------------------------------------
if (dop1(sat))
    doppler_pred_range1_R(sat,1) = ph1(sat) - dop1(sat);
end
if (dop2(sat))
    doppler_pred_range2_R(sat,1) = ph2(sat) - dop2(sat);
end

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
