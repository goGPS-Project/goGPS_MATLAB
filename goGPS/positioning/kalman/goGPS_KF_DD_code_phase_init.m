function [kalman_initialized] = goGPS_KF_DD_code_phase_init(XR0, XM, time_rx, pr1_R, pr1_M, ...
         ph1_R, ph1_M, dop1_R, dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
         dop2_R, dop2_M, snr_R, snr_M, Eph, SP3, iono, lambda, frequencies, p_rate, dtMdot, flag_IAR, flag_XR, flag_tropo, sbas)

% SYNTAX:
%   [kalman_initialized] = goGPS_KF_DD_code_phase_init(XR0, XM, time_rx, pr1_R, pr1_M, ...
%        ph1_R, ph1_M, dop1_R, dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%        dop2_R, dop2_M, snr_R, snr_M, Eph, SP3, iono, lambda, frequencies, p_rate, dtMdot, flag_IAR, flag_XR, flag_tropo, sbas);
%
% INPUT:
%   XR0 = rover approximate/apriori position (X,Y,Z)
%   XM  = master known position (X,Y,Z)
%   time_rx = GPS time
%   pr1_R  = ROVER-SATELLITE code pseudorange (carrier L1)
%   pr1_M  = MASTER-SATELLITE code pseudorange (carrier L1)
%   ph1_R  = ROVER-SATELLITE phase observation (carrier L1)
%   ph1_M  = MASTER-SATELLITE phase observation (carrier L1)
%   dop1_R = ROVER-SATELLITE Doppler observation (carrier L1)
%   dop1_M = MASTER-SATELLITE Doppler observation (carrier L1)
%   pr2_R  = ROVER-SATELLITE code pseudorange (carrier L2)
%   pr2_M  = MASTER-SATELLITE code pseudorange (carrier L2)
%   ph2_R  = ROVER-SATELLITE phase observation (carrier L2)
%   ph2_M  = MASTER-SATELLITE phase observation (carrier L2)
%   dop2_R = ROVER-SATELLITE Doppler observation (carrier L2)
%   dop2_M = MASTER-SATELLITE Doppler observation (carrier L2)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   Eph = satellites ephemerides
%   SP3 = structure containing precise ephemeris data
%   iono =  ionospheric parameters (vector of zeroes if not available)
%   lambda = wavelength matrix (depending on the enabled constellations)
%   frequencies = L1 carrier (phase=1), L2 carrier (phase=2), L1&L2 (phase=[1 2])
%   p_rate = processing interval [s]
%   dtMdot = master receiver clock drift
%   flag_IAR = boolean variable to enable/disable integer ambiguity resolution
%   flag_XR  = 0: unknown
%              1: approximated
%              2: fixed
%   flag_tropo = boolean variable to enable/disable tropospheric delay estimation
%   sbas = SBAS corrections
%
% OUTPUT:
%   kalman_initialized = flag to point out whether Kalman has been successfully initialized
%
% DESCRIPTION:
%   Kalman filter initialization with the computation of the ROVER
%   initial position (X,Y,Z).

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     Andrea Nardo, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------
% KALMAN FILTER PARAMETERS
%--------------------------------------------------------------------------------------------

global sigmaq0 sigmaq0_N sigmaq0_tropo zero_time
global cutoff snr_threshold cond_num_threshold o1 o2 o3 nN nT
global n_sys

global Xhat_t_t X_t1_t T I Cee conf_sat conf_cs pivot pivot_old interval
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global doppler_pred_range1_R doppler_pred_range2_R
global doppler_pred_range1_M doppler_pred_range2_M
global ratiotest mutest succ_rate fixed_solution

kalman_initialized = 0;

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1_R,1);

%topocentric coordinates initialization
azR = zeros(nSatTot,1);
elR = zeros(nSatTot,1);
distR = zeros(nSatTot,1);
azM = zeros(nSatTot,1);
elM = zeros(nSatTot,1);
distM = zeros(nSatTot,1);

%compute inter-frequency factors (for the ionospheric delay)
ionoFactor = goGNSS.getInterFreqIonoFactor(lambda);

%--------------------------------------------------------------------------------------------
% SELECTION SINGLE / DUAL FREQUENCY
%--------------------------------------------------------------------------------------------

%number of unknown phase ambiguities
if (length(frequencies) == 1)
    nN = nSatTot;
else
    nN = nSatTot*2;
end

%--------------------------------------------------------------------------------------------
% NUMBER OF TROPOSPHERIC PARAMETERS
%--------------------------------------------------------------------------------------------

nT = 2;

%--------------------------------------------------------------------------------------------
% DYNAMIC MODEL OF THE KALMAN FILTER
%--------------------------------------------------------------------------------------------

%zeroes vectors useful in the matrices definition
Z_nN_o1 = zeros(nN,o1);
Z_o1_nN = zeros(o1,nN);
Z_nT_nN = zeros(nT,nN);
Z_nN_nT = zeros(nN,nT);
Z_nT_o1 = zeros(nT,o1);
Z_o1_nT = zeros(o1,nT);
Z_o1_o1 = zeros(o1);

%T matrix construction - system dynamics
%position and velocity equations
T0 = eye(o1) + diag(ones(o1-1,1),1)*interval;

%second degree polynomial
% T0 = [1 1; 0 1];
%third degree polynomial
% T0 = [1 1 0; 0 1 1; 0 0 1]

%matrix structure of initial comb_N
N0 = eye(nN);

%matrix structure of tropospheric parameters
TT = eye(nT);

%system dynamics
%X(t+1)  = X(t) + Vx(t)
%Vx(t+1) = Vx(t)
%... <-- for the other two variables Y e Z
%comb_N(t+1) = comb_N(t)

T = [T0      Z_o1_o1 Z_o1_o1 Z_o1_nN Z_o1_nT;
     Z_o1_o1 T0      Z_o1_o1 Z_o1_nN Z_o1_nT;
     Z_o1_o1 Z_o1_o1 T0      Z_o1_nN Z_o1_nT;
     Z_nN_o1 Z_nN_o1 Z_nN_o1 N0      Z_nN_nT;
     Z_nT_o1 Z_nT_o1 Z_nT_o1 Z_nT_nN TT];

%construction of an identity matrix
I = eye(o3+nN+nT);

%--------------------------------------------------------------------------------------------
% SELECTION OF THE SATELLITES
%--------------------------------------------------------------------------------------------

if (length(frequencies) == 2)
    sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) & (pr2_R ~= 0) & (pr2_M ~= 0) );
    sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) & (ph1_R ~= 0) & (ph1_M ~= 0) & ...
                (pr2_R ~= 0) & (pr2_M ~= 0) & (ph2_R ~= 0) & (ph2_M ~= 0) );
else
    if (frequencies == 1)
        sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) );
        sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) & ...
                    (ph1_R ~= 0) & (ph1_M ~= 0) );
    else
        sat_pr = find( (pr2_R ~= 0) & (pr2_M ~= 0) );
        sat = find( (pr2_R ~= 0) & (pr2_M ~= 0) & ...
                    (ph2_R ~= 0) & (ph2_M ~= 0) );
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

%------------------------------------------------------------------------------------
% APPROXIMATE POSITION
%-----------------------------------------------------------------------------------

if ((sum(abs(XR0)) == 0) || isempty(XR0))
    %approximate position not available
    flag_XR = 0;
    XR0 = [];
% else
%     %approximate position available
%     flag_XR = 1;
end

%--------------------------------------------------------------------------------------------
% KALMAN FILTER INITIAL STATE
%--------------------------------------------------------------------------------------------

%zeroes vectors useful in the matrices definition
Z_om_1 = zeros(o1-1,1);
sigma2_N = zeros(nN,1);

min_nsat_LS = 3 + n_sys;

if (length(sat_pr) >= min_nsat_LS)

    sat_pr_old = sat_pr;

    if (frequencies(1) == 1)
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono1_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr1_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, sbas,  XM,  [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, frequencies, p_rate,      2, 0); %#ok<ASGLU>
        if (length(sat_pr_M) < min_nsat_LS); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono1_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, sbas,  XR0, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 1); %#ok<ASGLU>
    else
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono1_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr2_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, sbas,  XM,  [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, frequencies, p_rate,       2, 0); %#ok<ASGLU>
        if (length(sat_pr_M) < min_nsat_LS); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono1_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, sbas,  XR0, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 1); %#ok<ASGLU>
    end

    err_iono2_M = err_iono1_M .* ionoFactor(sat_pr_M,2);
    err_iono2_R = err_iono1_R .* ionoFactor(sat_pr_R,2);

%     if flag_XR==2
%         cov_XR=eye(3).*(0.01^2);  % da sistemare pi? accurata
%         cond_num=0;
%     end

    %keep only satellites that rover and master have in common
    [sat_pr, iR, iM] = intersect(sat_pr_R, sat_pr_M);
    err_tropo_R = err_tropo_R(iR);
    err_tropo_M = err_tropo_M(iM);
    err_iono1_R = err_iono1_R(iR);
    err_iono1_M = err_iono1_M(iM);
    err_iono2_R = err_iono2_R(iR);
    err_iono2_M = err_iono2_M(iM);

    %apply cutoffs also to phase satellites
    sat_removed = setdiff(sat_pr_old, sat_pr);
    sat(ismember(sat,sat_removed)) = [];

    for i = 1:size(sat_pr)
        if (nargin > 23 && ~isempty(dtMdot) && dop1_M(sat_pr(i)) == 0 && any(Eph(:)))
            [dop1_M(sat_pr(i)), dop2_M(sat_pr(i))] = doppler_shift_approx(XM, zeros(3,1), XS_tx(i,:)', VS_tx(i,:)', time_tx(i), dtMdot, sat_pr(i), Eph, lambda(sat_pr(i),:));
        end
    end

    %--------------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
    %--------------------------------------------------------------------------------------------

    %satellites configuration: code only (-1), both code and phase (+1);
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat_pr) = -1;
    conf_sat(sat) = +1;

    %cycle-slip configuration (no cycle-slip)
    conf_cs = zeros(nSatTot,1);

    %previous pivot
    pivot_old = 0;

    %current pivot
    if ~isempty(sat)
        [null_max_elR, pivot_index] = max(elR(sat)); %#ok<ASGLU>
        pivot = sat(pivot_index);
    else
        [null_max_elR, pivot_index] = max(elR(sat_pr)); %#ok<ASGLU>
        pivot = sat_pr(pivot_index);
    end

    %if at least min_nsat_LS satellites are available after the cutoffs, and if the
    % condition number in the least squares does not exceed the threshold
    if (size(sat_pr,1) >= min_nsat_LS && cond_num < cond_num_threshold)

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

%do not use least squares ambiguity estimation
% NOTE: LS amb. estimation is automatically switched off if the number of
% satellites with phase available is not sufficient
if (size(sat_pr,1) + size(sat,1) - 2 <= 3 + size(sat,1) - 1 || size(sat,1) <= min_nsat_LS)

    %ambiguity initialization: initialized value
    %if the satellite is visible, 0 if the satellite is not visible
    N1 = zeros(nSatTot,1);
    N2 = zeros(nSatTot,1);
    sigma2_N1 = zeros(nSatTot,1);
    sigma2_N2 = zeros(nSatTot,1);

    %computation of the phase double differences in order to estimate N
    if ~isempty(sat)
        [N1(sat), sigma2_N1(sat)] = amb_estimate_observ(pr1_R(sat), pr1_M(sat), ph1_R(sat), ph1_M(sat), pivot, sat, lambda(sat,1));
        [N2(sat), sigma2_N2(sat)] = amb_estimate_observ(pr2_R(sat), pr2_M(sat), ph2_R(sat), ph2_M(sat), pivot, sat, lambda(sat,2));
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

    %find the indices of the satellites with phase available
    [~, index] = intersect(sat_pr,sat);

    % if fixed apriori rover coordinates are available, estimation of
    % ambiguities from geometrical range
    if flag_XR==2
        %frequency 1
        [N1(sat)]=ambiguity_init_from_range(distR(sat), distM(sat), pivot_index, ph1_R(sat), ph1_M(sat), lambda(sat,1));
        cov_N1=diag(sum(sigma2_XR) ./ (lambda(sat).^2));
        N=N1;
        sigma2_N(sat) = diag(cov_N1);

    else

        %ROVER positioning improvement with code and phase double differences
        if ~isempty(sat)
            [     XR, N1(sat),      cov_XR, cov_N1, PDOP, HDOP, VDOP] = LS_DD_code_phase(XR, XM, XS(index,:), pr1_R(sat), ph1_R(sat), snr_R(sat), pr1_M(sat), ph1_M(sat), snr_M(sat), elR(sat), elM(sat), err_tropo_R(index), err_iono1_R(index), err_tropo_M(index), err_iono1_M(index), pivot_index, lambda(sat,1), 0);
            [null_XR, N2(sat), null_cov_XR, cov_N2]                   = LS_DD_code_phase(XR, XM, XS(index,:), pr2_R(sat), ph2_R(sat), snr_R(sat), pr2_M(sat), ph2_M(sat), snr_M(sat), elR(sat), elM(sat), err_tropo_R(index), err_iono2_R(index), err_tropo_M(index), err_iono2_M(index), pivot_index, lambda(sat,2), 0); %#ok<ASGLU>
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

        if (length(frequencies) == 2)
            N = [N1; N2];
            sigma2_N(sat) = diag(cov_N1);
            sigma2_N(sat+nSatTot) = diag(cov_N2);
        else
            if (frequencies == 1)
                N = N1;
                sigma2_N(sat) = diag(cov_N1);
            else
                N = N1;
                sigma2_N(sat) = diag(cov_N2);
            end
        end
    end
end

%a-priori tropospheric delay
if (flag_tropo)

    [week, sow] = time2weektow(time_rx + zero_time);
    date = gps2date(week, sow);
    [~, mjd] = date2jd(date);

    [phi_R, lam_R, h_R] = cart2geod(XR(1), XR(2), XR(3));
    [phi_M, lam_M, h_M] = cart2geod(XM(1), XM(2), XM(3));

    %ZTD = saast_dry(goGNSS.STD_PRES, H, phi) + saast_wet(goGNSS.STD_TEMP, H); %H here is orthometric

    [pressure_R, temperature_R, undu_R] = gpt(mjd, phi_R, lam_R, h_R); %#ok<ASGLU>
    ZWD_R = saast_wet(temperature_R, goGNSS.STD_HUMI, h_R - undu_R);

    [pressure_M, temperature_M, undu_M] = gpt(mjd, phi_M, lam_M, h_M); %#ok<ASGLU>
    ZWD_M = saast_wet(temperature_M, goGNSS.STD_HUMI, h_M - undu_M);
else
    ZWD_R = 0;
    ZWD_M = 0;
end

%initialization of the state vector
Xhat_t_t = [XR(1); Z_om_1; XR(2); Z_om_1; XR(3); Z_om_1; N; ZWD_R; ZWD_M];

%update at step t+1 X Vx Y Vy Z Vz N
X_t1_t = T*Xhat_t_t;

%--------------------------------------------------------------------------------------------
% INITIAL STATE COVARIANCE MATRIX
%--------------------------------------------------------------------------------------------

%initial state covariance matrix
Cee(:,:) = zeros(o3+nN+nT);
Cee(1,1) = sigma2_XR(1);
Cee(o1+1,o1+1) = sigma2_XR(2);
Cee(o2+1,o2+1) = sigma2_XR(3);
Cee(2:o1,2:o1) = sigmaq0 * eye(o1-1);
Cee(o1+2:o2,o1+2:o2) = sigmaq0 * eye(o1-1);
Cee(o2+2:o3,o2+2:o3) = sigmaq0 * eye(o1-1);
Cee(o3+1:o3+nN,o3+1:o3+nN) = diag(sigma2_N);
Cee(o3+nN+1:o3+nN+nT,o3+nN+1:o3+nN+nT) = sigmaq0_tropo * eye(nT);

%--------------------------------------------------------------------------------------------
% INTEGER AMBIGUITY SOLVING BY LAMBDA METHOD
%--------------------------------------------------------------------------------------------

sat_np = sat(sat ~= pivot);
if (flag_IAR && ~isempty(sat_np))
    %try to solve integer ambiguities
    [Xhat_t_t([1 o1+1 o2+1]), Xhat_t_t(o3+sat_np)] = lambdafix(Xhat_t_t([1 o1+1 o2+1]), Xhat_t_t(o3+sat_np), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), Cee(o3+sat_np,o3+sat_np), Cee([1 o1+1 o2+1],o3+sat_np));
else
    ratiotest = [ratiotest NaN];
    mutest    = [mutest NaN];
    succ_rate = [succ_rate NaN];
    fixed_solution = [fixed_solution 0];
end

%--------------------------------------------------------------------------------------------
% DOPPLER-BASED PREDICTION OF PHASE RANGES
%--------------------------------------------------------------------------------------------
if (dop1_R(sat))
    doppler_pred_range1_R(sat,1) = ph1_R(sat) - dop1_R(sat);
end
if (dop2_R(sat))
    doppler_pred_range2_R(sat,1) = ph2_R(sat) - dop2_R(sat);
end
if (dop1_M(sat))
    doppler_pred_range1_M(sat,1) = ph1_M(sat) - dop1_M(sat);
end
if (dop2_M(sat))
    doppler_pred_range2_M(sat,1) = ph2_M(sat) - dop2_M(sat);
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
