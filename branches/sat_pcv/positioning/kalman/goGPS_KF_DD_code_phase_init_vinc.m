function [kalman_initialized] = goGPS_KF_DD_code_phase_init_vinc(XR0, XM, time_rx, pr1_R, pr1_M, ...
         ph1_R, ph1_M, dop1_R, dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
         dop2_R, dop2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase, ref, dtMdot)

% SYNTAX:
%   [kalman_initialized] = goGPS_KF_DD_code_phase_init(XR0, XM, time_rx, pr1_R, pr1_M, ...
%        ph1_R, ph1_M, dop1_R, dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%        dop2_R, dop2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase, ref, dtMdot);
%
% INPUT:
%   XR0 = rover approximate position (X,Y,Z)
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
%   phase = carrier L1 (phase=1) carrier L2 (phase=2)
%   ref = reference trajectory
%   dtMdot = master receiver clock drift
%
% OUTPUT:
%   kalman_initialized = flag to point out whether Kalman has been successfully initialized
%
% DESCRIPTION:
%   Kalman filter initialization with the computation of the ROVER
%   initial position (X,Y,Z). Constrained path.

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

%--------------------------------------------------------------------------------------------
% KALMAN FILTER PARAMETERS
%--------------------------------------------------------------------------------------------

global sigmaq0 sigmaq_vel sigmaq0_N
global cutoff snr_threshold cond_num_threshold o1 nN
global s0 ax ay az

global Xhat_t_t X_t1_t Yhat_t_t Y_t1_t T I Cee conf_sat conf_cs pivot pivot_old interval
global azR elR distR azM elM distM
global PDOP HDOP VDOP
global doppler_pred_range1_R doppler_pred_range2_R
global doppler_pred_range1_M doppler_pred_range2_M
global n_sys

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
% SINGLE / DOUBLE FREQUENCY SELECTION
%--------------------------------------------------------------------------------------------

%number of unknown phase ambiguities
if (length(phase) == 1)
    nN = nSatTot;
else
    nN = nSatTot*2;
end

%--------------------------------------------------------------------------------------------
% CONSTRAINED PATH
%--------------------------------------------------------------------------------------------

%angular coefficients in the three dimensions of hortogonal axes estimate
ax = ref(2:end,1) - ref(1:end-1,1);
ay = ref(2:end,2) - ref(1:end-1,2);
az = ref(2:end,3) - ref(1:end-1,3);

%normalization on the segment distance
ad = sqrt(ax.^2 + ay.^2 + az.^2);
ax = ax ./ ad;
ay = ay ./ ad;
az = az ./ ad;

%curvilinear coordinate offset
s0 = [0; cumsum(ad)];

%--------------------------------------------------------------------------------------------
% DYNAMIC MODEL OF THE KALMAN FILTER
%--------------------------------------------------------------------------------------------

%zeroes vectors useful in the matrices definition
Z_nN_o1 = zeros(nN,o1);
Z_o1_nN = zeros(o1,nN);

%T matrix construction - system dynamics
%position and velocity equations
T0 = eye(o1) + diag(ones(o1-1,1),1)*interval;

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

T = [T0      Z_o1_nN;
     Z_nN_o1 N0];

%construction of an identity matrix
I = eye(o1+nN);

%model error covariance matrix
Cvv = zeros(o1+nN);
Cvv(o1,o1) = sigmaq_vel; %#ok<NASGU>

%--------------------------------------------------------------------------------------------
% SELECTION OF THE SATELLITES
%--------------------------------------------------------------------------------------------

if (length(phase) == 2)
    sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) & (pr2_R ~= 0) & (pr2_M ~= 0) );
    sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) & (ph1_R ~= 0) & (ph1_M ~= 0) & ...
                (pr2_R ~= 0) & (pr2_M ~= 0) & (ph2_R ~= 0) & (ph2_M ~= 0) );
else
    if (phase == 1)
        sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) );
        sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) & ...
                    (ph1_R ~= 0) & (ph1_M ~= 0) );       
    else
        sat_pr = find( (pr2_R ~= 0) & (pr2_M ~= 0) );
          sat = find( (pr2_R ~= 0) & (pr2_M ~= 0) & ...
                      (ph2_R ~= 0) & (ph2_M ~= 0) );
    end
end
sat_pr = sat_pr(ismember(sat_pr, Eph(30,:)));
sat = sat(ismember(sat, Eph(30,:)));

%only satellites with code and phase
%sat_pr = sat;

%------------------------------------------------------------------------------------
% APPROXIMATE POSITION
%-----------------------------------------------------------------------------------

if ((sum(abs(XR0)) == 0) || isempty(XR0))
    %approximate position not available
    flag_XR = 0;
    XR0 = [];
else
    %approximate position available
    flag_XR = 1;
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
    
    if (phase(1) == 1)
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono1_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr1_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, [], XM,  [],  [], sat_pr,   [], lambda(sat_pr,:),   cutoff, snr_threshold, phase,       2, 0); %#ok<NASGU,ASGLU>
        if (length(sat_pr_M) < min_nsat_LS); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono1_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, [], XR0, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, phase, flag_XR, 1); %#ok<ASGLU>
        
        err_iono2_M = err_iono1_M .* ionoFactor(sat_pr_M,2);
        err_iono2_R = err_iono1_R .* ionoFactor(sat_pr_R,2);
    else
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono2_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr2_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, [], XM,  [],  [], sat_pr,   [], lambda(sat_pr,:),   cutoff, snr_threshold, phase,       2, 0); %#ok<NASGU,ASGLU>
        if (length(sat_pr_M) < min_nsat_LS); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono2_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, [], XR0, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, phase, flag_XR, 1); %#ok<ASGLU>
        
        err_iono1_M = err_iono2_M ./ ionoFactor(sat_pr_M,2);
        err_iono1_R = err_iono2_R ./ ionoFactor(sat_pr_R,2);
    end
    
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
        
        %projection over the constrain
        bx = XR(1) - ref(1:end-1,1) + ax.*s0(1:end-1);
        by = XR(2) - ref(1:end-1,2) + ay.*s0(1:end-1);
        bz = XR(3) - ref(1:end-1,3) + az.*s0(1:end-1);
        
        s_R = (ax.*bx + ay.*by + az.*bz) ./ (ax.^2 + ay.^2 + az.^2);
        
        XR_proj(:,1) = ref(1:end-1,1) + ax .* (s_R - s0(1:end-1));
        XR_proj(:,2) = ref(1:end-1,2) + ay .* (s_R - s0(1:end-1));
        XR_proj(:,3) = ref(1:end-1,3) + az .* (s_R - s0(1:end-1));
        
        %minimum distance estimate
        d = sqrt((XR(1) - XR_proj(:,1)).^2 + ...
            (XR(2) - XR_proj(:,2)).^2 + ...
            (XR(3) - XR_proj(:,3)).^2);
        
        [dmin, i] = min(d); %#ok<ASGLU>
        
        %cartesian coordinates positioning
        if ((XR_proj(i,1) >= min(ref(i,1),ref(i+1,1))) && (XR_proj(i,1) <= max(ref(i,1),ref(i+1,1))) && ...
                (XR_proj(i,2) >= min(ref(i,2),ref(i+1,2))) && (XR_proj(i,2) <= max(ref(i,2),ref(i+1,2))) && ...
                (XR_proj(i,3) >= min(ref(i,3),ref(i+1,3))) && (XR_proj(i,3) <= max(ref(i,3),ref(i+1,3))))
            
            s_R = s_R(i);
            XR = XR_proj(i,:)';
        else
            d = sqrt((XR(1) - ref(:,1)).^2 + ...
                (XR(2) - ref(:,2)).^2 + ...
                (XR(3) - ref(:,3)).^2);
            
            [dmin, i] = min(d); %#ok<ASGLU>
            
            s_R = s0(i);
            XR = ref(i,:)';
        end
        
        %propagated error
        sigma2_s_R = (ax(i)^2*sigma2_XR(1) + ay(i)^2*sigma2_XR(2) + az(i)^2*sigma2_XR(3)) ./ (ax(i)^2 + ay(i)^2 + az(i)^2)^2;
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

    if (length(phase) == 2)
        N = [N1; N2];
        sigma2_N = [sigma2_N1; sigma2_N2];
    else
        if (phase == 1)
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

    %ROVER positioning improvement with code and phase double differences
    if ~isempty(sat)
        [     XR, N1(sat),      cov_XR, cov_N1, PDOP, HDOP, VDOP] = LS_DD_code_phase(XR, XM, XS(index,:), pr1_R(sat), ph1_R(sat), snr_R(sat), pr1_M(sat), ph1_M(sat), snr_M(sat), elR(sat), elM(sat), err_tropo_R(index), err_iono1_R(index), err_tropo_M(index), err_iono1_M(index), pivot_index, lambda(sat,1), 0);
        [null_XR, N2(sat), null_cov_XR, cov_N2] =                   LS_DD_code_phase(XR, XM, XS(index,:), pr2_R(sat), ph2_R(sat), snr_R(sat), pr2_M(sat), ph2_M(sat), snr_M(sat), elR(sat), elM(sat), err_tropo_R(index), err_iono2_R(index), err_tropo_M(index), err_iono2_M(index), pivot_index, lambda(sat,2), 0); %#ok<ASGLU>
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
    
    if (length(phase) == 2)
        N = [N1; N2];
        sigma2_N(sat) = diag(cov_N1);
        sigma2_N(sat+nSatTot) = diag(cov_N2);
    else
        if (phase == 1)
            N = N1;
            sigma2_N(sat) = diag(cov_N1);
        else
            N = N1;
            sigma2_N(sat) = diag(cov_N2);
        end
    end
end

%initialization of the state vector
Xhat_t_t = [s_R; Z_om_1; N];

%state update at step t+1 X Vx Y Vy Z Vz N
%estimation at step t, because the initial velocity is equal to 0
X_t1_t = T*Xhat_t_t;

%--------------------------------------------------------------------------------------------
% CARTESIAN COORDINATES ESTIMATE
%--------------------------------------------------------------------------------------------

%curvilinear coordinate localization
i = find((Xhat_t_t(1) >= s0(1:end-1)) & (Xhat_t_t(1) < s0(2:end)));

%cartesian coordinates estimate
Yhat_t_t(1,1) = ref(i,1) + ax(i) * (Xhat_t_t(1) - s0(i));
Yhat_t_t(2,1) = ref(i,2) + ay(i) * (Xhat_t_t(1) - s0(i));
Yhat_t_t(3,1) = ref(i,3) + az(i) * (Xhat_t_t(1) - s0(i));

%curvilinear coordinate localization
i = find((X_t1_t(1) >= s0(1:end-1)) & (X_t1_t(1) < s0(2:end)));

%cartesian coordinates estimate
Y_t1_t(1,1) = ref(i,1) + ax(i) * (X_t1_t(1) - s0(i));
Y_t1_t(1,2) = ref(i,2) + ay(i) * (X_t1_t(1) - s0(i));
Y_t1_t(1,3) = ref(i,3) + az(i) * (X_t1_t(1) - s0(i));

%--------------------------------------------------------------------------------------------
% INITIAL STATE COVARIANCE MATRIX
%--------------------------------------------------------------------------------------------

%initial state covariance matrix
Cee(:,:) = zeros(o1+nN);
Cee(1,1) = sigma2_s_R;
Cee(2:o1,2:o1) = sigmaq0 * eye(o1-1);
Cee(o1+1:o1+nN,o1+1:o1+nN) = diag(sigma2_N);

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

% %--------------------------------------------------------------------------------------------
% % INITIAL KALMAN FILTER DOP
% %--------------------------------------------------------------------------------------------
% 
% %covariance propagation
% Cee_XYZ = Cee([1 o1+1 o2+1],[1 o1+1 o2+1]);
% Cee_ENU = global2localCov(Cee_XYZ, Xhat_t_t([1 o1+1 o2+1]));
% 
% %KF DOP computation
% KPDOP = sqrt(Cee_XYZ(1,1) + Cee_XYZ(2,2) + Cee_XYZ(3,3));
% KHDOP = sqrt(Cee_ENU(1,1) + Cee_ENU(2,2));
% KVDOP = sqrt(Cee_ENU(3,3));

kalman_initialized = 1;
