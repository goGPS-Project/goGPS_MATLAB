function [kalman_initialized] = kalman_goGPS_vinc_init (pos_R, pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ...
              ph1_Rsat, ph1_Msat, dop1_Rsat, dop1_Msat, pr2_Rsat, pr2_Msat, ...
              ph2_Rsat, ph2_Msat, dop2_Rsat, dop2_Msat, snr_R, snr_M, phase, ref, dtMdot)

% SYNTAX:
%   [kalman_initialized] = kalman_goGPS_vinc_init (pos_R, pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ...
%        ph1_Rsat, ph1_Msat, dop1_Rsat, dop1_Msat, pr2_Rsat, pr2_Msat, ...
%        ph2_Rsat, ph2_Msat, dop2_Rsat, dop2_Msat, snr_R, snr_M, phase, ref, dtMdot);
%
% INPUT:
%   pos_R = rover approximate coordinates (X,Y,Z)
%   pos_M = master known coordinates (X,Y,Z)
%   time = GPS time
%   Eph = satellites ephemerides
%   iono = ionosphere parameters (vector of zeros if not available)
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
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   phase = carrier L1 (phase=1), carrier L2 (phase=2)
%   ref = reference line
%   dtMdot = master receiver clock drift
%
% OUTPUT:
%   kalman_initialized = flag to point out whether Kalman has been successfully initialized
%
% DESCRIPTION:
%   Kalman filter initialization, POVER initial position estimate included
%   (X,Y,Z). Constrained path.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
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

global sigmaq0 sigmaq_vel sigmaq0_N
global cutoff o1 nN
global s0 ax ay az

global Xhat_t_t X_t1_t Yhat_t_t Y_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP
global doppler_pred_range1_R doppler_pred_range2_R
global doppler_pred_range1_M doppler_pred_range2_M

kalman_initialized = 0;

%--------------------------------------------------------------------------------------------
% SINGLE / DOUBLE FREQUENCY SELECTION
%--------------------------------------------------------------------------------------------

%nN represents the number of unknown phase ambiguities
if (length(phase) == 1)
    nN = 32;
else
    nN = 64;
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
% DYNAMIC MODEL OF KALMAN FILTER
%--------------------------------------------------------------------------------------------

%vectors of zeros useful for matrix construction
Z_nN_o1 = zeros(nN,o1);
Z_o1_nN = zeros(o1,nN);

%T matrix construction - system dynamics
%position and velocity equations
T0 = eye(o1) + diag(ones(o1-1,1),1);

%second order polynomial
% T0 = [1 1; 0 1];
%third order polynomial
% T0 = [1 1 0; 0 1 1; 0 0 1]

%initial comb_N matrix structure
N0 = eye(nN);

%system dynamics
%X(t+1)  = X(t) + Vx(t)
%Vx(t+1) = Vx(t)
%... <-- for the other two variables Y and Z
%comb_N(t+1) = comb_N(t)
T = [T0      Z_o1_nN;
     Z_nN_o1 N0];

%identity matrix construction of 38 variables (6 for positions and
%velocities + 32 or 64 for the number of satellites)
%for the following computation
I = eye(o1+nN);

%model error covariance matrix
Cvv = zeros(o1+nN);
Cvv(o1,o1) = sigmaq_vel; %#ok<NASGU>

%--------------------------------------------------------------------------------------------
% SATELLITES SELECTION
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

%only code and phase satellites
%sat_pr = sat;

%--------------------------------------------------------------------------------------------
% INITIAL POSITION ESTIMATE WITH BANCROFT ALGORITHM
%--------------------------------------------------------------------------------------------

if (sum(pos_R) == 0)
    if (length(sat_pr) >= 4)
        [pos_R] = input_bancroft(pr1_Rsat(sat_pr), sat_pr, time(1), Eph);
    else
        return
    end
end

%--------------------------------------------------------------------------------------------
% SATELLITE POSITION CORRECTION (AND MASTER DOPPLER SHIFT COMPUTATION)
%--------------------------------------------------------------------------------------------
pos_S = zeros(3,size(sat_pr));
pos_S_ttime = zeros(3,32);
vel_S = zeros(3,32);
ttime = zeros(32,1);
for i = 1:size(sat_pr)
    
    i_sat = sat_pr(i);
    
    %satellite position (with clock error and Earth rotation corrections)
    [pos_S(:,i), dt_S, pos_S_ttime(:,i_sat), vel_S(:,i_sat), ttime(i_sat,1)] = sat_corr(Eph, i_sat, time, pr1_Rsat(i_sat)); %#ok<ASGLU>

    if (nargin > 20 & ~isempty(dtMdot) & dop1_Msat(i_sat) == 0)
        [dop1_Msat(i_sat), dop2_Msat(i_sat)] = doppler_shift_approx(pos_M, zeros(3,1), pos_S_ttime(:,i_sat), vel_S(:,i_sat), ttime(i_sat,1), dtMdot, i_sat, Eph);
    end
end

%------------------------------------------------------------------------------------
% SATELLITES ELEVATION, PIVOT AND CUT-OFF CHECK
%-----------------------------------------------------------------------------------

%initialization
azR = zeros(32,1);
elR = zeros(32,1);
distR = zeros(32,1);
azM = zeros(32,1);
elM = zeros(32,1);
distM = zeros(32,1);

%azimuth, elevation, ROVER-SATELLITE distance estimate
[azR(sat_pr), elR(sat_pr), distR(sat_pr)] = topocent(pos_R, pos_S');

%azimuth, elevation, MASTER-SATELLITE distance estimate
[azM(sat_pr), elM(sat_pr), distM(sat_pr)] = topocent(pos_M, pos_S');

%elevation cut-off
sat_cutoff = find(elR > cutoff);
sat_pr = intersect(sat_pr,sat_cutoff);
sat = intersect(sat,sat_cutoff);

%previous pivot
pivot_old = 0;

%current pivot
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

%satellites configuration: code only (-1), code and phase (+1);
conf_sat = zeros(32,1);
conf_sat(sat_pr) = -1;
conf_sat(sat) = +1;

%cycle-slips configuration (no cycle-slips)
conf_cs = zeros(32,1);

%--------------------------------------------------------------------------------------------
% KALMAN FILTER INITIAL STATE
%--------------------------------------------------------------------------------------------

%vectors of zeros useful for matrix declaration
Z_om_1 = zeros(o1-1,1);
sigmaq_comb_N = zeros(32,1);

if (length(sat_pr) >= 4)
    
    %ROVER positioning by means of code double differences
    if (phase(1) == 1)
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr1_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono); %#ok<NASGU>
    else
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr2_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono); %#ok<NASGU>
    end
    
    %re-estimate to obtain better accuracy
    if (phase(1) == 1)
        [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_double_diff(pos_R, pr1_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr1_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono);
    else
        [pos_R, cov_pos_R, PDOP, HDOP, VDOP] = code_double_diff(pos_R, pr2_Rsat(sat_pr), snr_R(sat_pr), pos_M, pr2_Msat(sat_pr), snr_M(sat_pr), time, sat_pr, pivot, Eph, iono);
    end
    
    if isempty(cov_pos_R) %if a covariance matrix estimate was not possible (iso-determined problem)
        cov_pos_R = sigmaq0 * eye(3);
    end
    sigmaq_pos_R = diag(cov_pos_R);
    
    %projection over the constrain
    bx = pos_R(1) - ref(1:end-1,1) + ax.*s0(1:end-1);
    by = pos_R(2) - ref(1:end-1,2) + ay.*s0(1:end-1);
    bz = pos_R(3) - ref(1:end-1,3) + az.*s0(1:end-1);
    
    s_R = (ax.*bx + ay.*by + az.*bz) ./ (ax.^2 + ay.^2 + az.^2);
    
    pos_R_proj(:,1) = ref(1:end-1,1) + ax .* (s_R - s0(1:end-1));
    pos_R_proj(:,2) = ref(1:end-1,2) + ay .* (s_R - s0(1:end-1));
    pos_R_proj(:,3) = ref(1:end-1,3) + az .* (s_R - s0(1:end-1));
    
    %minimum distance estimate
    d = sqrt((pos_R(1) - pos_R_proj(:,1)).^2 + ...
        (pos_R(2) - pos_R_proj(:,2)).^2 + ...
        (pos_R(3) - pos_R_proj(:,3)).^2);
    
    [dmin i] = min(d); %#ok<ASGLU>
    
    %cartesian coordinates positioning
    if ((pos_R_proj(i,1) >= min(ref(i,1),ref(i+1,1))) & (pos_R_proj(i,1) <= max(ref(i,1),ref(i+1,1))) & ...
            (pos_R_proj(i,2) >= min(ref(i,2),ref(i+1,2))) & (pos_R_proj(i,2) <= max(ref(i,2),ref(i+1,2))) & ...
            (pos_R_proj(i,3) >= min(ref(i,3),ref(i+1,3))) & (pos_R_proj(i,3) <= max(ref(i,3),ref(i+1,3))))
        
        s_R = s_R(i);
        pos_R = pos_R_proj(i,:);
        
    else
        
        d = sqrt((pos_R(1) - ref(:,1)).^2 + ...
            (pos_R(2) - ref(:,2)).^2 + ...
            (pos_R(3) - ref(:,3)).^2);
        
        [dmin i] = min(d); %#ok<ASGLU>
        
        s_R = s0(i);
        pos_R = ref(i,:);
        
    end
    
    %propagated error
    sigmaq_s_R = (ax(i)^2*sigmaq_pos_R(1) + ay(i)^2*sigmaq_pos_R(2) + az(i)^2*sigmaq_pos_R(3)) ./ (ax(i)^2 + ay(i)^2 + az(i)^2)^2;
    
    %do not use least squares ambiguity estimation
    % NOTE: LS amb. estimation is automatically switched off if the number of
    % satellites with phase available is not sufficient
    if (size(sat_pr,1) + size(sat,1) - 2 <= 3 + size(sat,1) - 1)
        
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
            [null_pos_R, null_cov_pos_R, comb_N1_stim(sat), cov_comb_N1_stim] = code_phase_double_diff(pos_R', pr1_Rsat(sat), ph1_Rsat(sat), snr_R(sat), pos_M, pr1_Msat(sat), ph1_Msat(sat), snr_M(sat), time, sat, pivot, Eph, 1, iono); %#ok<ASGLU>
            [null_pos_R, null_cov_pos_R, comb_N2_stim(sat), cov_comb_N2_stim] = code_phase_double_diff(pos_R', pr2_Rsat(sat), ph2_Rsat(sat), snr_R(sat), pos_M, pr2_Msat(sat), ph2_Msat(sat), snr_M(sat), time, sat, pivot, Eph, 2, iono); %#ok<ASGLU>
        end
        
        if isempty(cov_comb_N1_stim) %if it was not possible to compute the covariance matrix
            cov_comb_N1_stim = sigmaq0_N * eye(length(sat));
        end
        
        if isempty(cov_comb_N2_stim) %if it was not possible to compute the covariance matrix
            cov_comb_N2_stim = sigmaq0_N * eye(length(sat));
        end
        
        if (length(phase) == 2)
            comb_N_stim = [comb_N1_stim; comb_N2_stim];
            sigmaq_comb_N(sat) = diag(cov_comb_N1_stim);
            sigmaq_comb_N(sat+32) = diag(cov_comb_N2_stim);
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
else
    return
end

%initial point initialization, composed by 6(positions and velocities) +
%32 o 64 (N combinations) variables
Xhat_t_t = [s_R; Z_om_1; comb_N_stim];

%point estimate at the step t+1 X Vx Y Vy Z Vz comb_N
%estimate at the step t, because the initial velocities is equal to 0
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
Cee(1,1) = sigmaq_s_R;
Cee(2:o1,2:o1) = sigmaq0 * eye(o1-1);
Cee(o1+1:o1+nN,o1+1:o1+nN) = diag(sigmaq_comb_N);

%--------------------------------------------------------------------------------------------
% DOPPLER-BASED PREDICTION OF PHASE RANGES
%--------------------------------------------------------------------------------------------
if (dop1_Rsat(sat))
    doppler_pred_range1_R(sat,1) = ph1_Rsat(sat) - dop1_Rsat(sat);
end
if (dop2_Rsat(sat))
    doppler_pred_range2_R(sat,1) = ph2_Rsat(sat) - dop2_Rsat(sat);
end
if (dop1_Msat(sat))
    doppler_pred_range1_M(sat,1) = ph1_Msat(sat) - dop1_Msat(sat);
end
if (dop2_Msat(sat))
    doppler_pred_range2_M(sat,1) = ph2_Msat(sat) - dop2_Msat(sat);
end

kalman_initialized = 1;
