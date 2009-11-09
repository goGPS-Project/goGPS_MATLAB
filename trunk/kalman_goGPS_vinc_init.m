function kalman_goGPS_vinc_init (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ...
         ph1_Rsat, ph1_Msat, pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, phase, ref)

% SYNTAX:
%   kalman_goGPS_vinc_init (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ...
%   ph1_Rsat, ph1_Msat, pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, phase, ref);
%
% INPUT:
%   pos_M = Master given coordinates (X,Y,Z)
%   time = GPS time
%   Eph = satellites ephemerides
%   iono = ionosphere parameters (vector of zeros if not available)
%   pr1_Rsat = ROVER-SATELLITE code-pseudorange (carrier L1)
%   pr1_Msat = MASTER-SATELLITE code-pseudorange (carrier L1)
%   ph1_Rsat = ROVER-SATELLITE phase observations (carrier L1)
%   ph1_Msat = MASTER-SATELLITE phase observations (carrier L1)
%   pr2_Rsat = ROVER-SATELLITE code-pseudorange (carrier L2)
%   pr2_Msat = MASTER-SATELLITE code-pseudorange (carrier L2)
%   ph2_Rsat = ROVER-SATELLITE phase observations (carrier L2)
%   ph2_Msat = MASTER-SATELLITE phase observations (carrier L2)
%   phase = carrier L1 (phase=1), carrier L2 (phase=2)
%   ref = reference line
%
% DESCRIPTION:
%   Kalman filter initialization, POVER initial position estimate included
%   (X,Y,Z). Constrained path.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

global a f
global sigmaq0 sigmaq_vel sigmaq0_N
global cutoff o1 o2 o3 nN
global s0 ax ay az

global Xhat_t_t X_t1_t Yhat_t_t Y_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM

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

%vectors of zeros useful for matrix declaration
Z_nN_o1 = zeros(nN,o1);
Z_o1_nN = zeros(o1,nN);
Z_o1_o1 = zeros(o1);

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
Cvv(o1,o1) = sigmaq_vel;

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

if (length(sat_pr) >= 4)
    [pos_R, pos_SAT] = input_bancroft(pr1_Rsat(sat_pr), sat_pr, time(1), Eph);
else
    error('%d satellites are not enough to apply Bancroft algorithm\n', length(sat_pr));
end

pos_R = pos_R(1:3);
pos_SAT = pos_SAT(:,1:3);

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
[azR(sat_pr), elR(sat_pr), distR(sat_pr)] = topocent(pos_R, pos_SAT, a, f);

%azimuth, elevation, MASTER-SATELLITE distance estimate
[azM(sat_pr), elM(sat_pr), distM(sat_pr)] = topocent(pos_M, pos_SAT, a, f);

%elevation cut-off
sat_cutoff = find(elR > cutoff);
sat_pr = intersect(sat_pr,sat_cutoff);
sat = intersect(sat,sat_cutoff);

%previous pivot
pivot_old = 0;

%current pivot
if ~isempty(sat)
    [max_elR, i] = max(elR(sat));
    pivot = sat(i);
else
    [max_elR, i] = max(elR(sat_pr));
    pivot = sat_pr(i);
end
%pivot = find(elR == max(elR));

%--------------------------------------------------------------------------------------------
% SATELLITES CONFIGURATION
%--------------------------------------------------------------------------------------------

%satellites configuration: code only (-1), code and phase (+1);
conf_sat = zeros(32,1);
conf_sat(sat_pr) = -1;
conf_sat(sat) = +1;

%cycle-slips configuration (no cycle-slips)
conf_cs = zeros(32,1);

%number of satellites in view (not used)
nsat = size(sat_pr,1);

%--------------------------------------------------------------------------------------------
% KALMAN FILTER INITIAL STATE
%--------------------------------------------------------------------------------------------

%vectors of zeros useful for matrix declaration
Z_om_1 = zeros(o1-1,1);

%ROVER positioning by means of code double differences
if (phase(1) == 1)
    if (sum(abs(iono)) == 0) %if ionospheric parameters are not available; they are set equal to zero
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), time, sat_pr, pivot, Eph);
    else
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), time, sat_pr, pivot, Eph, iono);
    end
else
    if (sum(abs(iono)) == 0) %if ionospheric parameters are not available; they are set equal to zero
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), time, sat_pr, pivot, Eph);
    else
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), time, sat_pr, pivot, Eph, iono);
    end
end

%required re-estimate because only at the previous step a positioning
%precision in the order of centimeters can be obtained. 
if (phase(1) == 1)
    if (sum(abs(iono)) == 0) %if ionospheric parameters are not available; they are set equal to zero
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), time, sat_pr, pivot, Eph);
    else
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), time, sat_pr, pivot, Eph, iono);
    end
else
    if (sum(abs(iono)) == 0) %if ionospheric parameters are not available; they are set equal to zero
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), time, sat_pr, pivot, Eph);
    else
        [pos_R, cov_pos_R] = code_double_diff(pos_R, pr2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), time, sat_pr, pivot, Eph, iono);
    end
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

[dmin i] = min(d);

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

    [dmin i] = min(d);

    s_R = s0(i);
    pos_R = ref(i,:);

%     d1 = sqrt(sum((ref(i,:) - pos_R_proj(i,:)).^2));
%     d2 = sqrt(sum((ref(i+1,:) - pos_R_proj(i,:)).^2));
% 
%     if (d1 < d2)
%         s_R = s0(i);
%         pos_R = ref(i,:);
%     else
%         s_R = s0(i+1);
%         pos_R = ref(i+1,:);
%     end

end

%propagated error
sigmaq_s_R = (ax(i)^2*sigmaq_pos_R(1) + ay(i)^2*sigmaq_pos_R(2) + az(i)^2*sigmaq_pos_R(3)) ./ (ax(i)^2 + ay(i)^2 + az(i)^2)^2;

%satellites combinations initialization: initialized value
%if the satellite is in view, 0 if it is not
N1_stim = zeros(32,1);
N2_stim = zeros(32,1);
sigmaq_N1 = zeros(32,1);
sigmaq_N2 = zeros(32,1);

%phase double differences estimate, so to estimate the N values
if ~isempty(sat)
    [N1_stim(sat), sigmaq_N1(sat)] = amb_estimate_observ(pos_R, pos_M, pr1_Rsat(sat), pr1_Msat(sat), ph1_Rsat(sat), ph1_Msat(sat), Eph, time, pivot, sat, 1);
    [N2_stim(sat), sigmaq_N2(sat)] = amb_estimate_observ(pos_R, pos_M, pr2_Rsat(sat), pr2_Msat(sat), ph2_Rsat(sat), ph2_Msat(sat), Eph, time, pivot, sat, 2);
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

%initial point initialization, composed by 6(positions and velocities) + 
%32 o 64 (N combinations) variables
Xhat_t_t = [s_R; Z_om_1; N_stim];

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
%Cee(o1+1:o1+nN,o1+1:o1+nN) = sigmaq_N * eye(nN);
Cee(o1+1:o1+nN,o1+1:o1+nN) = sigmaq0_N * eye(nN);
