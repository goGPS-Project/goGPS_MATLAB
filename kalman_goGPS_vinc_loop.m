function [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop ...
         (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ph1_Rsat, ph1_Msat,...
         pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, snr_R, snr_M, phase, ref)

% SYNTAX:
%   [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop ...
%   (pos_M, time, Eph, iono, pr1_Rsat, pr1_Msat, ph1_Rsat, ph1_Msat,...
%   pr2_Rsat, pr2_Msat, ph2_Rsat, ph2_Msat, snr_R, snr_M, phase, ref);
%
% INPUT:
%   pos_M = Master given cooridnates (X,Y,Z)
%   time = GPS time
%   Eph = satellites ephemerides
%   iono = ionosphere parameters  (not used)
%   pr1_Rsat = ROVER-SATELLITE code-pseudorange (carrier L1)
%   pr1_Msat = MASTER-SATELLITE code-pseudorange (carrier L1)
%   ph1_Rsat = ROVER-SATELLITE phase observations (carrier L1)
%   ph1_Msat = MASTER-SATELLITE phase observations (carrier L1)
%   pr2_Rsat = ROVER-SATELLITE code-pseudorange (carrier L2)
%   pr2_Msat = MASTER-SATELLITE code-pseudorange (carrier L2)
%   ph2_Rsat = ROVER-SATELLITE phase observations (carrier L2)
%   ph2_Msat = MASTER-SATELLITE phase observations (carrier L2)
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

global lambda1 lambda2 f1 f2 a f

global sigmaq_vel sigmaq0_N
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph weights
global min_nsat cutoff o1 o2 o3 nN
global s0 ax ay az

global Xhat_t_t X_t1_t Yhat_t_t Y_t1_t T I Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM

%----------------------------------------------------------------------------------------
% INITIALIZATION
%----------------------------------------------------------------------------------------

%output variables for events warning (birth, death, etc)
check_on = 0;
check_off = 0;
check_pivot = 0;
check_cs = 0;

%azimuth, elevation, ROVER-satellite and MASTER-satellite distance
azR = zeros(32,1);
azM = zeros(32,1);
elR = zeros(32,1);
elM = zeros(32,1);
distR = zeros(32,1);
distM = zeros(32,1);

%----------------------------------------------------------------------------------------
% APPROXIMATE COORDINATES FOR OBSERVATION EQUATIONS LINEARIZATION
%----------------------------------------------------------------------------------------

%approximate coordinates X Y Z
X_app = Y_t1_t(1);
Y_app = Y_t1_t(2);
Z_app = Y_t1_t(3);

%----------------------------------------------------------------------------------------
% COVARIANCE MATRIX FOR THE MODEL ERROR
%----------------------------------------------------------------------------------------

%re-initialization of the model error Cvv matrix, since it is frequently
%modified when cycle slips or new satellite births occurr

Cvv = zeros(o1+nN);
Cvv(o1,o1) = sigmaq_vel;

%------------------------------------------------------------------------------------
% SATELLITES SELECTION
%------------------------------------------------------------------------------------

%extraction of the satellites in view
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

%------------------------------------------------------------------------------------
% SATELLITES ELEVATION, PIVOT AND CUT-OFF CHECK
%------------------------------------------------------------------------------------

j = 1;
bad_sat = [];

for i = 1:size(sat_pr)

    %satellite position correction (clock and rotation)
    Rot_X = sat_corr(Eph, sat_pr(i), time, pr1_Rsat(sat_pr(i)), Y_t1_t');

    %azimuth, elevation, ROVER-SATELLITE distance estimate
    [azR(sat_pr(i)), elR(sat_pr(i)), distR(sat_pr(i))] = topocent(Y_t1_t, Rot_X', a, f);

    %azimuth, elevation, MASTER-SATELLITE distance estimate
    [azM(sat_pr(i)), elM(sat_pr(i)), distM(sat_pr(i))] = topocent(pos_M, Rot_X', a, f);

    %elevation test
    if (elR(sat_pr(i)) < cutoff)
        bad_sat(j,1) = sat_pr(i);
        j = j + 1;
    end
end

%elevation cut-off
sat_pr(find(ismember(sat_pr,bad_sat) == 1)) = [];
sat(find(ismember(sat,bad_sat) == 1)) = [];

%previous pivot
if (pivot ~= 0)
    pivot_old = pivot;
end

%current pivot
if ~isempty(sat)
    [max_elR, i] = max(elR(sat));
    pivot = sat(i);
else
    [max_elR, i] = max(elR(sat_pr));
    pivot = sat_pr(i);
end
%pivot = find(elR == max(elR));

%----------------------------------------------------------------------------------------
% SATELLITES CONFIGURATION
%----------------------------------------------------------------------------------------

%previous satellites configuration (phase measurement only)
sat_old = find(conf_sat == 1);

%satellites configuration: code only (-1), code and phase (+1);
conf_sat = zeros(32,1);
conf_sat(sat_pr) = -1;
conf_sat(sat) = +1;

%cycle-slips configuration
conf_cs = zeros(32,1);

%total number of satellites in view
nsat = size(sat_pr,1);
n = nsat - 1;

%------------------------------------------------------------------------------------
% OBSERVATION EQUATIONS
%------------------------------------------------------------------------------------

%when the number of satellites in view is bigger or equal to min_sat
if (nsat >= min_nsat)

    %rows in which the phase observation is available
    p = find(ismember(setdiff(sat_pr,pivot),setdiff(sat,pivot))==1);

    %function that allows to estimate the Kalman filter parameters
    [alfa1, prstim1, ddc1, ddp1] = input_kalman_vinc (Y_t1_t', pr1_Rsat(sat_pr), ph1_Rsat(sat_pr), pos_M, pr1_Msat(sat_pr), ph1_Msat(sat_pr), time, sat_pr, pivot, Eph, 1);
    [alfa2, prstim2, ddc2, ddp2] = input_kalman_vinc (Y_t1_t', pr2_Rsat(sat_pr), ph2_Rsat(sat_pr), pos_M, pr2_Msat(sat_pr), ph2_Msat(sat_pr), time, sat_pr, pivot, Eph, 2);

    %vectors of zeros useful for the matrix declaration
    Z_1_nN = zeros(1,nN);
    Z_n_nN = zeros(n,nN);
    Z_n_om = zeros(n,o1-1);
    Z_1_om = zeros(1,o1-1);

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
    y0_cod1 = ddc1 - prstim1 + alfa1*X_t1_t(1);
    y0_cod2 = ddc2 - prstim2 + alfa2*X_t1_t(1);

    %Y0 vector estimate for the phase
    if ~isempty(p)
        y0_fas1 = ddp1(p) - prstim1(p) + alfa1(p)*X_t1_t(1);
        y0_fas2 = ddp2(p) - prstim2(p) + alfa2(p)*X_t1_t(1);
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
    % COVARIANCE MATRIX OF THE OBSERVATIONS
    %------------------------------------------------------------------------------------
 
    %vectors of zeros useful for the matrix declaration
    Z_n_n = zeros(n,n);
    
    %weight function parameters
    snr_a = 30;
    snr_0 = 10;
    snr_1 = 50;
    snr_A = 30;

    if (weights == 0)

        %covariance matrix of the noise measurement
        %matrix code-code or phase-phase construction base Q
        Q = 2*ones(n) + 2*eye(n);

    else
        if (weights == 1)

            %weight vectors (elevation)
            q_R = 1 ./ (sin(elR * pi/180).^2);
            q_M = 1 ./ (sin(elM * pi/180).^2);

        elseif (weights == 2)

		    %weight vectors (signal-to-noise ratio)
            q_R = 10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1);
            q_R(find(snr_R >= snr_1)) = 1;
            q_M = 10.^(-(snr_M-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_M-snr_1)+1);
            q_M(find(snr_M >= snr_1)) = 1;

        elseif (weights == 3)

            %weight vectors (elevation and signal-to-noise ratio)
            q_R = 1 ./ (sin(elR * pi/180).^2) .* (10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1));
            q_R(find(snr_R >= snr_1)) = 1;
            q_M = 1 ./ (sin(elM * pi/180).^2) .* (10.^(-(snr_M-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_M-snr_1)+1));
            q_M(find(snr_M >= snr_1)) = 1;

        end

        q_RP = q_R(pivot,1);                  % ROVER-PIVOT
        q_MP = q_M(pivot,1);                  % MASTER-PIVOT
        q_RS = q_R(sat_pr);                   % ROVER-generic satellite
        q_MS = q_M(sat_pr);                   % MASTER-generic satellite
        q_RS(find(sat_pr==pivot)) = [];       % ROVER-generic satellite (without pivot)
        q_MS(find(sat_pr==pivot)) = [];       % MASTER-generic satellite (without pivot)

        %covariance matrix of the noise measurement
        %matrix code-code or phase-phase construction base Q
        Q = (q_RP + q_MP) * ones(n) + diag(q_RS + q_MS);
    end

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
    % SATELLITES BIRTH/DEATH
    %------------------------------------------------------------------------------------

    %dying satellite quest
    if (length(sat) < length(sat_old))

        check_off = 1;

        % dead satellites saving
        sat_dead = setdiff(sat_old,sat);

        %for dying satellites it is fundamental to set its N combination 
        %with the PIVOT to 0. Furthermore it is convenient to raise
        %its uncertainty (not necessary because it is raised anyway
        %in caso of its future birth)
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

        %dying satellites display
        ['Satelliti morti al tempo ' num2str(time) ' ' num2str(sat_dead')];
    end

    %borning satellite quest
    if (length(sat) > length(sat_old))

        check_on = 1;

        %born satellites saving
        sat_born = setdiff(sat,sat_old);

        %for borning satellites it is fundamental to estimate its N 
        %combination with the PIVOT. Furthermore it is convenient to raise
        %its uncertainty because it is necessary to wait for about
        %10 epochs in order to have an estimate correct enough.
        sigmaq_s_R = diag(T*Cee*T');
        sigmaq_s_R = sigmaq_s_R(1);

        %curvilinear coordinate localization
        i = find((X_t1_t(1) >= s0(1:end-1)) & (X_t1_t(1) <= s0(2:end)));

        sigmaq_pos_R(1,1) = ax(i)^2 * sigmaq_s_R(1);
        sigmaq_pos_R(2,1) = ay(i)^2 * sigmaq_s_R(1);
        sigmaq_pos_R(3,1) = az(i)^2 * sigmaq_s_R(1);

        %N combination estimate according to the equality between
        %crossed phases and pseudoranges
        [comb_N1, sigmaq_N1] = amb_estimate_approx(Y_t1_t', pos_M, sigmaq_pos_R, pr1_Rsat(sat), pr1_Msat(sat), ph1_Rsat(sat), ph1_Msat(sat), Eph, time, pivot, sat, 1);
        [comb_N2, sigmaq_N2] = amb_estimate_approx(Y_t1_t', pos_M, sigmaq_pos_R, pr2_Rsat(sat), pr2_Msat(sat), ph2_Rsat(sat), ph2_Msat(sat), Eph, time, pivot, sat, 2);

        index = find(ismember(sat,sat_born) == 0);
        comb_N1(index) = [];
        comb_N2(index) = [];

        %estimated parameters saving
        if (length(phase) == 2)
            X_t1_t(o1+sat_born,1) = comb_N1;
            X_t1_t(o1+32+sat_born,1) = comb_N2;
            %Cvv(o1+sat_born,o1+sat_born) = sigmaq_N1 * eye(size(sat_born,1));
            %Cvv(o1+32+sat_born,o1+32+sat_born) = sigmaq_N2 * eye(size(sat_born,1));
            Cvv(o1+sat_born,o1+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            Cvv(o1+32+sat_born,o1+32+sat_born) = sigmaq0_N * eye(size(sat_born,1));
        else
            if (phase == 1)
                X_t1_t(o1+sat_born,1) = comb_N1;
                %Cvv(o1+sat_born,o1+sat_born) = sigmaq_N1 * eye(size(sat_born,1));
                Cvv(o1+sat_born,o1+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            else
                X_t1_t(o1+sat_born,1) = comb_N2;
                %Cvv(o1+sat_born,o1+sat_born) = sigmaq_N2 * eye(size(sat_born,1));
                Cvv(o1+sat_born,o1+sat_born) = sigmaq0_N * eye(size(sat_born,1));
            end
        end

        %borning satellites display
        ['Satelliti nati al tempo ' num2str(time) ' ' num2str(sat_born')]; 
    end

    %------------------------------------------------------------------------------------
    % PIVOT CHANGE
    %------------------------------------------------------------------------------------

    %PIVOT change quest
    if (pivot ~= pivot_old)

        check_pivot = 1;

        %matrix construction for updating the PIVOT change
        %sat is the vector with the current satellites in view
        %nsat is the sat vector size
        R = zeros(32);
        R(sat,sat) = eye(length(sat));
        R(sat,pivot) = -1;
        R(pivot_old,pivot_old) = 0;
        R(pivot,pivot) = 0;

        I0 = eye(o1);
        Z_32_o1 = zeros(32,o1);
        Z_o1_32 = zeros(o1,32);
        Z_32_32 = zeros(32,32);

        %full matrix construction
        %sat_old, sat
        if (length(phase) == 2)
            A = [I0 Z_o1_32 Z_o1_32; Z_32_o1 R Z_32_32; Z_32_o1 Z_32_32 R];
        else
            A = [I0 Z_o1_32; Z_32_o1 R];
        end

        %new estimated state matrix construction
        X_t1_t = A*X_t1_t;

        %re-estimate of the covariance matrix Cee at the previous epoch
        Cee = A*Cee*A';

        %pivot changes display
        ['Pivot change at time ' num2str(time) ' from ' num2str(pivot_old) ' to ' num2str(pivot)];
    end

    %------------------------------------------------------------------------------------
    % CYCLE-SLIPS
    %------------------------------------------------------------------------------------

    if ~isempty(sat)

        %verify cycle slip occurring or epoch under analysis lack.
        %The system state changes only for phase ambiguities.

        if (length(phase) == 2)
            [cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip_kalman(pos_M, Y_t1_t', X_t1_t(o1+1:o1+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
            [cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip_kalman(pos_M, Y_t1_t', X_t1_t(o1+33:o1+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
            %[cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip_observv(X_t1_t(o1+1:o1+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
            %[cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip_observv(X_t1_t(o1+33:o1+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
            %[cycle_slip_found1, N_slip1, sat_slip1] = cycle_slip(pos_M, Y_t1_t', X_t1_t(o1+1:o1+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, pivot, sat, time-1, time, 1);
            %[cycle_slip_found2, N_slip2, sat_slip2] = cycle_slip(pos_M, Y_t1_t', X_t1_t(o1+33:o1+64), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, pivot, sat, time-1, time, 2);

            if (cycle_slip_found1 == 1)
                check_cs = 1;
                conf_cs(sat_slip1) = 1;
                X_t1_t(o1+sat_slip1) = N_slip1;
                Cvv(o1+sat_slip1,o1+sat_slip1) = sigmaq0_N * eye(size(sat_slip1,1));
            end
            if (cycle_slip_found2 == 1)
                check_cs = 1;
                conf_cs(sat_slip2) = 1;
                X_t1_t(o1+32+sat_slip2) = N_slip2;
                Cvv(o1+32+sat_slip2,o1+32+sat_slip2) = sigmaq0_N * eye(size(sat_slip2,1));
            end
        else
            if (phase == 1)
                [cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman(pos_M, Y_t1_t', X_t1_t(o1+1:o1+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip_observv(X_t1_t(o1+1:o1+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, time, pivot, sat, 3, 1);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip(pos_M, Y_t1_t', X_t1_t(o1+1:o1+32), ph1_Rsat(sat), ph1_Msat(sat), pr1_Rsat(sat), pr1_Msat(sat), Eph, pivot, sat, time-1, time, 1);
            else
                [cycle_slip_found, N_slip, sat_slip] = cycle_slip_kalman(pos_M, Y_t1_t', X_t1_t(o1+1:o1+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip_observv(X_t1_t(o1+1:o1+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, time, pivot, sat, 3, 2);
                %[cycle_slip_found, N_slip, sat_slip] = cycle_slip(pos_M, Y_t1_t', X_t1_t(o1+1:o1+32), ph2_Rsat(sat), ph2_Msat(sat), pr2_Rsat(sat), pr2_Msat(sat), Eph, pivot, sat, time-1, time, 2);
            end
            if (cycle_slip_found == 1)
                check_cs = 1;
                conf_cs(sat_slip) = 1;
                X_t1_t(o1+sat_slip) = N_slip;
                Cvv(o1+sat_slip,o1+sat_slip) = sigmaq0_N * eye(size(sat_slip,1));
            end
        end
    end
else
    %to show that for any satellites configuration
    %the data have not been processed (motion with dynamic only)
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