function goGPS_LS_DD_code_phase(time_rx, XM, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase, flag_IAR)

% SYNTAX:
%   goGPS_LS_DD_code_phase(time_rx, XM, pr1_R, pr1_M, pr2_R, pr2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase, flag_IAR);
%
% INPUT:
%   time_rx = GPS reception time
%   XM    = MASTER position
%   pr1_R = ROVER code observations (L1 carrier)
%   pr1_M = MASTER code observations (L1 carrier)
%   pr2_R = ROVER code observations (L2 carrier)
%   pr2_M = MASTER code observations (L2 carrier)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   Eph   = satellite ephemeris
%   SP3   = structure containing precise ephemeris and clock
%   iono  = ionosphere parameters
%   lambda = wavelength matrix (depending on the enabled constellations)
%   phase  = L1 carrier (phase=1), L2 carrier (phase=2)
%   flag_IAR = boolean variable to enable/disable integer ambiguity resolution
%
% DESCRIPTION:
%   Computation of the receiver position (X,Y,Z).
%   Relative (double difference) positioning by least squares adjustment
%   on code and phase observations.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Hendy F. Suhandri
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
global cutoff snr_threshold cond_num_threshold o1 o2 o3

global Xhat_t_t Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP
global ratiotest mutest succ_rate fixed_solution
global n_sys

%number of unknown phase ambiguities
if (length(phase) == 1)
    nN = 32;
else
    nN = 64;
end

%compute inter-frequency factors (for the ionospheric delay)
ionoFactor = goGNSS.getInterFreqIonoFactor(lambda);

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1_R,1);

%covariance matrix initialization
cov_XR = [];
cov_N1 = [];
cov_N2 = [];

%topocentric coordinate initialization
azR   = zeros(nSatTot,1);
elR   = zeros(nSatTot,1);
distR = zeros(nSatTot,1);
azM   = zeros(nSatTot,1);
elM   = zeros(nSatTot,1);
distM = zeros(nSatTot,1);

%--------------------------------------------------------------------------------------------
% SATELLITE SELECTION
%--------------------------------------------------------------------------------------------

if (length(phase) == 2)
    sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) & ...
                   (pr2_R ~= 0) & (pr2_M ~= 0) );
    sat    = find( (ph1_R ~= 0) & (ph1_M ~= 0) & ...
                   (ph2_R ~= 0) & (ph2_M ~= 0) );
else
    if (phase == 1)
        sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) );
        sat    = find( (ph1_R ~= 0) & (ph1_M ~= 0) );
    else
        sat_pr = find( (pr2_R ~= 0) & (pr2_M ~= 0) );
        sat    = find( (ph2_R ~= 0) & (ph2_M ~= 0) );
    end
end
sat_pr = sat_pr(ismember(sat_pr, Eph(30,:)));
sat = sat(ismember(sat, Eph(30,:)));

N1 = zeros(nSatTot,1);
N2 = zeros(nSatTot,1);
Z_om_1 = zeros(o1-1,1);
sigma2_N = zeros(nN,1);

min_nsat_LS = 3 + n_sys;

if (size(sat_pr,1) >= min_nsat_LS)
    
    sat_pr_old = sat_pr;
    
    if (phase == 1)
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr1_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, [], XM, [],  [], sat_pr,   [], lambda(sat_pr,:),   cutoff, snr_threshold, phase, 2, 0); %#ok<NASGU,ASGLU>
        if (length(sat_pr_M) < min_nsat_LS); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, [], [], XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, phase, 0, 1); %#ok<ASGLU>
    else
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr2_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, [], XM, [],  [], sat_pr,   [], lambda(sat_pr,:),   cutoff, snr_threshold, phase, 2, 0); %#ok<NASGU,ASGLU>
        if (length(sat_pr_M) < min_nsat_LS); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, [], [], XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, phase, 0, 1); %#ok<ASGLU>
    end
    
    %keep only satellites that rover and master have in common
    [sat_pr, iR, iM] = intersect(sat_pr_R, sat_pr_M);
    XS = XS(iR,:);
    if (~isempty(err_tropo_R))
        err_tropo_R = err_tropo_R(iR);
        err_iono_R  = err_iono_R (iR);
        err_tropo_M = err_tropo_M(iM);
        err_iono_M  = err_iono_M (iM);
    end
    
    %apply cutoffs also to phase satellites
    sat_removed = setdiff(sat_pr_old, sat_pr);
    sat(ismember(sat,sat_removed)) = [];
    
    %--------------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
    %--------------------------------------------------------------------------------------------
    
    %satellite configuration
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat_pr,1) = -1;
    conf_sat(sat,1)    = +1;
    
    %no cycle-slips when working with code only
    conf_cs = zeros(nSatTot,1);
    
    %previous pivot
    pivot_old = 0;
    
    %actual pivot
    [null_max_elR, pivot_index] = max(elR(sat)); %#ok<ASGLU>
    pivot = sat(pivot_index);
    
    %--------------------------------------------------------------------------------------------
    % LEAST SQUARES SOLUTION
    %--------------------------------------------------------------------------------------------
    
    %if at least min_nsat_LS satellites are available after the cutoffs, and if the
    % condition number in the least squares does not exceed the threshold
    if (size(sat,1) >= min_nsat_LS && cond_num < cond_num_threshold)
        
        %loop is needed to improve the atmospheric error correction
        for i = 1 : 3

            if (phase == 1)
                [XR, N1(sat), cov_XR, cov_N1, PDOP, HDOP, VDOP] = LS_DD_code_phase(XR, XM, XS, pr1_R(sat), ph1_R(sat), snr_R(sat), pr1_M(sat), ph1_M(sat), snr_M(sat), elR(sat), elM(sat), err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda(sat,1), flag_IAR);
            else
                [XR, N2(sat), cov_XR, cov_N2, PDOP, HDOP, VDOP] = LS_DD_code_phase(XR, XM, XS, pr2_R(sat), ph2_R(sat), snr_R(sat), pr2_M(sat), ph2_M(sat), snr_M(sat), elR(sat), elM(sat), err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda(sat,2), flag_IAR);
            end
            
            if (i < 3)
                if (~isempty(ratiotest))
                    ratiotest(end) = [];
                    mutest(end) = [];
                end
                if (~isempty(fixed_solution))
                    fixed_solution(end) = [];
                    succ_rate(end) = [];
                end
            end
            
            [phiR, lamR, hR] = cart2geod(XR(1), XR(2), XR(3));
            [azR(azR ~= 0), elR(elR ~= 0), distR(distR ~= 0)] = topocent(XR, XS);
            
            err_tropo_R = tropo_error_correction(elR(elR ~= 0), hR);
            err_iono_R = iono_error_correction(phiR*180/pi, lamR*180/pi, azR(azR ~= 0), elR(elR ~= 0), time_rx, iono, []);
            
            %correct the ionospheric errors for different frequencies
            err_iono_R = ionoFactor(sat,phase).*err_iono_R;
        end
        
        if isempty(cov_N1) %if it was not possible to compute the covariance matrix
            cov_N1 = sigmaq0_N * eye(length(sat));
        end
        
        if isempty(cov_N2) %if it was not possible to compute the covariance matrix
            cov_N2 = sigmaq0_N * eye(length(sat));
        end
        
        if (length(phase) == 2)
            N = [N1; N2];
            sigma2_N(sat) = diag(cov_N1);
            %sigma2_N(sat) = (sigmaq_cod1 / lambda1^2) * ones(length(sat),1);
            sigma2_N(sat+nN) = diag(cov_N2);
            %sigma2_N(sat+nN) = (sigmaq_cod2 / lambda2^2) * ones(length(sat),1);
        else
            if (phase == 1)
                N = N1;
                sigma2_N(sat) = diag(cov_N1);
                %sigma2_N(sat) = (sigmaq_cod1 / lambda1^2) * ones(length(sat),1);
            else
                N = N2;
                sigma2_N(sat) = diag(cov_N2);
                %sigma2_N(sat) = (sigmaq_cod2 / lambda2^2) * ones(length(sat),1);
            end
        end
    else
        if (~isempty(Xhat_t_t))
            XR = Xhat_t_t([1,o1+1,o2+1]);
            N  = Xhat_t_t(o3+1:end);
            pivot = 0;

            fixed_solution = [fixed_solution 0];
            succ_rate = [succ_rate NaN];
        else
            return
        end
    end
       
else
    if (~isempty(Xhat_t_t))
        XR = Xhat_t_t([1,o1+1,o2+1]);
        N  = Xhat_t_t(o3+1:end);
        pivot = 0;

        fixed_solution = [fixed_solution 0];
        succ_rate = [succ_rate NaN];
    else
        return
    end
end

if isempty(cov_XR) %if it was not possible to compute the covariance matrix
    cov_XR = sigmaq0 * eye(3);
end
sigma2_XR = diag(cov_XR);

%initialization of the initial point with 6(positions and velocities) +
%nSatTot or nSatTotx2 (N combinations) variables
Xhat_t_t = [XR(1); Z_om_1; XR(2); Z_om_1; XR(3); Z_om_1; N];

%--------------------------------------------------------------------------------------------
% INITIAL STATE COVARIANCE MATRIX
%--------------------------------------------------------------------------------------------

%initial state covariance matrix
Cee(:,:) = zeros(o3+nN);
Cee(1,1) = sigma2_XR(1);
Cee(o1+1,o1+1) = sigma2_XR(2);
Cee(o2+1,o2+1) = sigma2_XR(3);
Cee(2:o1,2:o1) = sigmaq0 * eye(o1-1);
Cee(o1+2:o2,o1+2:o2) = sigmaq0 * eye(o1-1);
Cee(o2+2:o3,o2+2:o3) = sigmaq0 * eye(o1-1);
Cee(o3+1:o3+nN,o3+1:o3+nN) = diag(sigma2_N);
