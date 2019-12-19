function goGPS_BLK_DD_code_phase_static(time_rx, XR0, XM, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M, Eph, SP3, iono, lambda, frequencies, p_rate, antenna_PCV)

% SYNTAX:
%   goGPS_LS_DD_code_phase(time_rx, XR0, XM, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase, p_rate, antenna_PCV);
%
% INPUT:
%   time_rx = GPS reception time
%   XR0   = ROVER approximate position
%   XM    = MASTER position
%   pr1_R = ROVER code observations (L1 carrier)
%   pr1_M = MASTER code observations (L1 carrier)
%   pr2_R = ROVER code observations (L2 carrier)
%   pr2_M = MASTER code observations (L2 carrier)
%   ph1_R = ROVER phase observations (L1 carrier)
%   ph1_M = MASTER phase observations (L1 carrier)
%   ph2_R = ROVER phase observations (L2 carrier)
%   ph2_M = MASTER phase observations (L2 carrier)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   Eph   = satellite ephemeris
%   SP3   = structure containing precise ephemeris and clock
%   iono  = ionosphere parameters
%   lambda = wavelength matrix (depending on the enabled constellations)
%   phase  = L1 carrier (phase=1), L2 carrier (phase=2)
%   p_rate = processing interval [s]
%   antenna_PCV = antenna phase center variation
%
% DESCRIPTION:
%   Computation of the receiver position (X,Y,Z).
%   Relative (double difference) positioning by least squares adjustment
%   on code and phase observations.


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
%  Contributors:     Hendy F. Suhandri, ...
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

global cutoff snr_threshold cond_num_threshold

global conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP
global n_sys

global y0_epo A_epo b_epo Q_epo
y0_epo = [];
A_epo  = [];
b_epo  = [];
Q_epo  = [];

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1_R,1);

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

% Find sat in common, between master and rover
if (length(frequencies) == 2)
    sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) & ...
                   (pr2_R ~= 0) & (pr2_M ~= 0) );
    sat    = find( (ph1_R ~= 0) & (ph1_M ~= 0) & ...
                   (ph2_R ~= 0) & (ph2_M ~= 0) );
else
    if (frequencies == 1)
        sat_pr = find( (pr1_R ~= 0) & (pr1_M ~= 0) );
        sat    = find( (ph1_R ~= 0) & (ph1_M ~= 0) );
    else
        sat_pr = find( (pr2_R ~= 0) & (pr2_M ~= 0) );
        sat    = find( (ph2_R ~= 0) & (ph2_M ~= 0) );
    end
end
if (isempty(SP3))
    eph_avail = Eph(30,:);
else
    eph_avail = SP3.avail;
end
sat_pr = sat_pr(ismember(sat_pr, eph_avail));
sat = sat(ismember(sat, eph_avail));

min_nsat_LS = 3 + n_sys;

flag_XR = 2;

if (size(sat_pr,1) >= min_nsat_LS)

    sat_pr_old = sat_pr;

    if (frequencies == 1)
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr1_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, [],  XM, [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, frequencies, p_rate,       2, 0); %#ok<ASGLU>
        if (length(sat_pr_M) < min_nsat_LS); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, [], XR0, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 1); %#ok<ASGLU>
    else
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_pr_M, elM(sat_pr_M), azM(sat_pr_M), distM(sat_pr_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr2_M(sat_pr),   snr_M(sat_pr),   Eph, SP3, iono, [],  XM, [],  [], sat_pr,    [], lambda(sat_pr,:),   cutoff, snr_threshold, frequencies, p_rate,       2, 0); %#ok<ASGLU>
        if (length(sat_pr_M) < min_nsat_LS); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_pr_R, elR(sat_pr_R), azR(sat_pr_R), distR(sat_pr_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2_R(sat_pr_M), snr_R(sat_pr_M), Eph, SP3, iono, [], XR0, XS, dtS, sat_pr_M, sys, lambda(sat_pr_M,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 1); %#ok<ASGLU>
    end

    %keep only satellites that rover and master have in common
    [sat_pr, iR, iM] = intersect(sat_pr_R, sat_pr_M);
    XS = XS(iR,:);
    sys = sys(iR);
    if (~isempty(err_tropo_R))
        err_tropo_R = err_tropo_R(iR);
        err_iono_R  = err_iono_R (iR);
        err_tropo_M = err_tropo_M(iM);
        err_iono_M  = err_iono_M (iM);
    end

    %apply cutoffs also to phase satellites
    sat_removed = setdiff(sat_pr_old, sat_pr);
    sat(ismember(sat,sat_removed)) = [];

    % keep only satellites that rover and master have in common both in phase and code
    [sat_pr, i_pr] = intersect(sat_pr, sat);
    XS = XS(i_pr,:);
    sys = sys(i_pr);
    if (~isempty(err_tropo_R))
        err_tropo_R = err_tropo_R(i_pr);
        err_iono_R  = err_iono_R (i_pr);
        err_tropo_M = err_tropo_M(i_pr);
        err_iono_M  = err_iono_M (i_pr);
    end
    
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
    % PHASE CENTER VARIATIONS
    %--------------------------------------------------------------------------------------------

    %compute PCV: phase and code 1
    [~, index_ph]=intersect(sat_pr,sat);

    if (~isempty(antenna_PCV) && antenna_PCV(2).n_frequency ~= 0) % master
        index_master = 2;
        PCO1_M = PCO_correction(antenna_PCV(index_master), XR0, XS, sys, 1);
        PCV1_M = PCV_correction(antenna_PCV(index_master), 90-elM(sat_pr), azM(sat_pr), sys, 1);
        pr1_M(sat_pr) = pr1_M(sat_pr) - (PCO1_M + PCV1_M);
        ph1_M(sat)    = ph1_M(sat)    - (PCO1_M(index_ph) + PCV1_M(index_ph))./lambda(sat,1);

        if (length(frequencies) == 2 || frequencies(1) == 2)
            PCO2_M = PCO_correction(antenna_PCV(index_master), XR0, XS, sys, 2);
            PCV2_M = PCV_correction(antenna_PCV(index_master), 90-elM(sat_pr), azM(sat_pr), sys, 2);
            pr2_M(sat_pr) = pr2_M(sat_pr) - (PCO2_M + PCV2_M);
            ph2_M(sat)    = ph2_M(sat)    - (PCO2_M(index_ph) + PCV2_M(index_ph))./lambda(sat,2);
        end
    end

    if (~isempty(antenna_PCV) && antenna_PCV(1).n_frequency ~= 0) % rover
        index_rover = 1;
        PCO1_R = PCO_correction(antenna_PCV(index_rover), XR0, XS, sys, 1);
        PCV1_R = PCV_correction(antenna_PCV(index_rover), 90-elR(sat_pr), azR(sat_pr), sys, 1);
        pr1_R(sat_pr) = pr1_R(sat_pr) - (PCO1_R + PCV1_R);
        ph1_R(sat)    = ph1_R(sat)    - (PCO1_R(index_ph) + PCV1_R(index_ph))./lambda(sat,1);

        if (length(frequencies) == 2 || frequencies(1) == 2)
            PCO1_R = PCO_correction(antenna_PCV(index_rover), XR0, XS, sys, 2);
            PCV2_R = PCV_correction(antenna_PCV(index_rover), 90-elM(sat_pr), azM(sat_pr), sys, 2);
            pr2_R(sat_pr) = pr2_R(sat_pr) - (PCO1_R + PCV2_R);
            ph2_R(sat)    = ph2_R(sat)    - (PCO1_R(index_ph) + PCV2_R(index_ph))./lambda(sat,2);
        end
    end

    %--------------------------------------------------------------------------------------------
    % PREPARE INPUT FOR LEAST SQUARES BATCH
    %--------------------------------------------------------------------------------------------

    %if at least min_nsat_LS satellites are available after the cutoffs, and if the
    % condition number in the least squares does not exceed the threshold
    if (size(sat,1) >= min_nsat_LS && (isempty(cond_num) || cond_num < cond_num_threshold))

        if (frequencies == 1)
            [y0_epo, A_epo, b_epo, Q_epo] = input_BLK_DD_code_phase_static(XR0, XM, XS, pr1_R(sat), ph1_R(sat), snr_R(sat), pr1_M(sat), ph1_M(sat), snr_M(sat), elR(sat), elM(sat), err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda(sat,1));
        else
            [y0_epo, A_epo, b_epo, Q_epo] = input_BLK_DD_code_phase_static(XR0, XM, XS, pr2_R(sat), ph2_R(sat), snr_R(sat), pr2_M(sat), ph2_M(sat), snr_M(sat), elR(sat), elM(sat), err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda(sat,2));
        end

%         if (any(bad_obs))
%             pr_out = bad_obs(bad_obs <= size(sat,1));
%             ph_out = bad_obs(bad_obs >  size(sat,1));
%             only_pr = setdiff(pr_out, ph_out);
%             only_ph = setdiff(ph_out, pr_out);
%
%             conf_sat(sat(pr_out),1) = 0;
%             conf_sat(sat(ph_out),1) = 0;
%             conf_sat(sat(only_pr),1) = +2; %satellite with phase, but without code
%             conf_sat(sat(only_ph),1) = -1;
%         end
    else
        pivot = 0;
    end
else
    pivot = 0;
end
