function goGPS_LS_DD_code(time_rx, XM, pr1_R, pr1_M, pr2_R, pr2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase)

% SYNTAX:
%   goGPS_LS_DD_code(time_rx, XM, pr1_R, pr1_M, pr2_R, pr2_M, snr_R, snr_M, Eph, SP3, iono, lambda, phase);
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
%   SP3   = structure containing precise ephemeris data
%   iono  = ionosphere parameters
%   lambda = wavelength matrix (depending on the enabled constellations)
%   phase  = L1 carrier (phase=1), L2 carrier (phase=2)
%
% DESCRIPTION:
%   Computation of the receiver position (X,Y,Z).
%   Relative (double difference) positioning by least squares adjustment
%   on code observations.

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

global sigmaq0
global cutoff snr_threshold cond_num_threshold o1 o2 o3

global Xhat_t_t Cee conf_sat conf_cs pivot pivot_old
global azR elR distR azM elM distM
global PDOP HDOP VDOP

%covariance matrix initialization
cov_XR = [];

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1_R,1);

%compute inter-frequency factors (for the ionospheric delay)
ionoFactor = goGNSS.getInterFreqIonoFactor(lambda);

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
    sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) & ...
                (pr2_R ~= 0) & (pr2_M ~= 0) );
else
    if (phase == 1)
        sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) );
    else
        sat = find( (pr2_R ~= 0) & (pr2_M ~= 0) );
    end
end
if (isempty(SP3))
    eph_avail = Eph(30,:);
else
    eph_avail = SP3.avail;
end
sat = sat(ismember(sat, eph_avail));

min_nsat = 4;

if (size(sat,1) >= min_nsat)
    
    if (phase(1) == 1)
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_M, elM(sat_M), azM(sat_M), distM(sat_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr1_M(sat),   snr_M(sat),   Eph, SP3, iono, [], XM, [],  [], sat,   [], lambda(sat,:),   cutoff, snr_threshold, phase, 2, 0); %#ok<NASGU,ASGLU>
        if (length(sat_M) < min_nsat); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_R, elR(sat_R), azR(sat_R), distR(sat_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1_R(sat_M), snr_R(sat_M), Eph, SP3, iono, [], [], XS, dtS, sat_M, sys, lambda(sat_M,:), cutoff, snr_threshold, phase, 0, 1); %#ok<ASGLU>
    else
        [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_M, err_iono_M, sat_M, elM(sat_M), azM(sat_M), distM(sat_M), sys, cov_XM, var_dtM]                             = init_positioning(time_rx, pr2_M(sat),   snr_M(sat),   Eph, SP3, iono, [], XM, [],  [], sat,   [], lambda(sat,:),   cutoff, snr_threshold, phase, 2, 0); %#ok<NASGU,ASGLU>
        if (length(sat_M) < min_nsat); return; end
        [XR, dtR, XS, dtS,     ~,     ~,       ~, err_tropo_R, err_iono_R, sat_R, elR(sat_R), azR(sat_R), distR(sat_R), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2_R(sat_M), snr_R(sat_M), Eph, SP3, iono, [], [], XS, dtS, sat_M, sys, lambda(sat_M,:), cutoff, snr_threshold, phase, 0, 1); %#ok<ASGLU>
    end
    
    %keep only satellites that rover and master have in common
    [sat, iR, iM] = intersect(sat_R, sat_M);
    XS = XS(iR,:);
    if (~isempty(err_tropo_R))
        err_tropo_R = err_tropo_R(iR);
        err_iono_R  = err_iono_R (iR);
        err_tropo_M = err_tropo_M(iM);
        err_iono_M  = err_iono_M (iM);
    end
    
    %--------------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
    %--------------------------------------------------------------------------------------------
    
    %satellite configuration
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat,1) = +1;
    
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
    
    %if a sufficient number of satellites is available after the cutoffs, and if the 
    % condition number in the least squares does not exceed the threshold
    if (size(sat,1) >= min_nsat && cond_num < cond_num_threshold)
        
        %loop is needed to improve the atmospheric error correction
        for i = 1 : 3

            if (phase == 1)
                [XR, cov_XR] = LS_DD_code(XR, XS, pr1_R(sat), pr1_M(sat), snr_R(sat), snr_M(sat), elR(sat), elM(sat), distR(sat), distM(sat), err_tropo_R, err_tropo_M, err_iono_R, err_iono_M, pivot_index);
            else
                [XR, cov_XR] = LS_DD_code(XR, XS, pr2_R(sat), pr2_M(sat), snr_R(sat), snr_M(sat), elR(sat), elM(sat), distR(sat), distM(sat), err_tropo_R, err_tropo_M, err_iono_R, err_iono_M, pivot_index);
            end
            
            [phiR, lamR, hR] = cart2geod(XR(1), XR(2), XR(3));
            [azR(azR ~= 0), elR(elR ~= 0), distR(distR ~= 0)] = topocent(XR, XS);
            
            err_tropo_R = tropo_error_correction(elR(elR ~= 0), hR);
            err_iono_R = iono_error_correction(phiR*180/pi, lamR*180/pi, azR(azR ~= 0), elR(elR ~= 0), time_rx, iono, []);
            
            %correct the ionospheric errors for different frequencies
            err_iono_R = ionoFactor(sat,phase).*err_iono_R;
        end
    else
        if (~isempty(Xhat_t_t))
            XR = Xhat_t_t([1,o1+1,o2+1]);
            pivot = 0;
        else
            return
        end
    end
else
    if (~isempty(Xhat_t_t))
        XR = Xhat_t_t([1,o1+1,o2+1]);
        pivot = 0;
    else
        return
    end
end

if isempty(cov_XR) %if it was not possible to compute the covariance matrix
    cov_XR = sigmaq0 * eye(3);
end
sigma2_XR = diag(cov_XR);

%--------------------------------------------------------------------------------------------
% GOGPS OUTPUT SAVING
%--------------------------------------------------------------------------------------------

Xhat_t_t = zeros(o3,1);
Xhat_t_t(1)    = XR(1);
Xhat_t_t(o1+1) = XR(2);
Xhat_t_t(o2+1) = XR(3);
Cee(:,:) = zeros(o3);
Cee(1,1) = sigma2_XR(1);
Cee(o1+1,o1+1) = sigma2_XR(2);
Cee(o2+1,o2+1) = sigma2_XR(3);
