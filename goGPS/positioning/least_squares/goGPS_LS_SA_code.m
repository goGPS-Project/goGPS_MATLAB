function goGPS_LS_SA_code(time_rx, pr1, pr2, snr, Eph, SP3, iono, sbas, lambda, frequencies, p_rate, obs_comb, XR0)

% SYNTAX:
%   goGPS_LS_SA_code(time_rx, pr1, pr2, snr, Eph, SP3, iono, sbas, frequencies,  p_rate, obs_comb, XR0);
%
% INPUT:
%   time_rx = GPS reception time
%   pr1     = code observations (L1 carrier)
%   pr2     = code observations (L2 carrier)
%   snr     = signal-to-noise ratio
%   Eph     = satellite ephemeris
%   SP3     = structure containing precise ephemeris data
%   iono    = ionosphere parameters
%   sbas    = SBAS corrections
%   lambda  = wavelength matrix (depending on the enabled constellations)
%   frequencies = L1 carrier (phase=1), L2 carrier (phase=2)
%   p_rate  = processing interval [s]
%   obs_comb = observations combination (e.g. iono-free: obs_comb = 'IONO_FREE')
%   XR0     = approximated receiver position
%
% DESCRIPTION:
%   Computation of the receiver position (X,Y,Z).
%   Standalone code positioning by least squares adjustment.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
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

global sigmaq0
global cutoff snr_threshold cond_num_threshold o1 o2 o3

global Xhat_t_t Cee conf_sat conf_cs pivot pivot_old
global azR elR distR
global PDOP HDOP VDOP

global residuals_float outliers
residuals_float(:)=NaN;
outliers(:)=NaN;

%covariance matrix initialization
cov_XR = [];

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1,1);

%topocentric coordinate initialization
azR   = zeros(nSatTot,1);
elR   = zeros(nSatTot,1);
distR = zeros(nSatTot,1);

%iono-free coefficients
alpha1 = lambda(:,4);
alpha2 = lambda(:,5);

%available satellites
if (length(frequencies) == 2)
    sat = find( (pr1(:,1) ~= 0) & (pr2(:,1) ~= 0) );
else
    if (frequencies == 1)
        sat = find(pr1(:,1) ~= 0);
    else
        sat = find(pr2(:,1) ~= 0);
    end
end
if (isempty(SP3))
    eph_avail = Eph(30,:);
else
    eph_avail = SP3.avail;
end
sat = sat(ismember(sat, eph_avail));

%--------------------------------------------------------------------------------------------
% SBAS FAST CORRECTIONS
%--------------------------------------------------------------------------------------------

if (~isempty(sbas))
    %apply SBAS fast (pseudorange) corrections
    pr1(sat) = pr1(sat) + sbas.prc(sat)';
end

%--------------------------------------------------------------------------------------------
% POSITIONING
%--------------------------------------------------------------------------------------------

min_nsat = 4;

if any(XR0)
    flag_XR = 1;
%    flag_XR = 2;
else
    flag_XR = 0;
end

if (size(sat,1) >= min_nsat)

    if (frequencies(1) == 1)
        if (length(frequencies) < 2 || ~strcmp(obs_comb,'IONO_FREE'))
            [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, elR(sat), azR(sat), distR(sat), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, obs_outlier, ~, ~, residuals] = init_positioning(time_rx, pr1(sat), snr(sat), Eph, SP3, iono, sbas, XR0, [], [], sat, [], lambda(sat,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 0, 1); %#ok<ASGLU>
        else
            [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, elR(sat), azR(sat), distR(sat), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, obs_outlier, ~, ~, residuals] = init_positioning(time_rx, alpha1(sat).*pr1(sat) - alpha2(sat).*pr2(sat), snr(sat), Eph, SP3, zeros(8,1), sbas, XR0, [], [], sat, [], zeros(length(sat),2), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 0, 1); %#ok<ASGLU>
        end
    else
        [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, elR(sat), azR(sat), distR(sat), sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, obs_outlier, ~, ~, residuals] = init_positioning(time_rx, pr2(sat), snr(sat), Eph, SP3, iono, sbas, XR0, [], [], sat, [], lambda(sat,:), cutoff, snr_threshold, frequencies, p_rate, flag_XR, 0, 1); %#ok<ASGLU>
    end

    if ~isempty(dtR)
        residuals_float(residuals(:,2))=residuals(:,1);
        outliers(find(obs_outlier==1))=1;  %#ok<FNDSB>
    end

    %--------------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING
    %--------------------------------------------------------------------------------------------

    %satellite configuration
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat,1) = +1;

    %no cycle-slips when working with code only
    conf_cs = zeros(nSatTot,1);

    %previous pivot
    pivot_old = 0;

    %current pivot
    [null_max_elR, i] = max(elR(sat)); %#ok<ASGLU>
    pivot = sat(i);

    %if less than min_nsat satellites are available after the cutoffs, or if the
    % condition number in the least squares exceeds the threshold
    if (size(sat,1) < min_nsat | cond_num > cond_num_threshold)

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

%-------------------------------------------------------------------------------

Xhat_t_t = zeros(o3,1);
Xhat_t_t(1)    = XR(1);
Xhat_t_t(o1+1) = XR(2);
Xhat_t_t(o2+1) = XR(3);
Cee(:,:) = zeros(o3);
Cee(1,1) = sigma2_XR(1);
Cee(o1+1,o1+1) = sigma2_XR(2);
Cee(o2+1,o2+1) = sigma2_XR(3);
