function goGPS_LS_SA_code_MR(time_rx, pr1_R, pr2_R, snr_R, Eph, SP3, iono, sbas, lambda, phase)

% SYNTAX:
%   goGPS_LS_SA_code_MR(time_rx, pr1_R, pr2_R, snr_R, Eph, SP3, iono, sbas, lambda, phase);
%
% INPUT:
%   time_rx = GPS reception time
%   pr1_R = ROVER code observations (L1 carrier)   (nsat x nrec)
%   pr2_R = ROVER code observations (L2 carrier)   (nsat x nrec)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio  (nsat x nrec)
%   Eph   = satellite ephemeris
%   SP3   = structure containing precise ephemeris and clock
%   iono  = ionosphere parameters
%   sbas    = SBAS corrections
%   lambda  = wavelength matrix (depending on the enabled constellations)
%   phase  = L1 carrier (phase=1), L2 carrier (phase=2)
%
% DESCRIPTION:
%   Computation of the receiver position (X,Y,Z).
%   Average of multiple receiver position estimates (separate least squares
%   adjustment) on code observations.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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
global azR elR distR
global PDOP HDOP VDOP

%total number of satellite slots (depending on the constellations enabled)
nSatTot = size(pr1_R,1);

%number of rover receivers
nRov = size(pr1_R,2);

%covariance matrix initialization
cov_XR = [];

%initialization
XR = zeros(3,nRov);
dtR = zeros(nRov,1);

azR   = zeros(nSatTot,nRov);
elR   = zeros(nSatTot,nRov);
distR = zeros(nSatTot,nRov);

err_tropo_R = zeros(nSatTot,nRov);
err_iono_R  = zeros(nSatTot,nRov);

%--------------------------------------------------------------------------------------------
% SATELLITE SELECTION
%--------------------------------------------------------------------------------------------

if (length(phase) == 2)
    sat_pr = find( (prod(single(pr1_R ~= 0),2)) & ...
                   (prod(single(pr2_R ~= 0),2)) );
else
    if (phase == 1)
        sat_pr = find( (prod(single(pr1_R ~= 0),2)) );
    else
        sat_pr = find( (prod(single(pr2_R ~= 0),2)) );
    end
end
sat_pr = sat_pr(ismember(sat_pr, Eph(30,:)));

%--------------------------------------------------------------------------------------------
% SBAS FAST CORRECTIONS
%--------------------------------------------------------------------------------------------

if (~isempty(sbas))
    %apply SBAS fast (pseudorange) corrections
    pr1_R(sat,:) = pr1_R(sat,:) + sbas.prc(ones(nRov,1),sat)';
end

%--------------------------------------------------------------------------------------------
% POSITIONING
%--------------------------------------------------------------------------------------------

if (size(sat_pr,1) >= 4)
    
    sat_pr_R = (1 : nSatTot)';
    for r = 1 : nRov
        if (phase == 1)
            [XR(:,r), dtR(r,1), XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_R_tmp, err_iono_R_tmp, sat_pr_R_tmp, elR(sat_pr_R_tmp,r), azR(sat_pr_R_tmp,r), distR(sat_pr_R_tmp,r), is_GLO, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr1_R(sat_pr,r), snr_R(sat_pr,r), Eph, SP3, iono, sbas, [], [], [], sat_pr, lambda(sat_pr,:), cutoff, snr_threshold, phase, 0, 0); %#ok<ASGLU>
        else
            [XR(:,r), dtR(r,1), XS, dtS, XS_tx, VS_tx, time_tx, err_tropo_R_tmp, err_iono_R_tmp, sat_pr_R_tmp, elR(sat_pr_R_tmp,r), azR(sat_pr_R_tmp,r), distR(sat_pr_R_tmp,r), is_GLO, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pr2_R(sat_pr,r), snr_R(sat_pr,r), Eph, SP3, iono, sbas, [], [], [], sat_pr, lambda(sat_pr,:), cutoff, snr_threshold, phase, 0, 0); %#ok<ASGLU>
        end
        
        err_tropo_R(sat_pr_R_tmp,r) = err_tropo_R_tmp;
        err_iono_R(sat_pr_R_tmp,r) = err_iono_R_tmp;
        sat_pr_R = intersect(sat_pr_R, sat_pr_R_tmp);
    end
    
    %keep only satellites above the elevation cutoff/snr threshold
    sat_pr = sat_pr_R;
    
    %--------------------------------------------------------------------------------------------
    % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
    %--------------------------------------------------------------------------------------------
    
    %satellite configuration
    conf_sat = zeros(nSatTot,1);
    conf_sat(sat_pr,1) = +1;
    
    %no cycle-slips when working with code only
    conf_cs = zeros(nSatTot,1);
    
    %previous pivot
    pivot_old = 0;
    
    %actual pivot
    [~, pivot_index] = max(elR(sat_pr));
    pivot = sat_pr(pivot_index);
    
    %--------------------------------------------------------------------------------------------
    % AVERAGING OF MULTIPLE SOLUTIONS
    %--------------------------------------------------------------------------------------------

    XR = mean(XR,2);
       
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
