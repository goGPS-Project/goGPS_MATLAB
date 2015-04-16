function [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, el, az, dist, sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, obs_outlier, bad_epoch, var_SPP, residuals_obs, is_bias] = init_positioning(time_rx, pseudorange, snr, Eph, SP3, iono, sbas, XR0, XS0, dtS0, sat0, sys0, lambda, cutoff_el, cutoff_snr, phase, flag_XR, flag_XS, flag_OOLO)

% SYNTAX:
%   [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, el, az, dist, sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_sat, bad_epoch] = init_positioning(time_rx, pseudorange, snr, Eph, SP3, iono, sbas, XR0, XS0, dtS0, sat0, sys0, sys0, lambda, cutoff_el, cutoff_snr, phase, flag_XR, flag_XS, flag_OOLO);
%
% INPUT:
%   time_rx     = reception time
%   pseudorange = observed code pseudoranges
%   snr         = observed signal-to-noise ratio
% 
%   Eph         = ephemeris
%   SP3         = structure containing precise ephemeris data
%   iono        = ionosphere parameters (Klobuchar)
%   sbas        = SBAS corrections
%   XR0         = receiver position (=[] if not available)
%   XS0         = satellite positions (=[] if not available)
%   dtS0        = satellite clocks (=[] if not available)
%   sat0        = available satellite PRNs
%   sys0        = array with different values for different systems (used if flag_XS == 1)
%   lambda      = wavelength matrix (depending on the enabled constellations)
%   cutoff_el   = elevation cutoff
%   cutoff_snr  = signal-to-noise ratio cutoff
%   phase       = L1 carrier (phase=1), L2 carrier (phase=2)
%   flag_XR     = 0: unknown
%                 1: approximated
%                 2: fixed
%   flag_XS     = 0: unknown
%                 1: already estimated
%   flag_OOLO   = 0: outlier detection enabled
%                 1: outlier detection disabled
%
% OUTPUT:
%   XR          = receiver position (X,Y,Z)
%   dtR         = receiver clock error (scalar)
%   XS          = satellite position at transmission time in ECEF(time_rx) (X,Y,Z)
%   dtS         = satellite clock error (vector)
%   XS_tx       = satellite position at transmission time in ECEF(time_tx) (X,Y,Z)
%   VS_tx       = satellite velocity at transmission time in ECEF(time_tx) (X,Y,Z)
%   time_tx     = transmission time (vector)
%   err_tropo   = tropospheric error (vector)
%   err_iono    = ionospheric error (vector)
%   sat         = satellites available after cutoffs
%   el          = satellite elevation (vector)
%   az          = satellite azimuth (vector)
%   dist        = satellite-receiver geometric distance (vector)
%   sys         = array with different values for different systems
%   cov_XR      = receiver position error covariance matrix [3x3]
%   var_dtR     = receiver clock error variance (scalar)
%   PDOP        = position dilution of precision (scalar)
%   HDOP        = horizontal dilution of precision (scalar)
%   VDOP        = vertical dilution of precision (scalar)
%   cond_num    = condition number from the least squares N matrix (scalar) 
%   obs_outlier = vector with 0 or 1 for sats found as outlier
%   bad_epoch   = 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   var_SPP     = [code single point positioning a posteriori sigma, sum of
%                weighted squared residuals, redoundancy]
%   residuals_obs=vector with [residuals of all input observation, computed from the final estimates, correspondent sat id]
%   is_bias     = inter-systems bias (vector with all possibile systems)
%
% DESCRIPTION:
%   Compute initial receiver and satellite position and clocks using
%   Bancroft and least-squares iterative correction. Requires at least
%   four satellites available.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3 beta
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Stefano Caldera
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

%compute inter-frequency factors (for the ionospheric delay)
ionoFactor = goGNSS.getInterFreqIonoFactor(lambda);

bad_epoch=NaN;
var_SPP = NaN(1,3);
is_bias=NaN(6,1);

%----------------------------------------------------------------------------------------------
% FIRST ESTIMATE OF SATELLITE POSITIONS
%----------------------------------------------------------------------------------------------
%number of satellites
nsat = size(pseudorange,1);

%receiver clock error, troposphere and ionosphere initialization (not yet estimated)
dtR = 0;
err_tropo = zeros(nsat,1);
err_iono  = zeros(nsat,1);

%output variable initialization
cov_XR = [];
var_dtR = NaN;
PDOP = -9999;
HDOP = -9999;
VDOP = -9999;
cond_num = [];
obs_outlier = [];

if (flag_XS == 0)
    %satellite position and clock error
    [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, pseudorange, sat0, Eph, SP3, sbas, err_tropo, err_iono, dtR);

else
    XS  = XS0;
    dtS = dtS0;
    XS_tx = [];
    VS_tx = [];
    time_tx = [];
    no_eph = zeros(nsat,1);
    sys = sys0;
end

%if multi-system observations, then an inter-system bias parameter for each additional system must be estimated
num_sys  = length(unique(sys(sys ~= 0)));
min_nsat = 3 + num_sys;

% maximum RMS of code single point positioning to accept current epoch
SPP_threshold=4; %meters 
bad_sat=[];

%----------------------------------------------------------------------------------------------
% APPROXIMATE RECEIVER POSITION
%----------------------------------------------------------------------------------------------

if (isempty(XR0))
    XR0 = zeros(3,1);
end

if (flag_XR < 2 || ~any(XR0))
    
    index = find(no_eph == 0);
    
    nsat_avail = length(index);

    if (nsat_avail < min_nsat) %if available observations are not enough, return empty variables
        XR   = [];
        dtR  = [];
        az   = [];
        el   = [];
        dist = [];
        sat  = [];
        err_tropo = [];
        err_iono  = [];
        
        if (flag_XR == 2)
            cov_XR = zeros(3,3);
        end
        
        return
    end

    %iterative least-squares from XR0,i.e. given coordinates or the center of the Earth (i.e. [0; 0; 0])
    n_iter_max = 5;
    n_iter = 0;
    var_SPP(1) = Inf;
    while(var_SPP(1) > SPP_threshold^2 && n_iter < n_iter_max)
        [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP] = LS_SA_code(XR0, XS(index,:), pseudorange(index), zeros(nsat_avail,1), zeros(nsat_avail,1), zeros(nsat_avail,1), dtS(index), zeros(nsat_avail,1), zeros(nsat_avail,1), sys(index), SPP_threshold);       
        %bad_sat(sat(bad_obs))=1;
        XR0 = XR;
        n_iter = n_iter + 1;
    end
else
    XR = XR0; %known receiver coordinates
end

%----------------------------------------------------------------------------------------------
% ELEVATION CUTOFF, SNR CUTOFF AND REMOVAL OF SATELLITES WITHOUT EPHEMERIS
%----------------------------------------------------------------------------------------------

%satellite topocentric coordinates (azimuth, elevation, distance)
[az, el, dist] = topocent(XR, XS);

%elevation cutoff, SNR cutoff and removal of satellites without ephemeris
if (any(snr))
    index = find((el > cutoff_el) & ((snr ~= 0) & (snr > cutoff_snr)) & (no_eph == 0));
else
    index = find((el > cutoff_el) & (no_eph == 0));
end
sat   = sat0(index);
pseudorange = pseudorange(index);
snr  = snr(index);
el   = el(index);
az   = az(index);
dist = dist(index);
XS   = XS(index,:);
dtS  = dtS(index);
sys  = sys(index);
ionoFactor = ionoFactor(index,:);
nsat = size(pseudorange,1);
if (flag_XS == 1)
    XS0  = XS0(index,:);
end
residuals_obs=NaN(nsat,2);

%--------------------------------------------------------------------------------------------
% LEAST SQUARES SOLUTION
%--------------------------------------------------------------------------------------------

%if multi-system observations, then an inter-system bias parameter for each additional system must be estimated
num_sys  = length(unique(sys(sys ~= 0)));
min_nsat = 3 + num_sys;

%there are not enough satellites to improve the receiver position
if ((flag_XR == 1) && (nsat < min_nsat))
    flag_XR = 2;
end

%if the receiver position is fixed, only one satellite is required for estimating the receiver clock
if (flag_XR == 2)
    nsat_required = min_nsat - 3;
else
    nsat_required = min_nsat;
end

if (nsat >= nsat_required)

    XSold = zeros(nsat,3);      %old version of satellite positions
    n_iter = 0;                 %iteration counter
    
    if (flag_XS == 0)
        n_iter_max = 2;         %maximum number of iterations (for estimating satellite positions)
    else
        n_iter_max = 1;         %no iterations (satellite positions already estimated)
    end

    %threshold = sum of the coordinate variation over all satellites [m]
    threshold = 0.01;
    index_obs=1:length(pseudorange);

    while ((sqrt(sum((XS(:)-XSold(:)).^2)) > threshold) && (n_iter < n_iter_max))

        XSold = XS;             %save old version of satellite positions
        n_iter = n_iter+1;      %increase iteration counter

        %cartesian to geodetic conversion of ROVER coordinates
        [phiR, lamR, hR] = cart2geod(XR(1), XR(2), XR(3));

        %radians to degrees
        phiR = phiR * 180 / pi;
        lamR = lamR * 180 / pi;

        %computation of tropospheric errors
        err_tropo = tropo_error_correction(el, hR);

        %computation of ionospheric errors
        err_iono = iono_error_correction(phiR, lamR, az, el, time_rx, iono, sbas);      
        
        %correct the ionospheric errors for different frequencies
        err_iono = ionoFactor(:,phase).*err_iono;

        if (flag_XR < 2) %if unknown or approximate receiver position

            [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP, residuals_obs(index_obs,1), is_bias] = LS_SA_code(XR, XS, pseudorange, snr, el, dist, dtS, err_tropo, err_iono, sys, SPP_threshold);
            residuals_obs(index_obs,2)=sat;
            if (exist('flag_OOLO','var') && flag_OOLO==1)
                if ~isempty(bad_obs)
                    % add outlier satellite to obs_outlier
                    obs_outlier(sat(bad_obs))=1;
                    
                    % remove bad satellite from observations
                    XS(bad_obs,:)=[];
                    XSold(bad_obs,:)=[];
                    snr(bad_obs)=[];
                    sat(bad_obs)=[];
                    pseudorange(bad_obs)=[];
                    dtS(bad_obs)=[];
                    err_tropo(bad_obs)=[];
                    err_iono(bad_obs)=[];
                    sys(bad_obs)=[];
                    ionoFactor(bad_obs,:)=[];
                    index_obs(bad_obs)=[];
                    nsat = size(pseudorange,1);
                    if (flag_XS == 1)
                        XS0(bad_obs,:)  = [];
                    end
                    
                    %satellite topocentric coordinates (azimuth, elevation, distance)
                    [az, el, dist] = topocent(XR, XS);
                    
                    %elevation cutoff
                    index = find((el > cutoff_el));
                    
                    sat   = sat(index);
                    pseudorange = pseudorange(index);
                    snr  = snr(index);
                    el   = el(index);
                    az   = az(index);
                    dist = dist(index);
                    XS   = XS(index,:);
                    XSold= XSold(index,:);
                    dtS  = dtS(index);
                    err_tropo=err_tropo(index);
                    err_iono=err_iono(index);
                    sys  = sys(index);
                    ionoFactor = ionoFactor(index,:);
                    nsat = size(pseudorange,1);
                    if (flag_XS == 1)
                        XS0  = XS0(index,:);
                    end
                end
            end
            
        else
            [dtR, var_dtR, bad_obs, bad_epoch, var_SPP, residuals_obs(index_obs,1), is_bias] = LS_SA_code_clock(pseudorange, snr, el, dist, dtS, err_tropo, err_iono, sys, SPP_threshold);
            residuals_obs(index_obs,2)=sat;
            if (exist('flag_OOLO','var') && flag_OOLO==1)
                if ~isempty(bad_obs)
                    % add outlier satellite to obs_outlier
                    obs_outlier(sat(bad_obs))=1;
                    % remove bad satellite from observations
                    XS(bad_obs,:)=[];
                    XSold(bad_obs,:)=[];
                    snr(bad_obs)=[];
                    sat(bad_obs)=[];
                    pseudorange(bad_obs)=[];
                    dtS(bad_obs)=[];
                    err_tropo(bad_obs)=[];
                    err_iono(bad_obs)=[];
                    sys(bad_obs)=[];
                    ionoFactor(bad_obs,:)=[];
                    index_obs(bad_obs)=[];
                    nsat = size(pseudorange,1);
                    if (flag_XS == 1)
                        XS0(bad_obs,:)  = [];
                    end
                    
                    %satellite topocentric coordinates (azimuth, elevation, distance)
                    [az, el, dist] = topocent(XR, XS);
                    
                    %elevation cutoff
                    index = find((el > cutoff_el));
                    
                    sat   = sat(index);
                    pseudorange = pseudorange(index);
                    snr  = snr(index);
                    el   = el(index);
                    az   = az(index);
                    dist = dist(index);
                    XS   = XS(index,:);
                    XSold= XSold(index,:);
                    dtS  = dtS(index);
                    err_tropo=err_tropo(index);
                    err_iono=err_iono(index);
                    sys  = sys(index);
                    ionoFactor = ionoFactor(index,:);
                    nsat = size(pseudorange,1);
                    if (flag_XS == 1)
                        XS0  = XS0(index,:);
                    end
                end
            end
        end

        if (flag_XS == 0)
            %satellite position and clock error
            [XS, dtS, XS_tx, VS_tx, time_tx] = satellite_positions(time_rx, pseudorange, sat, Eph, SP3, sbas, err_tropo, err_iono, dtR);
        else
            XS  = XS0;
            dtS = dtS0;
        end

        %satellite topocentric coordinates (azimuth, elevation, distance)
        [az, el, dist] = topocent(XR, XS);
    end
    
else
    %empty variables
    
    dtR  = [];
    az   = [];
    el   = [];
    dist = [];
    sat  = [];
    err_tropo = [];
    err_iono  = [];
    
    if (flag_XR == 2)
        cov_XR = zeros(3,3);
    else
        XR   = [];
    end
    
end
