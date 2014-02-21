function [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, el, az, dist, sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pseudorange, snr, Eph, SP3, iono, sbas, XR0, XS0, dtS0, sat0, sys0, lambda, cutoff_el, cutoff_snr, phase, flag_XR, flag_XS)

% SYNTAX:
%   [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, el, az, dist, sys, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pseudorange, snr, Eph, SP3, iono, sbas, XR0, XS0, dtS0, sat0, sys0, sys0, lambda, cutoff_el, cutoff_snr, phase, flag_XR, flag_XS);
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
%
% OUTPUT:
%   XR        = receiver position (X,Y,Z)
%   dtR       = receiver clock error (scalar)
%   XS        = satellite position at transmission time in ECEF(time_rx) (X,Y,Z)
%   dtS       = satellite clock error (vector)
%   XS_tx     = satellite position at transmission time in ECEF(time_tx) (X,Y,Z)
%   VS_tx     = satellite velocity at transmission time in ECEF(time_tx) (X,Y,Z)
%   time_tx   = transmission time (vector)
%   err_tropo = tropospheric error (vector)
%   err_iono  = ionospheric error (vector)
%   sat       = satellites available after cutoffs
%   el        = satellite elevation (vector)
%   az        = satellite azimuth (vector)
%   dist      = satellite-receiver geometric distance (vector)
%   sys       = array with different values for different systems
%   cov_XR    = receiver position error covariance matrix [3x3]
%   var_dtR   = receiver clock error variance (scalar)
%   PDOP      = position dilution of precision (scalar)
%   HDOP      = horizontal dilution of precision (scalar)
%   VDOP      = vertical dilution of precision (scalar)
%   cond_num  = condition number from the least squares N matrix (scalar) 
%
% DESCRIPTION:
%   Compute initial receiver and satellite position and clocks using
%   Bancroft and least-squares iterative correction. Requires at least
%   four satellites available.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.2 beta
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

%compute inter-frequency factors (for the ionospheric delay)
ionoFactor = goGNSS.getInterFreqIonoFactor(lambda);

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
var_dtR = [];
PDOP = -9999;
HDOP = -9999;
VDOP = -9999;
cond_num = [];

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

%----------------------------------------------------------------------------------------------
% APPROXIMATE RECEIVER POSITION
%----------------------------------------------------------------------------------------------

if (flag_XR == 0)
    
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

    %iterative least-squares from the center of the Earth (i.e. [0; 0; 0])
    XR0 = zeros(3,1);
    for i = 1 : 3
        [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = LS_SA_code(XR0, XS(index,:), pseudorange(index), zeros(nsat_avail,1), zeros(nsat_avail,1), zeros(nsat_avail,1), dtS(index), zeros(nsat_avail,1), zeros(nsat_avail,1), sys(index));
        XR0 = XR;
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
            [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = LS_SA_code(XR, XS, pseudorange, snr, el, dist, dtS, err_tropo, err_iono, sys);
        else
            [dtR, var_dtR] = LS_SA_code_clock(pseudorange, snr, el, dist, dtS, err_tropo, err_iono, sys);
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
    dtR = [];
    err_tropo = [];
    err_iono  = [];

    if (flag_XR == 2)
        cov_XR = zeros(3,3);
    end
end
