function [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, el, az, dist, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pseudorange, snr, Eph, SP3_time, SP3_coor, SP3_clck, iono, sbas, XR0, XS0, dtS0, sat0, cutoff_el, cutoff_snr, flag_XR, flag_XS)

% SYNTAX:
%   [XR, dtR, XS, dtS, XS_tx, VS_tx, time_tx, err_tropo, err_iono, sat, el, az, dist, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = init_positioning(time_rx, pseudorange, snr, Eph, SP3_time, SP3_coor, SP3_clck, iono, sbas, XR0, XS0, dtS0, sat0, cutoff_el, cutoff_snr, flag_XR, flag_XS);
%
% INPUT:
%   time_rx     = reception time
%   pseudorange = observed code pseudoranges
%   snr         = observed signal-to-noise ratio
% 
%   Eph         = ephemeris
%   SP3_time    = precise ephemeris time
%   SP3_coor    = precise ephemeris coordinates
%   SP3_clck    = precise ephemeris clocks
%   iono        = ionosphere parameters (Klobuchar)
%   sbas        = SBAS corrections
%   XR0         = receiver position (=[] if not available)
%   XS0         = satellite positions (=[] if not available)
%   dtS0        = satellite clocks (=[] if not available)
%   sat0        = available satellite PRNs
%   cutoff_el   = elevation cutoff
%   cutoff_snr  = signal-to-noise ratio cutoff
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
%   redsat    = satellites available after cutoffs
%   el        = satellite elevation (vector)
%   az        = satellite azimuth (vector)
%   dist      = satellite-receiver geometric distance (vector)
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
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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

global v_light

%----------------------------------------------------------------------------------------------
% FIRST ESTIMATE OF SATELLITE POSITIONS
%----------------------------------------------------------------------------------------------

%number of satellites
nsat = size(pseudorange,1);

%receiver clock error, troposphere and ionosphere initialization (not yet estimated)
dtR = 0;
err_tropo = zeros(nsat,1);
err_iono  = zeros(nsat,1);

if (flag_XS == 0)
    %satellite position and clock error
    [XS, dtS, XS_tx, VS_tx, time_tx, no_eph] = satellite_positions(time_rx, pseudorange, sat0, Eph, SP3_time, SP3_coor, SP3_clck, sbas, err_tropo, err_iono, dtR);
else
    XS  = XS0;
    dtS = dtS0;
    XS_tx = [];
    VS_tx = [];
    time_tx = [];
    no_eph = zeros(nsat,1);
end

%----------------------------------------------------------------------------------------------
% APPROXIMATE RECEIVER POSITION BY BANCROFT ALGORITHM
%----------------------------------------------------------------------------------------------

if (flag_XR == 0)
    
    index = find(no_eph == 0);
    
    %NOTE: satellite selection may enhance the solution
    B = [XS(index,:), pseudorange(index) + v_light * dtS(index)]; %Bancroft matrix
    x = bancroft(B);                                         %estimated parameters
    XR = x(1:3);                                         %receiver coordinates [m]
    %dtR = x(4) / v_light;                               %receiver clock error [s]
else
    XR = XR0;                                          %known receiver coordinates
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
nsat = size(pseudorange,1);
if (flag_XS == 1)
    XS0  = XS0(index,:);
end

%--------------------------------------------------------------------------------------------
% LEAST SQUARES SOLUTION
%--------------------------------------------------------------------------------------------

%there are not enough satellites to improve the receiver position
if ((flag_XR == 1) & (nsat < 4))
    flag_XR = 2;
end

%if the receiver position is fixed, only one satellite is required for estimating the receiver clock
if (flag_XR == 2)
    nsat_required = 1;
else
    nsat_required = 4;
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

    while ((sqrt(sum((XS(:)-XSold(:)).^2)) > threshold) & (n_iter < n_iter_max))

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

        if (flag_XR < 2) %if unknown or approximate receiver position
            [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = LS_SA_code(XR, XS, pseudorange, snr, el, dist, dtS, err_tropo, err_iono);
        else
            [dtR, var_dtR] = LS_SA_code_clock(pseudorange, snr, el, dist, dtS, err_tropo, err_iono);
            cov_XR = [];
            PDOP = -9999;
            HDOP = -9999;
            VDOP = -9999;
            cond_num = [];
        end

        if (flag_XS == 0)
            %satellite position and clock error
            [XS, dtS, XS_tx, VS_tx, time_tx] = satellite_positions(time_rx, pseudorange, sat, Eph, SP3_time, SP3_coor, SP3_clck, sbas, err_tropo, err_iono, dtR);
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
    else
        cov_XR = [];
    end
    var_dtR = [];

    PDOP = -9999;
    HDOP = -9999;
    VDOP = -9999;
    cond_num = [];
end

