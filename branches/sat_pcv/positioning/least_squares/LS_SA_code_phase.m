function [XR, dtR, N_hat, cov_XR, var_dtR, cov_N, PDOP, HDOP, VDOP] = LS_SA_code_phase(XR_approx, XS, pr, ph, snr, elR, distR_approx, sat_pr, sat_ph, dtS, err_tropo, err_iono, sys, lambda)

% SYNTAX:
%   [XR, dtR, N_hat, cov_XR, var_dtR, cov_N, PDOP, HDOP, VDOP] = LS_SA_code_phase(XR_approx, XS, pr, ph, snr, elR, distR_approx, sat_pr, sat_ph, dtS, err_tropo, err_iono, sys, lambda);
%
% INPUT:
%   XR_approx    = receiver approximate position (X,Y,Z)
%   XS           = satellite position (X,Y,Z) (with both code and phase)
%   pr           = code observations
%   ph           = phase observations
%   snr          = signal-to-noise ratio
%   elR          = satellite elevation (vector)
%   distR_approx = approximate receiver-satellite distance (vector)
%   sat_pr       = available satellites
%   sat_ph       = available satellites with phase
%   dtS          = satellite clock error (vector)
%   err_tropo    = tropospheric error
%   err_iono     = ionospheric error
%   sys          = array with different values for different systems
%   lambda       = vector containing GNSS wavelengths for available satellites
%
% OUTPUT:
%   pos_R = estimated position (X,Y,Z)
%   cov_pos_R = covariance matrix of estimation errors (rover position)
%   N_stim = linear combination of ambiguity estimate
%   cov_N_stim = covariance matrix of estimation errors (combined ambiguity values)
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%
% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

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

%variable initialization
global sigmaq_cod1 sigmaq_ph

v_light = goGNSS.V_LIGHT;

%data indexes
[~, index] = intersect(sat_pr,sat_ph); %sat_ph is a subset of sat_pr

%number of observations (assuming that sat_ph is a subset of sat_pr)
nsat_pr = length(sat_pr);
nsat_ph = length(sat_ph);
n = nsat_pr + nsat_ph;

%number of unknown parameters
m = 4 + nsat_ph;

% %approximate receiver-satellite distance
% XR_mat = XR_approx(:,ones(n,1))';
% distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));

%design matrix (code)
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ...   %column for X coordinate
     (XR_approx(2) - XS(:,2)) ./ distR_approx, ...   %column for Y coordinate
     (XR_approx(3) - XS(:,3)) ./ distR_approx, ...   %column for Z coordinate
      zeros(nsat_pr,nsat_ph), ... %column for phase ambiguities   (here zero)
      ones(nsat_pr,1)];    %column for receiver clock delay (multiplied by c)

%design matrix (phase)
A = [A; (XR_approx(1) - XS(index,1)) ./ distR_approx(index), ... %column for X coordinate
        (XR_approx(2) - XS(index,2)) ./ distR_approx(index), ... %column for Y coordinate
        (XR_approx(3) - XS(index,3)) ./ distR_approx(index), ... %column for Z coordinate
         diag(-lambda(index)) .* eye(nsat_ph), ...               %column for phase ambiguities
         ones(nsat_ph,1)];             %column for receiver clock delay (multiplied by c)
     
%if multi-system observations, then estimate an inter-system bias parameter for each additional system
uni_sys = unique(sys(sys ~= 0));
num_sys = length(uni_sys);
ISB = zeros(n,1);
if (num_sys > 1)
    m = m + num_sys - 1;
    for s = 2 : num_sys
        ISB(sys == uni_sys(s)) = 1;
        A = [A, ISB];
        ISB = zeros(n,1);
    end
end

%known term vector
b_pr = distR_approx - v_light*dtS + err_tropo + err_iono; %code
b_ph = distR_approx - v_light*dtS + err_tropo - err_iono; %phase
b = [b_pr; b_ph(index)];

%observation vector
y0 = [pr; lambda(index).*ph(index)];

%observation noise covariance matrix
Q = zeros(n);
Q1 = cofactor_matrix_SA(elR, snr);
Q2 = Q1(index,index);
Q(1:nsat_pr,1:nsat_pr) = sigmaq_cod1 * Q1;
Q(nsat_pr+1:end,nsat_pr+1:end) = sigmaq_ph * Q2;

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x_hat = (N^-1)*A'*(Q^-1)*(y0-b);
XR = XR_approx + x_hat(1:3);

%estimated phase ambiguities
N_hat = x_hat(3+[1:nsat_ph]);

%estimated receiver clock
dtR = x_hat(3+nsat_ph+1) / v_light;

%estimation of the variance of the observation error
y_hat = A*x_hat + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat * (N^-1);
    cov_XR  = Cxx(1:3,1:3);
    cov_N   = Cxx(3+[1:nsat_ph],3+[1:nsat_ph]);
    var_dtR = Cxx(3+nsat_ph+1,3+nsat_ph+1) / v_light^2;
else
    cov_XR  = [];
    cov_N   = [];
    var_dtR = []; 
end

%DOP computation
if (nargout > 6)
    A(nsat_pr+1:end,:) = []; A(:,4:end-1) = [];
    cov_XYZ = (A'*A)^-1;
    cov_XYZ = cov_XYZ(1:3,1:3);
    cov_ENU = global2localCov(cov_XYZ, XR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end
