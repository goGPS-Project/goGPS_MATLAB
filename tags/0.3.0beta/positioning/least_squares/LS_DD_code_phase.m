function [XR, N_hat, cov_XR, cov_N, PDOP, HDOP, VDOP] = LS_DD_code_phase ...
         (XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, phase)

% SYNTAX:
%   [XR, N_hat, cov_XR, cov_N, PDOP, HDOP, VDOP] = LS_DD_code_phase ...
%   (XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, phase);
%
% INPUT:
%   XR_approx   = receiver approximate position (X,Y,Z)
%   XM          = master station position (X,Y,Z)
%   XS          = satellite position (X,Y,Z)
%   pr_R        = receiver code observations
%   ph_R        = receiver phase observations
%   pr_M        = master code observations
%   pr_M        = master phase observations
%   snr_R       = receiversignal-to-noise ratio
%   snr_M       = mastersignal-to-noise ratio
%   elR         = satellite elevation (vector)
%   elM         = satellite elevation (vector)
%   err_tropo_R = tropospheric error
%   err_tropo_M = tropospheric error
%   err_iono_R  = ionospheric error
%   err_iono_M  = ionospheric error
%   pivot_index  = index identifying the pivot satellite
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% OUTPUT:
%   XR = estimated position (X,Y,Z)
%   N_hat = linear combination of ambiguity estimate
%   cov_XR = covariance matrix of estimation errors (rover position)
%   cov_N = covariance matrix of estimation errors (combined ambiguity values)
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision

%
% DESCRIPTION:
%   Least squares solution using code and phase double differences.
%   Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
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

%variable initialization
global lambda1 lambda2
global sigmaq_cod1 sigmaq_ph

if (phase == 1)
    lambda = lambda1;
else
    lambda = lambda2;
end

%number of observations
n = 2*length(pr_R);

%number of unknown parameters
m = 3 + (n/2 - 1);

%approximate receiver-satellite distance
XR_mat = XR_approx(:,ones(n/2,1))';
XM_mat = XM(:,ones(n/2,1))';
distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
distM = sqrt(sum((XS-XM_mat).^2 ,2));

%design matrix (code)
A = [((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index,1)) / distR_approx(pivot_index)), ... %column for X coordinate
     ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index,2)) / distR_approx(pivot_index)), ... %column for Y coordinate
     ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index,3)) / distR_approx(pivot_index)), ... %column for Z coordinate
       zeros(n/2,n/2)]; %column for phase ambiguities   (here zero)
 
%design matrix (phase)
A = [A; ((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index,1)) / distR_approx(pivot_index)), ... %column for X coordinate
        ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index,2)) / distR_approx(pivot_index)), ... %column for Y coordinate
        ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index,3)) / distR_approx(pivot_index)), ... %column for Z coordinate
          -lambda * eye(n/2)]; %column for phase ambiguities

%known term vector
b_pr =        (distR_approx - distM)      - (distR_approx(pivot_index) - distM(pivot_index));       %approximate pseudorange DD
b_pr = b_pr + (err_tropo_R - err_tropo_M) - (err_tropo_R(pivot_index)  - err_tropo_M(pivot_index)); %tropospheric error DD
b_pr = b_pr + (err_iono_R  - err_iono_M)  - (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD (code)
b_ph = b_pr - (err_iono_R  - err_iono_M)  - (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD (phase)
b = [b_pr; b_ph];

%observation vector
y0_pr =         (pr_R - pr_M) - (pr_R(pivot_index) - pr_M(pivot_index));
y0_ph = lambda*((ph_R - ph_M) - (ph_R(pivot_index) - ph_M(pivot_index)));
y0 = [y0_pr; y0_ph];

%remove pivot-pivot lines
A( [pivot_index, pivot_index+n/2], :) = [];
A(            :, pivot_index+3)       = [];
b( [pivot_index, pivot_index+n/2])    = [];
y0([pivot_index, pivot_index+n/2])    = [];
n = n - 2;

%observation noise covariance matrix
Q = zeros(n);
Q1 = cofactor_matrix(elR, elM, snr_R, snr_M, pivot_index);
Q(1:n/2,1:n/2) = sigmaq_cod1 * Q1;
Q(n/2+1:end,n/2+1:end) = sigmaq_ph * Q1;

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x_hat = (N^-1)*A'*(Q^-1)*(y0-b);
XR = XR_approx + x_hat(1:3);

%estimated double difference ambiguities (without PIVOT)
N_hat_nopivot = x_hat(4:end);

%add a zero at PIVOT position
N_hat = zeros(n/2+1,1);
N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);

%estimation of the variance of the observation error
y_hat = A*x_hat + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat * (N^-1);

    %rover position covariance matrix
    cov_XR = Cxx(1:3,1:3);

    %combined ambiguity covariance matrix
    cov_N_nopivot = Cxx(4:end,4:end);

    %add one line and one column (zeros) at PIVOT position
    cov_N = zeros(n/2+1);
    cov_N(1:pivot_index-1,1:pivot_index-1)     = cov_N_nopivot(1:pivot_index-1,1:pivot_index-1);
    cov_N(pivot_index+1:end,pivot_index+1:end) = cov_N_nopivot(pivot_index:end,pivot_index:end);

else
    cov_XR = [];

    cov_N = [];
end

%DOP computation
if (nargout > 4)
    cov_XYZ = (A(1:n/2,1:3)'*A(1:n/2,1:3))^-1;
    cov_ENU = global2localCov(cov_XYZ, XR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end
