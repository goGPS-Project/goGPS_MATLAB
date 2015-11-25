function [xR, cov_XR] = LS_DD_code(XR_approx, XS, pr_R, pr_M, snr_R, snr_M, elR, elM, distR_approx, distM, err_tropo_R, err_tropo_M, err_iono_R, err_iono_M, pivot_index)

% SYNTAX:
%   [xR, cov_XR] = LS_DD_code(XR_approx, XS, pr_R, pr_M, snr_R, snr_M, elR, elM, distR_approx, distM, err_tropo_R, err_tropo_M, err_iono_R, err_iono_M, pivot_index);
%
% INPUT:
%   XR_approx    = receiver approximate position (X,Y,Z)
%   XS           = satellite position (X,Y,Z)
%   pr_R         = receiver code observations
%   pr_M         = master code observations
%   snr_R        = receiversignal-to-noise ratio
%   snr_M        = mastersignal-to-noise ratio
%   elR          = satellite elevation (vector)
%   elM          = satellite elevation (vector)
%   distR_approx = approximate receiver-satellite distance (vector)
%   distM        = master station-satellite distance (vector)
%   err_tropo_R  = tropospheric error
%   err_tropo_M  = tropospheric error
%   err_iono_R   = ionospheric error
%   err_iono_M   = ionospheric error
%   pivot_index  = index identifying the pivot satellite
%
% OUTPUT:
%   xR = estimated position (X,Y,Z)
%   cov_XR = estimated position error covariance matrix
%
% DESCRIPTION:
%   Least squares solution using code double differences.
%   Epoch-by-epoch solution.

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

%number of observations
n = length(pr_R);

%number of unknown parameters
m = 3;

% %approximate receiver-satellite distance
% XR_mat = XR_approx(:,ones(n,1))';
% XM_mat = XM(:,ones(n,1))';
% distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
% distM = sqrt(sum((XS-XM_mat).^2 ,2));

%design matrix
A = [((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index,1)) / distR_approx(pivot_index)), ... %column for X coordinate
     ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index,2)) / distR_approx(pivot_index)), ... %column for Y coordinate
     ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index,3)) / distR_approx(pivot_index))];    %column for Z coordinate

%known term vector
b =     (distR_approx - distM)      - (distR_approx(pivot_index) - distM(pivot_index));       %approximate pseudorange DD
b = b + (err_tropo_R - err_tropo_M) - (err_tropo_R(pivot_index)  - err_tropo_M(pivot_index)); %tropospheric error DD
b = b + (err_iono_R  - err_iono_M)  - (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD

%observation vector
y0 = (pr_R - pr_M) - (pr_R(pivot_index) - pr_M(pivot_index));

%remove pivot-pivot lines
A(pivot_index, :) = [];
b(pivot_index)    = [];
y0(pivot_index)   = [];
n = n - 1;

%observation covariance matrix
Q = cofactor_matrix(elR, elM, snr_R, snr_M, pivot_index);

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x  = (N^-1)*A'*(Q^-1)*(y0-b);
xR = XR_approx + x;

%estimation of the variance of the observation error
y_hat = A*x + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    cov_XR = sigma02_hat * (N^-1);
else
    cov_XR = [];
end
