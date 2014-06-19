function [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch] = LS_SA_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold)

% SYNTAX:
%   [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = LS_SA_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys);
%
% INPUT:
%   XR_approx    = receiver approximate position (X,Y,Z)
%   XS           = satellite position (X,Y,Z)
%   pr_R         = code observations
%   snr_R        = signal-to-noise ratio
%   elR          = satellite elevation (vector)
%   distR_approx = approximate receiver-satellite distance (vector)
%   dtS          = satellite clock error (vector)
%   err_tropo_RS = tropospheric error
%   err_iono_RS  = ionospheric error
%   sys          = array with different values for different systems
%   SPP_threshold= maximum RMS of code single point positioning to accept current epoch
%
% OUTPUT:
%   XR   = estimated position (X,Y,Z)
%   dtR  = receiver clock error (scalar)
%   cov_XR  = estimated position error covariance matrix
%   var_dtR = estimated clock error variance
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%   cond_num = condition number on the eigenvalues of the N matrix
%   bad_obs = index of observations found as outlier
%   bad_epoch = 1:bad epoch, 0:good epoch
% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.2 beta
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

v_light = goGNSS.V_LIGHT;

%number of observations
n = length(pr_R);

%number of unknown parameters
m = 4;

if (~any(distR_approx))
    %approximate receiver-satellite distance
    XR_mat = XR_approx(:,ones(n,1))';
    distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
end

%design matrix
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
     (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
     (XR_approx(3) - XS(:,3)) ./ distR_approx, ... %column for Z coordinate
      ones(n,1)];        %column for receiver clock delay (multiplied by c)
  
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
b = distR_approx - v_light*dtS + err_tropo_RS + err_iono_RS;

%observation vector
y0 = pr_R;

%observation covariance matrix
Q = cofactor_matrix_SA(elR, snr_R);
invQ=diag((diag(Q).^-1));

%normal matrix
N = (A'*(invQ)*A);

if nargin<10 || (n == m) || exist('SPP_threshold','var')==0
    %least squares solution
    x   = (N^-1)*A'*(invQ)*(y0-b);
    %estimation of the variance of the observation error
    y_hat = A*x + b;
    v_hat = y0 - y_hat;
    sigma02_hat = (v_hat'*(invQ)*v_hat) / (n-m);
    
    if n==m
        bad_epoch=1;
    else
        bad_epoch=0;
    end
    bad_obs=[];
else    
    %------------------------------------------------------------------------------------
    % OUTLIER DETECTION (OPTIMIZED LEAVE ONE OUT)
    %------------------------------------------------------------------------------------
    search_for_outlier=1;
    bad_obs=[];
    
    while search_for_outlier==1
        [index_outlier,x,sigma02_hat]=OLOO(A, y0-b, Q);
        if index_outlier~=0
            bad_obs=[bad_obs;index_outlier];
            %fprintf('\nOUTLIER FOUND! obs %d/%d\n',index_outlier,length(y0));
            A(index_outlier,:)=[];
            y0(index_outlier,:)=[];
            b(index_outlier,:)=[];
            Q(index_outlier,:)=[];
            Q(:,index_outlier)=[];
            invQ=diag((diag(Q).^-1));
            [n,m]=size(A);
        else
            search_for_outlier=0;
        end
    end
    N = (A'*(invQ)*A);
    
    if n>m
        bad_epoch=(sigma02_hat>SPP_threshold^2);
    else
        bad_epoch=1;
    end
end

XR  = XR_approx + x(1:3);
dtR = x(4) / v_light;

%computation of the condition number on the eigenvalues of N
N_min_eig = min(eig(N));
N_max_eig = max(eig(N));
cond_num = N_max_eig / N_min_eig;

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat * (N^-1);
    cov_XR  = Cxx(1:3,1:3);
    var_dtR = Cxx(4,4) / v_light;
else
    cov_XR  = [];
    var_dtR = []; 
end

%DOP computation
if (nargout > 4)
    cov_XYZ = (A(:,1:3)'*A(:,1:3))^-1;
    cov_ENU = global2localCov(cov_XYZ, XR);

    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end
