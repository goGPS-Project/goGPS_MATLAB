function [dtR, var_dtR, bad_obs, bad_epoch, sigma02_hat, residuals_obs, is_bias] = LS_SA_code_clock(pr_R, snr_R, elR, distR, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold)

% SYNTAX:
%   [dtR, var_dtR] = LS_SA_code_clock(pr_R, snr_R, elR, distR, dtS, err_tropo_RS, err_iono_RS, sys);
%
% INPUT:
%   pr_R  = code observations (vector)
%   snr_R = signal-to-noise ratio (vector)
%   elR   = satellite elevation (vector)
%   distR = satellite-receiver distance (vector)
%   dtS   = satellite clock error (vector)
%   err_tropo_RS = tropospheric error (vector)
%   err_iono_RS  = ionospheric error (vector)
%   sys          = array with different values for different systems
%
% OUTPUT:
%   dtR = receiver clock error (scalar)
%   var_dtR = estimate error variance (scalar)
%   bad_obs = vector with ids of observations found as outlier
%   bad_epoch = 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   sigma02_hat = [a posteriori sigma (SPP sigma), v_hat'*(invQ)*v_hat), n-m] 
%   residuals_obs = vector with residuals of all input observation, computed from the final estimates
%   is_bias = inter-systems bias (vector with all possibile systems)
%
% DESCRIPTION:
%   Receiver clock error estimation by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
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
sigma02_hat = NaN(1,3);
residuals_obs = NaN(length(pr_R),1);
is_bias=NaN(6,1);

%number of observations
n = length(pr_R);

%number of unknown parameters
m = 1;

%design matrix
A = ones(n,1);

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
b = distR - v_light*dtS + err_tropo_RS + err_iono_RS;

%observation vector
y0 = pr_R;

%observation covariance matrix
Q = cofactor_matrix_SA(elR, snr_R);
invQ=diag((diag(Q).^-1));

%normal matrix
N = (A'*(invQ)*A);


if nargin<8 || (n == m) || exist('SPP_threshold','var')==0
    %least squares solution
    x   = (N^-1)*A'*(invQ)*(y0-b);
    %estimation of the variance of the observation error
    y_hat = A*x + b;
    v_hat = y0 - y_hat;
    sigma02_hat(1,1) = (v_hat'*(invQ)*v_hat) / (n-m);
    sigma02_hat(1,2) = (v_hat'*(invQ)*v_hat);
    sigma02_hat(1,3) = n-m;
    residuals_obs=v_hat;
    if (num_sys > 1)
        is_bias(uni_sys(2:end))=x(end-(num_sys-2):end);
    end
    
    if n==m
        bad_epoch=-1;
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
    index_obs = 1:length(y0);
    A0=A;
    y00=y0-b;
    if (num_sys > 1)
        sys0=sys;
        uni_sys0 = uni_sys;
        pos_isbias=m-(num_sys-2):m;
    end
    while search_for_outlier==1
        [index_outlier,x,sigma02_hat(1,1),v_hat]=OLOO(A, y0-b, Q);
        if index_outlier~=0
            bad_obs=[bad_obs;index_obs(index_outlier)];
            if (num_sys > 1)
                sys(index_outlier)=[];
                uni_sys = unique(sys(sys ~= 0));
                num_sys = length(uni_sys);
                if length(uni_sys)<length(uni_sys0) % an i-s bias is not estimable anymore
                    A(:,1+find(uni_sys0==sys0(index_outlier))-1)=[];
                end
            end
            %fprintf('\nOUTLIER FOUND! obs %d/%d\n',index_outlier,length(y0));
            A(index_outlier,:)=[];
            y0(index_outlier,:)=[];
            b(index_outlier,:)=[];
            Q(index_outlier,:)=[];
            Q(:,index_outlier)=[];
            invQ=diag((diag(Q).^-1));
            index_obs(index_outlier)=[];
            [n,m]=size(A);
        else
            search_for_outlier=0;
        end
    end
    if (num_sys > 1)
        is_bias(uni_sys(2:end))=x(end-(num_sys-2):end);
    end
    N = (A'*(invQ)*A);
    residuals_obs = y00 - A0*x;
    sigma02_hat(1,2) = (v_hat'*(invQ)*v_hat);
    sigma02_hat(1,3) = n-m;
    sigma02_hat(1,1) = sigma02_hat(1,2)/sigma02_hat(1,3);
    if n>m
        bad_epoch=(sigma02_hat(1)>SPP_threshold^2);
    else
        bad_epoch=-1;
    end
    
end

dtR = x(1) / v_light;

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat(1)*(N^-1);
    var_dtR = Cxx(1,1) / v_light^2;
else
    var_dtR = NaN; 
end
