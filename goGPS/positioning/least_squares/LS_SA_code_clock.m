function [dtR, ISBs, var_dtR, var_ISBs, bad_obs, bad_epoch, sigma02_hat, residuals_obs, y0, b, A, Q] = LS_SA_code_clock(pr_R, snr_R, elR, distR, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold, flag_residual_thres)

% SYNTAX:
%   [dtR, ISBs, var_dtR, var_ISBs, bad_obs, bad_epoch, sigma02_hat, residuals_obs, y0, b, A, Q] = LS_SA_code_clock(pr_R, snr_R, elR, distR, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold, flag_residual_thres);
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
%   ISBs = estimated inter-system biases
%   var_dtR = estimate error variance (scalar)
%   var_ISBs = variance of estimation errors (inter-system biases)
%   bad_obs = vector with ids of observations found as outlier
%   bad_epoch = 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   sigma02_hat = [a posteriori sigma (SPP sigma), v_hat'*(invQ)*v_hat), n-m]
%   residuals_obs = vector with residuals of all input observation, computed from the final estimates
%   y0 = least-squares observation vector
%   A = least-squares design matrix
%   b = least-squares known term vector
%   Q = observation covariance matrix
%
% DESCRIPTION:
%   Receiver clock error estimation by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.


%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     Stefano Caldera, ...
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
global flag_outlier flag_outlier_OLOO max_code_residual

v_light = goGNSS.V_LIGHT;
sigma02_hat = NaN(1,3);
residuals_obs = NaN(length(pr_R),1); %#ok<NASGU>
ISBs = [];
var_ISBs = [];

%number of observations
n = length(pr_R);

%number of unknown parameters
m = 1;

%design matrix
A = ones(n,1);

%if multi-system observations, then estimate an inter-system bias parameter for each additional system
uni_sys = unique(sys(sys ~= 0),'stable');
num_sys = length(uni_sys);
ISB = zeros(n,1);
if (num_sys > 1)
    if (any(floor(uni_sys) ~= 2)) %not only GLONASS
        ref_clock = 1;
    else %only GLONASS
        ref_clock = 0;
        A = [];
        m = m - 1;
    end
    m = m + num_sys - ref_clock;
    for s = (1 + ref_clock) : num_sys
        ISB(sys == uni_sys(s)) = 1;
        A = [A, ISB]; %#ok<AGROW>
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

bad_epoch=0;
bad_obs=[];
A0=A;
y00=y0-b;
index_obs = 1:length(y00); index_obs = index_obs';

% arg 10 => flag_residual_thres
if flag_outlier && nargin >= 10 && flag_residual_thres
    %remove observations with residuals exceeding thresholds
    x = (N^-1)*A'*(invQ)*(y0-b);
    y_hat = A*x + b;
    v_hat = y0 - y_hat;
    index_outlier = find(abs(v_hat) > max_code_residual);
    if (~isempty(index_outlier))
        if (n-length(index_outlier) < m)
            bad_epoch = 1;
        else
            bad_obs=unique([bad_obs;index_obs(index_outlier)]);
            A(index_outlier,:)=[];
            y0(index_outlier,:)=[];
            b(index_outlier,:)=[];
            Q(index_outlier,:)=[];
            Q(:,index_outlier)=[];
            invQ=diag((diag(Q).^-1));
            N = (A'*(invQ)*A);
            index_obs(index_outlier) = [];
            [n,m]=size(A);
        end
    end
end

% arg 9 => SPP_threshold
if nargin<9 || (n == m)
    %least squares solution
    x   = (N^-1)*A'*(invQ)*(y0-b);
    %estimation of the variance of the observation error
    y_hat = A*x + b;
    v_hat = y0 - y_hat;
    sigma02_hat(1,1) = (v_hat'*(invQ)*v_hat) / (n-m);
    sigma02_hat(1,2) = (v_hat'*(invQ)*v_hat);
    sigma02_hat(1,3) = n-m;
    residuals_obs = y00 - A0*x;
    if n==m
        bad_epoch=-1;
    end
else
    %------------------------------------------------------------------------------------
    % OUTLIER DETECTION (OPTIMIZED LEAVE ONE OUT)
    %------------------------------------------------------------------------------------
    search_for_outlier = flag_outlier_OLOO;

    if (num_sys > 1)
        sys0=sys;
        uni_sys0 = uni_sys;
    end
    if (search_for_outlier == 1)
        while search_for_outlier==1
            [index_outlier,x,sigma02_hat(1,1),v_hat]=OLOO(A, y0-b, Q);
            if index_outlier~=0
                bad_obs=[bad_obs;index_obs(index_outlier)]; %#ok<AGROW>
                if (num_sys > 1)
                    sys(index_outlier)=[];
                    uni_sys = unique(sys(sys ~= 0));
                    num_sys = length(uni_sys);
                    if length(uni_sys)<length(uni_sys0) % an i-s bias is not estimable anymore
                        A(:,1+find(uni_sys0==sys0(index_outlier))-1)=[];
                        A0(:,4+find(uni_sys0==sys0(index_outlier))-1)=[];
                        sys0=sys;
                        uni_sys0 = uni_sys;
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
    else
        x = (N^-1)*A'*(invQ)*(y0-b);
        y_hat = A*x + b;
        v_hat = y0 - y_hat;
    end
    N = (A'*(invQ)*A);
    residuals_obs = y00 - A0*x;
    sigma02_hat(1,2) = (v_hat'*(invQ)*v_hat);
    sigma02_hat(1,3) = n-m;
    sigma02_hat(1,1) = sigma02_hat(1,2)/sigma02_hat(1,3);
    if n>m
        if (flag_outlier && sigma02_hat(1)>SPP_threshold^2)
            bad_epoch=1;
        end
    else
        bad_epoch=-1;
    end
end

dtR = x(1) / v_light;

if (num_sys > 1)
    ISBs = x(ref_clock+[1:num_sys-ref_clock]) / v_light;
end

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat(1)*(N^-1);
    var_dtR = Cxx(1,1) / v_light^2;
    if (num_sys > 1)
        var_ISBs = Cxx(2:end,2:end) / v_light^2;
    end
else
    var_dtR = NaN;
    var_ISBs = NaN;
end

%prepare A and ISBs variables for output (only when GLONASS is enabled)
if (any(floor(uni_sys) == 2))
    pos = find(floor(uni_sys) == 2);
    GLO_slots = single((uni_sys(pos) - 2)*100 + 1);
    GLO_IFBs = zeros(14,1);
    GLO_IFBs(GLO_slots) = ISBs(pos-ref_clock);
    ISBs = [GLO_IFBs; ISBs(pos(end)+1:end)];
    GLO_rows = find(floor(sys) == 2);
    if (ref_clock ~= 0)
        A = [A(:,1) zeros(n,14) A(:,1+pos(end)+1:end)];
    else
        A = zeros(n,14);
    end
    for r = 1 : length(GLO_rows)
        A(GLO_rows(r),1+GLO_slots(r)) = 1;
    end
end
