function [XR, dtR, ISBs, cov_XR, var_dtR, var_ISBs, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, sigma02_hat, residuals_obs, y0, b, A, Q] = LS_SA_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold, flag_residual_thres)

% SYNTAX:
%   [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, sigma02_hat, residuals_obs, y0, b, A, Q] = LS_SA_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold, flag_residual_thres);
%
% INPUT:
%   XR_approx     = receiver approximate position (X,Y,Z)
%   XS            = satellite position (X,Y,Z)
%   pr_R          = code observations
%   snr_R         = signal-to-noise ratio
%   elR           = satellite elevation (vector)
%   distR_approx  = approximate receiver-satellite distance (vector)
%   dtS           = satellite clock error (vector)
%   err_tropo_RS  = tropospheric error
%   err_iono_RS   = ionospheric error
%   sys           = array with different values for different systems
%   SPP_threshold = maximum RMS of code single point positioning to accept current epoch
%
% OUTPUT:
%   XR   = estimated position (X,Y,Z)
%   dtR  = receiver clock error (scalar)
%   ISBs = estimated inter-system biases
%   cov_XR  = estimated position error covariance matrix
%   var_dtR = estimated clock error variance
%   var_ISBs = variance of estimation errors (inter-system biases)
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%   cond_num = condition number on the eigenvalues of the N matrix
%   bad_obs = vector with ids of observations found as outlier
%   bad_epoch = 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   sigma02_hat = [a posteriori sigma (SPP sigma), v_hat'*(invQ)*v_hat), n-m]
%   residuals_obs = vector with residuals of all input observation, computed from the final estimates
%   y0 = least-squares observation vector
%   A = least-squares design matrix
%   b = least-squares known term vector
%   Q = observation covariance matrix

% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 3
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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

global flag_outlier max_code_residual

v_light = goGNSS.V_LIGHT;
sigma02_hat = NaN(1,3);
residuals_obs = NaN(length(pr_R),1); %#ok<NASGU>
ISBs = [];
var_ISBs = [];

%number of observations
n = length(pr_R);

%number of unknown parameters
m = 4;

if (~any(distR_approx))
    if (any(XR_approx))
        %receiver-satellite distance, corrected by the Shapiro delay
        [~, distR_approx] = relativistic_range_error_correction(XR_approx, XS);
    else
        %approximate receiver-satellite distance
        XR_mat = XR_approx(:,ones(n,1))';
        distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
    end
end

%design matrix
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
     (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
     (XR_approx(3) - XS(:,3)) ./ distR_approx, ... %column for Z coordinate
      ones(n,1)];        %column for receiver clock delay (multiplied by c)

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
b = distR_approx - v_light*dtS + err_tropo_RS + err_iono_RS;

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

if (flag_outlier && exist('flag_residual_thres','var') && flag_residual_thres == 1)
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

if nargin<10 || (n == m) || exist('SPP_threshold','var')==0
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
    search_for_outlier = flag_outlier;

    if (num_sys > 1)
        sys0=sys;
        uni_sys0 = uni_sys;
    end
    if (search_for_outlier == 1)
        while search_for_outlier==1
            [index_outlier,x,sigma02_hat(1,1),v_hat]=OLOO(A, y0-b, Q);
            if index_outlier~=0
                bad_obs=unique([bad_obs;index_obs(index_outlier)]);
                if (num_sys > 1)
                    sys(index_outlier)=[];
                    uni_sys = unique(sys(sys ~= 0));
                    num_sys = length(uni_sys);
                    if length(uni_sys)<length(uni_sys0) % an i-s bias is not estimable anymore
                        A(:,4+find(uni_sys0==sys0(index_outlier))-1)=[];
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
        if (sigma02_hat(1)>SPP_threshold^2)
            bad_epoch=1;
        end
    else
        bad_epoch=-1;
    end
end

XR  = XR_approx + x(1:3);
dtR = x(4) / v_light;

%estimated inter-system biases
if (num_sys > 1)
    ISBs = x(3+ref_clock+[1:num_sys-ref_clock]) / v_light;
end

%computation of the condition number on the eigenvalues of N
N_min_eig = min(eig(N));
N_max_eig = max(eig(N));
cond_num = N_max_eig / N_min_eig;

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat(1) * (N^-1);
    cov_XR  = Cxx(1:3,1:3);
    var_dtR = Cxx(4,4) / v_light^2;
    if (num_sys > 1)
        var_ISBs = Cxx(5:end,5:end) / v_light^2;
    end
else
    cov_XR  = [];
    var_dtR = [];
    var_ISBs = [];
end

%DOP computation
if (nargout > 4)
    cov_XYZ = (A'*A)^-1;
    cov_XYZ = cov_XYZ(1:3,1:3);
    cov_ENU = global2localCov(cov_XYZ, XR);

    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
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
        A = [A(:,1:3) zeros(n,14) A(:,4+pos(end)+1:end)];
    else
        A = zeros(n,14);
    end
    for r = 1 : length(GLO_rows)
        A(GLO_rows(r),4+GLO_slots(r)) = 1;
    end
end
