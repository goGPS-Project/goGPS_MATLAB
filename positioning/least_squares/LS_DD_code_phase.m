function [XR, N_hat, cov_XR, cov_N, PDOP, HDOP, VDOP] = LS_DD_code_phase ...
         (XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda, flag_LAMBDA, flag_Tykhon, flag_iter, flag_iter_zero)

% SYNTAX:
%   [XR, N_hat, cov_XR, cov_N, PDOP, HDOP, VDOP] = LS_DD_code_phase ...
%   (XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda, flag_LAMBDA, flag_Tykhon, flag_iter, flag_iter_zero);
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
%   pivot_index = index identifying the pivot satellite
%   lambda      = vector containing GNSS wavelengths for available satellites
%   flag_LAMBDA    = flag to enable/disable ambiguity resolution by LAMBDA 
%   flag_Tykhon    = flag to enable/disable Tykhonov-Phillips regularization (experimental)
%   flag_iter      = flag to enable/disable iterative least squares (experimental)
%   flag_iter_zero = flag to enable/disable initialization to zero for iterative least squares (experimental)
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
%   Epoch-by-epoch solution, optionally with Tykhonov-Phillips regularization.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Hendy F. Suhandri
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

if (nargin < 19)
    flag_Tykhon = 0;    %<--- set 1 to force Tykhonov (experimental)
end

if (nargin < 20)
    flag_iter = 0;      %<--- set 1 to force iterative least squares (experimental)
end

if (nargin < 21)
    flag_iter_zero = 0; %<--- set 1 to force zero-initialization for iterative least squares (experimental)
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
        diag(-lambda) .* eye(n/2)]; %column for phase ambiguities

%known term vector
b    =     (distR_approx - distM)      - (distR_approx(pivot_index) - distM(pivot_index));       %approximate pseudorange DD
b    = b + (err_tropo_R - err_tropo_M) - (err_tropo_R(pivot_index)  - err_tropo_M(pivot_index)); %tropospheric error DD
b_pr = b + (err_iono_R  - err_iono_M)  - (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD (code)
b_ph = b - (err_iono_R  - err_iono_M)  + (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD (phase)
b = [b_pr; b_ph];

%observation vector
y0_pr =         (pr_R - pr_M) - (pr_R(pivot_index) - pr_M(pivot_index));
y0_ph = lambda.*((ph_R - ph_M) - (ph_R(pivot_index) - ph_M(pivot_index)));
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

%estimation of the variance of the observation error
y_hat = A*x_hat + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

if (~flag_LAMBDA)
    %apply least squares solution
    XR = XR_approx + x_hat(1:3);
    
    %estimated double difference ambiguities (without PIVOT)
    N_hat_nopivot = x_hat(4:end);
    
    %add a zero at PIVOT position
    N_hat = zeros(n/2+1,1);
    N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
    N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
    
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
    
else %apply LAMBDA
    
    if (flag_Tykhon)
        %Processing with Tykhonov-Phillips regularization
        [x_hat, bias, lbd] = tykhonov_regularization(x_hat, y0, b, A, Q); %#ok<ASGLU>
        
        %computation of the condition number on the eigenvalues of N
        N_min_eig = min(eig(N + lbd*eye(m)));
        N_max_eig = max(eig(N + lbd*eye(m)));
        cond_num = ceil(log10(N_max_eig / N_min_eig)); %#ok<NASGU>
    else
        %computation of the condition number on the eigenvalues of N
        N_min_eig = min(eig(N));
        N_max_eig = max(eig(N));
        cond_num = ceil(log10(N_max_eig / N_min_eig)); %#ok<NASGU>
    end
    
    if (~flag_Tykhon)
        Qxx = (N^-1); %vcm of parameter computed from inverse of normal matrix
    else
        Qxx = ((N + lbd*eye(m))^-1)*N*((N + lbd*eye(m))^-1); %vcm of parameter
    end
    
    if (flag_iter)
        
        if (flag_iter_zero)
            %force initial value to zeros (otherwise use x_hat and Qxx from previous steps)
            x_hat = zeros(m,1);
            Qxx   = eye(m);
        end
        
        %iterative least squares (float) solution
        M = ((A'*Q^-1*A) + eye(m))^-1; %normal matrix for iterative LS
        U = A'*Q^-1*(y0-b); %part A'*P*Y of normal equation
        dummy = zeros(m);
        x_n = ones(m,1);
        tol = 1e-4; %threshold of iterative process
        cnt = 1;
        
        while (norm(x_n) > tol)
            x_old = x_hat;
            dummy = dummy + M^(cnt); %compute iterative normal matrix
            x_hat = (dummy)*U + M^(cnt)*x_hat; %estimate parameter
            Qxx   = (dummy)*A'*Q^-1*A*(dummy)' + M^(cnt)*Qxx*M^(cnt)'; %estimate vcm of parameter
            
            x_n = x_hat - x_old;
            cnt = cnt + 1;
        end
    end
    
    if (n > m && sigmaq_ph ~= 1e30)
        
        %covariance matrix
        Cxx = sigma02_hat * (N^-1);
        
        [U] = chol(Cxx); %compute cholesky decomp. for identical comp of vcm purpose
        Cxx = U'*U; %find back the vcm of parameter, now the off diag. comp. are identical :)
        
        %%this part for future statistical test of data reliability :)
        %     Q_hat = A*Cxx*A';
        %     Qv_hat= Q - Q_hat;
        %     alpha = 0.05;
        %     [identify,BNR_xhat, BNR_y] = test_statistic(Q_hat, Q, A, alpha, v_hat,Qv_hat)

        Cbb = Cxx(1:3,1:3); %vcm block of baseline component
        Cba = Cxx(1:3,4:end);% correlation of baselines and ambiguities comp.
        Cab = Cxx(4:end,1:3);
        Caa = Cxx(4:end,4:end);%vcm block of ambiguity component
        
        bhat = x_hat(1:3);
        ahat = x_hat(4:end);
        
        [bcheck, acheck, Qzhat] = lambdafix(bhat, ahat, Cbb, Caa, Cba);
        
        XR     = (XR_approx + bcheck);
        cov_XR = Cbb - (Cba * Caa^-1 * Cab) + (Cba * Caa^-1*Qzhat*(Caa^-1)'*Cab);
        
        %ambiguity estimation
        %afl = x_hat(4:end); %extract float ambiguity from float least-squares
        %[afix, sqnorm, Qahat, Z, D, L] = lambda_routine2(afl,Caa); %using LAMBDA routines

        %Integer ambiguity validation
        %----------------------------
        %O1 = (afl - afix(:,1))' * (Caa)^-1 * (afl - afix(:,1)); %nominator component
        %O2 = (afl - afix(:,2))' * (Caa)^-1 * (afl - afix(:,2)); %denominator component
        %O  = O1 / O2; %ratio test
        
        %if O <= 0.5%0.787%0.6    %if below threshold, use the LAMBDA result
        %    z = afix(:,1);
        %    posType = 1; %fixed
        %else                     %if above threshold, reject the LAMBDA result
        %    z = afl;
        %    posType = 0; %float
        %end
        
        %if (~flag_Tykhon)
        %    bias = zeros(size(afix,1),1); %no bias
        %else
        %    bias = bias(4:end); %use bias component from ambiguity
        %end
        %[up_bound, lo_bound] = success_rate(D,L,bias); %compute success rate of ambiguity
        
        %Computing definite coordinate
        %-----------------------------
        %bfl    = x_hat(1:3); %float baseline component
        %bdef   = bfl - Cba * Caa^-1 * (afl - z); %compute the conditional baseline component
        %cov_XR = Cbb - (Cba * Caa^-1 * Cab) + (Cba * Caa^-1*Qahat*(Caa^-1)'*Cab); %compute its vcm
        %XR     = (XR_approx + bdef); %define the definitive coordinate
        %std    = sqrt(diag(cov_XR)); %standard deviation of baseline
        %bl     = norm(XM - XR); %baseline length
        
        %estimated double difference ambiguities (without PIVOT)
        N_hat_nopivot = acheck;
        
        %add a zero at PIVOT position
        N_hat = zeros(n/2+1,1);
        N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
        N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
        
        %combined ambiguity covariance matrix
        cov_N_nopivot = Qzhat;
        
        %add one line and one column (zeros) at PIVOT position
        cov_N = zeros(n/2+1);
        cov_N(1:pivot_index-1,1:pivot_index-1)     = cov_N_nopivot(1:pivot_index-1,1:pivot_index-1);
        cov_N(pivot_index+1:end,pivot_index+1:end) = cov_N_nopivot(pivot_index:end,pivot_index:end);
        
    else
        XR = XR_approx + x_hat(1:3);
        N_hat_nopivot = x_hat(4:end);
        
        %add a zero at PIVOT position
        N_hat = zeros(n/2+1,1);
        N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
        N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
        
        cov_XR = [];
        cov_N = [];
    end
end

%DOP computation
if (nargout > 4)
    cov_XYZ = (A(1:n/2,1:3)'*A(1:n/2,1:3))^-1;
    cov_ENU = global2localCov(cov_XYZ, XR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end
