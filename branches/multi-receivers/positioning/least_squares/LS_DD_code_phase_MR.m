function [Xb, N_hat, cov_XR, cov_N, PDOP, HDOP, VDOP, up_bound, lo_bound, posType, attitude, XR] = LS_DD_code_phase_MR ...
         (Xb_approx, XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, phase, attitude, geometry, flag_Tykhon)

% SYNTAX:
%   [XR, N_hat, cov_XR, cov_N, PDOP, HDOP, VDOP, up_bound, lo_bound, posType] = LS_DD_code_phase ...
%   (XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, phase, flag_Tykhon);
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
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%   attitude    = apriori attitude angles (raw, pitch, yaw)
%   geometry    = instrument geometry
%   flag_Tykhon = flag to enable/disable Tykhonov-Phillips regularization

%
% OUTPUT:
%   Xb = estimated position of barycenter (X,Y,Z)
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
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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
global lambda1 lambda2
global sigmaq_cod1 sigmaq_ph

if (phase == 1)
    lambda = lambda1;
else
    lambda = lambda2;
end

if (nargin < 21)
    flag_Tykhon = 0;
end


%number of observations
nRec=size(pr_R,3);
n = 2*size(pr_R,1)*(nRec);

%number of unknown parameters
m = 3 + (n/2/nRec) + 3; % must be corrected later for the pivot pivot ambiguity


%approximate receiver-satellite distance
distR_approx=NaN(n/2/(nRec),nRec);

XM_mat = XM(:,ones(n/2/(nRec),1))';
distM = sqrt(sum((XS-XM_mat).^2 ,2));
for i=1:nRec
    XR_approx_i=XR_approx(:,1,i);
    XR_mat = XR_approx_i(:,ones(n/2/(nRec),1))';
    distR_approx(:,i) = sqrt(sum((XS-XR_mat).^2 ,2));
end

%Xb_approx=mean([XR_approx(:,1,1), XR_approx(:,1,2),XR_approx(:,1,3)],2);
[phi_b_apriori, lam_b_apriori, h_b_apriori] = cart2geod(Xb_approx(1), Xb_approx(2), Xb_approx(3));


%% compute diff --------- must be put outside to avoid the recomputation every epoch
syms s_Xb s_Yb s_Zb s_phi_b s_lam_b s_roll s_pitch s_yaw s_x s_y s_z s_XS s_YS s_ZS s_XS_Piv s_YS_Piv s_ZS_Piv s_XM s_YM s_ZM
syms r11 r12 r13 r21 r22 r23 r31 r32 r33

% rotation from local to global
s_Rgl=[-sin(s_lam_b) cos(s_lam_b) 0; ...
    -sin(s_phi_b)*cos(s_lam_b) -sin(s_phi_b)*sin(s_lam_b) cos(s_phi_b); ...
    cos(s_phi_b)*cos(s_lam_b) cos(s_phi_b)*sin(s_lam_b) sin(s_phi_b)];


% rotation from instrumental to local
s_Rli=[cos(s_roll)*cos(s_pitch) -sin(s_pitch)*cos(s_yaw)+cos(s_roll)*sin(s_pitch)+sin(s_yaw) sin(s_roll)*sin(s_yaw)+cos(s_roll)*sin(s_pitch)*cos(s_yaw); ...
    sin(s_roll)+cos(s_pitch) sin(s_roll)*sin(s_pitch)*sin(s_yaw) sin(s_roll)*sin(s_pitch)*cos(s_yaw)-cos(s_roll)*sin(s_yaw); ...
    -sin(s_pitch) cos(s_pitch)*sin(s_yaw) cos(s_pitch)*cos(s_yaw)];

% rotation from instrumental to local
%s_Rli=[r11 r12 r13; r21 r22 r23; r31 r32 r33];


s_X=[s_Xb; s_Yb; s_Zb]+s_Rgl*s_Rli*[s_x;s_y;s_z];

s_PR_DD=sqrt((s_X(1)-s_XS)^2+(s_X(2)-s_YS)^2 + (s_X(3)-s_ZS)^2) - sqrt((s_XM-s_XS)^2+(s_YM-s_YS)^2 + (s_ZM-s_ZS)^2) - ...
    (sqrt((s_X(1)-s_XS_Piv)^2+(s_X(2)-s_YS_Piv)^2 + (s_X(3)-s_ZS_Piv)^2) - sqrt((s_XM-s_XS_Piv)^2+(s_YM-s_YS_Piv)^2 + (s_ZM-s_ZS_Piv)^2));


s_Ai=[diff(s_PR_DD,s_Xb) diff(s_PR_DD,s_Yb) diff(s_PR_DD,s_Zb) diff(s_PR_DD,s_roll) diff(s_PR_DD,s_pitch) diff(s_PR_DD,s_yaw)];
%s_Ai=[diff(s_PR_DD,s_Xb) diff(s_PR_DD,s_Yb) diff(s_PR_DD,s_Zb) diff(s_PR_DD,r11) diff(s_PR_DD,r12) diff(s_PR_DD,r13) diff(s_PR_DD,r21) diff(s_PR_DD,r22) diff(s_PR_DD,r23) diff(s_PR_DD,r31) diff(s_PR_DD,r32) diff(s_PR_DD,r33)];




F_Ai=inline(s_Ai); %% <-- colonne della matrice disegno corrispondenti alle incognite geometriche
F_PR_DD=inline(s_PR_DD); %% <-- parte geometrica dei termini noti
F_s_X=inline(s_X); %% <-- per calcolare le coordinate aggiornate di ogni ricevitore a partire da quelle nuove del punto 'baricentro'
%% 




A=zeros(n,m);
b=[];
y0=[];
%design matrix , known term vector, observation vector (code)
for i=1:nRec
    for j=1:size(pr_R,1)
        A(j+(i-1)*size(pr_R,1),1:6)=F_Ai(XS(j,1), XS(pivot_index,1), Xb_approx(1), XS(j,2), XS(pivot_index,2), Xb_approx(2), XS(j,3), XS(pivot_index,3), Xb_approx(3), lam_b_apriori, phi_b_apriori, attitude(2), attitude(1), geometry(1,i), geometry(2,i), attitude(3), geometry(3,i));
        b=[b;F_PR_DD(XM(1),XS(j,1), XS(pivot_index,1), Xb_approx(1), XM(2), XS(j,2), XS(pivot_index,2), Xb_approx(2), XM(3), XS(j,3), XS(pivot_index,3), Xb_approx(3), lam_b_apriori, phi_b_apriori, attitude(2), attitude(1), geometry(1,i), geometry(2,i), attitude(3), geometry(3,i)) + ...
            (err_tropo_R(j,i) - err_tropo_M(j,1)) - (err_tropo_R(pivot_index,i)  - err_tropo_M(pivot_index,1)) + ...
            (err_iono_R(j,i)  - err_iono_M(j,1))  - (err_iono_R(pivot_index,i)   - err_iono_M(pivot_index,1))];
        y0=[y0; (pr_R(j,1,i) - pr_M(j,1)) - (pr_R(pivot_index,1,i) - pr_M(pivot_index,1))];

    
    end
end

%design matrix , known term vector, observation vector (phase)
for i=1:nRec
    for j=1:size(pr_R,1)
        A(nRec*size(pr_R,1)+j+(i-1)*size(pr_R,1),1:6)=F_Ai(XS(j,1), XS(pivot_index,1), Xb_approx(1), XS(j,2), XS(pivot_index,2), Xb_approx(2), XS(j,3), XS(pivot_index,3), Xb_approx(3), lam_b_apriori, phi_b_apriori, attitude(2), attitude(1), geometry(1,i), geometry(2,i), attitude(3), geometry(3,i));
        b=[b;F_PR_DD(XM(1),XS(j,1), XS(pivot_index,1), Xb_approx(1), XM(2), XS(j,2), XS(pivot_index,2), Xb_approx(2), XM(3), XS(j,3), XS(pivot_index,3), Xb_approx(3), lam_b_apriori, phi_b_apriori, attitude(2), attitude(1), geometry(1,i), geometry(2,i), attitude(3), geometry(3,i)) + ...
            (err_tropo_R(j,i) - err_tropo_M(j,1)) - (err_tropo_R(pivot_index,i)  - err_tropo_M(pivot_index,1)) - ...
            (err_iono_R(j,i)  - err_iono_M(j,1))  + (err_iono_R(pivot_index,i)   - err_iono_M(pivot_index,1))];
        y0=[y0; lambda*((ph_R(j,1,i) - ph_M(j,1)) - (ph_R(pivot_index,1,i) - ph_M(pivot_index,1)))];
        A(nRec*size(pr_R,1)+j+(i-1)*size(pr_R,1),6+j)=-lambda;
    
    end
end


%remove pivot-pivot lines
A(pivot_index:size(pr_R,1):end, :) = [];
A(            :, pivot_index+6)       = [];
b(pivot_index:size(pr_R,1):end)    = [];
y0(pivot_index:size(pr_R,1):end)    = [];
[n,m]=size(A);


%observation noise covariance matrix
Q = zeros(n);
for i=1:nRec    
    Q1 = cofactor_matrix(elR(:,i), elM, snr_R(:,1,i), snr_M, pivot_index);
    Q((i-1)*(size(pr_R,1)-1)+1 : (i-1)*(size(pr_R,1)-1)+ (size(pr_R,1)-1), (i-1)*(size(pr_R,1)-1)+1 : (i-1)*(size(pr_R,1)-1)+ (size(pr_R,1)-1)) = sigmaq_cod1 * Q1;
    Q(n/2+(i-1)*(size(pr_R,1)-1)+1: n/2+(i-1)*(size(pr_R,1)-1)+size(pr_R,1)-1,n/2+(i-1)*(size(pr_R,1)-1)+1:n/2+(i-1)*(size(pr_R,1)-1)+size(pr_R,1)-1) = sigmaq_ph * Q1;
end

Q=eye(n);
      
%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x_hat = (N^-1)*A'*(Q^-1)*(y0-b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Xb=Xb_approx+x_hat(1:3);
[phi_b, lam_b, h_b] = cart2geod(Xb(1), Xb(2), Xb(3));
attitude= attitude+x_hat(4:6);

for i=1:nRec    
    XR(1:3,1,i)=F_s_X(Xb(1),Xb(2),Xb(3),lam_b,phi_b,attitude(2),attitude(1),geometry(1,i),geometry(2,i),attitude(3),geometry(3,i));
end





N_hat=x_hat(7:end);
cov_XR=[];
cov_N=[];
PDOP=[];
HDOP=[];
VDOP=[];
up_bound=[];
lo_bound=[];
posType=[];













% 
% 
% %design matrix (code)
% A = [((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index,1)) / distR_approx(pivot_index)), ... %column for X coordinate
%      ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index,2)) / distR_approx(pivot_index)), ... %column for Y coordinate
%      ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index,3)) / distR_approx(pivot_index)), ... %column for Z coordinate
%      zeros(n/2,n/2)]; %column for phase ambiguities   (here zero)
% 
% %design matrix (phase)
% A = [A; ((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index,1)) / distR_approx(pivot_index)), ... %column for X coordinate
%         ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index,2)) / distR_approx(pivot_index)), ... %column for Y coordinate
%         ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index,3)) / distR_approx(pivot_index)), ... %column for Z coordinate
%         -lambda * eye(n/2)]; %column for phase ambiguities
% 
% %known term vector
% b_pr =        (distR_approx - distM)      - (distR_approx(pivot_index) - distM(pivot_index));       %approximate pseudorange DD
% b_pr = b_pr + (err_tropo_R - err_tropo_M) - (err_tropo_R(pivot_index)  - err_tropo_M(pivot_index)); %tropospheric error DD
% b_pr = b_pr + (err_iono_R  - err_iono_M)  - (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD (code)
% b_ph = b_pr - (err_iono_R  - err_iono_M)  - (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD (phase)
% b = [b_pr; b_ph];
% 
% %observation vector
% y0_pr =         (pr_R - pr_M) - (pr_R(pivot_index) - pr_M(pivot_index));
% y0_ph = lambda*((ph_R - ph_M) - (ph_R(pivot_index) - ph_M(pivot_index)));
% y0 = [y0_pr; y0_ph];
% 
% %remove pivot-pivot lines
% A( [pivot_index, pivot_index+n/2], :) = [];
% A(            :, pivot_index+3)       = [];
% b( [pivot_index, pivot_index+n/2])    = [];
% y0([pivot_index, pivot_index+n/2])    = [];
% n = n - 2;
% 
% %observation noise covariance matrix
% Q = zeros(n);
% Q1 = cofactor_matrix(elR, elM, snr_R, snr_M, pivot_index);
% Q(1:n/2,1:n/2) = sigmaq_cod1 * Q1;
% Q(n/2+1:end,n/2+1:end) = sigmaq_ph * Q1;
% 
% %normal matrix
% N = (A'*(Q^-1)*A);
% 
% %least squares solution
% x_hat = (N^-1)*A'*(Q^-1)*(y0-b);

if (flag_Tykhon)
    %Processing with Tykhonov-Phillips regularization
    [x_hat, bias, lbd] = tykhonov_regularization(x_hat, y0, b, A, Q);
    
    %computation of the condition number on the eigenvalues of N
    N_min_eig = min(eig(N + lbd*eye(m)));
    N_max_eig = max(eig(N + lbd*eye(m)));
    cond_num = ceil(log10(N_max_eig / N_min_eig));
else
    %computation of the condition number on the eigenvalues of N
    N_min_eig = min(eig(N));
    N_max_eig = max(eig(N));
    cond_num = ceil(log10(N_max_eig / N_min_eig));
end

%covariance matrix of the estimation error
% if (n > m)
%     if (~flag_Tykhon)
%         Qxx = (N^-1);   %vcm of parameter computed from inverse of normal matrix
%     else
%         Qxx = ((N + lbd*eye(m))^-1)*N*((N + lbd*eye(m))^-1); %vcm of parameter
%     end
%     [U] = chol(Qxx);    %compute cholesky decomp. for identical comp of vcm purpose
%     Cxx = U'*U; %find back the vcm of parameter, now the off diag. comp. are identical :)
%     
%     %%this part for future statistical test of data reliability :)
%     %     Q_hat = A*Cxx*A';
%     %     Qv_hat= Q - Q_hat;
%     %     alpha = 0.05;
%     %     [identify,BNR_xhat, BNR_y] = test_statistic(Q_hat, Q, A, alpha, v_hat,Qv_hat)
%     
%     %rover position covariance matrix
%     %     cov_XR = Cxx(1:3,1:3);
%     Cbb = Cxx(1:3,1:3); %vcm block of baseline component
%     Cba = Cxx(1:3,4:end);% correlation of baselines and ambiguities comp.
%     Cab = Cxx(4:end,1:3);
%     Caa = Cxx(4:end,4:end);%vcm block of ambiguity component 
%     
%     %ambiguity estimation
%     afl = x_hat(4:end); %extract float ambiguity from float least-squares
%     [afix, sqnorm, Qahat, Z, D, L] = lambda_routine2(afl,Caa); %using LAMBDA routines
%     
%     %Integer ambiguity validation
%     %----------------------------
%     O1 = (afl - afix(:,1))' * (Caa)^-1 * (afl - afix(:,1)); %nominator component
%     O2 = (afl - afix(:,2))' * (Caa)^-1 * (afl - afix(:,2)); %denominator component
%     O  = O1 / O2; %ratio test
%     
%     if O <= 0.5%0.787%0.6    %if below threshold, use the LAMBDA result
%         z = afix(:,1);
%         posType = 1; %fixed
%     else                     %if above threshold, reject the LAMBDA result
%         z = afl;
%         posType = 0; %float
%     end
%     
%     if (~flag_Tykhon)
%         bias = check_bias(D); %no bias detected
%     else
%         bias = check_bias(D,bias(4:end)); %use bias component from ambiguity
%     end
%     [up_bound, lo_bound] = success_rate(D,L,bias); %compute success rate of ambiguity
%     
%     %Computing definite coordinate
%     %-----------------------------
%     bfl    = x_hat(1:3); %float baseline component
%     bdef   = bfl - Cba * Caa^-1 * (afl - z); %compute the conditional baseline component
%     cov_XR = Cbb - (Cba * Caa^-1 * Cab) + (Cba * Caa^-1*Qahat*(Caa^-1)'*Cab); %compute its vcm
%     XR     = (XR_approx + bdef); %define the definitive coordinate
%     %std    = sqrt(diag(cov_XR)); %standard deviation of baseline
%     %bl     = norm(XM - XR); %baseline length
%     
%     %estimated double difference ambiguities (without PIVOT)
%     N_hat_nopivot = z;
%     
%     %add a zero at PIVOT position
%     N_hat = zeros(n/2+1,1);
%     N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
%     N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
%     
%     %combined ambiguity covariance matrix
%     cov_N_nopivot = Qahat;
%     
%     %add one line and one column (zeros) at PIVOT position
%     cov_N = zeros(n/2+1);
%     cov_N(1:pivot_index-1,1:pivot_index-1)     = cov_N_nopivot(1:pivot_index-1,1:pivot_index-1);
%     cov_N(pivot_index+1:end,pivot_index+1:end) = cov_N_nopivot(pivot_index:end,pivot_index:end);
%     
% else
%     cov_XR = [];
%     
%     cov_N = [];
% end
% 
% %DOP computation
% if (nargout > 4)
%     cov_XYZ = (A(1:n/2,1:3)'*A(1:n/2,1:3))^-1;
%     cov_ENU = global2localCov(cov_XYZ, XR);
%     
%     PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
%     HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
%     VDOP = sqrt(cov_ENU(3,3));
% end
