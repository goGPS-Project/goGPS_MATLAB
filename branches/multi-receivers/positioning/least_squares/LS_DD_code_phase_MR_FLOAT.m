function [Xb, N_hat, cov_Xb, cov_N, cov_ATT, attitude, XR, PDOP, HDOP, VDOP] = LS_DD_code_phase_MR_FLOAT ...
    (Xb_approx, XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, phase, attitude, geometry, flag_Tykhon,F_Ai, F_PR_DD, F_s_X,threshold)

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
global lambda1 lambda2 O
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
nRec=size(pr_R,2);
n = 2*size(pr_R,1)*(nRec);

%number of unknown parameters
m = 3 + (n/2) + 3; % must be corrected later for the pivot pivot ambiguity


%approximate receiver-satellite distance
distR_approx=NaN(n/2/(nRec),nRec);

XM_mat = XM(:,ones(n/2/(nRec),1))';
distM = sqrt(sum((XS-XM_mat).^2 ,2));
for i=1:nRec
    XR_approx_i=XR_approx(:,i);
    XR_mat = XR_approx_i(:,ones(n/2/(nRec),1))';
    distR_approx(:,i) = sqrt(sum((XS-XR_mat).^2 ,2));
end

%Xb_approx=mean([XR_approx(:,1,1), XR_approx(:,1,2),XR_approx(:,1,3)],2);
[phi_b_apriori, lam_b_apriori, h_b_apriori] = cart2geod(Xb_approx(1), Xb_approx(2), Xb_approx(3));

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
        y0=[y0; (pr_R(j,i) - pr_M(j,1)) - (pr_R(pivot_index,i) - pr_M(pivot_index,1))];
        
        
    end
end

%design matrix , known term vector, observation vector (phase)
for i=1:nRec
    for j=1:size(pr_R,1)
        A(nRec*size(pr_R,1)+j+(i-1)*size(pr_R,1),1:6)=F_Ai(XS(j,1), XS(pivot_index,1), Xb_approx(1), XS(j,2), XS(pivot_index,2), Xb_approx(2), XS(j,3), XS(pivot_index,3), Xb_approx(3), lam_b_apriori, phi_b_apriori, attitude(2), attitude(1), geometry(1,i), geometry(2,i), attitude(3), geometry(3,i));
        b=[b;F_PR_DD(XM(1),XS(j,1), XS(pivot_index,1), Xb_approx(1), XM(2), XS(j,2), XS(pivot_index,2), Xb_approx(2), XM(3), XS(j,3), XS(pivot_index,3), Xb_approx(3), lam_b_apriori, phi_b_apriori, attitude(2), attitude(1), geometry(1,i), geometry(2,i), attitude(3), geometry(3,i)) + ...
            (err_tropo_R(j,i) - err_tropo_M(j,1)) - (err_tropo_R(pivot_index,i)  - err_tropo_M(pivot_index,1)) - ...
            (err_iono_R(j,i)  - err_iono_M(j,1))  + (err_iono_R(pivot_index,i)   - err_iono_M(pivot_index,1))];
        y0=[y0; lambda*((ph_R(j,i) - ph_M(j,1)) - (ph_R(pivot_index,i) - ph_M(pivot_index,1)))];
        A(nRec*size(pr_R,1)+j+(i-1)*size(pr_R,1),6+(i-1)*size(pr_R,1)+j)=-lambda;
        
    end
end


%remove pivot-pivot lines
A(pivot_index:size(pr_R,1):end, :) = [];
A(            :, 6+pivot_index:size(pr_R,1):end)= [];
b(pivot_index:size(pr_R,1):end)    = [];
y0(pivot_index:size(pr_R,1):end)    = [];
[n,m]=size(A);

%observation noise covariance matrix
Q = zeros(n);
for i=1:nRec
    Q1 = cofactor_matrix(elR(:,i), elM, snr_R(:,i), snr_M, pivot_index);
    Q((i-1)*(size(pr_R,1)-1)+1 : (i-1)*(size(pr_R,1)-1)+ (size(pr_R,1)-1), (i-1)*(size(pr_R,1)-1)+1 : (i-1)*(size(pr_R,1)-1)+ (size(pr_R,1)-1)) = sigmaq_cod1 * Q1;
    Q(n/2+(i-1)*(size(pr_R,1)-1)+1: n/2+(i-1)*(size(pr_R,1)-1)+size(pr_R,1)-1,n/2+(i-1)*(size(pr_R,1)-1)+1:n/2+(i-1)*(size(pr_R,1)-1)+size(pr_R,1)-1) = sigmaq_ph * Q1;
end

% % add further constraints: I assume roll and pitch to be near zero!
% 
% y0=[y0;-attitude(1);-attitude(2)];
% A(n+1,4)=1;
% A(n+2,5)=1;
% b(n+1)=0;
% b(n+2)=0;
% 
% Q(n+1,n+1)=1E-3;
% Q(n+2,n+2)=1E-3;
% [n,m]=size(A);


%normal matrix
N = (A'*(Q^-1)*A);
%least squares solution
x_hat = (N^-1)*A'*(Q^-1)*(y0-b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xb=Xb_approx+x_hat(1:3);
N_hat=x_hat(7:end);
attitude=attitude+x_hat(4:6);


[phi_b, lam_b, h_b] = cart2geod(Xb(1), Xb(2), Xb(3));
for i=1:nRec
    XR(1:3,1,i)=F_s_X(Xb(1),Xb(2),Xb(3),lam_b,phi_b,attitude(2),attitude(1),geometry(1,i),geometry(2,i),attitude(3),geometry(3,i));
end

Qxx = (N^-1);  
[U] = chol(Qxx);    %compute cholesky decomp. for identical comp of vcm purpose
Cxx = U'*U; %find back the vcm of parameter, now the off diag. comp. are identical :)


cov_Xb = Cxx(1:3,1:3); %rover position covariance matrix
cov_ATT = Cxx(4:6,4:6);
cov_N = Cxx(7:end,7:end);


%DOP computation
if (nargout > 7)
    for i=1:nRec
        cov_XYZ = (A((i-1)*length(pr_R)+1:i*length(pr_R),1:3)'*A((i-1)*length(pr_R)+1:i*length(pr_R),1:3))^-1;
        cov_ENU = global2localCov(cov_XYZ, XR(1:3,1,i));
        PDOP(1,i) = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
        HDOP(1,i) = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
        VDOP(1,i) = sqrt(cov_ENU(3,3));
    end
end



