function [R, N_hat, cov_R, cov_N] = LS_DD_code_phase_MR(XR_approx, F, XM, XS, ...
    pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, ...
    err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda, flag_IAR)
%
% function [R, N_hat, cov_R, cov_N] = LS_DD_code_phase_MR(XR_approx, F, XM, XS, ...
%               pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, ...
%               err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda, flag_IAR);
%
%--------------------------------
%%%%   UPDATE STATUS : 26/02/2014
%--------------------------------
%
% INPUT:
%   XR_approx    = rover approximate position [X;Y;Z] x nrec
%   F            = matrix of rover antennae in body coord. frame [X;Y;Z] x nrec
%   XM           = master station position (X,Y,Z)
%   XS           = satellite position (X,Y,Z)
%   pr_R         = rover code observations
%   ph_R         = rover phase observations
%   snr_R        = rover signal-to-noise ratio
%   pr_M         = master code observations
%   ph_M         = master phase observations
%   snr_M        = master signal-to-noise ratio
%   elR          = rover satellite elevation
%   elM          = master satellite elevation
%   err_tropo_R  = rover tropospheric error
%   err_iono_R   = rover ionospheric error
%   err_tropo_M  = tropospheric error
%   err_iono_M   = ionospheric error
%   pivot_index  = index identifying the pivot satellite
%   lambda       = vector containing GNSS wavelengths for available satellites
%   flag_IAR     = flag to enable/disable integer ambiguity resolution

% OUTPUT:
%   R            = 9-DCM components
%   N_hat        = linear combination of ambiguity estimate
%   cov_XR       = covariance matrix of estimation errors (9-DCM components)
%   cov_N        = covariance matrix of estimation errors (combined ambiguity values)

%
% DESCRIPTION:
%   Multibaseline least squares solution using code and phase double differences.
%   Epoch-by-epoch solution

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
%
% Modified for fixed-multibaseline purpose by Hendy F. Suhandri
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

global sigmaq_cod1 sigmaq_ph

%number of baselines (n)
n = size(XR_approx,2);

%number of (pesudo) DD-obs. in one baseline (m) (because the number of the pivot satellite is still
%included for matrix structure reasons)
m = 2*length(pr_R); %--->INPUT: pr_M = numb. of CA obs. in Master

%number of unknown parameters (u)
%q is span of matrix F, (F = matrix baselines in body system)
if n < 1
    error('LS_DD_code_phase_MR: no available baselines.')
elseif n == 1
    q = 1;
elseif n == 2
    q = 2;
else %n > 2
    q = 3;
end
q;

u = 3*q + (m/2-1)*n;

%approximate receiver-satellite distance
XM_mat  = XM(:,ones(m/2,1))';
distM   = sqrt(sum((XS-XM_mat).^2 ,2));
XR_mat  = zeros(m/2,3,n);
distR_approx = zeros(m/2,n);
for r = 1 : n
    XR_mat(:,:,r) = XR_approx(:,r*ones(m/2,1))';
    distR_approx(:,r) = sqrt(sum((XS-XR_mat(:,:,r)).^2 ,2));
end
% XS
% XR_mat
% XR_approx
% distR_approx

%Coordinate transformation from Cartesion to Geodetic
[phi, lam, h] = cart2geod(XM_mat(1,1), XM_mat(1,2), XM_mat(1,3));

% %Coordinate transformation from ECEF to Local level (NED)
% R_el = [-sin(phi)*cos(lam) -sin(phi)*sin(lam) -cos(phi) ;
%         -sin(lam) cos(lam) 0 ;
%         cos(phi)*cos(lam) cos(phi)*sin(lam) -sin(phi)];

%Coordinate transformation from Local level (NED) to ECEF
R_le = [-sin(phi)*cos(lam) -sin(lam) -cos(phi)*cos(lam);
        -sin(phi)*sin(lam) cos(lam) -cos(phi)*sin(lam);
        cos(phi) 0 -sin(phi)];

%line-of-sight matrix (code)
% G = [((XR_approx(1,1) - XS(:,1)) ./ distR_approx(:,1)) - ((XR_approx(1,1) - XS(pivot_index,1)) / distR_approx(pivot_index,1)), ... %column for X coordinate
%      ((XR_approx(2,1) - XS(:,2)) ./ distR_approx(:,1)) - ((XR_approx(2,1) - XS(pivot_index,2)) / distR_approx(pivot_index,1)), ... %column for Y coordinate
%      ((XR_approx(3,1) - XS(:,3)) ./ distR_approx(:,1)) - ((XR_approx(3,1) - XS(pivot_index,3)) / distR_approx(pivot_index,1))]    %column for Z coordinate];

format long g

G = [-((XS(:,1) - XR_approx(1,1)) ./ distR_approx(:,1)) + ((XS(pivot_index,1) - XR_approx(1,1)) / distR_approx(pivot_index,1)), ... %column for X coordinate
    -((XS(:,2) - XR_approx(2,1)) ./ distR_approx(:,1)) + ((XS(pivot_index,2) - XR_approx(2,1)) / distR_approx(pivot_index,1)), ... %column for Y coordinate
    -((XS(:,3) - XR_approx(3,1)) ./ distR_approx(:,1)) + ((XS(pivot_index,3) - XR_approx(3,1)) / distR_approx(pivot_index,1))];    %column for Z coordinate];


%line-of-sight matrix (code and phase)
G = [G; G];

%wavelength matrix (code and phase)
A = [zeros(m/2,m/2); diag(lambda) .* eye(m/2)];

%known approximated term matrix

%b    =      distR_approx - repmat(distR_approx(pivot_index,:),m/2,1) - repmat(distM,1,n) + repmat(distM(pivot_index),m/2,n)
b    =     (distR_approx - repmat(distM,1,n))       - repmat(distR_approx(pivot_index,:) - repmat(distM(pivot_index),1,n),m/2,1);       %approximate pseudorange DD
b    = b + (err_tropo_R  - repmat(err_tropo_M,1,n)) - repmat(err_tropo_R(pivot_index,:)  - repmat(err_tropo_M(pivot_index),1,n),m/2,1); %tropospheric error DD
b_pr = b + (err_iono_R   - repmat(err_iono_M,1,n))  - repmat(err_iono_R(pivot_index,:)   - repmat(err_iono_M(pivot_index),1,n),m/2,1);  %ionoshperic error DD (code)
b_ph = b - (err_iono_R   - repmat(err_iono_M,1,n))  + repmat(err_iono_R(pivot_index,:)   - repmat(err_iono_M(pivot_index),1,n),m/2,1);  %ionoshperic error DD (phase)

b = [b_pr; b_ph];

%observation matrix
y0_pr = (pr_R - repmat(pr_M,1,n)) - repmat(pr_R(pivot_index,:) - repmat(pr_M(pivot_index),1,n),m/2,1);
y0_ph = repmat(lambda,1,n).*((ph_R - repmat(ph_M,1,n)) - repmat(ph_R(pivot_index,:) - repmat(ph_M(pivot_index),1,n),m/2,1));

y0 = [y0_pr; y0_ph];

%vectorize known term and observation matrices
b  = reshape(b,  numel(b), 1);
y0 = reshape(y0, numel(b), 1);

%pivot indices for vectorized matrices
pivot_index_vec = zeros(n,1);
for r = 1 : n
    pivot_index_vec(r) = pivot_index+(m/2)*(r-1);
end
pivot_index_vec = [pivot_index_vec; pivot_index_vec+(m/2*n)];

%remove pivot-pivot lines
A( [pivot_index, pivot_index+m/2],:)  = [];
A( :, pivot_index)  = [];
G( [pivot_index, pivot_index+m/2], :)  = [];
b(  pivot_index_vec) = [];
y0( pivot_index_vec) = [];
m = (m - 2);

%"observations minus known terms"
%y = vec2mat((y0 - b),n);
y = reshape((y0 - b),m,n);

%observation noise covariance matrix
Q  = zeros(m);
Q1 = cofactor_matrix(elR(:,1), elM, snr_R(:,1), snr_M, pivot_index);
Q(1:m/2,1:m/2) = sigmaq_cod1 * Q1;
Q(m/2+1:end,m/2+1:end) = sigmaq_ph * Q1;
Q;
P  = 0.5*(eye(n)+ones(n));
Qy = kron(P,Q);

%antenna baselines in body coord. frame
for k = 1 : length(F(1,:))
    for i = 1 : length(F(:,1))
        if (F(i,k) ~= 0)
            Fn(i,k) = F(i,k);
        end
    end
end

Fn;

%attitude-based float solution
Gp = (eye(m) - A*(A'*Q^-1*A)^-1*A'*Q^-1)*G;
Rfloat = R_le'*((Gp'*Q^-1*Gp)^-1*Gp'*Q^-1*y*P^-1*Fn'*(Fn*P^-1*Fn')^-1); %floated-Rotation Matrix 
Zfloat = (A'*Q^-1*A)^-1*A'*Q^-1*(y - G*Rfloat*Fn); %floated-ambiguities
Q_vecR = kron((Fn*(P^-1)*Fn')^-1,(Gp'*(Q^-1)*Gp)^-1);
Q_vecZ = (kron(P^-1,A'*Q^-1*A) - kron(P^-1*Fn',A'*Q^-1*G)*(kron(Fn*P^-1*Fn',G'*Q^-1*G))^-1*kron(Fn*P^-1,G'*Q^-1*A))^-1;

%float rotation angles for single baseline case
if n  == 1
    yaw   = atand(Rfloat(2,1)/Rfloat(1,1))
    pitch = atand(-Rfloat(3,1)/sqrt(Rfloat(1,1)^2+Rfloat(2,1)^2))

elseif n = 2
    yaw   = atand(Rfloat(2,1)/Rfloat(1,1))
    pitch = atand(-Rfloat(3,1)/sqrt(Rfloat(1,1)^2+Rfloat(2,1)^2))
    roll  = asind(Rfloat(3,2)/cosd(pitch))
    
else
    yaw   = atand(Rfloat(2,1)/Rfloat(1,1))
    pitch = atand(-Rfloat(3,1)/sqrt(Rfloat(1,1)^2+Rfloat(2,1)^2))
    roll  = atand(Rfloat(3,2))/Rfloat(3,3)
end

if (~flag_IAR)
    %apply least squares solution
    R = Rfloat;
%
%     %estimated double difference ambiguities (without PIVOT)
%     N_hat_nopivot = vecZfloat;
%
%     %add a zero at PIVOT position
%     N_hat = zeros((m/2)*n+1,1);
%     N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
%     N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
%
%     %covariance matrix of the estimation error
%     if (m*n > u)
%
%         %rover position covariance matrix
%         cov_R = QvecR;
%
%         %combined ambiguity covariance matrix
%         cov_N_nopivot = QvecZ;
%
%         %add one line and one column (zeros) at PIVOT position
%         cov_N = zeros((m/2)*n+1);
%         cov_N(1:pivot_index-1,1:pivot_index-1)     = cov_N_nopivot(1:pivot_index-1,1:pivot_index-1);
%         cov_N(pivot_index+1:end,pivot_index+1:end) = cov_N_nopivot(pivot_index:end,pivot_index:end);
%     else
%         cov_R = [];
%         cov_N = [];
%     end
%
% else %apply LAMBDA
%
%     if (m*n > u && sigmaq_ph ~= 1e30)
%
%         [vecZcheck, QvecZhat] = lambdafix(vecZfloat, QvecZ);
%
%         R     = ((G'*(Q^-1)*G)^-1)*G'*(Q^-1)*(y-A*vecZcheck)*(P^-1)*Fn'*((Fn*(P^-1)*Fn')^-1);
%         cov_R = kron((Fn*(P^-1)*Fn')^-1,(G'*(Q^-1)*G)^-1); % = QRR - QRZr*(QZrZr^-1)*QZrR;
%
%
%         %estimated double difference ambiguities (without PIVOT)
%         N_hat_nopivot = vecZcheck;
%
%         %add a zero at PIVOT position
%         N_hat = zeros(m*n/2+1,1);
%         N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
%         N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
%
%         %combined ambiguity covariance matrix
%         cov_N_nopivot = QvecZhat;
%
%         %add one line and one column (zeros) at PIVOT position
%         cov_N = zeros(n/2+1);
%         cov_N(1:pivot_index-1,1:pivot_index-1)     = cov_N_nopivot(1:pivot_index-1,1:pivot_index-1);
%         cov_N(pivot_index+1:end,pivot_index+1:end) = cov_N_nopivot(pivot_index:end,pivot_index:end);
%
%     else
%         R = Rfloat;
%         N_hat_nopivot = vecZfloat;
%
%         %add a zero at PIVOT position
%         N_hat = zeros((m/2)*n+1,1);
%         N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
%         N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
%
%         cov_R = [];
%         cov_N = [];
%     end
end
