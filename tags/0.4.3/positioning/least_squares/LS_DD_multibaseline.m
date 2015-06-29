function [R, N_hat, cov_R, cov_N] = LS_DD_multibaseline(Ant, F, XS, XM, ...
              XR1_approx, XR2_approx, XR3_approx, XR4_approx, pr_R1, pr_R2, ...
              pr_R3, pr_R4, ph_R1, ph_R2, ph_R3, ph_R4, pr_M, ph_M, snr_R1, ...
              snr_M, elR1, elM, err_tropo_R1, err_tropo_R2, err_tropo_R3, ...
              err_tropo_R4, err_iono_R1, err_iono_R2, err_iono_R3, err_iono_R4, ...
              err_tropo_M, err_iono_M, pivot_index, lambda)
%
% function [R, N_hat, cov_R, cov_N] = LS_DD_multibaseline(Ant, F, XS, XM, ...
%               XR1_approx, XR2_approx, XR3_approx, XR4_approx, pr_R1, pr_R2, ...
%               pr_R3, pr_R4, ph_R1, ph_R2, ph_R3, ph_R4, pr_M, ph_M, snr_R1, ...
%               snr_M, elR1, elM, err_tropo_R1, err_tropo_R2, err_tropo_R3, ...
%               err_tropo_R4, err_iono_R1, err_iono_R2, err_iono_R3, err_iono_R4, ...
%               err_tropo_M, err_iono_M, pivot_index, lambda)
%
% INPUT:
%   Ant          = number of antenna
%   F            = matrix of antennae in body coord. frame
%   XS           = satellite position (X,Y,Z)
%   XM           = master station position (X,Y,Z)
%   XR1_approx   = rover-1 approximate position (X,Y,Z)
%   XR2_approx   = rover-2 approximate position (X,Y,Z)
%   XR3_approx   = rover-3 approximate position (X,Y,Z)
%   XR4_approx   = rover-4 approximate position (X,Y,Z)
%   pr_R1        = rover-1 code observations
%   pr_R2        = rover-2 code observations
%   pr_R3        = rover-3 code observations
%   pr_R4        = rover-4 code observations
%   ph_R1        = rover-1 phase observations
%   ph_R2        = rover-2 phase observations
%   ph_R3        = rover-3 phase observations
%   ph_R4        = rover-4 phase observations
%   pr_M         = master code observations?
%   ph_M         = master phase observations
%   snr_R1       = rover-1 signal-to-noise ratio
%   snr_M        = mastersignal-to-noise ratio
%   elR1         = rover-1 satellite elevation (vector)
%   elM          = master satellite elevation (vector)
%   err_tropo_R1 = rover-1 tropospheric error
%   err_tropo_R2 = rover-2 tropospheric error
%   err_tropo_R3 = rover-3 tropospheric error
%   err_tropo_R4 = rover-4 tropospheric error
%   err_iono_R1  = rover-1 ionospheric error
%   err_iono_R2  = rover-2 ionospheric error
%   err_iono_R3  = rover-3 ionospheric error
%   err_iono_R4  = rover-4 ionospheric error
%   err_tropo_M  = tropospheric error
%   err_iono_M   = ionospheric error
%   pivot_index  = index identifying the pivot satellite
%   lambda       = vector containing GNSS wavelengths for available satellites

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

% %variable initialization
% global sigmaq_cod1 sigmaq_ph

%number of baseline (n)
n = Ant-1; %---> INPUT: ant = numbers of antenna      

%number of DD-obs. in one baseline (m)
m = 2*(length(pr_M)); %--->INPUT: pr_M = numb. of CA obs. in Master

%number of unknown parameters (u)
 %q is span of matrix F, (F=matrix baselines in body syst.)
if n < 1;
    error('no baseline!')
elseif n == 1;
    q = 1;
elseif n == 2;
    q = 2;
else n > 2;
    q = 3;
end
q

u = 3*q + (m/2 - 1)*n


%approximate receiver-satellite distance
XM_mat  = XM(:,ones(m/2,1))';
distM = sqrt(sum((XS-XM_mat).^2 ,2));
XR1_mat = XR1_approx(:,ones(m/2,1))';
XR2_mat = XR2_approx(:,ones(m/2,1))';
XR3_mat = XR3_approx(:,ones(m/2,1))';
XR4_mat = XR4_approx(:,ones(m/2,1))';
distR1_approx = sqrt(sum((XS-XR1_mat).^2 ,2));
distR2_approx = sqrt(sum((XS-XR1_mat).^2 ,2));
distR3_approx = sqrt(sum((XS-XR1_mat).^2 ,2));
distR4_approx = sqrt(sum((XS-XR1_mat).^2 ,2));

%collecting all approximate and observation value vectors
pr_R         = [pr_R1, pr_R2, pr_R3, pr_R4];
ph_R         = [ph_R1, ph_R2, ph_R3, ph_R4];
distR_approx = [distR1_approx1, distR2_approx, distR3_approx, distR4_approx];
err_tropo_R  = [err_tropo_R1, err_tropo_R2, err_tropo_R3, err_tropo_R4];
err_iono_R   = [err_iono_R1, err_iono_R2, err_iono_R3, err_iono_R4];


%line-of-sign matrix (code)
G = [((XR1_approx(1) - XS(:,1)) ./ distR1_approx) - ((XR1_approx(1) - XS(pivot_index,1)) / distR1_approx(pivot_index)), ... %column for X coordinate
     ((XR1_approx(2) - XS(:,2)) ./ distR1_approx) - ((XR1_approx(2) - XS(pivot_index,2)) / distR1_approx(pivot_index)), ... %column for Y coordinate
     ((XR1_approx(3) - XS(:,3)) ./ distR1_approx) - ((XR1_approx(3) - XS(pivot_index,3)) / distR1_approx(pivot_index))];    %column for Z coordinate]; 

%line-of-sign matrix (phase)
G = [G; ((XR1_approx(1) - XS(:,1)) ./ distR1_approx) - ((XR1_approx(1) - XS(pivot_index,1)) / distR1_approx(pivot_index)), ... %column for X coordinate
        ((XR1_approx(2) - XS(:,2)) ./ distR1_approx) - ((XR1_approx(2) - XS(pivot_index,2)) / distR1_approx(pivot_index)), ... %column for Y coordinate
        ((XR1_approx(3) - XS(:,3)) ./ distR1_approx) - ((XR1_approx(3) - XS(pivot_index,3)) / distR1_approx(pivot_index))];    %column for Z coordinate

%wavelength matrix (code & carrier)
A = [zeros(m/2,m/2); diag(-lambda) .* eye(m/2)];

%known term vector Rover 1
b   =     (distR_approx - repmat(distM,1,4)) - (distR_approx(pivot_index) - repmat(distM(pivot_index),1,4));       %approximate pseudorange DD
b   = b + (err_tropo_R - repmat(err_tropo_M,1,4)) - (err_tropo_R(pivot_index) - repmat(err_tropo_M(pivot_index),1,4)); %tropospheric error DD
b_pr= b + (err_iono_R  - repmat(err_iono_M,1,4))  - (err_iono_R(pivot_index) - repmat(err_iono_M(pivot_index),1,4));  %ionoshperic error DD (code)
b_ph= b - (err_iono_R  - repmat(err_iono_M,1,4))  + (err_iono_R(pivot_index) - repmat(err_iono_M(pivot_index),1,4));  %ionoshperic error DD (phase)

b = cat(1,b_pr, b_ph);

%observation vector
y0_pr =         (pr_R - repmat(pr_M,1,4)) - (pr_R(pivot_index) - repmat(pr_M(pivot_index),1,4));
y0_ph = lambda.*((ph_R - repmat(ph_M,1,4)) - (ph_R(pivot_index) - repmat(ph_M(pivot_index),1,4)));
y0 = cat(1,y0_pr, y0_ph);

%"observations minus known terms"
y = (y0 - b);

%remove pivot-pivot lines
A( [pivot_index, pivot_index+m/2], :)  = [];
G(            :, pivot_index+3)        = [];
b( [pivot_index, pivot_index+m/2]) = [];
y0([pivot_index, pivot_index+m/2]) = [];
m = (m - 2);
 
%observation noise covariance matrix
Q  = zeros(m);
Q1 = cofactor_matrix(elR1, elM, snr_R1, snr_M, pivot_index);
Q(1:m/2,1:m/2) = sigmaq_cod1 * Q1;
Q(m/2+1:end,m/2+1:end) = sigmaq_ph * Q1;

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
Gp = (eye(m) - A*((A'*(Q^-1)*A)^-1)*A'*(Q^-1))*G;
Pfn= Fn'*((Fn*(P^-1)*Fn')^-1)*Fn*(P^-1);
Pg = G*((G'*(Q^-1)*A)^-1)*A'*(Q^-1);
% Ap = (eye(m) - G*((G'*(Q^-1)*G)^-1)*G'*(Q^-1))*A;

Rfloat    = ((Gp'*(Q^-1)*Gp)^-1)*Gp'*(Q^-1)*y*(P^-1)*Fn'*((Fn*(P^-1)*Fn')^-1);
vecZfloat = (kron(eye(n),((A'*(Q^-1)*A)^-1)*A'*(Q^-1)))*(y(:) - (kron(Fn',G))*Rfloat(:));

QvecR = kron((Fn*(P^-1)*Fn')^-1,(Gp*(Q^-1)*Gp)^-1);
%QZrR  = kron(-F'*(F*(P^-1)*F)^-1,((A'*(Q^-1)*A)^-1)*A'*(Q^-1)*G*((Gp'*(Q^-1)*Gp)^-1));
%QRZr  = QZrR';
QvecZ = ((kron(P^-1,A'*(Q^-1)))*(eye(m*n)-kron(Pfn,Pg))*(kron(eye(n),A)))^-1; 
 
if (~flag_LAMBDA)
    %apply least squares solution
    R = Rfloat;
    
    %estimated double difference ambiguities (without PIVOT)
    N_hat_nopivot = vecZfloat;
    
    %add a zero at PIVOT position
    N_hat = zeros((m/2)*n+1,1);
    N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
    N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
    
    %covariance matrix of the estimation error
    if (m*n > u)
        
        %rover position covariance matrix
        cov_R = QvecR;
        
        %combined ambiguity covariance matrix
        cov_N_nopivot = QvecZ;
        
        %add one line and one column (zeros) at PIVOT position
        cov_N = zeros((m/2)*n+1);
        cov_N(1:pivot_index-1,1:pivot_index-1)     = cov_N_nopivot(1:pivot_index-1,1:pivot_index-1);
        cov_N(pivot_index+1:end,pivot_index+1:end) = cov_N_nopivot(pivot_index:end,pivot_index:end);
    else
        cov_R = [];
        cov_N = [];
    end
    
else %apply LAMBDA
        
    if (m*n > u && sigmaq_ph ~= 1e30)
                        
        [vecZcheck, QvecZhat] = lambdafix(vecZfloat, QvecZ);
        
        R     = ((G'*(Q^-1)*G)^-1)*G'*(Q^-1)*(y-A*vecZcheck)*(P^-1)*Fn'*((Fn*(P^-1)*Fn')^-1);
        cov_R = kron((Fn*(P^-1)*Fn')^-1,(G'*(Q^-1)*G)^-1); % = QRR - QRZr*(QZrZr^-1)*QZrR;
        
        
        %estimated double difference ambiguities (without PIVOT)
        N_hat_nopivot = vecZcheck;
        
        %add a zero at PIVOT position
        N_hat = zeros(m*n/2+1,1);
        N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
        N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
        
        %combined ambiguity covariance matrix
        cov_N_nopivot = QvecZhat;
        
        %add one line and one column (zeros) at PIVOT position
        cov_N = zeros(n/2+1);
        cov_N(1:pivot_index-1,1:pivot_index-1)     = cov_N_nopivot(1:pivot_index-1,1:pivot_index-1);
        cov_N(pivot_index+1:end,pivot_index+1:end) = cov_N_nopivot(pivot_index:end,pivot_index:end);
        
    else
        R = Rfloat;
        N_hat_nopivot = vecZfloat;
        
        %add a zero at PIVOT position
        N_hat = zeros((m/2)*n+1,1);
        N_hat(1:pivot_index-1)   = N_hat_nopivot(1:pivot_index-1);
        N_hat(pivot_index+1:end) = N_hat_nopivot(pivot_index:end);
        
        cov_R = [];
        cov_N = [];
    end
end
 

