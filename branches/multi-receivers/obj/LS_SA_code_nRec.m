function [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = LS_SA_code_nRec(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sat, geometry)

% SYNTAX:
%   [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num] = LS_SA_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS);
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
%   sat = index of usable satellites
%   geometry = receiver coordinates in object RF
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
%
% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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
keyboard

%% compute distances from object coordinates
n_receivers=size(pr_R,2);
if n_receivers>1 
    
    % graph
    baselines=[ones(n_receivers-1,1),(2:n_receivers)'];
    
    % distances
    distance_3D=(sum((geometry(:,2:end)-repmat(geometry(:,1),1,n_receivers-1)).^2,1)).^.5;
    distance_2D=(sum((geometry(1:2,2:end)-repmat(geometry(1:2,1),1,n_receivers-1)).^2,1)).^.5;
    
  
end


%% temp
distance_3D=[distance_3D, sqrt(sum((geometry(:,2)-geometry(:,3)).^2))];
baselines=[baselines;[2 3]];
%%
% y0: observation vector
y0 = NaN(sum(sat(:))+length(distance_3D),1);

%known term vector
b = NaN(sum(sat(:))+length(distance_3D),1);

%observation covariance matrix
Q = zeros(sum(sat(:))+length(distance_3D));

% A: design matrix.. 
A = zeros(sum(sat(:)),4*n_receivers);
n_row_A=0;


% pseuduoranges of all the receivers from all the visible satellites and known distances between receivers
for i=1:n_receivers
    n_sat_i= sum(sat(:,i));
    index_sat_i=find(sat(:,i)==1);    
    
    %satellite-receivers geometry
    A(n_row_A+1:n_row_A+n_sat_i,(i-1)*3+1:i*3)=[(XR_approx(1,i) - XS(index_sat_i,1)) ./ distR_approx(index_sat_i,i), ... %column for X coordinate
     (XR_approx(2,i) - XS(index_sat_i,2)) ./ distR_approx(index_sat_i,i), ... %column for Y coordinate
     (XR_approx(3,i) - XS(index_sat_i,3)) ./ distR_approx(index_sat_i,i)];
 
    %receiver clocks
    A(n_row_A+1:n_row_A+n_sat_i,3*n_receivers+i)=ones(n_sat_i,1);
    
    y0(n_row_A+1:n_row_A+n_sat_i,1) = pr_R(index_sat_i,i);
    b(n_row_A+1:n_row_A+n_sat_i,1) = distR_approx(index_sat_i,i) - goGNSS.V_LIGHT.*dtS(index_sat_i) + err_tropo_RS(index_sat_i,i) + err_iono_RS(index_sat_i,i);
    
    Q(n_row_A+1:n_row_A+n_sat_i,n_row_A+1:n_row_A+n_sat_i)= cofactor_matrix_SA(elR(index_sat_i,i), snr_R(index_sat_i,i));
    
    n_row_A=n_row_A+n_sat_i;
end

% constraints on distances
for i=1:size(distance_3D,2);
    % indexes of the two receivers
    rec_i = baselines(i,1);
    rec_j = baselines(i,2);
    dist_ij_approx = sqrt(sum((XR_approx(:,rec_i)-XR_approx(:,rec_j)).^2));
    A(n_row_A+1,(rec_i-1)*3+1:rec_i*3) = [-(XR_approx(1,rec_i)-XR_approx(1,rec_j))/dist_ij_approx ...
                                          -(XR_approx(2,rec_i)-XR_approx(2,rec_j))/dist_ij_approx ...
                                          -(XR_approx(3,rec_i)-XR_approx(3,rec_j))/dist_ij_approx];
    A(n_row_A+1,(rec_j-1)*3+1:rec_j*3) = [(XR_approx(1,rec_i)-XR_approx(1,rec_j))/dist_ij_approx ...
                                          (XR_approx(2,rec_i)-XR_approx(2,rec_j))/dist_ij_approx ...
                                          (XR_approx(3,rec_i)-XR_approx(3,rec_j))/dist_ij_approx];
    
    y0(n_row_A+1,1)= distance_3D(i);                                 
    b(n_row_A+1,1) = 0;
         
    Q(n_row_A+1,n_row_A+1) = 10^-2;
    n_row_A = n_row_A+1;
end


%Least-Squares adjustement
%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x   = (N^-1)*A'*(Q^-1)*(y0-b);
XR  = XR_approx + reshape(x(1:n_receivers*3),3,n_receivers);
dtR = x(3*n_receivers+1:end) ./ goGNSS.V_LIGHT;

% number of observations (n) and unknown (m)
[n,m] = size(A);

%estimation of the variance of the observation error
y_hat = A*x + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%computation of the condition number on the eigenvalues of N
N_min_eig = min(eig(N));
N_max_eig = max(eig(N));
cond_num = N_max_eig / N_min_eig;

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat * (N^-1);
    cov_XR  = Cxx(1:3*n_receivers,1:3*n_receivers);
    var_dtR = Cxx(3*n_receivers+1:end,3*n_receivers+1:end);
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
