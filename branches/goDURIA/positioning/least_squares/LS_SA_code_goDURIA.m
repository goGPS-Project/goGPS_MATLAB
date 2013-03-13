function [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP] = LS_SA_goDURIA(XR_approx_t0,XR_approx_t1, XS_t0,XS_t1,ph_t0,ph_t1, snr_t0,snr_t1, elR, distR_approx_t0,distR_approx_t1, sat_pr, sat_ph, dtS_t0, dtS_t1, err_tropo_t0,err_tropo_t1, err_iono_t0,err_iono_t1, phase)
                                                                                       
% SYNTAX:
%   [XR, dtR, N_hat, cov_XR, var_dtR, cov_N, PDOP, HDOP, VDOP] = LS_SA_code_phase(XR_approx, XS, pr, ph, snr, elR, distR_approx, sat_pr, sat_ph, dtS, err_tropo, err_iono, phase);
%
% INPUT:
%   XR_approx    = receiver approximate position (X,Y,Z)
%   XS           = satellite position (X,Y,Z) (with both code and phase)
%   pr           = code observations
%   ph           = phase observations
%   snr          = signal-to-noise ratio
%   elR          = satellite elevation (vector)
%   distR_approx = approximate receiver-satellite distance (vector)
%   sat_pr       = available satellites
%   sat_ph       = available satellites with phase
%   dtS          = satellite clock error (vector)
%   err_tropo    = tropospheric error
%   err_iono     = ionospheric error
%   phase        = L1 carrier (phase=1), L2 carrier (phase=2)
%
% OUTPUT:
%   pos_R = estimated position (X,Y,Z)
%   cov_pos_R = covariance matrix of estimation errors (rover position)
%   N_stim = linear combination of ambiguity estimate
%   cov_N_stim = covariance matrix of estimation errors (combined ambiguity values)
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%
% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
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

%variable initialization
global v_light
global lambda1 lambda2
global sigmaq_cod1 sigmaq_ph

if (phase == 1)
    lambda = lambda1;
else
    lambda = lambda2;
end

%data indexes
[~, index] = intersect(sat_pr,sat_ph); %sat_ph is a subset of sat_pr

%number of observations (assuming that sat_ph is a subset of sat_pr)
nsat_pr = length(sat_pr);
nsat_ph = length(sat_ph);
n = nsat_ph;

%number of unknown parameters
m = 4;

% %approximate receiver-satellite distance
% XR_mat = XR_approx(:,ones(n,1))';
% distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));

%design matrix (code)



A = [];


for i = 1 : n
%design matrix (phase)

    eij_approx_t0=((XR_approx_t0 - XS_t0(i,:)'))./distR_approx_t0(i);
    eij_approx_t1=((XR_approx_t1 - XS_t1(i,:)'))./distR_approx_t1(i);
    eij_approx=(eij_approx_t0+eij_approx_t1)./2;

    A = [A; eij_approx(1) eij_approx(2) eij_approx(3) 1];

%known term vector
end


    try
    b = [distR_approx_t0'-distR_approx_t1' - v_light*(dtS_t0-dtS_t1)]; %phase
    catch
        distR_approx_t0';
    end
   

%observation vector
    y0 = [ lambda*(ph_t0(index)-ph_t1(index))];

%observation noise covariance matrix
Q = zeros(n);
Q1 = cofactor_matrix_SA(elR, snr_t0);
Q2 = Q1(index,index);
%Q(1:nsat_pr,1:nsat_pr) = sigmaq_cod1 * Q1;
Q(1:nsat_pr,1:nsat_pr) = sigmaq_ph * Q2;

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x_hat = (N^-1)*A'*(Q^-1)*(y0-b);
XR = XR_approx_t0-XR_approx_t1 + x_hat(1:3);

%estimated phase ambiguities


%estimated receiver clock
dtR = x_hat(end) ;

%estimation of the variance of the observation error
y_hat = A*x_hat + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat * (N^-1);
    cov_XR  = Cxx(1:3,1:3);
    cov_N   = Cxx(4:end-1,4:end-1);
    var_dtR = Cxx(end,end);
else
    cov_XR  = [];
    cov_N   = [];
    var_dtR = []; 
end

%DOP computation
if (nargout > 6)
    cov_XYZ = (A(1:nsat_pr,1:3)'*A(1:nsat_pr,1:3))^-1;
    cov_ENU = global2localCov(cov_XYZ, XR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end
