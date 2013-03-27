function [XR, cov_XR, PDOP, HDOP, VDOP] = LS_DD_phase_variometric(XR_approx_t0, XR_approx_t1, XS_t0, XS_t1, ph_t0, ph_t1, snr, elR, distR_approx_t0, distR_approx_t1, sat_ph, dtS_t0, dtS_t1, err_tropo_t0, err_tropo_t1, err_iono_t0, err_iono_t1, pivot_index, phase)
                                                                                       
% SYNTAX:
%   [XR, cov_XR, PDOP, HDOP, VDOP] = LS_DD_phase_variometric(XR_approx_t0, XR_approx_t1, XS_t0, XS_t1, ph_t0, ph_t1, snr, elR, distR_approx_t0, distR_approx_t1, sat_ph, dtS_t0, dtS_t1, err_tropo_t0, err_tropo_t1, err_iono_t0, err_iono_t1, pivot_index, phase);
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
global lambda1 lambda2
global sigmaq_ph

if (phase == 1)
    lambda = lambda1;
else
    lambda = lambda2;
end

%number of observations
nsat_ph = length(sat_ph);
n = nsat_ph;

%number of unknown parameters
m = 3;

% %approximate receiver-satellite distance
% XR_mat = XR_approx(:,ones(n,1))';
% distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));

%design matrix (phase)
A = [];

for i = 1 : n
    
    if (i ~= pivot_index)
        eij_approx_t0 = ((XR_approx_t0 - XS_t0(i,:)')./distR_approx_t0(i)) - ((XR_approx_t0 - XS_t0(pivot_index,:)')./distR_approx_t0(pivot_index));
        eij_approx_t1 = ((XR_approx_t1 - XS_t1(i,:)')./distR_approx_t1(i)) - ((XR_approx_t1 - XS_t1(pivot_index,:)')./distR_approx_t1(pivot_index));
        eij_approx    = (eij_approx_t0 + eij_approx_t1)./2;
        
        A = [A; eij_approx(1) eij_approx(2) eij_approx(3)];
    end
end

%known term vector
% try
    b = (distR_approx_t0' - distR_approx_t1') - (distR_approx_t0(pivot_index) - distR_approx_t1(pivot_index));
% catch
%     distR_approx_t0'; %DEBUG
% end
   
%observation vector
y0 = lambda*((ph_t0 - ph_t1) - (ph_t0(pivot_index) - ph_t1(pivot_index)));

%remove pivot-pivot lines
b(pivot_index) = [];
y0(pivot_index) = [];
elR(pivot_index) = [];
snr(pivot_index) = [];
n = n - 1;

%observation noise covariance matrix
Q1 = cofactor_matrix_SA(elR, snr);
Q  = sigmaq_ph * Q1;

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x_hat = (N^-1)*A'*(Q^-1)*(y0-b);
XR = XR_approx_t0 - XR_approx_t1 + x_hat;

%estimation of the variance of the observation error
y_hat = A*x_hat + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat * (N^-1);
    cov_XR  = Cxx(1:3,1:3);
else
    cov_XR  = [];
end

%DOP computation
if (nargout > 2)
    cov_XYZ = (A'*A)^-1;
    cov_ENU = global2localCov(cov_XYZ, XR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end
