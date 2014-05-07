function [XR, XR_a, dtR, dtR_a, cov_XR, cov_XR_a, var_dtR, var_dtR_a] = LS_SA_phase_variometric_a(XR_approx_t0, XR_approx_t1,XR_approx_t2, XS_t0, XS_t1,XS_t2, ph_t0, ph_t1,ph_t2, snr, elR, distR_approx_t0, distR_approx_t1,distR_approx_t2, sat_ph, dtS_t0, dtS_t1,dtS_t2, err_tropo_t0, err_tropo_t1, err_tropo_t2, err_iono_t0, err_iono_t1,err_iono_t2, phase, sys, lambda)
                                                                                    
% SYNTAX:
%   [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP] = LS_SA_phase_variometric(XR_approx_t0, XR_approx_t1, XS_t0, XS_t1, ph_t0, ph_t1, snr, elR, distR_approx_t0, distR_approx_t1, sat_ph, dtS_t0, dtS_t1, err_tropo_t0, err_tropo_t1, err_iono_t0, err_iono_t1, phase, sys, lambda);
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
%   sys          = array with different values for different systems
%   lambda       = vector containing GNSS wavelengths for available satellites
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
%                           goGPS v0.4.2 beta
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
global sigmaq_ph

v_light = goGNSS.V_LIGHT;

%number of observations
nsat_ph = length(sat_ph);
n = nsat_ph;

%number of unknown parameters
m = 4;

% %approximate receiver-satellite distance
% XR_mat = XR_approx(:,ones(n,1))';
% distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));

%design matrix (phase)
A = [];

for i = 1 : n

    eij_approx_t0=((XR_approx_t0 - XS_t0(i,:)'))./distR_approx_t0(i);
    eij_approx_t1=((XR_approx_t1 - XS_t1(i,:)'))./distR_approx_t1(i);
    eij_approx=(eij_approx_t0 + eij_approx_t1)./2;
    
    A = [A; eij_approx(1) eij_approx(2) eij_approx(3) 1];
end

%if multi-system observations, then estimate an inter-system bias parameter for each additional system
% uni_sys = unique(sys(sys ~= 0));
% num_sys = length(uni_sys);
% ISB = zeros(n,1);
% if (num_sys > 1)
%     m = m + num_sys - 1;
%     for s = 2 : num_sys
%         ISB(sys == uni_sys(s)) = 1;
%         A = [A, ISB];
%         ISB = zeros(n,1);
%     end
% end

%known term vector
try
    b = distR_approx_t1' - distR_approx_t0' - v_light*(dtS_t1 - dtS_t0); %velocity estimation
    b_a = distR_approx_t2' - 2.*distR_approx_t1' + distR_approx_t0' - v_light*(dtS_t2 - 2.*dtS_t1 + dtS_t0); %acceleration estimation
    b = b + (err_tropo_t1 - err_tropo_t0) - (err_iono_t1 - err_iono_t0);
    b_a = b_a + (err_tropo_t2 - 2.*err_tropo_t1 + err_tropo_t0) - (err_iono_t2 - 2.*err_iono_t1 + err_iono_t0);
catch
    size(distR_approx_t0') %DEBUG
    size(distR_approx_t1')
    size(dtS_t0)
    size(dtS_t1)
end
   
%observation vector
y0 = lambda.*(ph_t1 - ph_t0); %velocity
y0_a = lambda.*(ph_t2 - 2.*ph_t1 + ph_t0); %acceleration
%observation noise covariance matrix
Q1 = cofactor_matrix_SA(elR, snr);
Q  = sigmaq_ph * Q1;

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x_hat = (N^-1)*A'*(Q^-1)*(y0-b); %velocity
x_hat_a = (N^-1)*A'*(Q^-1)*(y0_a-b_a); %acceleration
XR = XR_approx_t1-XR_approx_t0 + x_hat(1:3);
XR_a = XR_approx_t2-2*XR_approx_t1+XR_approx_t0 + x_hat_a(1:3);
%XR = x_hat;
%XR_a = x_hat_a;

%estimated receiver clock drift
dtR = x_hat(4)/v_light;
dtR_a = x_hat_a(4)/v_light;

%estimation of the variance of the observation error velocity
y_hat = A*x_hat + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%estimation of the variance of the observation error accelerations
y_hat_a = A*x_hat_a + b_a;
v_hat_a = y0_a - y_hat_a;
sigma02_hat_a = (v_hat_a'*(Q^-1)*v_hat_a) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat * (N^-1);
    cov_XR  = Cxx(1:3,1:3);
    var_dtR = Cxx(end,end);
    Cxx_a = sigma02_hat_a * (N^-1);
    cov_XR_a  = Cxx_a(1:3,1:3);
    var_dtR_a = Cxx_a(end,end);
else
    cov_XR  = [];
    var_dtR = []; 
    cov_XR_a  = [];
    var_dtR_a = [];
end

%DOP computation
% if (nargout > 6)
%     cov_XYZ = (A(1:nsat_ph,1:3)'*A(1:nsat_ph,1:3))^-1;
%     cov_ENU = global2localCov(cov_XYZ, XR);
%     
%     PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
%     HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
%     VDOP = sqrt(cov_ENU(3,3));
% end
