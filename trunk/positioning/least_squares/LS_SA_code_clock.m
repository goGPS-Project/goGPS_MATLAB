function [dtR, var_dtR] = LS_SA_code_clock(pr_R, snr_R, elR, distR, dtS, err_tropo_RS, err_iono_RS, sys)

% SYNTAX:
%   [dtR, var_dtR] = LS_SA_code_clock(pr_R, snr_R, elR, distR, dtS, err_tropo_RS, err_iono_RS, sys);
%
% INPUT:
%   pr_R  = code observations (vector)
%   snr_R = signal-to-noise ratio (vector)
%   elR   = satellite elevation (vector)
%   distR = satellite-receiver distance (vector)
%   dtS   = satellite clock error (vector)
%   err_tropo_RS = tropospheric error (vector)
%   err_iono_RS  = ionospheric error (vector)
%   sys          = array with different values for different systems
%
% OUTPUT:
%   dtR = receiver clock error (scalar)
%   var_dtR = estimate error variance (scalar)
%
% DESCRIPTION:
%   Receiver clock error estimation by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.2 beta
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

v_light = goGNSS.V_LIGHT;

%number of observations
n = length(pr_R);

%number of unknown parameters
m = 1;

%design matrix
A = ones(n,1);

%if multi-system observations, then estimate an inter-system bias parameter for each additional system
uni_sys = unique(sys(sys ~= 0));
num_sys = length(uni_sys);
ISB = zeros(n,1);
if (num_sys > 1)
    m = m + num_sys - 1;
    for s = 2 : num_sys
        ISB(sys == uni_sys(s)) = 1;
        A = [A, ISB];
        ISB = zeros(n,1);
    end
end

%known term vector
b = distR - v_light*dtS + err_tropo_RS + err_iono_RS;

%observation vector
y0 = pr_R;

%observation covariance matrix
Q = cofactor_matrix_SA(elR, snr_R);

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x   = (N^-1)*A'*(Q^-1)*(y0-b);
dtR = x(1) / v_light;

%estimation of the variance of the observation error
y_hat = A*x + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat*(N^-1);
    var_dtR = Cxx(1,1) / v_light;
else
    var_dtR = []; 
end
