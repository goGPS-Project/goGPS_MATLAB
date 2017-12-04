function [dtR, ISBs, N_hat, var_dtR, var_ISBs, cov_N] = LS_SA_code_phase_clock(pr, ph, snr, elR, distR_approx, sat_pr, sat_ph, dtS, err_tropo, err_iono, phwindup, sys, lambda)

% SYNTAX:
%   [dtR, ISBs, N_hat, var_dtR, var_ISBs, cov_N] = LS_SA_code_phase(pr, ph, snr, elR, distR_approx, sat_pr, sat_ph, dtS, err_tropo, err_iono, phwindup, sys, lambda);
%
% INPUT:
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
%   phwindup     = phase wind-up
%   sys          = array with different values for different systems
%   lambda       = vector containing GNSS wavelengths for available satellites
%
% OUTPUT:
%   dtR = estimated receiver clock
%   ISBs = estimated inter-system biases
%   N_hat = linear combination of ambiguity estimate
%   var_dtR = variance of estimation errors (receiver clock)
%   var_ISBs = variance of estimation errors (inter-system biases)
%   cov_N = covariance matrix of estimation errors (ambiguity values)
%
% DESCRIPTION:
%   Receiver clock estimation by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.


%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------%variable initialization
global sigmaq_cod1 sigmaq_ph

v_light = goGNSS.V_LIGHT;
ISBs = [];
var_ISBs = [];

%data indexes
[~, index] = intersect(sat_pr,sat_ph); %sat_ph is a subset of sat_pr

%number of observations (assuming that sat_ph is a subset of sat_pr)
nsat_pr = length(sat_pr);
nsat_ph = length(sat_ph);
n = nsat_pr + nsat_ph;

%number of unknown parameters
m = 1 + nsat_ph;

% %approximate receiver-satellite distance
% XR_mat = XR_approx(:,ones(n,1))';
% distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));

%design matrix (code)
A = [ zeros(nsat_pr,nsat_ph), ... %column for phase ambiguities   (here zero)
      ones(nsat_pr,1)];    %column for receiver clock delay (multiplied by c)

%design matrix (phase)
A = [A;  diag(-lambda(index)) .* eye(nsat_ph), ...          %column for phase ambiguities
         ones(nsat_ph,1)];             %column for receiver clock delay (multiplied by c)

%if multi-system observations, then estimate an inter-system bias parameter for each additional system
uni_sys = unique(sys(sys ~= 0));
num_sys = length(uni_sys);
ISB = zeros(n,1);
sys = [sys; sys];
if (num_sys > 1)
    m = m + num_sys - 1;
    for s = 2 : num_sys
        ISB(sys == uni_sys(s)) = 1;
        A = [A, ISB];
        ISB = zeros(n,1);
    end
end

%known term vector
b_pr = distR_approx - v_light*dtS + err_tropo + err_iono; %code
b_ph = distR_approx - v_light*dtS + err_tropo - err_iono + lambda.*phwindup; %phase
b = [b_pr; b_ph(index)];

%observation vector
y0 = [pr; lambda(index).*ph(index)];

%observation noise covariance matrix
Q = zeros(n);
Q1 = cofactor_matrix_SA(elR, snr);
Q2 = Q1(index,index);
Q(1:nsat_pr,1:nsat_pr) = sigmaq_cod1 * Q1;
Q(nsat_pr+1:end,nsat_pr+1:end) = sigmaq_ph * Q2;

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x_hat = (N^-1)*A'*(Q^-1)*(y0-b);

%estimated phase ambiguities
N_hat = x_hat(1:nsat_ph);

%estimated receiver clock
dtR = x_hat(nsat_ph+1) / v_light;

%estimated inter-system biases
if (num_sys > 1)
    ISBs = x_hat(nsat_ph+1+[1:num_sys-1]) / v_light;
end

%estimation of the variance of the observation error
y_hat = A*x_hat + b;
v_hat = y0 - y_hat;
sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat * (N^-1);
    cov_N   = Cxx(1:nsat_ph,1:nsat_ph);
    var_dtR = Cxx(nsat_ph+1,nsat_ph+1) / v_light^2;
    if (num_sys > 1)
        var_ISBs = Cxx(nsat_ph+1+[1:num_sys-1],nsat_ph+1+[1:num_sys-1]) / v_light^2;
    end
else
    cov_N   = [];
    var_dtR = [];
end
