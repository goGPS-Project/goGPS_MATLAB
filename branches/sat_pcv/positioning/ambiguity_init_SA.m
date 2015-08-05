function [N_stim_slip, N_stim_born, dtR] = ambiguity_init_SA(XR_approx, XS, dtS, pr, ph, snr, ...
    elR, sat_pr, sat_ph, sat_slip, sat_born, distR_approx, err_tropo, err_iono, sys, lambda, N_kalman, Cee_N_kalman)

% SYNTAX:
%   [N_stim_slip, N_stim_born, dt_R] = ambiguity_init_SA(XR_approx, XS, dtS, pr, ph, snr, ...
%    elR, sat, sat_slip, sat_born, distR_approx, err_tropo, err_iono, sys, lambda, N_kalman, Cee_N_kalman);
%
% INPUT:
%   XR_approx = receiver approximate position (X,Y,Z)
%   XS = satellite positions (X,Y,Z)
%   dtS = satellite clock error
%   pr = code pseudorange
%   ph = phase measurement
%   snr = signal-to-noise ratio
%   elR = satellite elevation
%   sat_pr = available satellites
%   sat_ph = available satellites with phase
%   sat_slip = slipped satellites
%   sat_born = new satellites
%   distR_approx = approximate range
%   err_tropo = tropospheric error
%   err_iono = ionospheric error
%   sys = array with different values for different systems
%   lambda = vector containing GNSS wavelengths for available satellites
%   N_kalman = phase ambiguities estimated by Kalman filter  *** same size as ph ***
%   Cee_N_kalman = phase ambiguities estimated error  *** same size as ph ***
%
% OUTPUT:
%   N_stim_slip = phase ambiguity estimation for slipped satellites
%   N_stim_born = phase ambiguity estimation for new satellites
%   dtR = receiver clock error estimation
%
% DESCRIPTION:
%   This function estimates the phase ambiguities and the 
%   receiver clock error using a least-squares adjustment.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
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

global sigmaq_cod1 sigmaq_ph
% global clock_delay_thresh

v_light = goGNSS.V_LIGHT;

%remove zeros
index_zero_pr = (pr == 0);
index_zero_ph = (ph == 0);
pr(index_zero_pr) = [];
ph(index_zero_ph) = [];
lambda(index_zero_ph) = [];

%number of observations (assuming that sat_ph is a subset of sat_pr)
nsat_pr = length(pr);
nsat_ph = length(ph);
n = nsat_pr + nsat_ph;

%number of slipped satellites
nsat_slip = size(sat_slip,1);

%number of new satellites
nsat_born = size(sat_born,1);

%merge new and slipped satellites
sat_amb = [sat_slip; sat_born];
nsat_amb = nsat_slip + nsat_born;

%data indexes
[~, index] = intersect(sat_pr,sat_ph);        %sat_ph is a subset of sat_pr
[~, index_slip] = intersect(sat_ph,sat_slip); %sat_slip is a subset of sat_ph
[~, index_born] = intersect(sat_ph,sat_born); %sat_born is a subset of sat_ph
[~, index_noamb] = setdiff(sat_ph,sat_amb);   %satellites for which the ambiguity is already available
index_amb = [index_slip; index_born];         %satellites for which the ambiguity needs to be estimated
                                              % NOTE: 'intersect' would sort the values, so it can't be used here

%ambiguity columns in design matrix (lambda positions)
A_amb = zeros(nsat_ph,nsat_amb);
for i = 1:nsat_amb
    A_amb(index_amb(i),i) = -lambda(index_amb(i));
end

%design matrix (code)
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ...        %column for X coordinate
     (XR_approx(2) - XS(:,2)) ./ distR_approx, ...        %column for Y coordinate
     (XR_approx(3) - XS(:,3)) ./ distR_approx, ...        %column for Z coordinate
      zeros(nsat_pr,nsat_amb), ...     %column for phase ambiguities   (here zero)
      ones(nsat_pr,1)];         %column for receiver clock delay (multiplied by c)

%design matrix (phase)
A = [A; (XR_approx(1) - XS(index,1)) ./ distR_approx(index), ...  %column for X coordinate
        (XR_approx(2) - XS(index,2)) ./ distR_approx(index), ...  %column for Y coordinate
        (XR_approx(3) - XS(index,3)) ./ distR_approx(index), ...  %column for Z coordinate
         A_amb, ...                                          %column for phase ambiguities
         ones(nsat_ph,1)];              %column for receiver clock delay (multiplied by c)

%if multi-system observations, then estimate an inter-system bias parameter for each additional system
uni_sys = unique(sys(sys ~= 0));
num_sys = length(uni_sys);
ISB = zeros(n,1);
if (num_sys > 1)
    for s = 2 : num_sys
        ISB(sys == uni_sys(s)) = 1;
        A = [A, ISB];
        ISB = zeros(n,1);
    end
end

%known term vector
b_pr = distR_approx - v_light*dtS + err_tropo + err_iono; %code
b_ph = distR_approx - v_light*dtS + err_tropo - err_iono; %phase
b = [b_pr; b_ph(index)];

%observation vector
ph(index_noamb) = ph(index_noamb) + N_kalman(index_noamb);
y0 = [pr; lambda.*ph];

%observation noise covariance matrix
Q = zeros(n);
Q1 = cofactor_matrix_SA(elR, snr);
Q2 = Q1(index,index);

Q(1:nsat_pr,1:nsat_pr) = sigmaq_cod1 * Q1;
% if (nargin >= 18)
%     %ambiguity estimation error is taken into account (TO BE FIXED: not properly scaled
%     %with respect to input code and phase variances)
%     Q(nsat_pr+1:end,nsat_pr+1:end) = (sigmaq_ph * eye(nsat_ph) + lambda.^2.*Cee_N_kalman) .* Q2;
% else
    Q(nsat_pr+1:end,nsat_pr+1:end) = sigmaq_ph * Q2;
% end

% A_cod = A(1:nsat,:);
% Q_cod = Q(1:nsat,1:nsat);
% y0_cod = y0(1:nsat);
% b_cod = b(1:nsat);

%least squares solution using only code
% x_cod = ((A_cod'*Q_cod^-1*A_cod)^-1)*A_cod'*Q_cod^-1*(y0_cod-b_cod);

%normal matrix
N = (A'*(Q^-1)*A);

%least squares solution
x  = (N^-1)*A'*(Q^-1)*(y0-b);

%test on differences between code and code+phase receiver clock delay estimation
% while (abs(x_cod(4) - x(4)) > clock_delay_thresh)
%     
%     %delete phase observation with maximum error variance
%     [null_m, i] = max(diag(Q2));
%     Q2(i,:) = [];
%     Q2(:,i) = [];
%     
%     Q(nsat + i, :) = [];
%     Q(:, nsat + i) = [];
%     
%     A(nsat + i, :) = [];
%     y0(nsat + i) = [];
%     b(nsat + i) = [];
%     
%     %least squares solution
%     x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);
% end

if (nsat_slip ~= 0)
    N_stim_slip = x(4 : 4 + nsat_slip - 1);
else
    N_stim_slip = [];
end
if (nsat_born ~= 0)
    N_stim_born = x(4 + nsat_slip : 4 + nsat_amb - 1);
else
    N_stim_born = [];
end

dtR = x(end) / v_light;
