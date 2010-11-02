function [N_stim_slip, N_stim_born, dtR] = amb_estimate_LS_SA(posR_app, posS, dtS, pr_R, ph_R, snr_R, ...
    el_R, sat, sat_slip, sat_born, prRS_app, err_tropo_RS, err_iono_RS, phase, N_kalman, Cee_N_kalman)

% SYNTAX:
%   [N_stim_slip, N_stim_born, dt_R] = amb_estimate_LS_SA(posR_app, posS, dtS, pr_R, ph_R, snr_R, ...
%    el_R, sat, sat_slip, sat_born, prRS_app, err_tropo_RS, err_iono_RS, phase,
%    N_kalman, Cee_N_kalman);
%
% INPUT:
%   posR_app = receiver position (X,Y,Z)
%   posS = satellite positions (X,Y,Z)
%   dtS = satellite clock error
%   pr_R = ROVER-SATELLITE code pseudorange
%   ph_R = ROVER-SATELLITE phase measurement
%   snr_R = ROVER signal-to-noise ratio
%   el_R = ROVER satellite elevation
%   sat = available satellites
%   sat_slip = slipped satellites
%   sat_born = new satellites
%   prRS_app = approximate pseudorange
%   err_tropoRS = tropospheric error
%   err_ionoRS = ionospheric error
%   phase = GPS frequency selector
%   N_kalman = phase ambiguities estimated by Kalman filter
%   Cee_N_kalman = phase ambiguities estimated error
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
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
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
% global clock_delay_thresh

if (phase == 1)
    lambda = lambda1;
else
    lambda = lambda2;
end

%number of visible satellites
nsat = size(sat,1);

%number of slipped satellites
nsat_slip = size(sat_slip,1);

%number of new satellites
nsat_born = size(sat_born,1);

%merge new and slipped satellites
sat_amb = [sat_slip; sat_born];
nsat_amb = nsat_slip + nsat_born;

A = [];
b = [];
y0 = [];
tr = [];
io = [];

%CODE
for i = 1 : nsat

    %observed code measurement
    prRS_obs = pr_R(i);

    %design matrix computation
    A = [A; ((posR_app(1) - posS(1,i)) / prRS_app(i)) ...
            ((posR_app(2) - posS(2,i)) / prRS_app(i)) ...
            ((posR_app(3) - posS(3,i)) / prRS_app(i)) ...
            zeros(1,nsat_amb) 1];

    %approximate pseudoranges
    b = [b; prRS_app(i) - v_light*dtS(i)];

    %observed pseudoranges
    y0 = [y0; prRS_obs];

    %save tropospheric errors
    tr = [tr; err_tropo_RS(i)];

    %save ionospheric errors
    io = [io; err_iono_RS(i)];
end

% index vector of satellites with phase observations
p = [];

%PHASE
for i = 1 : nsat
    if (ph_R(i) ~= 0)

        %index vector
        p = [p; i];
        
        %observed phase measurement
        phRS_obs = ph_R(i);
        
        pos = find(sat(i) == sat_amb);
        
        %ambiguity vector in design matrix (lambda position)
        N_row = zeros(1, nsat_amb);
        if (~isempty(pos))
            N_row(pos) = -lambda;
        end
        
        %design matrix computation
        A = [A; ((posR_app(1) - posS(1,i)) / prRS_app(i)) ...
            ((posR_app(2) - posS(2,i)) / prRS_app(i)) ...
            ((posR_app(3) - posS(3,i)) / prRS_app(i)) ...
            N_row 1];
        
        %approximate pseudoranges
        b = [b; prRS_app(i) - v_light*dtS(i)];
        
        %observed pseudoranges
        if (~isempty(pos))
            y0 = [y0; lambda * phRS_obs];
        else
            y0 = [y0; lambda * (phRS_obs + N_kalman(i))];
        end
        
        %save tropospheric errors
        tr = [tr; err_tropo_RS(i)];
        
        %save ionospheric errors
        io = [io; -err_iono_RS(i)];
    end
end

%correction of the b known term
b = b + tr + io;

%number of observations
n = length(y0);

%number of unknown parameters
% m = 4 + nsat;

%observation noise covariance matrix
Q = zeros(n);
Q1 = cofactor_matrix_SA(el_R, snr_R, sat);
Q2 = Q1(p,p);

Q(1:nsat,1:nsat) = sigmaq_cod1 * Q1;
if (nargin == 11)
    %ambiguity estimation error is taken into account (TO BE FIXED: not properly scaled
    %with respect to input code and phase variances)
    Q(nsat+1:end,nsat+1:end) = (sigmaq_ph * eye(n - nsat) + lambda^2*Cee_N_kalman(p,p)) .* Q2;
else
    Q(nsat+1:end,nsat+1:end) = sigmaq_ph * Q2;
end

% A_cod = A(1:nsat,:);
% Q_cod = Q(1:nsat,1:nsat);
% y0_cod = y0(1:nsat);
% b_cod = b(1:nsat);

%least squares solution using only code
% x_cod = ((A_cod'*Q_cod^-1*A_cod)^-1)*A_cod'*Q_cod^-1*(y0_cod-b_cod);

%least squares solution
x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);

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
