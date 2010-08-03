function [A, b, err_iono_RS] = input_kalman_SA(pos_R_app, pr_R, ph_R, snr_R, sat, time, Eph, phase, iono, N_kalman, Cee_N_kalman)

% SYNTAX:
%   [A, b, err_iono_RS] = input_kalman_SA(pos_R_app, pr_R, ph_R, snr_R,
%   sat, time, Eph, phase, iono, N_kalman, Cee_N_kalman);
%
% INPUT:
%   posR = receiver position (X,Y,Z)
%   pr_Rsat = ROVER-SATELLITE code pseudorange
%   snr_R = ROVER signal-to-noise ratio
%   time = GPS time
%   sat = configuration of visible satellites
%   Eph = ephemerides matrix
%   iono = ionospheric parameters
%
% OUTPUT:
%   A = parameters obtained from the linearization of the observation equation,
%       e.g. (xR-xS)/prRS)
%   b = approximated code range
%   err_iono_RS = ionospheric error (for all the satellites)
%
% DESCRIPTION:
%   This function computes the parameters needed to apply the Kalman filter.
%   Transition matrix that link state variables to GPS observations.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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
global clock_delay_thresh

if (phase == 1)
    lambda = lambda1;
else
    lambda = lambda2;
end

%number of visible satellites
nsat = size(sat,1);

%ROVER-satellite elevation initialization
elR = zeros(nsat,1);

%cartesian to geodetic conversion of ROVER coordinates
[phiR, lamR, hR] = cart2geod(pos_R_app(1), pos_R_app(2), pos_R_app(3));

%radians to degrees
phiR = phiR * 180 / pi;
lamR = lamR * 180 / pi;

for i = 1 : nsat
    %satellite position (with clock error and Earth rotation corrections)
    [posS(i,:) dtS(i)] = sat_corr(Eph, sat(i), time, pr_R(i), pos_R_app);

    %computation of the satellite azimuth and elevation
    [azR, elR(i)] = topocent(pos_R_app, posS(i,:));
    
    %computation of ROVER-SATELLITE approximated pseudorange
    prRS_app(i) = sqrt(sum((pos_R_app - posS(i,:)').^2));
    
    %computation of tropospheric errors
    err_tropo_RS(i) = err_tropo(elR(i), hR);
    
    %computation of ionospheric errors
    err_iono_RS(i) = err_iono(iono, phiR, lamR, azR, elR(i), time);
end

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
    A = [A; ((pos_R_app(1) - posS(i,1)) / prRS_app(i)) ...
            ((pos_R_app(2) - posS(i,2)) / prRS_app(i)) ...
            ((pos_R_app(3) - posS(i,3)) / prRS_app(i)) 1];

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
        
        %design matrix computation
        A = [A; ((pos_R_app(1) - posS(i,1)) / prRS_app(i)) ...
            ((pos_R_app(2) - posS(i,2)) / prRS_app(i)) ...
            ((pos_R_app(3) - posS(i,3)) / prRS_app(i)) 1];
        
        %approximate pseudoranges
        b = [b; prRS_app(i) - v_light*dtS(i)];
        
        %observed pseudoranges
        y0 = [y0; lambda * (phRS_obs + N_kalman(i))];
        
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
Q1 = cofactor_matrix_SA(elR, snr_R, sat);
Q2 = Q1(p,p);

Q(1:nsat,1:nsat) = sigmaq_cod1 * Q1;
if (nargin == 11)
    %ambiguity estimation error is taken into account (TO BE FIXED: not properly scaled
    %with respect to input code and phase variances)
    Q(nsat+1:end,nsat+1:end) = (sigmaq_ph * eye(n - nsat) + lambda^2*Cee_N_kalman(p,p)) .* Q2;
else
    Q(nsat+1:end,nsat+1:end) = sigmaq_ph * Q2;
end

A_cod = A(1:nsat,:);
Q_cod = Q(1:nsat,1:nsat);
y0_cod = y0(1:nsat);
b_cod = b(1:nsat);

%least squares solution using only code
x_cod = ((A_cod'*Q_cod^-1*A_cod)^-1)*A_cod'*Q_cod^-1*(y0_cod-b_cod);

%least squares solution
x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);

%test on differences between code and code+phase receiver clock delay estimation
while (abs(x_cod(4) - x(4)) > clock_delay_thresh)
    
    %delete phase observation with maximum error variance
    [null_m, i] = max(diag(Q2));
    Q2(i,:) = [];
    Q2(:,i) = [];
    
    Q(nsat + i, :) = [];
    Q(:, nsat + i) = [];
    
    A(nsat + i, :) = [];
    y0(nsat + i) = [];
    b(nsat + i) = [];
    
    %least squares solution
    x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);
end

A = A(1:nsat, 1:3);

b = b(1:nsat) + x(4);
