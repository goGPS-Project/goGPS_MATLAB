function [pos_R, cov_pos_R, N_stim, cov_N_stim, PDOP, HDOP, VDOP] = code_phase_SA(pos_R_app, pr_R, ph_R, snr_R, sat, time, Eph, phase, iono)

% SYNTAX:
%   [pos_R, cov_pos_R, N_stim, cov_N_stim, PDOP, HDOP, VDOP] = code_phse_SA(pos_R_app, pr_R, ph_R, snr_R, sat, time, Eph, phase, iono);
%
% INPUT:
%   pos_R_app = ROVER position (X,Y,Z)
%   pr_R = ROVER code observations (L1 carrier)
%   ph_R = ROVER-SATELLITE phase observation
%   snr_R = ROVER signal-to-noise ratio
%   sat = visible satellite configuration
%   time = GPS time
%   Eph = ephemerides
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%   iono = ionosphere parameters
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
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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
    [posS(i,:) dtS(i)] = sat_corr(Eph, sat(i), time, pr_R(i));

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
            ((pos_R_app(3) - posS(i,3)) / prRS_app(i)) ...
            zeros(1, nsat) 1];

    %approximate pseudoranges
    b = [b; prRS_app(i) - v_light*dtS(i)];

    %observed pseudoranges
    y0 = [y0; prRS_obs];

    %save tropospheric errors
    tr = [tr; err_tropo_RS(i)];
        
    %save ionospheric errors
    io = [io; err_iono_RS(i)];
end

%PHASE
for i = 1 : nsat

    %observed phase measurement
    phRS_obs = ph_R(i);

    %ambiguity vector in design matrix (lambda position)
    N_row = zeros(1, nsat);
    N_row(i) = -lambda;
    
    %design matrix computation
    A = [A; ((pos_R_app(1) - posS(i,1)) / prRS_app(i)) ...
            ((pos_R_app(2) - posS(i,2)) / prRS_app(i)) ...
            ((pos_R_app(3) - posS(i,3)) / prRS_app(i)) ...
            N_row 1];

    %approximate pseudoranges
    b = [b; prRS_app(i) - v_light*dtS(i)];

    %observed pseudoranges
    y0 = [y0; phRS_obs*lambda];

    %save tropospheric errors
    tr = [tr; err_tropo_RS(i)];

    %save ionospheric errors
    io = [io; -err_iono_RS(i)];
end

%correction of the b known term
b = b + tr + io;

%number of observations
n = length(y0);

%number of unknown parameters
m = 4 + nsat;

%observation noise covariance matrix
Q = zeros(n);
Q1 = cofactor_matrix_SA(elR, snr_R, sat);
Q(1:n/2,1:n/2) = sigmaq_cod1 * Q1;
Q(n/2+1:end,n/2+1:end) = sigmaq_ph * Q1;

%parameter vector
xR = [pos_R_app; zeros(nsat,1); 0];

%least squares solution
x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);

try
    xR = xR + x;
catch
    %DEBUG
    fprintf('Size of parameter vector: %d\n', size(xR,1));
    fprintf('Size of observation vector : %d\n', size(y0,1));
    fprintf('Number of satellites: %d\n', nsat);
    xR = xR + x;
end

%estimation of the variance of the observation error
y_stim = A*x + b;
v_stim = y0 - y_stim;
sigma0q_stim = (v_stim'* Q^-1 * v_stim) / (n-m);

%estimated rover position
pos_R = xR(1:3);

%estimated combined ambiguity values
N_stim = xR(4:end-1);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma0q_stim * ((A'*Q^-1*A)^-1);

    %rover position covariance matrix
    cov_pos_R = Cxx(1:3,1:3);

    %combined ambiguity covariance matrix
    cov_N_stim = Cxx(4:end-1,4:end-1);

else
    cov_pos_R = [];

    cov_N_stim = [];
end

%DOP computation
if (nargout > 4)
    cov_XYZ = (A(1:n/2,1:3)'*A(1:n/2,1:3))^-1;
    cov_ENU = global2localCov(cov_XYZ, xR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end