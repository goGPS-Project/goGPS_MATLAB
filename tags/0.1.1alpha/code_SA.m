function [xR, Cxx, PDOP, HDOP, VDOP, A] = code_SA(posR, pr1_R, snr_R, sat, time, Eph, iono)

% SYNTAX:
%   [xR, Cxx, PDOP, HDOP, VDOP, A] = code_SA(posR, pr1_R, , snr_R, sat, time, sat, Eph, iono);
%
% INPUT:
%   posR = ROVER position (X,Y,Z)
%   pr1_R = ROVER code observations (L1 carrier)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   sat = visible satellite configuration
%   time = GPS time
%   Eph = ephemerides
%   iono = ionosphere parameters
%
% OUTPUT:
%   xR = estimated position (X,Y,Z)
%   Cxx = estimate error covariance matrix
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%   A = design matrix
%
% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.1 alpha
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

global v_light

%number of visible satellites
nsat = size(sat,1);

%ROVER-satellite elevation initialization
elR = zeros(nsat,1);

%cartesian to geodetic conversion of ROVER coordinates
[phiR, lamR, hR] = cart2geod(posR(1), posR(2), posR(3));

%radians to degrees
phiR = phiR * 180 / pi;
lamR = lamR * 180 / pi;

A = [];
b = [];
y0 = [];
tr = [];
io = [];

for i = 1 : nsat
    
    %satellite position (with clock error and Earth rotation corrections)
    [posS dtS] = sat_corr(Eph, sat(i), time, pr1_R(i), posR);

    %computation of the satellite azimuth and elevation
    [azR, elR(i), distR] = topocent(posR, posS'); %#ok<NASGU>
    
    %computation of ROVER-SATELLITE approximated pseudorange
    prRS_app = sqrt(sum((posR - posS).^2));
    
    %observed code pseudorange
    prRS_obs = pr1_R(i);
    
    %design matrix computation
    A = [A; ((posR(1) - posS(1)) / prRS_app) ...
            ((posR(2) - posS(2)) / prRS_app) ...
            ((posR(3) - posS(3)) / prRS_app) 1];
    
    %approximate pseudoranges
    b = [b; prRS_app];
    
    %observed pseudoranges
    y0 = [y0; prRS_obs + v_light*dtS];
    
    %computation of tropospheric errors
    err_tropo_RS = err_tropo(elR(i), hR);
    
    %save tropospheric errors
    tr = [tr; err_tropo_RS];
    
    %if ionospheric parameters are available
    if (nargin == 7)
        
        %computation of ionospheric errors
        err_iono_RS = err_iono(iono, phiR, lamR, azR, elR(i), time);
        
        %save ionospheric errors
        io = [io; err_iono_RS];
    end
end

%correction of the b known term
b = b + tr;
if (nargin == 7)
   b = b + io;
end

%number of observations
n = length(y0);

%number of unknown parameters
m = 4;

%observation covariance matrix
Q = cofactor_matrix_SA(elR, snr_R, sat);

%least squares solution
x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);
xR = posR + x(1:3);

%estimation of the variance of the observation error
y_stim = A*x + b;
v_stim = y0 - y_stim;
sigma0q_stim = (v_stim'* v_stim) / (n-m);

%covariance matrix of the estimation error
Cxx = sigma0q_stim * ((A'*Q^-1*A)^-1);

%DOP computation
if (nargout > 2)
    cov_XYZ = (A(:,1:3)'*A(:,1:3))^-1;
    cov_ENU = global2localCov(cov_XYZ, xR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end