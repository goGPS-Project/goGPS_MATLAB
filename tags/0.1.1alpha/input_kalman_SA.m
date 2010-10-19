function [A, b, err_iono_RS] = input_kalman_SA(posR, pr_Rsat, snr_R, time, sat, Eph, iono)

% SYNTAX:
%   [A, b, err_iono_RS] = input_kalman_SA(posR, pr_Rsat, time, sat, Eph, iono)
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

%variable initialization
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
err_iono_RS = 0;

for i = 1 : nsat
    
    %observed code measurement
    prRS_obs = pr_Rsat(i);
    
    %satellite position (with clock error and Earth rotation corrections)
    [posS dtS] = sat_corr(Eph, sat(i), time, pr_Rsat(i), posR);
    
    %computation of the satellite azimuth and elevation
    [azR, elR(i)] = topocent(posR, posS');
    
    %computation of the ROVER-SATELLITE approximated pseudoranges
    prRS_app = sqrt(sum((posR - posS).^2));
    
    %computation of tropospheric errors
    err_tropo_RS = err_tropo(elR(i), hR);
    
    %if ionospheric parameters are available
    if (nargin == 7)
        
        %computation of ionospheric errors
        err_iono_RS(i) = err_iono(iono, phiR, lamR, azR, elR(i), time);
    end
    
    %construction of the transition matrix
    A = [A; ((posR(1) - posS(1)) / prRS_app) ...
            ((posR(2) - posS(2)) / prRS_app) ...
            ((posR(3) - posS(3)) / prRS_app) 1];

    %approximate pseudoranges
    b = [b; prRS_app - v_light*dtS + err_tropo_RS + err_iono_RS(i)];
    
    %observed pseudoranges
    y0 = [y0; prRS_obs];
end

%observation noise covariance matrix
Q = cofactor_matrix_SA(elR, snr_R, sat);

%least squares solution
x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);

A = A(:,1:3);

b = b + x(4);
