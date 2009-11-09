function [xR, Cxx, comb_pr_app, comb_pr_obs, A] = code_double_diff ...
         (posR, pr1_R, posM, pr1_M, time, sat, pivot, Eph, iono)

% SYNTAX:
%   [xR, Cxx, comb_pr_app, comb_pr_obs, A] = code_double_diff ...
%   (posR, pr1_R, posM, pr1_M, time, sat, pivot, Eph, iono);
%
% INPUT:
%   posR = ROVER position (X,Y,Z)
%   pr1_R = ROVER-SATELLITE code pseudorange
%   posM = MASTER position (X,Y,Z)
%   pr1_M = MASTER-SATELLITE code pseudorange
%   time = GPS time
%   sat = visible satellite configuration
%   pivot = pivot satellite
%   cutoff = cut-off on satellite elevation [degrees]
%   Eph = ephemerides matrix
%   iono = ionospheric parameters
%
% OUTPUT:
%   xR = estimated position (X,Y,Z)
%   Cxx = covariance matrix of estimation errors
%   comb_pr_app = crossed approximated pseudoranges
%                 (useful to verify that computations are done correctly)
%   comb_pr_obs = crossed observed pseudoranges
%                 (useful to verify that computations are done correctly)
%   A = design matrix, useful to be passed to the function for the computation
%       of phase double differences, without doing again all the computations
%
% DESCRIPTION:
%   Least squares solution using code double differences.
%   Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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
global a f
global cutoff_init

%rover position coordinates X Y Z
X_R = posR(1);
Y_R = posR(2);
Z_R = posR(3);

%conversion from cartesian to geodetic coordinates
[phi_R, lam_R, h_R] = cart2geod(X_R, Y_R, Z_R);

%master position coordinates X Y Z
X_M = posM(1);
Y_M = posM(2);
Z_M = posM(3);

%conversion from cartesian to geodetic coordinates
[phi_M, lam_M, h_M] = cart2geod(X_M, Y_M, Z_M);

%number of visible satellites
nsat = size(sat,1);

%PIVOT satellite index
i = find(pivot == sat);

%PIVOT position (with clock error and Earth rotation corrections)
posP = sat_corr(Eph, sat(i), time, pr1_R(i), posR);

%computation of ROVER-PIVOT and MASTER-PIVOT approximated pseudoranges
prRP_app = sqrt(sum((posR - posP).^2));
prMP_app = sqrt(sum((posM - posP).^2));

%observed code pseudorange
prRP_obs = pr1_R(i);
prMP_obs = pr1_M(i);

%atmospheric error computation
if (nargin == 9)

   %ROVER-PIVOT and MASTER-PIVOT tropospheric error computation
   [azR, elR, distR] = topocent(posR, posP', a, f);
   [azM, elM, distM] = topocent(posM, posP', a, f);
   err_tropo_RP = err_tropo(elR, h_R);
   err_tropo_MP = err_tropo(elM, h_M);

   %ROVER-PIVOT and MASTER-PIVOT ionospheric error computation
   [phiR, lamR, hR] = cart2geod(posR(1), posR(2), posR(3));
   [phiM, lamM, hM] = cart2geod(posM(1), posM(2), posM(3));
   phiR = phiR * 180 / pi;
   lamR = lamR * 180 / pi;
   phiM = phiM * 180 / pi;
   lamM = lamM * 180 / pi;
   err_iono_RP = err_iono(iono, phiR, lamR, azR, elR, time);
   err_iono_MP = err_iono(iono, phiM, lamM, azM, elM, time);
end

A = [];
tr = [];
io = [];
comb_pr_app = [];
comb_pr_obs = [];

%computation of all linear combinations between PIVOT and other satellites
for i = 1 : nsat
    if (sat(i) ~= pivot)

        %satellite position (with clock error and Earth rotation corrections)
        posS = sat_corr(Eph, sat(i), time, pr1_R(i), posR);

        %computation of the satellite azimuth and elevation
        [azR, elR, distR] = topocent(posR, posS', a, f);
        [azM, elM, distM] = topocent(posM, posS', a, f);

        %cut-off threshold to eliminate too low satellite observations
        if (abs(elR) > cutoff_init)

            %computation of ROVER-PIVOT and MASTER-PIVOT approximated pseudoranges
            prRS_app = sqrt(sum((posR - posS).^2));
            prMS_app = sqrt(sum((posM - posS).^2));

            %observed code pseudoranges
            prRS_obs = pr1_R(i);
            prMS_obs = pr1_M(i);

            %design matrix computation
            A = [A; (((posR(1) - posS(1)) / prRS_app) - ((posR(1) - posP(1)) / prRP_app)) ...
                    (((posR(2) - posS(2)) / prRS_app) - ((posR(2) - posP(2)) / prRP_app)) ...
                    (((posR(3) - posS(3)) / prRS_app) - ((posR(3) - posP(3)) / prRP_app))];

            %computation of crossed approximated pseudoranges
            comb_pr_app = [comb_pr_app; (prRS_app - prMS_app) - (prRP_app - prMP_app)];

            %computation of crossed observed pseudoranges
            comb_pr_obs = [comb_pr_obs; (prRS_obs - prMS_obs) - (prRP_obs - prMP_obs)];

            %computation of crossed atmospheric errors
            if (nargin == 9)

                %computation of tropospheric errors
                err_tropo_RS = err_tropo(elR, h_R);
                err_tropo_MS = err_tropo(elM, h_M);

                %computation of ionospheric errors
                err_iono_RS = err_iono(iono, phiR, lamR, azR, elR, time);
                err_iono_MS = err_iono(iono, phiM, lamM, azM, elM, time);

                %computation of crossed tropospheric errors
                tr = [tr; (err_tropo_RS - err_tropo_MS) - (err_tropo_RP - err_tropo_MP)];

                %computation of crossed ionospheric errors
                io = [io; (err_iono_RS - err_iono_MS) - (err_iono_RP - err_iono_MP)];
            end
        end
    end
end

%vector of the b known term
b = comb_pr_app;

%correction of the b known term
if (nargin == 9)
   b = b + tr - io;
end

%observation vector
y0 = comb_pr_obs;

%number of observations
n = length(y0);

%number of unknown parameters
m = 3;

%observation covariance matrix
Q = 2*ones(n) + 2*eye(n);

%least squares solution
x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);
xR = posR + x;

%estimation of the variance of the observation error
y_stim = A*x + b;
v_stim = y0 - y_stim;
sigma0q_stim = (v_stim' * Q^-1 * v_stim) / (n-m);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma0q_stim * ((A'*Q^-1*A)^-1);
else
    Cxx = [];
end