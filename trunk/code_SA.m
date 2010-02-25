function [xR, Cxx, A] = code_SA(posR, pr1_R, time, Eph, iono)

% SYNTAX:
%   [xR, Cxx, A] = code_SA(posR, pr1_R, time, sat, Eph, iono);
%
% INPUT:
%   posR = ROVER position (X,Y,Z)
%   pr1_R = ROVER code observations (L1 carrier)
%   time = GPS time
%   Eph = ephemerides
%   iono = ionosphere parameters
%
% OUTPUT:
%   xR = estimated position (X,Y,Z)
%   Cxx = estimate error covariance matrix
%   A = design matrix
%
% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
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
global cutoff

%number of visible satellites
sat = find(pr1_R ~= 0);

%cartesian to geodetic conversion of ROVER coordinates
[phiR, lamR, hR] = cart2geod(posR(1), posR(2), posR(3));

A = [];
b = [];
y0 = [];
tr = [];
io = [];

for i = 1 : length(sat)

    %satellite position (with clock error and Earth rotation corrections)
    [posS dtS] = sat_corr(Eph, sat(i), time, pr1_R(sat(i)), posR);

    %computation of the satellite azimuth and elevation
    [azR, elR, distR] = topocent(posR, posS'); %#ok<NASGU>

    %cut-off threshold to eliminate too low satellite observations
    if (elR > cutoff)

        %computation of ROVER-SATELLITE approximated pseudorange
        prRS_app = sqrt(sum((posR - posS).^2));

        %observed code pseudorange
        prRS_obs = pr1_R(sat(i));

        %design matrix computation
        A = [A; ((posR(1) - posS(1)) / prRS_app) ...
                ((posR(2) - posS(2)) / prRS_app) ...
                ((posR(3) - posS(3)) / prRS_app) 1];

        %approximate pseudoranges
        b = [b; prRS_app];

        %observed pseudoranges
        y0 = [y0; prRS_obs + v_light*dtS];

        %computation of atmospheric errors
        if (nargin == 5)

            %computation of tropospheric errors
            err_tropo_RS = err_tropo(elR, hR);

            %computation of ionospheric errors
            err_iono_RS = err_iono(iono, phiR, lamR, azR, elR, time);

            %save tropospheric errors
            tr = [tr; err_tropo_RS];

            %save ionospheric errors
            io = [io; err_iono_RS];
        end
    end
end

%correction of the b known term
if (nargin == 5)
   b = b + tr - io;
end

%number of observations
n = length(y0);

%number of unknown parameters
m = 4;

%least squares solution
x = ((A'*A)^-1)*A'*(y0-b);
xR = posR + x(1:3);

%estimation of the variance of the observation error
y_stim = A*x + b;
v_stim = y0 - y_stim;
sigma0q_stim = (v_stim'* v_stim) / (n-m);

%covariance matrix of the estimation error
Cxx = sigma0q_stim * ((A'*A)^-1);
