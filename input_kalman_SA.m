function [A, prstim_pr1, prstim_ph1, prstim_pr2, prstim_ph2] = input_kalman_SA(posR_app, posS, prRS_app, dtR, dtS, err_tropo_RS, err_iono_RS)

% SYNTAX:
%   [A, prstim1, prstim2] = input_kalman_SA(posR_app, posS, prRS_app, dtR, dtS, err_tropo_RS, err_iono_RS);
%
% INPUT:
%   posR_app = receiver position (X,Y,Z)
%   posS = satellite position (X,Y,Z)
%   prRS_app = approximated pseudorange
%
% OUTPUT:
%   A = parameters obtained from the linearization of the observation equation,
%       e.g. (xR-xS)/prRS)
%   prstim_pr1 = approximate P1 pseudorange (corrected by clock-tropo-iono delays)
%   prstim_ph1 = approximate L1 pseudorange (corrected by clock-tropo-iono delays)
%   prstim_pr2 = approximate P2 pseudorange (corrected by clock-tropo-iono delays)
%   prstim_ph2 = approximate L2 pseudorange (corrected by clock-tropo-iono delays)
%
% DESCRIPTION:
%   This function computes the parameters needed to apply the Kalman filter
%   transition matrix, that links state variables to GPS observations.

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

%number of visible satellites
nsat = size(posS,1);

A = [];

for i = 1 : nsat

    %design matrix computation
    A = [A; ((posR_app(1) - posS(i,1)) / prRS_app(i)) ...
            ((posR_app(2) - posS(i,2)) / prRS_app(i)) ...
            ((posR_app(3) - posS(i,3)) / prRS_app(i))];
end

prstim_pr1 = prRS_app + v_light*(dtR - dtS) + err_tropo_RS + err_iono_RS;
prstim_ph1 = prRS_app + v_light*(dtR - dtS) + err_tropo_RS - err_iono_RS;
prstim_pr2 = prRS_app + v_light*(dtR - dtS) + err_tropo_RS + (lambda2/lambda1)^2 * err_iono_RS;
prstim_ph2 = prRS_app + v_light*(dtR - dtS) + err_tropo_RS - (lambda2/lambda1)^2 * err_iono_RS;