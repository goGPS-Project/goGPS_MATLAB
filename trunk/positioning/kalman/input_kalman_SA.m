function [A, prstim_pr1, prstim_ph1, prstim_pr2, prstim_ph2] = input_kalman_SA(XR_approx, XS, distR_approx, dtR, dtS, err_tropo, err_iono, lambda)

% SYNTAX:
%   [A, prstim_pr1, prstim_ph1, prstim_pr2, prstim_ph2] = input_kalman_SA(XR_approx, XS, distR_approx, dtR, dtS, err_tropo, err_iono, lambda);
%
% INPUT:
%   XR_approx = receiver approximate position (X,Y,Z)
%   XS = satellite position (X,Y,Z)
%   distR_approx = receiver-satellite approximate range
%   dtR = receiver clock error
%   dtS = satellite clock error
%   err_tropo = tropospheric error
%   err_iono = ionospheric error
%   lambda = matrix containing GNSS wavelengths for available satellites
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
%                           goGPS v0.4.1 beta
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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
v_light = goGNSS.V_LIGHT;

lambda1 = lambda(:,1);
lambda2 = lambda(:,2);

%design matrix
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
     (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
     (XR_approx(3) - XS(:,3)) ./ distR_approx];    %column for Z coordinate

prstim_pr1 = distR_approx + v_light*(dtR - dtS) + err_tropo + err_iono;
prstim_ph1 = distR_approx + v_light*(dtR - dtS) + err_tropo - err_iono;
prstim_pr2 = distR_approx + v_light*(dtR - dtS) + err_tropo + (lambda2./lambda1).^2 .* err_iono;
prstim_ph2 = distR_approx + v_light*(dtR - dtS) + err_tropo - (lambda2./lambda1).^2 .* err_iono;
