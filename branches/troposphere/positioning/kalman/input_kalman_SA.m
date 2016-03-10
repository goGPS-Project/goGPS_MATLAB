function [A, prapp_pr1, prapp_ph1, prapp_pr2, prapp_ph2, probs_prIF, probs_phIF, prapp_prIF, prapp_phIF] = input_kalman_SA(XR_approx, XS, pr1, ph1, pr2, ph2, distR_approx, dtR, dtS, err_tropo, err_iono1, err_iono2, lambda)

% SYNTAX:
%   [A, prapp_pr1, prapp_ph1, prapp_pr2, prapp_ph2, probs_prIF, probs_phIF, prapp_prIF, prapp_phIF] = input_kalman_SA(XR_approx, XS, pr1, ph1, pr2, ph2, distR_approx, dtR, dtS, err_tropo, err_iono1, err_iono2, lambda);
%
% INPUT:
%   XR_approx = receiver approximate position (X,Y,Z)
%   XS = satellite position (X,Y,Z)
%   pr1 = code pseudorange (carrier L1)
%   ph1 = phase observation (carrier L1)
%   pr2 = code pseudorange (carrier L2)
%   ph2 = phase observation (carrier L2)
%   distR_approx = receiver-satellite approximate range
%   dtR = receiver clock error
%   dtS = satellite clock error
%   err_tropo = tropospheric error
%   err_iono1 = ionospheric error (L1 carrier)
%   err_iono2 = ionospheric error (L2 carrier)
%   lambda = matrix containing GNSS wavelengths for available satellites
%
% OUTPUT:
%   A = parameters obtained from the linearization of the observation equation,
%       e.g. (xR-xS)/prRS)
%   prapp_pr1 = approximate P1 pseudorange (corrected by clock-tropo-iono delays)
%   prapp_ph1 = approximate L1 pseudorange (corrected by clock-tropo-iono delays)
%   prapp_pr2 = approximate P2 pseudorange (corrected by clock-tropo-iono delays)
%   prapp_ph2 = approximate L2 pseudorange (corrected by clock-tropo-iono delays)
%   probs_prIF = observed code (iono-free combination)
%   probs_phIF = observed phase (iono-free combination)
%   prapp_prIF = approximate code (iono-free combination)
%   prapp_phIF = approximate phase (iono-free combination)
%
% DESCRIPTION:
%   This function computes the parameters needed to apply the Kalman filter
%   transition matrix, that links state variables to GPS observations.

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

%variable initialization
v_light = goGNSS.V_LIGHT;

%design matrix
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
     (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
     (XR_approx(3) - XS(:,3)) ./ distR_approx];    %column for Z coordinate

%approximate pseudoranges
prapp_pr  = distR_approx + v_light*(dtR - dtS) + err_tropo;
prapp_pr1 = prapp_pr + err_iono1;
prapp_ph1 = prapp_pr - err_iono1;
prapp_pr2 = prapp_pr + err_iono2;
prapp_ph2 = prapp_pr - err_iono2;

%observed iono-free combinations
alpha1 = (goGNSS.F1^2/(goGNSS.F1^2 - goGNSS.F2^2));
alpha2 = (goGNSS.F2^2/(goGNSS.F1^2 - goGNSS.F2^2));
probs_prIF  = alpha1 *                pr1 - alpha2 *                pr2; %observed pseudorange (iono-free code)
probs_phIF  = alpha1 * lambda(:,1) .* ph1 - alpha2 * lambda(:,2) .* ph2; %observed pseudorange (iono-free phase)

%approximate iono-free combinations (alpha1 - alpha2 = 1)
prapp_prIF = prapp_pr;
prapp_phIF = prapp_pr;
