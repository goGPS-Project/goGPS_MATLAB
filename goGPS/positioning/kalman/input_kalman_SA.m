function [A, prapp_pr1, prapp_ph1, prapp_pr2, prapp_ph2, probs_prIF, probs_phIF, prapp_prIF, prapp_phIF] = input_kalman_SA(XR_approx, XS, pr1, ph1, pr2, ph2, distR_approx, dtS, err_tropo, err_iono1, err_iono2, phwindup, lambda, PCV_S)

% SYNTAX:
%   [A, prapp_pr1, prapp_ph1, prapp_pr2, prapp_ph2, probs_prIF, probs_phIF, prapp_prIF, prapp_phIF] = input_kalman_SA(XR_approx, XS, pr1, ph1, pr2, ph2, distR_approx, dtS, err_tropo, err_iono1, err_iono2, phwindup, lambda, PCV_S);
%
% INPUT:
%   XR_approx = receiver approximate position (X,Y,Z)
%   XS = satellite position (X,Y,Z)
%   pr1 = code pseudorange (carrier L1)
%   ph1 = phase observation (carrier L1)
%   pr2 = code pseudorange (carrier L2)
%   ph2 = phase observation (carrier L2)
%   distR_approx = receiver-satellite approximate range
%   dtS = satellite clock error
%   err_tropo = tropospheric error
%   err_iono1 = ionospheric error (L1 carrier)
%   err_iono2 = ionospheric error (L2 carrier)
%   phwindup = phase wind-up
%   lambda = matrix containing GNSS wavelengths for available satellites
%   PCV_S = satellite PCV corrections (iono-free only)
%
% OUTPUT:
%   A = parameters obtained from the linearization of the observation equation, e.g. (xR-xS)/prRS)
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

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

%variable initialization
v_light = Core_Utils.V_LIGHT;

%design matrix
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
     (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
     (XR_approx(3) - XS(:,3)) ./ distR_approx];    %column for Z coordinate

%approximate pseudoranges
prapp_pr  = distR_approx - v_light*dtS + err_tropo;
prapp_pr1 = prapp_pr + err_iono1;
prapp_ph1 = prapp_pr - err_iono1 + lambda(:,1).*phwindup;
prapp_pr2 = prapp_pr + err_iono2;
prapp_ph2 = prapp_pr - err_iono2 + lambda(:,2).*phwindup;

%observed iono-free combinations
alpha1 = lambda(:,4);
alpha2 = lambda(:,5);
alphat = lambda(:,6);
alphan = lambda(:,7);
probs_prIF  = alpha1 .* pr1 - alpha2 .* pr2; %observed pseudorange (iono-free code)
probs_phIF  = alphat .* ph1 - alphan .* ph2; %observed pseudorange (iono-free phase)

%approximate iono-free combinations (alpha1 - alpha2 = 1)
prapp_prIF = prapp_pr + PCV_S;
prapp_phIF = prapp_pr + (alpha1 .* lambda(:,1) .* phwindup - alpha2 .* lambda(:,2) .* phwindup) + PCV_S;
