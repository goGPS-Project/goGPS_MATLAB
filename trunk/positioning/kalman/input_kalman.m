function [A, probs_pr1, probs_ph1, prapp_pr1, prapp_ph1, probs_pr2, probs_ph2, prapp_pr2, prapp_ph2] = input_kalman(XR_approx, XS, pr1_R, ph1_R, pr1_M, ph1_M, pr2_R, ph2_R, pr2_M, ph2_M, err_tropo_R, err_iono1_R, err_iono2_R, err_tropo_M, err_iono1_M, err_iono2_M, distR_approx, distM, sat, pivot, lambda)

% SYNTAX:
%   [A, probs_pr1, probs_ph1, prapp_pr1, prapp_ph1, probs_pr2, probs_ph2, prapp_pr2, prapp_ph2] = input_kalman(XR_approx, XS, pr1_R, ph1_R, pr1_M, ph1_M, pr2_R, ph2_R, pr2_M, ph2_M, err_tropo_R, err_iono1_R, err_iono2_R, err_tropo_M, err_iono1_M, err_iono2_M, distR_approx, distM, sat, pivot, lambda);
%
% INPUT:
%   XR_approx = receiver approximate position (X,Y,Z)
%   XS = satellite position (X,Y,Z)
%   pr1_R = ROVER-SATELLITE code pseudorange (carrier L1)
%   ph1_R = ROVER-SATELLITE phase observations (carrier L1)
%   pr1_M = MASTER-SATELLITE code pseudorange (carrier L1)
%   ph1_M = MASTER-SATELLITE phase observations (carrier L1)
%   pr2_R = ROVER-SATELLITE code pseudorange (carrier L2)
%   ph2_R = ROVER-SATELLITE phase observations (carrier L2)
%   pr2_M = MASTER-SATELLITE code pseudorange (carrier L2)
%   ph2_M = MASTER-SATELLITE phase observations (carrier L2)
%   err_tropo_R = ROVER-SATELLITE tropospheric error
%   err_iono1_R  = ROVER-SATELLITE ionospheric error (carrier L1)
%   err_iono2_R  = ROVER-SATELLITE ionospheric error (carrier L2)
%   err_tropo_M = MASTER-SATELLITE tropospheric error
%   err_iono1_M  = MASTER-SATELLITE ionospheric error (carrier L1)
%   err_iono2_R  = ROVER-SATELLITE ionospheric error (carrier L2)
%   distR_approx = ROVER-SATELLITE approximate range
%   distM = MASTER-SATELLITE range
%   sat = configuration of visible satellites
%   pivot = pivot satellite
%   lambda = matrix containing GNSS wavelengths for available satellites
%
% OUTPUT:
%   A = parameters obtained from the linearization of the observation equation,
%       e.g. ((xR-xS)/prRS)-((xR-xP)/prRP)
%   probs_pr1 = observed code double differences (carrier L1)
%   probs_ph1 = observed phase double differences (carrier L1)
%   prapp_pr1 = approximate code double differences (carrier L1)
%   prapp_ph1 = approximate phase double differences (carrier L1)
%   probs_pr2 = observed code double differences (carrier L2)
%   probs_ph2 = observed phase double differences (carrier L2)
%   prapp_pr2 = approximate code double differences (carrier L2)
%   prapp_ph2 = approximate phase double differences (carrier L2)
%
% DESCRIPTION:
%   This function computes the parameters needed to apply the Kalman filter.
%   Transition matrix that link state variables to GPS observations.

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
lambda1 = lambda(:,1);
lambda2 = lambda(:,2);

%pivot search
pivot_index = find(pivot == sat);

%design matrix
A = [((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index,1)) / distR_approx(pivot_index)), ... %column for X coordinate
     ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index,2)) / distR_approx(pivot_index)), ... %column for Y coordinate
     ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index,3)) / distR_approx(pivot_index))];    %column for Z coordinate

%observed pseudoranges
probs_pr1  = (pr1_R - pr1_M) - (pr1_R(pivot_index) - pr1_M(pivot_index));  %observed pseudorange DD (L1 code)
probs_pr2  = (pr2_R - pr2_M) - (pr2_R(pivot_index) - pr2_M(pivot_index));  %observed pseudorange DD (L2 code)
probs_ph1  = (lambda1 .* ph1_R - lambda1 .* ph1_M) - (lambda1(pivot_index) * ph1_R(pivot_index) - lambda1(pivot_index) * ph1_M(pivot_index)); %observed pseudorange DD (L1 phase)
probs_ph2  = (lambda2 .* ph2_R - lambda2 .* ph2_M) - (lambda2(pivot_index) * ph2_R(pivot_index) - lambda2(pivot_index) * ph2_M(pivot_index)); %observed pseudorange DD (L2 phase)

%approximate pseudoranges
prapp_pr  =            (distR_approx - distM)      - (distR_approx(pivot_index) - distM(pivot_index));       %approximate pseudorange DD
prapp_pr  = prapp_pr + (err_tropo_R - err_tropo_M) - (err_tropo_R(pivot_index)  - err_tropo_M(pivot_index)); %tropospheric error DD
prapp_pr1 = prapp_pr + (err_iono1_R - err_iono1_M) - (err_iono1_R(pivot_index)  - err_iono1_M(pivot_index)); %ionoshperic error DD (L1 code)
prapp_ph1 = prapp_pr - (err_iono1_R - err_iono1_M) + (err_iono1_R(pivot_index)  - err_iono1_M(pivot_index)); %ionoshperic error DD (L1 phase)
prapp_pr2 = prapp_pr + (err_iono2_R - err_iono2_M) - (err_iono2_R(pivot_index)  - err_iono2_M(pivot_index)); %ionoshperic error DD (L2 code)
prapp_ph2 = prapp_pr - (err_iono2_R - err_iono2_M) + (err_iono2_R(pivot_index)  - err_iono2_M(pivot_index)); %ionoshperic error DD (L2 phase)

%remove pivot-pivot lines
A(pivot_index, :)      = [];
probs_pr1(pivot_index) = [];
probs_ph1(pivot_index) = [];
probs_pr2(pivot_index) = [];
probs_ph2(pivot_index) = [];
prapp_pr1(pivot_index) = [];
prapp_ph1(pivot_index) = [];
prapp_pr2(pivot_index) = [];
prapp_ph2(pivot_index) = [];
