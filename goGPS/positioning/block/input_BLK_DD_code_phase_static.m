function [y0, A, b, Q] = input_BLK_DD_code_phase_static ...
         (XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda)

% SYNTAX:
%   [y0, A, b, Q] = input_LS_DD_code_phase_batch ...
%   (XR_approx, XM, XS, pr_R, ph_R, snr_R, pr_M, ph_M, snr_M, elR, elM, err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, lambda);
%
% INPUT:
%   XR_approx   = receiver approximate position (X,Y,Z)
%   XM          = master station position (X,Y,Z)
%   XS          = satellite position (X,Y,Z)
%   pr_R        = receiver code observations
%   ph_R        = receiver phase observations
%   pr_M        = master code observations
%   pr_M        = master phase observations
%   snr_R       = receiversignal-to-noise ratio
%   snr_M       = mastersignal-to-noise ratio
%   elR         = satellite elevation (vector)
%   elM         = satellite elevation (vector)
%   err_tropo_R = tropospheric error
%   err_tropo_M = tropospheric error
%   err_iono_R  = ionospheric error
%   err_iono_M  = ionospheric error
%   pivot_index = index identifying the pivot satellite
%   lambda      = vector containing GNSS wavelengths for available satellites
%
% OUTPUT:
%   y0 = observation vector
%   A = design matrix
%   b = known term vector
%   Q = observation covariance matrix
%
% DESCRIPTION:
%   Function that prepares the input matrices for the least squares batch solution.

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
%  Contributors:     Hendy F. Suhandri, ...
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
global sigmaq_cod1 sigmaq_ph

%number of observations
n = 2*length(pr_R);

%approximate receiver-satellite distance
XR_mat = XR_approx(:,ones(n/2,1))';
XM_mat = XM(:,ones(n/2,1))';
distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
distM = sqrt(sum((XS-XM_mat).^2 ,2));

%design matrix (code)
A = [((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index,1)) / distR_approx(pivot_index)), ... %column for X coordinate
     ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index,2)) / distR_approx(pivot_index)), ... %column for Y coordinate
     ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index,3)) / distR_approx(pivot_index))];    %column for Z coordinate

%design matrix (phase)
A = [A; ((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index,1)) / distR_approx(pivot_index)), ... %column for X coordinate
        ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index,2)) / distR_approx(pivot_index)), ... %column for Y coordinate
        ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index,3)) / distR_approx(pivot_index))];    %column for Z coordinate

%known term vector
b    =     (distR_approx - distM)      - (distR_approx(pivot_index) - distM(pivot_index));       %approximate pseudorange DD
b    = b + (err_tropo_R - err_tropo_M) - (err_tropo_R(pivot_index)  - err_tropo_M(pivot_index)); %tropospheric error DD
b_pr = b + (err_iono_R  - err_iono_M)  - (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD (code)
b_ph = b - (err_iono_R  - err_iono_M)  + (err_iono_R(pivot_index)   - err_iono_M(pivot_index));  %ionoshperic error DD (phase)
b = [b_pr; b_ph];

%observation vector
y0_pr =         (pr_R - pr_M) - (pr_R(pivot_index) - pr_M(pivot_index));
y0_ph = lambda.*((ph_R - ph_M) - (ph_R(pivot_index) - ph_M(pivot_index)));
y0 = [y0_pr; y0_ph];

%remove pivot-pivot lines
A( [pivot_index, pivot_index+n/2], :) = [];
b( [pivot_index, pivot_index+n/2])    = [];
y0([pivot_index, pivot_index+n/2])    = [];
n = n - 2;
% n0 = n;

%observation noise covariance matrix
Q = zeros(n);
Q1 = cofactor_matrix(elR, elM, snr_R, snr_M, pivot_index);
Q(1:n/2,1:n/2) = sigmaq_cod1 * Q1;
Q(n/2+1:end,n/2+1:end) = sigmaq_ph * Q1;

%keep only phase observations
y0 = y0(n/2+1:end);
b = b(n/2+1:end);
A = A(n/2+1:end,:);
Q = Q(n/2+1:end,n/2+1:end);
