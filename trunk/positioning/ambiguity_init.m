function [N_stim_slip, N_stim_born, var_N_slip, var_N_born] = ambiguity_init(XR_approx, XS, pr_R, pr_M, ...
    ph_R, ph_M, snr_R, snr_M, elR, elM, sat_pr, sat_ph, sat_slip, sat_born, distR_approx, distM, ...
    err_tropo_R, err_tropo_M, err_iono_R, err_iono_M, pivot, lambda, N_kalman, Cee_N_kalman, sigmaq_pos_R)

% SYNTAX:
%   [N_stim_slip, N_stim_born] = ambiguity_init(XR_approx, XS, pr_R, pr_M, ...
%    ph_R, ph_M, snr_R, snr_M, elR, elM, sat_pr, sat_ph, sat_slip, sat_born, distR_approx, distM, ...
%    err_tropo_R, err_tropo_M, err_iono_R, err_iono_M, pivot, lambda, N_kalman, Cee_N_kalman);
%
% INPUT:
%   XR_approx = receiver approximate position (X,Y,Z)
%   XS = satellite positions (X,Y,Z)
%   pr_R = ROVER-SATELLITE code pseudorange
%   pr_M = MASTER-SATELLITE code pseudorange
%   ph_R = ROVER-SATELLITE phase measurement
%   ph_M = MASTER-SATELLITE phase measurement
%   snr_R = ROVER signal-to-noise ratio
%   snr_M = MASTER signal-to-noise ratio
%   el_R = ROVER satellite elevation
%   el_M = MASTER satellite elevation
%   sat_pr = available satellites
%   sat_ph = available satellites with phase
%   sat_slip = slipped satellites
%   sat_born = new satellites
%   distR_approx = ROVER-SATELLITE approximate range
%   distM = MASTER-SATELLITE approximate range
%   err_tropo_R = ROVER-SATELLITE tropospheric error
%   err_tropo_M = MASTER-SATELLITE tropospheric error
%   err_iono_R = ROVER-SATELLITE ionospheric error
%   err_iono_M = MASTER-SATELLITE ionospheric error
%   pivot = pivot satellite
%   lambda = vector containing GNSS wavelengths for available satellites
%   N_kalman = phase ambiguities estimated by Kalman filter
%   Cee_N_kalman = phase ambiguities estimated error
%
% OUTPUT:
%   N_stim_slip = phase ambiguity estimation for slipped satellites
%   N_stim_born = phase ambiguity estimation for new satellites
%   var_N_slip = error variance of phase ambiguity estimation for slipped satellites
%   var_N_born = error variance of phase ambiguity estimation for born satellites
%
% DESCRIPTION:
%   This function estimates phase ambiguities using a least-squares
%   adjustment.

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

global sigmaq_cod1 sigmaq_ph
global amb_restart_method

N_stim_slip = [];
N_stim_born = [];
var_N_slip=[];
var_N_born=[];

%number of slipped satellites
nsat_slip = length(sat_slip);

%number of new satellites
nsat_born = length(sat_born);

%merge new and slipped satellites
sat_amb = [sat_slip; sat_born];
nsat_amb = nsat_slip + nsat_born;

%data indexes
[~, index] = intersect(sat_pr,sat_ph);        %sat_ph is a subset of sat_pr
[~, index_slip] = intersect(sat_ph,sat_slip); %sat_slip is a subset of sat_ph
[~, index_born] = intersect(sat_ph,sat_born); %sat_born is a subset of sat_ph
[~, index_noamb] = setdiff(sat_ph,sat_amb);   %satellites for which the ambiguity is already available
index_amb = [index_slip; index_born];         %satellites for which the ambiguity needs to be estimated
                                              % NOTE: 'intersect' would sort the values, so it can't be used here

%keep only available phase observations
ph_R = ph_R(index);
ph_M = ph_M(index);
lambda = lambda(index);

%remove zeros
index_zero_pr = or(pr_R==0,pr_M==0);
index_zero_ph = or(ph_R==0,ph_M==0);
pr_R(index_zero_pr) = [];
pr_M(index_zero_pr) = [];
ph_R(index_zero_ph) = [];
ph_M(index_zero_ph) = [];
lambda(index_zero_ph) = [];

%number of observations (assuming that sat_ph is a subset of sat_pr)
nsat_pr = length(pr_R);
nsat_ph = length(ph_R);

%pivot indexes
pivot_index_pr = find(pivot == sat_pr);
pivot_index_ph = find(pivot == sat_ph);

%if the selected method is observed code - phase comparison
if (amb_restart_method == 0)
    
    %observed code double differences
    comb_pr = (pr_R - pr_M) - (pr_R(pivot_index_pr) - pr_M(pivot_index_pr));
    
    %observed phase double differences
    comb_ph = (ph_R - ph_M) - (ph_R(pivot_index_ph) - ph_M(pivot_index_ph));
    
    %linear combination of initial ambiguity estimate
    N_stim = comb_pr(index) ./ lambda - comb_ph;
    sigmaq_N_stim = 4*sigmaq_cod1 ./ lambda(1).^2;
    
    %new ambiguity for slipped satellites
    N_stim_slip = N_stim(index_slip);
    var_N_slip = ones(nsat_slip,1)*sigmaq_N_stim;
    
    %new ambiguity for new satellite
    N_stim_born = N_stim(index_born);
    var_N_born = ones(nsat_born,1)*sigmaq_N_stim;
    
%if the number of observations is not sufficient to apply least squares adjustment
%or if the selected method is Kalman-estimated code - phase comparison
elseif (nsat_pr + nsat_ph - 2 <= 3 + nsat_amb) || fix(nsat_ph/2) < nsat_amb || (amb_restart_method == 1)
    
    %KEPT AS A REFERENCE: it should be used in the calling functions and
    %passed as an argument
    %sigmaq_pos_R = diag(T*Cee*T');
    %sigmaq_pos_R = sigmaq_pos_R([1,o1+1,o2+1]);
    
    %observed code double differences
    comb_pr = (distR_approx - distM) - (distR_approx(pivot_index_pr) - distM(pivot_index_pr));
    
    %observed phase double differences
    comb_ph = (ph_R - ph_M) - (ph_R(pivot_index_ph) - ph_M(pivot_index_ph));
    
    %linear combination of initial ambiguity estimate
    N_stim = comb_pr(index) ./ lambda - comb_ph;
    sigmaq_N_stim = sum(sigmaq_pos_R) ./ lambda(1).^2;
    
    %new ambiguity for slipped satellites
    N_stim_slip = N_stim(index_slip);
    var_N_slip = ones(nsat_slip,1)*sigmaq_N_stim;
    
    %new ambiguity for new satellite
    N_stim_born = N_stim(index_born);
    var_N_born = ones(nsat_born,1)*sigmaq_N_stim;

%if the selected method is Least squares adjustment
else
    
    %ambiguity columns in design matrix (lambda positions)
    A_amb = zeros(nsat_ph,nsat_amb);
    for i = 1:nsat_amb
        A_amb(index_amb(i),i) = -lambda(index_amb(i));
    end
    
    %design matrix (code)
    A = [((XR_approx(1) - XS(:,1)) ./ distR_approx) - ((XR_approx(1) - XS(pivot_index_pr,1)) / distR_approx(pivot_index_pr)), ... %column for X coordinate
         ((XR_approx(2) - XS(:,2)) ./ distR_approx) - ((XR_approx(2) - XS(pivot_index_pr,2)) / distR_approx(pivot_index_pr)), ... %column for Y coordinate
         ((XR_approx(3) - XS(:,3)) ./ distR_approx) - ((XR_approx(3) - XS(pivot_index_pr,3)) / distR_approx(pivot_index_pr)), ... %column for Z coordinate
           zeros(nsat_pr, nsat_amb)]; %column for phase ambiguities   (here zero)
    
    %design matrix (phase)
    A = [A; ((XR_approx(1) - XS(index,1)) ./ distR_approx(index)) - ((XR_approx(1) - XS(pivot_index_pr,1)) / distR_approx(pivot_index_pr)), ... %column for X coordinate
            ((XR_approx(2) - XS(index,2)) ./ distR_approx(index)) - ((XR_approx(2) - XS(pivot_index_pr,2)) / distR_approx(pivot_index_pr)), ... %column for Y coordinate
            ((XR_approx(3) - XS(index,3)) ./ distR_approx(index)) - ((XR_approx(3) - XS(pivot_index_pr,3)) / distR_approx(pivot_index_pr)), ... %column for Z coordinate
              A_amb]; %column for phase ambiguities

    %observed pseudoranges
    probs_pr  = (pr_R - pr_M) - (pr_R(pivot_index_pr) - pr_M(pivot_index_pr));                                                                       %observed pseudorange DD (code)
    probs_ph  = (lambda .* ph_R - lambda .* ph_M) - (lambda(pivot_index_ph) * ph_R(pivot_index_ph) - lambda(pivot_index_ph) * ph_M(pivot_index_ph)); %observed pseudorange DD (phase)
    
    %approximate pseudoranges
    prapp    =         (distR_approx - distM)      - (distR_approx(pivot_index_pr) - distM(pivot_index_pr));       %approximate pseudorange DD
    prapp    = prapp + (err_tropo_R - err_tropo_M) - (err_tropo_R(pivot_index_pr)  - err_tropo_M(pivot_index_pr)); %tropospheric error DD
    prapp_pr = prapp + (err_iono_R  - err_iono_M)  - (err_iono_R(pivot_index_pr)   - err_iono_M(pivot_index_pr));  %ionoshperic error DD
    prapp_ph = prapp - (err_iono_R  - err_iono_M)  + (err_iono_R(pivot_index_pr)   - err_iono_M(pivot_index_pr));  %ionoshperic error DD
    
    %remove pivot-pivot lines
    A(pivot_index_pr, :) = [];
    A(nsat_pr - 1 + pivot_index_ph, :) = [];
    probs_pr(pivot_index_pr) = [];
    probs_ph(pivot_index_ph) = [];
    prapp_pr(pivot_index_pr) = [];
    prapp_ph(pivot_index_pr) = [];
    N_kalman(pivot_index_pr) = [];
    lambda(pivot_index_ph)   = [];

    %update indexes
    pos = find(index == pivot_index_pr);
    index(pos) = [];
    index(pos:end) = index(pos:end) - 1;
    pos = find(index_noamb == pivot_index_ph);
    index_noamb(pos) = [];
    index_noamb(pos:end) = index_noamb(pos:end) - 1;
    
    %known term vector
    b = [prapp_pr; prapp_ph(index)];
    
    %observation vector
    N_kalman = N_kalman(index);
    probs_ph(index_noamb) = probs_ph(index_noamb) + lambda(index_noamb).*N_kalman(index_noamb);
    y0 = [probs_pr; probs_ph];

    %number of observations
    n = length(y0);
    
    %observation noise covariance matrix
    Q = zeros(n);
    Q1 = cofactor_matrix(elR, elM, snr_R, snr_M, pivot_index_pr);
    Q2 = Q1(index,index);
    
    Q(1:nsat_pr-1,1:nsat_pr-1) = sigmaq_cod1 * Q1;
%     if (nargin >= 24)
%         Cee_N_kalman(pivot_index_pr,:) = [];
%         Cee_N_kalman(:,pivot_index_pr) = [];
%         %ambiguity estimation error is taken into account (TO BE FIXED: not properly scaled with respect to input code and phase variances)
%         Q(nsat_pr:end,nsat_pr:end) = (sigmaq_ph * eye(nsat_ph - 1) + lambda(1)^2*Cee_N_kalman(index,index)) .* Q2;
%     else
        Q(nsat_pr:end,nsat_pr:end) = sigmaq_ph * Q2;
%     end

    %normal matrix
    N = (A'*(Q^-1)*A);
    
    %least squares solution
    x  = (N^-1)*A'*(Q^-1)*(y0-b);
    
    y_hat = A*x + b;
    v_hat = y0 - y_hat;
    m=size(A,2);
    if n-m > 0
        sigma02_hat = (v_hat'*(Q^-1)*v_hat) / (n-m);    
        Cxx = sigma02_hat * (N^-1);
        cov_XR  = Cxx(1:3,1:3);
        cov_N = Cxx(4:end,4:end);
        var_N = diag(cov_N);
              
        if (nsat_slip ~= 0)
            N_stim_slip = x(4 : 4 + nsat_slip - 1);
            var_N_slip = var_N(1:nsat_slip);
        end
        if (nsat_born ~= 0)
            N_stim_born = x(4 + nsat_slip : 4 + nsat_amb - 1);
            var_N_born = var_N(nsat_slip+1:end);
        end
    end
        
        
        
end
