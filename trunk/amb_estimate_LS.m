function [N_stim_slip, N_stim_born] = amb_estimate_LS(posR_app, posS, pr_R, pr_M, ...
    ph_R, ph_M, snr_R, snr_M, elR, elM, sat_pr, sat, sat_slip, sat_born, prRS_app, prMS_app, ...
    err_tropo_RS, err_tropo_MS, err_iono_RS, err_iono_MS, pivot, phase, N_kalman, Cee_N_kalman)

% SYNTAX:
%   [N_stim_slip, N_stim_born] = amb_estimate_LS(posR_app, posS, pr_R, pr_M, ...
%   ph_R, ph_M, snr_R, snr_M, elR, elM, sat_pr, sat, sat_slip, sat_born, prRS_app, prMS_app, ...
%   err_tropo_RS, err_tropo_MS, err_iono_RS, err_iono_MS, pivot, phase, N_kalman, Cee_N_kalman);
%
% INPUT:
%   posR_app = receiver position (X,Y,Z)
%   posS = satellite positions (X,Y,Z)
%   dtS = satellite clock error
%   pr_R = ROVER-SATELLITE code pseudorange
%   ph_R = ROVER-SATELLITE phase measurement
%   snr_R = ROVER signal-to-noise ratio
%   el_R = ROVER satellite elevation
%   sat_pr = available satellites
%   sat = available satellites with phase
%   sat_slip = slipped satellites
%   sat_born = new satellites
%   prRS_app = approximate pseudorange
%   err_tropoRS = tropospheric error
%   err_ionoRS = ionospheric error
%   pivot = ID of pivot satellite
%   phase = GPS frequency selector
%   N_kalman = phase ambiguities estimated by Kalman filter
%   Cee_N_kalman = phase ambiguities estimated error
%
% OUTPUT:
%   N_stim_slip = phase ambiguity estimation for slipped satellites
%   N_stim_born = phase ambiguity estimation for new satellites
%
% DESCRIPTION:
%   This function estimates phase ambiguities using a least-squares
%   adjustment.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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
global lambda1 lambda2
global sigmaq_cod1 sigmaq_ph

if (phase == 1)
    lambda = lambda1;
else
    lambda = lambda2;
end

%number of visible satellites
nsat_pr = size(sat_pr,1);

%number of visible satellites with phase
nsat = size(sat,1);

%PIVOT satellite index
j = find(pivot == sat_pr);

%PIVOT satellite position
posP = posS(:,j);

%ROVER-PIVOT and MASTER-PIVOT observed code measurements
prRP_obs = pr_R(j);
prMP_obs = pr_M(j);

%ROVER-PIVOT and MASTER-PIVOT observed phase measurements
phRP_obs = ph_R(j);
phMP_obs = ph_M(j);

%ROVER-PIVOT and MASTER-PIVOT approximate code measurements
prRP_app = prRS_app(j);
prMP_app = prMS_app(j);

%ROVER-PIVOT and MASTER-PIVOT tropospheric errors
err_tropo_RP = err_tropo_RS(j);
err_tropo_MP = err_tropo_MS(j);

%ROVER-PIVOT and MASTER-PIVOT ionospheric errors
err_iono_RP = err_iono_RS(j);
err_iono_MP = err_iono_MS(j);

%number of slipped satellites
nsat_slip = size(sat_slip,1);

%number of new satellites
nsat_born = size(sat_born,1);

%merge new and slipped satellites
sat_amb = [sat_slip; sat_born];
nsat_amb = nsat_slip + nsat_born;

%if the number of observations is not sufficient to apply least squares adjustment
if (nsat_pr + nsat - 2 <= 3 + nsat - 1)
    %observed code double differences
    comb_pr = (pr_Rsat(sat_slip) - pr_Msat(sat_slip)) - (pr_RP - pr_MP);
    
    %observed phase double differences
    comb_ph = (ph_Rsat(sat_slip) - ph_Msat(sat_slip)) - (ph_RP - ph_MP);
    
    %linear combination of ambiguities
    N_stim_slip = ((comb_pr - comb_ph * lambda)) / lambda;
    %sigmaq_N_stim_slip = 4*sigmaq_cod1 / lambda^2;
    
    %observed code double differences
    comb_pr = (pr_Rsat(sat_born) - pr_Msat(sat_born)) - (pr_RP - pr_MP);
    
    %observed phase double differences
    comb_ph = (ph_Rsat(sat_born) - ph_Msat(sat_born)) - (ph_RP - ph_MP);
    
    %linear combination of ambiguities
    N_stim_born = ((comb_pr - comb_ph * lambda)) / lambda;
    %sigmaq_N_stim_born = 4*sigmaq_cod1 / lambda^2;
else
    A = [];
    tr = [];
    io = [];
    comb_pr_app = [];
    comb_pr_obs = [];
    
    %CODE
    for i = 1 : nsat_pr
        if (sat_pr(i) ~= pivot)
            
            %observed code measurement
            prRS_obs = pr_R(i);
            prMS_obs = pr_M(i);
            
            %design matrix computation
            A = [A; (((posR_app(1) - posS(1,i)) / prRS_app(i)) - ((posR_app(1) - posP(1)) / prRP_app)) ...
                (((posR_app(2) - posS(2,i)) / prRS_app(i)) - ((posR_app(2) - posP(2)) / prRP_app)) ...
                (((posR_app(3) - posS(3,i)) / prRS_app(i)) - ((posR_app(3) - posP(3)) / prRP_app)) ...
                zeros(1, nsat_amb)];
            
            %computation of crossed approximate pseudoranges
            comb_pr_app = [comb_pr_app; (prRS_app(i) - prMS_app(i)) - (prRP_app - prMP_app)];
            
            %computation of crossed observed code pseudoranges
            comb_pr_obs = [comb_pr_obs; (prRS_obs - prMS_obs) - (prRP_obs - prMP_obs)];
            
            %computation of crossed tropospheric errors
            tr = [tr; (err_tropo_RS(i) - err_tropo_MS(i)) - (err_tropo_RP - err_tropo_MP)];
            
            %computation of crossed ionospheric errors
            io = [io; (err_iono_RS(i) - err_iono_MS(i)) - (err_iono_RP - err_iono_MP)];
        end
    end
    
    % index vector of satellites with phase observations (over nsat_pr)
    p = [];
    % index vector of satellites with phase observations (over nsat_amb)
    r = [];
    
    %PHASE
    for i = 1 : nsat_pr
        if (ph_R(i) ~= 0) & (ph_M(i) ~= 0) & (sat_pr(i) ~= pivot)
            
            %index vectors
            p = [p; i];
            if (i < j)
                r = [r; i];
            else
                r = [r; i-1];
            end
            
            %observed phase measurement
            phRS_obs = ph_R(i);
            phMS_obs = ph_M(i);
            
            pos = find(sat_pr(i) == sat_amb);
            
            %ambiguity vector in design matrix (lambda position)
            N_row = zeros(1, nsat_amb);
            if (~isempty(pos))
                N_row(pos) = -lambda;
            end
            
            %design matrix computation
            A = [A; (((posR_app(1) - posS(1,i)) / prRS_app(i)) - ((posR_app(1) - posP(1)) / prRP_app)) ...
                (((posR_app(2) - posS(2,i)) / prRS_app(i)) - ((posR_app(2) - posP(2)) / prRP_app)) ...
                (((posR_app(3) - posS(3,i)) / prRS_app(i)) - ((posR_app(3) - posP(3)) / prRP_app)) ...
                N_row];
            
            %computation of crossed approximated pseudoranges
            comb_pr_app = [comb_pr_app; (prRS_app(i) - prMS_app(i)) - (prRP_app - prMP_app)];
            
            %computation of crossed observed phase ranges
            if (~isempty(pos))
                comb_pr_obs = [comb_pr_obs; lambda * ((phRS_obs - phMS_obs) - (phRP_obs - phMP_obs))];
            else
                comb_pr_obs = [comb_pr_obs; lambda * ((phRS_obs - phMS_obs) - (phRP_obs - phMP_obs) + N_kalman(i))];
            end
            
            %computation of crossed tropospheric errors
            tr = [tr; (err_tropo_RS(i) - err_tropo_MS(i)) - (err_tropo_RP - err_tropo_MP)];
            
            %computation of crossed ionospheric errors
            io = [io; -((err_iono_RS(i) - err_iono_MS(i)) - (err_iono_RP - err_iono_MP))];
        end
    end
    
    %vector of the b known term
    b = comb_pr_app;
    
    %correction of the b known term
    b = b + tr + io;
    
    %observation vector
    y0 = comb_pr_obs;
    
    %number of observations
    n = length(y0);
    
    %number of unknown parameters
    % m = 3 + nsat_pr-1;
    
    %observation noise covariance matrix
    Q = zeros(n);
    Q1 = cofactor_matrix(elR, elM, snr_R, snr_M, sat_pr, pivot);
    Q2 = Q1(r,r);
    
    Q(1:nsat_pr-1,1:nsat_pr-1) = sigmaq_cod1 * Q1;
    if (nargin == 24)
        %ambiguity estimation error is taken into account (TO BE FIXED: not properly scaled
        %with respect to input code and phase variances)
        Q(nsat_pr:end,nsat_pr:end) = (sigmaq_ph * eye(n - (nsat_pr - 1)) + lambda^2*Cee_N_kalman(p,p)) .* Q2;
    else
        Q(nsat_pr:end,nsat_pr:end) = sigmaq_ph * Q2;
    end
    % sat_slip
    % Q
    % pause
    % A_cod = A(1:nsat_pr-1,:);
    % Q_cod = Q(1:nsat_pr-1,1:nsat_pr-1);
    % y0_cod = y0(1:nsat_pr-1);
    % b_cod = b(1:nsat_pr-1);
    
    %least squares solution using only code
    % x_cod = ((A_cod'*Q_cod^-1*A_cod)^-1)*A_cod'*Q_cod^-1*(y0_cod-b_cod);
    
    %least squares solution
    x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);
    
    %test on differences between code and code+phase receiver clock delay estimation
    % while (abs(x_cod(4) - x(4)) > clock_delay_thresh)
    %
    %     %delete phase observation with maximum error variance
    %     [null_m, i] = max(diag(Q2));
    %     Q2(i,:) = [];
    %     Q2(:,i) = [];
    %
    %     Q(nsat_pr + i, :) = [];
    %     Q(:, nsat_pr + i) = [];
    %
    %     A(nsat_pr + i, :) = [];
    %     y0(nsat_pr + i) = [];
    %     b(nsat_pr + i) = [];
    %
    %     %least squares solution
    %     x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);
    % end
    
    if (nsat_slip ~= 0)
        N_stim_slip = x(4 : 4 + nsat_slip - 1);
    else
        N_stim_slip = [];
    end
    if (nsat_born ~= 0)
        N_stim_born = x(4 + nsat_slip : 4 + nsat_amb - 1);
    else
        N_stim_born = [];
    end
end
