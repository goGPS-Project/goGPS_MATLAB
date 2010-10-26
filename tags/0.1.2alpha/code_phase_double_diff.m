function [pos_R, cov_pos_R, N_stim, cov_N_stim, PDOP, HDOP, VDOP] = code_phase_double_diff ...
         (pos_R_app, pr_R, ph_R, snr_R, pos_M, pr_M, ph_M, snr_M, time, sat, pivot, Eph, phase, iono)

% SYNTAX:
%   [pos_R, cov_pos_R, N_stim, cov_N_stim, PDOP, HDOP, VDOP] = code_phase_double_diff ...
%   (pos_R_app, pr_R, ph_R, snr_R, pos_M, pr_M, ph_M, snr_M, time, sat, pivot, Eph, phase, iono);
%
% INPUT:
%   pos_R_app = ROVER position (X,Y,Z)
%   pr_R = ROVER-SATELLITE code pseudorange
%   ph_R = ROVER-SATELLITE phase observation
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   pos_M = MASTER position (X,Y,Z)
%   pr_M = MASTER-SATELLITE code pseudorange
%   ph_M = MASTER-SATELLITE phase observation
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   time = GPS time
%   sat = visible satellite configuration
%   pivot = pivot satellite
%   Eph = ephemerides matrix
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%   iono = ionospheric parameters
%
% OUTPUT:
%   pos_R = estimated position (X,Y,Z)
%   cov_pos_R = covariance matrix of estimation errors (rover position)
%   N_stim = linear combination of ambiguity estimate
%   cov_N_stim = covariance matrix of estimation errors (combined ambiguity values)
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision

%
% DESCRIPTION:
%   Least squares solution using code and phase double differences.
%   Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
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

%rover position coordinates X Y Z
X_R = pos_R_app(1);
Y_R = pos_R_app(2);
Z_R = pos_R_app(3);

%conversion from cartesian to geodetic coordinates
[phiR, lamR, hR] = cart2geod(X_R, Y_R, Z_R);

%master position coordinates X Y Z
X_M = pos_M(1);
Y_M = pos_M(2);
Z_M = pos_M(3);

%conversion from cartesian to geodetic coordinates
[phiM, lamM, hM] = cart2geod(X_M, Y_M, Z_M);

%radians to degrees
phiR = phiR * 180 / pi;
lamR = lamR * 180 / pi;
phiM = phiM * 180 / pi;
lamM = lamM * 180 / pi;

%number of visible satellites
nsat = size(sat,1);

%PIVOT satellite index
i = find(pivot == sat);

%PIVOT position (with clock error and Earth rotation corrections)
posP = sat_corr(Eph, sat(i), time, pr_R(i), pos_R_app);

%computation of ROVER-PIVOT and MASTER-PIVOT approximated pseudoranges
prRP_app = sqrt(sum((pos_R_app - posP).^2));
prMP_app = sqrt(sum((pos_M - posP).^2));

%ROVER-PIVOT and MASTER-PIVOT observed code measurements
prRP_obs = pr_R(i);
prMP_obs = pr_M(i);

%ROVER-PIVOT and MASTER-PIVOT observed phase measurements
phRP_obs = ph_R(i);
phMP_obs = ph_M(i);

%ROVER-satellite and MASTER-satellite elevation initialization
elR = zeros(nsat,1);
elM = zeros(nsat,1);

%computation of PIVOT azimuth and elevation
[azR, elR(i)] = topocent(pos_R_app, posP');
[azM, elM(i)] = topocent(pos_M, posP');

%ROVER-PIVOT and MASTER-PIVOT tropospheric error computation
err_tropo_RP = err_tropo(elR(i), hR);
err_tropo_MP = err_tropo(elM(i), hM);

%if ionospheric parameters are available
if (nargin == 14)

   %ROVER-PIVOT and MASTER-PIVOT ionospheric error computation
   err_iono_RP = err_iono(iono, phiR, lamR, azR, elR(i), time);
   err_iono_MP = err_iono(iono, phiM, lamM, azM, elM(i), time);
end

for j = 1 : nsat
    if (sat(j) ~= pivot)

        %satellite position (with clock error and Earth rotation corrections)
        posS(j,:) = sat_corr(Eph, sat(j), time, pr_R(j), pos_R_app);

        %computation of the satellite azimuth and elevation
        [azR, elR(j)] = topocent(pos_R_app, posS(j,:));
        [azM, elM(j)] = topocent(pos_M, posS(j,:));

        %computation of ROVER-SATELLITE and MASTER-SATELLITE approximated pseudoranges
        prRS_app(j) = sqrt(sum((pos_R_app - posS(j,:)').^2));
        prMS_app(j) = sqrt(sum((pos_M - posS(j,:)').^2));
        
        %computation of tropospheric errors
        err_tropo_RS(j) = err_tropo(elR(j), hR);
        err_tropo_MS(j) = err_tropo(elM(j), hM);

        %if ionospheric parameters are available
        if (nargin == 14)

            %computation of ionospheric errors
            err_iono_RS(j) = err_iono(iono, phiR, lamR, azR, elR(j), time);
            err_iono_MS(j) = err_iono(iono, phiM, lamM, azM, elM(j), time);
        end

    end
end

A = [];
tr = [];
io = [];
comb_pr_app = [];
comb_pr_obs = [];

%computation of all linear combinations between PIVOT and other satellites (CODE)
for j = 1 : nsat
    if (sat(j) ~= pivot)

        %observed code measurement
        prRS_obs = pr_R(j);
        prMS_obs = pr_M(j);

        %design matrix computation
        A = [A; (((pos_R_app(1) - posS(j,1)) / prRS_app(j)) - ((pos_R_app(1) - posP(1)) / prRP_app)) ...
                (((pos_R_app(2) - posS(j,2)) / prRS_app(j)) - ((pos_R_app(2) - posP(2)) / prRP_app)) ...
                (((pos_R_app(3) - posS(j,3)) / prRS_app(j)) - ((pos_R_app(3) - posP(3)) / prRP_app)) ...
                zeros(1, nsat-1)];

        %computation of crossed approximate pseudoranges
        comb_pr_app = [comb_pr_app; (prRS_app(j) - prMS_app(j)) - (prRP_app - prMP_app)];

        %computation of crossed observed code pseudoranges
        comb_pr_obs = [comb_pr_obs; (prRS_obs - prMS_obs) - (prRP_obs - prMP_obs)];
        
        %computation of crossed tropospheric errors
        tr = [tr; (err_tropo_RS(j) - err_tropo_MS(j)) - (err_tropo_RP - err_tropo_MP)];

        %if ionospheric parameters are available
        if (nargin == 14)

            %computation of crossed ionospheric errors
            io = [io; (err_iono_RS(j) - err_iono_MS(j)) - (err_iono_RP - err_iono_MP)];
        end

    end
end

%computation of all linear combinations between PIVOT and other satellites (PHASE)

%index for identifying the correct slot for ambiguities (excluding PIVOT)
k = 1;

for j = 1 : nsat
    if (sat(j) ~= pivot)

        %observed phase measurement
        phRS_obs = ph_R(j);
        phMS_obs = ph_M(j);

        %ambiguity vector in design matrix (lambda position)
        N_row = zeros(1, nsat-1);
        N_row(k) = -lambda;

        %design matrix computation
        A = [A; (((pos_R_app(1) - posS(j,1)) / prRS_app(j)) - ((pos_R_app(1) - posP(1)) / prRP_app)) ...
                (((pos_R_app(2) - posS(j,2)) / prRS_app(j)) - ((pos_R_app(2) - posP(2)) / prRP_app)) ...
                (((pos_R_app(3) - posS(j,3)) / prRS_app(j)) - ((pos_R_app(3) - posP(3)) / prRP_app)) ...
                N_row];

        %computation of crossed approximated pseudoranges
        comb_pr_app = [comb_pr_app; (prRS_app(j) - prMS_app(j)) - (prRP_app - prMP_app)];

        %computation of crossed observed phase ranges
        comb_pr_obs = [comb_pr_obs; ((phRS_obs - phMS_obs) - (phRP_obs - phMP_obs)) * lambda];
        
        %computation of crossed tropospheric errors
        tr = [tr; (err_tropo_RS(j) - err_tropo_MS(j)) - (err_tropo_RP - err_tropo_MP)];

        %if ionospheric parameters are available
        if (nargin == 14)

            %computation of crossed ionospheric errors
            io = [io; -((err_iono_RS(j) - err_iono_MS(j)) - (err_iono_RP - err_iono_MP))];
        end
        k = k + 1;
    end
end

%vector of the b known term
b = comb_pr_app;

%correction of the b known term
b = b + tr;
if (nargin == 14)
   b = b + io;
end

%observation vector
y0 = comb_pr_obs;

%number of observations
n = length(y0);

%number of unknown parameters
m = 3 + (nsat - 1);

%observation noise covariance matrix
Q = zeros(n, n);
Q1 = cofactor_matrix(elR, elM, snr_R, snr_M, sat, pivot);
Q(1:n/2,1:n/2) = sigmaq_cod1 * Q1;
Q(n/2+1:end,n/2+1:end) = sigmaq_ph * Q1;

%parameter vector
xR = [pos_R_app; zeros(nsat-1,1)];

%least squares solution
x = ((A'*Q^-1*A)^-1)*A'*Q^-1*(y0-b);

try
    xR = xR + x;
catch
    %DEBUG
    fprintf('Size of parameter vector: %d\n', size(xR,1));
    fprintf('Size of observation vector : %d\n', size(y0,1));
    fprintf('Number of satellites (with pivot): %d\n', nsat);
    xR = xR + x;
end

%estimation of the variance of the observation error
y_stim = A*x + b;
v_stim = y0 - y_stim;
sigma0q_stim = (v_stim' * Q^-1 * v_stim) / (n-m);

%estimated rover position
pos_R = xR(1:3);

%estimated combined ambiguity values (without PIVOT)
N_stim_nopivot = xR(4:end);

%add a zero at PIVOT position
N_stim = zeros(nsat,1);
N_stim(1:i-1) = N_stim_nopivot(1:i-1);
N_stim(i+1:end) = N_stim_nopivot(i:end);

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma0q_stim * ((A'*Q^-1*A)^-1);

    %rover position covariance matrix
    cov_pos_R = Cxx(1:3,1:3);

    %combined ambiguity covariance matrix
    cov_N_stim_nopivot = Cxx(4:end,4:end);

    %add one line and one column (zeros) at PIVOT position
    cov_N_stim = zeros(nsat,nsat);
    cov_N_stim(1:i-1,1:i-1) = cov_N_stim_nopivot(1:i-1,1:i-1);
    cov_N_stim(i+1:end,i+1:end) = cov_N_stim_nopivot(i:end,i:end);

else
    cov_pos_R = [];

    cov_N_stim = [];
end

%DOP computation
if (nargout > 4)
    cov_XYZ = (A(1:n/2,1:3)'*A(1:n/2,1:3))^-1;
    cov_ENU = global2localCov(cov_XYZ, xR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end