function [A, prstim_pr, prstim_ph, ddc, ddp, A0] = input_kalman_vinc(posR_app, posS, ...
          prRS_app, prMS_app, pr_Rsat, ph_Rsat, pr_Msat, ph_Msat, ...
          err_tropo_RS, err_iono_RS, err_tropo_MS, err_iono_MS, ...
          sat, pivot, phase)

% SYNTAX:
%   [A, prstim_pr, prstim_ph, ddc, ddp, A0] = input_kalman_vinc(posR_app, posS, ...
%         prRS_app, prMS_app, pr_Rsat, ph_Rsat, pr_Msat, ph_Msat, ...
%         err_tropo_RS, err_iono_RS, err_tropo_MS, err_iono_MS, ...
%         sat, pivot, phase);
%
% INPUT:
%   posR_app = receiver position (X,Y,Z)
%   posS = satellite posizion (X,Y,Z)
%   prRS_app = ROVER-SATELLITE approximate pseudorange
%   prMS_app = MASTER-SATELLITE approximate pseudorange
%   pr_Rsat = ROVER-SATELLITE code pseudorange
%   ph_Rsat = ROVER-SATELLITE phase observations
%   pr_Msat = MASTER-SATELLITE code pseudorange
%   ph_Msat = MASTER-SATELLITE phase observations
%   err_tropoRS = ROVER-SATELLITE tropospheric error
%   err_ionoRS = ROVER-SATELLITE ionospheric error
%   err_tropoMS = MASTER-SATELLITE tropospheric error
%   err_ionoMS = MASTER-SATELLITE ionospheric error
%   sat = configuration of visible satellites
%   pivot = pivot satellite
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% OUTPUT:
%   A = parameters obtained from the linearization of the observation equation,
%       projected on the constraint
%   prstim_pr = approximated code double differences
%   prstim_ph = approximated phase double differences
%   ddc = observed code double differences
%   ddp = observed phase double differences
%   A0 = parameters obtained from the linearization of the observation
%        equation (useful for DOP computation)
%
% DESCRIPTION:
%   This function computes the parameters needed to apply the Kalman filter.
%   Transition matrix that link state variables to GPS observations.
%   Constrained path.

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
global lambda1;
global lambda2;
global X_t1_t
global ax ay az s0

%number of visible satellites
nsat = length(prRS_app);

%PIVOT search
i = find(pivot == sat);

%PIVOT satellite position
posP = posS(:,i);

%ROVER-PIVOT and MASTER-PIVOT approximate pseudoranges
prRP_app = prRS_app(i);
prMP_app = prMS_app(i);

%ROVER-PIVOT and MASTER-PIVOT code observations
prRP = pr_Rsat(i);
prMP = pr_Msat(i);

%ROVER-PIVOT and MASTER-PIVOT phase observations
if (phase == 1)
    phRP = lambda1 * ph_Rsat(i);
    phMP = lambda1 * ph_Msat(i);
else
    phRP = lambda2 * ph_Rsat(i);
    phMP = lambda2 * ph_Msat(i);
end

%ROVER-PIVOT and MASTER-PIVOT tropospheric errors
err_tropo_RP = err_tropo_RS(i);
err_tropo_MP = err_tropo_MS(i);

%ROVER-PIVOT and MASTER-PIVOT ionospheric errors
err_iono_RP = err_iono_RS(i);
err_iono_MP = err_iono_MS(i);

A = [];
A0 = [];
ddc_app = [];
ddc = [];
ddp = [];
tr = [];
io = [];

%curvilinear coordinate localization
j = find((X_t1_t(1) >= s0(1:end-1)) & (X_t1_t(1) < s0(2:end)));

%computation of all the PIVOT-SATELLITE linear combinations
for i = 1 : nsat
    if (sat(i) ~= pivot)

        %construction of the transition matrix (only satellite-receiver geometry)
        A0 = [A0; (((posR_app(1) - posS(1,i)) / prRS_app(i)) - ((posR_app(1) - posP(1)) / prRP_app)) ...
                  (((posR_app(2) - posS(2,i)) / prRS_app(i)) - ((posR_app(2) - posP(2)) / prRP_app)) ...
                  (((posR_app(3) - posS(3,i)) / prRS_app(i)) - ((posR_app(3) - posP(3)) / prRP_app))];

        %construction of the transition matrix (projected on the constraint)
        A = [A; ax(j)*A0(end,1) + ay(j)*A0(end,2) + az(j)*A0(end,3)];

        %computation of the estimated code double differences
        ddc_app = [ddc_app; (prRS_app(i) - prMS_app(i)) - (prRP_app - prMP_app)];

        %computation of the observed code double differences
        ddc = [ddc; (pr_Rsat(i) - pr_Msat(i)) - (prRP - prMP)];

        %computation of the observed phase double differences
        if (phase == 1)
            ddp = [ddp; (lambda1 * ph_Rsat(i) - lambda1 * ph_Msat(i)) - (phRP - phMP)];
        else
            ddp = [ddp; (lambda2 * ph_Rsat(i) - lambda2 * ph_Msat(i)) - (phRP - phMP)];
        end
        
        %computation of tropospheric residuals
        tr = [tr; (err_tropo_RS(i) - err_tropo_MS(i)) - (err_tropo_RP - err_tropo_MP)];
        
        %computation of ionospheric residuals
        io = [io; (err_iono_RS(i) - err_iono_MS(i)) - (err_iono_RP - err_iono_MP)];
    end
end

if (phase == 1)
    prstim_pr = ddc_app + tr + io;
    prstim_ph = ddc_app + tr - io;
else
    prstim_pr = ddc_app + tr + io*(lambda2/lambda1)^2;
    prstim_ph = ddc_app + tr - io*(lambda2/lambda1)^2;
end