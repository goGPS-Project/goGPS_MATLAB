function [N_stim, sigmaq_N_stim] = amb_estimate_approx_SA(pos_R, sigmaq_pos_R, ...
         pr_Rsat, ph_Rsat, Eph, time, sat, phase)

% SYNTAX:
%   [N_stim, sigmaq_N_stim] = amb_estimate_approx_SA(pos_R, sigmaq_pos_R, ...
%   pr_Rsat, ph_Rsat, Eph, time, sat, phase);
%
% INPUT:
%   pos_R = ROVER assessed position (X,Y,Z)
%   sigmaq_pos_R = rounded ROVER position variance
%   pr_Rsat = ROVER-SATELLITE code-pseudorange
%   ph_Rsat = ROVER-SATELLITE phase-pseudorange
%   Eph = ephemerides matrix
%   time = GPS time
%   sat = configuration of satellites in view
%   phase = carrier L1 (phase=1), carrier L2 (phase=2)
%
% OUTPUT:
%   N_stim = linear combination of ambiguity estimate
%   sigmaq_N_stim = assessed variances of combined ambiguity
%
% DESCRIPTION:
%   Estimation of phase ambiguities (and of their error variance) by using
%   both phase observations and satellite-receiver distance, based on the
%   ROVER approximate position.

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

%variables initialization
global lambda1
global lambda2
global v_light

err_iono_RS = 0;

%cartesian to geodetic conversion of ROVER coordinates
[phiR, lamR, hR] = cart2geod(pos_R(1), pos_R(2), pos_R(3));

%radians to degrees
phiR = phiR * 180 / pi;
lamR = lamR * 180 / pi;

%loop on all used satellites
for m = 1 : size(sat,1)
    
    %new satellites position correction (clock and Earth rotation)
    [pos_S dtS]= sat_corr(Eph, sat(m), time, pr_Rsat(m));

    %computation of the satellite azimuth and elevation
    [azR, elR] = topocent(pos_R, pos_S'); %#ok<ASGLU>
    
    %computation of tropospheric errors
    err_tropo_RS = err_tropo(elR, hR);
    
    %if ionospheric parameters are available
    %if (nargin == 7)
        
        %computation of ionospheric errors
        %err_iono_RS = err_iono(iono, phiR, lamR, azR, elR, time);
    %end
    
    %ROVER,MASTER-SATELLITES pseudorange estimate
    pr_stim_Rsat(m,1) = sqrt(sum((pos_R - pos_S).^2)) - v_light*dtS + err_tropo_RS + err_iono_RS;
end

%linear combination of initial ambiguity estimate
if (phase == 1)
    N_stim = ((pr_stim_Rsat - ph_Rsat * lambda1)) / lambda1;
    sigmaq_N_stim = sum(sigmaq_pos_R) / lambda1^2;
else
    N_stim = ((pr_stim_Rsat - ph_Rsat * lambda2)) / lambda2;
    sigmaq_N_stim = sum(sigmaq_pos_R) / lambda2^2;
end
