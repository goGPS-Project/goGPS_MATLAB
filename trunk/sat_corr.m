function [X, tcorr] = sat_corr(Eph, sat, time, pr_Rsat, pos_R)

% SYNTAX:
%   [X, tcorr] = sat_corr(Eph, sat, time, pr_Rsat, pos_R);
%
% INPUT:
%   Eph = satellite ephemerides matrix
%   sat = satellite id number
%   time = GPS time
%   pr_Rsat = ROVER-SATELLITE code pseudorange
%   pos_R = rover position
%
% OUTPUT:
%   X = corrected satellite position
%   tcorr = correction due to satellite clock error
%
% DESCRIPTION:
%   Correction of the satellite position by taking into account
%   the satellite clock error and the Earth rotation.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) Kai Borre
% Kai Borre and C.C. Goad 11-24-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------

global v_light

%--------------------------------------------------------------------------------------------
% CLOCK ERROR CORRECTION
%--------------------------------------------------------------------------------------------

col_Eph = find_eph(Eph, sat, time);
k = col_Eph;
tx_RAW = time - pr_Rsat / v_light;
t0m = Eph(21,k);
dt = check_t(tx_RAW - t0m);
tcorr = (Eph(2,k) * dt + Eph(20,k)) * dt + Eph(19,k);
tx_GPS = tx_RAW - tcorr;
dt = check_t(tx_GPS - t0m);
tcorr = (Eph(2,k) * dt + Eph(20,k)) * dt + Eph(19,k);
tx_GPS = tx_RAW - tcorr;

%N = 10;
%for i = 1 : N
%   dt = check_t(tx_GPS - t0m);
%   tcorr = (Eph(2,k) * dt + Eph(20,k)) * dt + Eph(19,k);
%   tx_GPS = tx_GPS - tcorr;
%end

%position with original GPS time
%X = sat_pos(time, Eph(:,k));

%computation of clock-corrected satellite position
X = sat_pos(tx_GPS, Eph(:,k));

%if the receiver position is not known, return
if (nargin == 4)
    return
end

%--------------------------------------------------------------------------------------------
% EARTH ROTATION ERROR CORRECTION
%--------------------------------------------------------------------------------------------

rho2 = (X(1) - pos_R(1))^2 + (X(2) - pos_R(2))^2 + (X(3) - pos_R(3))^2;
traveltime = sqrt(rho2) / v_light;

%computation of rotation-corrected satellite position
X = e_r_corr(traveltime, X);
