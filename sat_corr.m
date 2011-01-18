function [Xcorr, tcorr, X, V] = sat_corr(Eph, sat, time, pr_Rsat, pos_R)

% SYNTAX:
%   [Xcorr, tcorr, X, V] = sat_corr(Eph, sat, time, pr_Rsat, pos_R);
%
% INPUT:
%   Eph = satellite ephemerides matrix
%   sat = satellite id number
%   time = GPS time
%   pr_Rsat = ROVER-SATELLITE code pseudorange
%   pos_R = rover position
%
% OUTPUT:
%   Xcorr = corrected satellite position
%   tcorr = correction due to satellite clock error
%   X = satellite position at transmission time
%   V = satellite velocity
%
% DESCRIPTION:
%   Correction of the satellite position by taking into account
%   the satellite clock error and the Earth rotation.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
%
% Partially based on scripts found in EASY suite by Kai Borre
%----------------------------------------------------------------------------------------------

global v_light

%--------------------------------------------------------------------------------------------
% CLOCK ERROR CORRECTION
%--------------------------------------------------------------------------------------------

Xcorr = [];
tcorr = [];
X = [];
V = [];

k = find_eph(Eph, sat, time);

if (isempty(k))
    return
end

af2   = Eph(2,k);
roota = Eph(4,k);
ecc   = Eph(6,k);
af0   = Eph(19,k);
af1   = Eph(20,k);
tom   = Eph(21,k);
tgd   = Eph(28,k);

tx_RAW = time - pr_Rsat / v_light;

%relativistic correction term computation
Ek = ecc_anomaly(tx_RAW, Eph(:,k));
dtr = -4.442807633e-10 * ecc * roota * sin(Ek);

dt = check_t(tx_RAW - tom);
tcorr = (af2 * dt + af1) * dt + af0 + dtr - tgd;
tx_GPS = tx_RAW - tcorr;
dt = check_t(tx_GPS - tom);
tcorr = (af2 * dt + af1) * dt + af0 + dtr - tgd;
tx_GPS = tx_RAW - tcorr;

N = 10;
for i = 1 : N
  dt = check_t(tx_GPS - tom);
  tcorr = (af2 * dt + af1) * dt + af0 + dtr - tgd;
  tx_GPS = tx_GPS - tcorr;
end

%computation of clock-corrected satellite position (and velocity)
if (nargout > 2)
    [X, V] = sat_pos(tx_GPS, Eph(:,k));
    %position and velocity with original GPS time
    %[X, V] = sat_pos(time, Eph(:,k));
else
    X = sat_pos(tx_GPS, Eph(:,k));
    %position with original GPS time
    %X = sat_pos(time, Eph(:,k));
end

%if the receiver position is not known, return
% if (nargin == 4)
%     Xcorr = X;
%     return
% end

%--------------------------------------------------------------------------------------------
% EARTH ROTATION ERROR CORRECTION
%--------------------------------------------------------------------------------------------

% rho2 = (X(1) - pos_R(1))^2 + (X(2) - pos_R(2))^2 + (X(3) - pos_R(3))^2;
% traveltime = sqrt(rho2) / v_light;
traveltime = time - tx_GPS;

%computation of rotation-corrected satellite position
Xcorr = e_r_corr(traveltime, X);

