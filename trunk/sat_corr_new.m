function [Xcorr, tcorr, X, V, tx_GPS] = sat_corr(Eph, sat, time, pseudorange)

% SYNTAX:
%   [Xcorr, tcorr, X, V, tx_GPS] = sat_corr(Eph, sat, time, pseudorange);
%
% INPUT:
%   Eph = satellite ephemerides matrix
%   sat = satellite id number
%   time = GPS time
%   pseudorange = code pseudorange
%
% OUTPUT:
%   Xcorr = corrected satellite position
%   tcorr = correction due to satellite clock error
%   X = satellite position at transmission time
%   V = satellite velocity
%   tx_GPS = clock-corrected transmission time
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
tgd   = Eph(28,k); %This correction term is only for the benefit of "single-frequency" (L1 P(Y) or L2 P(Y)) users

tx_RAW = time - pseudorange / v_light;

%relativistic correction term computation
Ek = ecc_anomaly(tx_RAW, Eph(:,k));
dtr = -4.442807633e-10 * ecc * roota * sin(Ek);

dt = check_t(tx_RAW - tom);
tcorr = (af2 * dt + af1) * dt + af0 + dtr - tgd;
tx_GPS = tx_RAW - tcorr;
dt = check_t(tx_GPS - tom);
tcorr = (af2 * dt + af1) * dt + af0 + dtr - tgd;
tx_GPS = tx_RAW - tcorr;

% N = 10;
% for i = 1 : N
%   dt = check_t(tx_GPS - tom);
%   tcorr = (af2 * dt + af1) * dt + af0 + dtr - tgd;
%   tx_GPS = tx_GPS - tcorr;
% end

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

traveltime = time - tx_GPS;

%computation of rotation-corrected satellite position
Xcorr = e_r_corr(traveltime, X);

