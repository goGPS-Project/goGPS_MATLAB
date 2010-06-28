function [X_sat_rot] = e_r_corr(traveltime, X_sat)

% SYNTAX:
%   [X_sat_rot] = e_r_corr(traveltime, X_sat);
%
% INPUT:
%   traveltime = signal travel time
%   X_sat = satellite position
%
% OUTPUT:
%   X_sat_rot = corrected satellite position
%
% DESCRIPTION:
%   Correct satellite position according to Earth rotation
%   during signal travel time.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.1 alpha
%
% Copyright (C) Kai Borre
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------

global Omegae_dot

%find the rotation angle
omegatau   = Omegae_dot * traveltime;

%build a rotation matrix
R3 = [ cos(omegatau)    sin(omegatau)   0;
      -sin(omegatau)    cos(omegatau)   0;
       0                0               1];

%apply the rotation
X_sat_rot = R3 * X_sat;