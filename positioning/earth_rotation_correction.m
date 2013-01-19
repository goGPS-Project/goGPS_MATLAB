function [Xsat_rot] = earth_rotation_correction(traveltime, Xsat)

% SYNTAX:
%   [Xsat_rot] = earth_rotation_correction(traveltime, Xsat);
%
% INPUT:
%   traveltime = signal travel time
%   Xsat = satellite position
%
% OUTPUT:
%   Xsat_rot = corrected satellite position
%
% DESCRIPTION:
%   Correct satellite position according to Earth rotation
%   during signal travel time.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
%
% Copyright (C) Kai Borre
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

global Omegae_dot

%Xsat to column
Xsat = Xsat(:);

%find the rotation angle
omegatau = Omegae_dot * traveltime;

%build a rotation matrix
R3 = [ cos(omegatau)    sin(omegatau)   0;
      -sin(omegatau)    cos(omegatau)   0;
       0                0               1];

%apply the rotation
Xsat_rot = R3 * Xsat;
