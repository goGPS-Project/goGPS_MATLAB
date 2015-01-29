function [Xsat_rot] = earth_rotation_correction(traveltime, Xsat, Omegae_dot)

% SYNTAX:
%   [Xsat_rot] = earth_rotation_correction(traveltime, Xsat, Omegae_dot);
%
% INPUT:
%   traveltime = signal travel time
%   Xsat = satellite position
%   Omegae_dot = angular velocity of the Earth rotation [rad/s]
%
% OUTPUT:
%   Xsat_rot = corrected satellite position
%
% DESCRIPTION:
%   Correct satellite position according to Earth rotation
%   during signal travel time.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) Kai Borre
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

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
