function [Xsat_rot] = earth_rotation_correction(traveltime, Xsat, Omegae_dot)

% SYNTAX:
%   [Xsat_rot] = earth_rotation_correction(traveltime, Xsat, Omegae_dot);
%
% INPUT:
%   traveltime = signal travel time
%   Xsat = satellite position at transmission time (time_tx) in ECEF (X,Y,Z)
%   Omegae_dot = angular velocity of the Earth rotation [rad/s]
%
% OUTPUT:
%   Xsat_rot = corrected satellite position
%
% DESCRIPTION:
%   Correct satellite position according to Earth rotation
%   during signal travel time.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) Kai Borre
%  Written by:       (C) Kai Borre
%  Contributors:     Mirko Reguzzoni, Eugenio Realini, 2009
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

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
