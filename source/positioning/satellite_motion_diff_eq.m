function [pos_dot, vel_dot] = satellite_motion_diff_eq(pos, vel, acc, a, GM, J2, Omegae_dot)

% SYNTAX:
%   [pos_dot, vel_dot] = satellite_motion_diff_eq(pos, vel, acc, a, GM, J2, Omegae_dot);
%
% INPUT:
%   pos = satellite position (XYZ)
%   vel = satellite velocity
%   acc = acceleration due to lunar-solar gravitational perturbation
%   a   = ellipsoid semi-major axis [m]
%   GM  = gravitational constant (mass of Earth) [m^3/s^2]
%   J2  = second zonal harmonic of the geopotential
%   Omegae_dot = angular velocity of the Earth rotation [rad/s]
%                (if provided, an Earth-fixed system will be used)
%
% OUTPUT:
%   pos_dot = differential position
%   vel_dot = differential velocity
%
% DESCRIPTION:
%   Computation of the differential equations of perturbed orbital motion.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

if (nargin < 6)
    %do the computation in an inertial system
    Omegae_dot = 0;
end

%differential position
pos_dot = vel;

%renaming variables for better readability
%
%position
x = pos(1);
y = pos(2);
z = pos(3);
%velocity
vx = vel(1);
vy = vel(2);
%acceleration (i.e. perturbation)
ax = acc(1);
ay = acc(2);
az = acc(3);

%parameters
r = sqrt(x^2 + y^2 + z^2);
g = -GM/r^3;
h = J2*1.5*(a/r)^2;
k = 5*z^2/r^2;

%differential velocity
vel_dot = zeros(size(vel));
vel_dot(1) = g*x*(1 - h*(k - 1)) + ax + Omegae_dot^2*x + 2*Omegae_dot*vy;
vel_dot(2) = g*y*(1 - h*(k - 1)) + ay + Omegae_dot^2*y - 2*Omegae_dot*vx;
vel_dot(3) = g*z*(1 - h*(k - 3)) + az;
