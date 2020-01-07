function [eclipsed] = check_eclipse_condition(time, XS, SP3, sat, p_rate)

% SYNTAX:
%   [eclipsed] = check_eclipse_condition(time, XS, SP3, sat, p_rate);
%
% INPUT:
%   time     = GPS time
%   XS       = satellite position (X,Y,Z)
%   SP3      = structure containing precise ephemeris data
%   p_rate   = processing interval [s]
%   sat      = satellite PRN
%
% OUTPUT:
%   eclipsed = boolean value to define satellite eclipse condition (0: OK, 1: eclipsed)
%
% DESCRIPTION:
%   Check if the input satellite is under eclipse condition.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     Andrea Gatti...
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

eclipsed = 0;

t_sun = SP3.t_sun;
X_sun = SP3.X_sun;

%[~, q] = min(abs(t_sun - time));
% speed improvement of the above line
% supposing t_sun regularly sampled
q = round((time - t_sun(1)) / p_rate) + 1;

X_sun = X_sun(:,q);

%satellite geocentric position
XS_n = norm(XS);
XS_u = XS / XS_n;

%sun geocentric position
X_sun_n = norm(X_sun);
X_sun_u = X_sun / X_sun_n;

%satellite-sun angle
%cosPhi = dot(XS_u, X_sun_u);
% speed improvement of the above line
cosPhi = sum(conj(XS_u').*X_sun_u);


%threshold to detect noon/midnight maneuvers
if (~isempty(strfind(SP3.satType{sat},'BLOCK IIA')))
    t = 4.9*pi/180; % maximum yaw rate of 0.098 deg/sec (Kouba, 2009)
elseif (~isempty(strfind(SP3.satType{sat},'BLOCK IIR')))
    t = 2.6*pi/180; % maximum yaw rate of 0.2 deg/sec (Kouba, 2009)
elseif (~isempty(strfind(SP3.satType{sat},'BLOCK IIF')))
    t = 4.35*pi/180; % maximum yaw rate of 0.11 deg/sec (Dilssner, 2010)
else
    t = 0; %ignore noon/midnight maneuvers for other constellations (TBD)
end

%shadow crossing affects only BLOCK IIA satellites
shadowCrossing = cosPhi < 0 && XS_n*sqrt(1 - cosPhi^2) < goGNSS.ELL_A_GPS;
if (shadowCrossing && ~isempty(strfind(SP3.satType{sat},'BLOCK IIA')))
    eclipsed = 1;
end

%noon/midnight maneuvers affect all satellites
noonMidnightTurn = acos(abs(cosPhi)) < t;
if (noonMidnightTurn)
    eclipsed = 3;
end
