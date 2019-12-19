function [latpp, lonpp, fpp] = iono_pierce_point(latR, lonR, azS, elS)

% SYNTAX:
%   [latpp, lonpp, fpp] = iono_pierce_point(latR, lonR, azS, elS);
%
% INPUT:
%   latR = receiver position latitude  [rad]
%   lonR = receiver position longitude [rad]
%   azS  = satellite azimuth [rad]
%   elS  = satellite elevation [rad]
%
% OUTPUT:
%   latpp = ionosphere piercing point latitude  [rad]
%   lonpp = ionosphere piercing point longitude [rad]
%   fpp   = slant factor (mapping function)
%
% DESCRIPTION:
%   Computation of the ionosphere piercing point (IPP).

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Antonio Herrera Olmo, 2012
%  Contributors:     Eugenio Realini, 2013
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

R  = 6378.1363; %Earth radius [km]
hI = 350;       %ionosphere thin shell height [km]

k = (R/(R+hI))*cos(elS);
psipp = (pi/2) - elS - asin(k);

%set azimuth from -180 to 180
azS = mod((azS+pi),2*pi)-pi;

%latitude of the ionosphere piercing point
latpp = asin(sin(latR)*cos(psipp) + cos(latR)*sin(psipp)*cos(azS));

%longitude of the ionosphere piercing point
if ((latpp >  70*pi/180) & (tan(psipp)*cos(azS)      > tan((pi/2) - latR))) | ...
   ((latpp < -70*pi/180) & (tan(psipp)*cos(azS + pi) > tan((pi/2) + latR)))

    lonpp = lonR + pi - asin(sin(psipp)*sin(azS)/cos(latpp));
else
    lonpp = lonR + asin(sin(psipp)*sin(azS)/cos(latpp));
end

%slant (obliquity) factor
fpp = (1-(k)^2)^(-1/2);
