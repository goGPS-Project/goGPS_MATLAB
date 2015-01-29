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
%   Computation of the satellite position (X,Y,Z) and velocity by means
%   of its ephemerides.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Antonio Herrera Olmo, 2012
% Adapted by Eugenio Realini, 2013
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

R  = 6378.1363; %Earth radius [km]
hI = 350;       %ionosphere thin shell height [km]

k = (R/(R+hI))*cos(elS);
phipp = (pi/2) - elS - asin(k);

%latitude of the ionosphere piercing point
latpp = asin(sin(latR)*cos(phipp) + cos(latR)*sin(phipp)*cos(azS));

%longitude of the ionosphere piercing point
if ((latpp >  70*pi/180) & (tan(phipp)*cos(azS)      > tan((pi/2) - latR))) | ...
   ((latpp < -70*pi/180) & (tan(phipp)*cos(azS + pi) > tan((pi/2) + latR)))
    
    lonpp = lonR + pi - asin((sin(phipp)*sin(azS))/(cos(latpp)));
else
    lonpp = lonR + asin((sin(phipp)*sin(azS))/(cos(latpp)));
end

%slant (obliquity) factor
fpp = (1-(k)^2)^(-1/2);
