function [pos_S_SP3, vel_S, dtr] = interpolate_SP3_coor(time, SP3_time, SP3_coor)

% SYNTAX:
%   [pos_S_SP3, vel_S, dtr] = interpolate_SP3_coor(time, SP3_time, SP3_coor);
%
% INPUT:
%   time = interpolation time (GPS time)
%   SP3_time = precise ephemeris epochs (GPS time)
%   SP3_coor = satellite coordinates  [m]
%   sat = satellite PRN
%   tcorr = interpolated SP3 clock correction
%
% OUTPUT:
%   pos_S_SP3 = interpolated satellite coordinates
%   vel_S = satellite velocity
%   dtr = relativistic correction term
%
% DESCRIPTION:
%   SP3 (precise ephemeris) coordinates 1-second interpolation by Lagrange
%   polynomials. Satellite velocity computation. Relativistic correction.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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

global v_light

%number of seconds in a quarter of an hour
quarter_sec = 900;

%find the SP3 epoch closest to the interpolation time
[~, p] = min(abs(SP3_time - time));

b = SP3_time(p) - time;

pos_S_SP3 = zeros(3,1);

%extract the SP3 coordinates
SP3_X = []; SP3_Y = []; SP3_Z = [];
for i = -3 : +3
    SP3_X = [SP3_X SP3_coor(1,p+i)];
    SP3_Y = [SP3_Y SP3_coor(2,p+i)];
    SP3_Z = [SP3_Z SP3_coor(3,p+i)];
end

x = 1 : 7;

% %Lagrange interpolation (coordinates)
% u = 4 - b/quarter_sec;
% pos_S_SP3(1,1) = LagrangeInter(x, SP3_X, u);
% pos_S_SP3(2,1) = LagrangeInter(x, SP3_Y, u);
% pos_S_SP3(3,1) = LagrangeInter(x, SP3_Z, u);
% 
% %interpolate on the previous and following seconds for computing velocity
% u = [4 - b/quarter_sec - 1/quarter_sec, 4 - b/quarter_sec + 1/quarter_sec];
% pos_S_v(1,1:2) = LagrangeInter(x, SP3_X, u);
% pos_S_v(2,1:2) = LagrangeInter(x, SP3_Y, u);
% pos_S_v(3,1:2) = LagrangeInter(x, SP3_Z, u);

%Lagrange interpolation (coordinates)
s = 1/quarter_sec;
t = 4 - b/quarter_sec;
u = t - s : s : t + s;

LI_SP3_X = LagrangeInter(x, SP3_X, u);
LI_SP3_Y = LagrangeInter(x, SP3_Y, u);
LI_SP3_Z = LagrangeInter(x, SP3_Z, u);

pos_S_SP3(1,1) = LI_SP3_X(2);
pos_S_SP3(2,1) = LI_SP3_Y(2);
pos_S_SP3(3,1) = LI_SP3_Z(2);

pos_S_v(1,1:2) = [LI_SP3_X(1) LI_SP3_X(3)];
pos_S_v(2,1:2) = [LI_SP3_Y(1) LI_SP3_Y(3)];
pos_S_v(3,1:2) = [LI_SP3_Z(1) LI_SP3_Z(3)];

%compute velocity
vel_S = (pos_S_v(:,2) - pos_S_v(:,1)) / 2;

%compute the relativistic correction term for the satellite clock
dtr = -2*dot(pos_S_SP3,vel_S)/(v_light^2);

% %Lagrange interpolation (coordinates)
% u = 5 - b/quarter_sec - dtr/quarter_sec;
% pos_S_SP3(1,1) = LagrangeInter(x, SP3_X, u);
% pos_S_SP3(2,1) = LagrangeInter(x, SP3_Y, u);
% pos_S_SP3(3,1) = LagrangeInter(x, SP3_Z, u);
