function [pos_S, vel_S] = interpolate_SP3_coord(time, SP3, sat)

% SYNTAX:
%   [pos_S, vel_S] = interpolate_SP3_coord(time, SP3, sat);
%
% INPUT:
%   time = interpolation time (GPS time, continuous since 6-1-1980)
%   SP3  = structure containing precise ephemeris data
%   sat = satellite PRN
%
% OUTPUT:
%   pos_S = interpolated satellite coordinates
%   vel_S = satellite velocity
%
% DESCRIPTION:
%   SP3 (precise ephemeris) coordinates 1-second interpolation by Lagrange
%   polynomials. Satellite velocity computation. Relativistic correction.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

SP3_time  = SP3.time;
SP3_coord = SP3.coord(:, sat, :);

%degree of interpolation polynomial (Lagrange)
n = 10;

%number of seconds in a quarter of an hour
quarter_sec = 900;

%find the SP3 epoch closest to the interpolation time
[~, p] = min(abs(SP3_time - time));

b = SP3_time(p) - time;

pos_S = zeros(3,1);

%extract the SP3 coordinates
SP3_X = []; SP3_Y = []; SP3_Z = [];
for i = -n/2 : n/2
    SP3_X = [SP3_X SP3_coord(1,p+i)];
    SP3_Y = [SP3_Y SP3_coord(2,p+i)];
    SP3_Z = [SP3_Z SP3_coord(3,p+i)];
end

x = 1 : n+1;

% %Lagrange interpolation (coordinates)
% u = 4 - b/quarter_sec;
% pos_S(1,1) = LagrangeInter(x, SP3_X, u);
% pos_S(2,1) = LagrangeInter(x, SP3_Y, u);
% pos_S(3,1) = LagrangeInter(x, SP3_Z, u);
% 
% %interpolate on the previous and following seconds for computing velocity
% u = [4 - b/quarter_sec - 1/quarter_sec, 4 - b/quarter_sec + 1/quarter_sec];
% pos_S_v(1,1:2) = LagrangeInter(x, SP3_X, u);
% pos_S_v(2,1:2) = LagrangeInter(x, SP3_Y, u);
% pos_S_v(3,1:2) = LagrangeInter(x, SP3_Z, u);

%Lagrange interpolation (coordinates)
s = 1/quarter_sec;
t = n/2+1 - b/quarter_sec;
u = t - s : s : t + s;

LI_SP3_X = LagrangeInter(x, SP3_X, u);
LI_SP3_Y = LagrangeInter(x, SP3_Y, u);
LI_SP3_Z = LagrangeInter(x, SP3_Z, u);

pos_S(1,1) = LI_SP3_X(2);
pos_S(2,1) = LI_SP3_Y(2);
pos_S(3,1) = LI_SP3_Z(2);

pos_S_v(1,1:2) = [LI_SP3_X(1) LI_SP3_X(3)];
pos_S_v(2,1:2) = [LI_SP3_Y(1) LI_SP3_Y(3)];
pos_S_v(3,1:2) = [LI_SP3_Z(1) LI_SP3_Z(3)];

%compute velocity
vel_S = (pos_S_v(:,2) - pos_S_v(:,1)) / 2;

% if (nargout > 2)
%     %compute the relativistic correction term for the satellite clock
%     dtrel = -2*dot(pos_S,vel_S)/(v_light^2);
% end

% %Lagrange interpolation (coordinates)
% u = 5 - b/quarter_sec - dtr/quarter_sec;
% pos_S(1,1) = LagrangeInter(x, SP3_X, u);
% pos_S(2,1) = LagrangeInter(x, SP3_Y, u);
% pos_S(3,1) = LagrangeInter(x, SP3_Z, u);
