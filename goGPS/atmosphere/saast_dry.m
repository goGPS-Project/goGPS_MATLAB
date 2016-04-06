function [ZHD] = saast_dry(P, H, lat)

% SYNTAX:
%   [ZHD] = saast_dry(P, H, lat);
%
% INPUT:
%   P = atmospheric pressure [hPa]
%   H = orthometric height [m]
%   lat = latitude [deg]
%
% OUTPUT:
%   ZHD = Zenith Hydrostatic Delay
%
% DESCRIPTION:
%   Zenith Hydrostatic Delay (ZHD) computation by Saastamoinen model.

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

%ZHD (Saastamoinen model)
ZHD = 0.0022768*P*(1+0.00266*cosd(2*lat)+0.00000028*H);
