function [ZWD] = saast_wet(T, H, h)

% SYNTAX:
%   [ZWD] = saast_wet(T, H, h);
%
% INPUT:
%   T = air temperature
%   H = humidity
%   h = orthometric height
%
% OUTPUT:
%   ZWD = Zenith Wet Delay
%
% DESCRIPTION:
%   Zenith Wet Delay (ZWD) computation by Saastamoinen model.

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

% Convert C -> K
T = T + 273.15;

%height correction
H = H * exp(-0.0006396*h);

% Convert humidity
H = H./100;

c = -37.2465 + 0.213166*T - 2.56908 * (10^-4) * (T.^2);
e = H .* exp(c);

%ZWD (Saastamoinen model)
ZWD = 0.0022768*(((1255/T)+0.05)*e);
