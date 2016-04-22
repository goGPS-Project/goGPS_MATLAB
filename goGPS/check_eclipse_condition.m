function [eclipsed] = check_eclipse_condition(time, XS, SP3)

% SYNTAX:
%   [eclipsed] = check_eclipse_condition(time, XS, SP3);
%
% INPUT:
%   time = GPS time
%   XS   = satellite position (X,Y,Z)
%   SP3  = structure containing precise ephemeris data
%
% OUTPUT:
%   eclipsed = boolean value to define satellite eclipse condition (0: OK, 1: eclipsed)
%
% DESCRIPTION:
%   Check if the input satellite is under eclipse condition.

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

t_sun = SP3.t_sun;
X_sun = SP3.X_sun;

[~, q] = min(abs(t_sun - time));
X_sun = X_sun(:,q);

XS_n = norm(XS);
X_sun_n = norm(X_sun);

cosPhi = dot(XS, X_sun)/(XS_n*X_sun_n);

if (cosPhi < 0 && XS_n*sqrt(1 - cosPhi^2) < goGNSS.ELL_A_GPS)
    eclipsed = 1;
else
    eclipsed = 0;
end
