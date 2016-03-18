function [stidecorr] = solid_earth_tide_correction(time, XR, XS, SP3, phi)

% SYNTAX:
%   [stidecorr] = solid_earth_tide_correction(time, XR, XS, SP3, phi);
%
% INPUT:
%   time = GPS time
%   XR   = receiver position  (X,Y,Z)
%   XS   = satellite position (X,Y,Z)
%   SP3  = structure containing precise ephemeris data
%   phi  = receiver latitude (rad)
%
% OUTPUT:
%   stidecorr = solid Earth tide correction terms (along the satellite-receiver line-of-sight)
%
% DESCRIPTION:
%   Computation of the solid Earth tide displacement terms.

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

if (nargin < 5)
    phi = cart2geod(XR(1,1), XR(2,1), XR(3,1));
end

%Sun and Moon position
t_sun  = SP3.t_sun;
X_sun  = SP3.X_sun;
X_moon = SP3.X_moon;
[~, q] = min(abs(t_sun - time));
X_sun  = X_sun(:,q);
X_moon = X_moon(:,q);

%receiver geocentric position
XR_n = norm(XR);
XR_u = XR / XR_n;

%sun geocentric position
X_sun_n = norm(X_sun);
X_sun_u = X_sun / X_sun_n;

%moon geocentric position
X_moon_n = norm(X_moon);
X_moon_u = X_moon / X_moon_n;

%latitude dependence
p = (3*sin(phi)^2-1)/2;

%nominal degree 2 Love number
H2 = 0.6078 - 0.0006*p;
%nominal degree 2 Shida number
L2 = 0.0847 + 0.0002*p;

%gravitational parameters
GE = goGNSS.GM_GAL; %Earth
GS = GE*332946.0; %Sun
GM = GE*0.01230002; %Moon

%Earth equatorial radius
R = 6378136.6;

%solid Earth tide displacement
Vsun  = dot(X_sun_u,XR_u);
Vmoon = dot(X_moon_u,XR_u);
r_sun  = (GS*R^4)/(GE*X_sun_n^3) *(H2*XR_u*(1.5*Vsun^2  - 0.5) + 3*L2*Vsun *(X_sun_u  - Vsun *XR_u));
r_moon = (GM*R^4)/(GE*X_moon_n^3)*(H2*XR_u*(1.5*Vmoon^2 - 0.5) + 3*L2*Vmoon*(X_moon_u - Vmoon*XR_u));
r = r_sun + r_moon;

%displacement along the receiver-satellite line-of-sight
stidecorr = zeros(size(XS,1),1);
for s = 1 : size(XS,1)
    LOS  = XS(s,:)' - XR;
    LOSu = LOS / norm(LOS);
    stidecorr(s,1) = -dot(r,LOSu);
end
