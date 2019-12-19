function [stidecorr] = solid_earth_tide_correction(time, XR, XS, SP3, p_rate, phiC, lam)

% SYNTAX:
%   [stidecorr] = solid_earth_tide_correction(time, XR, XS, SP3, phiC, lam);
%
% INPUT:
%   time = GPS time
%   XR   = receiver position  (X,Y,Z)
%   XS   = satellite position (X,Y,Z)
%   SP3  = structure containing precise ephemeris data
%   p_rate   = processing interval [s]
%   phiC = receiver geocentric latitude (rad)
%   lam  = receiver longitude (rad)
%
% OUTPUT:
%   stidecorr = solid Earth tide correction terms (along the satellite-receiver line-of-sight)
%
% DESCRIPTION:
%   Computation of the solid Earth tide displacement terms.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
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

if (nargin < 6)
    [~, lam, ~, phiC] = cart2geod(XR(1,1), XR(2,1), XR(3,1));
end
%north (b) and radial (c) local unit vectors
b = [-sin(phiC)*cos(lam); -sin(phiC)*sin(lam); cos(phiC)];
c = [+cos(phiC)*cos(lam); +cos(phiC)*sin(lam); sin(phiC)];

%Sun and Moon position
t_sun  = SP3.t_sun;
X_sun  = SP3.X_sun;
X_moon = SP3.X_moon;
%[~, q] = min(abs(t_sun - time));
% speed improvement of the above line
% supposing t_sun regularly sampled
q = round((time - t_sun(1)) / p_rate) + 1;
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
p = (3*sin(phiC)^2-1)/2;

%gravitational parameters
GE = goGNSS.GM_GAL; %Earth
GS = GE*332946.0; %Sun
GM = GE*0.01230002; %Moon

%Earth equatorial radius
R = 6378136.6;

%nominal degree 2 Love number
H2 = 0.6078 - 0.0006*p;
%nominal degree 2 Shida number
L2 = 0.0847 + 0.0002*p;

%solid Earth tide displacement (degree 2)
Vsun  = sum(conj(X_sun_u) .* XR_u);
Vmoon = sum(conj(X_moon_u) .* XR_u);
r_sun2  = (GS*R^4)/(GE*X_sun_n^3) *(H2*XR_u*(1.5*Vsun^2  - 0.5) + 3*L2*Vsun *(X_sun_u  - Vsun *XR_u));
r_moon2 = (GM*R^4)/(GE*X_moon_n^3)*(H2*XR_u*(1.5*Vmoon^2 - 0.5) + 3*L2*Vmoon*(X_moon_u - Vmoon*XR_u));
r = r_sun2 + r_moon2;

%nominal degree 3 Love number
H3 = 0.292;
%nominal degree 3 Shida number
L3 = 0.015;

%solid Earth tide displacement (degree 3)
r_sun3  = (GS*R^5)/(GE*X_sun_n^4) *(H3*XR_u*(2.5*Vsun^3  - 1.5*Vsun)  +   L3*(7.5*Vsun^2  - 1.5)*(X_sun_u  - Vsun *XR_u));
r_moon3 = (GM*R^5)/(GE*X_moon_n^4)*(H3*XR_u*(2.5*Vmoon^3 - 1.5*Vmoon) +   L3*(7.5*Vmoon^2 - 1.5)*(X_moon_u - Vmoon*XR_u));
r = r + r_sun3 + r_moon3;

%from "conventional tide free" to "mean tide"
radial = (-0.1206 + 0.0001*p)*p;
north  = (-0.0252 + 0.0001*p)*sin(2*phiC);
r = r + radial*c + north*b;

%displacement along the receiver-satellite line-of-sight
stidecorr = zeros(size(XS,1),1);
for s = 1 : size(XS,1)
    LOS  = XR - XS(s,:)';
    LOSu = LOS / norm(LOS);
    %stidecorr(s,1) = dot(r,LOSu);
    stidecorr(s,1) = sum(conj(r).*LOSu);
end
