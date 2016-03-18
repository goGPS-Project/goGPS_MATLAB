function [phwindup] = phase_windup_correction(time, XR, XS, SP3, phwindup)

% SYNTAX:
%   [phwindup] = phase_windup_correction(time, XR, XS, SP3, phwindup);
%
% INPUT:
%   time = GPS time
%   XR   = receiver position  (X,Y,Z)
%   XS   = satellite position (X,Y,Z)
%   SP3  = structure containing precise ephemeris data
%   phwindup = phase wind-up (previous value)
%
% OUTPUT:
%   phwindup = phase wind-up (updated value)
%
% DESCRIPTION:
%   Computation of the phase wind-up terms.

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

%east (a) and north (b) local unit vectors
[phi, lam] = cart2geod(XR(1,1), XR(2,1), XR(3,1));
a = [-sin(lam); cos(lam); 0];
b = [-sin(phi)*cos(lam); -sin(phi)*sin(lam); cos(phi)];
    
for s = 1 : size(XS,1)
    %satellite-fixed local unit vectors
    [i, j, k] = satellite_fixed_frame(time, XS(s,:)', SP3);

    %receiver and satellites effective dipole vectors
    Dr = a - k*dot(k,a) + cross(k,b);
    Ds = i - k*dot(k,i) - cross(k,j);
    
    %phase wind-up computation
    psi = dot(k, cross(Ds, Dr));
    arg = dot(Ds,Dr)/(norm(Ds)*norm(Dr));
    if (arg < -1)
        arg = -1;
    elseif (arg > 1)
        arg = 1;
    end
    dPhi = sign(psi)*acos(arg)/(2*pi);
    if (isempty(phwindup) || phwindup(s,1) == 0)
        N = 0;
    else
        N = round(phwindup(s,1) - dPhi);
    end
    phwindup(s,1) = dPhi + N;
end
