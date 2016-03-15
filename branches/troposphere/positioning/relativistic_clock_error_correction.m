function [corr] = relativistic_clock_error_correction(time, Eph, XS, VS)

% SYNTAX:
%   [corr] = relativistic_clock_error_correction(time, Eph, XS, VS);
%
% INPUT:
%   time = GPS time
%   Eph = satellite ephemeris vector
%   XS = satellite position (X,Y,Z)
%   VS = satellite velocity (X,Y,Z)
%
% OUTPUT:
%   corr = relativistic clock error correction term
%
% DESCRIPTION:
%   Computation of the relativistic clock error correction term. From the
%   Interface Specification document revision E (IS-GPS-200E), page 86.

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

if (Eph(4)~=0) %if not using SP3 ephemeris
    roota = Eph(4);
    ecc   = Eph(6);
    
    Ek = ecc_anomaly(time, Eph);
    corr = -4.442807633e-10 * ecc * roota * sin(Ek);
else
    corr = -2*dot(XS,VS)/(goGNSS.V_LIGHT^2);
end
