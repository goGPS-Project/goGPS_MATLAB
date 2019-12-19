function [corr] = relativistic_clock_error_correction(time, Eph, SP3, XS, VS)

% SYNTAX:
%   [corr] = relativistic_clock_error_correction(time, Eph, SP3, XS, VS);
%
% INPUT:
%   time = GPS time
%   Eph = satellite ephemeris vector
%   SP3 = struct with SP3 data
%   XS = satellite position (X,Y,Z)
%   VS = satellite velocity (X,Y,Z)
%
% OUTPUT:
%   corr = relativistic clock error correction term
%
% DESCRIPTION:
%   Computation of the relativistic clock error correction term. From the
%   Interface Specification document revision E (IS-GPS-200E), page 86.

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
%  Contributors:     ...
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

if (isempty(SP3)) %if not using SP3 ephemeris
    roota = Eph(4);
    ecc   = Eph(6);

    Ek = ecc_anomaly(time, Eph);
    corr = -4.442807633e-10 * ecc * roota * sin(Ek);
else
    % corr = -2 * dot(XS,VS) / (Core_Utils.V_LIGHT^2); % slower
    corr = -2 * sum(conj(XS) .* VS) / (Core_Utils.V_LIGHT ^ 2); % faster
end
