function [corr] = sat_clock_error_correction(time, Eph)

% SYNTAX:
%   [corr] = sat_clock_error_correction(time, Eph);
%
% INPUT:
%   time = GPS time
%   Eph = satellite ephemeris vector
%
% OUTPUT:
%   corr = satellite clock correction term
%
% DESCRIPTION:
%   Computation of the satellite clock error correction term. From the
%   Interface Specification document revision E (IS-GPS-200E), page 86.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
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

%if GLONASS
if (strcmp(char(Eph(31)),'R'))

    TauN     = Eph(2);
    GammaN   = Eph(3);
    ref_toe  = Eph(32);

    dt = check_t(time - ref_toe);
    corr = -TauN + GammaN*dt;

else %if GPS/Galileo/QZSS/BeiDou

    af2 = Eph(2);
    af0 = Eph(19);
    af1 = Eph(20);
    ref_toc = Eph(33);

    %consider BeiDou time (BDT) for BeiDou satellites
    if (strcmp(char(Eph(31)),'C'))
        time = time - 14;
    end

    dt = check_t(time - ref_toc);
    corr = (af2 * dt + af1).* dt + af0;
end
