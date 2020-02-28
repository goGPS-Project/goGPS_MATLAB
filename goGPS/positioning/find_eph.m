function icol = find_eph(Eph, sat, time, override_dtmax)

% SYNTAX:
%   icol = find_eph(Eph, sat, time);
%
% INPUT:
%   Eph = ephemerides matrix
%   sat = satellite index
%   time = GPS time
%
% OUTPUT:
%   icol = column index for the selected satellite
%
% DESCRIPTION:
%   Selection of the column corresponding to the specified satellite
%   (with respect to the specified GPS time) in the ephemerides matrix.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       (C) Kai Borre
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

isat = find(Eph(30,:) == sat);

n = size(isat,2);
if (n == 0)
    icol = [];
    return
end
icol = isat(1);

delta = 0;

%consider BeiDou time (BDT) for BeiDou satellites
if (strcmp(char(Eph(31)),'C'))
    delta = 14;
    time = time - delta;
end

time_eph = Eph(32,icol);
dtmin = time_eph - time;
for t = isat

    time_eph = Eph(32,t);
    dt = time_eph - time;
    if (abs(dt) < abs(dtmin))
        icol = t;
        dtmin = dt;
    end
end

if nargin == 4
    dtmax = override_dtmax;
else
    %maximum interval from ephemeris reference time
    fit_interval = Eph(29,icol);
    if (fit_interval ~= 0)
        dtmax = fit_interval*3600/2;
    else
        switch (char(Eph(31,icol)))
            case 'R' %GLONASS
                dtmax = 950; %900 + 50 to account for leap seconds difference
            case 'J' %QZSS
                dtmax = 3600;
            otherwise
                dtmax = 7200;
        end
    end
end

%if (fix(abs(dtmin)) - delta > dtmax)
if (fix(abs(dtmin)) > dtmax)
    icol = [];
    return
end

%check satellite health
%the second and third conditions are temporary (QZSS and Galileo health flag is kept on for tests)
if (Eph(27,icol) ~= 0 && ~strcmp(char(Eph(31,icol)),'J') && ~strcmp(char(Eph(31,icol)),'E'))
    icol = [];
    return
end
