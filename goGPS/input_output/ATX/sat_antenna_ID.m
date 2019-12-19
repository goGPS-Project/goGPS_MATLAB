function [antmod] = sat_antenna_ID(constellations)

% SYNTAX:
%   [antmod] = sat_antenna_ID(constellations);
%
% INPUT:
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   antmod = cell-array containing satellite antenna model strings
%
% DESCRIPTION:
%   Create a cell-array with satellite antenna ID strings.

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

antmod = cell(constellations.nEnabledSat,1);

if (constellations.GPS.enabled)
    for p = 1 : constellations.GPS.numSat
        antmod{constellations.GPS.indexes(p)} = ['G' num2str(constellations.GPS.PRN(p),'%02d')];
    end
end

if (constellations.GLONASS.enabled)
    for p = 1 : constellations.GLONASS.numSat
        antmod{constellations.GLONASS.indexes(p)} = ['R' num2str(constellations.GLONASS.PRN(p),'%02d')];
    end
end

if (constellations.Galileo.enabled)
    for p = 1 : constellations.Galileo.numSat
        antmod{constellations.Galileo.indexes(p)} = ['E' num2str(constellations.Galileo.PRN(p),'%02d')];
    end
end

if (constellations.BeiDou.enabled)
    for p = 1 : constellations.BeiDou.numSat
        antmod{constellations.BeiDou.indexes(p)} = ['C' num2str(constellations.BeiDou.PRN(p),'%02d')];
    end
end

if (constellations.QZSS.enabled)
    for p = 1 : constellations.QZSS.numSat
        antmod{constellations.QZSS.indexes(p)} = ['J' num2str(constellations.QZSS.PRN(p)-192,'%02d')];
    end
end
