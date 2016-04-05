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

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.2 beta
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
        antmod{constellations.QZSS.indexes(p)} = ['J' num2str(constellations.QZSS.PRN(p),'%02d')];
    end
end
