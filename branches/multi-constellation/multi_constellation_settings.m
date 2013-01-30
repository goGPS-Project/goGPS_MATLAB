% SYNTAX:
%   multi_constellation_settings;
%
% DESCRIPTION:
%   Multi-constellation settings and initialization.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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

constellations.GPS     = struct('numSat', 32, 'enabled', 1, 'indexes', []);
constellations.GLONASS = struct('numSat', 30, 'enabled', 0, 'indexes', []);
constellations.Galileo = struct('numSat', 30, 'enabled', 0, 'indexes', []);
constellations.BeiDou  = struct('numSat', 30, 'enabled', 0, 'indexes', []);
constellations.QZSS    = struct('numSat',  4, 'enabled', 0, 'indexes', []);

nSatTot = 0; %total number of satellites used given the enabled constellations
q = 0;                          %counter for enabled constellations

systems = fieldnames(constellations); 
for i = 1:numel(systems)
    if(constellations.(systems{i}).enabled)
        nSatTot = nSatTot + constellations.(systems{i}).numSat;
        q = q + 1;
        if (q == 1)
            indexes_tmp = [1 : constellations.(systems{i}).numSat];
        else
            indexes_tmp = [indexes_tmp(end) + 1 : indexes_tmp(end) + constellations.(systems{i}).numSat];
        end
        constellations.(systems{i}).indexes = indexes_tmp;
    end
end

constellations.nEnabledSat = nSatTot;
