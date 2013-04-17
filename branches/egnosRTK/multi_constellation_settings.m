function [constellations] = multi_constellation_settings(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)

% SYNTAX:
%   [constellations] = multi_constellation_settings(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
%
% INPUT:
%   GPS_flag = boolean flag for enabling/disabling GPS usage
%   GLO_flag = boolean flag for enabling/disabling GLONASS usage
%   GAL_flag = boolean flag for enabling/disabling Galileo usage
%   BDS_flag = boolean flag for enabling/disabling BeiDou usage
%   QZS_flag = boolean flag for enabling/disabling QZSS usage
%   SBS_flag = boolean flag for enabling/disabling SBAS usage (for ranging)
%
% OUTPUT:
%   constellations = struct with multi-constellation settings
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

GPS_PRN = [1:32];
GLO_PRN = [1:24];
GAL_PRN = [1:30];
BDS_PRN = [1:37];
QZS_PRN = [193:196];
SBS_PRN = 0; %SBAS ranging not supported yet

constellations.GPS     = struct('numSat', numel(GPS_PRN), 'enabled', GPS_flag, 'indexes', 0, 'PRN', GPS_PRN);
constellations.GLONASS = struct('numSat', numel(GLO_PRN), 'enabled', GLO_flag, 'indexes', 0, 'PRN', GLO_PRN);
constellations.Galileo = struct('numSat', numel(GAL_PRN), 'enabled', GAL_flag, 'indexes', 0, 'PRN', GAL_PRN);
constellations.BeiDou  = struct('numSat', numel(BDS_PRN), 'enabled', BDS_flag, 'indexes', 0, 'PRN', BDS_PRN);
constellations.QZSS    = struct('numSat', numel(QZS_PRN), 'enabled', QZS_flag, 'indexes', 0, 'PRN', QZS_PRN);
constellations.SBAS    = struct('numSat', numel(SBS_PRN), 'enabled', 0,        'indexes', 0, 'PRN', SBS_PRN); %SBAS ranging not supported yet

nSatTot = 0; %total number of satellites used given the enabled constellations
q = 0;       %counter for enabled constellations

systems = fieldnames(constellations);
constellations.indexes = [];
constellations.PRN = [];
for i = 1 : numel(systems)
    if(constellations.(systems{i}).enabled)
        nSatTot = nSatTot + constellations.(systems{i}).numSat;
        q = q + 1;
        if (q == 1)
            indexes_tmp = [1 : constellations.(systems{i}).numSat];
        else
            indexes_tmp = [indexes_tmp(end) + 1 : indexes_tmp(end) + constellations.(systems{i}).numSat];
        end
        constellations.(systems{i}).indexes = indexes_tmp;
        constellations.indexes = [constellations.indexes, indexes_tmp];
        constellations.PRN = [constellations.PRN, constellations.(systems{i}).PRN];
    end
end

constellations.nEnabledSat = nSatTot;
