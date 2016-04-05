function nmeastring = NMEA_PGGPK_gen(sat, KPDOP, KHDOP, KVDOP, mode)

% SYNTAX:
%   nmeastring = NMEA_PGGPK_gen(sat, KPDOP, KHDOP, KVDOP, mode);
%
% INPUT:
%   sat = list of active satellites
%   KPDOP = Kalman filter position dilution of precision
%   KHDOP = Kalman filter horizontal dilution of precision
%   KVDOP = Kalman filter vertical dilution of precision
%   mode = parameter to specifiy "fix" method
%          'S': satellite positioning enabled
%          'D': Kalman filter dynamics only
%
% OUTPUT:
%   nmeastring = $PGGPK goGPS proprietary sentence (NMEA)
%
% DESCRIPTION:
%   Returns a $PGGPK goGPS proprietary sentence in NMEA 0183 format.

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

%number of satellites
nsat = size(sat,1);

%-----------------------------------------------------------------------------------------------
% COMPOSITION OF THE NMEA SENTENCE
%-----------------------------------------------------------------------------------------------

nmeastring = sprintf('$PGGPK,%c', mode); 

for i = 1 : 12
    if (i <= nsat)
        nmeastring = [nmeastring sprintf(',%d', sat(i))];
    else
        nmeastring = [nmeastring ','];
    end
end

nmeastring = [nmeastring sprintf(',%.2f,%.2f,%.2f', KPDOP, KHDOP, KVDOP)];

%checksum computation
checksum = NMEA_checksum(nmeastring);
nmeastring = [nmeastring '*' checksum];
