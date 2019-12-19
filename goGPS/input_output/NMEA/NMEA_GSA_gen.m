function nmeastring = NMEA_GSA_gen(sat, PDOP, HDOP, VDOP, sel, mode)

% SYNTAX:
%   nmeastring = NMEA_GSA_gen(sat, PDOP, HDOP, VDOP, sel, mode);
%
% INPUT:
%   sat = list of active satellites
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%   sel  = selection mode
%          'M' = Manual: forced to operate in 2D or 3D mode
%          'A' = Automatic
%   mode = parameter to specifiy "fix" method
%          '1': no-fix
%          '2': 2D fix
%          '3': 3D fix
%
% OUTPUT:
%   nmeastring = $GPGSA sentence (NMEA)
%
% DESCRIPTION:
%   Returns a $GPGSA sentence in NMEA 0183 format.

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

%number of satellites
nsat = size(sat,1);

%-----------------------------------------------------------------------------------------------
% COMPOSITION OF THE NMEA SENTENCE
%-----------------------------------------------------------------------------------------------

nmeastring = sprintf('$GPGSA,%c,%c', sel, mode);

for i = 1 : 12
    if (i <= nsat)
        nmeastring = [nmeastring sprintf(',%d', sat(i))];
    else
        nmeastring = [nmeastring ','];
    end
end

nmeastring = [nmeastring sprintf(',%.2f,%.2f,%.2f', PDOP, HDOP, VDOP)];

%checksum computation
checksum = NMEA_checksum(nmeastring);
nmeastring = [nmeastring '*' checksum];
