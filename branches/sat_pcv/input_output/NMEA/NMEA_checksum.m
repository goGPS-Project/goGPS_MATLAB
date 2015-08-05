function checksum = NMEA_checksum(nmeastring)

% SYNTAX:
%   checksum = NMEA_checksum(nmeastring);
%
% INPUT:
%   nmeastring = NMEA sentence without checksum
%
% OUTPUT:
%   checksum = checksum to be appended to NMEA sentence
%
% DESCRIPTION:
%   Checksum computation as required by NMEA 0183 format.

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

checksum = uint8(0);
nmeastring_d = uint8(nmeastring);

for count = 2:length(nmeastring)                     % checksum computation ignores $ at start
    checksum = bitxor(checksum,nmeastring_d(count)); % checksum computation
end

% convert checksum to hex value
checksum = dec2hex(checksum, 2);
