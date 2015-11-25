function nmeastring = NMEA_PORZA_gen(port, baudrate, protocol)

% SYNTAX:
%   nmeastring = NMEA_PORZA_gen(port, baudrate, protocol);
%
% INPUT:
%   port = receiver port to be configured (0:current port; 1:COM1; 2:COM2)
%   baudrate = baud rate to be set
%   protocol = protocol to be enabled on the specified port:
%                0: disable
%                1: NMEA 0183
%                2: RTCM-104 (differential corrections reception)
%                3: BINR
%                4: BINR2
%
% OUTPUT:
%   nmeastring = $PORZA sentence (NMEA)
%
% DESCRIPTION:
%   Returns a $PORZA sentence (NVS proprietary) in NMEA 0183 format.

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

%-----------------------------------------------------------------------------------------------
% COMPOSITION OF THE NMEA SENTENCE
%-----------------------------------------------------------------------------------------------

nmeastring = sprintf('$PORZA,%d,%d,%d',port,baudrate,protocol);

%-----------------------------------------------------------------------------------------------
% CHECKSUM COMPUTATION
%-----------------------------------------------------------------------------------------------

checksum = NMEA_checksum(nmeastring);
nmeastring = [nmeastring '*' checksum];
