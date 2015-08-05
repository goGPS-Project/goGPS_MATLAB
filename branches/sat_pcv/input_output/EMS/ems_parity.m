function [crc, parity] = ems_parity(msg)

% SYNTAX:
%   [crc, parity] = ems_parity(msg);
%
% INPUT:
%   msg = EMS message
%
% OUTPUT:
%   crc    = CRC value read from the input message
%   parity = parity value (computed by CRC-24Q QualComm algorithm)
%
% DESCRIPTION:
%   EMS parity check tool.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
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

h2d = hex2dec(msg');         %hexadecimal to decimal conversion
d2b = dec2bin(h2d, 4);       %4-bit binary number 
s   = reshape(d2b', 1, 256); %binary string that represents all the 256 bits of the message

%CRC data
crc = fbin2dec(s(227:250));

%portion of message used for parity checking
data = s(1:226);

%compute CRC-24Q QualComm algorithm
[parity] = crc24q(data);
parity = fbin2dec(parity);
