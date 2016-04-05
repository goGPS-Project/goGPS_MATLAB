function [parity] = crc24q(msg)

% SYNTAX:
%   [parity] = crc24q(msg);
%
% INPUT:
%   msg = binary message
%
% OUTPUT:
%   parity = crc parity (24 bits)
%
% DESCRIPTION:
%   Applies CRC-24Q QualComm algorithm.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% ('rtcm3torinex.c', by Dirk Stoecker, BKG Ntrip Client (BNC) Version 1.6.1
%  was used as a reference)
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

parity = uint32(0);

%check the length of the input string, in case make it splittable byte-wise
remainder = rem(length(msg),8);
if (remainder ~= 0)
    fill = char(ones(1,8-remainder)*48); %fill string of zeroes
    msg = [fill msg];
end

Nbits = length(msg);

%pre-allocate to increase speed
Nbytes = Nbits / 8;
bytes = cell(1,Nbytes);
k = 1;
for j = 1 : 8 : Nbits
    bytes{k} = msg(j:j+7);
    k = k + 1;
end
%call 'fbin2dec' and 'bitshift' only once (to optimize speed)
bytes = bitshift(fbin2dec(bytes), 16);
bytes = uint32(bytes);
for i = 1 : Nbytes
    parity = bitxor(parity, bytes(i));
    for j = 1 : 8
        parity = bitshift(parity, 1);
        if bitand(parity, 16777216)
            parity = bitxor(parity, 25578747);
        end
    end
end

parity = dec2bin(parity, 24);
