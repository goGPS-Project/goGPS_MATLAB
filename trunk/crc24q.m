function [parity] = crc24q(msg)

% SYNTAX:
%   [parity] = crc24q(msg)
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
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
% ('rtcm3torinex.c', by Dirk Stöcker, BKG Ntrip Client (BNC) Version 1.6.1
%  was used as a reference)
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

parity = 0;

for i = 1 : 8 : length(msg)
    parity = bitxor(parity, bitshift(bin2dec(msg(i:i+7)), 16));
    for j = 1 : 8
        parity = bitshift(parity, 1);
        if bitand(parity, 16777216)
            parity = bitxor(parity, 25578747);
        end
    end
end

parity = dec2bin(parity, 24);
