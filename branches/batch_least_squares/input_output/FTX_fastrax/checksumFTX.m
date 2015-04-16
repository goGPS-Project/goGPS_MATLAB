function checksum = checksumFTX(bit_msg, len)

% SYNTAX:
%   checksum = checksumFTX(bit_msg, len);
%
% INPUT:
%   bit_msg = binary message
%   len = payload length
%
% DESCRIPTION:
%   Fastrax checksum computation.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Ivan Reguzzoni
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
%--------------------------------------------------------------------------

% initialization
checksum = 0;
pos = 1;

for i = 0 : (len-1)
    tmp = (checksum + 1) * (fbin2dec([bit_msg(pos+8:pos+15) bit_msg(pos:pos+7)])+i);
    sh_tmp = bitshift(tmp,-16);
    tmp = bitxor(tmp, sh_tmp);
    checksum = bitxor(checksum, tmp);
    % optimization
    if (checksum > 65535)
        checksum = checksum - (floor(checksum/65536))*65536;
%         checksum = dec2bin(checksum,16);
%         checksum = fbin2dec(checksum(end-15:end));
    end
    pos = pos+16;
end

checksum = dec2bin(checksum,16);
checksum = [checksum(9:16) checksum(1:8)];

% ---- Test ----
% dec2hex(fbin2dec(checksum))
