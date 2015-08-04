function [data] = decode_subframe_3(msg)

% SYNTAX:
%   [data] = decode_subframe_3(msg);
%
% INPUT:
%   msg = binary message containing subframe data
%
% OUTPUT:
%   data = array that contains information from subframe 3
%          1) Cic
%          2) omega0
%          3) Cis
%          4) i0
%          5) Crc
%          6) omega
%          7) omegadot
%          8) IODE3
%          9) IDOT
%
% DESCRIPTION:
%   Decode GPS subframe 3 (WORDS 3-10), according to IS-GPS-200D.

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

pos = 1;
data = [];

if (length(msg) == 256)     %32-bit WORDs (u-blox)
    little_endian = 1;
    delta = 32;
    s = 9; %starting bit
elseif (length(msg) == 192) %24-bit WORDs (SkyTraq)
    little_endian = 0;
    delta = 24;
    s = 1; %starting bit
else
    return
end

%Subframe 3
SF3D0 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 3 - Word 0 (=Word 3 in IS-GPS-200D)
SF3D1 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 3 - Word 1 (=Word 4 in IS-GPS-200D)
SF3D2 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 3 - Word 2 (=Word 5 in IS-GPS-200D)
SF3D3 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 3 - Word 3 (=Word 6 in IS-GPS-200D)
SF3D4 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 3 - Word 4 (=Word 7 in IS-GPS-200D)
SF3D5 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 3 - Word 5 (=Word 8 in IS-GPS-200D)
SF3D6 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 3 - Word 6 (=Word 9 in IS-GPS-200D)
SF3D7 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 3 - Word 7 (=Word 10 in IS-GPS-200D)

if (little_endian)
    SF3D0 = fliplr(reshape(SF3D0,8,[])); SF3D0 = SF3D0(:)';               %
    SF3D1 = fliplr(reshape(SF3D1,8,[])); SF3D1 = SF3D1(:)';               %
    SF3D2 = fliplr(reshape(SF3D2,8,[])); SF3D2 = SF3D2(:)';               %
    SF3D3 = fliplr(reshape(SF3D3,8,[])); SF3D3 = SF3D3(:)';               % byte order inversion (little endian)
    SF3D4 = fliplr(reshape(SF3D4,8,[])); SF3D4 = SF3D4(:)';               %
    SF3D5 = fliplr(reshape(SF3D5,8,[])); SF3D5 = SF3D5(:)';               %
    SF3D6 = fliplr(reshape(SF3D6,8,[])); SF3D6 = SF3D6(:)';               %
    SF3D7 = fliplr(reshape(SF3D7,8,[])); SF3D7 = SF3D7(:)';               %
end

data(1) = twos_complement(SF3D0(s:s+15)) * (2^-29);                          % Cic
data(2) = twos_complement([SF3D0(s+16:s+23) SF3D1(s:s+23)]) * pi * (2^-31);  % omega0
data(3) = twos_complement(SF3D2(s:s+15)) * (2^-29);                          % Cis
data(4) = twos_complement([SF3D2(s+16:s+23) SF3D3(s:s+23)]) * pi * (2^-31);  % i0
data(5) = twos_complement(SF3D4(s:s+15)) * (2^-5);                           % Crc
data(6) = twos_complement([SF3D4(s+16:s+23) SF3D5(s:s+23)]) * pi * (2^-31);  % omega
data(7) = twos_complement(SF3D6(s:s+23)) * pi * (2^-43);                     % omegadot
data(8) = fbin2dec(SF3D7(s:s+7));                                            % IODE3
data(9) = twos_complement(SF3D7(s+8:s+21)) * pi * (2^-43);                   % IDOT
