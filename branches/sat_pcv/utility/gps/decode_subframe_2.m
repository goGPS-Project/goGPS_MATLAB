function [data] = decode_subframe_2(msg)

% SYNTAX:
%   [data] = decode_subframe_2(msg);
%
% INPUT:
%   msg = binary message containing subframe data
%
% OUTPUT:
%   data = array that contains information from subframe 2
%          1) IODE2
%          2) Crs
%          3) delta_n
%          4) M0
%          5) Cuc
%          6) e
%          7) Cus
%          8) root_A
%          9) toe
%         10) fit_int
%
% DESCRIPTION:
%   Decode GPS subframe 2 (WORDS 3-10), according to IS-GPS-200D.

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

%Subframe 2
SF2D0 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 2 - Word 0 (=Word 3 in IS-GPS-200D)
SF2D1 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 2 - Word 1 (=Word 4 in IS-GPS-200D)
SF2D2 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 2 - Word 2 (=Word 5 in IS-GPS-200D)
SF2D3 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 2 - Word 3 (=Word 6 in IS-GPS-200D)
SF2D4 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 2 - Word 4 (=Word 7 in IS-GPS-200D)
SF2D5 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 2 - Word 5 (=Word 8 in IS-GPS-200D)
SF2D6 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 2 - Word 6 (=Word 9 in IS-GPS-200D)
SF2D7 = msg(pos:pos+delta-1); pos = pos + delta;          % subframe 2 - Word 7 (=Word 10 in IS-GPS-200D)

if (little_endian)
    SF2D0 = fliplr(reshape(SF2D0,8,[])); SF2D0 = SF2D0(:)';               %
    SF2D1 = fliplr(reshape(SF2D1,8,[])); SF2D1 = SF2D1(:)';               %
    SF2D2 = fliplr(reshape(SF2D2,8,[])); SF2D2 = SF2D2(:)';               %
    SF2D3 = fliplr(reshape(SF2D3,8,[])); SF2D3 = SF2D3(:)';               % byte order inversion (little endian)
    SF2D4 = fliplr(reshape(SF2D4,8,[])); SF2D4 = SF2D4(:)';               %
    SF2D5 = fliplr(reshape(SF2D5,8,[])); SF2D5 = SF2D5(:)';               %
    SF2D6 = fliplr(reshape(SF2D6,8,[])); SF2D6 = SF2D6(:)';               %
    SF2D7 = fliplr(reshape(SF2D7,8,[])); SF2D7 = SF2D7(:)';               %
end

data(1) = fbin2dec(SF2D0(s:s+7));                                            % IODE2
data(2) = twos_complement(SF2D0(s+8:s+23)) * (2^-5);                         % Crs
data(3) = twos_complement(SF2D1(s:s+15)) * pi * (2^-43);                     % delta_n
data(4) = twos_complement([SF2D1(s+16:s+23) SF2D2(s:s+23)]) * pi * (2^-31);  % M0
data(5) = twos_complement(SF2D3(s:s+15)) * (2^-29);                          % Cuc
data(6) = fbin2dec([SF2D3(s+16:s+23) SF2D4(s:s+23)]) * (2^-33);              % e
data(7) = twos_complement(SF2D5(s:s+15)) * (2^-29);                          % Cus
data(8) = fbin2dec([SF2D5(s+16:s+23) SF2D6(s:s+23)]) * (2^-19);              % root_A
data(9) = fbin2dec(SF2D7(s:s+15)) * (2^4);                                   % toe
data(10) = fbin2dec(SF2D7(s+16));                                            % fit_int
