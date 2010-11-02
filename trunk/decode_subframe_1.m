function [data] = decode_subframe_1(msg)

% SYNTAX:
%   [data] = decode_subframe_1(msg);
%
% INPUT:
%   msg = binary message containing subframe data
%
% OUTPUT:
%   data = array that contains information from subframe 1
%          1) weekno
%          2) code_on_L2
%          3) svaccur
%          4) svhealth
%          5) IODC
%          6) L2flag
%          7) tgd
%          8) toc
%          9) af2
%         10) af1
%         11) af0
%
% DESCRIPTION:
%   Decode GPS subframe 1 (WORDS 3-10), according to IS-GPS-200D.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
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

%Subframe 1
SF1D0 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 0 (=Word 3 in IS-GPS-200D)
SF1D1 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 1 (=Word 4 in IS-GPS-200D)
SF1D2 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 2 (=Word 5 in IS-GPS-200D)
SF1D3 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 3 (=Word 6 in IS-GPS-200D)
SF1D4 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 4 (=Word 7 in IS-GPS-200D)
SF1D5 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 5 (=Word 8 in IS-GPS-200D)
SF1D6 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 6 (=Word 9 in IS-GPS-200D)
SF1D7 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 7 (=Word 10 in IS-GPS-200D)

SF1D0 = fliplr(reshape(SF1D0,8,[])); SF1D0 = SF1D0(:)';               % byte order inversion (little endian)
SF1D1 = fliplr(reshape(SF1D1,8,[])); SF1D1 = SF1D1(:)';               %
SF1D4 = fliplr(reshape(SF1D4,8,[])); SF1D4 = SF1D4(:)';               %
SF1D5 = fliplr(reshape(SF1D5,8,[])); SF1D5 = SF1D5(:)';               %
SF1D6 = fliplr(reshape(SF1D6,8,[])); SF1D6 = SF1D6(:)';               %
SF1D7 = fliplr(reshape(SF1D7,8,[])); SF1D7 = SF1D7(:)';               %

data(1) = fbin2dec(SF1D0(9:18));                          % weekno
data(2) = fbin2dec(SF1D0(19:20));                         % code_on_L2
data(3) = fbin2dec(SF1D0(21:24));                         % svaccur
data(4) = fbin2dec(SF1D0(25:30));                         % svhealth
IODC_2MSBs = SF1D0(31:32);                           % IODC_2MSBs
IODC_8LSBs = SF1D5(9:16);                            % IODC_8LSBs
data(5) = fbin2dec([IODC_2MSBs IODC_8LSBs]);              % IODC
data(6) = fbin2dec(SF1D1(9));                             % L2flag
data(7) = twos_complement(SF1D4(25:32)) * (2^-31);   % tgd
data(8) = fbin2dec(SF1D5(17:32)) * (2^4);                 % toc
data(9) = twos_complement(SF1D6(9:16)) * (2^-55);    % af2
data(10) = twos_complement(SF1D6(17:32)) * (2^-43);  % af1
data(11) = twos_complement(SF1D7(9:30)) * (2^-31);   % af0
