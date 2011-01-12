function [data] = decode_subframe_4(msg)

% SYNTAX:
%   [data] = decode_subframe_4(msg);
%
% INPUT:
%   msg = binary message containing subframe data
%
% OUTPUT:
%   data = array that contains information from subframe 4
%          1) a0
%          2) a1
%          3) a2
%          4) a3
%          5) b0
%          6) b1
%          7) b2
%          8) b3
%          9) leap seconds
%
% DESCRIPTION:
%   Decode GPS subframe 4 (WORDS 3-10), according to IS-GPS-200D.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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

%Subframe 4
SF4D0 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 0 (=Word 3 in IS-GPS-200D)
SF4D1 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 1 (=Word 4 in IS-GPS-200D)
SF4D2 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 2 (=Word 5 in IS-GPS-200D)
SF4D3 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 3 (=Word 6 in IS-GPS-200D)
SF4D4 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 4 (=Word 7 in IS-GPS-200D)
SF4D5 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 5 (=Word 8 in IS-GPS-200D)
SF4D6 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 6 (=Word 9 in IS-GPS-200D)
SF4D7 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 7 (=Word 10 in IS-GPS-200D)

SF4D0 = fliplr(reshape(SF4D0,8,[])); SF4D0 = SF4D0(:)';               %
SF4D1 = fliplr(reshape(SF4D1,8,[])); SF4D1 = SF4D1(:)';               %
SF4D2 = fliplr(reshape(SF4D2,8,[])); SF4D2 = SF4D2(:)';               % byte order inversion (little endian)
SF4D6 = fliplr(reshape(SF4D6,8,[])); SF4D6 = SF4D6(:)';               %

data(1) = twos_complement(SF4D0(17:24)) * (2^-30);   % a0
data(2) = twos_complement(SF4D0(25:32)) * (2^-27);   % a1
data(3) = twos_complement(SF4D1(9:16))  * (2^-24);   % a2
data(4) = twos_complement(SF4D1(17:24)) * (2^-24);   % a3
data(5) = twos_complement(SF4D1(25:32)) * (2^11);    % b0
data(6) = twos_complement(SF4D2(9:16))  * (2^14);    % b1
data(7) = twos_complement(SF4D2(17:24)) * (2^16);    % b2
data(8) = twos_complement(SF4D2(25:32)) * (2^16);    % b3
data(9) = twos_complement(SF4D6(9:16));              % leap seconds
