function [data] = decode_1029(msg)

% SYNTAX:
%   [data] = decode_1029(msg)
%
% INPUT:
%   msg = binary message received from the master station
%
% OUTPUT:
%   data = cell-array that contains the 1029 packet information
%          1.1)  DF002 = message number = 1029
%          2.1) DF003 = Reference Station ID
%          2.2) DF051 = Modified Julian Day (MJD) Number
%          2.3) DF052 = Seconds of Day (UTC)
%          2.4) DF138 = Number of Characters to Follow - This represents the number of fully formed Unicode characters in the message text. It is not necessarily the number of bytes that are needed to represent the characters as UTF-8. Note that for some messages it may not be possible to utilize the full range of this field, e.g. where many characters require 3 or 4 byte representations and together will exceed 255 code units.
%          2.5) DF139 = Number of UTF-8 Code Units (N) - The length of the message is limited by this field, or possibly by DF+1.
%          2.6) DF140 = UTF-8 Character Code Units
%
% DESCRIPTION:
%   RTCM format 1029 message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Ivan Reguzzoni
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

%   RTCM format 1029 message decoding.

%message pointer initialization
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(25,1);

%Message Number DF002 uint12 12 1029
DF002 = fbin2dec(msg(pos:pos+11)); pos = pos + 12;

%Reference Station ID DF003 uint12 12
DF003 = fbin2dec(msg(pos:pos+11)); pos = pos + 12;

%Modified Julian Day (MJD) Number DF051 uint16 16
DF051 = fbin2dec(msg(pos:pos+15)); pos = pos + 16;

%Seconds of Day (UTC) DF052 uint17 17
DF052 = fbin2dec(msg(pos:pos+16)); pos = pos + 17;

%Number of Characters to Follow DF138 uint7 7 This represents the number of fully formed Unicode characters in the message text. It is not necessarily the number of bytes that are needed to represent the characters as UTF-8. Note that for some messages it may not be possible to utilize the full range of this field, e.g. where many characters require 3 or 4 byte representations and together will exceed 255 code units.
DF138 = fbin2dec(msg(pos:pos+6)); pos = pos + 7;

%Number of UTF-8 Code Units (N) DF139 uint8 8 The length of the message is limited by this field, or possibly by DF+1.
DF139 = fbin2dec(msg(pos:pos+7)); pos = pos + 8;

%UTF-8 Character Code Units DF140 utf8(N) 8*N
DF140 = fbin2dec(msg(pos:pos+8*DF139-1)); pos = pos + 8*DF139;

%output data save
data{1} = DF002;
data{2}(1) = DF003;
data{2}(2) = DF051;
data{2}(3) = DF052;
data{2}(4) = DF138;
data{2}(5) = DF139;
data{2}(6) = DF140;
