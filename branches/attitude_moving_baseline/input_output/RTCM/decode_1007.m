function [data] = decode_1007(msg)

% SYNTAX:
%   [data] = decode_1007(msg)
%
% INPUT:
%   msg = binary message received from the master station
%
% OUTPUT:
%   data = cell-array that contains the 1007 packet information
%          1.1) DF002 = message number = 1007
%          2.1) DF003 = reference station id
%          2.2) DF031 = antenna setup
%          3.1) DF030 = antenna description (number of characters in DF029)
%
% DESCRIPTION:
%   RTCM format 1007 message decoding.

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

%message pointer initialization
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(2,1);
data{3} = '';

%message number = 1007
DF002 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%reference station id
DF003 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%char number for the antenna description
DF029 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

%antenna description
DF030 = char(fbin2dec(reshape(msg(pos:pos+8*DF029-1),8,[])'))';  pos = pos + 8*DF029;

%antenna setup
DF031 = fbin2dec(msg(pos:pos+7));

%--------------------------------------------------------------------------------------------

%output data save
data{1} = DF002;
data{2}(1) = DF003;
data{2}(2) = DF031;
data{3} = DF030;