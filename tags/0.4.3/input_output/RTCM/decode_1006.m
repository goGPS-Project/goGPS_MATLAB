function [data] = decode_1006(msg)

% SYNTAX:
%   [data] = decode_1006(msg)
%
% INPUT:
%   msg = binary message received from the master station
%
% OUTPUT:
%   data = cell-array that contains the 1006 packet information
%          1.1)  DF002 = message number = 1006
%          2.1)  DF003 = reference station id
%          2.2)  DF021 = ITRF reference frame year (default = 0)
%          2.3)  DF022 = is GPS service supported? YES=1, NO=0
%          2.4)  DF023 = is GLONASS service supported? YES=1, NO=0
%          2.5)  DF024 = is Galileo service supported? YES=1, NO=0
%          2.6)  DF141 = reference station type (real=0, virtual=1)
%          2.7)  DF142 = observations acquired at the same instant YES=1, NO=0
%          2.8)  (DF025 * 0.0001) = antenna X coordinate in meters (ITRF)
%          2.9)  (DF026 * 0.0001) = antenna Y coordinate in meters (ITRF)
%          2.10) (DF027 * 0.0001) = antenna Z coordinate in meters (ITRF)
%          2.11) (DF028 * 0.0001) = antenna height (with respect to the marker)
%
% DESCRIPTION:
%   RTCM format 1006 message decoding.

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
data{2} = zeros(11,1);

%message number = 1006
DF002 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%reference station id
DF003 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%ITRF reference frame year
DF021 = fbin2dec(msg(pos:pos+5));  pos = pos + 6;

%GPS service index (YES=1, NO=0)
DF022 = fbin2dec(msg(pos));  pos = pos + 1;

%GLONASS service index (YES=1, NO=0)
DF023 = fbin2dec(msg(pos));  pos = pos + 1;

%Galileo service index (YES=1, NO=0)
DF024 = fbin2dec(msg(pos));  pos = pos + 1;

%master station type (real=0, virtual=1)
DF141 = fbin2dec(msg(pos));  pos = pos + 1;

%antenna X coordinate
DF025 = twos_complement(msg(pos:pos+37));  pos = pos + 38;

%single oscillator index (YES=1, NO=0)
DF142 = fbin2dec(msg(pos));  pos = pos + 1;

%skip reserved field
pos = pos + 1;

%antenna Y coordinate
DF026 = twos_complement(msg(pos:pos+37));  pos = pos + 38;

%skip reserved field
pos = pos + 2;

%antenna Z coordinate
DF027 = twos_complement(msg(pos:pos+37));  pos = pos + 38;

%antenna height
DF028 = fbin2dec(msg(pos:pos+15));

%--------------------------------------------------------------------------------------------

%output data save
data{1} = DF002;
data{2}(1)  = DF003;
data{2}(2)  = DF021;
data{2}(3)  = DF022;
data{2}(4)  = DF023;
data{2}(5)  = DF024;
data{2}(6)  = DF141;
data{2}(7)  = DF142;
data{2}(8)  = (DF025 * 0.0001);
data{2}(9)  = (DF026 * 0.0001);
data{2}(10) = (DF027 * 0.0001);
data{2}(11) = (DF028 * 0.0001);