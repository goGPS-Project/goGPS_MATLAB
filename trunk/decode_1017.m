function [data] = decode_1017(msg)

% SYNTAX:
%   [data] = decode_1017(msg)
%
% INPUT:
%   msg = binary message received from the master station
%
% OUTPUT:
%   data = cell-array that contains the 1017 packet information
%          1.1)  DF002 = message number = 1017
%          2.2)  DF059 = Network ID
%          2.3)  DF072 = Subnetwork ID
%          2.4)  DF065 = GPS epoch time (TOW) (seconds)
%          2.5)  DF066 = GPS multiple message indicator (0 = last message, 1 = message with same message number and epoch time transmitted in sequence) 
%          2.6)  DF060 = Master reference station ID
%          2.7)  DF061 = Auxiliary reference station ID
%          2.8)  DF067 = # of GPS sats
%          3.1)  DF074 = Ambiguity status flag (0 = reserved, 1 = Correct Integer Ambiguity Level for L1 and L2, 2 = Correct Integer Ambiguity Level for L1-L2 widelane, 3 = Uncertain Integer Ambiguity Level)  
%          3.2)  DF075 = GPS non sync count (cycle slip counter)
%          3.3)  DF070 = GPS geometric carrier phase correction difference (mm)
%          3.4)  DF071 = GPS IODE
%          3.5)  DF069 = GPS ionospheric carrier phase correction (mm)
%
% DESCRIPTION:
%   RTCM format 1017 message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**, Sara Lucca*
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

%message pointer initialization
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(7,1);
data{3} = zeros(32,5);

%message number = 1017
DF002 = bin2dec(msg(pos:pos+11));  pos = pos + 12;

%Network ID
DF059 = bin2dec(msg(pos:pos+7));  pos = pos + 8;

%Subnetwork ID
DF072 = bin2dec(msg(pos:pos+3));  pos = pos + 4;

%GPS TOW (s)
DF065 = bin2dec(msg(pos:pos+22))*0.1;  pos = pos + 23;

%GPS multiple message indicator
DF066 = bin2dec(msg(pos));  pos = pos + 1;

%Master reference station ID
DF060 = bin2dec(msg(pos:pos+11));  pos = pos + 12;

%Auxiliary reference station ID
DF061 = bin2dec(msg(pos:pos+11));  pos = pos + 12;

%# of GPS Sats
DF067 = bin2dec(msg(pos:pos+3));  pos = pos + 4;

%output data save
data{1} = DF002;
data{2}(1) = DF059;
data{2}(2) = DF072;
data{2}(3) = DF065;
data{2}(4) = DF066;
data{2}(5) = DF060;
data{2}(6) = DF061;
data{2}(7) = DF067;

%data decoding for each satellite
for i = 1 : DF067

    %analyzed satellite number
    SV = bin2dec(msg(pos:pos+5));  pos = pos + 6;

    %GPS ambiguity status flag
    DF074 = bin2dec(msg(pos:pos+1));  pos = pos + 2;

    %GPS non sync count
    DF075 = bin2dec(msg(pos:pos+2));  pos = pos + 3;

    %GPS geometry carrier phase correction difference
    DF070 = twos_complement(msg(pos:pos+16))*0.5;  pos = pos + 17;

    %GPS IODE
    DF071 = bin2dec(msg(pos:pos+7));  pos = pos + 8;

    %GPS ionospheric carrier phase correction difference
    DF069 = twos_complement(msg(pos:pos+16))*0.5;  pos = pos + 17;
   
    %---------------------------------------------------------

    %output data save
    data{3}(SV,1)  = DF074;
    data{3}(SV,2)  = DF075;
    data{3}(SV,3)  = DF070;
    data{3}(SV,4)  = DF071;
    data{3}(SV,5)  = DF069;
   
end
