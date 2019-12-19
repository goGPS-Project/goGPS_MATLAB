function [data] = decode_1014(msg)

% SYNTAX:
%   [data] = decode_1014(msg)
%
% INPUT:
%   msg = binary message received from the master station
%
% OUTPUT:
%   data = cell-array that contains the 1014 packet information
%          1.1)  DF002 = message number = 1014
%          2.2)  DF059 = Network ID
%          2.3)  DF072 = Subnetwork ID
%          2.4)  DF058 = Number of auxiliary stations transmitted
%          2.5)  DF060 = Master reference station ID
%          2.6)  DF061 = Auxiliary reference station ID
%          2.7)  DF062 = Aux-Master delta latitude (degrees)
%          2.8)  DF063 = Aux-Master delta longitude (degrees)
%          2.8)  DF064 = Aux-Master delta height (mm)
%
% DESCRIPTION:
%   RTCM format 1014 message decoding.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     Sara Lucca, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

%message pointer initialization
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(8,1);

%message number = 1014
DF002 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%Network ID
DF059 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

%Subnetwork ID
DF072 = fbin2dec(msg(pos:pos+3));  pos = pos + 4;

%Number of auxiliary stations transmitted
DF058 = fbin2dec(msg(pos:pos+4));  pos = pos + 5;

%Master reference station ID
DF060 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%Auxiliary reference station ID
DF061 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%Aux-Master delta latitude
DF062 = twos_complement(msg(pos:pos+19)) * 25 * 1e-6;  pos = pos + 20;

%Aux-Master delta longitude
DF063 = twos_complement(msg(pos:pos+20)) * 25 * 1e-6;  pos = pos + 21;

%Aux-Master delta height (mm)
DF064 = twos_complement(msg(pos:pos+22))*1;

%output data save
data{1} = DF002;
data{2}(1) = DF059;
data{2}(2) = DF072;
data{2}(3) = DF058;
data{2}(4) = DF060;
data{2}(5) = DF061;
data{2}(6) = DF062;
data{2}(7) = DF063;
data{2}(8) = DF064;
