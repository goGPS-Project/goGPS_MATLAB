function [data] = decode_1019(msg, constellations)

% SYNTAX:
%   [data] = decode_1019(msg, constellations)
%
% INPUT:
%   msg = binary message received from the master station
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the 1019 packet information
%          1.1) DF002 = message number = 1019
%          2.1) DF009 = GPS satellite id
%          2.24)DF076 = GPS week number
%          2.26)DF077 = GPS SV ACCURACY
%          2.23)DF078 = GPS CODE ON L2 (00 = reserved, 01 = P code, 10 = C/A code, 11 = L2C
%          2.13)DF079 = GPS IDOT
%          2.22)DF071 = GPS IODE
%          2.21)DF081 = GPS toc
%          2.2) DF082 = GPS af2
%          2.20)DF083 = GPS af1
%          2.19)DF084 = GPS af0
%               DF085 = GPS IODC
%          2.11)DF086 = GPS Crs
%          2.5) DF087 = GPS delta-N
%          2.3) DF088 = GPS M0
%          2.8) DF089 = GPS Cuc
%          2.6) DF090 = GPS eccentricity
%          2.9) DF091 = GPS Cus
%          2.4) DF092 = GPS root A
%          2.18)DF093 = GPS toe
%          2.14)DF094 = GPS Cic
%          2.16)DF095 = GPS omega0
%          2.15)DF096 = GPS Cis
%          2.12)DF097 = GPS i0
%          2.10)DF098 = GPS Crc
%          2.7) DF099 = GPS omega
%          2.17)DF100 = GPS OMEGADOT
%          2.28)DF101 = GPS tGD
%          2.27)DF102 = SV health
%          2.25)DF103 = L2 P data flag
%          2.29)DF137 = Fit interval
%          2.30)multi-constellation satellite index (here only GPS is assumed)
%
% DESCRIPTION:
%   RTCM format 1019 message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Sara Lucca
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
data{2} = zeros(33,1);

%message number
DF002 = fbin2dec(msg(pos:pos+11)); pos = pos + 12;

%analyzed satellite number
DF009 = fbin2dec(msg(pos:pos+5)); pos = pos + 6;

%GPS week number
DF076 = fbin2dec(msg(pos:pos+9)); pos = pos + 10;

%GPS SV accuracy
DF077 = fbin2dec(msg(pos:pos+3)); pos = pos + 4;

%GPS code on L2
DF078 = fbin2dec(msg(pos:pos+1)); pos = pos + 2;

%GPS IDOT
DF079 = twos_complement(msg(pos:pos+13)) * pi * (2^-43); pos = pos + 14;

%GPS IODE
DF071 = fbin2dec(msg(pos:pos+7)); pos = pos + 8;

%GPS toc
DF081 = fbin2dec(msg(pos:pos+15))*(2^4); pos = pos + 16;

%GPS af2
DF082 = twos_complement(msg(pos:pos+7)) * (2^-55); pos = pos + 8;

%GPS af1
DF083 = twos_complement(msg(pos:pos+15)) * (2^-43); pos = pos + 16;

%af0
DF084 = twos_complement(msg(pos:pos+21)) * (2^-31); pos = pos + 22;

%GPS IODC
DF085 = fbin2dec(msg(pos:pos+9)); pos = pos + 10;

%GPS Crs
DF086 = twos_complement(msg(pos:pos+15)) * (2^-5); pos = pos + 16;

%GPS delta-n
DF087 = twos_complement(msg(pos:pos+15)) * pi * (2^-43); pos = pos + 16;

%GPS M0
DF088 = twos_complement(msg(pos:pos+31)) * pi * (2^-31); pos = pos + 32;

%GPS Cuc
DF089 = twos_complement(msg(pos:pos+15)) * (2^-29); pos = pos + 16;

%GPS ecc
DF090 = fbin2dec(msg(pos:pos+31))*(2^-33); pos = pos + 32;

%GPS Cus
DF091 = twos_complement(msg(pos:pos+15)) * (2^-29); pos = pos + 16;

%GPS rootA
DF092 = fbin2dec(msg(pos:pos+31))*(2^-19); pos = pos + 32;

%GPS toe
DF093 = fbin2dec(msg(pos:pos+15))*(2^4); pos = pos + 16;

%GPS Cic
DF094 = twos_complement(msg(pos:pos+15)) * (2^-29); pos = pos + 16;

%GPS omega0
DF095 = twos_complement(msg(pos:pos+31)) * pi * (2^-31); pos = pos + 32;

%GPS Cis
DF096 = twos_complement(msg(pos:pos+15)) * (2^-29); pos = pos + 16;

%GPS i0
DF097 = twos_complement(msg(pos:pos+31)) * pi*(2^-31); pos = pos + 32;

%GPS Crc
DF098 = twos_complement(msg(pos:pos+15)) * (2^-5); pos = pos + 16;

%GPS omega
DF099 = twos_complement(msg(pos:pos+31)) * pi * (2^-31); pos = pos + 32;

%GPS omegadot
DF100= twos_complement(msg(pos:pos+23)) * pi * (2^-43); pos = pos + 24;

%GPS tGD
DF101 = twos_complement(msg(pos:pos+7)) * (2^-31); pos = pos + 8;

%GPS SV health
DF102 = fbin2dec(msg(pos:pos+5)); pos = pos + 6;

%GPS L2 P data flag
DF103 = fbin2dec(msg(pos)); pos = pos + 1;

%GPS fit interval
DF137 = fbin2dec (msg(pos));

%force GPS system
System = int8('G');

%------------------------------------------------
%output data save (if IODC == IODE)
if (DF085 == DF071)
    data{1} = DF002;
    data{2}(1) = DF009;
    data{2}(2) = DF082;
    data{2}(3) = DF088;
    data{2}(4) = DF092;
    data{2}(5) = DF087;
    data{2}(6) = DF090;
    data{2}(7) = DF099;
    data{2}(8) = DF089;
    data{2}(9) = DF091;
    data{2}(10) = DF098;
    data{2}(11) = DF086;
    data{2}(12) = DF097;
    data{2}(13) = DF079;
    data{2}(14) = DF094;
    data{2}(15) = DF096;
    data{2}(16) = DF095;
    data{2}(17) = DF100;
    data{2}(18) = DF093;
    data{2}(19) = DF084;
    data{2}(20) = DF083;
    data{2}(21) = DF081;
    data{2}(22) = DF071;
    data{2}(23) = DF078;
    data{2}(24) = DF076;
    data{2}(25) = DF103;
    data{2}(26) = DF077;
    data{2}(27) = DF102;
    data{2}(28) = DF101;
    data{2}(29) = DF137;
    data{2}(30) = constellations.GPS.indexes(DF009);
    data{2}(31) = int8('G');
    data{2}(32) = weektow2time(DF076, DF093, 'G');
    data{2}(33) = weektow2time(DF076, DF081, 'G');
end