function [data] = decode_1020(msg, constellations)

% SYNTAX:
%   [data] = decode_1020(msg, constellations)
%
% INPUT:
%   msg = binary message received from the master station
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the 1020 packet information
%          1.1) DF002 = message number = 1020
%          2.1) DF038 = GLONASS satellite ID
%          2.12) (DF040 - 7) * 0.5625 + 1602.0 = frequency on L1
%          2.13) (DF040 - 7) * 0.4375 + 1246.0 = frequency on L2
%               DF104 = GLONASS almanac health - Cn (0=non-operability, 1=operability of n-satellite)
%               DF105 = GLONASS almanac health availability indicator (0=not available, 1=available)
%               DF106_1 = GLONASS P1 (DF106) - time interval between two adjacent values of tb (minutes)(00 = 0 min, 01 = 30 min, 10 = 45 min, 11=60 min)
%          2.4) DF107 = GLONASS tk - time referenced to the beginning of the frame within the current day (sec)
%          2.8) DF108 = GLONASS MSB of Bn word (1 = malfunctioning of given satellite)
%               DF109 = GLONASS P2 - flag of oddness (1) or evenness (0) of the value of tb
%               DF110 = GLONASS tb - time to which GLONASS navigation data are referenced (sec)
%          2.5) DF111 = GLONASS ECEF-X component of satellite velocity vector in PZ-90 datum (Km/s)
%          2.6) DF112 = GLONASS ECEF-X component of satellite coordinates in PZ-90 datum (Km)
%          2.7) DF113 = GLONASS ECEF-X component of satellite acceleration in PZ-90 datum (Km/s^2)
%          2.9) DF114 = GLONASS ECEF-Y component of satellite velocity vector in PZ-90 datum (Km/s)
%          2.10)DF115 = GLONASS ECEF-Y component of satellite coordinates in PZ-90 datum (Km)
%          2.11)DF116 = GLONASS ECEF-Y component of satellite acceleration in PZ-90 datum (Km/s^2)
%          2.14)DF117 = GLONASS ECEF-Z component of satellite velocity vector in PZ-90 datum (Km/s)
%          2.15)DF118 = GLONASS ECEF-Z component of satellite coordinates in PZ-90 datum (Km)
%          2.16)DF119 = GLONASS ECEF-Z component of satellite acceleration in PZ-90 datum (Km/s^2)
%               DF120 = GLONASS P3 - number of satellites for which almanac is transmitted within given frame (0 = 4, 1 = 5)
%          2.3) DF121 = GLONASS (gamma-n) - relative deviation of predicted satellite carrier frequency from nominal value
%               DF122 = GLONASS-M P - modification flag for the satellite transmitting given navigation signal (00 = GLONASS, 01 = GLONASS-M)
%               DF123 = GLONASS ln - satellite health flag - (1 = malfunctioning)
%          2.2) DF124 = GLONASS tau-n - correction to the satellite time relative to GLONASS system time
%               DF125 = GLONASS delta-tau -time difference between navigation RF signal transmitted in L2 sub-band and navigation RF signal transmitted in L1 sub-band
%          2.17)DF126 = GLONASS En - The age of GLONASS navigation data
%               DF127 = GLONASS-M P4 - flag of ephemeris parameters updating (1 = updated)
%               DF128 = GLONASS-M Ft predicted satellite user range accuracy at time tb
%               DF129 = GLONASS Nt calendar number of day within four-year interval starting from the 1-st of January in a leap year.
%               DF130 = GLONASS-M M - Type of GLONASS satellite (01 = GLONASS-M, 00 = GLONASS)
%               DF131 = GLONASS additional data (0 - no data)
%               DF132 = GLONASS calendar number of day within the four-year period to which ?c is referenced.
%               DF133 = GLONASS tau-c - Difference between GLONASS system time and UTC(SU) (referred to Na).
%               DF134 = GLONASS-M N4 = GLONASS year interval number starting from 1996
%               DF135 = GLONASS-M tauGPS = correction to GPS system time relative to GLONASS system time
%               DF136 = GLONASS-M ln = satellite health flag - (1 = malfunctioning)
%
% DESCRIPTION:
%   RTCM format 1020 message decoding.

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

if (~constellations.GLONASS.enabled)
    return
end

%message number = 1020
DF002 = fbin2dec(msg(pos:pos+11)); pos = pos + 12;

%GLONASS satellite ID
DF038 = fbin2dec(msg(pos:pos+5)); pos = pos + 6;

%GLONASS satellite frequency channel number
%0 --> channel number -07
%1 --> channel number -06
%...
%20 --> channel number +13
DF040 = fbin2dec(msg(pos:pos+4)); pos = pos + 5;
DF040_1 = (DF040 - 7) * 0.5625 + 1602.0;
DF040_2 = (DF040 - 7) * 0.4375 + 1246.0;

%GLONASS almanac health (Cn word)
DF104 = fbin2dec(msg(pos)); pos = pos + 1; %#ok<*NASGU>

%GLONASS almanac health availability indicator (DF105)
DF105 = fbin2dec(msg(pos)); pos = pos + 1;

%GLONASS P1
DF106_1 = (msg(pos:pos+1)); pos = pos + 2;

switch DF106_1
    %0 minutes
    case '00'
        DF106 = 0;
    %30 minutes
    case '01'
        DF106 = 30;
    %45 minutes
    case '10'
        DF106 = 45;
    %60 minutes
    case '11'
        DF106 = 60;
end

%GLONASS tk
%number of thirty-second intervals
DF107_1 = fbin2dec(msg(pos:pos)); pos = pos + 1;

%integer number of minutes
DF107_2 = fbin2dec(fliplr(msg(pos:pos+5))); pos = pos + 6;

%integer number of hours elapsed since the beginning of current day
DF107_3 = fbin2dec(fliplr(msg(pos:pos+4))); pos = pos + 5;

%time referenced to the beginning of the frame within the current day (sec)
DF107 = DF107_1 * 30 + DF107_2 * 60 + DF107_3 * 60 * 60;

%GLONASS MSB of Bn word
DF108 = fbin2dec(msg(pos)); pos = pos + 1;

%GLONASS P2
DF109 = fbin2dec(msg(pos)); pos = pos + 1;

%GLONASS tb (sec)
DF110 = fbin2dec(msg(pos:pos+6)) * 15 * 60; pos = pos + 7;

%GLONASS ECEF-X component of satellite velocity vector
%sign 0 = +, 1 = -
DF111_1 = fbin2dec(msg(pos)); pos=pos+1;
%velocity component
DF111_2 = fbin2dec(msg(pos:pos+22))*(2^-20); pos = pos + 23;
if (DF111_1 == 0)
    DF111 = DF111_2;
    else
    DF111 = DF111_2*(-1);
end

% GLONASS ECEF-X component of satellite coordinates
DF112_1 = fbin2dec(msg(pos)); pos=pos + 1;
%coordinate component
DF112_2 = fbin2dec(msg(pos:pos+25)) * (2^-11); pos = pos + 26;
if (DF112_1 == 0)
    DF112 = DF112_2;
    else
    DF112 = DF112_2*(-1);
end

%GLONASS ECEF-X component of satellite acceleration vector
DF113_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF113_2 = fbin2dec(msg(pos:pos+3)) * (2^-30); pos = pos + 4;
if (DF113_1 == 0)
    DF113 = DF113_2;
    else
    DF113 = DF113_2*(-1);
end

%GLONASS ECEF-Y component of satellite velocity vector
DF114_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF114_2 = fbin2dec(msg(pos:pos+22)) * (2^-20); pos = pos + 23;
if (DF114_1 == 0)
    DF114 = DF114_2;
    else
    DF114 = DF114_2*(-1);
end

%GLONASS ECEF-Y component of satellite coordinates
DF115_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF115_2 = fbin2dec(msg(pos:pos+25)) * (2^-11); pos = pos + 26;
if (DF115_1 == 0)
    DF115 = DF115_2;
    else
    DF115 = DF115_2*(-1);
end

%GLONASS ECEF-Y component of satellite acceleration vector
DF116_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF116_2 = fbin2dec(msg(pos:pos+3)) * (2^-30); pos = pos + 4;
if (DF116_1 == 0)
    DF116 = DF116_2;
    else
    DF116 = DF116_2*(-1);
end

%GLONASS ECEF-Z component of satellite velocity vector
DF117_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF117_2 = fbin2dec(msg(pos:pos+22)) * (2^-20); pos = pos + 23;
if (DF117_1 == 0)
    DF117 = DF117_2;
    else
    DF117 = DF117_2*(-1);
end

%GLONASS ECEF-Z component of satellite coordinates
DF118_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF118_2 = fbin2dec(msg(pos:pos+25)) * (2^-11); pos = pos + 26;
if (DF118_1 == 0)
    DF118 = DF118_2;
    else
    DF118 = DF118_2*(-1);
end

%GLONASS ECEF-Z component of satellite acceleration vector
DF119_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF119_2 = fbin2dec(msg(pos:pos+3)) * (2^-30); pos = pos + 4;
if (DF119_1 == 0)
    DF119 = DF119_2;
    else
    DF119 = DF119_2*(-1);
end

%GLONASS P3
DF120_1 = fbin2dec(msg(pos)); pos = pos + 1;
if (DF120_1 == 0)
    DF120 = 4;
else DF120 = 5;
end

%GLONASS gamma-n
DF121_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF121_2 = fbin2dec(msg(pos:pos+9)) * (2^-40); pos = pos + 10;
if (DF121_1 == 0)
    DF121 = DF121_2;
    else
    DF121 = DF121_2*(-1);
end

%GLONASS M-P
DF122 = msg(pos:pos+1); pos = pos + 2;

%GLONASS l-n
DF123 = fbin2dec(msg(pos)); pos = pos + 1;

%GLONASS tau-n
DF124_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF124_2 = fbin2dec(msg(pos:pos+20)) * (2^-30); pos = pos + 21;
if (DF124_1 == 0)
    DF124 = DF124_2;
    else
    DF124 = DF124_2*(-1);
end

%GLONASS delta-tau
DF125_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF125_2 = fbin2dec(msg(pos:pos+3)) * (2^-30); pos = pos + 4;
if (DF125_1 == 0)
    DF125 = DF125_2;
    else
    DF125 = DF125_2*(-1);
end

%GLONASS En (day)
DF126 = fbin2dec(msg(pos:pos+4)); pos = pos + 5;

%GLONASS-M P4
DF127 = fbin2dec(msg(pos)); pos = pos + 1;

%GLONASS-M Ft
DF128 = fbin2dec(msg(pos:pos+3)); pos = pos + 4;

%GLONASS-M Nt (day)
DF129 = fbin2dec(msg(pos:pos+10)); pos = pos + 11;

%GLONASS-M M
DF130 = msg(pos:pos+1); pos = pos + 2;

%GLONASS additional data
DF131 = fbin2dec(msg(pos)); pos = pos + 1;

%GLONASS N-a
DF132 = fbin2dec(msg(pos:pos+10)); pos = pos + 11;

%GLONASS tau-c (sec)
DF133_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF133_2 = fbin2dec(msg(pos:pos+30)) * (2^-31); pos = pos + 31;
if (DF133_1 == 0)
    DF133 = DF133_2;
    else
    DF133 = DF133_2*(-1);
end

%GLONASS-M N4
DF134 = fbin2dec(msg(pos:pos+4)) * 4; pos = pos + 5;

%GLONASS-M tauGPS
DF135_1 = fbin2dec(msg(pos)); pos = pos + 1;
DF135_2 = fbin2dec(msg(pos:pos+20)) * (2^-31); pos = pos + 21;
if (DF135_1 == 0)
    DF135 = DF135_2;
    else
    DF135 = DF135_2*(-1);
end

%GLONASS-M ln
DF136 = fbin2dec(msg(pos)); pos = pos + 1;

%output data save
data{1} = DF002;
data{2}(1) = DF038;
data{2}(2) = DF124;
data{2}(3) = DF121;
data{2}(4) = DF107;
data{2}(5) = DF111;
data{2}(6) = DF112;
data{2}(7) = DF113;
data{2}(8) = DF108;
data{2}(9) = DF114;
data{2}(10) = DF115;
data{2}(11) = DF116;
data{2}(12) = DF040_1;
data{2}(13) = DF040_2;
data{2}(14) = DF117;
data{2}(15) = DF118;
data{2}(16) = DF119;
data{2}(17) = DF126;
data{2}(18) = NaN; %toe not available, taken care of by the caller
data{2}(19) = 0;
data{2}(20) = 0;
data{2}(21) = 0;
data{2}(22) = 0;
data{2}(23) = 0;
data{2}(24) = 0; %weekno not available, taken care of by the caller
data{2}(25) = 0;
data{2}(26) = 0;
data{2}(27) = 0; %sv health not available
data{2}(28) = 0;
data{2}(29) = 0;
data{2}(30) = constellations.GLONASS.indexes(DF038);
data{2}(31) = int8('R');
data{2}(32) = 0; %continuous toe taken care of by the caller
data{2}(33) = 0;
