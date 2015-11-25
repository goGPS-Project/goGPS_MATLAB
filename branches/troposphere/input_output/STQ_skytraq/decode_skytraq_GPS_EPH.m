function [data] = decode_skytraq_GPS_EPH(msg, constellations)

% SYNTAX:
%   [data] = decode_skytraq_GPS_EPH(msg, constellations);
%
% INPUT:
%   msg = message transmitted by the SkyTraq receiver
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the GPS Ephemeris data information
%          1.1) message class-id (GPS-EPH)
%          2.1) GPS satellite id
%          2.2) GPS af2
%          2.3) GPS M0
%          2.4) GPS root A
%          2.5) GPS delta-N
%          2.6) GPS eccentricity
%          2.7) GPS omega
%          2.8) GPS Cuc
%          2.9) GPS Cus
%          2.10)GPS Crc
%          2.11)GPS Crs
%          2.12)GPS i0
%          2.13)GPS IDOT
%          2.14)GPS Cic
%          2.15)GPS Cis
%          2.16)GPS omega0
%          2.17)GPS omegadot
%          2.18)GPS toe
%          2.19)GPS af0
%          2.20)GPS af1
%          2.21)GPS toc
%          2.22)GPS IODE
%          2.23)GPS codes;
%          2.24)GPS weekno;
%          2.25)GPS L2flag;
%          2.26)GPS svaccur;
%          2.27)GPS svhealth;
%          2.28)GPS tgd;
%          2.29)GPS fit_int;
%          2.30)multi-constellation satellite index (here only GPS is assumed)
%
% DESCRIPTION:
%   GPS Ephemeris data binary message decoding.

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

if (nargin < 2 || isempty(constellations))
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

%first message initial index
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(33,1);

%output data save
data{1} = 'GPS_EPH';

%satellite ID decoding (2 byte)
PRN = fbin2dec(msg(pos:pos+15));  pos = pos + 16;

%reserved (1 byte)
pos = pos + 8;

%skip WORD 2
pos = pos + 24;

%Subframe 1
[subframe_1_data]  = decode_subframe_1(msg(pos:pos+191)); pos = pos + 192;

weekno     = subframe_1_data(1);
code_on_L2 = subframe_1_data(2);
svaccur    = subframe_1_data(3);
svhealth   = subframe_1_data(4);
IODC       = subframe_1_data(5);
L2flag     = subframe_1_data(6);
tgd        = subframe_1_data(7);
toc        = subframe_1_data(8);
af2        = subframe_1_data(9);
af1        = subframe_1_data(10);
af0        = subframe_1_data(11);

%reserved (1 byte)
pos = pos + 8;

%skip WORD 2
pos = pos + 24;

%Subframe 2
[subframe_2_data] = decode_subframe_2(msg(pos:pos+191)); pos = pos + 192;

IODE2   = subframe_2_data(1);
Crs     = subframe_2_data(2);
delta_n = subframe_2_data(3);
M0      = subframe_2_data(4);
Cuc     = subframe_2_data(5);
e       = subframe_2_data(6);
Cus     = subframe_2_data(7);
root_A  = subframe_2_data(8);
toe     = subframe_2_data(9);
fit_int = subframe_2_data(10);

%reserved (1 byte)
pos = pos + 8;

%skip WORD 2
pos = pos + 24;

%Subframe 3
[subframe_3_data] = decode_subframe_3(msg(pos:pos+191)); pos = pos + 192;

Cic      = subframe_3_data(1);
omega0   = subframe_3_data(2);
Cis      = subframe_3_data(3);
i0       = subframe_3_data(4);
Crc      = subframe_3_data(5);
omega    = subframe_3_data(6);
omegadot = subframe_3_data(7);
IODE3    = subframe_3_data(8);
IDOT     = subframe_3_data(9);

%output and reorder ephemerides data (if IODC == IODE)
if ((IODC == IODE2) && (IODC == IODE3) && constellations.GPS.enabled)
    data{2}(1) = PRN;
    data{2}(2) = af2;
    data{2}(3) = M0;
    data{2}(4) = root_A;
    data{2}(5) = delta_n;
    data{2}(6) = e;
    data{2}(7) = omega;
    data{2}(8) = Cuc;
    data{2}(9) = Cus;
    data{2}(10) = Crc;
    data{2}(11) = Crs;
    data{2}(12) = i0;
    data{2}(13) = IDOT;
    data{2}(14) = Cic;
    data{2}(15) = Cis;
    data{2}(16) = omega0;
    data{2}(17) = omegadot;
    data{2}(18) = toe;
    data{2}(19) = af0;
    data{2}(20) = af1;
    data{2}(21) = toc;
    data{2}(22) = IODE3;
    data{2}(23) = code_on_L2;
    data{2}(24) = weekno;
    data{2}(25) = L2flag;
    data{2}(26) = svaccur;
    data{2}(27) = svhealth;
    data{2}(28) = tgd;
    data{2}(29) = fit_int;
    data{2}(30) = constellations.GPS.indexes(PRN);
    data{2}(31) = int8('G');
    data{2}(32) = weektow2time(weekno, toe, 'G');
    data{2}(33) = weektow2time(weekno, toc, 'G');
end

% Check, no data --> delete header to improve performance
if (sum(data{2,1}(2:29)) == 0)
    data{1} = '';
end
