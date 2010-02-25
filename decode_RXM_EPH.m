function [data] = decode_RXM_EPH(msg)

% SYNTAX:
%   [data] = decode_RXM_EPH(msg)
%
% INPUT:
%   msg = message transmitted by the u-blox receiver
%
% OUTPUT:
%   data = cell-array that contains the RXM-EPH packet information
%          1.1) message class-id (RXM-EPH)
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
%
% DESCRIPTION:
%   RXM-EPH binary message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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

% first message initial index
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(21,1);

%output data save
data{1} = 'RXM-EPH';

%satellite SVN (4 bytes)
SVN = msg(pos:pos+31); pos = pos + 32;
SVN = fliplr(reshape(SVN,8,[]));                  % byte order inversion (little endian)
SVN = SVN(:)';
SVN = bin2dec((SVN(1:32)));

%Hand-Over Word
HOW = msg(pos:pos+31); pos = pos + 32;
HOW = fliplr(reshape(HOW,8,[]));                  % byte order inversion (little endian)
HOW = HOW(:)';

HOW = bin2dec(HOW(28:30)); %#ok<*NASGU>

%Subframe 1
SF1D0 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 0 (=Word 3 in ICD-GPS-200)
SF1D1 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 1 (=Word 4 in ICD-GPS-200)
SF1D2 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 2 (=Word 5 in ICD-GPS-200)
SF1D3 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 3 (=Word 6 in ICD-GPS-200)
SF1D4 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 4 (=Word 7 in ICD-GPS-200)
SF1D5 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 5 (=Word 8 in ICD-GPS-200)
SF1D6 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 6 (=Word 9 in ICD-GPS-200)
SF1D7 = msg(pos:pos+31); pos = pos + 32;          % subframe 1 - Word 7 (=Word 10 in ICD-GPS-200)

SF1D5 = fliplr(reshape(SF1D5,8,[])); SF1D5 = SF1D5(:)';               % byte order inversion (little endian)
SF1D6 = fliplr(reshape(SF1D6,8,[])); SF1D6 = SF1D6(:)';               %
SF1D7 = fliplr(reshape(SF1D7,8,[])); SF1D7 = SF1D7(:)';               %

IODC_8LSBs = bin2dec(SF1D5(9:16));
toc = bin2dec(SF1D5(17:32)) * (2^4);
af2 = twos_complement(SF1D6(9:16)) * (2^-55);
af1 = twos_complement(SF1D6(17:32)) * (2^-43);
af0 = twos_complement(SF1D7(9:30)) * (2^-31);

%Subframe 2
SF2D0 = msg(pos:pos+31); pos = pos + 32;          % subframe 2 - Word 0 (=Word 3 in ICD-GPS-200)
SF2D1 = msg(pos:pos+31); pos = pos + 32;          % subframe 2 - Word 1 (=Word 4 in ICD-GPS-200)
SF2D2 = msg(pos:pos+31); pos = pos + 32;          % subframe 2 - Word 2 (=Word 5 in ICD-GPS-200)
SF2D3 = msg(pos:pos+31); pos = pos + 32;          % subframe 2 - Word 3 (=Word 6 in ICD-GPS-200)
SF2D4 = msg(pos:pos+31); pos = pos + 32;          % subframe 2 - Word 4 (=Word 7 in ICD-GPS-200)
SF2D5 = msg(pos:pos+31); pos = pos + 32;          % subframe 2 - Word 5 (=Word 8 in ICD-GPS-200)
SF2D6 = msg(pos:pos+31); pos = pos + 32;          % subframe 2 - Word 6 (=Word 9 in ICD-GPS-200)
SF2D7 = msg(pos:pos+31); pos = pos + 32;          % subframe 2 - Word 7 (=Word 10 in ICD-GPS-200)

SF2D0 = fliplr(reshape(SF2D0,8,[])); SF2D0 = SF2D0(:)';               %
SF2D1 = fliplr(reshape(SF2D1,8,[])); SF2D1 = SF2D1(:)';               %
SF2D2 = fliplr(reshape(SF2D2,8,[])); SF2D2 = SF2D2(:)';               %
SF2D3 = fliplr(reshape(SF2D3,8,[])); SF2D3 = SF2D3(:)';               % byte order inversion (little endian)
SF2D4 = fliplr(reshape(SF2D4,8,[])); SF2D4 = SF2D4(:)';               %
SF2D5 = fliplr(reshape(SF2D5,8,[])); SF2D5 = SF2D5(:)';               %
SF2D6 = fliplr(reshape(SF2D6,8,[])); SF2D6 = SF2D6(:)';               %
SF2D7 = fliplr(reshape(SF2D7,8,[])); SF2D7 = SF2D7(:)';               %

IODE2 = bin2dec(SF2D0(9:16));
Crs = twos_complement(SF2D0(17:32)) * (2^-5);
delta_n = twos_complement(SF2D1(9:24)) * pi * (2^-43);
M0 = twos_complement([SF2D1(25:32) SF2D2(9:32)]) * pi * (2^-31);
Cuc = twos_complement(SF2D3(9:24)) * (2^-29);
e = bin2dec([SF2D3(25:32) SF2D4(9:32)]) * (2^-33);
Cus = twos_complement(SF2D5(9:24)) * (2^-29);
root_A = bin2dec([SF2D5(25:32) SF2D6(9:32)]) * (2^-19);
toe = bin2dec(SF2D7(9:24)) * (2^4);

%Subframe 3
SF3D0 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 0 (=Word 3 in ICD-GPS-200)
SF3D1 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 1 (=Word 4 in ICD-GPS-200)
SF3D2 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 2 (=Word 5 in ICD-GPS-200)
SF3D3 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 3 (=Word 6 in ICD-GPS-200)
SF3D4 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 4 (=Word 7 in ICD-GPS-200)
SF3D5 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 5 (=Word 8 in ICD-GPS-200)
SF3D6 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 6 (=Word 9 in ICD-GPS-200)
SF3D7 = msg(pos:pos+31); pos = pos + 32;          % subframe 3 - Word 7 (=Word 10 in ICD-GPS-200)

SF3D0 = fliplr(reshape(SF3D0,8,[])); SF3D0 = SF3D0(:)';               %
SF3D1 = fliplr(reshape(SF3D1,8,[])); SF3D1 = SF3D1(:)';               %
SF3D2 = fliplr(reshape(SF3D2,8,[])); SF3D2 = SF3D2(:)';               %
SF3D3 = fliplr(reshape(SF3D3,8,[])); SF3D3 = SF3D3(:)';               % byte order inversion (little endian)
SF3D4 = fliplr(reshape(SF3D4,8,[])); SF3D4 = SF3D4(:)';               %
SF3D5 = fliplr(reshape(SF3D5,8,[])); SF3D5 = SF3D5(:)';               %
SF3D6 = fliplr(reshape(SF3D6,8,[])); SF3D6 = SF3D6(:)';               %
SF3D7 = fliplr(reshape(SF3D7,8,[])); SF3D7 = SF3D7(:)';               %

Cic = twos_complement(SF3D0(9:24)) * (2^-29);
omega0 = twos_complement([SF3D0(25:32) SF3D1(9:32)]) * pi * (2^-31);
Cis = twos_complement(SF3D2(9:24)) * (2^-29);
i0 = twos_complement([SF3D2(25:32) SF3D3(9:32)]) * pi * (2^-31);
Crc = twos_complement(SF3D4(9:24)) * (2^-5);
omega = twos_complement([SF3D4(25:32) SF3D5(9:32)]) * pi * (2^-31);
omegadot = twos_complement(SF3D6(9:32)) * pi * (2^-43);
IODE3 = bin2dec(SF3D7(9:16));
IDOT = twos_complement(SF3D7(17:30)) * pi * (2^-43);

if (IODC_8LSBs == IODE2) & (IODC_8LSBs == IODE3)
    data{2}(1) = SVN;
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
end