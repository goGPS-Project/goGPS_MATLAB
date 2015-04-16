function [data] = decode_RXM_SFRB(msg, constellations)

% SYNTAX:
%   [data] = decode_RXM_SFRB(msg, constellations);
%
% INPUT:
%   msg = message transmitted by the u-blox receiver
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the RXM-SFRB packet information
%          1.1) message class-id (RXM-SFRB)
%          2.1) ionosphere parameter (a0)
%          2.2) ionosphere parameter (a1)
%          2.3) ionosphere parameter (a2)
%          2.4) ionosphere parameter (a3)
%          2.5) ionosphere parameter (b0)
%          2.6) ionosphere parameter (b1)
%          2.7) ionosphere parameter (b2)
%          2.8) ionosphere parameter (b3)
%          2.9) leap seconds
%          3.1) GPS satellite id
%          3.2) GPS af2
%          3.3) GPS M0
%          3.4) GPS root A
%          3.5) GPS delta-N
%          3.6) GPS eccentricity
%          3.7) GPS omega
%          3.8) GPS Cuc
%          3.9) GPS Cus
%          3.10)GPS Crc
%          3.11)GPS Crs
%          3.12)GPS i0
%          3.13)GPS IDOT
%          3.14)GPS Cic
%          3.15)GPS Cis
%          3.16)GPS omega0
%          3.17)GPS omegadot
%          3.18)GPS toe
%          3.19)GPS af0
%          3.20)GPS af1
%          3.21)GPS toc
%          3.22)GPS IODE
%          3.23)GPS codes;
%          3.24)GPS weekno;
%          3.25)GPS L2flag;
%          3.26)GPS svaccur;
%          3.27)GPS svhealth;
%          3.28)GPS tgd;
%          3.29)GPS fit_int;
%          3.30)multi-constellation satellite index (here only GPS is assumed)
%
% DESCRIPTION:
%   RXM-SFRB binary message decoding (OBSOLETE).

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

% first message initial index
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(9,1);
data{3} = zeros(33,1);

%output data save
data{1} = 'RXM-SFRB';

%channel number (1 byte)
CHN = fbin2dec(msg(pos:pos+7));  pos = pos + 8; %#ok<NASGU>

%satellite ID (1 byte)
PRN = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

%if GPS satellite
if (PRN <= 32)

    %check HOW to see which subframe is it
    pos = pos + 32;
    HOW = msg(pos:pos+31); pos = pos + 32;
    HOW = fliplr(reshape(HOW,8,[]));                  % byte order inversion (little endian)
    HOW = HOW(:)';
    %subframe ID
    SFID = fbin2dec(HOW(28:30));
    
    switch SFID
        case 1
            %Subframe 1
            [subframe_1_data]  = decode_subframe_1(msg(pos:pos+255));
            
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
            
        case 2
            %Subframe 2
            [subframe_2_data] = decode_subframe_2(msg(pos:pos+255));
            
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
            
        case 3
            %Subframe 3
            [subframe_3_data] = decode_subframe_3(msg(pos:pos+255));
            
            Cic      = subframe_3_data(1);
            omega0   = subframe_3_data(2);
            Cis      = subframe_3_data(3);
            i0       = subframe_3_data(4);
            Crc      = subframe_3_data(5);
            omega    = subframe_3_data(6);
            omegadot = subframe_3_data(7);
            IODE3    = subframe_3_data(8);
            IDOT     = subframe_3_data(9);
        case 4
            %check if it is page 18 of subframe 4
            WORD3 = msg(pos:pos+31);
            WORD3 = fliplr(reshape(WORD3,8,[])); WORD3 = WORD3(:)'; % byte order inversion (little endian)

            if (fbin2dec(WORD3(11:16)) == 56) %SVID "56" <--> page "18"
                %Subframe 4
                [subframe_4_data] = decode_subframe_4(msg(pos:pos+255));

                data{2}(1) = subframe_4_data(1); % a0
                data{2}(2) = subframe_4_data(2); % a1
                data{2}(3) = subframe_4_data(3); % a2
                data{2}(4) = subframe_4_data(4); % a3
                data{2}(5) = subframe_4_data(5); % b0
                data{2}(6) = subframe_4_data(6); % b1
                data{2}(7) = subframe_4_data(7); % b2
                data{2}(8) = subframe_4_data(8); % b3
                data{2}(9) = subframe_4_data(9); % leap seconds
            end
    end

    %output and reorder ephemerides data (if IODC == IODE)
    if ((IODC == IODE2) && (IODC == IODE3) && constellations.GPS.enabled)
        data{3}(1) = PRN;
        data{3}(2) = af2;
        data{3}(3) = M0;
        data{3}(4) = root_A;
        data{3}(5) = delta_n;
        data{3}(6) = e;
        data{3}(7) = omega;
        data{3}(8) = Cuc;
        data{3}(9) = Cus;
        data{3}(10) = Crc;
        data{3}(11) = Crs;
        data{3}(12) = i0;
        data{3}(13) = IDOT;
        data{3}(14) = Cic;
        data{3}(15) = Cis;
        data{3}(16) = omega0;
        data{3}(17) = omegadot;
        data{3}(18) = toe;
        data{3}(19) = af0;
        data{3}(20) = af1;
        data{3}(21) = toc;
        data{3}(22) = IODE3;
        data{3}(23) = code_on_L2;
        data{3}(24) = weekno;
        data{3}(25) = L2flag;
        data{3}(26) = svaccur;
        data{3}(27) = svhealth;
        data{3}(28) = tgd;
        data{3}(29) = fit_int;
        data{3}(30) = constellations.GPS.indexes(PRN);
        data{2}(31) = int8('G');
        data{2}(32) = weektow2time(weekno, toe, 'G');
        data{2}(33) = weektow2time(weekno, toc, 'G');
    end
end
