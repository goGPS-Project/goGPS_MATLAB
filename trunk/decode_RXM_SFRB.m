function [data] = decode_RXM_SFRB(msg)

% SYNTAX:
%   [data] = decode_RXM_SFRB(msg);
%
% INPUT:
%   msg = message transmitted by the u-blox receiver
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
%          2.9) leap seconds [NOT USED]
%
% DESCRIPTION:
%   RXM-SFRB binary message decoding (only subframe 4).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
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
data{2} = zeros(9,1);

%output data save
data{1} = 'RXM-SFRB';

%channel number (1 byte)
CHN = fbin2dec(msg(pos:pos+7));  pos = pos + 8; %#ok<NASGU>

%satellite ID (1 byte)
SVN = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

%if GPS satellite
if (SVN <= 32)

    %check HOW to see which subframe is it
    pos = pos + 32;
    HOW = msg(pos:pos+31); pos = pos + 32;
    HOW = fliplr(reshape(HOW,8,[]));                  % byte order inversion (little endian)
    HOW = HOW(:)';
    %subframe ID
    SFID = fbin2dec(HOW(28:30));
    
    switch SFID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% WARNING: subframes 1, 2, 3 read from RXM/AID-EPH messages %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         case 1
%             %Subframe 1
%             [subframe_1_data]  = decode_subframe_1(msg(pos:pos+255));
%             
%             weekno     = subframe_1_data(1);
%             code_on_L2 = subframe_1_data(2);
%             svaccur    = subframe_1_data(3);
%             svhealth   = subframe_1_data(4);
%             IODC       = subframe_1_data(5);
%             L2flag     = subframe_1_data(6);
%             tgd        = subframe_1_data(7);
%             toc        = subframe_1_data(8);
%             af2        = subframe_1_data(9);
%             af1        = subframe_1_data(10);
%             af0        = subframe_1_data(11);
            
%         case 2
%             %Subframe 2
%             [subframe_2_data] = decode_subframe_2(msg(pos:pos+255));
%             
%             IODE2   = subframe_2_data(1);
%             Crs     = subframe_2_data(2);
%             delta_n = subframe_2_data(3);
%             M0      = subframe_2_data(4);
%             Cuc     = subframe_2_data(5);
%             e       = subframe_2_data(6);
%             Cus     = subframe_2_data(7);
%             root_A  = subframe_2_data(8);
%             toe     = subframe_2_data(9);
%             fit_int = subframe_2_data(10);
            
%         case 3
%             %Subframe 3
%             [subframe_3_data] = decode_subframe_3(msg(pos:pos+255));
%             
%             Cic      = subframe_3_data(1);
%             omega0   = subframe_3_data(2);
%             Cis      = subframe_3_data(3);
%             i0       = subframe_3_data(4);
%             Crc      = subframe_3_data(5);
%             omega    = subframe_3_data(6);
%             omegadot = subframe_3_data(7);
%             IODE3    = subframe_3_data(8);
%             IDOT     = subframe_3_data(9);
            
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

%     %output and reorder ephemerides data (if IODC == IODE)
%     if (IODC == IODE2) & (IODC == IODE3)
%         data{2}(1) = SVN;
%         data{2}(2) = af2;
%         data{2}(3) = M0;
%         data{2}(4) = root_A;
%         data{2}(5) = delta_n;
%         data{2}(6) = e;
%         data{2}(7) = omega;
%         data{2}(8) = Cuc;
%         data{2}(9) = Cus;
%         data{2}(10) = Crc;
%         data{2}(11) = Crs;
%         data{2}(12) = i0;
%         data{2}(13) = IDOT;
%         data{2}(14) = Cic;
%         data{2}(15) = Cis;
%         data{2}(16) = omega0;
%         data{2}(17) = omegadot;
%         data{2}(18) = toe;
%         data{2}(19) = af0;
%         data{2}(20) = af1;
%         data{2}(21) = toc;
%         data{2}(22) = IODE3;
%         data{2}(23) = code_on_L2;
%         data{2}(24) = weekno;
%         data{2}(25) = L2flag;
%         data{2}(26) = svaccur;
%         data{2}(27) = svhealth;
%         data{2}(28) = tgd;
%         data{2}(29) = fit_int;
%     end
end