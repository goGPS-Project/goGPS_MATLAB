function [data] = decode_NAV_SBAS(msg)

% SYNTAX:
%   [data] = decode_NAV_SBAS(msg)
%
% INPUT:
%   msg = message transmitted by the u-blox receiver
%
% OUTPUT:
%   data = cell-array that contains the NAV_SBAS packet information
%          1.1) message class-id (NAV_SBAS)
%          2.1) TOW  = week time (in seconds)
%          2.2) GEO  = PRN Number of the GEO where correction and integrity data is used from
%          2.3) MODE = SBAS Mode: 0 Disabled, 1 Enabled Integrity, 3 Enabled Testmode 
%          2.4) SYS  = SBAS System (WAAS/EGNOS/MSAS) : -1 Unknown, 0 WAAS, 1 EGNOS, 2 MSAS, 16 GPS 
%          2.5) SERV = SBAS Services available: bit0 Ranging, bit1 Corrections, bit2 Integrity, bit3 Testmode 
%          2.6) CNT  = Number of SV data following
%          2.7) RES  = reserved field (not used) [3 bytes, 2.7)RES1, 2.8)RES2, 2.9)RES3]
%          3.1) SVID = SV id
%          3.2) FLAGS= Flags for this SV 
%          3.3) UDRE = Monitoring status
%          3.4) SYSn = System (WAAS/EGNOS/MSAS) same as SYS
%          3.5) SERVn= Services available same as SERVICE
%          3.6) RES0n= reserved field (not used) 
%          3.7) PRC  = Pseudo Range correction (in m) 
%          3.8) RES1n= reserved field (not used) 
%          3.9) IC   = Ionosphere correction (in m)
%
% DESCRIPTION:
%   NAV_SBAS binary message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini, Antonio Herrera
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
data{3} = zeros(32,9);

%output data save
data{1} = 'NAV-SBAS';

% week time decoding (4 byte) [in ms]
TOW1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW3 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW4 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW = TOW1 + (TOW2 * 2^8) + (TOW3 * 2^16) + (TOW4 * 2^24);  % little endian
TOW = TOW / 1000;

% PRN Number of the GEO where correction and integrity data is used from (1 byte)
GEO = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

% MODE = SBAS Mode: 0 Disabled, 1 Enabled Integrity, 3 Enabled Testmode (1 byte)
MODE = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

% SBAS System (WAAS/EGNOS/MSAS) : -1 Unknown, 0 WAAS, 1 EGNOS, 2 MSAS, 16 GPS (1 byte)
SYS = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

% SBAS Services available: bit0 Ranging, bit1 Corrections, bit2 Integrity, bit3 Testmode (1 byte)    ******** Info bit
SERV = fbin2dec(msg(pos:pos+7));  pos = pos + 8;  % i nº dec indica que hay.

% Number of SV data following (1 byte)
CNT = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

% reserved field (3 byte)
RES1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
RES2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
RES3 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

%output data save
data{2}(1) = TOW;
data{2}(2) = GEO;
data{2}(3) = MODE;
data{2}(4) = SYS;
data{2}(5) = SERV;
data{2}(6) = CNT;
data{2}(7) = RES1;
data{2}(8) = RES2;
data{2}(9) = RES3;

%read the measurements of every satellite
for j = 1 : CNT

    % satellite number decoding (1 byte)
    SVID = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    % Flags for this SV (1 byte)
    FLAGS = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    % Monitoring status (1 byte)
    UDRE = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    % System (WAAS/EGNOS/MSAS): : 0 WAAS, 1 EGNOS, 2 MSAS, 16 GPS (1 byte)
    SYSn = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    % SBAS Services available: bit0 Ranging, bit1 Corrections, bit2 Integrity, bit3 Testmode (1 byte)    ******** Info bit
    SERVn = fbin2dec(msg(pos:pos+7));  pos = pos + 8; % i nº dec indica que hay.
    
    % Reserved field (not used) [U1] (1 byte)
    RES0n = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    % Pseudo Range correction [in cm] --> [m] (2 byte)
    PRC1 = msg(pos:pos+7);  pos = pos + 8;
    PRC2 = msg(pos:pos+7);  pos = pos + 8;
    PRC  = twos_complement([PRC2 PRC1]);          % little endian
    PRC  = PRC / 100;
    
    % Reserved field (not used) [I2] (2 byte)
    RES1n1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    RES1n2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    RES1n  = RES1n1 + (RES1n2 * 2^8);        % little endian
    
    % Ionosphere correction [in cm] --> [m] (2 byte)
    IC1 = msg(pos:pos+7);  pos = pos + 8;
    IC2 = msg(pos:pos+7);  pos = pos + 8;
    IC  = twos_complement([IC2 IC1]);          % little endian
    IC  = IC / 100;

    % data output save
    data{3}(SVID,1) = SVID;
    data{3}(SVID,2) = FLAGS;
    data{3}(SVID,3) = UDRE;
    data{3}(SVID,4) = SYSn;
    data{3}(SVID,5) = SERVn;
    data{3}(SVID,6) = RES0n;
    data{3}(SVID,7) = PRC;
    data{3}(SVID,8) = RES1n;
    data{3}(SVID,9) = IC;

end