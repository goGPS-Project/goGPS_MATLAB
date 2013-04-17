function [out] = ublox_CFG_SBAS(serialObj, action)

% SYNTAX:
%   [out] = ublox_CFG_SBAS(serialObj, action)
%
% INPUT:
%   serialObj = serial Object identifier
%   action = 'default': load defaults settings
%            
%
% OUTPUT:
%   out = outcome of the request
%
% DESCRIPTION:
%   SBAS u-blox receiver configurations.

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

header1 = 'B5';                                % header (hexadecimal value)
header2 = '62';                                % header (hexadecimal value)
ID1 = '06';                                    % CFG (hexadecimal value)
ID2 = '16';                                    % SBAS (hexadecimal value)

codeHEX = [header1; header2; ID1; ID2];        % initial hexadecimal poll message
codeDEC = hex2dec(codeHEX);                    % conversion to decimal

LEN1 = 0;                                      % payload length
LEN2 = 8;                                     % (little-endian)
codeDEC = [codeDEC; LEN2; LEN1];


mode_mask = ['00'];                              % 1 Byte
use_mask  = ['00'];                              % 1 Byte
nch_SBAS  = ['00'];                              % 1 Byte, Max nº channel SBAS  range:0-3
res_byte  = ['00'];                              % 1 Byte, Reserved for future use.
scan_mask = ['00';'00';'00';'00'];               % 4 Byte, If all Bits are set to zero, auto-scan (i.e. all valid PRNs) are searched


% Define Action 

% SBAS Mode. (Bitmask) yes
% Bit 0: SBAS Enabled (1) / Disabled (0)
% Bit 1: SBAS Testbed: Use data anyhow (1) / Ignore data when in Test Mode (SBAS Msg 0) siempre 0
% Bits 2-7: reserved for future use 


% SBAS Usage (Bitmask)
% Bit 0: Use SBAS GEOs as a ranging source (for navigation) 1
% Bit 1: Use SBAS Differential Corrections 1
% Bit 2: Use SBAS Integrity Information   0 at the moment

% Max nº channel SBAS  (valid range:0-3) 3

% Which SBAS PRN numbers to search for (Bitmask)
% If all Bits are set to zero, auto-scan (i.e. all valid PRNs) are searched.
% Every bit corresponds to a PRN number
% Bit 0: PRN120
% Bit 1: PRN121
% ....
% Bit 18: PRN138
% Bits 19-31: reserved. set to zero 


% Action Default
if (strcmp(action, 'default'))
    
mode_mask = ['01'];                 % SBAS Enabled and Ignore data when in Test Mode
% use_mask  = ['03'];                 % SBAS GEO and SBAS Differential Corrections
% use_mask  = ['00'];                 % Don´t use
use_mask  = ['02'];                 % Corrections only
nch_SBAS  = ['03'];                 % Max nº channel SBAS: 3 (valid range:0-3)
res_byte  = ['00'];                 % Reserved for future use.
scan_mask = ['00';'00';'00';'00'];  % All valid PRNs

end


codeDEC = [codeDEC; hex2dec(mode_mask); hex2dec(use_mask); hex2dec(nch_SBAS); hex2dec(res_byte); hex2dec(scan_mask)];

% checksum
CK_A = 0; CK_B = 0;
for i = 3 : length(codeDEC)
    CK_A = CK_A + codeDEC(i);
    CK_B = CK_B + CK_A;
end

CK_A = mod(CK_A,256);
CK_B = mod(CK_B,256);

codeDEC = [codeDEC; CK_A; CK_B];

%serial port checking
reply_1 = get(serialObj,'BytesAvailable');

if (reply_1 ~= 0)
    %clear the serial port (data not decoded)
    reply = fread(serialObj, reply_1, 'uint8'); %#ok<NASGU>
end

% send message
try
    fwrite(serialObj, codeDEC, 'uint8', 'async');
catch
    stopasync(serialObj);
    fwrite(serialObj, codeDEC, 'uint8', 'async');
end

[out] = ublox_check_ACK(serialObj, ID1, ID2);