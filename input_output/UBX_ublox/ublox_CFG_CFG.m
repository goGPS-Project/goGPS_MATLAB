function [out] = ublox_CFG_CFG(serialObj, action)

% SYNTAX:
%   [out] = ublox_CFG_CFG(serialObj, action)
%
% INPUT:
%   serialObj = serial Object identifier
%   action = 'clear': load factory defaults to active settings
%            'save' : save active settings to non-volatile memory
%            'load' : load settings from non-volatile memory to active settings
%
% OUTPUT:
%   out = outcome of the request
%
% DESCRIPTION:
%   Clear, save or load u-blox receiver configurations.

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

header1 = 'B5';                                % header (hexadecimal value)
header2 = '62';                                % header (hexadecimal value)
ID1 = '06';                                    % CFG (hexadecimal value)
ID2 = '09';                                    % CFG (hexadecimal value)

codeHEX = [header1; header2; ID1; ID2];        % initial hexadecimal poll message
codeDEC = hex2dec(codeHEX);                    % conversion to decimal

LEN1 = 0;                                      % payload length
LEN2 = 13;                                     % (little-endian)
codeDEC = [codeDEC; LEN2; LEN1];

clear_mask = ['00';'00';'00';'00'];
save_mask  = ['00';'00';'00';'00'];
load_mask  = ['00';'00';'00';'00'];

if (strcmp(action, 'clear'))

    clear_mask = ['FF';'FF';'00';'00'];
    load_mask  = ['FF';'FF';'00';'00'];

elseif (strcmp(action, 'save'))

    save_mask = ['FF';'FF';'00';'00'];

elseif (strcmp(action, 'load'))

    load_mask = ['FF';'FF';'00';'00'];

end

%devices = '17';                                % devices: BBR, FLASH, I2C-EEPROM, SPI-FLASH
devices = '07';                                % devices: BBR, FLASH, I2C-EEPROM

codeDEC = [codeDEC; hex2dec(clear_mask); hex2dec(save_mask); hex2dec(load_mask); hex2dec(devices)];

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
    pause(0.01);
end

% send message
try
    fwrite(serialObj, codeDEC, 'uint8', 'async');
    pause(0.01);
catch
    stopasync(serialObj);
    fwrite(serialObj, codeDEC, 'uint8', 'async');
    pause(0.01);
end

[out] = ublox_check_ACK(serialObj, ID1, ID2);
pause(0.01);
