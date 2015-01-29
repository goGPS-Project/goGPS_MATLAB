function [out] = skytraq_message_format(serialObj, format, memory)

% SYNTAX:
%   [out] = skytraq_message_format(serialObj, format, memory);
%
% INPUT:
%   serialObj = serial Object identifier
%   format = parameter to specify whether to set NMEA or binary output
%            (accepted values: 'NOOUT', 'NMEA', 'BIN'); default = 'BIN'
%   memory = parameter to specify whether to update only the SRAM or
%            SRAM+FLASH memory (accpeted values: 'SRAM', 'FLASH'); default = 'SRAM'
%
% OUTPUT:
%   out = receiver reply
%
% DESCRIPTION:
%   Set NMEA or binary format output on SkyTraq receivers.

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

header1 = 'A0';                            % header init (hexadecimal value)
header2 = 'A1';                            % header init (hexadecimal value)

header3 = '0D';                            % header close (hexadecimal value)
header4 = '0A';                            % header close (hexadecimal value)

ID  = 9;
payload_length = 3;

%default values
message_body1 = 2; %output binary data
message_body2 = 0; %update only SRAM

if (nargin > 1)
    if (strcmp(format, 'NOOUT'))
        message_body1 = 0;
    elseif (strcmp(format, 'NMEA'))
        message_body1 = 1;
    end
    if (strcmp(memory, 'FLASH'))
        message_body2 = 1;
    end
end
message_body = [message_body1; message_body2];

codeHEX = [header1; header2];
HDR = hex2dec(codeHEX);

codeDEC = [HDR; 0; payload_length; ID; message_body];

% checksum
CS = 0;
CS = bitxor(CS,ID);
CS = bitxor(CS,message_body1);
CS = bitxor(CS,message_body2);

codeHEX = [header3; header4];
HDR = hex2dec(codeHEX);

codeDEC = [codeDEC; CS; HDR];

% %serial port checking
% reply_1 = get(serialObj, 'BytesAvailable');
%
% if (reply_1 ~= 0)
%     %clear the serial port (data not decoded)
%     reply = fread(serialObj, reply_1, 'uint8');
% end

% send message
% try
    fwrite(serialObj, codeDEC, 'uint8', 'async');
% catch
%     stopasync(serialObj);
%     fwrite(serialObj, codeDEC, 'uint8', 'async');
% end

[out] = skytraq_check_ACK(serialObj, ID);