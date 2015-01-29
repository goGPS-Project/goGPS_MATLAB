function [out] = skytraq_binary_output_rate(serialObj, rate)

% SYNTAX:
%   [out] = skytraq_binary_output_rate(serialObj, rate);
%
% INPUT:
%   serialObj = serial Object identifier
%   rate = binary measurement output rate for Meas_time / Raw_meas / SV_CH_Status
%          (supported output rate configurations 1 / 2 / 4 / 5 / 10 / 20 Hz)
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

ID  = 18;
payload_length = 8;

%default value
message_body1 = 0; %measurement rate

%fixed values (settings hard-coded for goGPS; can be changed in the future...)
message_body2 = 1; %MEAS_TIME enabled
message_body3 = 1; %RAW_MEAS enabled
message_body4 = 0; %SV_CH_Status disabled
message_body5 = 0; %RCV_State disabled
message_body6 = 0; %Subframe disabled
message_body7 = 0; %update only SRAM

if (nargin > 1)
    switch rate
        case 1
            message_body1 = 0;
        case 2
            message_body1 = 1;
        case 4
            message_body1 = 2;
        case 5
            message_body1 = 3;
        case 10
            message_body1 = 4;
        case 20
            message_body1 = 5;
    end
end
message_body = [message_body1; message_body2; message_body3; message_body4; message_body5; message_body6; message_body7];

codeHEX = [header1; header2];
HDR = hex2dec(codeHEX);

codeDEC = [HDR; 0; payload_length; ID; message_body];

% checksum
CS = 0;
CS = bitxor(CS,ID);
for i = 1 : 7
    CS = bitxor(CS,message_body(i));
end

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