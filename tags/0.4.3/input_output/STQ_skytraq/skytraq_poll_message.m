function skytraq_poll_message(serialObj, MsgID, parameter)

% SYNTAX:
%   skytraq_poll_message(serialObj, MsgID, parameter);
%
% INPUT:
%   serialObj = serial Object identifier
%   MsgIDLab  = SkyTraq message ID ('11' = almanac; '30' = ephemeris)
%   parameter = parameter sent to poll all satellites or one specific
%               satellite data (0 = all; 1-32 specific satellite)
%
% DESCRIPTION:
%   Poll SkyTraq alamanac or ephemeris.

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

payload_length = 2;
ID  = hex2dec(MsgID);
message_body = parameter;

codeHEX = [header1; header2];
HDR = hex2dec(codeHEX);

codeDEC = [HDR; 0; payload_length; ID; message_body];

% checksum
CS = 0;
CS = bitxor(CS,ID);
CS = bitxor(CS,message_body);

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