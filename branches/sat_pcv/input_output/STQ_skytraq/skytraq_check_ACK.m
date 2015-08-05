function [out] = skytraq_check_ACK(serialObj, MsgID)

% SYNTAX:
%   [out] = skytraq_check_ACK(serialObj, MsgID);
%
% INPUT:
%   serialObj = serial Object identifier
%   MsgID  = message ID of the requested message (hex)
%
% OUTPUT:
%   out = acknowledge outcome
%
% DESCRIPTION:
%   Check acknowledge reply after polling/sending SkyTraq messages.

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

out = 0;

header1 = 'A0';                            % header init (hexadecimal value)
header2 = 'A1';                            % header init (hexadecimal value)

header3 = '0D';                            % header close (hexadecimal value)
header4 = '0A';                            % header close (hexadecimal value)

ID  = 131;
payload_length = 2;
message_body = MsgID;

codeHEX = [header1; header2];
HDR = hex2dec(codeHEX);

ACK_DEC = [HDR; 0; payload_length; ID; message_body];

% checksum
CS = 0;
CS = bitxor(CS,ID);
CS = bitxor(CS,message_body);

codeHEX = [header3; header4];
HDR = hex2dec(codeHEX);

ACK_DEC = [ACK_DEC; CS; HDR];

%time acquisition
start_time = toc;

%maximum waiting time
dtMax = 0.5;

reply_1 = 0;
reply_2 = 0;

while (reply_1 ~= reply_2) | (reply_1 == 0)
    
    %time acquisition
    current_time = toc;
    
    %check if maximum waiting time is expired
    if (current_time - start_time > dtMax)
        return
    end
    
    % serial port checking
    reply_1 = get(serialObj, 'BytesAvailable');
    pause(0.1);
    reply_2 = get(serialObj, 'BytesAvailable');
end

reply = fread(serialObj, reply_1, 'uint8');

% search for acknowledge in reply
index = strfind(reply',ACK_DEC(1:6)');

if (~isempty(index) & length(reply(index:end)) >= 9)
    % extract acknowledge message from reply
    reply = reply(index:index+8);

    if (reply == ACK_DEC) %#ok<BDSCI>
        out = 1;
    end
end