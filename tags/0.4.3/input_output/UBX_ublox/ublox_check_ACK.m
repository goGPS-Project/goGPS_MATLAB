function [out] = ublox_check_ACK(serialObj, Class, MsgID)

% SYNTAX:
%   [out] = ublox_check_ACK(serialObj, Class, MsgID);
%
% INPUT:
%   serialObj = serial Object identifier
%   Class  = u-blox message class (hexadecimal value)
%   MsgID  = u-blox message ID (hexadecimal value)
%
% OUTPUT:
%   out = acknowledge outcome
%
% DESCRIPTION:
%   Check acknowledge reply after polling u-blox messages.

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

% ACK-ACK message (without checksum)
ACK_HEX = ['B5'; '62'; '05'; '01'; '02'; '00'; Class; MsgID];
ACK_DEC = hex2dec(ACK_HEX);

% checksum
CK_A = 0; CK_B = 0;
for i = 3 : length(ACK_DEC)
    CK_A = CK_A + ACK_DEC(i);
    CK_B = CK_B + CK_A;
end

CK_A = mod(CK_A,256);
CK_B = mod(CK_B,256);

ACK = [ACK_DEC; CK_A; CK_B];

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
index = strfind(reply',ACK_DEC(1:4)');

if (~isempty(index) & length(reply(index:end)) >= 10)
    % extract acknowledge message from reply
    reply = reply(index:index+9);

    if (reply == ACK) %#ok<BDSCI>
        out = 1;
    end
end