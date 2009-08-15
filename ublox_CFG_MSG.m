function [reply] = ublox_CFG_MSG(serialObj, ClassLab, MsgIDLab, mode)

% SYNTAX:
%   [reply] = ublox_CFG_MSG(serialObj, Class, MsgID, mode)
%
% INPUT:
%   serialObj = serial Object identifier
%   ClassLab  = u-blox message class (label - e.g. 'RXM')
%   MsgIDLab  = u-blox message ID (label - e.g. 'RAW')
%   mode = 0: disable periodic message
%          1: enable periodic message
%
% OUTPUT:
%   reply = receiver reply
%
% DESCRIPTION:
%   Poll, enable or disable u-blox messages.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
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

[Class, MsgID] = ublox_UBX_codes(ClassLab, MsgIDLab);

header1 = 'B5';                            % header (hexadecimal value)
header2 = '62';                            % header (hexadecimal value)

ID1 = '06';                                % CFG (hexadecimal value)
ID2 = '01';                                % MSG (hexadecimal value)

codeHEX = [header1; header2; ID1; ID2];    % initial hexadecimal poll message
codeDEC = hex2dec(codeHEX);                % conversion to decimal

LEN1 = 0;                                  % payload length
LEN2 = 3;                                  % (enable/disable a periodic message)
codeDEC = [codeDEC; LEN2; LEN1];           % (little-endian)
if (mode == 1)
    Rate = '01';                           % enable periodic message
else
    Rate = '00';                           % disable periodic message
end
codeHEX = [Class; MsgID; Rate];
codeDEC = [codeDEC; hex2dec(codeHEX)];

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
    reply = fread(serialObj, reply_1, 'uint8');
end

% send message
fwrite(serialObj, codeDEC, 'uint8', 'async');

reply_1 = 0;
reply_2 = 0;

while (reply_1 ~= reply_2) | (reply_1 == 0)
    % serial port checking
    reply_1 = get(serialObj, 'BytesAvailable');
    pause(0.05);
    reply_2 = get(serialObj, 'BytesAvailable');
end

reply = fread(serialObj, 10, 'uint8');                              % read Acknowledge

ACK_HEX = ['B5'; '62'; '05'; '01'; '02'; '00'; ID1; ID2];           % ACK-ACK message (without checksum)
ACK_DEC = hex2dec(ACK_HEX);                                         % conversion to decimal

% checksum
CK_A = 0; CK_B = 0;
for i = 3 : length(ACK_DEC)
    CK_A = CK_A + ACK_DEC(i);
    CK_B = CK_B + CK_A;
end

CK_A = mod(CK_A,256);
CK_B = mod(CK_B,256);

ACK = [ACK_DEC; CK_A; CK_B];

if (mode == 1) % periodic message enable request

    if (reply == ACK)
        %fprintf('u-blox %s-%s message enabled.\n', ClassLab, MsgIDLab);
        % positive reply
        reply = 1;
    else
        %fprintf('It was not possible to enable u-blox %s-%s messages.\n', ClassLab, MsgIDLab);
        % negative reply
        reply = 0;
    end

else

    if (reply == ACK)
        %fprintf('u-blox %s-%s message disabled.\n', ClassLab, MsgIDLab);
        % positive reply
        reply = 1;
    else
        %fprintf('It was not possible to disable u-blox %s-%s messages.\n', ClassLab, MsgIDLab);
        % negative reply
        reply = 0;
    end
end
