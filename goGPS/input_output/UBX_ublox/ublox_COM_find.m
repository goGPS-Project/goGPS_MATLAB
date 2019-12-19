function [COMPort] = ublox_COM_find()

% SYNTAX:
%   [COMPort] = ublox_COM_find()
%
% OUTPUT:
%   COMPort = COM port on which the u-blox receiver has been found
%
% DESCRIPTION:
%   Scans all the COM ports and tries to detect if an u-blox receiver is connected.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

header1 = 'B5';      % header (hexadecimal value)
header2 = '62';      % header (hexadecimal value)

Class = '02';
MsgId = '10';

COMPort = [];

% start detection
fprintf('Detecting u-blox COM port...\n');

serialInfo = instrhwinfo('serial');
for i = 1 : length(serialInfo.AvailableSerialPorts)

    s = serial(serialInfo.AvailableSerialPorts(i),'BaudRate',57600);
    set(s,'InputBufferSize',16384);
    if ~ismac
        set(s,'FlowControl','hardware');
    else
        set(s,'FlowControl','software');
    end
    try
        fopen(s);
    catch
        break
    end

    ublox_poll_message(s, 'RXM', 'RAW', 0);

    bytes_1 = 0;
    bytes_2 = 0;

    while (bytes_1 ~= bytes_2) | (bytes_1 == 0)

        bytes_1 = get(s, 'BytesAvailable');
        pause(0.5);
        bytes_2 = get(s, 'BytesAvailable');
    end

    reply = fread(s, bytes_1, 'uint8');
    replyBIN = dec2bin(reply,8);
    replyBIN = replyBIN';
    replyBIN = replyBIN(:)';

    codeHEX = [header1 header2 Class MsgId];
    codeBIN = dec2bin(hex2dec(codeHEX),32);

    pos = strfind(replyBIN, codeBIN);

    if (~isempty(pos))
        LEN = fbin2dec(replyBIN(pos(1)+32:pos(1)+39)) + (fbin2dec(replyBIN(pos(1)+40:pos(1)+47)) * 2^8);

        if (LEN ~= 0)
            COMPort = serialInfo.AvailableSerialPorts(i);

            % successful detection
            fprintf('u-blox receiver detected on port %s.\n', COMPort{1});

            fclose(s);
            delete(s)
            clear s
            break
        end
    end

    fclose(s);
    delete(s)
    clear s
end

if isempty(COMPort)
    % unsuccessful detection
    fprintf('u-blox receiver could not be detected.\nYou can try to connect the device to the USB port BEFORE starting MATLAB.\n');
    answer = input('Do you want to specify the COM port manually? (Y/N): ','s');
    if (answer == 'Y') | (answer == 'y')
        COMPort = input('Please specify the COM port (e.g. COM6): ','s');
    end
end
