function [out] = nvs_message_format(serialObj, format, baudrate)

% SYNTAX:
%   [out] = nvs_message_format(serialObj, format, baudrate);
%
% INPUT:
%   serialObj = serial Object identifier
%   format = parameter to specify whether to set NMEA or binary output
%            (accepted values: 'NMEA', 'BIN'); default = 'BIN'
%   baudrate = requested baudrate
%
% OUTPUT:
%   out = receiver reply
%
% DESCRIPTION:
%   Set NMEA or BINR binary format output on NVS receivers.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.1 beta
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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

%default values
if (nargin < 2)
    format = 'BIN';
end
if (nargin < 3)
    baudrate = 115200;
end

%BINR data
DLE = '10'; %beginning of message
ETX = '03'; %end of message
BINR_0Bh = '0B';

%standard values for checking replies
check_NMEA = uint8('PORZA')';
check_BINR = hex2dec(['10';'50']);

%check the status of the current port (i.e. whether it's sending NMEA or BINR)
nmeastring = '$CCGPQ,PORZA';
checksum = NMEA_checksum(nmeastring);
nmeastring = [nmeastring '*' checksum];
codeDEC_NMEA = [uint8(nmeastring)'; 13; 10];
codeDEC_BINR = hex2dec([DLE; BINR_0Bh; DLE; ETX]);

flag_nmea = [];

for i = 1 : 2
    if (i == 1)
        codeDEC = codeDEC_NMEA;
        check   = check_NMEA;
        set(serialObj,'Parity','none');
    else
        codeDEC = codeDEC_BINR;
        check   = check_BINR;
        set(serialObj,'Parity','odd');
    end

    flushinput(serialObj);
    
    fwrite(serialObj, codeDEC, 'uint8', 'async');

    reply = get_reply(serialObj);

    %check reply
    if (~isempty(strfind(reply',check')))
        if (i == 1)
            flag_nmea = 1;
        else %i = 2
            flag_nmea = 0;
        end
        break
    end
end

if (isempty(flag_nmea))
    set(serialObj,'Parity','none');
    flag_nmea = 1;
    check = check_NMEA;
end

%set user-requested format
switch (format)
    case 'NMEA'
        if (flag_nmea)
            protocol = 1;
        else
            protocol = '02';
        end
    otherwise %BIN
        if (flag_nmea)
            protocol = 3;
        else
            protocol = '04';
        end
end

if (flag_nmea) %NMEA command
    port = 0; %current port
    NMEA_cmd = NMEA_PORZA_gen(port, baudrate, protocol);
    codeDEC = [uint8(NMEA_cmd)'; 13; 10];
else %BINR command
    port = '00'; %current port
    baudHEX = dec2hex(baudrate,8);
    codeDEC = hex2dec([DLE; BINR_0Bh; port; baudHEX(7:8); baudHEX(5:6); baudHEX(3:4); baudHEX(1:2); protocol; DLE; ETX]);
end

flushinput(serialObj);

% send message
% try
    fwrite(serialObj, codeDEC, 'uint8', 'async');
% catch
%     stopasync(serialObj);
%     fwrite(serialObj, codeDEC, 'uint8', 'async');
% end

[reply] = get_reply(serialObj);

%check reply
if (~isempty(strfind(reply',check')))
    out = 1;
    switch (format)
        case 'NMEA'
            set(serialObj,'Parity','none');
        otherwise %BINR
            set(serialObj,'Parity','odd');
    end
else
    out = 0;
end

end

%--------------------------------------------------------------------------

function [reply] = get_reply(serialObj)

reply = [];

%time acquisition
start_time = toc;

%maximum waiting time
dtMax = 12;

reply_1 = 0;
reply_2 = 0;

while (reply_1 ~= reply_2) || (reply_1 == 0)
    
    %time acquisition
    current_time = toc;
    
    %check if maximum waiting time is expired
    if (current_time - start_time > dtMax)
        return
    end
    
    % serial port checking
    reply_1 = get(serialObj, 'BytesAvailable');
    pause(3);
    reply_2 = get(serialObj, 'BytesAvailable');
    
    fprintf('.')
end

reply = fread(serialObj, reply_1, 'uint8');

end
