function [out] = fastrax_check_ACK(serialObj, tran_id, msg_id, payload)

% SYNTAX:
%   [out] = fastrax_check_ACK(serialObj, tran_id, msg_id, payload);
%
% INPUT:
%   serialObj = serial Object identifier
%   tran_id = transaction ID of the request [dec]
%   msg_id  = ID of the message that was sent [hex]
%   payload = payload of the message that was sent [dec array]
%
% OUTPUT:
%   out = acknowledge outcome
%
% DESCRIPTION:
%   Check acknowledge reply after sending Fastrax messages.

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

header1 = '3C';     % header (hexadecimal value) - 60 (decimal value) - ascii '<'
header2 = '21';     % header (hexadecimal value) - 33 (decimal value) - ascii '!'

source = '20';      % source information
destination = '42'; % destination information

% switch Response bit from 0 to 1
tran_id = dec2bin(tran_id,8);
tran_id(1) = '1';
tran_id = fbin2dec(tran_id);

payload(2)  = 64; % modify payload

pay_size = size(payload,1)/2; % payload size

payload_bin = dec2bin(payload,8)';
payload_bin = payload_bin(:)';

checksum = checksumFTX(payload_bin,pay_size); % checksum

end1 = '3E';

ACK_HEX = [header1; header2; source; destination; msg_id];
ACK_DEC = hex2dec(ACK_HEX);
ACK_DEC = [ACK_DEC; tran_id; pay_size; payload; fbin2dec(checksum(1:8)); fbin2dec(checksum(9:16)); hex2dec(end1)];

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

if (~isempty(index) & length(reply(index:end)) >= 10+pay_size*2)
    % extract acknowledge message from reply
    reply = reply(index:index+9+pay_size*2);

    if (reply == ACK_DEC) %#ok<BDSCI>
        out = 1;
    end
end