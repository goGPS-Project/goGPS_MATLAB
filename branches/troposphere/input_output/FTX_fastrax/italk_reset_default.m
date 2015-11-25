function italk_reset_default(serialObj, tran_id)

% SYNTAX:
%   italk_reset_default(serialObj, tran_id);
%
% INPUT:
%   serialObj = serial Object identifier
%   tran_id = transaction ID
%
% DESCRIPTION:
%   ITALK configuration message.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Ivan Reguzzoni
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

%----------------------------------------------------------------------------------------------
% FTX MESSAGE HEADER
%----------------------------------------------------------------------------------------------

header1 = '3C';      % header (hexadecimal value) - 60 (decimal value) - ascii '<'
header2 = '21';      % header (hexadecimal value) - 33 (decimal value) - ascii '!'

source = '42';       % source information
destination = '20';  % destination information

msg_id   = '47';     % message id

payload  = ['03';'00';'00';'00';'00';'00']; % payload (reset all parameters to default values)

pay_size = size(payload,1)/2;               % payload size

payload_bin = dec2bin(hex2dec(payload),8)';
payload_bin = payload_bin(:)';

checksum = checksumFTX(payload_bin,pay_size);    % checksum

end1 = '3E';

codeHEX = [header1; header2; source; destination; msg_id];
codeDEC = hex2dec(codeHEX);
codeDEC = [codeDEC; tran_id; pay_size; hex2dec(payload); fbin2dec(checksum(1:8)); fbin2dec(checksum(9:16)); hex2dec(end1) ];

%serial port checking
reply_1 = get(serialObj,'BytesAvailable');

if (reply_1 ~= 0)
    %clear the serial port (data not decoded)
    reply = fread(serialObj, reply_1, 'uint8'); %#ok<NASGU>
end

% send message
fwrite(serialObj, codeDEC, 'uint8', 'async');

% wait for the asynchronous writing to finish
pause(0.5)
