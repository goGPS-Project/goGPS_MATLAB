function [out] = italk_PPS_MEAS_MS(serialObj, tran_id, offset)

% SYNTAX:
%   [out] = italk_PPS_MEAS_MS(serialObj, tran_id, offset);
%
% INPUT:
%   serialObj = serial Object identifier
%   tran_id = transaction ID
%   offset = measurement time offset from GPS TOW [ms]
%
% OUTPUT:
%   out = receiver reply
%
% DESCRIPTION:
%   ITALK configuration message.
%
%   "Target measurement time offset [ms] : If the PPS_SYNC_TRACK flag is
%   YES, then the tracking measurements are aligned with the GPS tow with
%   an offset specified by this parameter. A value of 500 would result in
%   measurements every GPS tow + 0.5 seconds."
%   (http://isuite.fastrax.fi/sdk/341/system/sys_parameters.html)

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

%time offset in hexadecimal format
offsetHEX = dec2hex(offset,4);

payload  = ['01';'00';'01';'00';'44';'24';offsetHEX(3:4);offsetHEX(1:2)]; % payload (set offset measurement offset from GPS TOW)

pay_size = size(payload,1)/2;               % payload size

payload_bin = dec2bin(hex2dec(payload),8)';
payload_bin = payload_bin(:)';
payload_dec = hex2dec(payload);

checksum = checksumFTX(payload_bin,pay_size);    % checksum

end1 = '3E';

codeHEX = [header1; header2; source; destination; msg_id];
codeDEC = hex2dec(codeHEX);
codeDEC = [codeDEC; tran_id; pay_size; payload_dec; fbin2dec(checksum(1:8)); fbin2dec(checksum(9:16)); hex2dec(end1) ];

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

[out] = fastrax_check_ACK(serialObj, tran_id, msg_id, payload_dec);
