function [out] = nvs_raw_output_rate(serialObj, rate)

% SYNTAX:
%   [out] = nvs_raw_output_rate(serialObj, rate);
%
% INPUT:
%   serialObj = serial Object identifier
%   rate = raw measurement output rate
%          (supported output rate configurations 1 / 2 / 5 / 10 Hz)
%
% OUTPUT:
%   out = receiver reply
%
% DESCRIPTION:
%   Enable raw data output at the specified rate on NVS receivers.

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

%--------------------------------------------------------------------------
% NAVIGATION RATE
%--------------------------------------------------------------------------

%BINR format
DLE = '10'; %beginning of message
ETX = '03'; %end of message
DType = '02'; %data type (2 == navigation rate)
BINR_D7h = 'D7';
BINR_E7h = 'E7';
rateHEX = dec2hex(rate,2); %in Hz

codeHEX = [DLE; BINR_D7h; DType; rateHEX; DLE; ETX];
codeDEC = hex2dec(codeHEX);

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

check = hex2dec([DLE; BINR_E7h]);

[out1] = nvs_check_reply(serialObj, check);

%--------------------------------------------------------------------------
% RAW DATA RATE
%--------------------------------------------------------------------------

%BINR format
DLE = '10'; %beginning of message
ETX = '03'; %end of message
BINR_F4h = 'F4';
BINR_70h = '70';
rateHEX = dec2hex((1/rate)/0.1,2); %number of intervals of 100 ms

codeHEX = [DLE; BINR_F4h; rateHEX; DLE; ETX];
codeDEC = hex2dec(codeHEX);

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

check = hex2dec([DLE; BINR_70h]);

[out2] = nvs_check_reply(serialObj, check);

%--------------------------------------------------------------------------

if (out1 && out2)
    out = 1;
else
    out = 0;
end
