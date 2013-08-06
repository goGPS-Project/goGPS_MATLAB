function [data] = decode_nvs(msg, wait_dlg)

% SYNTAX:
%   [data] = decode_nvs(msg, wait_dlg);
%
% INPUT:
%   msg = binary message received by the NVS receiver
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   data = cell-array that contains the decoded NVS messages
%          (message class and id are in the first cell-array field)
%
% DESCRIPTION:
%   NVS messages decoding (also in sequence).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.1 beta
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Daisuke Yoshida
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

warning off

% output variable initialization
data = cell(0);

%----------------------------------------------------------------------------------------------
% NVS MESSAGE HEADER
%----------------------------------------------------------------------------------------------

codeHEX = '10';                          % <DLE> header (hexadecimal value)
codeBIN = dec2bin(hex2dec(codeHEX),8);   % <DLE> header (binary value)

pos_all = findstr(msg, codeBIN);          % key sequence initial indexes

pos = 1;

for i = 1 : length(pos_all)
    next_byte = dec2hex(bin2dec(msg(pos_all(i)+7)));
    if (strcmp(next_byte,'F5') || ...
        strcmp(next_byte,'F7') || ...
        strcmp(next_byte,'4A') || ...
        strcmp(next_byte,'62'))
        
        pos = pos_all(i);
        
        break
    end
end

codeBIN_HDR = codeBIN;

%----------------------------------------------------------------------------------------------
% NVS MESSAGE FOOTER
%----------------------------------------------------------------------------------------------

header1 = '10';      % header (hexadecimal value)
header2 = '03';      % header (hexadecimal value)

codeHEX = [header1 header2];              % initial hexadecimal stream
codeBIN = dec2bin(hex2dec(codeHEX),16);   % initial binary stream

pos_FTR = findstr(msg, codeBIN);          % message initial index

i = 1;
while (i <= length(pos_FTR) && pos_FTR(i)+23 <= length(msg))
    if (~strcmp(msg(pos_FTR(i)+16:pos_FTR(i)+23),codeBIN_HDR))
        pos_FTR(i) = [];
    end
    i = i + 1;
end

%----------------------------------------------------------------------------------------------
% MESSAGE DECODING LOOP
%----------------------------------------------------------------------------------------------

% counter initialization
i = 1;

if (nargin == 2)
    waitbar(0,wait_dlg,'Decoding rover stream...')
end

while (pos + 7 < length(msg) && i <= length(pos_FTR))
    
    if (nargin == 2)
        waitbar(pos/length(msg),wait_dlg)
    end
    
    % check if there is an NVS header (<DLE>)
    if (strcmp(msg(pos:pos+7),codeBIN_HDR))
        
        % skip the NVS header (8 bit)
        pos = pos + 8;
        
        % message id (1 byte)
        id = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
        id = dec2hex(id,2);
        
        % position of the last bit of the data message (i.e. before <DLE><ETX>)
        data_msg_end = pos_FTR(i) - 1;

        % counter increment
        i = i + 1;
        
        switch id
            % RAW (raw measurement)
            case 'F5', [data(:,i)] = decode_F5h(msg(pos:data_msg_end));
                
            % HUI (sat. Health / UTC / Ionosphere)
            case '4A', [data(:,i)] = decode_4Ah(msg(pos:data_msg_end));
                
            % EPH (ephemeris)
            case 'F7', [data(:,i)] = decode_F7h(msg(pos:data_msg_end));
        end
        
        % skip the message body
        pos = data_msg_end + 1;
        
        % skip the 2 closing bytes (i.e. <DLE><ETX>)
        pos = pos + 16;
    else
        break
    end
end
