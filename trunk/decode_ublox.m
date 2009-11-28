function [data, nmea_string] = decode_ublox(msg)

% SYNTAX:
%   [data, nmea_string] = decode_ublox(msg)
%
% INPUT:
%   msg = binary message received by the u-blox receiver
%
% OUTPUT:
%   data = cell-array that contains the decoded u-blox messages
%          (message class and id are in the first cell-array field)
%   nmea_string = string containing all NMEA sentences found in input msg
%
% DESCRIPTION:
%   u-blox UBX messages decoding (also in sequence).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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
% UBX MESSAGE HEADER
%----------------------------------------------------------------------------------------------

header1 = 'B5';      % header (hexadecimal value)
header2 = '62';      % header (hexadecimal value)

codeHEX = [header1 header2];              % initial hexadecimal stream
codeBIN = dec2bin(hex2dec(codeHEX),16);   % initial binary stream

pos_UBX = findstr(msg, codeBIN);          % message initial index

%----------------------------------------------------------------------------------------------
% NMEA MESSAGE HEADER
%----------------------------------------------------------------------------------------------

headerNMEA1 = '24';                      % NMEA header ($)
headerNMEA2 = '47';                      % NMEA header (G)
headerNMEA3 = '50';                      % NMEA header (P)

codeHEX_NMEA = [headerNMEA1 headerNMEA2 headerNMEA3];      % initial hexadecimal stream
codeBIN_NMEA = dec2bin(hex2dec(codeHEX_NMEA),24);          % initial binary stream

pos_NMEA = findstr(msg, codeBIN_NMEA);   % NMEA message initial index

%----------------------------------------------------------------------------------------------
% MESSAGE IDENTIFICATION
%----------------------------------------------------------------------------------------------

% output variable initialization
data = cell(0);

nmea_string = '';

% find the index of the first message, if any
if (~isempty(pos_UBX) & ~isempty(pos_NMEA))
    
    if (pos_UBX(1) < pos_NMEA(1))
        pos = pos_UBX(1);
    else
        pos = pos_NMEA(1);
    end
    
elseif (~isempty(pos_UBX))
    pos = pos_UBX(1);

elseif (~isempty(pos_NMEA))
    pos = pos_NMEA(1);

else
    return
end
    
    % counter initialization
    i = 0;

    while (pos + 15 <= length(msg))

        % check if there is an UBX header
        if (strcmp(msg(pos:pos+15),codeBIN))

            % counter increment
            i = i + 1;

            % skip the u-blox header (16 bit)
            pos = pos + 16;

            if (pos + 31 <= length(msg))

                % message class (1 byte)
                class = bin2dec(msg(pos:pos+7));  pos = pos + 8;
                class = dec2hex(class,2);

                % message id (1 byte)
                id = bin2dec(msg(pos:pos+7));  pos = pos + 8;
                id = dec2hex(id,2);

                % payload length (2 bytes)
                LEN1 = bin2dec(msg(pos:pos+7));  pos = pos + 8;
                LEN2 = bin2dec(msg(pos:pos+7));  pos = pos + 8;
                LEN = LEN1 + (LEN2 * 2^8);      % little endian
                clear LEN1 LEN2

                if (pos + 8*LEN - 1 <= length(msg))

                    % message identification
                    switch class

                        % RXM (receiver manager)
                        case '02'
                            switch id
                                % RAW (raw measurement)
                                case '10', [data(:,i)] = decode_RXM_RAW(msg(pos:pos+8*LEN-1));
                                % EPH (ephemerides)
                                case '31'
                                    if (LEN == 104) %(ephemerides available)
                                        [data(:,i)] = decode_RXM_EPH(msg(pos:pos+8*LEN-1));
                                    end
                            end
                    end

                end

                pos = pos + 8*LEN;

            else
                break
            end

            % skip the 2 checksum bytes
            pos = pos + 16;

        % check if a NMEA message is starting
        elseif (pos + 23 <= length(msg)) & (strcmp(msg(pos:pos+23),codeBIN_NMEA))

            % save the NMEA message (search for <CR>)
            while (pos + 7 <= length(msg)) & (bin2dec(msg(pos:pos+7)) ~= 13)
                nmea_string = [nmea_string char(bin2dec(msg(pos:pos+7)))];
                pos = pos + 8;
            end

            % save just <LF> (without <CR>, otherwise MATLAB fails in interpreting it)
            pos = pos + 8;
            nmea_string = [nmea_string char(bin2dec(msg(pos:pos+7)))];
            pos = pos + 8;

        else
            break
        end

    end
end