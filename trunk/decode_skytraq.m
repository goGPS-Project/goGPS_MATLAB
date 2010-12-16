function [data, nmea_sentences] = decode_skytraq(msg, wait_dlg)

% SYNTAX:
%   [data, nmea_sentences] = decode_skytraq(msg, wait_dlg);
%
% INPUT:
%   msg = binary message received by the SkyTraq receiver
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   data = cell-array that contains the decoded SkyTraq messages
%          (message class and id are in the first cell-array field)
%   nmea_sentences = cell-array containing all NMEA sentences found in input msg
%
% DESCRIPTION:
%   SkyTraq binary messages decoding (also in sequence).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
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
% MESSAGE HEADER
%----------------------------------------------------------------------------------------------

header1 = 'A0';      % header (hexadecimal value)
header2 = 'A1';      % header (hexadecimal value)

codeHEX = [header1 header2];              % initial hexadecimal stream
codeBIN = dec2bin(hex2dec(codeHEX),16);   % initial binary stream

pos_HDR = findstr(msg, codeBIN);          % message initial index

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
% MESSAGE STARTING POINT
%----------------------------------------------------------------------------------------------

% output variable initialization
data = cell(0);
nmea_sentences = cell(0);
nmea_counter = 1;
nmea_string = '';

% find the index of the first message, if any
if (~isempty(pos_HDR) & ~isempty(pos_NMEA))

    if (pos_HDR(1) < pos_NMEA(1))
        pos = pos_HDR(1);
    else
        pos = pos_NMEA(1);
    end

elseif (~isempty(pos_HDR))
    pos = pos_HDR(1);

elseif (~isempty(pos_NMEA))
    pos = pos_NMEA(1);

else
    return
end

%----------------------------------------------------------------------------------------------
% MESSAGE DECODING LOOP
%----------------------------------------------------------------------------------------------

% counter initialization
i = 0;

if (nargin == 2)
    waitbar(0,wait_dlg,'Decoding rover stream...')
end

while (pos + 15 <= length(msg))

    if (nargin == 2)
        waitbar(pos/length(msg),wait_dlg)
    end

    % check if there is an header
    if (strcmp(msg(pos:pos+15),codeBIN))

        % counter increment
        i = i + 1;

        % skip the header (16 bit)
        pos = pos + 16;

        if (pos + 23 <= length(msg))
            
            % payload length (2 bytes)
            LEN1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
            LEN2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
            LEN = LEN2 + (LEN1 * 2^8);      % big endian

            if (LEN ~= 0)
                if (pos + 8*LEN + 23 <= length(msg))
                    
                    % checksum
                    CS = 0;
                    slices = cell(1,LEN);                     % pre-allocate to increase speed
                    k = 1;
                    for j = pos : 8 : (pos + 8*LEN - 1)
                        slices{k} = msg(j:j+7);
                        k = k + 1;
                    end
                    slices = fbin2dec(slices);                %call 'fbin2dec' only once (to optimize speed)
                    for k = 1 : LEN
                        CS = bitxor(CS,slices(k));
                    end
                    CS_rec = fbin2dec(msg(pos + 8*LEN:pos + 8*LEN + 7));
                    
                    % if checksum matches
                    if (CS == CS_rec)
                        
                        % message id (1 byte)
                        id = fbin2dec(msg(pos:pos+7));
                        id = dec2hex(id,2);
                        
                        % message identification
                        switch id
                            % MEAS_TIME (Measurement time information)
                            case 'DC', fprintf('TIME\n');
                                
                            % RAW_MEAS (Raw channel measurements)
                            case 'DD', fprintf('RAW DATA\n');
                                
                            % SUBFRAME (Subframe buffer data)
                            case 'E0', fprintf('SUBFRAME\n');
                        end
                    else
                        %fprintf('Checksum error!\n');
                    end
                    
                    % skip the payload
                    pos = pos + 8*LEN;
                    
                    % skip the checksum byte and "end of sequence" bytes
                    pos = pos + 24;
                else
                    break
                end
            end
        else
            break
        end

    % check if a NMEA sentence is starting
    elseif (pos + 23 <= length(msg)) & (strcmp(msg(pos:pos+23),codeBIN_NMEA))
        
        % search for <CR><LF>
        % The maximum number of characters for a valid NMEA 0183 sentence
        % is 82, but in order not to miss invalid length NMEA sentences
        % (i.e. not standard), a maximum of 100 characters is used.
        % Thus the search for the end delimiter is restricted within
        % 100*8 = 800 bits or the end of the message whichever comes first.
        if ((length(msg)-pos)<799)
            pos_ENDNMEA = findstr(msg(pos:end),[dec2bin(13,8) dec2bin(10,8)]);
        else
            pos_ENDNMEA = findstr(msg(pos:pos+799),[dec2bin(13,8) dec2bin(10,8)]);
        end
        
        if ~isempty(pos_ENDNMEA)
            % save the NMEA sentence
            while (~strcmp(msg(pos:pos+7),'00001101'))
                nmea_string = [nmea_string char(fbin2dec(msg(pos:pos+7)))];
                pos = pos + 8;
            end
            
            % save just <LF> (without <CR>, otherwise MATLAB fails in interpreting it)
            pos = pos + 8;
            nmea_string = [nmea_string char(fbin2dec(msg(pos:pos+7)))];
            pos = pos + 8;
            
            nmea_sentences{nmea_counter,1} = nmea_string;
            nmea_counter = nmea_counter + 1;
            nmea_string = '';
        else
            % if a NMEA sentence is started but its end is not available,
            % just jump over the header and continue
            pos = pos + 24;
        end

    % check if there are other packages
    else

        % find the index of the first message, if any
        pos_HDR = findstr(msg(pos:end),codeBIN);
        if ~isempty(pos_HDR)
            pos = pos + pos_HDR(1) - 1;
        else
            break
        end
    end

end