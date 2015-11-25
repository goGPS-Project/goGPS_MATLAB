function [data, nmea_sentences] = decode_ublox(msg, constellations, wait_dlg)

% SYNTAX:
%   [data, nmea_sentences] = decode_ublox(msg, constellations, wait_dlg);
%
% INPUT:
%   msg = binary message received by the u-blox receiver
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   data = cell-array that contains the decoded u-blox messages
%          (message class and id are in the first cell-array field)
%   nmea_sentences = cell-array containing all NMEA sentences found in input msg
%
% DESCRIPTION:
%   u-blox UBX messages decoding (also in sequence).

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

if (nargin < 2 || isempty(constellations))
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

%----------------------------------------------------------------------------------------------
% UBX MESSAGE HEADER
%----------------------------------------------------------------------------------------------

header1 = 'B5';      % header (hexadecimal value)
header2 = '62';      % header (hexadecimal value)

codeHEX = [header1 header2];              % initial hexadecimal stream
codeBIN = dec2bin(hex2dec(codeHEX),16);   % initial binary stream

pos_UBX = strfind(msg, codeBIN);          % message initial index

%----------------------------------------------------------------------------------------------
% NMEA MESSAGE HEADER
%----------------------------------------------------------------------------------------------

headerNMEA1 = '24';                      % NMEA header ($)
headerNMEA2 = '47';                      % NMEA header (G)
headerNMEA3 = '50';                      % NMEA header (P)

codeHEX_NMEA = [headerNMEA1 headerNMEA2 headerNMEA3];      % initial hexadecimal stream
codeBIN_NMEA = dec2bin(hex2dec(codeHEX_NMEA),24);          % initial binary stream

pos_NMEA = strfind(msg, codeBIN_NMEA);   % NMEA message initial index

%----------------------------------------------------------------------------------------------
% MESSAGE STARTING POINT
%----------------------------------------------------------------------------------------------

% output variable initialization
data = cell(0);
nmea_sentences = cell(0);
nmea_counter = 1;
nmea_string = '';

% find the index of the first message, if any
if (~isempty(pos_UBX) && ~isempty(pos_NMEA))

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

%----------------------------------------------------------------------------------------------
% MESSAGE DECODING LOOP
%----------------------------------------------------------------------------------------------

% counter initialization
i = 0;

if (nargin == 3)
    waitbar(0,wait_dlg,'Decoding rover stream...')
end

while (pos + 15 <= length(msg))

    if (nargin == 3)
        waitbar(pos/length(msg),wait_dlg)
    end

    % check if there is an UBX header
    if (strcmp(msg(pos:pos+15),codeBIN))

        % counter increment
        i = i + 1;

        % skip the u-blox header (16 bit)
        pos = pos + 16;

        if (pos + 31 <= length(msg))

            % message class (1 byte)
            class = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
            class = dec2hex(class,2);

            % message id (1 byte)
            id = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
            id = dec2hex(id,2);

            %to detect truncated messages (before computing LEN, because
            %sometimes the LEN field is truncated as well)
            pos_nxt = pos_UBX(find(pos_UBX>pos,1));
            pos_rem = pos_nxt-pos;
            
            % payload length (2 bytes)
            LEN1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
            LEN2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
            LEN = LEN1 + (LEN2 * 2^8);      % little endian

            if (LEN ~= 0)
                if (pos + 8*LEN + 15 <= length(msg))
                    
                    % checksum
                    CK_A = 0; CK_B = 0;
                    Nslices = floor((8*LEN + 31) / 8 + 1);    %pre-allocate to
                    slices = cell(Nslices,1);                 %increase speed
                    k = 1;
                    for j = (pos - 32) : 8 : (pos + 8*LEN - 1)
                        slices{k,1} = msg(j:j+7);
                        k = k + 1;
                    end
                    slices = fbin2dec(slices);         %call 'fbin2dec' only once (to optimize speed)
                    for r = 1 : k-1
                        CK_A = CK_A + slices(r,1);
                        CK_B = CK_B + CK_A;
                    end
                    CK_A = mod(CK_A,256);
                    CK_B = mod(CK_B,256);
                    CK_A_rec = fbin2dec(msg(pos + 8*LEN:pos + 8*LEN + 7));
                    CK_B_rec = fbin2dec(msg(pos + 8*LEN + 8:pos + 8*LEN + 15));

                    % if checksum matches
                    if (CK_A == CK_A_rec) && (CK_B == CK_B_rec)

                        % message identification
                        switch class
                            
                            % RXM (receiver manager)
                            case '02'
                                switch id
                                    % RAW (raw measurement)
                                    case '10', [data(:,i)] = decode_RXM_RAW(msg(pos:pos+8*LEN-1), constellations);
                                        
                                    % SFRB (subframe buffer)
                                    %case '11', [data(:,i)] = decode_RXM_SFRB(msg(pos:pos+8*LEN-1), constellations);
                                        
                                    % EPH (*OBSOLETE* ephemeris)
                                    case '31'
                                       if (LEN == 104) %(ephemeris available)
                                           [data(:,i)] = decode_RXM_EPH(msg(pos:pos+8*LEN-1), constellations);
                                       end
                                end
                                
                            % AID (aiding messages)
                            case '0B'
                                switch id
                                    % HUI (sat. Health / UTC / Ionosphere)
                                    case '02', [data(:,i)] = decode_AID_HUI(msg(pos:pos+8*LEN-1));
                                        
                                    % EPH (ephemeris)
                                    case '31'
                                        if (LEN == 104) %(ephemeris available)
                                            [data(:,i)] = decode_AID_EPH(msg(pos:pos+8*LEN-1), constellations);
                                        end
                                end
                        end
                        
                    else
                        warning('UBXDecoder:ChecksumError','checksum error (class: 0x%s, id: 0x%s))\n',class,id);
                        
                        %skip truncated messages
                        % +4 to include the two length bytes and the two checksum bytes
                        if (~isempty(pos_rem) && ~mod(pos_rem,8) && 8*(LEN+4) > pos_rem)
                            warning('UBXDecoder:TruncatedMessage','truncated UBX message detected and skipped (class: 0x%s, id: 0x%s))\n',class,id);
                            pos = pos_nxt;
                            continue
                        end
                    end
                    
                    % skip the message body
                    pos = pos + 8*LEN;
                    
                    % skip the 2 checksum bytes
                    pos = pos + 16;
                else
                    break
                end
            end
        else
            break
        end

    % check if a NMEA sentence is starting
    elseif (pos + 23 <= length(msg)) && (strcmp(msg(pos:pos+23),codeBIN_NMEA))
        
        % search for <CR><LF>
        % The maximum number of characters for a valid NMEA 0183 sentence
        % is 82, but in order not to miss invalid length NMEA sentences
        % (i.e. not standard), a maximum of 100 characters is used.
        % Thus the search for the end delimiter is restricted within
        % 100*8 = 800 bits or the end of the message whichever comes first.
        if ((length(msg)-pos)<799)
            pos_ENDNMEA = strfind(msg(pos:end),[dec2bin(13,8) dec2bin(10,8)]);
        else
            pos_ENDNMEA = strfind(msg(pos:pos+799),[dec2bin(13,8) dec2bin(10,8)]);
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
        pos = pos_UBX(find(pos_UBX>pos,1));
        if (isempty(pos)), break, end;
    end
end
