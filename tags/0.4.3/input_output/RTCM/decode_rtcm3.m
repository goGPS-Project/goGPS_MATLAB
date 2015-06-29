function [data] = decode_rtcm3(msg, constellations, wait_dlg)

% SYNTAX:
%   [data] = decode_rtcm3(msg, constellations, wait_dlg);
%
% INPUT:
%   msg = binary message received from the master station
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   data = cell-array that contains the decoded RTCM messages
%          (packet number is in the first cell-array field)
%
% DESCRIPTION:
%   RTCM 3.1 binary messages decoding (also in sequence).

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

if (nargin < 2 || isempty(constellations))
    [constellations] = goGNSS.initConstellation(1, 1, 0, 0, 0, 0);
end

%----------------------------------------------------------------------------------------------
% MESSAGE ID CODE
%----------------------------------------------------------------------------------------------

preamble = '11010011';      % FIXED transport layer header (8 bit)
%reserved = '000000';        % reserved field (6 bit). It could change in the future!

%codeBIN = [preamble reserved];      % binary initial stream

pos_all = strfind(msg, preamble);    % message initial index

%----------------------------------------------------------------------------------------------
% MESSAGE IDENTIFICATION
%----------------------------------------------------------------------------------------------

% output variable initialization
data = cell(0);

% if there is at least one RTCM messagge
if ~isempty(pos_all)

    % counter initialization
    i = 0;
    
    if (nargin == 3)
        waitbar(0,wait_dlg,'Decoding master stream...')
    end

    % pointer initialization
    pos0 = pos_all(1);

    while ~isempty(pos0)
        
        % skip the "preamble" (8 bit) and "reserved" (6 bit) fields
        pos = pos0 + 14;

        if (pos + 9 <= length(msg))
            
            if (nargin == 3)
                waitbar(pos/length(msg),wait_dlg)
            end

            % message length (10 bit)
            LEN = fbin2dec(msg(pos:pos+9));  pos = pos + 10;

            if (pos + 8*LEN + 23 <= length(msg))

                % CRC check
                CRC_comp = crc24q(msg(pos0 : pos + 8*LEN - 1));
                CRC_read = msg(pos + 8*LEN : pos + 8*LEN + 23);

                if (strcmp(CRC_comp, CRC_read))

                    % counter increment
                    i = i + 1;

                    % message number (12 bit)
                    DF002 = fbin2dec(msg(pos:pos+11));

                    % message identification
                    switch DF002

                        % GPS observations on L1 carrier
                        case 1001
                            [data(:,i)] = decode_1001(msg(pos:pos+8*LEN-1), constellations);

                        % GPS observations on L1 carrier
                        case 1002
                            [data(:,i)] = decode_1002(msg(pos:pos+8*LEN-1), constellations);

                        % GPS observations on L1/L2 carrier
                        case 1003
                            [data(:,i)] = decode_1003(msg(pos:pos+8*LEN-1), constellations);

                        % GPS observations on L1/L2 carrier
                        case 1004
                            [data(:,i)] = decode_1004(msg(pos:pos+8*LEN-1), constellations);

                        % master station coordinates
                        case 1005
                            [data(:,i)] = decode_1005(msg(pos:pos+8*LEN-1));

                        % master station coordinates + antenna height
                        case 1006
                            [data(:,i)] = decode_1006(msg(pos:pos+8*LEN-1));

                        % antenna description
                        case 1007
                            [data(:,i)] = decode_1007(msg(pos:pos+8*LEN-1));

                        % antenna description + serial number
                        case 1008
                            [data(:,i)] = decode_1008(msg(pos:pos+8*LEN-1));

                        % GLONASS observations on L1 carrier
                        case 1009
                            [data(:,i)] = decode_1009(msg(pos:pos+8*LEN-1), constellations);

                        % GLONASS observations on L1 carrier
                        case 1010
                            [data(:,i)] = decode_1010(msg(pos:pos+8*LEN-1), constellations);

                        % GLONASS observations on L1/L2 carrier
                        case 1011
                            [data(:,i)] = decode_1011(msg(pos:pos+8*LEN-1), constellations);

                        % GLONASS observations on L1/L2 carrier
                        case 1012
                            [data(:,i)] = decode_1012(msg(pos:pos+8*LEN-1), constellations);

                        % system parameters
                        case 1013
                            [data(:,i)] = decode_1013(msg(pos:pos+8*LEN-1), constellations);

                        % auxiliary network information
                        case 1014
                            [data(:,i)] = decode_1014(msg(pos:pos+8*LEN-1));

                        % ionospheric correction differences
                        case 1015
                            [data(:,i)] = decode_1015(msg(pos:pos+8*LEN-1), constellations);

                        % geometric correction differences
                        case 1016
                            [data(:,i)] = decode_1016(msg(pos:pos+8*LEN-1), constellations);

                        % combined ionospheric and geometric correction differences
                        case 1017
                            [data(:,i)] = decode_1017(msg(pos:pos+8*LEN-1), constellations);

                        %GPS ephemerides
                        case 1019
                            [data(:,i)] = decode_1019(msg(pos:pos+8*LEN-1), constellations);

                        %GLONASS ephemerides
                        case 1020
                            [data(:,i)] = decode_1020(msg(pos:pos+8*LEN-1), constellations);

                        % not implemented
                        case 1029
                            [data(:,i)] = decode_1029(msg(pos:pos+8*LEN-1));
                            
                    end

                    % move pointer
                    pos = pos + 8*LEN;

                    % skip the 24 bits of CRC
                    pos = pos + 24;

                end
            end
        end

        % message ending byte
        while (mod(pos,8) ~= 1)
            pos = pos + 1;
        end

        % pointer position
        pos0 = pos_all(find(pos_all>=pos,1));

    end
end