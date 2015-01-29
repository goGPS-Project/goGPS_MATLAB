function [data] = decode_rtcm2(msg, constellations, time_GPS)

% SYNTAX:
%   [data] = decode_rtcm2(msg, constellations, time_GPS)
%
% INPUT:
%   msg = binary message received from the master station
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%   time_GPS = current GPS time (in case message 14 is not available)
%
% OUTPUT:
%   data = cell-array that contains the decoded RTCM messages
%          (packet number is in the first cell-array field)
%
% DESCRIPTION:
%   RTCM 2.3 binary messages decoding (also in sequence).

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
    [constellations] = goGNSS.initConstellation(1, 1, 0, 0, 0, 0);
end

%----------------------------------------------------------------------------------------------
% RTCM2 HEADER DETECTION
%----------------------------------------------------------------------------------------------

preamble = '01100110';           %RTCM2 preamble
preamble_inv = '10011001';       %inverted RTCM2 preamble

check = ['0' preamble];          %preamble is good only if preceding bit is '0'
check_inv = ['1' preamble_inv];  %inverted preamble is good only if preceding bit is '1'

pos_preamble = strfind(msg, check) + 1;
pos_preamble_inv = strfind(msg, check_inv) + 1;

%----------------------------------------------------------------------------------------------
% FIRST MESSAGE SYNCHRONIZATION AND DECODING
%----------------------------------------------------------------------------------------------

%discard potential preambles for which the two preceding bits are not available
pos_preamble(pos_preamble < 3) = [];
pos_preamble_inv(pos_preamble_inv < 3) = [];

pos = [];

if (~isempty(pos_preamble) | ~isempty(pos_preamble_inv))
    sync = 0;
    i = 1;
    j = 1;

    while (~sync & i <= length(pos_preamble) & j <= length(pos_preamble_inv) &...
            pos_preamble(i)+29 < length(msg) & pos_preamble_inv(i)+29 < length(msg))

        parity = 0;
        parity_inv = 0;

        if (~isempty(pos_preamble(i)))
            [parity, decoded_word] = gps_parity(msg(pos_preamble(i)-2:pos_preamble(i)-1), msg(pos_preamble(i):pos_preamble(i)+29));
            i = i + 1;
        end
        if (~isempty(pos_preamble_inv(j)))
            [parity_inv, decoded_word_inv] = gps_parity(msg(pos_preamble_inv(j)-2:pos_preamble_inv(j)-1), msg(pos_preamble_inv(j):pos_preamble_inv(j)+29));
            j = j + 1;
        end

        if (parity & ~parity_inv)
            pos = pos_preamble(i-1);
            sync = 1;
        elseif (~parity & parity_inv)
            pos = pos_preamble_inv(j-1);
            parity = parity_inv;
            decoded_word = decoded_word_inv;
            sync = 1;
        elseif (parity & parity_inv)
            if (pos_preamble(i-1) < pos_preamble_inv(j-1))
                pos = pos_preamble(i-1);
                sync = 1;
            else
                pos = pos_preamble_inv(j-1);
                parity = parity_inv;
                decoded_word = decoded_word_inv;
                sync = 1;
            end
        end
    end

end

%----------------------------------------------------------------------------------------------
% MESSAGE IDENTIFICATION & DECODING LOOP
%----------------------------------------------------------------------------------------------

% output variable initialization
data = cell(0);

% if there is at least one RTCM messagge
if ~isempty(pos)

    % counter initialization
    i = 0;

    while (pos+59 <= length(msg) & (strcmp(msg(pos-1:pos+7),check) | strcmp(msg(pos-1:pos+7),check_inv)))

        % counter increment
        i = i + 1;

        % read the first word of the 2-word header
        if (i > 1)
            [parity, decoded_word] = gps_parity(msg(pos-2:pos-1), msg(pos:pos+29));
        end

        if (parity)
            % message type (6 bit)
            type = fbin2dec(decoded_word(9:14));
            % station number (10 bit)
            station = fbin2dec(decoded_word(15:24)); %#ok<NASGU>

            pos = pos + 30;

            % read the second word of the 2-word header
            [parity, decoded_word] = gps_parity(msg(pos-2:pos-1), msg(pos:pos+29));

            if (parity)
                % modified Z-count (13 bit)
                modz = fbin2dec(decoded_word(1:13));
                % seq (3 bit)
                seq = fbin2dec(decoded_word(14:16)); %#ok<NASGU>
                % number of words following
                n_words = fbin2dec(decoded_word(17:21));
                % reference station health
                health = fbin2dec(decoded_word(22:24)); %#ok<NASGU>

                pos = pos + 30;

                if (pos+n_words*30-1 <= length(msg))

                    % message identification
                    switch type

                        % GPS reference station parameters
                        case 3
                            [data(:,i)] = decode_3(msg(pos-2:pos+n_words*30-1), n_words);

                        % RTK uncorrected carrier phases
                        case 18
                            [data(:,i)] = decode_18(msg(pos-2:pos+n_words*30-1), n_words, modz, constellations);
                            if(nargin == 2)
                                sec_of_hour = data{2,i}(2);
                                sec_of_week = time_GPS - mod(time_GPS,3600) + sec_of_hour;
                                data{2,i}(2) = sec_of_week;
                            end

                        % RTK uncorrected pseudoranges
                        case 19
                            [data(:,i)] = decode_19(msg(pos-2:pos+n_words*30-1), n_words, modz, constellations);
                            if(nargin == 2)
                                sec_of_hour = data{2,i}(2);
                                sec_of_week = time_GPS - mod(time_GPS,3600) + sec_of_hour;
                                data{2,i}(2) = sec_of_week;
                            end
                    end

                    pos = pos + n_words*30;
                end
            else
                pos = pos + 30;
            end
        end

    end
end
