function [data] = decode_3(msg, n_words)

% SYNTAX:
%   [data] = decode_3(msg, n_words)
%
% INPUT:
%   msg = binary message received from the master station
%   n_words = number of words composing the message
%
% OUTPUT:
%   data = cell-array that contains the message '3' information
%          1.1) message type
%          2.1) master station antenna X coordinate in meters (ECEF)
%          2.2) master station antenna Y coordinate in meters (ECEF)
%          2.3) master station antenna Z coordinate in meters (ECEF)
%
% DESCRIPTION:
%   RTCM 2 format, message '3' decoding.

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

%message pointer initialization (after preceding bytes)
pos = 3;

%variable initialization
parity = zeros(4,1);

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(3,1);

if (pos + n_words*30-1 <= length(msg))

    for i = 1 : n_words
        [parity(i), decoded_word(i,1:24)] = gps_parity(msg(pos-2:pos-1), msg(pos:pos+29));
        pos = pos + 30;
    end

    if (parity)

        %message type = 3
        type = 3;

        %antenna X coordinate
        Xcoord = twos_complement([decoded_word(1,:) decoded_word(2,1:8)]);

        %antenna Y coordinate
        Ycoord = twos_complement([decoded_word(2,9:24) decoded_word(3,1:16)]);

        %antenna Z coordinate
        Zcoord = twos_complement([decoded_word(3,17:24) decoded_word(4,1:24)]);

        %output data save
        data{1} = type;
        data{2}(1) = Xcoord * 0.01;
        data{2}(2) = Ycoord * 0.01;
        data{2}(3) = Zcoord * 0.01;

    end
end