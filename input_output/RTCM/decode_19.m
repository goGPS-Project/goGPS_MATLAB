function [data] = decode_19(msg, n_words, modz, constellations)

% SYNTAX:
%   [data] = decode_19(msg, n_words, modz, constellations)
%
% INPUT:
%   msg = binary message received from the master station
%   n_words = number of words composing the message
%   modz = modified Z-count
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the message '18' information
%          1.1) message type
%          2.1) Frequency indicator (L1 or L2)
%          2.2) time of measurement
%          3.1) Multiple message indicator
%          3.2) C/A - P code indicator
%          3.3) GPS/GLONASS satellite constellation indicator
%          3.4) Satellite ID
%          3.5) Data quality
%          3.6) Cumulative loss of continuity indicator
%          3.7) Carrier phase
%
% DESCRIPTION:
%   RTCM 2 format, message '19' decoding.

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

%message pointer initialization
pos = 3;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(2,1);
data{3} = zeros(constellations.nEnabledSat,7);

if (~constellations.GPS.enabled)
    return
end

if (pos + n_words*30-1 <= length(msg))

    for i = 1 : n_words
        [parity(i), decoded_word(i,1:24)] = gps_parity(msg(pos-2:pos-1), msg(pos:pos+29));
        pos = pos + 30;
    end

    if (parity)

        %message type = 19
        type = 19;

        %frequency indicator
        F = fbin2dec(decoded_word(1,1:2));

        %time of measurement (seconds of the hour)
        time = fbin2dec(decoded_word(1,5:24));

        %output data save
        data{1} = type;
        data{2}(1) = F;
        data{2}(2) = 0.6*modz + time*1e-6;

        for i = 1 : (n_words-1)/2

            %multiple message indicator
            MMI = fbin2dec(decoded_word(i*2,1));

            %C/A - P code indicator
            CI = fbin2dec(decoded_word(i*2,2));

            %GPS/GLONASS satellite constellation indicator
            sys = fbin2dec(decoded_word(i*2,3));

            %satellite ID
            SV = fbin2dec(decoded_word(i*2,4:8));
            if (SV == 0); SV=32; end

            %data quality
            DQ = fbin2dec(decoded_word(i*2,9:11));

            %cumulative loss of continuity indicator
            loss = fbin2dec(decoded_word(i*2,12:16));

            %carrier phase
            pr = fbin2dec([decoded_word(i*2,17:24) decoded_word(i*2+1,1:24)]);

            %data output save
            data{3}(SV,1) = MMI;
            data{3}(SV,2) = CI;
            data{3}(SV,3) = sys;
            data{3}(SV,4) = SV;
            data{3}(SV,5) = DQ;
            data{3}(SV,6) = loss;
            data{3}(SV,7) = pr * 0.02;
        end
    end
end