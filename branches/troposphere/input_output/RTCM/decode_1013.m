function [data] = decode_1013(msg, constellations)

% SYNTAX:
%   [data] = decode_1013(msg, constellations)
%
% INPUT:
%   msg = binary message received from the master station
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the 1013 packet information
%          1.1)  DF002 = message number = 1013
%          2.1)  DF003 = reference station id
%          2.2)  DF051 = Modified Julian Day (MJD)
%          2.3)  DF052 = Seconds of Day (UTC)
%          2.4)  DF053 = No. of Message ID Announcements to Follow
%          2.5)  DF054 = Leap Seconds, GPS-UTC
%          3.1)  DF055 = Message ID
%          3.2)  DF056 = Message SyncFlag (0 = Asynchronous, 1 = Synchronous)
%          3.3)  DF057 = Message transmission interval
%
% DESCRIPTION:
%   RTCM format 1013 message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Sara Lucca
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
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(5,1);
data{3} = zeros(constellations.nEnabledSat,3);

%message number = 1013
DF002 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%reference station id
DF003 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%Modified Julian Day (MJD)
DF051 = fbin2dec(msg(pos:pos+15));  pos = pos + 16;

%UTC (sec)
DF052 = fbin2dec(msg(pos:pos+16));  pos = pos + 17;

%Number of Message ID Announcements to Follow
DF053 = fbin2dec(msg(pos:pos+4));  pos = pos + 5;

%Leap Seconds,GPS-UTC
DF054 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

%output data save
data{1} = DF002;
data{2}(1) = DF003;
data{2}(2) = DF051;
data{2}(3) = DF052;
data{2}(4) = DF053;
data{2}(5) = DF054;

%-------------------------------------------------

%satellite number
NMF = data{2}(4);

%data decoding for every satellite
for i = 1 : NMF

    %message number
    DF055 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

    %Sync flag
    DF056 = fbin2dec(msg(pos));  pos = pos + 1;

    %Transmission interval
    DF057 = fbin2dec(msg(pos:pos+15)) * 0.1;  pos = pos + 16;

    %------------------------------------------------
    
    % assign constellation-specific indexes
    idx = [];
    if (i <= 24 && constellations.GLONASS.enabled)
        idx = constellations.GLONASS.indexes(i);
    end
    
    %output data save
    data{3}(idx,1)  = DF055;
    data{3}(idx,2)  = DF056;
    data{3}(idx,3)  = DF057;

end