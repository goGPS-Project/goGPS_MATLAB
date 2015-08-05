function [data] = decode_1012(msg, constellations)

% SYNTAX:
%   [data] = decode_1012(msg, constellations)
%
% INPUT:
%   msg = binary message received from the master station
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the 1012 packet information
%          1.1)  DF002 = message number = 1012
%          2.1)  DF003 = reference station id
%          2.2)  (DF034 / 1000) = UTC time + 3 hours, in seconds (GLONASS epoch)
%          2.3)  DF005 = are there other messages of the same epoch? YES=1, NO=0
%          2.4)  DF035 = number of visible satellites
%          2.5)  DF036 = phase-smoothed code? YES=1, NO=0
%          2.6)  DF037 = smoothing window
%          3.1)  DF039 = code type vector on L1: C/A=0, P=1
%          3.2)  (DF041 * 0.02) + (DF044 * 599584.92) = code observation vector on L1
%          3.3)  (code observation L1 + (DF042*0.0005)) / lambda1 = phase observation vector on L1
%          3.4)  DF043 = how long L1 has been locked? index vector (cycle-slip=0)
%          3.5)  (DF045 * 0.25) = signal-to-noise ratio vector in dBHz on L1(from 0 to 63.75 dBHz)
%          3.6)  DF046 = vector code type on L2: C/A=0, P=1, reserved=2,3
%          3.7)  ((DF041 + DF047) * 0.02) + (DF044 * 599584.92) = code observation vector on L2
%          3.8)  (code observation L1 + (DF048*0.0005)) / lambda2 = phase observation vector on L2
%          3.9)  DF049 = how long L2 has been locked? index vector (cycle-slip=0)
%          3.10) (DF050 * 0.25) = signal-to-noise ratio vector in dBHz on L2(from 0 to 63.75 dBHz)
%          3.11) (DF040 - 7) * 0.5625 + 1602.0 = frequency vector on L1
%          3.12) (DF040 - 7) * 0.4375 + 1246.0 = frequency vector on L2
%
% DESCRIPTION:
%   RTCM format 1012 message decoding.

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
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(6,1);
data{3} = zeros(constellations.nEnabledSat,12);

if (~constellations.GLONASS.enabled)
    return
end

%message number = 1012
DF002 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%reference station id
DF003 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%TOD = time of day in milliseconds
DF034 = fbin2dec(msg(pos:pos+26));  pos = pos + 27;

%other synchronous RTCM messages flag (YES=1, NO=0)
DF005 = fbin2dec(msg(pos));  pos = pos + 1;

%number of visible satellites
DF035 = fbin2dec(msg(pos:pos+4));  pos = pos + 5;

%phase-smoothed code flag (YES=1, NO=0)
DF036 = fbin2dec(msg(pos));  pos = pos + 1;

%smoothing window
DF037 = fbin2dec(msg(pos:pos+2));  pos = pos + 3;

%output data save
data{1} = DF002;
data{2}(1) = DF003;
data{2}(2) = DF034 / 1000;
data{2}(3) = DF005;
data{2}(4) = DF035;
data{2}(5) = DF036;
data{2}(6) = DF037;

%-------------------------------------------------

%satellite number
NSV = data{2}(4);

%data decoding for every satellite
for i = 1 : NSV

    %analyzed satellite number
    SV = fbin2dec(msg(pos:pos+5));  pos = pos + 6;

    %if GLONASS satellite (known slot)
    if (SV >= 1 & SV <= 24)

        %code type (C/A=0, P=1)
        DF039 = fbin2dec(msg(pos));  pos = pos + 1;

        %frequency indicator
        DF040 = fbin2dec(msg(pos:pos+4));  pos = pos + 5;

        %L1 pseudorange
        DF041 = fbin2dec(msg(pos:pos+24));  pos = pos + 25;

        %L1 phaserange - L1 pseudorange
        DF042 = twos_complement(msg(pos:pos+19));  pos = pos + 20;

        %L1 lock-time index (see Table 4.3-2 on RTCM manual)
        DF043 = fbin2dec(msg(pos:pos+6));  pos = pos + 7;

        %L1 pseudorange integer ambiguity
        DF044 = fbin2dec(msg(pos:pos+6));  pos = pos + 7;

        %L1-CNR (carrier-to-noise ratio) apart from resolution
        DF045 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

        %---------------------------------------------------------

        %L2 code type (C/A=0, P=1,2,3)
        DF046 = fbin2dec(msg(pos:pos+1));  pos = pos + 2;

        %L2-L1 pseudorange
        DF047 = twos_complement(msg(pos:pos+13));  pos = pos + 14;

        %L2 phaserange - L1 pseudorange
        DF048 = twos_complement(msg(pos:pos+19));  pos = pos + 20;

        %L2lock-time index (see Table 4.3-2 on RTCM manual)
        DF049 = fbin2dec(msg(pos:pos+6));  pos = pos + 7;

        %L2-CNR (carrier-to-noise ratio) apart from resolution
        DF050 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

        %---------------------------------------------------------
        
        % assign constellation-specific indexes
        if (constellations.GLONASS.enabled)

            idx = constellations.GLONASS.indexes(SV);
            
            %L1 carrier frequency [MHz]
            data{3}(idx,11) = (DF040 - 7) * 0.5625 + 1602.0;
            
            %L2 carrier frequency [MHz]
            data{3}(idx,12) = (DF040 - 7) * 0.4375 + 1246.0;
            
            %debugging
            %v_light / (data{3}(idx,11) * 1e6)
            %v_light / (data{3}(idx,12) * 1e6)
            
            %---------------------------------------------------------
            
            %output data save
            data{3}(idx,1)  = DF039;
            data{3}(idx,2)  = (DF041 * 0.02)  + (DF044 * 599584.92);
            data{3}(idx,3)  = (data{3}(idx,2) + (DF042 * 0.0005)) * data{3}(idx,11) * 1e6 / goGNSS.V_LIGHT;
            data{3}(idx,4)  = DF043;
            data{3}(idx,5)  = DF045 * 0.25;
            data{3}(idx,6)  = DF046;
            data{3}(idx,7)  = (data{3}(idx,2) + (DF047 * 0.02));
            data{3}(idx,8)  = (data{3}(idx,2) + (DF048 * 0.0005)) * data{3}(idx,12) * 1e6 / goGNSS.V_LIGHT;
            data{3}(idx,9)  = DF049;
            data{3}(idx,10) = DF050 * 0.25;
        end

    else %SBAS satellites

        %do not store SBAS satellite information
        pos = pos + 124;

    end

end