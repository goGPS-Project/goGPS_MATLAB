function [data] = decode_1009(msg, constellations)

% SYNTAX:
%   [data] = decode_1009(msg, constellations)
%
% INPUT:
%   msg = binary message received from the master station
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the 1009 packet information
%          1.1) DF002 = message number = 1009
%          2.1) DF003 = reference station id
%          2.2) (DF034 / 1000) = UTC time + 3 hours, in seconds (GLONASS epoch)
%          2.3) DF005 = are there other messages of the same epoch? YES=1, NO=0
%          2.4) DF035 = number of visible satellites
%          2.5) DF036 = phase-smoothed code? YES=1, NO=0
%          2.6) DF037 = smoothing window
%          3.1) DF039 = code type vector on L1: C/A=0, P=1
%          3.2) (DF041 * 0.02) = code observation vector on L1
%          3.3) (code observation L1 + (DF042*0.0005)) / lambda1 = phase observation vector on L1
%          3.4) DF043 = how long L1 has been locked? index vector (cycle-slip=0)
%          3.5) (DF040 - 7) * 0.5625 + 1602.0 = frequencies vector on L1
%
% DESCRIPTION:
%   RTCM format 1009 message decoding.

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

%message number = 1009
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

%number of satellites
NSV = data{2}(4);

%satellite data decoding
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

        %---------------------------------------------------------
        
        % assign constellation-specific indexes
        if (constellations.GLONASS.enabled)
            idx = constellations.GLONASS.indexes(SV);
            
            %carrier L1 frequency [MHz]
            data{3}(idx,5) = (DF040 - 7) * 0.5625 + 1602.0;
            
            %debugging
            %v_light / (data{3}(SV,6) * 1e6)
            
            %---------------------------------------------------------
            
            %output data save
            data{3}(idx,1)  = DF039;
            data{3}(idx,2)  = (DF041 * 0.02);
            data{3}(idx,3)  = (data{3}(idx,2) + (DF042 * 0.0005)) * data{3}(idx,6) * 1e6 / goGNSS.V_LIGHT;
            data{3}(idx,4)  = DF043;
        end

    else %SBAS satellites

        %do not store SBAS satellite information
        pos = pos + 58;

    end

end