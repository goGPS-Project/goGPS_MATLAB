function [data] = decode_1004(msg, constellations)

% SYNTAX:
%   [data] = decode_1004(msg, constellations)
%
% INPUT:
%   msg = binary message received from the master station
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the 1004 packet information
%          1.1)  DF002 = message number = 1004
%          2.1)  DF003 = reference station id
%          2.2)  (DF004 / 1000) = time of week in seconds (GPS epoch)
%          2.3)  DF005 = are there other messages of the same epoch? YES=1, NO=0
%          2.4)  DF006 = number of visible satellites
%          2.5)  DF007 = phase-smoothed code? YES SI=1, NO=0
%          2.6)  DF008 = smoothing window
%          3.1)  DF010 = code type vector on L1: C/A=0, P=1
%          3.2)  (DF011 * 0.02) + (DF014 * 299792.458) = code observation vector on L1
%          3.3)  (code observation L1 + (DF012*0.0005)) / lambda1 = phase observation vector on L1
%          3.4)  DF013 = how long L1 has been locked? index vector (cycle-slip=0)
%          3.5)  (DF015 * 0.25) = signal-to-noise ratio vector on L1 in dBHz (from 0 to 63.75 dBHz)
%          3.6)  DF016 = code type vector on  L2: C/A=0, P_direct=1, P_cross_corr=2, P_corr=3
%          3.7)  ((DF011 + DF017) * 0.02) + (DF014 * 299792.458) = code observation vector on L2
%          3.8)  (code observation L1 + (DF018*0.0005)) / lambda2 = phase observation vector on L2
%          3.9)  DF019 = how long L2 has been locked? index vector (cycle-slip=0)
%          3.10) (DF020 * 0.25) = signal-to-noise ratio vector on L2 in dBHz (from 0 to 63.75 dBHz)
%
% DESCRIPTION:
%   RTCM format 1004 message decoding.

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

%retrieve GPS L1/L2 wavelengths
lambda(1) = goGNSS.getWavelength(goGNSS.ID_GPS, 1);
lambda(2) = goGNSS.getWavelength(goGNSS.ID_GPS, 2);

%message pointer initialization
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(6,1);
data{3} = zeros(constellations.nEnabledSat,10);

if (~constellations.GPS.enabled)
    return
end

%message number = 1004
DF002 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%reference station id
DF003 = fbin2dec(msg(pos:pos+11));  pos = pos + 12;

%TOW = time of week in milliseconds
DF004 = fbin2dec(msg(pos:pos+29));  pos = pos + 30;

%other synchronous RTCM messages flag (YES=1, NO=0)
DF005 = fbin2dec(msg(pos));  pos = pos + 1;

%number of visible satellites
DF006 = fbin2dec(msg(pos:pos+4));  pos = pos + 5;

%phase-smoothed code flag (YES=1, NO=0)
DF007 = fbin2dec(msg(pos));  pos = pos + 1;

%smoothing window
DF008 = fbin2dec(msg(pos:pos+2));  pos = pos + 3;

%output data save
data{1} = DF002;
data{2}(1) = DF003;
data{2}(2) = DF004 / 1000;
data{2}(3) = DF005;
data{2}(4) = DF006;
data{2}(5) = DF007;
data{2}(6) = DF008;

%-------------------------------------------------

%number of satellites
NSV = data{2}(4);

%data decoding for each satellite
for i = 1 : NSV

    %analyzed satellite number
    SV = fbin2dec(msg(pos:pos+5));  pos = pos + 6;

    %if GPS satellite
    if (SV >= 1 & SV <= 32)

        %L1 code type (C/A=0, P=1)
        DF010 = fbin2dec(msg(pos));  pos = pos + 1;

        %L1 pseudorange
        DF011 = fbin2dec(msg(pos:pos+23));  pos = pos + 24;

        %L1 phaserange - L1 pseudorange
        DF012 = twos_complement(msg(pos:pos+19));  pos = pos + 20;

        %L1 lock-time index (see Table 4.3-2 on RTCM manual)
        DF013 = fbin2dec(msg(pos:pos+6));  pos = pos + 7;

        %L1 pseudorange initial ambiguity
        DF014 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

        %L1-CNR (carrier-to-noise ratio): integer to be multiplied by the resolution
        DF015 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

        %---------------------------------------------------------

        %L2 code type (C/A=0, P=1,2,3)
        DF016 = fbin2dec(msg(pos:pos+1));  pos = pos + 2;

        %L2-L1 pseudorange
        DF017 = twos_complement(msg(pos:pos+13));  pos = pos + 14;

        %L2 phaserange - L1 pseudorange
        DF018 = twos_complement(msg(pos:pos+19));  pos = pos + 20;

        %L2 lock-time index (see Table 4.3-2 on RTCM specifications)
        DF019 = fbin2dec(msg(pos:pos+6));  pos = pos + 7;

        %L2-CNR (carrier-to-noise ratio): integer to be multiplied by the resolution
        DF020 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

        %---------------------------------------------------------
        
        % assign constellation-specific indexes
        idx = [];
        if (SV <= 32 && constellations.GPS.enabled)
            idx = constellations.GPS.indexes(SV);
        end

        %output data save
        data{3}(idx,1)  = DF010;
        data{3}(idx,2)  = (DF011 * 0.02) + (DF014 * 299792.458);
        data{3}(idx,3)  = (data{3}(SV,2) + (DF012*0.0005)) / lambda(1);
        data{3}(idx,4)  = DF013;
        data{3}(idx,5)  = DF015 * 0.25;
        data{3}(idx,6)  = DF016;
        data{3}(idx,7)  = (data{3}(SV,2) + (DF017 * 0.02));
        data{3}(idx,8)  = (data{3}(SV,2) + (DF018*0.0005)) / lambda(2);
        data{3}(idx,9)  = DF019;
        data{3}(idx,10) = DF020 * 0.25;

    else %SBAS satellites

        %do not store SBAS satellite information
        pos = pos + 119;

    end
end
