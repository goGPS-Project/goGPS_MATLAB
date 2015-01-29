function [data] = decode_F5h(msg, constellations)

% SYNTAX:
%   [data] = decode_F5h(msg, constellations)
%
% INPUT:
%   msg = message transmitted by the NVS receiver
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the RXM-RAW packet information
%          1.1) message id (F5h)
%          2.1) TOW_GPS = GPS time-of-week (in seconds)
%          2.2) WEEK    = GPS week
%          2.3) NSV     = number of visible satellites
%          3.1) CPM     = phase measurements (in cycles)
%          3.2) PRM     = pseudorange measurements (C/A code in meters)
%          3.3) DOM     = doppler measurements (in Hertz)
%          3.4) SV      = space vehicle number
%          3.5) RDF     = raw data flags
%          3.6) CNO     = signal-to-noise ratio (in dbHz)
%
% DESCRIPTION:
%   BINR F5h binary message decoding.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Daisuke Yoshida
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

%first message initial index
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(3,1);
data{3} = zeros(constellations.nEnabledSat,6);

%output data save
data{1} = 'F5h';

%------------------------------------------------

%check the minimum allowed length of the message
if (length(msg) < 216)
    return
end

%calculate the number of visible satellites
NSV = (length(msg) - 216)/240; % 27*8 bits = 216, 30*8 bits = 240

%------------------------------------------------

%time of measurement (UTC) in ms (8 byte)
TOW_UTC = msg(pos:pos+63); pos = pos + 64;
TOW_UTC = fliplr(reshape(TOW_UTC,8,[]));                  % byte order inversion (little endian)
TOW_UTC = TOW_UTC(:)';

%floating point value decoding (double floating point)
sign = str2num(TOW_UTC(1));
esp  = fbin2dec(TOW_UTC(2:12));
mant = fbin2dec(TOW_UTC(13:64)) / 2^52;
TOW_UTC = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

%------------------------------------------------

%GPS week decoding (2 byte)
WEEK1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
WEEK2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
WEEK = WEEK1 + (WEEK2 * 2^8);        % little endian

%------------------------------------------------

%GPS-UTC time shift in ms (8 byte)
GPS_time_shift = msg(pos:pos+63); pos = pos + 64;
GPS_time_shift = fliplr(reshape(GPS_time_shift,8,[]));                  % byte order inversion (little endian)
GPS_time_shift = GPS_time_shift(:)';

%floating point value decoding (double floating point)
sign = str2num(GPS_time_shift(1));
esp  = fbin2dec(GPS_time_shift(2:12));
mant = fbin2dec(GPS_time_shift(13:64)) / 2^52;
GPS_time_shift = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

%------------------------------------------------

%GLONASS-UTC time shift in ms (8 byte)
GLONASS_time_shift = msg(pos:pos+63); pos = pos + 64;
GLONASS_time_shift = fliplr(reshape(GLONASS_time_shift,8,[]));          % byte order inversion (little endian)
GLONASS_time_shift = GLONASS_time_shift(:)';

%floating point value decoding (double floating point)
sign = str2num(GLONASS_time_shift(1));
esp  = fbin2dec(GLONASS_time_shift(2:12));
mant = fbin2dec(GLONASS_time_shift(13:64)) / 2^52;
GLONASS_time_shift = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

%------------------------------------------------

%receiver time scale correction in ms (1 byte)
Time_Correction = twos_complement(msg(pos:pos+7)); pos = pos + 8;

%------------------------------------------------

%compute GPS and GLONASS TOWs (in seconds)
TOW_GPS = (TOW_UTC + GPS_time_shift) / 1e3;
TOW_GLO = (TOW_UTC + GLONASS_time_shift) / 1e3;

%output data save
data{2}(1) = TOW_GPS;
data{2}(2) = WEEK;
data{2}(3) = NSV;

nGLO = 0;
nGPS = 0;
nQZS = 0;
nGAL = 0;

%read the measurements of every satellite
for j = 1 : NSV
    
    %signal type 01=GLONASS, 02=GPS/QZSS, 04=SBAS, 08=Galileo
    signal_type = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    %------------------------------------------------
    
    %satellite PRN
    SID = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    %------------------------------------------------
    
    %carrier number for GLONASS
    num4GLONASS = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    %------------------------------------------------
    
    %signal-to-noise ratio (in dBHz)
    CNO = fbin2dec(msg(pos:pos+7)); pos = pos + 8;
    
    %------------------------------------------------
    
    %L1 phase measurement (in cycles)
    L1field = msg(pos:pos+63);
    pos = pos + 64;
    
    %byte order inversion (little endian)
    L1field = fliplr(reshape(L1field,8,[]));
    L1field = L1field(:)';
    
    %floating point value decoding (double floating point)
    sign = fbin2dec(L1field(1));
    esp  = fbin2dec(L1field(2:12));
    mant = fbin2dec(L1field(13:64)) / 2^52;
    L1 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);%
    
    %------------------------------------------------
    
    %C/A pseudorange measurement (in ms)
    C1field = msg(pos:pos+63);
    pos = pos + 64;
    
    %byte order inversion (little endian)
    C1field = fliplr(reshape(C1field,8,[]));
    C1field = C1field(:)';
    
    %floating point value decoding (double floating point)
    sign = fbin2dec(C1field(1));
    esp  = fbin2dec(C1field(2:12));
    mant = fbin2dec(C1field(13:64)) / 2^52;
    C1 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);
    
    %------------------------------------------------
    
    %doppler measurements decoding (in dB-Hz)
    D1field = msg(pos:pos+63);
    pos = pos + 64;
    
    %byte order inversion (little endian)
    D1field = fliplr(reshape(D1field,8,[]));
    D1field = D1field(:)';
    
    %floating point value decoding (single floating point)
    sign = fbin2dec(D1field(1));
    esp  = fbin2dec(D1field(2:12));
    mant = fbin2dec(D1field(13:64)) / 2^52;
    D1 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);
    if (abs(D1) > 1e5), D1 = 0; end
    
    %------------------------------------------------
    
    %raw data flags
    RDF = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
    
    %------------------------------------------------
    
    %reserved
    pos = pos + 8;
    
    %assign constellation-specific indexes
    idx = [];
    if (SID && signal_type == 1 && constellations.GLONASS.enabled && SID <= constellations.GLONASS.numSat)
        
        idx = constellations.GLONASS.indexes(SID);
        nGLO = nGLO + 1;
        
    elseif (SID && signal_type == 2 && constellations.GPS.enabled && SID <= constellations.GPS.numSat)
        
        idx = constellations.GPS.indexes(SID);
        nGPS = nGPS + 1;
        
    elseif (SID && signal_type == 2 && constellations.QZSS.enabled && SID == 33)
        
        SID = SID-32;
        idx = constellations.QZSS.indexes(SID);
        nQZS = nQZS + 1;
        
    elseif (SID && signal_type == 8 && constellations.Galileo.enabled && SID <= constellations.Galileo.numSat)
        
        idx = constellations.Galileo.indexes(SID);
        nGAL = nGAL + 1;
    end
    
    %phase, code and doppler measure save
    CPM = L1;
    PRM = C1*goGNSS.V_LIGHT*1e-3; %in meters
    DOM = D1;
    
    if (~isempty(idx))
        %data output save
        data{3}(idx,1) = CPM;
        data{3}(idx,2) = PRM;
        data{3}(idx,3) = DOM;
        data{3}(idx,4) = SID;
        data{3}(idx,5) = RDF;
        data{3}(idx,6) = CNO;
    end
end

if (any(data{3}(:,2) < 0) || any(data{3}(:,2) > 60e6) || NSV <= 0 || nGPS+nGLO+nGAL+nQZS == 0)

    %discard everything if there is at least one anomalous pseudorange
    % (i.e. if pseudorange ambiguity not solved yet)
    data = cell(3,1);
    data{1} = 0;
    data{2} = zeros(3,1);
    data{3} = zeros(constellations.nEnabledSat,6);

%     %discard satellites with anomalous pseudoranges
%     pos1 = find(data{3}(:,2) < 0);
%     pos2 = find(data{3}(:,2) > 60e6);
%     pos  = union(pos1,pos2);
%     data{3}(pos,:) = zeros(length(pos),6);
end
