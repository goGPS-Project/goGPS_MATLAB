function [data] = decode_RXM_RAW(msg, constellations)

% SYNTAX:
%   [data] = decode_RXM_RAW(msg, constellations)
%
% INPUT:
%   msg = message transmitted by the u-blox receiver
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the RXM-RAW packet information
%          1.1) message class-id (RXM-RAW)
%          2.1) TOW  = week time (in seconds)
%          2.2) WEEK = GPS week
%          2.3) NSV  = number of visible satellites
%          2.4) RES  = reserved field (not used)
%          3.1) CPM  = phase measurements (in cycles)
%          3.2) PRM  = pseudorange measurements (C/A code in meters)
%          3.3) DOM  = doppler measurements (in Hertz)
%          3.4) SV   = space vehicle number
%          3.5) MQI  = measurement quality index
%          3.6) CNO  = signal-to-noise ratio (in dbHz)
%          3.7) LLI  = loss of lock indicator
%
% DESCRIPTION:
%   RXM-RAW binary message decoding.

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

% first message initial index
pos = 1;

%output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(4,1);
data{3} = zeros(constellations.nEnabledSat,7);

%output data save
data{1} = 'RXM-RAW';

% week time decoding (4 byte)
TOW1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW3 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW4 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
TOW = TOW1 + (TOW2 * 2^8) + (TOW3 * 2^16) + (TOW4 * 2^24);  % little endian
TOW = TOW / 1000;

% GPS week decoding (2 byte)
WEEK1 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
WEEK2 = fbin2dec(msg(pos:pos+7));  pos = pos + 8;
WEEK = WEEK1 + (WEEK2 * 2^8);        % little endian

% number of visible satellites (1 byte)
NSV = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

% reserved field (1 byte)
RES = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

%output data save
data{2}(1) = TOW;
data{2}(2) = WEEK;
data{2}(3) = NSV;
data{2}(4) = RES;

%read the measurements of every satellite
for j = 1 : NSV

    % L1 phase measurement decoding (in cycles)
    L1field = msg(pos:pos+63);
    pos = pos + 64;

    % byte order inversion (little endian)
    L1field = fliplr(reshape(L1field,8,[]));
    L1field = L1field(:)';

    % floating point value decoding (double floating point)
    sign = fbin2dec(L1field(1));
    esp  = fbin2dec(L1field(2:12));
    mant = fbin2dec(L1field(13:64)) / 2^52;
    L1 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------

    % C/A pseudorange measurement decoding (in meters)
    C1field = msg(pos:pos+63);
    pos = pos + 64;

    % byte order inversion (little endian)
    C1field = fliplr(reshape(C1field,8,[]));
    C1field = C1field(:)';

    % floating point value decoding (double floating point)
    sign = fbin2dec(C1field(1));
    esp  = fbin2dec(C1field(2:12));
    mant = fbin2dec(C1field(13:64)) / 2^52;
    C1 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------

    % doppler measurements decoding (in Hz)
    D1field = msg(pos:pos+31);
    pos = pos + 32;

    % byte order inversion (little endian)
    D1field = fliplr(reshape(D1field,8,[]));
    D1field = D1field(:)';

    % floating point value decoding (single floating point)
    sign = fbin2dec(D1field(1));
    esp  = fbin2dec(D1field(2:9));
    mant = fbin2dec(D1field(10:32)) / 2^23;
    D1 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

    %------------------------------------------------
    
    % satellite number decoding
    SV = fbin2dec(msg(pos:pos+7));
    pos = pos + 8;
    
    % quality index decoding
    MQI = fbin2dec(msg(pos:pos+7));
    pos = pos + 8;
    
    % signal-to-noise ratio decoding (in dBHz)
    CNO = fbin2dec(msg(pos:pos+7));
    pos = pos + 8;
    
    % signal loss index decoding
    LLI = fbin2dec(msg(pos:pos+7));
    pos = pos + 8;

    % assign constellation-specific indexes
    idx = [];
    if (SV <= 32 && constellations.GPS.enabled)
        idx = constellations.GPS.indexes(SV);
    end
        
    % phase, code and doppler measure save
    CPM = L1;
    PRM = C1;
    DOM = D1;
    
    % data output save
    data{3}(idx,1) = CPM;
    data{3}(idx,2) = PRM;
    data{3}(idx,3) = DOM;
    data{3}(idx,4) = SV;
    data{3}(idx,5) = MQI;
    data{3}(idx,6) = CNO;
    data{3}(idx,7) = LLI;
end
