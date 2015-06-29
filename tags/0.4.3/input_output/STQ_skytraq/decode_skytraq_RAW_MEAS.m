function [data] = decode_skytraq_RAW_MEAS(msg, constellations)

% SYNTAX:
%   [data] = decode_skytraq_RAW_MEAS(msg, constellations);
%
% INPUT:
%   msg = message transmitted by the SkyTraq receiver
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   data = cell-array that contains the RAW_MEAS packet information
%          1.1) message class-id (RAW_MEAS)
%          2.1) IOD   = issue of data (0-255)
%          2.2) NMEAS = number of measurements
%          3.1) PRN   = satellite PRN
%          3.2) CN0   = signal-to-noise ratio (in dbHz)
%          3.3) PRM   = pseudorange measurements (C/A code in meters)
%          3.4) CPM   = carrier phase measurements (in cycles)
%          3.5) DOM   = doppler measurements (in Hertz)
%          3.6) channel indicator
%
% DESCRIPTION:
%   RAW_MEAS binary message decoding.

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
data{2} = zeros(2,1);
data{3} = zeros(constellations.nEnabledSat,6);

%output data save
data{1} = 'RAW_MEAS';

% IOD decoding (1 byte)
IOD = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

% number of measurements decoding (1 byte)
NMEAS = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

%output data save
data{2}(1) = IOD;
data{2}(2) = NMEAS;

%read the measurements of every satellite
for j = 1 : NMEAS
    
    % satellite PRN decoding
    PRN = fbin2dec(msg(pos:pos+7));
    pos = pos + 8;

    % signal-to-noise ratio decoding (in dBHz)
    CN0 = fbin2dec(msg(pos:pos+7));
    pos = pos + 8;

    %------------------------------------------------

    % C/A pseudorange measurement decoding (in meters)
    C1field = msg(pos:pos+63);
    pos = pos + 64;
    
    % floating point value decoding (double floating point)
    sign = fbin2dec(C1field(1));
    esp  = fbin2dec(C1field(2:12));
    mant = fbin2dec(C1field(13:64)) / 2^52;
    C1 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);
    
    %------------------------------------------------

    % L1 phase measurement decoding (in cycles)
    L1field = msg(pos:pos+63);
    pos = pos + 64;

    % floating point value decoding (double floating point)
    sign = fbin2dec(L1field(1));
    esp  = fbin2dec(L1field(2:12));
    mant = fbin2dec(L1field(13:64)) / 2^52;
    L1 = (-1)^sign * (2^(esp - 1023)) * (1 + mant);

    %------------------------------------------------

    % doppler measurements decoding (in Hz)
    D1field = msg(pos:pos+31);
    pos = pos + 32;

    % floating point value decoding (single floating point)
    sign = fbin2dec(D1field(1));
    esp  = fbin2dec(D1field(2:9));
    mant = fbin2dec(D1field(10:32)) / 2^23;
    D1 = (-1)^sign * (2^(esp - 127)) * (1 + mant);

    %------------------------------------------------

    % channel indicator
    channel_indicator = fbin2dec(msg(pos:pos+7));
    pos = pos + 8;

    % assign constellation-specific indexes
    idx = [];
    if (PRN <= 32)
        idx = constellations.GPS.indexes(PRN);
    end
    
    % phase, code and doppler measure save
    CPM = L1;
    PRM = C1;
    DOM = D1;
    
    % data output save
    data{3}(idx,1) = PRN;
    data{3}(idx,2) = CN0;
    data{3}(idx,3) = PRM;
    data{3}(idx,4) = CPM;
    data{3}(idx,5) = DOM;
    data{3}(idx,6) = channel_indicator;
end
