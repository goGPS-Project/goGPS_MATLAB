function [data] = decode_skytraq_MEAS_TIME(msg)

% SYNTAX:
%   [data] = decode_skytraq_MEAS_TIME(msg);
%
% INPUT:
%   msg = message transmitted by the SkyTraq receiver
%
% OUTPUT:
%   data = cell-array that contains the MEAS_TIME message information
%          1.1) message class-id (MEAS_TIME)
%          2.1) IOD  = issue of data (0-255)
%          2.2) WN = GPS week number
%          2.3) TOW  = week time (in seconds)
%          2.4) measurement period
%
% DESCRIPTION:
%   MEAS_TIME binary message decoding.

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

% first message initial index
pos = 1;

% output variable initialization
data = cell(3,1);
data{1} = 0;
data{2} = zeros(4,1);

% output data save
data{1} = 'MEAS_TIME';

% IOD decoding (1 byte)
IOD = fbin2dec(msg(pos:pos+7));   pos = pos + 8;

% GPS week decoding (2 bytes)
WN = fbin2dec(msg(pos:pos+15));   pos = pos + 16;

% time of week decoding (4 bytes)
TOW = fbin2dec(msg(pos:pos+31));  pos = pos + 32;
TOW = TOW / 1000;

% measurement period (2 bytes)
meas_period = fbin2dec(msg(pos:pos+15));

% output data save
data{2}(1) = IOD;
data{2}(2) = WN;
data{2}(3) = TOW;
data{2}(4) = meas_period;
