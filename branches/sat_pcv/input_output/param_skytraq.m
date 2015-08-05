function [out] = param_skytraq

% SYNTAX:
%   [out] = param_skytraq
%
% OUTPUT:
%   out = data vector
%
% DESCRIPTION:
%   Read skytraq receiver informations.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Ivan Reguzzoni
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

% inizialization
out = cell(6,2);

% String Name
Str_name = 'STQ';
%Str_name = 'SkyTraq';

% BaudRate - Serial
BaudRate = 115200;

% Buffer size - USB
Buffer_size = 16384;

% Minimun number of bytes used to synchronize data
Min_bytes = 0;

% String name in decoding phase (NA = Not Available)
TimeRaw_id = 'NA';        % message with both timing and observations
Eph_id     = 'GPS_EPH';   % message with ephemeris
Hui_id     = 'NA';        % message with ionosphere parameters
Time_id    = 'MEAS_TIME'; % message with timing
Raw_id     = 'RAW_MEAS';  % message with observations
Track_id   = 'NA';        % message with tracking data

out(1,1) = {Str_name};
out(2,1) = {BaudRate};
out(3,1) = {Buffer_size};
out(4,1) = {Min_bytes};

out(1,2) = {TimeRaw_id};
out(2,2) = {Eph_id};
out(3,2) = {Hui_id};
out(4,2) = {Time_id};
out(5,2) = {Raw_id};
out(6,2) = {Track_id};
