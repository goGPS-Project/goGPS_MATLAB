function [time] = weektow2time(week, sow, sys)

% SYNTAX:
%   [time] = weektow2time(week, sow, sys);
%
% INPUT:
%   week = GNSS week
%   sow  = GNSS seconds-of-week
%   sys  = GNSS system identifier
%
% OUTPUT:
%   time = GPS time (continuous since 6-1-1980)
%
% DESCRIPTION:
%   Conversion from GNSS time in week, seconds-of-week format to GPS time
%   in continuous format (similar to datenum).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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

if (~strcmp(sys, 'C'))
    time = week*7*86400 + sow;
else
    GPS_BDS_week = 1356; %GPS week number on 1st January 2006 (start of BeiDou time)
    time = (GPS_BDS_week + week)*7*86400 + sow;
end
