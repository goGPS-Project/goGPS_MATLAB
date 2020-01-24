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

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------


time = week*7*86400 + sow;

BDS_mask = (sys == 'C');
if (any(BDS_mask))
    GPS_BDS_week = 1356; %GPS week number on 1st January 2006 (start of BeiDou time)
    time(BDS_mask) = (GPS_BDS_week + week(BDS_mask))*7*86400 + sow(BDS_mask);
end
