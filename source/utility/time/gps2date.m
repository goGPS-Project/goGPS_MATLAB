function [date, doy, dow] = gps2date(gps_week, gps_sow)

% SYNTAX:
%   [date, doy, dow] = gps2date(gps_week, gps_sow);
%
% INPUT:
%   gps_week = GPS week
%   gps_sow  = GPS seconds of week
%
% OUTPUT:
%   date = date [year month day hour min sec]
%   doy  = day of year
%   dow  = day of week
%
% DESCRIPTION:
%   Conversion from GPS time to calendar date and day of year (DOY).

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
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

gps_start_datenum = 723186; %This is datenum([1980,1,6,0,0,0])

gps_dow = fix(gps_sow/86400);                             %day of week
date = datevec(gps_start_datenum + 7*gps_week + gps_dow); %calendar date up to days
gps_sod = gps_sow - gps_dow*86400;                        %seconds of day
date(:,4) = floor(gps_sod/3600);                          %hours
date(:,5) = floor(gps_sod/60 - date(:,4)*60);             %minutes
date(:,6) = gps_sod - date(:,4)*3600 - date(:,5)*60;      %seconds

%day of year (DOY)
if (nargout > 1)
    doy = date2doy(datenum(date));
    doy = floor(doy);
end

%day of week
if (nargout > 2)
    dow = gps_dow;
end
