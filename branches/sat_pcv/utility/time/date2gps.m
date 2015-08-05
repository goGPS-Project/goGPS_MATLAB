function [gps_week, gps_sow, gps_dow] = date2gps(date)

% SYNTAX:
%   [gps_week, gps_sow, gps_dow] = date2gps(date);
%
% INPUT:
%   date = date [year, month, day, hour, min, sec]
%
% OUTPUT:
%   gps_week = GPS week
%   gps_sow  = GPS seconds of week
%   gps_dow  = GPS day of week
%
% DESCRIPTION:
%   Conversion from calendar date to GPS time.

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

gps_start_datenum = 723186; %This is datenum([1980,1,6,0,0,0])

%number of days since the beginning of GPS time
deltat   = (datenum([date(:,1), date(:,2), date(:,3)]) - gps_start_datenum);

gps_week = floor(deltat/7);            %GPS week
gps_dow  = floor(deltat - gps_week*7); %GPS day of week
gps_sow  = (deltat - gps_week*7)*86400; 
gps_sow = gps_sow + date(:,4)*3600 + date(:,5)*60 + date(:,6); %GPS seconds of week

% %alternative way, using the Julian day
% jd = date2jd(date);
% [gps_week, gps_dow, gps_sow] = jd2gps(jd);
