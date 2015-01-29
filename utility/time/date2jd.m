function [jd, mjd] = date2jd(date)

% SYNTAX:
%   [jd, mjd] = date2jd(date);
%
% INPUT:
%   date = date [year, month, day, hour, min, sec]
%
% OUTPUT:
%   jd  = julian day
%   mjd = modified julian day
%
% DESCRIPTION:
%   Conversion from date to julian day and modified julian day.

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

year  = date(:,1);
month = date(:,2);
day   = date(:,3);

pos = find(month <= 2);
year(pos)  = year(pos) - 1;
month(pos) = month(pos) + 12;

%julian day
jd = floor(365.25*(year+4716)) + floor(30.6001*(month+1)) + day - 1537.5;

%modified julian day
mjd = jd - 2400000.5;
