function [time] = weektow2time(week, sow)

% SYNTAX:
%   [time] = weektow2time(week, sow);
%
% INPUT:
%   week = GPS week
%   sow  = GPS seconds-of-week
%
% OUTPUT:
%   time = GPS time (continuous since 6-1-1980)
%
% DESCRIPTION:
%   Conversion from GPS time in week, seconds-of-week format to GPS time in
%   continuous format (similar to datenum).

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

time = week*7*86400 + sow;
