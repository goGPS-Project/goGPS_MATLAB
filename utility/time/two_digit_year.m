function [year_out] = two_digit_year(year_in)

% SYNTAX:
%   [year_out] = two_digit_year(year_in);
%
% INPUT:
%   year_in = four-digit year
%
% OUTPUT:
%   year_out = two-digit year
%
% DESCRIPTION:
%   Conversion of four-digit year values to two-digit year values.

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

year_out = zeros(size(year_in));

year_out(year_in >= 2000) = year_in(year_in >= 2000) - 2000;
year_out(year_in <  2000) = year_in(year_in <  2000) - 1900;
