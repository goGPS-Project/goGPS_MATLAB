function [year_out] = four_digit_year(year_in)

% SYNTAX:
%   [year_out] = four_digit_year(year_in);
%
% INPUT:
%   year_in = two-digit year
%
% OUTPUT:
%   year_out = four-digit year
%
% DESCRIPTION:
%   Conversion of two-digit year values to four-digit year values
%   (following RINEX2 standard).

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

year_out(year_in > 79)  = year_in(year_in > 79)  + 1900;
year_out(year_in <= 79) = year_in(year_in <= 79) + 2000;
