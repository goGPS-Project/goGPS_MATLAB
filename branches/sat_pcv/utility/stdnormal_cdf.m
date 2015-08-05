function cdf = stdnormal_cdf (x)
% For each component of x, compute the CDF of the standard normal
% distribution at x.

% Description: CDF of the standard normal distribution

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 1995-2011 Kurt Hornik
% Copyright (C) 2012-2014 Mirko Reguzzoni, Eugenio Realini
%
% Adapted from Octave.
% Author: KH <Kurt.Hornik@wu-wien.ac.at>
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

if (nargin ~= 1)
    error('Requires one input arguments.');
end

if (numel(x) == 0)
    error ('stdnormal_cdf: X must not be empty');
end

cdf = erfc (x / (-sqrt(2))) / 2;
