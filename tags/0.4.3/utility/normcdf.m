function cdf = normcdf (x, m, s)
% For each element of x, compute the cumulative distribution
% function (CDF) at x of the normal distribution with mean
% m and standard deviation s.
% Default values are m = 0, s = 1.

% Description: CDF of the normal distribution

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 1995-2011 Kurt Hornik
% Copyright (C) 2012-2014 Mirko Reguzzoni, Eugenio Realini
%
% Adapted from Octave.
% Author: TT <Teresa.Twaroch@ci.tuwien.ac.at>
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

if (~ ((nargin == 1) || (nargin == 3)))
    error('Requires one or three input arguments.');
end

if (nargin == 1)
    m = 0;
    s = 1;
end

if (~isscalar (m) || ~isscalar (s))
    [retval, x, m, s] = common_size (x, m, s);
    if (retval > 0)
        error ('normcdf: X, M and S must be of common size or scalar');
    end
end

sz = size (x);
cdf = zeros (sz);

if (isscalar (m) && isscalar(s))
    if (find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf)))
        cdf = NaN (sz);
    else
        cdf =  stdnormal_cdf ((x - m) ./ s);
    end
else
    k = find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf));
    if (any (k))
        cdf(k) = NaN;
    end
    
    k = find (~isinf (m) & ~isnan (m) & (s >= 0) & (s < Inf));
    if (any (k))
        cdf(k) = stdnormal_cdf ((x(k) - m(k)) ./ s(k));
    end
end

cdf((s == 0) & (x == m)) = 0.5;
