function inv = tinv (x, n)
% For each probability value x, compute the inverse of the
% cumulative distribution function (CDF) of the t (Student)
% distribution with degrees of freedom n.  This function is
% analogous to looking in a table for the t-value of a single-tailed
% distribution.

% For very large n, the "correct" formula does not really work well,
% and the quantiles of the standard normal distribution are used
% directly.

% Description: Quantile function of the t distribution

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

if (nargin ~= 2)
    error('Requires two input arguments.');
end

if (~isscalar (n))
    [retval, x, n] = common_size (x, n);
    if (retval > 0)
        error ('tinv: X and N must be of common size or scalar');
    end
end

inv = zeros (size (x));

k = find ((x < 0) | (x > 1) | isnan (x) | ~(n > 0));
if (any (k))
    inv(k) = NaN;
end

k = find ((x == 0) & (n > 0));
if (any (k))
    inv(k) = -Inf;
end

k = find ((x == 1) & (n > 0));
if (any (k))
    inv(k) = Inf;
end

k = find ((x > 0) & (x < 1) & (n > 0) & (n < 10000));
if (any (k))
    if (isscalar (n))
        inv(k) = (sign (x(k) - 1/2) .* sqrt (n .* (1 ./ betainv (2*min (x(k), 1 - x(k)), n/2, 1/2) - 1)));
    else
        inv(k) = (sign (x(k) - 1/2) .* sqrt (n(k) .* (1 ./ betainv (2*min (x(k), 1 - x(k)), n(k)/2, 1/2) - 1)));
    end
end

% For large n, use the quantiles of the standard normal
k = find ((x > 0) & (x < 1) & (n >= 10000));
if (any (k))
    inv(k) = stdnormal_inv (x(k));
end
