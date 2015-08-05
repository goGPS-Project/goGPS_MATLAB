function inv = finv (x, m, n)
% For each component of x, compute the quantile (the inverse of
% the CDF) at x of the F distribution with parameters m and
% n.

% Description: Quantile function of the F distribution

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

if (nargin ~= 3)
    error('Requires three input arguments.');
end

if (~isscalar (m) || ~isscalar (n))
    [retval, x, m, n] = common_size (x, m, n);
    if (retval > 0)
        error ('finv: X, M and N must be of common size or scalar');
    end
end

sz = size (x);
inv = zeros (sz);

k = find ((x < 0) | (x > 1) | isnan (x) | ~(m > 0) | ~(n > 0));
if (any (k))
    inv(k) = NaN;
end

k = find ((x == 1) & (m > 0) & (n > 0));
if (any (k))
    inv(k) = Inf;
end

k = find ((x > 0) & (x < 1) & (m > 0) & (n > 0));
if (any (k))
    if (isscalar (m) && isscalar (n))
        inv(k) = ((1 ./ betainv (1 - x(k), n / 2, m / 2) - 1) .* n ./ m);
    else
        inv(k) = ((1 ./ betainv (1 - x(k), n(k) / 2, m(k) / 2) - 1) .* n(k) ./ m(k));
    end
end
