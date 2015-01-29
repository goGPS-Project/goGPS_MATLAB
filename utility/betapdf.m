function pdf = betapdf (x, a, b)
% For each element of x, returns the PDF at x of the beta
% distribution with parameters a and b.

% Description: PDF of the Beta distribution

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 1995-2011 Kurt Hornik
% Copyright (C) 2010 Christos Dimitrakakis
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Adapted from Octave.
% Author: KH <Kurt.Hornik@wu-wien.ac.at>, CD <christos.dimitrakakis@gmail.com>
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

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('betapdf: X, A and B must be of common size or scalar');
    end
end

sz = size (x);
pdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    pdf (k) = NaN;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0) & ((a ~= 1) | (b ~= 1)));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = exp ((a - 1) .* log (x(k)) + (b - 1) .* log (1 - x(k)) + gammaln(a + b) - gammaln(a) - gammaln(b));
    else
        pdf(k) = exp ((a(k) - 1) .* log (x(k)) + (b(k) - 1) .* log (1 - x(k)) + gammaln(a(k) + b(k)) - gammaln(a(k)) - gammaln(b(k)));
    end
end

% Most important special cases when the density is finite.
k = find ((x == 0) & (a == 1) & (b > 0) & (b ~= 1));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = exp(gammaln(a + b) - gammaln(a) - gammaln(b));
    else
        pdf(k) = exp(gammaln(a(k) + b(k)) - gammaln(a(k)) - gammaln(b(k)));
    end
end

k = find ((x == 1) & (b == 1) & (a > 0) & (a ~= 1));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = exp(gammaln(a + b) - gammaln(a) - gammaln(b));
    else
        pdf(k) = exp(gammaln(a(k) + b(k)) - gammaln(a(k)) - gammaln(b(k)));
    end
end

k = find ((x >= 0) & (x <= 1) & (a == 1) & (b == 1));
if (any (k))
    pdf(k) = 1;
end

% Other special case when the density at the boundary is infinite.
k = find ((x == 0) & (a < 1));
if (any (k))
    pdf(k) = Inf;
end

k = find ((x == 1) & (b < 1));
if (any (k))
    pdf(k) = Inf;
end
