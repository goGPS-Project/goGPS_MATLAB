function cdf = normcdf (x, m, s)

% NORMCDF  CDF of the normal distribution
%  CDF = normcdf(X, M, S) computes the cumulative distribution
%  function (CDF) at X of the normal distribution with mean M
%  and standard deviation S.
%
%  CDF = normcdf(X) is equivalent to CDF = normcdf(X, 0, 1)

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/normcdf.m
% Original author: TT <Teresa.Twaroch@ci.tuwien.ac.at>

% Copyright (C) 1995, 1996, 1997, 2005, 2006, 2007 Kurt Hornik
% Copyright (C) 2008 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

  if (~ ((nargin == 1) || (nargin == 3)))
    error ('normcdf: you must give one or three arguments');
  end

  if (nargin == 1)
    m = 0;
    s = 1;
  end

  if (~isscalar (m) || ~isscalar (s))
    [retval, x, m, s] = common_size (x, m, s);
    if (retval > 0)
      error ('normcdf: x, m and s must be of common size or scalar');
    end
  end

  sz = size (x);
  cdf = zeros (sz);

  if (isscalar (m) && isscalar(s))
    if (find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf)))
      cdf = NaN * ones (sz);
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

end
