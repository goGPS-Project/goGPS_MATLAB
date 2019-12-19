function inv = stdnormal_inv (x)
% For each component of x, compute the quantile (the
% inverse of the CDF) at x of the standard normal distribution.

% Description: Quantile function of the standard normal distribution

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
% Copyright (C) 1995-2011 Kurt Hornik
% Copyright (C) 2012-2014 Mirko Reguzzoni, Eugenio Realini
%  Written by:       KH <Kurt.Hornik@wu-wien.ac.at>
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

if (nargin ~= 1)
    error('Requires one input arguments.');
end

inv = sqrt (2) * erfinv (2 * x - 1);
