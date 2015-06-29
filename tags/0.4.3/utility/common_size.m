function [errorcode, varargout] = common_size (varargin)

% COMMON_SIZE  Checks that all inputs are either scalar or of common size
%  [ERR, Y1, ...] = common_size(X1, ...)
%  Determine if all input arguments are either scalar or of common
%  size.  If so, ERR is zero, and YI is a matrix of the
%  common size with all entries equal to XI if this is a scalar or
%  XI otherwise.  If the inputs cannot be brought to a common size,
%  errorcode is 1, and YI is XI.
%
%  Example:
%   [errorcode, a, b] = common_size([1 2; 3 4], 5)
%      >> errorcode = 0
%      >> a = [ 1, 2; 3, 4 ]
%      >> b = [ 5, 5; 5, 5 ]

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 1995-2007 Kurt Hornik
% Copyright (C) 2008 Dynare Team
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Adapted from Octave, Dynare.
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

if (nargin < 2)
    error ('common_size: only makes sense if nargin >= 2');
end

len = 2;
for i = 1 : nargin
    sz =  size (varargin{i});
    if (length (sz) < len)
        s(i,:) = [sz, ones(1,len - length(sz))];
    else
        if (length (sz) > len)
            if (i > 1)
                s = [s, ones(size(s,1), length(sz) - len)];
            end
            len = length (sz);
        end
        s(i,:) = sz;
    end
end

m = max (s);
if (any (any ((s ~= 1)') & any ((s ~= ones (nargin, 1) * m)')))
    errorcode = 1;
    varargout = varargin;
else
    errorcode = 0;
    for i = 1 : nargin
        varargout{i} = varargin{i};
        if (prod (s(i,:)) == 1)
            varargout{i} = varargout{i} * ones (m);
        end
    end
end
