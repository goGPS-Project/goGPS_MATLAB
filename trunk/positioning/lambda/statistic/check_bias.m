function bias = check_bias(D, b)

% SYNTAX:
%   bias = check_bias(D, b);
%
% INPUT:
%   D = diagonal of Qa matrix
%   b = bias matrix, if available (optional)
%
% OUTPUT:
%   bias = bias vector
%
% DESCRIPTION:
%   Check bias availability. If bias is available the value is returned, 
%   otherwise it is a zero-vector. Minimum input one-variable, maximum input
%   two-variables.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Hendy F. Suhandri
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

%'narginchk' function introduced since MATLAB 7.13 (R2011b)
%
% minargs=1; maxargs=2;
% narginchk(minargs, maxargs)
% 
% nVararg = nargin;
% 
% r = length(varargin{1});
% if nVararg == maxargs
%     bias = varargin{2};
%     
% elseif nVararg == minargs
%     bias = zeros(r,1);
% end

if (nargin > 1)
    bias = b;
else
    bias = zeros(length(D));
end

