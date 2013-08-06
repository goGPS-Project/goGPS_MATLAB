function [out] = remove_double_10h(msg)

% SYNTAX:
%   [out] = remove_double_10h(data_in);
%
% INPUT:
%   msg = BINR binary message that may contain double 10h values
%
% OUTPUT:
%   out = BINR binary message with double 10h values compressed into one
%
% DESCRIPTION:
%   NVS BINR message data repeats 10h values twice, to distinguish them
%   from header/footer bytes. This function compresses them back to a
%   single value.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.1 beta
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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

%check the presence of <DLE><DLE>; in case remove the first
j = 1;
while (j < length(msg))
    if (msg(j) == 16 && msg(j+1) == 16)
        msg(j) = [];
    end
    j = j + 1;
end
out = msg;
