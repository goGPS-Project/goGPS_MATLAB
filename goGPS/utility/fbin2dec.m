function x=fbin2dec(s)
%FBIN2DEC (fast bin2dec) Convert binary string to decimal integer.
%   X = FBIN2DEC(B) interprets the binary string B and returns in X the
%   equivalent decimal number. It is a stripped version of "bin2dec", with
%   a minimal check on input.
%
%   If B is a character array, or a cell array of strings, each row is
%   interpreted as a binary string.
%
%   Example
%       fbin2dec('010111') returns 23

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
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

% handle input
s = char(s);

[m,n] = size(s);

% Convert to numbers
v = s - '0';
twos = pow2(n-1:-1:0);
x = sum(v .* twos(ones(m,1),:),2);
