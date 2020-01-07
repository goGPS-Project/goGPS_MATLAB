function [out, flag] = decode_FP32(in)

% SYNTAX:
%   [out, flag] = decode_FP32(in);
%
% INPUT:
%   in = section of the binary message transmitted by the NVS receiver
%        (string of ones and zeros)
%
% OUTPUT:
%   out  = decoded value
%   flag = validity flag for the decoded value (0: not valid; 1: valid)
%
% DESCRIPTION:
%   Single precision floating point value decoding from BINR format.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giuliano Sironi 2011
%  Contributors:     Giuliano Sironi 2011, ...
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

flag = 1;

sign = str2num(in(1));
esp  = fbin2dec(in(2:9));
mant = fbin2dec(in(10:32)) / 2^23;

if (esp > 0 && esp < 255)
    out = (-1)^sign * (2^(esp - 127)) * (1 + mant);
elseif (esp == 0 && mant ~= 0)
    out = (-1)^sign * (2^(-126)) * mant;
elseif (esp == 0 && mant == 0)
    out = (-1)^sign * 0;
elseif (esp == 255 && mant == 0)
    out = (-1)^sign * Inf;
    flag = 0;
elseif (esp == 255 && mant ~= 0)
    out = NaN;
    flag = 0;
end
