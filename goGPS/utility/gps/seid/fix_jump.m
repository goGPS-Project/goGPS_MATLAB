function [fix_L2] = fix_jump(til_L2,PRN, threshold)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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


fix_L2 = til_L2(PRN,:);

notnan_L2 = fix_L2(~isnan(til_L2(PRN,:)));
fix_idx=find(abs(diff(notnan_L2)) > threshold);

L2_fixed = notnan_L2;
for i = 1 : length(fix_idx)
    L2_fixed(fix_idx(i) + 1:end) = L2_fixed(fix_idx(i) + 1:end) - (L2_fixed(fix_idx(i) + 1) - L2_fixed(fix_idx(i)));
end

fix_L2(~isnan(til_L2(PRN,:))) = L2_fixed;
