% SYNTAX:
%    flag = flagExpand(flag, expand_size)
%
% DESCRIPTION:
%    expand the flag array
%    if flag is a matrix this flagging expansion will work column by column
%
% INPUT:
%   flag          [n_obs x n_arrays]
%   expand_size   n_epochs with flags to activate at the border of a flagged interval

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti
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

function flag = flagExpand(flag, expand_size)
    % compute a moving window median to filter the data in input
    
    for c = 1 : size(flag, 2)
        flag(:, c) = conv(single(flag(:, c)), ones(2 * expand_size + 1, 1)', 'same') > 0;
    end
end
