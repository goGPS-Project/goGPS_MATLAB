function [flag_intervals] = getOutliers(flags, split_point)
% INPUT
%   flags           flag array (as logical)
%   split_point     points with forced splits (as logical)
%
% SYNTAX
%   [flag_intervals] = getOutliers(flags, <split_point>)
%
% DESCRIPTION
%   Returns start and end of flagged intervals
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 2
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
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

    % convert flags into array
    if isstruct(flags)
        flag_array = int8(struct2flagVec(flags, max(flags.pos + 1)));
        % add padding to avoid problem with flags on the borders
    else
        flag_array = flags;
    end
    flag_array = flag_array(:);
    flag_array = [0; flag_array; 0];
    flag_array(flag_array ~= 0) = 1;
    diff_tmp = flag_array(1:end-1) - flag_array(2 : end);
    clear flagArray;
    if nargin == 2 && ~isempty(split_point)
        split_point = split_point & flags;
        d_tmp = diff_tmp; 
        d_tmp(1 : end-1) = (d_tmp(1 : end-1) - 2 * split_point);
        f_i1 = find(d_tmp < 0);
        f_i2 = find([split_point; false] | d_tmp > 0);
        for i = 1 : numel(f_i1)
            if f_i1(i) == f_i2(i)
                f_i2(i) = [];
            end
        end
        flag_intervals = [f_i1, f_i2 - 1];
    else
        flag_intervals = [find(diff_tmp < 0), find(diff_tmp > 0) - 1];
    end
end

function [flag_array] = struct2flagVec(flags, maxSize)
    flag_array = false(maxSize,1);
    flag_array(flags.pos) = flags.val;
end
