function exportImg(file_name, data_map, c_map)
% Save on a file a matrix provided in "data"
%
% SYNTAX:
%   exportImg(file_name, data_map, <col_map>)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti ...
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

    narginchk(2,3);
    
    if islogical(data_map)
        imwrite(data_map, file_name);
    else
        if nargin < 3
            c_map = jet(256);
        end
        
        % Normalize data
        im_alpha = isnan(data_map);
        data_map = data_map - min(data_map(:));
        data_map = max(2, ceil(data_map ./ max(data_map(:)) .* (size(c_map, 1) -1)) + 1);
        data_map(im_alpha) = 1;
        c_map = [[0, 0, 0]; c_map];
        
        im = reshape([c_map(data_map(:), 1); c_map(data_map(:), 2); c_map(data_map(:), 3)], size(data_map, 1), size(data_map, 2), 3);
        imwrite(im, file_name, 'Alpha', 1-im_alpha);
    end
end
