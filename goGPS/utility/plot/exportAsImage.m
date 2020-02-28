function rgb_img = exportAsImage(file_name, data, cmap, min_max, x0, y0, x1, y1)
% SYNTAX:
%   rgb_img = exportAsImage(file_name, data, cmap, min_max, x0, y0, x1, y1)
%   rgb_img = exportAsImage(file_name, data, cmap, min_max, x1, y1)
%
% INPUT:
%   file_name
%   data      [n x m] image
%   min_max   limits of the image (empty to avoid saturation)
%   x0, y0    meshgrid of original coordinate ( or size of the original image)
%   x1, y1    meshgrid of final exported coordinate ( or size of the final exported image)
%
% OUTPUT:
%   rgb_img = [n x m x 3] RGB image
%
% DESCRIPTION:
%   Convert an grayscale image to a rgb image

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
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

    narginchk(3,8);

    if (nargin == 3) || isempty(min_max)
        min_max = [min(data(:)) max(data(:))];
    end

    if (nargin == 6)
        x1 = x0;
        y1 = y0;
        x0 = size(data,1);
        y0 = size(data,2);
    end
    if (nargin == 6) || (nargin == 8)
        if (numel(x0) == 1)
            x0 = 1 : x0;
            y0 = 1 : y0;
            x1 = 1 : x1;
            y1 = 1 : y1;
            [x0, y0] = meshgrid(x0, y0);
            [x1, y1] = meshgrid(x1, y1);
        end
        tic; data = interp2(x0, y0, data, x1, y1, 'spline'); toc;
    end
    clear x0 x1 y0 y1;
    rgb_img = img2rgb(data, cmap, min_max);
    imwrite(rgb_img, file_name);
end
