function rgb_img = img2rgb(img, cmap, min_max)
% SYNTAX:
%   rgb_img = img2rgb(img, minmaxm, cmap);
%
% INPUT:
%   img       image as double to be converted in rgb [n x m]
%   cmap      colormap
%   minmax    limits of the image (empty to avoid saturation)
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
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
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

    narginchk(2,3);

    if (nargin == 2)
        min_max = [min(img(:)) max(img(:))];
    end

    if diff(min_max) > 0
        id = (max(min_max(1), min(min_max(2), img)) - min_max(1)) / diff(min_max) ;
    else
        id = img * 0;
    end

    id  = floor(id * (size(cmap, 1) - 1) + 1);

    % Extract r,g,b components
    r = zeros(size(img)); r(:) = cmap(id,1);
    g = zeros(size(img)); g(:) = cmap(id,2);
    b = zeros(size(img)); b(:) = cmap(id,3);

    rgb_img = uint8(cat(3, r, g, b) * 255);
end
