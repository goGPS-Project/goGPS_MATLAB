% SYNTAX
%    [grid_s] = circConv2(grid, conv_x_pixels, conv_y_pixels)
%
% DESCRIPTION
%    Perform a convolution of x by y pixel, replicating the map to limit border effects
%
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
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

function [grid_s] = circConv2(grid, conv_x_pixels, conv_y_pixels)

if (nargin == 2)
    if numel(conv_x_pixels) > 2
        kernel = conv_x_pixels;
        conv_x_pixels = size(kernel,2);
        conv_y_pixels = size(kernel,1);
    elseif numel(conv_x_pixels) == 2
        conv_y_pixels = conv_x_pixels(2);
        conv_x_pixels = conv_x_pixels(1);
        kernel = ones(conv_y_pixels, conv_x_pixels) / (conv_y_pixels * conv_x_pixels);        
    else
        conv_y_pixels = conv_x_pixels(1);
        conv_x_pixels = conv_x_pixels(1);
        kernel = ones(conv_x_pixels, conv_x_pixels) / (conv_x_pixels * conv_x_pixels);
    end
else
    %    if Conv_x_pixels == Conv_y_pixels;
    %        kernel = fspecial('disk', Conv_y_pixels);
    %    else
    %        kernel = fspecial('average', [conv_y_pixels, conv_x_pixels]); % with image processing toolbox
    kernel = ones(conv_y_pixels, conv_x_pixels) / (conv_y_pixels * conv_x_pixels);
    %    end
end
grid_pad = [flipud(grid(1 : conv_y_pixels, :)); grid ; flipud(grid(end-conv_y_pixels : end, :))];
grid_pad = [fliplr(grid_pad(:, 1 : conv_x_pixels)) grid_pad fliplr(grid_pad(:, end - conv_x_pixels : end))];

GridPadS = conv2(grid_pad, kernel, 'same'); % old method (no image toolbox)
%GridPadS = imfilter(GridPad,kernel,'replicate'); % new method

grid_s = GridPadS(conv_y_pixels+1:end-conv_y_pixels-1,conv_x_pixels+1:end-conv_x_pixels-1);
