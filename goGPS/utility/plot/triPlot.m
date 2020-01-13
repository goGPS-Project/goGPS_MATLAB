% triPlot - Triangle plot 
%
% INPUT
%   x             list of triangles centroids (x-axes)    [n x 1]
%   y             list of triangles centroids (y-axes)    [n x 1]
%   ang           ang [0-360] from Up -> clockwise        [n x 1]
%   modulus       list of triangle magnitude              [n x 1]
%   max_size      maximum size of the triangle            [1 x 1]
%                 w.r.t. y axes (if scales are not equal) 
%                 default = 1
%   marker_scale  scale value                             [1 x 1]
%                 default = max_size / max(modulus)
%
% SYNTAX
%   triPlot(x, y, ang, modulus, <max_size>, <marker_scale>);
%
% DISCLAMER
%   Set axis ratio BEFOR using triPlot
%
% EXAMPLE
%   ang = 0:25:720;
%   x = (1 : numel(ang));
%   y = ones(1, numel(point_x)); 
%   modulus = ang;
%   max_size = 1;
%   marker_scale = max_size / max(modulus);
% 
%   figure;
%   triPlot(x, -2 * y, -90*1+0*ang, modulus, max_size, marker_scale);
%   triPlot(x, 0 * y, ang, fliplr(modulus), max_size, marker_scale);
%   triPlot(x, 2 * y, 90*1+0*ang, modulus, max_size, marker_scale);
%   axis equal
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
function triPlot(x, y, ang, modulus, max_size, marker_scale)
    if nargin < 5
        max_size = 1;
    end
    if nargin < 6
        marker_scale = max_size / max(modulus);
    end
    
    % Type conversion
    x = double(x);
    y = double(y);
    ang = double(ang);
    modulus = double(modulus);
    max_size = double(max_size);
    marker_scale = double(marker_scale);    
    
    min_size = max_size * 0.6;
        
    % Normalize point scale
    scale_factor = double(max(min_size, min(max_size, (modulus .* marker_scale))));
    ang = -ang;
    
    % Define arrow
    ax = gca;
    
    plot_scale_factor = ax.PlotBoxAspectRatio(2) * diff(xlim) / diff(ylim);
    
    x_tri = [-0.6 0 0.6 0];
    y_tri = [-sqrt(1 - x_tri(1)^2) 1 -sqrt(1 - x_tri(3)^2) -0.5];
    
    % Radius 1 -> 0.5
    x_tri = x_tri ./ 2;
    y_tri = y_tri ./ 2;
    
    hold on;
    % Fix aspect ratio of the axes    
    for i = 1 : numel(x)
        col = modulus(i);
        pos_x = x_tri .* scale_factor(i) + x(i);
        pos_y = y_tri .* scale_factor(i) + y(i);
        ht = patch(pos_x, pos_y, col); hold on;
        ht.Vertices(:,3) = 1000;
        rotate(ht, [0, 0, 1], double(ang(i)), [x(i), y(i) 1]);
        
        % Rescale vertexes if the ratio of the axis is not ok
        offset = repmat([x(i) y(i) mean(ht.Vertices(:,3))], numel(x_tri), 1);
        vrtx = ht.Vertices - offset; % translate to center of axis
        vrtx(:, 1) = vrtx(:, 1) .* plot_scale_factor;   % scale
        ht.Vertices = vrtx + offset; % translate back
    end
end

