% arrow - Arrow plot 
%
% INPUT
%   x0                arrow start X
%   y0                arrow start Y
%   vx                X modulus
%   vy                Y modulus
%   h_sscale          arrow scale (default 1)
%   color             color of the arrow (default [1 1 1])
%   <flag_oversize>   if true (default false) it will add withe dashed line on the arrow
%
%
% OUTPUT
%   hl            list of handles (one line, one patch)
%   hp            list of handle  (the point)
%
% SYNTAX
%   hl = arrow;
%
% DISCLAMER
%   Set axis ratio BEFORE using arrow, uses triPlot
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
%  Contributors:     Andrea Gatti, ...
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

function [hl, hp] = arrow(x0, y0, vx, vy, h_scale, color, flag_oversize, flag_dot)
    if nargin < 5 || isempty(h_scale)
        h_scale = 10;
    end
    if nargin < 6 || isempty(color)
        color = [0.6 0.6 0.6];
    end
    if nargin < 7 || isempty(flag_oversize)
        flag_oversize = false;
    end
    if nargin < 8 || isempty(flag_dot)
        flag_dot = false;
    end
    
    hl = []; % handle to the arrow head and body
    hp = []; % handle to the point
    
    if any([vx vy])
        line_width = 4;
        
        % new vertex
        
        ax = gca;
        old_state = ax.Units;
        ax.Units = 'points';
        hl = [];
        try
            ax.PlotBox;
        catch
            hl = line(x0 + [0 vx], y0 + [0 vy], 'LineWidth', line_width, 'Color', color); drawnow; % draw the full line
            cb = colorbar; drawnow; ax.PlotBox; delete(cb); % force PlotBox property creation
        end
        scale_x = diff(xlim) / ax.PlotBox(3); % points to X-axis coordinates
        scale_y = diff(ylim) / ax.PlotBox(4); % points to Y-axis coordinates
                
        aspect_ratio = pbaspect;
        az = cart2sph(vx * (aspect_ratio(1) * diff(ylim)), -vy * (aspect_ratio(2) * diff(xlim)), 0) * 180/pi + 90;
        aw_size = 1;
        
        red_x = - sind(az) * 0.5 * h_scale * scale_x;
        red_y = - cosd(az) * 0.5 * h_scale * scale_y;
        
        x1 = x0 + vx + red_x;
        y1 = y0 + vy + red_y;
        
        if ~isempty(hl)
            delete(hl);
        end
        hl = line(x0 + [0 vx + red_x], y0 + [0 vy + red_y], 'LineWidth', line_width, 'Color', color);
        hl = [hl hl];
        if (flag_oversize)
            hl = line(x0 + [0 vx + red_x], y0 + [0 vy + red_y], 'LineWidth', line_width, 'Color', [1 1 1], 'LineStyle', '--');
            hl = [hl hl];
        end
        
        if flag_dot
            hp = plot(x0, y0,'o', 'MarkerSize', 10, 'MarkerEdgeColor', max(0, Core_UI.LBLUE - 0.1), 'MarkerFaceColor', min(1, Core_UI.LBLUE + 0.1), 'LineWidth', 2); hold on;
        end
                
        tr = triPlot(x1, y1, az, aw_size, aw_size, h_scale);
        tr.FaceColor = color;
        tr.EdgeColor = 'none';
        ax.Units = old_state;
                
        hl = [hl tr];
    end
end
