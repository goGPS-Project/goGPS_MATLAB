function h = polarImagesc(az_grid, decl_grid, data, plot_bg)
% Similar to polarImagesc but for maps
% SYNTAX:
%    polarImagesc(az, decl, point_size, color, <flag_plot_bg>)
%
% INPUT
%   az      azimuth      [rad]
%   decl    declination  [rad]
%   color   data field for scatter
%
% OUTPUT
%   h       handle to the scattered points
%
% SEE ALSO
%   polarScatter

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2020 Andrea Gatti, Giulio Tagliaferro, Eugenio Realini
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
%----------------------------------------------------------------------------------------------


%%% INTERNAL PARAMETER
    scale = 1;
    %%%

    decl_n = decl_grid/(pi/2)*scale;
    x = repmat(sin(az_grid(:)'), numel(decl_n), 1) .* repmat(decl_n(:), 1, numel(az_grid));
    y = repmat(cos(az_grid(:)'), numel(decl_n), 1) .* repmat(decl_n(:), 1, numel(az_grid));
    %scatter(x(data~=0),y(data~=0),20,data(data~=0),'filled')
    
    dataInterp = scatteredInterpolant(x(:), y(:), data(:), 'linear' );
    x = -1 : 0.005 : 1;
    y = x;
    [x_mg, y_mg] = meshgrid(x, y);
    polar_data = nan(numel(x), numel(y));
    id_ok = hypot(x_mg, y_mg) < 1;
    polar_data(id_ok) = dataInterp(x_mg(id_ok), y_mg(id_ok));
    h = imagesc(x, y, polar_data);
    h.AlphaData = id_ok; 
    h.AlphaData(isnan(polar_data)) = 0; 
    set(gca,'Ydir','normal');
    
    is_hold = ishold();
    if nargin < 4
        plot_bg = ~is_hold;
    end
    if plot_bg
        hold on
        %plot parallel
        az_l = [0:pi/200:2*pi];
        d_step = 15/180*pi;
        decl_s = ([0:d_step:pi/2]/(pi/2))*scale;
        for d = decl_s
            x = cos(az_l).*d;
            y = sin(az_l).*d;
            plot(x,y,'color',[0.6 0.6 0.6]);
            text(cos(80/180*pi)*d,sin(80/180*pi)*d,sprintf('%d',round(d*90)),'HorizontalAlignment','center', 'FontWeight', 'bold');            
        end
        %plot meridian
        az_step = 30/180 *pi;
        az_s = [0:az_step:2*pi];
        decl_l = ([0 1])*scale;
        for a = az_s
            x = cos(a).*decl_l;
            y = sin(a).*decl_l;
            plot(x,y,'color',[0.6 0.6 0.6]);
            if abs(a-2*pi) > 0.0001
                text(cos(a)*1.1,sin(a)*1.1,sprintf('%d', mod(round((2*pi - a + pi/2) / pi * 180), 360)), 'HorizontalAlignment','center', 'FontWeight', 'bold');
            end
        end
        axis equal
        % xlim([-2 2])
        % ylim([-2 2])
        axis off
        set(gcf,'color','w');
        if ~is_hold
            hold off
        end
        xlim([-1.15 1.15]); ylim([-1.15 1.15]);
    end
end
