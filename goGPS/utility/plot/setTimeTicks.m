function setTimeTicks(num,format,ax)
% SYNTAX:
%   setTimeTicks(num,format);
%   setTimeTicks(ax, num,format);
%
% EXAMPLE:
%   setTimeTicks(4,'dd/mm/yyyy HH:MM');
%
% INPUT:
%   num      = number of tick to be shown
%   format   = time format => see datetick for more informations
%
% DESCRIPTION:
%   display a number of tick = num
%   showing the UTC date in the format required
%   the X axis should contain time values in UTC format
%
% REMARKS:
%  setTimeTicks uses the 'tag' and 'userdata' properties of the figure.
%  zoom and pan callbacks are used
%  (whenever another callback is present, it is appended - works with qplot)
%
% SEE ALSO:
%   resetTimeTicks qplot

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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

    if nargin < 1
        num = 4;
    end
    
    if nargin < 2 || isempty(format)
        format = 'auto';
    end

    if (nargin == 3)
        h = num;
        num = format;
        format = ax;
    else
        h = gca;
    end
    ttData.num = num;
    ttData.format = format;

    hZoom = zoom(h);
    hPan = pan(ancestor(h,'figure'));

    %     ttData.hZoom = get(hZoom,'ActionPostCallback');
    %     ttData.hPan = get(hPan,'ActionPostCallback');

    set(h, 'userdata', ttData);

    set(hZoom, 'ActionPostCallback', @zoomCallback);
    set(hPan , 'ActionPostCallback', @panCallback);

    resetTimeTicks(h, num, format);
end

function zoomCallback(obj, ev)
    resetTimeTicks(ev.Axes, ev.Axes.UserData.num, ev.Axes.UserData.format);
%     if ~isempty(obj.UserData.hZoom)
%         obj.UserData.hZoom(obj,ev);
%     end
end

function panCallback(obj, ev)
    resetTimeTicks(ev.Axes, ev.Axes.UserData.num, ev.Axes.UserData.format);
%     if ~isempty(obj.UserData.hPan)
%         obj.UserData.hPan(obj,ev);
%     end
end
