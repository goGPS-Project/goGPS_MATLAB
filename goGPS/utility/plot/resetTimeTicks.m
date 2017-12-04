function resetTimeTicks(h, num, format)
% reset function for setTimeTicks
% SEE ALSO: setTimeTicks

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
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

    if (nargin == 1)
        num = 4;
        format = 'dd/mm/yyyy HH:MMPM';
    end
    if isa(h, 'matlab.ui.Figure')
        h = get(h,'Children');
    end
    for i = 1 : numel(h)
        try
            ax = axis(h(i));
            step=(ax(2)-ax(1))/(num);
            tickPos = (ax(1)+(step/2)):step:(ax(2)-(step/2));
            set(h(i), 'XTick', tickPos);
            datetick(h(i),'x',format,'keepticks');
            axis(h(i),ax);
        catch
        end
    end
end
