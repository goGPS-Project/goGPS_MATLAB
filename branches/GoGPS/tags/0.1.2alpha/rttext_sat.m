function rttext_sat (t, az, el, snr, obs, pivot)

% SYNTAX:
%   rttext_sat (t, az, el, obs, pivot);
%
% INPUT:
%   t   = survey time (t=1,2,...)
%   az  = azimuth     [degrees]
%   el  = elevation   [degrees]
%   snr = signal-to-noise ratio
%   obs = observation type
%            0 = not used
%           +1 = code & phase
%           -1 = only code
%   pivot = pivot satellite
%
% DESCRIPTION:
%   Real time textual display of satellite data.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
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

global satid

% location on the screen
subplot(2,3,[3 6])

%----------------------------------------------------------------------------------------------
% TEXT CONTAINER
%----------------------------------------------------------------------------------------------

if (t == 1)
    % turn axes invisible
    axis off;
    tx = text(0.0,1,'Satellite configuration');
    set(tx,'FontWeight','Bold');
    text(0.0,0.95,'PRN  ELEV   AZIM     SNR');
end

%----------------------------------------------------------------------------------------------
% DISPLAY SATELLITE CONFIGURATION
%----------------------------------------------------------------------------------------------

vert_pos = 0.95;

for i = 1 : 32

    if (el(i) > 0)

        sat_string = [sprintf('%02d',i) '     ' sprintf('%05.2f',el(i)) '    ' sprintf('%06.2f',az(i)) '   ' sprintf('%d',snr(i))];

        vert_pos = vert_pos - 0.05;

        if (satid(i) == 0)
            satid(i) = text(0.0,vert_pos,sat_string);
        else
            set(satid(i),'String',sat_string)
            set(satid(i),'Position',[0 vert_pos 0])
        end

        if (obs(i) == 0); set(satid(i),'Color',[0.6 0.6 0.6]); end
        if (obs(i) == 1); set(satid(i),'Color','b'); end
        if (obs(i) == -1); set(satid(i),'Color','g'); end
        if (i == pivot); set(satid(i), 'Color', 'm'); end

    elseif (satid(i) > 0)

        delete(satid(i))
        satid(i) = 0;

        for j = i+1 : 32
            if (satid(j) > 0)
                pos = get(satid(j),'Position');
                set(satid(j),'Position',[0 pos(2)+0.05 0]);
            end
        end

    end
end
