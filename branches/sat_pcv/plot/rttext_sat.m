function rttext_sat (t, az, el, snr, obs, pivot, Eph, SP3)

% SYNTAX:
%   rttext_sat (t, az, el, obs, pivot, Eph, SP3);
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
%   Eph = matrix containing 33 navigation parameters for each satellite
%   SP3 = structure containing precise ephemeris data
%
% DESCRIPTION:
%   Real time textual display of satellite data.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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
    text(0.0,0.95,'PRN   ELEV   AZIM     SNR');
end

%----------------------------------------------------------------------------------------------
% DISPLAY SATELLITE CONFIGURATION
%----------------------------------------------------------------------------------------------

vert_pos = 0.95;

num_sat = numel(el);

sys = zeros(num_sat,1);
prn = zeros(num_sat,1);

sat = 1:numel(el);
sat = sat'.*abs(obs);

if (isempty(SP3))
    eph_avail = Eph(30,:);
    sys_avail = Eph(31,:);
    prn_avail = Eph(1,:);
else
    eph_avail = SP3.avail';
    sys_avail = SP3.sys';
    prn_avail = SP3.prn';
end

[~, idx1, idx2] = intersect(eph_avail, sat);

sys(idx2) = char(sys_avail(idx1));
prn(idx2) = prn_avail(idx1);

for i = 1 : num_sat

    if (el(i) > 0)

        sat_string = [sys(i) sprintf('%02d',prn(i)) '     ' sprintf('%04.1f',el(i)) '    ' sprintf('%05.1f',az(i)) '    ' sprintf('%04.1f',snr(i))];

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

        for j = i+1 : num_sat
            if (satid(j) > 0)
                pos = get(satid(j),'Position');
                set(satid(j),'Position',[0 pos(2)+0.05 0]);
            end
        end

    end
end
