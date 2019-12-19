function [ipp_lat, ipp_lon, elS] = IPP_satspec(elev_series, azim_series, commontime, stations_idx, PRN, pos_RM)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
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


n_sta=length(elev_series);

ipp_lat = NaN(n_sta,length(commontime));
ipp_lon = NaN(n_sta,length(commontime));

for s = 1 : n_sta
    %extract azimuth
    azS = azim_series{s}(PRN,stations_idx(s,:));
    pos = find(azS > 180);
    azS(pos) = azS(pos) - 360;

    %extract elevation
    elS = elev_series{s}(PRN,stations_idx(s,:));

    azS = azS*pi/180;
    elS = elS*pi/180;

    [latR, lonR] = cart2geod(pos_RM(1,1,s), pos_RM(2,1,s), pos_RM(3,1,s));

    for e = 1 : length(commontime)
        [latpp, lonpp] = iono_pierce_point(latR, lonR, azS(e), elS(e));
        ipp_lat(s,e) = latpp*180/pi;
        ipp_lon(s,e) = lonpp*180/pi;
    end
end
