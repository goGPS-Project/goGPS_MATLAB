function [igp4, tv] = sel_igp(latpp, lonpp, igp, lat_igp, lon_igp)

% SYNTAX:
%   [igp4, tv] = sel_igp(latpp, lonpp, igp, lat_igp, lon_igp);
%
% INPUT:
%   latpp = ionosphere pierce point latitude [rad]
%   lonpp = ionosphere pierce point longitude [rad]
%   igp   = IGP IDs (vector)
%   lat_igp = IGP latitude  [deg] (vector)
%   lon_igp = IGP longitude [deg] (vector)
%
% OUTPUT:
%   igp4 = IDs of the IGPs of the cell that contains the pierce point
%   tv   = igp4 latitude,longitude coordinates [deg] (matrix)
%
% DESCRIPTION:
%   Selection of the 4 IGP surroudning the ionosphere pierce point.
%
%   WARNING: this is a simplified version of the procedure described in the
%            RTCA229C appendix; valid only for latitudes between 0 and 60
%            degrees (central-southern Europe).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
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

%WARNING: this is a simplified version of the procedure described in the
%         RTCA229C appendix; valid only for latitudes between 0 and 60 deg

v_lat = [0 : 5 : 60];
v_lon = [-180 : 5 : 180];

latpp = latpp * 180/pi;
lonpp = lonpp * 180/pi;

%find the lower latitude of the closest nodes
r_lati = find(v_lat <= latpp,1,'last');
lati = v_lat(r_lati);
%find the higher latitude of the closest nodes
r_lats = find(v_lat > latpp,1,'first');
lats = v_lat(r_lats);

%find the lower longitude of the closest nodes 
r_loni = find(v_lon <= lonpp,1,'last');
loni = v_lon(r_loni);
%find the higher longitude of the closest nodes
r_lons = find(v_lon > lonpp,1,'first');
lons = v_lon(r_lons);

%coordinates of the 4 points (Appendix A, pag. 41)
tv_1 = [lats lons];
tv_2 = [lats loni];
tv_3 = [lati loni];
tv_4 = [lati lons];

tv = [tv_1; tv_2; tv_3; tv_4];

%lat/lon matrix
m_ll = [lat_igp' lon_igp'];

%initialization
igp4 = NaN(1,4);

for i = 1 : 4

    [val, cc, r] = intersect(tv(i,:), m_ll, 'rows'); %#ok<ASGLU>
    %NOTE: it may be needed to use a triangular cell instead of a square
    %cell (see procedure on pag. 37, app. A, RTCA 229C)
    if (~isempty(r))

        igp4(i) = igp(r);
    end
    %igp4(i) = igp(r);
end
