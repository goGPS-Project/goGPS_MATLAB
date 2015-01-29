function [dist,proj] = ref_2d_projection(ref,EAST,NORTH)

% SYNTAX:
%   [dist,proj] = ref_2d_projection(ref,EAST,NORTH)
%
% INPUT:
%   ref = reference path (X,Y,Z coordinates of the vertices)
%   EAST = estimated trajectory in UTM coordinates (EAST)
%   NORTH = estimated trajectory in UTM coordinates (NORTH)
%
% OUTPUT:
%   dist = 2D distance of each estimated point from the reference
%   proj = projected trajectory
%
% DESCRIPTION:
%   2D projection on a reference path.
%   At the moment working only for adjacency matrix as in the following
%   example:
%
%   adj_mat = [ 0 1 0 0 0 1
%               1 0 1 0 0 0
%               0 1 0 1 0 0
%               0 0 1 0 1 0
%               0 0 0 1 0 1
%               1 0 0 0 1 0 ];

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

%close the path by connecting the first and the last points
ref = [ref;ref(1,:)];

%computation of the angular coefficients
ax = ref(2:end,1) - ref(1:end-1,1);
ay = ref(2:end,2) - ref(1:end-1,2);

%normalization on the segment distance
ad = sqrt(ax.^2 + ay.^2);
ax = ax ./ ad;
ay = ay ./ ad;

%offset of the curvilinear coordinate
s0 = [0; cumsum(ad)];

for j = 1 : length(EAST)
    pos_R(1) = EAST(j);
    pos_R(2) = NORTH(j);

    d0 = sqrt((pos_R(1) - ref(:,1)).^2 + ...
              (pos_R(2) - ref(:,2)).^2);

    [dmin0 i0] = min(d0);

    %projection on the reference path
    bx = pos_R(1) - ref(1:end-1,1) + ax.*s0(1:end-1);
    by = pos_R(2) - ref(1:end-1,2) + ay.*s0(1:end-1);

    s_R = (ax.*bx + ay.*by) ./ (ax.^2 + ay.^2);

    pos_R_proj(:,1) = ref(1:end-1,1) + ax .* (s_R - s0(1:end-1));
    pos_R_proj(:,2) = ref(1:end-1,2) + ay .* (s_R - s0(1:end-1));

    %computation of the minimum distance
    d = sqrt((pos_R(1) - pos_R_proj(:,1)).^2 + ...
        (pos_R(2) - pos_R_proj(:,2)).^2);

    [dmin i] = min(d);

    %position in cartesian coordinates
    while (dmin < dmin0) & ((pos_R_proj(i,1) < min(ref(i,1),ref(i+1,1))) | (pos_R_proj(i,1) > max(ref(i,1),ref(i+1,1))) | ...
            (pos_R_proj(i,2) < min(ref(i,2),ref(i+1,2))) | (pos_R_proj(i,2) > max(ref(i,2),ref(i+1,2))))

        d(i) = 9e99;
        [dmin i] = min(d);

    end

    if dmin0 < dmin
        dist(j,1) = dmin0;

        proj(j,1) = ref(i0,1);
        proj(j,2) = ref(i0,2);

    else
        dist(j,1) = dmin;

        proj(j,1) = pos_R_proj(i,1);
        proj(j,2) = pos_R_proj(i,2);
    end
end