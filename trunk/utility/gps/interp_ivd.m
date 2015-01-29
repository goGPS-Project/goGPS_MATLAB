function [ivd_pp] = interp_ivd(igp4, igp, r_ivd, latpp, lonpp, tv)

% SYNTAX:
%   [ivd_pp] = interp_ivd(igp4, igp, r_ivd, latpp, lonpp, tv);
%
% INPUT:
%   igp4  = IDs of the IGPs of the cell that contains the pierce point
%   igp   = IGP IDs (vector)
%   r_ivd = ionosphere vertical delays (vector)
%   latpp = ionosphere pierce point latitude [rad]
%   lonpp = ionosphere pierce point longitude [rad]
%   tv    = igp4 latitude,longitude coordinates [deg] (matrix)
%
% OUTPUT:
%   ivd_pp = ionophere vertical delay at the piercing point
%
% DESCRIPTION:
%   Interpolation of the ionosphere vertical delay at the piercing point.
%   Valid for piercing points between -85 and +85 degrees.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

latpp = latpp * 180/pi;
lonpp = lonpp * 180/pi;

ivd4 = NaN(1,4);

for i = 1 : 4 
    c = find(igp == igp4(i));
    ivd4(i) = r_ivd(c); %#ok<FNDSB>
end

%check how many valid ivd values are available
%flag = sum(isfinite(ivd4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%WARNING: the interpolation might be done on 4 points (square cell) or 3
%points (triangular cell); for now only the 4 points case is implemented.
%If only 3 points are available, their mean is used as a virtual fourth
%point (TO BE FIXED)
mean_ivd4 = mean(ivd4(isfinite(ivd4)));
ivd4(~isfinite(ivd4)) = mean_ivd4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%boundaries
lati = tv(3,1); %lower latitude
loni = tv(3,2); %lower longitude
lats = tv(1,1); %higher latitude
lons = tv(1,2); %higher longitude

d_latpp = latpp - lati;
d_lonpp = lonpp - loni;

xpp = d_lonpp / (lons - loni);
ypp = d_latpp / (lats - lati);

%weight functions
w1 = xpp * ypp;
w2 = (1 - xpp) * ypp;
w3 = (1 - xpp) * (1 - ypp);
w4 = xpp * (1 - ypp);

w = [w1 w2 w3 w4];

ivd_pp = sum(ivd4 .* w);
