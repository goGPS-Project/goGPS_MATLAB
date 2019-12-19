function  [Z,residual] = plane_fitting(x, y, z, xi, yi)

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

%planar fitting (least squares)
X = [ones(size(x)) x y]; %design matrix
a = X\z;                    %parameters

Z0 = X*a;                    %fitted values
% plot3(x,y,Z0,'go','MarkerSize',7)
%plot fitting plane
min_x1 = min(x);
max_x1 = max(x);
min_x2 = min(y);
max_x2 = max(y);
step = 10;
residual = z - Z0;
x = min_x1 : (max_x1 - min_x1)/step : max_x1;
y = min_x2 : (max_x2 - min_x2)/step : max_x2;
[X1,X2] = meshgrid(x,y);

Z = zeros(size(X1));

%evaluate the plane at the requested coordinates
X = [ones(size(xi)) xi yi];
Z = X*a;


