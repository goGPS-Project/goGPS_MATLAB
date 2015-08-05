function [interp_value] = grid_bilin_interp(X_approx, Y_approx, grid, ncols, nrows, cellsize, Xll, Yll, nodata)

% SYNTAX:
%   [interp_value] = grid_bilin_interp(X_approx, Y_approx, grid, ncols, nrows, cellsize, Xll, Yll, nodata);
%
% INPUT:
%   X_approx = X coordinate of the interpolation point
%   Y_approx = Y coordinate of the interpolation point
%   grid = matrix containing the grid
%   ncols = number of columns of the grid
%   nrows = number of rows of the grid
%   cellsize = ground size of a cell
%   Xll = X coordinate of the center of the lower left cell
%   Yll = Y coordinate of the center of the lower left cell
%   nodata = value used for cells not containing data
%
% OUTPUT:
%   interp_value = interpolated value
%
% DESCRIPTION:
%   Function that applies a bilinear interpolation of the four nearest nodes
%   of a georeferenced grid in correspondence of a point of given coordinates.

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

%preparation of the grid axes
X = (Xll : cellsize : Xll + (ncols - 1) * cellsize)';
Y = (Yll : cellsize : Yll + (nrows - 1) * cellsize)';

if (X_approx <= X(1) | X_approx >= X(end) | Y_approx <= Y(1) | Y_approx >= Y(end))
    interp_value = nodata;
    return
end

%detection of the grid node nearest to the interpolation point
[mX, posX] = min(abs(X - X_approx));
[mY, posY] = min(abs(Y - Y_approx));

%definition of the four grid nodes that sorround the interpolation point
% (i,j) image coordinates (upper-left origin)
% (X,Y) ground coordinates (bottom-left origin)
if (X(posX) > X_approx) | (mX ==0)
    j_left = posX - 1;
    j_right = posX;
    X_left = X(posX - 1);
else
    j_left = posX;
    j_right = posX + 1;
    X_left = X(posX);
end

if (Y(posY) > Y_approx) | (mY ==0)
    i_up = nrows + 1 - posY;
    i_down = i_up + 1;
    Y_down = Y(posY - 1);
else
    i_down = nrows + 1 - posY;
    i_up = i_down - 1;
    Y_down = Y(posY);
end

%if one of the interp_value values of the four sorrounding points is a nodata value, do not interpolate and return nodata value
if (grid(i_up,j_left) == nodata | grid(i_up,j_right) == nodata | grid(i_down,j_left) == nodata | grid(i_down,j_right) == nodata)
    interp_value = nodata;
    return
end

%computation of the parameters of the bilinear function
%f(X, Y) = a*X*Y + b*X + c*Y + d

A = [0 0 cellsize 1; ...
     cellsize^2 cellsize cellsize 1; ...
     0 0 0 1; ...
     0 cellsize 0 1];

B = [grid(i_up,j_left);...
     grid(i_up,j_right);...
     grid(i_down,j_left);...
     grid(i_down,j_right)];

bilin_param = A\B;

i_approx = Y_approx - Y_down;
j_approx = X_approx - X_left;

%computation of the interpolated value
interp_value = bilin_param(1) * j_approx * i_approx + bilin_param(2) * j_approx + bilin_param(3) * i_approx + bilin_param(4);
interp_value = double(interp_value);
