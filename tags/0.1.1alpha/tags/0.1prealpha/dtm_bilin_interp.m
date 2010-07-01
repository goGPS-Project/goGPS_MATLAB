function [height] = dtm_bilin_interp(E_approx, N_approx, dtm, ncols, nrows, cellsize, Ell, Nll, nodata)

% SYNTAX:
%   [height] = dtm_bilin_interp(E_approx, N_approx, dtm, ncols, nrows, cellsize, Ell, Nll, nodata);
%
% INPUT:
%   E_approx = EAST coordinate of the interpolation point
%   N_approx = NORTH coordinate of the interpolation point
%   dtm = matrix containing the dtm grid
%   ncols = number of columns of the grid
%   nrows = number of rows of the grid
%   cellsize = ground size of a cell
%   xll = X coordinate of the center of the lower left cell
%   yll = Y coordinate of the center of the lower left cell
%   nodata = value used for cells not containing data
%
% OUTPUT:
%   height = interpolated height value
%
% DESCRIPTION:
%   Function that applies a bilinear interpolation of the four nearest nodes
%   of a DTM in correspondence of a point of given coordinates.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
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

%preparation of the axes of the georeferenced grid
E = (Ell : cellsize : Ell + (ncols - 1) * cellsize)';
N = (Nll : cellsize : Nll + (nrows - 1) * cellsize)';

if (E_approx <= E(1) | E_approx >= E(end) | N_approx <= N(1) | N_approx >= N(end))
    height = nodata;
    return
end;

%detection of the grid node nearest to the interpolation point
[mE, posE] = min(abs(E - E_approx));
[mN, posN] = min(abs(N - N_approx));

%definition of the four grid nodes that sorround the interpolation point
% (i,j) image coordinates (upper-left origin)
% (E,N) ground coordinates (bottom-left origin)
if (mE ~=0)
    if (E(posE) > E_approx)
        j_left = posE - 1;
        j_right = posE;
        E_left = E(posE - 1);
    else
        j_left = posE;
        j_right = posE + 1;
        E_left = E(posE);
    end;
else
    height = nodata;
    return
end;
if (mN ~=0)
    if (N(posN) > N_approx)
        i_up = nrows + 1 - posN;
        i_down = i_up + 1;
        N_down = N(posN - 1);
    else
        i_down = nrows + 1 - posN;
        i_up = i_down - 1;
        N_down = N(posN);
    end;
else
    height = nodata;
    return
end;

%if one of the height values of the four sorrounding points is a nodata value, do not interpolate and return nodata value
if (dtm(i_up,j_left) == nodata | dtm(i_up,j_right) == nodata | dtm(i_down,j_left) == nodata | dtm(i_down,j_right) == nodata)
    height = nodata;
    return
end

%computation of the parameters of the bilinear function
%f(E, N) = a*E*N + b*E + c*N + d

A = [0 0 cellsize 1; ...
     cellsize^2 cellsize cellsize 1; ...
     0 0 0 1; ...
     0 cellsize 0 1];

B = [dtm(i_up,j_left);...
     dtm(i_up,j_right);...
     dtm(i_down,j_left);...
     dtm(i_down,j_right)];

bilin_param = linsolve(A,B);

i_approx = N_approx - N_down;
j_approx = E_approx - E_left;

%computation of the interpolated value
height = bilin_param(1) * j_approx * i_approx + bilin_param(2) * j_approx + bilin_param(3) * i_approx + bilin_param(4);
height = double(height);
