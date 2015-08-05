function invA = cholinv(A)

% SYNTAX:
%   invA = cholinv(A);
%
% INPUT:
%   A = positive definite matrix to be inverted
%
% OUTPUT:
%   invA = inverse of matrix A
%
% DESCRIPTION:
%   Inverse of a positive definite matrix by using Cholesky decomposition.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Andrea Nardo
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

% compute cholesky decomposition
n = size(A,1);
U    = chol(A);
invU = U\eye(n);
%L    = inv(L);
invA = invU * invU';
