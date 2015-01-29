function [Cyy] = global2localCov(Cxx, X)

% SYNTAX:
%   [Cyy] = global2localCov(Cxx, X);
%
% INPUT:
%   Cxx = input covariance matrices
%   X   = position vectors
%
% OUTPUT:
%   Cyy = output covariance matrices
%
% DESCRIPTION:
%   Covariance propagation from Earth-fixed reference frame to local-level reference frame

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

%initialize new covariance matrices
Cyy = zeros(size(Cxx));

for i = 1 : size(X,2)

    %geodetic coordinates
    [phi, lam, h] = cart2geod(X(1,i), X(2,i), X(3,i)); %#ok<NASGU>

    %rotation matrix from global to local reference system
    R = [-sin(lam) cos(lam) 0;
         -sin(phi)*cos(lam) -sin(phi)*sin(lam) cos(phi);
         +cos(phi)*cos(lam) +cos(phi)*sin(lam) sin(phi)];

    %covariance propagation
    Cyy(:,:,i) = R * Cxx(:,:,i) * R';
end

%----------------------------------------------------------------------------------------------