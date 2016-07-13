function [ph_GF] = compute_geometry_free(ph1, ph2, lambda, err_iono)

% SYNTAX:
%   [ph_GF] = compute_geometry_free(ph1, ph2, lambda, err_iono);
%
% INPUT:
%   ph1 = carrier phase observation on the first frequency [cycles]
%   ph2 = carrier phase observation on the second frequency [cycles]
%   lambda = wavelength matrix (depending on the enabled constellations) [m]
%   err_iono = ionospheric error correction [m]
%
% OUTPUT:
%   ph_GF = geometry free observable
%
% DESCRIPTION:
%   Compute the geometry free observable.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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

ph_GF = zeros(size(ph1));
for s = 1 : size(ph1,1)
    if (any(ph1(s,:)) && any(ph2(s,:)))
        index_1 = find(ph1(s,:) ~= 0);
        index_2 = find(ph2(s,:) ~= 0);
        index = intersect(index_1, index_2);
        if (any(err_iono(s,:)))
            index_3 = find(err_iono(s,:) ~= 0);
            index = intersect(index,   index_3);
            p = polyfit(index,err_iono(s,index),3);
            err_iono_fit = polyval(p,index);
            corr = ((goGNSS.F1^2-goGNSS.F2^2)/goGNSS.F2^2)*err_iono_fit;
        else
            corr = zeros(size(err_iono(s,index)));
        end
        ph_GF(s,index) = (lambda(s,1)*ph1(s,index) - lambda(s,2)*ph2(s,index)) - corr;
    end
end
