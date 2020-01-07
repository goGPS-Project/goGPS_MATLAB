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

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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

ph_GF = zeros(size(ph1));
for s = 1 : size(ph1,1)
    if (any(ph1(s,:)) && any(ph2(s,:)))
        freq = Core_Utils.V_LIGHT ./ lambda(s,1:2);
        index_1 = find(ph1(s,:) ~= 0);
        index_2 = find(ph2(s,:) ~= 0);
        index = intersect(index_1, index_2);
        if (any(err_iono(s,:)))
            index_3 = find(err_iono(s,:) ~= 0);
            index = intersect(index,   index_3);
            [p,~,mu] = polyfit(index,err_iono(s,index),3);
            err_iono_fit = polyval(p,index,[],mu);
            corr = ((freq(1,1)^2-freq(1,2)^2)/freq(1,2)^2)*err_iono_fit;
        else
            corr = zeros(size(err_iono(s,index)));
        end
        ph_GF(s,index) = (lambda(s,1)*ph1(s,index) - lambda(s,2)*ph2(s,index)) - corr;
    end
end
