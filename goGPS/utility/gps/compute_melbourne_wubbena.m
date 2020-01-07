function [ph_MW, pr_NL, ph_WL] = compute_melbourne_wubbena(ph1, ph2, pr1, pr2, lambda)

% SYNTAX:
%   [ph_MW, pr_NL, ph_WL] = compute_melbourne_wubbena(ph1, ph2, pr1, pr2, lambda);
%
% INPUT:
%   ph1 = carrier phase observation on the first frequency [cycles]
%   ph2 = carrier phase observation on the second frequency [cycles]
%   pr1 = code observation on the first frequency [m]
%   pr2 = code observation on the second frequency [m]
%   lambda = wavelength matrix (depending on the enabled constellations) [m]
%
% OUTPUT:
%   ph_MW = Melbourne-Wubbena observable
%   pr_NL = narrow-lane observable
%   ph_WL = wide-lane observable
%
% DESCRIPTION:
%   Compute the Melbourne-Wubbena observable.

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

ph_WL = zeros(size(ph1));
pr_NL = zeros(size(ph1));
ph_MW = zeros(size(ph1));
for s = 1 : size(ph1,1)
    if (any(ph1(s,:)) && any(ph2(s,:)))
        freq = Core_Utils.V_LIGHT ./ lambda(s,1:2);
        index_1 = find(ph1(s,:) ~= 0);
        index_2 = find(ph2(s,:) ~= 0);
        index = intersect(index_1, index_2);
        ph_WL(s,index) = (freq(1,1)*lambda(s,1)*ph1(s,index) - freq(1,2)*lambda(s,2)*ph2(s,index))/(freq(1,1) - freq(1,2));
        pr_NL(s,index) = (freq(1,1)*pr1(s,index) + freq(1,2)*pr2(s,index))/(freq(1,1) + freq(1,2));
        ph_MW(s,index) = ph_WL(s,index) - pr_NL(s,index);
    end
end
