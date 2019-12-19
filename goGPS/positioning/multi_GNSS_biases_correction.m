function [pr1, ph1, pr2, ph2] = multi_GNSS_biases_correction(time_ref, pr1, ph1, pr2, ph2, ISBs, Eph, constellations, lambda)

% SYNTAX:
%   [pr1, ph1, pr2, ph2] = multi_GNSS_biases_correction(time_ref, pr1, ph1, pr2, ph2, ISBs, Eph, constellations, lambda);
%
% INPUT:
%   time_ref = GPS reference time
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   ISBs = inter-system/-frequency biases (as output by pre_processing.m)
%   Eph = matrix containing 33 ephemerides for each satellite
%   constellations = struct with multi-constellation settings
%   lambda  = wavelength matrix (depending on the enabled constellations)
%
% OUTPUT:
%   pr1 = processed code observation (L1 carrier)
%   ph1 = processed phase observation (L1 carrier)
%   pr2 = processed code observation (L2 carrier)
%   ph2 = processed phase observation (L2 carrier)
%
% DESCRIPTION:
%   Correction of code and phase observations based on inter-system and inter-frequency biases.

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

if (~isempty(ISBs))
    [~, idx] = unique(constellations.systems,'stable');
    offset = 0;
    for c = 2 : length(idx)
        if (c~=length(idx))
            e = idx(c+1)-1;
        else
            e = size(pr1,1);
        end
        if (strcmp(constellations.systems(idx(c)), 'R'))
            for r = 0 : constellations.GLONASS.numSat - 1
                k = find_eph(Eph, idx(c)+r, time_ref(round(length(time_ref)/2)), 86400);
                if (~isempty(k) && all(lambda(idx(c)+r,:)))
                    slot = Eph(15,k) + 8;
                    pr1(idx(c)+r,:) = pr1(idx(c)+r,:) - ISBs(c-2+slot)*Core_Utils.V_LIGHT*(pr1(idx(c)+r,:)~=0);
                    pr2(idx(c)+r,:) = pr2(idx(c)+r,:) - ISBs(c-2+slot)*Core_Utils.V_LIGHT*(pr2(idx(c)+r,:)~=0);
                    ph1(idx(c)+r,:) = ph1(idx(c)+r,:) - ISBs(c-2+slot)*Core_Utils.V_LIGHT*(ph1(idx(c)+r,:)~=0)/lambda(idx(c)+r,1);
                    ph2(idx(c)+r,:) = ph2(idx(c)+r,:) - ISBs(c-2+slot)*Core_Utils.V_LIGHT*(ph2(idx(c)+r,:)~=0)/lambda(idx(c)+r,2);
                else
                    pr1(idx(c)+r,:) = 0;
                    pr2(idx(c)+r,:) = 0;
                    ph1(idx(c)+r,:) = 0;
                    ph2(idx(c)+r,:) = 0;
                end
            end
            offset = slot;
        else
            pr1(idx(c):e,:) = pr1(idx(c):e,:) - ISBs(c-1+offset)*Core_Utils.V_LIGHT*(pr1(idx(c):e,:)~=0);
            pr2(idx(c):e,:) = pr2(idx(c):e,:) - ISBs(c-1+offset)*Core_Utils.V_LIGHT*(pr2(idx(c):e,:)~=0);
            ph1(idx(c):e,:) = ph1(idx(c):e,:) - ISBs(c-1+offset)*Core_Utils.V_LIGHT*(ph1(idx(c):e,:)~=0)/max(lambda(idx(c):e,1));
            ph2(idx(c):e,:) = ph2(idx(c):e,:) - ISBs(c-1+offset)*Core_Utils.V_LIGHT*(ph2(idx(c):e,:)~=0)/max(lambda(idx(c):e,2));
        end
    end
end

