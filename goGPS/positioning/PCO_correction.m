function PCO_corr = PCO_correction(antenna_PCV, XR, XS, sys, frequency)

% SYNTAX:
%   PCO_corr = PCO_correction(antenna_PCV, XR, XS, sys, frequency);
%
% INPUT:
%   antenna_PCV = antenna PCO/PCV struct
%   XR = zenithal angle to be interpolated (deg, array)
%   XS = azimutal angle to be interpolated (deg, array)
%   sys = system identification code to match correct constellation (1: GPS, 2: GLONASS, ...)
%   frequency = frequency of the observations (i.e.: 1 for L1, 2 for L2, ...)
%
% OUTPUT:
%   PCO_corr = value of the PCO correction projected along the receiver-satellite line-of-sight
%
% DESCRIPTION:
%   Function that computes the PCO corrections for the satellites

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
%  Contributors:     Stefano Caldera, ...
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

PCO_corr   = zeros(size(sys));

%extract constellations that have to be analyzed
constellation=unique(floor(sys));

%compute PCO correction: loop on every requested constellation
for j = 1 : length(constellation)

    %current system
    sys_i=constellation(j);

    %verify if corrections for this sys+frequency are available
    sysfreq_i = sys_i*10 + frequency;
    index_freq = find(antenna_PCV.sysfreq == sysfreq_i, 1);

    %if frequency is not available, use G01/G02 instead (as specified in the IGS ATX file)
    if (isempty(index_freq))
        sysfreq_i = 1*10 + frequency;
        index_freq = find(antenna_PCV.sysfreq == sysfreq_i, 1);
    end

    if ~isempty(index_freq) % corrections are available

        %index of input satellites belonging to constellation sys_i
        index_sat=find(sys==sys_i);
        XS_i = XS(index_sat,:);

        %correction along the receiver-satellite line-of-sight
        XRcorr = local2globalPos(antenna_PCV.offset(1,:,index_freq)', XR);
        corrXYZ = XRcorr - XR;
        PCO_corr_i = zeros(size(index_sat));
        for s = 1 : size(XS_i,1)
            LOS  = XR - XS_i(s,:)';
            LOSu = LOS / norm(LOS);
            PCO_corr_i(s,1) = dot(corrXYZ,LOSu);
        end

        PCO_corr(index_sat) = PCO_corr(index_sat) + PCO_corr_i;
    else
        % corrections not available
    end
end
