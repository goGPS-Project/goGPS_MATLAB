function PCV_corr = PCV_correction(antenna_PCV, zen, azi, sys, frequency)

% SYNTAX:
%   PCV_corr = PCV_correction(antenna_PCV, zen, azi, sys, frequency);
%
% INPUT:
%   antenna_PCV = antenna PCV struct
%   zen = zenithal angle to be interpolated (deg, array)
%   azi = azimutal angle to be interpolated (deg, array)
%   sys = system identification code to match correct constellation (1: GPS, 2: GLONASS, ...)
%   frequency = frequency of the observations (i.e.: 1 for L1, 2 for L2, ...)
%
% OUTPUT:
%   PCV_corr = interpolated value of the PCV correction
%
% DESCRIPTION:
%   Function that computes the PCV corrections for the array of input
%   az/el of satellite with system in sys and for the specified frequency

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

PCV_corr=zeros(size(zen));

sys = floor(sys);

% extract constellations that have to be analyzed
constellation=unique(sys);

% compute PCV correction: loop on every required constellation
for j=1:length(constellation)
    % current system
    sys_i=constellation(j);

    % verify if corrections for this sys+frequency are available
    sysfreq_i=sys_i*10+frequency;
    index_freq=find(antenna_PCV.sysfreq==sysfreq_i, 1);

    if ~isempty(index_freq)
        % corrections are available

        % index of input ele/az belonging to satellites of sys_i
        index_sat=find(sys==sys_i);
        zen_i=zen(index_sat);
        azi_i=azi(index_sat);

        % compute azi-zen depending corrections
        % verify if only 'noazi' is available

        if isnan(antenna_PCV.tablePCV(1,1,index_freq))
            % compute noazi corrections
            % only L1
            tableNOAZI=antenna_PCV.tableNOAZI(:,:,index_freq);
            table_zen=antenna_PCV.tablePCV_zen;

            % detection of the grid node nearest to the interpolation point
            [~, posX] = min(abs(repmat(table_zen,length(zen_i),1) - repmat(zen_i,1,length(table_zen))),[],2);
            posX(zen_i-table_zen(posX)'<0)= posX(zen_i-table_zen(posX)'<0)-1;   % get the lower node
            mX=zen_i-table_zen(posX)';

            % special case if posX corresponds to the last node of the table
            index_end=find(posX==length(table_zen));
            posX(index_end)=posX(index_end)-1;
            mX(index_end)=mX(index_end)+antenna_PCV.dzen;

            % interpolation
            PCV_corr_i=(tableNOAZI(posX+1)'-tableNOAZI(posX)')/antenna_PCV.dzen.*mX+tableNOAZI(posX)';

            % keep only the corrections included in the range from zen1 to zen2
            PCV_corr_i(zen_i > antenna_PCV.zen2) = 0;
        else
            % compute azi-zen depending corrections
            tablePCV=antenna_PCV.tablePCV(:,:,index_freq);
            table_zen=antenna_PCV.tablePCV_zen;
            table_azi=antenna_PCV.tablePCV_azi;

            % reverse the direction in azimut (so that azimuth goes high-low
            tablePCV=flipud(tablePCV);
            table_azi=fliplr(table_azi);

            %detection of the grid node nearest to the interpolation point
            [~, posX] = min(abs(repmat(table_zen,length(zen_i),1) - repmat(zen_i,1,length(table_zen))),[],2);
            [~, posY] = min(abs(repmat(table_azi,length(azi_i),1) - repmat(azi_i,1,length(table_azi))),[],2);
            posX(zen_i-table_zen(posX)'<0)= posX(zen_i-table_zen(posX)'<0)-1;   % get the lower node
            posY(azi_i-table_azi(posY)'<0)= posY(azi_i-table_azi(posY)'<0)+1;   % get the lower node
            mX=zen_i-table_zen(posX)';
            mY=azi_i-table_azi(posY)';

            % special case if posX or posY correspond to the last node of the table
            index_end=find(posX==length(table_zen));
            posX(index_end)=posX(index_end)-1;
            mX(index_end)=mX(index_end)+antenna_PCV.dzen;

            index_end=find(posY==1);
            posY(index_end)=posY(index_end)+1;
            mY(index_end)=mY(index_end)+antenna_PCV.dazi;

            % bilinear interpolation
            PCV_corr_i= tablePCV(sub2ind(size(tablePCV),posY,posX)) ./ (antenna_PCV.dzen*antenna_PCV.dazi) .* (antenna_PCV.dzen-mX) .* (antenna_PCV.dazi-mY) + ...
                tablePCV(sub2ind(size(tablePCV),posY,posX+1)) ./ (antenna_PCV.dzen*antenna_PCV.dazi) .* mX .* (antenna_PCV.dazi-mY) + ...
                tablePCV(sub2ind(size(tablePCV),posY-1,posX)) ./ (antenna_PCV.dzen*antenna_PCV.dazi) .* (antenna_PCV.dzen-mX) .* mY + ...
                tablePCV(sub2ind(size(tablePCV),posY-1,posX+1)) ./ (antenna_PCV.dzen*antenna_PCV.dazi) .* mX .* mY;
        end

        PCV_corr(index_sat)=PCV_corr(index_sat)+PCV_corr_i;
    else
        % corrections not available
    end
end
