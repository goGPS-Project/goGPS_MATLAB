function [dtR, dtRdot] = clock_error(pr1, Eph, iono, snr1, time, posR)

% SYNTAX:
%   [dtR, dtRdot] = clock_error(pr1, Eph, iono, snr1, time, posR);
%
% INPUT:
%   pr1 = code observation (L1 carrier)
%   Eph = matrix containing 29 ephemerides for each satellite
%   iono = matrix containing ionosphere parameters
%   snr1 = signal-to-noise ratio
%   time = GPS time
%   posR = receiver (approximate) position
%
% OUTPUT:
%   dtR = receiver clock error
%   dtRdot receiver clock drift
%
% DESCRIPTION:
%   Compute receiver clock error and drift for each epoch.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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

global cutoff rec_clock_error

%number of epochs
nEpochs = length(time);

%receiver clock error
dtR = zeros(nEpochs,1);

%receiver clock drift
dtRdot = zeros(nEpochs-1,1);

if (posR == 0)
    %find available satellites
    sat = find(pr1(:,1) ~= 0);

    %find corresponding ephemeris
    Eph_t = rt_find_eph (Eph, time(1));
    
    posR = input_bancroft(pr1(sat,1), sat, time(1), Eph_t);

    posR(1,1:nEpochs) = posR(1);
    posR(2,1:nEpochs) = posR(2);
    posR(3,1:nEpochs) = posR(3);
end

%----------------------------------------------------------------------------------------------
% RECEIVER CLOCK ERROR AND DRIFT
%----------------------------------------------------------------------------------------------

for i = 1 : nEpochs

    %find available satellites
    sat = find(pr1(:,i) ~= 0);
    
    %find corresponding ephemeris
    Eph_t = rt_find_eph (Eph, time(i));
    
    posS = zeros(size(sat,1),3);
    
    for j = 1 : size(sat,1)
        %satellite position computation
        [posS(j,:)] = sat_corr(Eph_t, sat(j), time(i), pr1(sat(j),i));
    end
    
    if (size(sat,1) >= 5)

        %initialization
        azR = zeros(32,1);
        elR = zeros(32,1);
        distR = zeros(32,1);

        %satellite azimuth, elevation, ROVER-SATELLITE distance
        [azR(sat), elR(sat), distR(sat)] = topocent(posR(:,i), posS); %#ok<NASGU,ASGLU>

        %elevation cut-off
        sat_cutoff = find(elR > cutoff);
        sat = intersect(sat,sat_cutoff);

        %estimate receiver clock error
        code_SA(posR(:,i), pr1(sat,i), snr1(sat,i), sat, time(i), Eph_t, iono);

        %store receiver clock error in an array
        dtR(i) = rec_clock_error;

        if (i > 1)
            %receiver clock drift
            if (dtR(i) ~= 0 && dtR(i-1) ~= 0)
                dtRdot(i-1) = (dtR(i) - dtR(i-1));
            end
        end
    else
        if (i > 2)
            dtR(i) = dtR(i-1) + (dtR(i-1) - dtR(i-2));
            dtRdot(i-1) = (dtR(i) - dtR(i-1));
        end
    end
end

%----------------------------------------------------------------------------------------------
% RECEIVER CLOCK DISCONTINUITIES
%----------------------------------------------------------------------------------------------

%check if there is any discontinuity in the clock drift
clock_thresh = 1e-4;
disc = find(abs(dtRdot-mean(dtRdot)) > clock_thresh);

%remove discontinuities from the clock drift
for i = 1 : length(disc)
    if (disc(i) < 5)
        dtRdot(disc(i)) = median([dtRdot(1) dtRdot(10)]);
    elseif (disc(i) <= nEpochs-6)
        dtRdot(disc(i)) = median([dtRdot(disc(i)-4) dtRdot(disc(i)+5)]);
    elseif (disc(i) > nEpochs-6)
        dtRdot(disc(i)) = median([dtRdot(nEpochs-10) dtRdot(nEpochs-1)]);
    end
end

dtRdot(end+1) = dtRdot(end);