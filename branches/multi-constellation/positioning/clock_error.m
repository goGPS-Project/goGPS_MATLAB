function [dtR, dtRdot] = clock_error(XR0, time_rx, pr, snr, Eph, SP3, iono, sbas)

% SYNTAX:
%   [dtR, dtRdot] = clock_error(XR0, time_rx, pr, snr, Eph, SP3, iono, sbas);
%
% INPUT:
%   XR0      = receiver (approximate) position
%   time_rx  = GPS time
%   pr       = code observation (L1 carrier)
%   snr      = signal-to-noise ratio
%   Eph      = matrix containing 30 navigation parameters for each satellite
%   SP3      = structure containing precise ephemeris data
%   iono     = ionosphere parameters (Klobuchar)
%   sbas     = SBAS corrections
%
% OUTPUT:
%   dtR = receiver clock error
%   dtRdot receiver clock drift
%
% DESCRIPTION:
%   Compute receiver clock error and drift for each epoch.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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

global cutoff snr_threshold

%number of epochs
nEpochs = length(time_rx);

%receiver clock error
dtR = zeros(nEpochs,1);

%receiver clock drift
dtRdot = zeros(nEpochs-1,1);

%------------------------------------------------------------------------------------
% APPROXIMATE POSITION
%-----------------------------------------------------------------------------------

if ((sum(abs(XR0)) == 0) | isempty(XR0))
    %approximate position not available
    flag_XR = 0;
else
    %approximate position available
    flag_XR = 1;
end

for i = 1 : nEpochs
    
    %--------------------------------------------------------------------------------------------
    % SATELLITE AND EPHEMERIS SELECTION
    %--------------------------------------------------------------------------------------------
    
    sat_pr = find(pr(:,i) ~= 0);
    nsat = length(sat_pr);

    Eph_t = rt_find_eph (Eph, time_rx(i), nsat);
    sbas_t = find_sbas(sbas, i);
    
    %----------------------------------------------------------------------------------------------
    % RECEIVER POSITION AND CLOCK ERROR
    %----------------------------------------------------------------------------------------------

    if (length(sat_pr) >= 4)
        [XR, dtR(i)] = init_positioning(time_rx(i), pr(sat_pr,i), snr(sat_pr,i), Eph_t, SP3, iono, sbas_t, XR0(:,i), [], [], sat_pr, cutoff, snr_threshold, flag_XR, 0); %#ok<ASGLU>
        
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
