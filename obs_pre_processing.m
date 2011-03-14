function [pr1, ph1, pr2, ph2, dop1, dop2, dtR, dtRdot] = obs_pre_processing(pr1, ph1, pr2, ...
          ph2, dop1, dop2, Eph, iono, snr1, time, posR)

% SYNTAX:
%   [pr1, ph1, pr2, ph2, dop1, dop2, dtR, dtRdot] = obs_pre_processing(pr1, ph1, pr2, ...
%    ph2, dop1, dop2, Eph, iono, snr1, time, posR);
%
% INPUT:
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   dop1 = Doppler observation (L1 carrier)
%   dop2 = Doppler observation (L2 carrier)
%   Eph = matrix containing 29 ephemerides for each satellite
%   iono = matrix containing ionosphere parameters
%   snr1 = signal-to-noise ratio
%   time = GPS time
%   posR = receiver (approximate) position
%
% OUTPUT:
%   pr1 = processed code observation (L1 carrier)
%   ph1 = processed phase observation (L1 carrier)
%   pr2 = processed code observation (L2 carrier)
%   ph2 = processed phase observation (L2 carrier)
%   dop1 = Doppler observation / integrated shift (L1 carrier)
%   dop2 = Doppler observation / integrated shift (L2 carrier)
%   dtR = receiver clock error
%   dtRdot receiver clock drift
%
% DESCRIPTION:
%   Script that checks and in case fixes observation discontinuities due to
%   receiver clock adjustments and cycle slips. Estimated receiver clock
%   error and drift for each epoch can also be returned.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
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

global lambda1 lambda2 v_light
global cutoff rec_clock_error

%detect epochs without data (GPS)
delepochs = find(time == 0);
if ~isempty(delepochs)
    fprintf('%d epochs without data\n', length(delepochs));
end

%number of epochs
nEpochs = length(time);

%receiver clock error
dtR = zeros(nEpochs,1);

%receiver clock drift
dtRdot = zeros(nEpochs-1,1);

%integrated Doppler shift
doppler_int1 = zeros(32,nEpochs);
doppler_int2 = zeros(32,nEpochs);

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

fprintf('Computing receiver clock error and drift...\n');

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

        %store receiver clock error in an array for later use
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
% RECEIVER CLOCK DISCONTINUITIES (DETECTION)
%----------------------------------------------------------------------------------------------

fprintf('Detecting clock drift discontinuities... ');

%check if there is any discontinuity in the clock drift
clock_thresh = 1e-4;
disc = find(abs(dtRdot-mean(dtRdot)) > clock_thresh);

%remove deleted epochs from discontinuities, because otherwise
%Doppler-predicted observations are going to be used on the wrong epochs;
%if an epoch is missing the discontinuity cannot be fixed
for i = 1 : length(delepochs)
	disc = setdiff(disc,(delepochs(i)-i));
end

% %display the clock error
% figure
% plot(dtR)
% %display the absolute value of the clock drift
% figure
% plot(abs(dtRdot))
% hold on
% plot([1, nEpochs], [clock_thresh clock_thresh],'r');

fprintf('%d detected.\n', length(disc));

%----------------------------------------------------------------------------------------------
% DOPPLER SHIFT
%----------------------------------------------------------------------------------------------

fprintf('Computing Doppler shift...\n');

%remove the discontinuities in order to compute a proper moving average later on
for i = 1 : length(disc)
    if (disc(i) < 5)
        dtRdot(disc(i)) = median([dtRdot(1) dtRdot(10)]);
    elseif (disc(i) <= nEpochs-6)
        dtRdot(disc(i)) = median([dtRdot(disc(i)-4) dtRdot(disc(i)+5)]);
    elseif (disc(i) > nEpochs-6)
        dtRdot(disc(i)) = median([dtRdot(nEpochs-10) dtRdot(nEpochs-1)]);
    end
end

for i = 1 : nEpochs
    
    %find available satellites
    sat = find(pr1(:,i) ~= 0);
    
    %find corresponding ephemeris
    Eph_t = rt_find_eph (Eph, time(i));
    
    %initialization
    posS = zeros(3,32);
    dtS  = zeros(32,1);
    posS_ttime = zeros(3,32);
    velS = zeros(3,32);
    ttime = zeros(32,1);
    
    for j = 1 : size(sat,1)
        %satellite position computation
        [posS(:,sat(j)), dtS(sat(j)), posS_ttime(:,sat(j)), velS(:,sat(j)), ttime(sat(j))] = sat_corr(Eph_t, sat(j), time(i), pr1(sat(j),i));
    end
    
    %get an estimation of the receiver clock drift by moving average
    if (i < 50)
        avg_dtRdot = mean(dtRdot(1:100));
    elseif (i <= nEpochs-51)
        avg_dtRdot = mean(dtRdot(i-49:i+50));
    elseif (i > nEpochs-51)
        avg_dtRdot = mean(dtRdot(nEpochs-100:nEpochs-1));
    end
    
    %Doppler shift computation using the estimated receiver clock drift
    for j = 1 : size(sat,1)
        [doppler_int1(sat(j),i), doppler_int2(sat(j),i)] = doppler_shift_approx(posR(:,i), zeros(3,1), posS_ttime(:,sat(j)), velS(:,sat(j)), ttime(sat(j)), avg_dtRdot, sat(j), Eph_t);
    end
end

% %DISPLAY FOR DEBUG
% for i = 2 : nEpochs
%     %find available satellites
%     sat = find(pr1(:,i) ~= 0);
%
%     %display GPS time
%     time(i)
%
%     %display difference between observed phase and Doppler-predicted phase
%     diff = ph1(sat,i) - (ph1(sat,i-1) - doppler_int1(sat,i-1))
%
%     pos = find(abs(diff) > 1);
%     if (~isempty(pos))
%         pause
%     end
% end

%----------------------------------------------------------------------------------------------
% RECEIVER CLOCK DISCONTINUITIES (FIXING)
%----------------------------------------------------------------------------------------------

if ~isempty(disc)

    if (isempty(delepochs))
        fprintf('Fixing %d clock discontinuities...\n', length(disc));
    else
        fprintf('Fixing %d clock discontinuities (missing epochs will not be fixed)...\n', length(disc));
    end

    for i = 1 : length(disc)
        
        %correct phase for clock shift effect
        for s = 1 : 32
            if (ph1(s,disc(i)+1) & ph1(s,disc(i)))
                if (dop1(s,disc(i)) ~= 0)
                    doppler1 = dop1(s,disc(i));
                else
                    doppler1 = doppler_int1(s,disc(i));
                end
                if (dop2(s,disc(i)) ~= 0)
                    doppler2 = dop2(s,disc(i));
                else
                    doppler2 = doppler_int2(s,disc(i));
                end
                diff1_ph = ph1(s,disc(i)+1) - (ph1(s,disc(i)) - doppler1);
                diff2_ph = ph2(s,disc(i)+1) - (ph2(s,disc(i)) - doppler2);
                for j = disc(i)+1 : nEpochs
                    if (ph1(s,j) ~= 0)
                        ph1(s,j) = ph1(s,j) - diff1_ph;
                    end
                    if (ph2(s,j) ~= 0)
                        ph2(s,j) = ph2(s,j) - diff2_ph;
                    end
                end
            end
        end
    end
end

%----------------------------------------------------------------------------------------------
% OBSERVATION CORRECTION
%----------------------------------------------------------------------------------------------

%correct pseudorange for receiver clock error and compute phase range
%from accumulated phase if necessary
for s = 1 : 32
    for i = 1 : nEpochs
        if (pr1(s,i) ~= 0)
            pr1(s,i) = pr1(s,i) - v_light*dtR(i);
            if (ph1(s,i) ~= 0 & abs(ph1(s,i)) < 6e7)
                ambig = 2^23;
                n = floor((pr1(s,i)/lambda1-ph1(s,i)) / ambig + 0.5 );
                ph1(s,i) = ph1(s,i) + n*ambig;
            end
        end
        if (pr2(s,i) ~= 0)
            pr2(s,i) = pr2(s,i) - v_light*dtR(i);
            if (ph2(s,i) ~= 0 & abs(ph2(s,i)) < 6e7)
                ambig = 2^23;
                n = floor((pr2(s,i)/lambda2-ph2(s,i)) / ambig + 0.5 );
                ph2(s,i) = ph2(s,i) + n*ambig;
            end
        end
    end
end

%----------------------------------------------------------------------------------------------
% CYCLE SLIPS (DETECTION AND FIXING)
%----------------------------------------------------------------------------------------------

fprintf('Detecting cycle slips... ');

%cycle-slip threshold
cs_thresh = 10;

flag_cs = 0;

for s = 1 : 32
    
    % L1
    index = find(ph1(s,:) ~= 0)';
    if ~isempty(index)
        ph = ph1(s,:);
        if (sum(dop1(s,:)) ~= 0)
            dp = dop1(s,:);
        else
            dp = doppler_int1(s,:);
        end
        diff_ph = zeros(nEpochs-1,1);
        for j = 2 : nEpochs
            if (ph(j) ~= 0 & ph(j-1) ~= 0 & dp(j-1) ~= 0)
                diff_ph(j) = ph(j) - (ph(j-1) - dp(j-1));
            end
        end

        %check if there is any cycle slip for satellite s
        cs = find(abs(diff_ph) > cs_thresh);
        
%         figure
%         plot(abs(diff_ph))
%         hold on
%         plot([1, nEpochs-1], [cs_thresh cs_thresh],'r');
        
        if ~isempty(cs)
            
            if (flag_cs == 0)
                fprintf('\n');
            end
            fprintf('SAT %d: fixing %d cycle slips on L1...\n', s, length(cs));
 
            for j = 1 : length(cs)
                for k = cs(j) : nEpochs
                    if (ph1(s,k) ~= 0)
                        ph1(s,k) = ph1(s,k) - diff_ph(cs(j));
                    end
                end
            end
            
            flag_cs = 1;
        end
    end
    
    % L2
    index = find(ph2(s,:) ~= 0)';
    if ~isempty(index)
        ph = ph2(s,:);
        if (sum(dop2(s,:)) ~= 0)
            dp = dop2(s,:);
        else
            dp = doppler_int2(s,:);
        end
        diff_ph = zeros(nEpochs-1,1);
        for j = 2 : nEpochs
            if (ph(j) ~= 0 & ph(j-1) ~= 0 & dp(j-1) ~= 0)
                diff_ph(j) = ph(j) - (ph(j-1) - dp(j-1));
            end
        end

        %check if there is any cycle slip for satellite s
        cs = find(abs(diff_ph) > cs_thresh);
        
        if ~isempty(cs)
            
            if (flag_cs == 0)
                fprintf('\n');
            end
            fprintf('SAT %d: fixing %d cycle slips on L2...\n', s, length(cs));
            
            for j = 1 : length(cs)
                for k = cs(j) : nEpochs
                    if (ph2(s,k) ~= 0)
                        ph2(s,k) = ph2(s,k) - diff_ph(cs(j));
                    end
                end
            end
            
            flag_cs = 1;
        end
    end
end
if (flag_cs == 0)
    fprintf('none detected.\n');
end

% select Doppler shift to be returned
if (dop1 == 0)
    dop1 = doppler_int1;
end
if (dop2 == 0)
    dop2 = doppler_int2;
end
