%function RINEX_fixing_tool(filename_obs, filename_nav)

% SYNTAX:
%   RINEX_fixing_tool(filename_obs, filename_nav);
%
% INPUT:
%   filename_obs = input RINEX observation file
%   filename_nav = input RINEX navigation file
%
% OUTPUT:
%
% DESCRIPTION:
%   Script that checks and in case fixes observation discontinuities in
%   RINEX files due to receiver clock shifts and cycle slips.

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

%filename_obs = '../../test/switzerland/static-manno-home/static-switzerland-manno-home_rover.obs';
%filename_nav = '../../test/switzerland/static-manno-home/static-switzerland-manno-home_rover.nav';
%filename_obs = '../../test/switzerland/static-manno-cryms/static-switzerland-manno-cryms_rover.obs';
%filename_obs = '../../test/switzerland/static-manno-cryms/VirA052L.11o';
%filename_nav = '../../test/switzerland/static-manno-cryms/static-switzerland-manno-cryms_rover.nav';
%filename_obs = '../../test/italy/goPilastro/raw/como_pillar_rover.obs';
filename_obs = '../../test/italy/goPilastro/raw/como_pillar_master.10o';
filename_nav = '../../test/italy/goPilastro/raw/como_pillar_rover.nav';

global_init;

fprintf('Reading RINEX files...\n');

%load RINEX files
[pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
    dop1_R, dop1_M, dop2_R, dop2_M, Eph_R, Eph_M, iono_R, iono_M, snr1_R, snr1_M, ...
    snr2_R, snr2_M, pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
    dop1_RR, dop1_MR, dop2_RR, dop2_MR, Eph_RR, Eph_MR, snr_RR, snr_MR, ...
    time_GPS, time_R, time_M, date, posR] = load_RINEX(filename_obs, filename_nav); %#ok<ASGLU>

%detect epochs without data (GPS)
delepochs = find(time_R == 0);
if ~isempty(delepochs)
    fprintf('%d epochs without data\n', length(delepochs));
end

%number of epochs
nEpochs = length(time_R);

%remove satellites without ephemerides (GPS)
delsat = setdiff(1:32,unique(Eph_R(1,:)));
pr1_R(delsat,:) = 0;
pr2_R(delsat,:) = 0;
ph1_R(delsat,:) = 0;
ph2_R(delsat,:) = 0;
snr1_R(delsat,:) = 0;
snr2_R(delsat,:) = 0;

%receiver clock error
dtR = zeros(nEpochs,1);

%receiver clock drift
dtRdot = zeros(nEpochs-1,1);

%Doppler shift
doppler_app1 = zeros(32,nEpochs);
doppler_app2 = zeros(32,nEpochs);

if (posR == 0)
    %find available satellites
    sat = find(pr1_R(:,1) ~= 0);

    %find corresponding ephemeris
    Eph_t = rt_find_eph (Eph_R, time_R(1));
    
    [posR] = input_bancroft(pr1_R(sat,1), sat, time_R(1), Eph_t);
end

%----------------------------------------------------------------------------------------------
% RECEIVER CLOCK ERROR AND DRIFT
%----------------------------------------------------------------------------------------------

fprintf('Computing receiver clock error and drift...\n');

for i = 1 : nEpochs

    %find available satellites
    sat = find(pr1_R(:,i) ~= 0);
    
    %find corresponding ephemeris
    Eph_t = rt_find_eph (Eph_R, time_R(i));
    
    posS = zeros(size(sat,1),3);
    
    for j = 1 : size(sat,1)
        %satellite position computation
        [posS(j,:)] = sat_corr(Eph_t, sat(j), time_R(i), pr1_R(sat(j),i));
    end
    
    if (size(sat,1) >= 5)

        %initialization
        azR = zeros(32,1);
        elR = zeros(32,1);
        distR = zeros(32,1);

        %satellite azimuth, elevation, ROVER-SATELLITE distance
        [azR(sat), elR(sat), distR(sat)] = topocent(posR, posS);

        %elevation cut-off
        sat_cutoff = find(elR > cutoff);
        sat = intersect(sat,sat_cutoff);

        %estimate receiver clock error
        code_SA(posR, pr1_R(sat,i), snr_R(sat,i), sat, time_R(i), Eph_t, iono_R);

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
disc = find((abs(dtRdot-mean(dtRdot))) > clock_thresh);

%remove deleted epochs from discontinuities, because otherwise
%Doppler-predicted observations are going to be used on the wrong epochs;
%if an epoch is missing the discontinuity cannot be fixed
for i = 1 : length(delepochs)
	disc = setdiff(disc,(delepochs(i)-i));
end

%display the clock error
figure
plot(dtR)
%display the clock drift
figure
plot(dtRdot)
hold on
plot([1, nEpochs], [clock_thresh clock_thresh],'r');

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
    sat = find(pr1_R(:,i) ~= 0);
    
    %find corresponding ephemeris
    Eph_t = rt_find_eph (Eph_R, time_R(i));
    
    %initialization
    posS = zeros(3,32);
    dtS  = zeros(32,1);
    posS_ttime = zeros(3,32);
    velS = zeros(3,32);
    ttime = zeros(32,1);
    
    %read the receiver clock error from the array
    rec_clock_error = dtR(i);
    
    for j = 1 : size(sat,1)
        %satellite position computation
        [posS(:,sat(j)), dtS(sat(j)), posS_ttime(:,sat(j)), velS(:,sat(j)), ttime(sat(j))] = sat_corr(Eph_t, sat(j), time_R(i), pr1_R(sat(j),i));
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
        [doppler_app1(sat(j),i), doppler_app2(sat(j),i)] = doppler_shift_approx(posR, zeros(3,1), posS_ttime(:,sat(j)), velS(:,sat(j)), ttime(sat(j)), avg_dtRdot, sat(j), Eph_t);
    end
end

% %DISPLAY FOR DEBUG
% for i = 2 : nEpochs
%     %find available satellites
%     sat = find(pr1_R(:,i) ~= 0);
%
%     %display GPS time
%     time_R(i)
%
%     %display difference between observed phase and Doppler-predicted phase
%     diff = ph1_R(sat,i) - (ph1_R(sat,i-1) - doppler_app1(sat,i-1))
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
    
    %correct pseudorange for receiver clock error
    for s = 1 : 32
        for i = 1 : nEpochs
            if (pr1_R(s,i) ~= 0)
                pr1_R(s,i) = pr1_R(s,i) - v_light*dtR(i);
            end
        end
    end
    
    for i = 1 : length(disc)
        
%         %check the presence of (and in case correct) any code "ambiguity" shift
%         % (useful to correct also those satellites with a signal outage during a clock discontinuity)
%         for s = 1 : 32
%             if (pr1_R(s,disc(i)+1) & pr1_R(s,disc(i)))
%                 if (abs(pr1_R(s,disc(i)+1)-pr1_R(s,disc(i))) > 2e5)
%                     diff = v_light * (dtR(disc(i)+1) - dtR(disc(i)));
%                     %if (pr1_R(s,disc(i)+1)<pr1_R(s,disc(i)))
%                         %diff = -diff;
%                     %end
%                     %n = round(abs(pr1_R(s,disc(i)+1)-pr1_R(s,disc(i)))/diff);
%                     for j = disc(i)+1 : nEpochs
%                         pos = find(pr1_R(:,j) ~= 0);
%                         pr1_R(pos,j) = pr1_R(pos,j) - diff;%*n;
%                         pos = find(pr2_R(:,j) ~= 0);
%                         pr2_R(pos,j) = pr2_R(pos,j) - diff;%*n;
%                     end
%                 end
%                 break
%             end
%         end

        %correct code and phase for clock shift effect
        for s = 1 : 32
            if (ph1_R(s,disc(i)+1) & ph1_R(s,disc(i)))
                if (dop1_R(s,disc(i)) ~= 0)
                    doppler1 = dop1_R(s,disc(i));
                else
                    doppler1 = doppler_app1(s,disc(i));
                end
                if (dop2_R(s,disc(i)) ~= 0)
                    doppler2 = dop2_R(s,disc(i));
                else
                    doppler2 = doppler_app2(s,disc(i));
                end
                %diff1_pr = pr1_R(s,disc(i)+1) - (pr1_R(s,disc(i)) - doppler1*lambda1);
                %diff2_pr = pr2_R(s,disc(i)+1) - (pr2_R(s,disc(i)) - doppler2*lambda2);
                diff1_ph = ph1_R(s,disc(i)+1) - (ph1_R(s,disc(i)) - doppler1);
                diff2_ph = ph2_R(s,disc(i)+1) - (ph2_R(s,disc(i)) - doppler2);
                for j = disc(i)+1 : nEpochs
%                     if (pr1_R(s,j) ~= 0)
%                         pr1_R(s,j) = pr1_R(s,j) - diff1_pr;
%                     end
%                     if (pr2_R(s,j) ~= 0)
%                         pr2_R(s,j) = pr2_R(s,j) - diff2_pr;
%                     end
                    if (ph1_R(s,j) ~= 0)
                        ph1_R(s,j) = ph1_R(s,j) - diff1_ph;
                    end
                    if (ph2_R(s,j) ~= 0)
                        ph2_R(s,j) = ph2_R(s,j) - diff2_ph;
                    end
                end
            end
        end
    end
end

%----------------------------------------------------------------------------------------------
% CYCLE SLIPS (DETECTION AND FIXING)
%----------------------------------------------------------------------------------------------

fprintf('Detecting cycle slips... ');

%cycle-slip threshold
cs_thresh = 3;

flag_cs = 0;

for s = 1 : 32
    
    % L1
    index = find(ph1_R(s,:) ~= 0)';
    if ~isempty(index)
        ph = ph1_R(s,:);
        if (sum(dop1_R(s,:)) ~= 0)
            dp = dop1_R(s,:);
        else
            dp = doppler_app1(s,:);
        end
        diff_ph = zeros(nEpochs-1,1);
        for j = 2 : nEpochs
            if (ph(j) ~= 0 & ph(j-1) ~= 0 & dp(j-1) ~= 0)
                diff_ph(j) = ph(j) - (ph(j-1) - dp(j-1));
            end
        end

        %check if there is any cycle slip for satellite s
        cs = find(abs(diff_ph) > cs_thresh);
        
        figure
        plot(abs(diff_ph))
        hold on
        plot([1, nEpochs-1], [cs_thresh cs_thresh],'r');
        
        if ~isempty(cs)
            
            if (flag_cs == 0)
                fprintf('\n');
            end
            fprintf('SAT %d: fixing %d cycle slips on L1...\n', s, length(cs));
 
            for j = 1 : length(cs)
                for k = cs(j) : nEpochs
                    if (ph1_R(s,k) ~= 0)
                        ph1_R(s,k) = ph1_R(s,k) - diff_ph(cs(j));
                    end
                end
            end
            
            flag_cs = 1;
        end
    end
    
    % L2
    index = find(ph2_R(s,:) ~= 0)';
    if ~isempty(index)
        ph = ph2_R(s,:);
        if (sum(dop2_R(s,:)) ~= 0)
            dp = dop2_R(s,:);
        else
            dp = doppler_app2(s,:);
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
            fprintf('SAT %d: fixing %d cycle slips on L1...\n', s, length(cs));
            
            for j = 1 : length(cs)
                for k = cs(j) : nEpochs
                    if (ph2_R(s,k) ~= 0)
                        ph2_R(s,k) = ph2_R(s,k) - diff_ph(cs(j));
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

if (~isempty(disc) | flag_cs)

    %----------------------------------------------------------------------------------------------
    % WRITE FIXED RINEX OBSERVATION FILE
    %----------------------------------------------------------------------------------------------
    
    %displaying
    fprintf(['Writing: ' filename_obs '.fixed\n']);
    
    %create RINEX observation file
    fid_obs = fopen([filename_obs '.fixed'],'wt');
    
    %write header
    fprintf(fid_obs,'     2.10           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE\n');
    fprintf(fid_obs,'goGPS                                                       PGM / RUN BY / DATE \n');
    fprintf(fid_obs,'                                                            MARKER NAME         \n');
    fprintf(fid_obs,'                                                            OBSERVER / AGENCY   \n');
    fprintf(fid_obs,'                                                            REC # / TYPE / VERS \n');
    fprintf(fid_obs,'                                                            ANT # / TYPE        \n');
    fprintf(fid_obs,'%14.4f%14.4f%14.4f                  APPROX POSITION XYZ \n', posR(1), posR(2), posR(3));
    fprintf(fid_obs,'        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N\n');
    fprintf(fid_obs,'     1     1     0                                          WAVELENGTH FACT L1/2\n');
    fprintf(fid_obs,'     6    C1    P2    L1    L2    S1    S2                  # / TYPES OF OBSERV \n');
    fprintf(fid_obs,'     1.000                                                  INTERVAL            \n');
    fprintf(fid_obs,'%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS   \n', ...
        date(1,1), date(1,2), date(1,3), date(1,4), date(1,5), date(1,6));
    fprintf(fid_obs,'                                                            END OF HEADER       \n');
    
    %-------------------------------------------------------------------------------
    
    if (nargin == 3)
        waitbar(0,wait_dlg,'Writing rover observation file...')
    end
    
    %write data
    for i = 1 : nEpochs
        if (nargin == 3)
            waitbar(i/nEpochs,wait_dlg)
        end
        
        sat = find(pr1_R(:,i) ~= 0);
        n = length(sat);
        
        %if no observations are available, do not write anything
        if (n > 0)
            fprintf(fid_obs,' %02d %2d %2d %2d %2d %10.7f  0 %2d', ...
                date(i,1), date(i,2), date(i,3), date(i,4), date(i,5), date(i,6), n);
            if (n>12)
                for j = 1 : 12
                    fprintf(fid_obs,'G%02d',sat(j));
                end
                fprintf(fid_obs,'\n');
                fprintf(fid_obs,'%32s','');
                for j = 13 : n
                    fprintf(fid_obs,'G%02d',sat(j));
                end
            else
                for j = 1 : n
                    fprintf(fid_obs,'G%02d',sat(j));
                end
            end
            fprintf(fid_obs,'\n');
            for j = 1 : n
                if (abs(pr1_R(sat(j),i)) > 1e-100)
                    fprintf(fid_obs,'%14.3f %1d',pr1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
                else
                    fprintf(fid_obs,'                ');
                end
                if (abs(pr2_R(sat(j),i)) > 1e-100)
                    fprintf(fid_obs,'%14.3f %1d',pr2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
                else
                    fprintf(fid_obs,'                ');
                end
                if (abs(ph1_R(sat(j),i)) > 1e-100)
                    fprintf(fid_obs,'%14.3f %1d',ph1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
                else
                    fprintf(fid_obs,'                ');
                end
                if (abs(ph2_R(sat(j),i)) > 1e-100)
                    fprintf(fid_obs,'%14.3f %1d',ph2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
                else
                    fprintf(fid_obs,'                ');
                end
                if (abs(snr1_R(sat(j),i)) > 1e-100)
                    fprintf(fid_obs,'%14.3f %1d',snr1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
                else
                    fprintf(fid_obs,'                ');
                end
                fprintf(fid_obs,'\n');
                if (abs(snr2_R(sat(j),i)) > 1e-100)
                    fprintf(fid_obs,'%14.3f %1d',snr2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
                else
                    fprintf(fid_obs,'                ');
                end
                fprintf(fid_obs,'\n');
            end
        end
    end
    
    %close RINEX observation file
    fclose(fid_obs);
else
    if (isempty(delepochs))
        fprintf('No fixing required.\n');
    else
        fprintf('No fixing required (missing epochs are not fixed).\n');
    end
end
