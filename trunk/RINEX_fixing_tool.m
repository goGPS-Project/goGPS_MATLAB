function RINEX_fixing_tool(filename_obs, filename_nav)

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
%   RINEX files due to receiver clock shifts.

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

global_init;

%load RINEX files
[pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
    dop1_R, dop1_M, dop2_R, dop2_M, Eph_R, Eph_M, iono_R, iono_M, snr1_R, snr1_M, ...
    snr2_R, snr2_M, pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
    dop1_RR, dop1_MR, dop2_RR, dop2_MR, Eph_RR, Eph_MR, snr_RR, snr_MR, ...
    time_GPS, time_R, time_M, date, posR] = load_RINEX(filename_obs, filename_nav); %#ok<ASGLU>

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

    if (size(sat,1) >= 4)
        %estimate receiver clock error by Bancroft algorithm
        input_bancroft(pr1_R(sat,i), sat, time_R(i), Eph_t);
    end
    
    %store receiver clock error in an array for later use
    dtR(i) = rec_clock_error; %#ok<NODEF>
    
    for j = 1 : size(sat,1)
        %satellite position computation
        [posS(:,sat(j)), dtS(sat(j)), posS_ttime(:,sat(j)), velS(:,sat(j)), ttime(sat(j))] = sat_corr(Eph_t, sat(j), time_R(i), pr1_R(sat(j),i));
    end

    if (i > 1)
        %receiver clock drift
        if (dtR(i) ~= 0 && dtR(i-1) ~= 0)
            dtRdot(i-1) = (dtR(i) - dtR(i-1));
        end
    end
end

%check if there is any discontinuity in the clock drift
sigma_dtRdot = std(dtRdot);
disc = find(abs(dtRdot) > 3*sigma_dtRdot);

%remove the discontinuities in order to compute a proper moving average later on
for i = 1 : length(disc)
    dtRdot(disc(i)) = median([dtRdot(disc(i)-5) dtRdot(disc(i)+5)]);
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
    rec_clock_error = dtR(i); %#ok<NASGU>
    
    for j = 1 : size(sat,1)
        %satellite position computation
        [posS(:,sat(j)), dtS(sat(j)), posS_ttime(:,sat(j)), velS(:,sat(j)), ttime(sat(j))] = sat_corr(Eph_t, sat(j), time_R(i), pr1_R(sat(j),i));
    end
    
    %get an estimation of the receiver clock drifting by moving average
    if (i < 50)
        avg_dtRdot = mean(dtRdot(1:50));
    elseif (i < nEpochs-51)
        avg_dtRdot = mean(dtRdot(i-49:i+50));
    elseif (i > nEpochs-51)
        avg_dtRdot = mean(dtRdot(nEpochs-51:nEpochs-1));
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

if ~isempty(disc)
    for i = 1 : length(disc)
        %check the presence of (and in case correct) any code "ambiguity" shift
        for s = 1 : 32
            if (pr1_R(s,disc(i)+1) & pr1_R(s,disc(i)))
                if (abs(pr1_R(s,disc(i)+1)-pr1_R(s,disc(i))) > 2e5)
                    diff = v_light * 0.001;
                    if (pr1_R(s,disc(i)+1)<pr1_R(s,disc(i)))
                        diff = -diff;
                    end
                    n = round(abs(pr1_R(s,disc(i)+1)-pr1_R(s,disc(i)))/diff);
                    for j = disc(i)+1 : nEpochs
                        pos = find(pr1_R(:,j) ~= 0);
                        pr1_R(pos,j) = pr1_R(pos,j) - diff*n;
                        pos = find(pr2_R(:,j) ~= 0);
                        pr2_R(pos,j) = pr2_R(pos,j) - diff*n;
                    end
                end
                break
            end
        end

        %correct code and phase for clock shift effect
        for s = 1 : 32
            if (ph1_R(s,disc(i)+1) & ph1_R(s,disc(i)))
                diff1 = ph1_R(s,disc(i)+1) - (ph1_R(s,disc(i)) - doppler_app1(s,disc(i)));
                diff2 = ph2_R(s,disc(i)+1) - (ph2_R(s,disc(i)) - doppler_app2(s,disc(i)));
                for j = disc(i)+1 : nEpochs
                    if (pr1_R(s,j) ~= 0)
                        pr1_R(s,j) = pr1_R(s,j) - diff1*lambda1;
                    end
                    if (pr2_R(s,j) ~= 0)
                        pr2_R(s,j) = pr2_R(s,j) - diff2*lambda2;
                    end
                    if (ph1_R(s,j) ~= 0)
                        ph1_R(s,j) = ph1_R(s,j) - diff1;
                    end
                    if (ph2_R(s,j) ~= 0)
                        ph2_R(s,j) = ph2_R(s,j) - diff2;
                    end
                end
            end
        end
    end
end

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
            fprintf(fid_obs,'%14.3f %1d',pr1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
            fprintf(fid_obs,'%14.3f %1d',pr2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
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
            fprintf(fid_obs,'%14.3f %1d',snr1_R(sat(j),i),floor(snr1_R(sat(j),i)/6));
            fprintf(fid_obs,'\n');
            fprintf(fid_obs,'%14.3f %1d',snr2_R(sat(j),i),floor(snr2_R(sat(j),i)/6));
            fprintf(fid_obs,'\n');
        end
    end
end

%close RINEX observation file
fclose(fid_obs);