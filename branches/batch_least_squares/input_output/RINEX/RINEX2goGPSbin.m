function RINEX2goGPSbin(filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav, filerootOUT, wait_dlg)

% SYNTAX:
%   RINEX2goGPSbin(filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav, filerootOUT, wait_dlg);
%
% INPUT:
%   filename_R_obs = rover observation RINEX file
%   filename_R_nav = rover navigation RINEX file
%   filename_M_obs = master observation RINEX file
%   filename_M_nav = master navigation RINEX file
%   filerootOUT = output file root
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%
% DESCRIPTION:
%   File conversion from rover and master RINEX files to goGPS binary
%   format (*_obs_* and *_eph_* files).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

if (isempty(dir(filename_R_nav)) && isempty(dir(filename_M_nav)))
    if (nargin >= 6)
        msgbox('Navigation data (ephemerides) are required to create goGPS binary data.');
    else
        fprintf('Navigation data (ephemerides) are required to create goGPS binary data.\n');
    end
    return
elseif (isempty(dir(filename_R_nav)))
    filename_nav = filename_M_nav;
elseif (isempty(dir(filename_M_nav)))
    filename_nav = filename_R_nav;
else
    filename_nav = filename_M_nav;
end

%load multi-constellation settings and initialize 'constellations' struct
global goIni;
GPS_flag = goIni.getData('Constellations','GPS');
GLO_flag = goIni.getData('Constellations','GLONASS');
GAL_flag = goIni.getData('Constellations','Galileo');
BDS_flag = goIni.getData('Constellations','BeiDou');
QZS_flag = goIni.getData('Constellations','QZSS');
SBS_flag = goIni.getData('Constellations','SBAS');
[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
nSatTot = constellations.nEnabledSat;
if (nSatTot == 0)
    fprintf('No constellations selected, setting default: GPS-only processing\n');
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
    nSatTot = constellations.nEnabledSat;
end

if (~isempty(dir(filename_R_obs)))
    
    %ROVER RINEX files reading
    if (nargin >= 6)

        [pr1_R, ~, ph1_R, ~, ~, ~, ~, ~, dop1_R, ~, ~, ~, snr1_R, ~, ~, ~,  ...
         ~, time_R, ~, week_R, ~, ~, ~, ~, ~, Eph_R, iono_R, interval_R] = ...
         load_RINEX(filename_nav, filename_R_obs, [], constellations, 0, wait_dlg);
    else
        [pr1_R, ~, ph1_R, ~, ~, ~, ~, ~, dop1_R, ~, ~, ~, snr1_R, ~, ~, ~,  ...
         ~, time_R, ~, week_R, ~, ~, ~, ~, ~, Eph_R, iono_R, interval_R] = ...
         load_RINEX(filename_nav, filename_R_obs, [], constellations, 0);
    end
    
    %TEMP
    snr_R = snr1_R;
else
    if (nargin >= 6)
        msgbox('Rover data are required to create goGPS binary data.');
    else
        fprintf('Rover data are required to create goGPS binary data.\n');
    end
    
    return
end

time_GPS = time_R;
epochs  = length(time_GPS);
num_sat = size(pr1_R,1);

%MASTER variable initialization
time_M = zeros(epochs,1);
pr1_M  = zeros(num_sat,epochs);
ph1_M  = zeros(num_sat,epochs);
dop1_M = zeros(num_sat,epochs);
snr_M  = zeros(num_sat,epochs);
pos_M  = zeros(3,epochs);
Eph_M  = zeros(33,num_sat,epochs);
iono   = zeros(8,epochs);

if (~isempty(dir(filename_M_obs)))
    
    %MASTER RINEX files reading
    if (nargin >= 6)
        
        [pr1_M, ~, ph1_M, ~, ~, ~, ~, ~, dop1_M, ~, ~, ~, snr1_M, ~, ~, ~,  ...
         ~, time_M, ~, ~, ~, ~, ~, pos_M, ~, Eph_M, iono_M, interval_M] = ...
         load_RINEX(filename_nav, filename_M_obs, [], constellations, 0, wait_dlg);
    else
        [pr1_M, ~, ph1_M, ~, ~, ~, ~, ~, dop1_M, ~, ~, ~, snr1_M, ~, ~, ~,  ...
         ~, time_M, ~, ~, ~, ~, ~, pos_M, ~, Eph_M, iono_M, interval_M] = ...
         load_RINEX(filename_nav, filename_M_obs, [], constellations, 0);
    end
    
    %TEMP
    snr_M = snr1_M;

    %--------------------------------------------------------------------------
    
    if (nargin >= 6)
        waitbar(0.33,wait_dlg,'Synchronizing data...')
    end

    %round time values for synchronizing rover and master epochs
    %interval_R = median(time_R(2:end) - time_R(1:end-1));
    roundtime_R = roundmod(time_R,interval_R);
    
    %interval_M = median(time_M(2:end) - time_M(1:end-1));
    roundtime_M = roundmod(time_M,interval_M);
    
    if ~isempty(time_R) && ~isempty(time_M)
        
        if (size(pos_M,2) == 1)
            pos_M(1,1:length(time_M)) = pos_M(1);
            pos_M(2,1:length(time_M)) = pos_M(2);
            pos_M(3,1:length(time_M)) = pos_M(3);
        end
        
        %head synchronization
        if (roundtime_R(1) < roundtime_M(1))
            pos = find(roundtime_R < roundtime_M(1));
            roundtime_R(pos) = [];                     %GPS time (rounded)
            time_R(pos)    = [];                       %GPS time
            week_R(pos)    = [];                       %GPS week
            pr1_R(:,pos)   = [];                       %code observations
            ph1_R(:,pos)   = [];                       %phase observations
            snr_R(:,pos)   = [];                       %signal-to-noise ratio
            dop1_R(:,pos)  = [];                       %doppler measurement
            iono(:,pos) = [];                          %ionosphere parameters
        end
        
        if (roundtime_M(1) < roundtime_R(1))
            pos = find(roundtime_M < roundtime_R(1));
            roundtime_M(pos) = [];                     %GPS time (rounded)
            time_M(pos)    = [];                       %GPS time
            pr1_M(:,pos)   = [];                       %code observations
            ph1_M(:,pos)   = [];                       %phase observations
            snr_M(:,pos)   = [];                       %signal-to-noise ratio
            pos_M(:,pos)   = [];                       %master station position
        end
        
        %tail synchronization
        if (roundtime_R(end) > roundtime_M(end))
            pos = find(roundtime_R > roundtime_M(end));
            roundtime_R(pos) = [];                     %GPS time (rounded)
            time_R(pos)    = [];                       %GPS time
            week_R(pos)    = [];                       %GPS week
            pr1_R(:,pos)   = [];                       %code observations
            ph1_R(:,pos)   = [];                       %phase observations
            snr_R(:,pos)   = [];                       %signal-to-noise ratio
            dop1_R(:,pos)  = [];                       %doppler measurement
            iono(:,pos) = [];                          %ionosphere parameters
        end
        
        if (roundtime_M(end) > roundtime_R(end))
            pos = find(roundtime_M > roundtime_R(end));
            roundtime_M(pos) = [];                     %GPS time (rounded)
            time_M(pos)    = [];                       %GPS time
            pr1_M(:,pos)   = [];                       %code observations
            ph1_M(:,pos)   = [];                       %phase observations
            snr_M(:,pos)   = [];                       %signal-to-noise ratio
            pos_M(:,pos)   = [];                       %master station position
        end
        
    end
    
    %-------------------------------------------------------------------------------
    
    if (nargin >= 6)
        waitbar(0.66,wait_dlg)
    end
    
    %signal losses
    time_GPS = union(roundtime_R,roundtime_M);           %overall reference time
    
    if ~isempty(time_GPS)
        
        interval = median(time_GPS(2:end) - time_GPS(1:end-1));
        
        time_GPS = (time_GPS(1) : interval : time_GPS(end))';   %GPS time without interruptions
        
        if ~isempty(time_R)
            
            newtime_R = setdiff(time_GPS, roundtime_R);  %ROVER missing epochs
            for i = 1 : length(newtime_R)
                
                pos = find(roundtime_R == newtime_R(i) - interval);  %position before the "holes"
                
                time_R = [time_R(1:pos);  newtime_R(i);  time_R(pos+1:end)];
                week_R = [week_R(1:pos);  week_R(pos);   week_R(pos+1:end)]; %does not take into account week change (TBD)
                pr1_R  = [pr1_R(:,1:pos)  zeros(num_sat,1)    pr1_R(:,pos+1:end)];
                ph1_R  = [ph1_R(:,1:pos)  zeros(num_sat,1)    ph1_R(:,pos+1:end)];
                dop1_R = [dop1_R(:,1:pos) zeros(num_sat,1)    dop1_R(:,pos+1:end)];
                snr_R  = [snr_R(:,1:pos)  zeros(num_sat,1)    snr_R(:,pos+1:end)];
                iono   = [iono(:,1:pos)   zeros(8,1)     iono(:,pos+1:end)];
                
                roundtime_R = roundmod(time_R,interval_R);
            end
        else
            time_R = time_GPS;
            week_R = zeros(1,length(time_GPS));
            pr1_R  = zeros(num_sat,length(time_GPS));
            ph1_R  = zeros(num_sat,length(time_GPS));
            dop1_R = zeros(num_sat,length(time_GPS));
            snr_R  = zeros(num_sat,length(time_GPS));
            iono   = zeros(8,length(time_GPS));
        end
        
        if ~isempty(time_M)
            
            newtime_M = setdiff(time_GPS, roundtime_M);  %MASTER missing epochs
            for i = 1 : length(newtime_M)
                
                pos = find(roundtime_M == newtime_M(i) - interval);  %position before the "holes"
                
                time_M = [time_M(1:pos);  newtime_M(i);  time_M(pos+1:end)];
                pr1_M  = [pr1_M(:,1:pos)  zeros(num_sat,1)    pr1_M(:,pos+1:end)];
                ph1_M  = [ph1_M(:,1:pos)  zeros(num_sat,1)    ph1_M(:,pos+1:end)];
                snr_M  = [snr_M(:,1:pos)  zeros(num_sat,1)    snr_M(:,pos+1:end)];
                pos_M  = [pos_M(:,1:pos)  zeros(3,1)     pos_M(:,pos+1:end)];
                
                roundtime_M = roundmod(time_M,interval_M);
            end
        else
            time_M = time_GPS;
            pr1_M  = zeros(num_sat,length(time_GPS));
            ph1_M  = zeros(num_sat,length(time_GPS));
            snr_M  = zeros(num_sat,length(time_GPS));
            pos_M  = zeros(3,length(time_GPS));
        end
    end
    
    if (nargin >= 6)
        waitbar(1,wait_dlg)
    end
end

epochs = length(time_GPS);

%--------------------------------------------------------------------------

%select ephemerides source
if (~Eph_M)
    Eph_tmp = Eph_R;
else
    Eph_tmp = Eph_M;
end

%reconstruct ephemerides complete set
if (nargin >= 6)
    waitbar(0,wait_dlg,'Reconstruct complete ephemerides set...')
end
Eph = zeros(33,num_sat,epochs);
for i = 1 : epochs
    Eph(:,:,i) = rt_find_eph(Eph_tmp, time_GPS(i), nSatTot);
    if (nargin >= 6)
        waitbar(i/epochs,wait_dlg)
    end
end

%select ionosphere parameters source
if (isempty(dir(filename_M_obs)))
    iono(1,1:epochs) = iono_R(1);
    iono(2,1:epochs) = iono_R(2);
    iono(3,1:epochs) = iono_R(3);
    iono(4,1:epochs) = iono_R(4);
    iono(5,1:epochs) = iono_R(5);
    iono(6,1:epochs) = iono_R(6);
    iono(7,1:epochs) = iono_R(7);
    iono(8,1:epochs) = iono_R(8);
else
    iono(1,1:epochs) = iono_M(1);
    iono(2,1:epochs) = iono_M(2);
    iono(3,1:epochs) = iono_M(3);
    iono(4,1:epochs) = iono_M(4);
    iono(5,1:epochs) = iono_M(5);
    iono(6,1:epochs) = iono_M(6);
    iono(7,1:epochs) = iono_M(7);
    iono(8,1:epochs) = iono_M(8);
end

%master position
pos_M(1,1:epochs) = pos_M(1);
pos_M(2,1:epochs) = pos_M(2);
pos_M(3,1:epochs) = pos_M(3);

%remove epochs in which not all ephemerides are available or available satellites are < 4
satEph = find(sum(abs(Eph(:,:,1))) ~= 0);
if (~isempty(dir(filename_M_obs)))
    satObs = find( (pr1_R(:,1) ~= 0) & (pr1_M(:,1) ~= 0));
else
    satObs = find(pr1_R(:,1) ~= 0);
end

while (length(satEph) < length(satObs)) || (length(satObs) < 4)

    time_GPS(1) = [];
    week_R(1)   = [];
    time_R(1)   = [];
    time_M(1)   = [];
    pr1_R(:,1)  = [];
    pr1_M(:,1)  = [];
    ph1_R(:,1)  = [];
    ph1_M(:,1)  = [];
    dop1_R(:,1) = [];
    dop1_M(:,1) = [];
    snr_R(:,1)  = [];
    snr_M(:,1)  = [];
    pos_M(:,1)  = [];
    Eph(:,:,1)  = [];
    iono(:,1)   = [];
    
    if (~isempty(dir(filename_M_obs)))
        satObs = find( (pr1_R(:,1) ~= 0) & (pr1_M(:,1) ~= 0));
    else
        satObs = find(pr1_R(:,1) ~= 0);
    end
    satEph = find(sum(abs(Eph(:,:,1)))~=0);
end

epochs = length(time_GPS);

%remove observations without ephemerides
if (nargin >= 6)
    waitbar(0,wait_dlg,'Remove observations without ephemerides...')
end
for i = 1 : epochs
    satEph = find(sum(abs(Eph(:,:,i)))~=0);
    delsat = setdiff(1:num_sat,satEph);
    pr1_R(delsat,i)  = 0;
    pr1_M(delsat,i)  = 0;
    ph1_R(delsat,i)  = 0;
    ph1_M(delsat,i)  = 0;
    dop1_R(delsat,i) = 0;
    dop1_M(delsat,i) = 0;
    snr_R(delsat,i)  = 0;
    snr_M(delsat,i)  = 0;
    
    if (nargin >= 6)
        waitbar(i/epochs,wait_dlg)
    end
end

%complete/partial path
tMin = 1;
tMax = 1e30;
tMin = max(tMin,1);
tMax = min(tMax,epochs);
time_GPS = time_GPS(tMin:tMax);
week_R = week_R(tMin:tMax);
time_R = time_R(tMin:tMax);
time_M = time_M(tMin:tMax);
pr1_R = pr1_R(:,tMin:tMax);
pr1_M = pr1_M(:,tMin:tMax);
ph1_R = ph1_R(:,tMin:tMax);
ph1_M = ph1_M(:,tMin:tMax);
dop1_R = dop1_R(:,tMin:tMax);
dop1_M = dop1_M(:,tMin:tMax);
snr_R = snr_R(:,tMin:tMax);
snr_M = snr_M(:,tMin:tMax);
pos_M = pos_M(:,tMin:tMax);
Eph = Eph(:,:,tMin:tMax);
iono = iono(:,tMin:tMax);

%do not overwrite existing files
i = 1;
j = length(filerootOUT);
while (~isempty(dir([filerootOUT '_obs*.bin'])) || ...
        ~isempty(dir([filerootOUT '_eph*.bin'])) )
    
    filerootOUT(j+1:j+4) = ['_' num2str(i,'%03d')];
    i = i + 1;
end

%open output files
fid_obs = fopen([filerootOUT '_obs_000.bin'],'w+');
fid_eph = fopen([filerootOUT '_eph_000.bin'],'w+');

%write number of satellites
fwrite(fid_obs, num_sat, 'int8');
fwrite(fid_eph, num_sat, 'int8');

%"file hour" variable
hour = 0;

if (nargin >= 6)
    waitbar(0,wait_dlg,'Writing goGPS binary data...')
end

%write output files
for t = 1 : length(time_GPS)
    
    if (nargin >= 6)
        waitbar(t/epochs,wait_dlg)
    end
    
    %-------------------------------------
    % file management
    %-------------------------------------
    
    if (floor(t/3600) > hour)
        
        hour = floor(t/3600);
        hour_str = num2str(hour,'%03d');
        
        fclose(fid_obs);
        fclose(fid_eph);
        
        fid_obs    = fopen([filerootOUT '_obs_' hour_str '.bin'],'w+');
        fid_eph    = fopen([filerootOUT '_eph_' hour_str '.bin'],'w+');
        
        fwrite(fid_obs, num_sat, 'int8');
        fwrite(fid_eph, num_sat, 'int8');
    end
    
    Eph_t = Eph(:,:,t);
    fwrite(fid_obs, [time_GPS(t); time_M(t); time_R(t); week_R(t); pr1_M(:,t); pr1_R(:,t); ph1_M(:,t); ph1_R(:,t); dop1_R(:,t); snr_M(:,t); snr_R(:,t); pos_M(:,t); iono(:,t)], 'double');
    fwrite(fid_eph, [time_GPS(t); Eph_t(:)], 'double');
end

%close files
fclose(fid_obs);
fclose(fid_eph);
