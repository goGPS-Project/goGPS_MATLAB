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
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
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

if (isempty(dir(filename_R_nav)) & isempty(dir(filename_M_nav)))
    if (nargin == 6)
        msgbox('Navigation data (ephemerides) are required to create goGPS binary data.');
    else
        fprintf('Navigation data (ephemerides) are required to create goGPS binary data.\n');
    end
    return
elseif (isempty(dir(filename_R_nav)))
    filename_R_nav = filename_M_nav;
elseif (isempty(dir(filename_M_nav)))
    filename_M_nav = filename_R_nav;
end

if (~isempty(dir(filename_R_obs)))
    
    %ROVER RINEX files reading
    if (nargin == 6)
        [pr1_R, ph1_R, pr2_R, ph2_R, ...
            Eph_R, iono_R, snr_R, ...
            pr1_RR, ph1_RR, pr2_RR, ph2_RR, ...
            Eph_RR, snr_RR, time_R, date] = ...
            load_RINEX_SA(filename_R_obs, filename_R_nav, wait_dlg); %#ok<ASGLU>
    else
        [pr1_R, ph1_R, pr2_R, ph2_R, ...
            Eph_R, iono_R, snr_R, ...
            pr1_RR, ph1_RR, pr2_RR, ph2_RR, ...
            Eph_RR, snr_RR, time_R, date] = ...
            load_RINEX_SA(filename_R_obs, filename_R_nav); %#ok<ASGLU>
    end
    
    %GPS week number
    date(:,1) = date(:,1) + 2000;
    week_R = floor((datenum(date) - datenum([1980,1,6,0,0,0]))/7);
else
    if (nargin == 6)
        msgbox('Rover data are required to create goGPS binary data.');
    else
        fprintf('Rover data are required to create goGPS binary data.\n');
    end
    
    return
end

time_GPS = time_R;
epochs = length(time_GPS);

%MASTER variable initialization
time_M = zeros(epochs);
pr1_M = zeros(32,epochs);
ph1_M = zeros(32,epochs);
snr_M = zeros(32,epochs);
pos_M = zeros(3,epochs);
Eph_M = zeros(29,32,epochs);
iono = zeros(8,epochs);

if (~isempty(dir(filename_M_obs)))
    
    %MASTER RINEX files reading
    if (nargin == 6)
        [pr1_M, ph1_M, pr2_M, ph2_M, ...
            Eph_M, iono_M, snr_M, ...
            pr1_MR, ph1_MR, pr2_MR, ph2_MR, ...
            Eph_MR, snr_MR, time_M, date, pos_M] = ...
            load_RINEX_SA(filename_M_obs, filename_M_nav, wait_dlg, time_GPS(end)); %#ok<ASGLU>
    else
        [pr1_M, ph1_M, pr2_M, ph2_M, ...
            Eph_M, iono_M, snr_M, ...
            pr1_MR, ph1_MR, pr2_MR, ph2_MR, ...
            Eph_MR, snr_MR, time_M, date, pos_M] = ...
            load_RINEX_SA(filename_M_obs, filename_M_nav); %#ok<ASGLU>
    end

    %--------------------------------------------------------------------------
    
    if (nargin == 2)
        waitbar(0.33,wait_dlg,'Synchronizing data...')
    end

    if ~isempty(time_R) & ~isempty(time_M)

        if (size(pos_M,2) == 1)
            pos_M(1,1:length(time_M)) = pos_M(1);
            pos_M(2,1:length(time_M)) = pos_M(2);
            pos_M(3,1:length(time_M)) = pos_M(3);
        end

        %head synchronization
        if (time_R(1) < time_M(1))
            pos = find(time_R < time_M(1));
            time_R(pos)    = [];                       %GPS time
            week_R(pos)    = [];                       %GPS week
            pr1_R(:,pos)   = [];                       %code observations
            ph1_R(:,pos)   = [];                       %phase observations
            snr_R(:,pos)   = [];                       %signal-to-noise ratio
            iono(:,pos)    = [];                       %ionosphere parameters
        end

        if (time_M(1) < time_R(1))
            pos = find(time_M < time_R(1));
            time_M(pos)    = [];                       %GPS time
            pr1_M(:,pos)   = [];                       %code observations
            ph1_M(:,pos)   = [];                       %phase observations
            snr_M(:,pos)   = [];                       %signal-to-noise ratio
            pos_M(:,pos)   = [];                       %master station position
        end
        
        %tail synchronization
        if (time_R(end) > time_M(end))
            pos = find(time_R > time_M(end));
            time_R(pos)    = [];                       %GPS time
            week_R(pos)    = [];                       %GPS week
            pr1_R(:,pos)   = [];                       %code observations
            ph1_R(:,pos)   = [];                       %phase observations
            snr_R(:,pos)   = [];                       %signal-to-noise ratio
            iono(:,pos)    = [];                       %ionosphere parameters
        end
        
        if (time_M(end) > time_R(end))
            pos = find(time_M > time_R(end));
            time_M(pos)    = [];                       %GPS time
            pr1_M(:,pos)   = [];                       %code observations
            ph1_M(:,pos)   = [];                       %phase observations
            snr_M(:,pos)   = [];                       %signal-to-noise ratio
            pos_M(:,pos)   = [];                       %master station position
        end
        
    end
    
    %-------------------------------------------------------------------------------
    
    if (nargin == 2)
        waitbar(0.66,wait_dlg)
    end
    
    %signal losses
    time_GPS = union(time_R,time_M);                     %overall reference time
    if ~isempty(time_GPS)
        
        time_GPS = (time_GPS(1) : 1 : time_GPS(end))';   %GPS time without interruptions
        
        if ~isempty(time_R)
            
            newtime_R = setdiff(time_GPS, time_R);       %ROVER missing epochs
            for i = 1 : length(newtime_R)
                
                pos = find(time_R == newtime_R(i) - 1);  %position before the "holes"

                time_R = [time_R(1:pos)   newtime_R(i)   time_R(pos+1:end)];
                week_R = [week_R(1:pos);  0;             week_R(pos+1:end)];
                pr1_R  = [pr1_R(:,1:pos)  zeros(32,1)    pr1_R(:,pos+1:end)];
                ph1_R  = [ph1_R(:,1:pos)  zeros(32,1)    ph1_R(:,pos+1:end)];
                snr_R  = [snr_R(:,1:pos)  zeros(32,1)    snr_R(:,pos+1:end)];
                iono   = [iono(:,1:pos)   zeros(8,1)     iono(:,pos+1:end)];
            end
        else
            time_R = time_GPS;
            week_R = zeros(1,length(time_GPS));
            pr1_R  = zeros(32,length(time_GPS));
            ph1_R  = zeros(32,length(time_GPS));
            snr_R  = zeros(32,length(time_GPS));
            iono   = zeros(8,length(time_GPS));
        end
        
        if ~isempty(time_M)
            
            newtime_M = setdiff(time_GPS, time_M);       %MASTER missing epochs
            for i = 1 : length(newtime_M)
                
                pos = find(time_M == newtime_M(i) - 1);  %position before the "holes"
                
                time_M = [time_M(1:pos)   newtime_M(i)   time_M(pos+1:end)];
                pr1_M  = [pr1_M(:,1:pos)  zeros(32,1)    pr1_M(:,pos+1:end)];
                ph1_M  = [ph1_M(:,1:pos)  zeros(32,1)    ph1_M(:,pos+1:end)];
                snr_M  = [snr_M(:,1:pos)  zeros(32,1)    snr_M(:,pos+1:end)];
                pos_M  = [pos_M(:,1:pos)  zeros(3,1)     pos_M(:,pos+1:end)];
            end
        else
            time_M = time_GPS;
            pr1_M  = zeros(32,length(time_GPS));
            ph1_M  = zeros(32,length(time_GPS));
            snr_M  = zeros(32,length(time_GPS));
            pos_M  = zeros(3,length(time_GPS));
        end
    end
    
    if (nargin == 2)
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
for i = 1 : epochs
    Eph(:,:,i) = rt_find_eph(Eph_tmp, time_GPS(i));
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

%remove epochs in which not all ephemerides are available or available
%satellites are < 4
satEph = find(sum(abs(Eph(:,:,1))) ~= 0);
if (~isempty(dir(filename_M_obs)))
    satObs = find( (pr1_R(:,1) ~= 0) & (pr1_M(:,1) ~= 0));
else
    satObs = find(pr1_R(:,1) ~= 0);
end
while (length(satEph) < length(satObs)) | (length(satObs) < 4)

    time_GPS(1) = [];
    week_R(1)   = [];
    time_R(1)   = [];
    time_M(1)   = [];
    pr1_R(:,1)  = [];
    pr1_M(:,1)  = [];
    ph1_R(:,1)  = [];
    ph1_M(:,1)  = [];
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
for i = 1 : epochs
    satEph = find(sum(abs(Eph(:,:,i)))~=0);
    delsat = setdiff(1:32,satEph);
    pr1_R(delsat,i) = 0;
    pr1_M(delsat,i) = 0;
    ph1_R(delsat,i) = 0;
    ph1_M(delsat,i) = 0;
    snr_R(delsat,i) = 0;
    snr_M(delsat,i) = 0;
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
snr_R = snr_R(:,tMin:tMax);
snr_M = snr_M(:,tMin:tMax);
pos_M = pos_M(:,tMin:tMax);
Eph = Eph(:,:,tMin:tMax);
iono = iono(:,tMin:tMax);

%do not overwrite existing files
i = 1;
j = length(filerootOUT);
while (~isempty(dir([filerootOUT '_obs*.bin'])) | ...
        ~isempty(dir([filerootOUT '_eph*.bin'])) )
    
    filerootOUT(j+1:j+3) = ['_' num2str(i,'%02d')];
    i = i + 1;
end

%open output files
fid_obs = fopen([filerootOUT '_obs_00.bin'],'w+');
fid_eph = fopen([filerootOUT '_eph_00.bin'],'w+');

%"file hour" variable
hour = 0;

if (nargin == 6)
    waitbar(0,wait_dlg,'Writing goGPS binary data...')
end

%write output files
for t = 1 : length(time_GPS)
    
    if (nargin == 6)
        waitbar(t/epochs,wait_dlg)
    end
    
    %-------------------------------------
    % file management
    %-------------------------------------
    
    if (floor(t/3600) > hour)
        
        hour = floor(t/3600);
        hour_str = num2str(hour,'%02d');
        
        fclose(fid_obs);
        fclose(fid_eph);
        
        fid_obs    = fopen([filerootOUT '_obs_'    hour_str '.bin'],'w+');
        fid_eph    = fopen([filerootOUT '_eph_'    hour_str '.bin'],'w+');
        
    end
    
    Eph_t = Eph(:,:,t);
    fwrite(fid_obs, [time_GPS(t); time_M(t); time_R(t); week_R(t); pr1_M(:,t); pr1_R(:,t); ph1_M(:,t); ph1_R(:,t); snr_M(:,t); snr_R(:,t); pos_M(:,t); iono(:,t)], 'double');
    fwrite(fid_eph, [time_GPS(t); Eph_t(:)], 'double');
end

%close files
fclose(fid_obs);
fclose(fid_eph);