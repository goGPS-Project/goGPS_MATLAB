function [pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, ...
          time_ref, time, week, date, pos, interval, antoff] = ...
          load_RINEX_obs(filename, constellations, wait_dlg)

% SYNTAX:
%   [pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, ...
%    time_ref, time, week, date, pos, interval, antoff] = ...
%    load_RINEX_obs(filename, constellations, wait_dlg);
%
% INPUT:
%   filename = RINEX observation file(s)
%   constellations = struct with multi-constellation settings
%                   (see 'multi_constellation_settings.m' - empty if not available)
%   wait_dlg = optional handler to waitbar figure (optional)
%
% OUTPUT:
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   dop1 = Doppler observation (L1 carrier)
%   dop2 = Doppler observation (L2 carrier)
%   snr1 = signal-to-noise ratio (L1 carrier)
%   snr2 = signal-to-noise ratio (L2 carrier)
%   time = receiver seconds-of-week
%   week = GPS week
%   date = date (year,month,day,hour,minute,second)
%   pos = rover approximate position
%   interval = observation time interval [s]
%   antoff = antenna offset [m]
%
% DESCRIPTION:
%   Parses RINEX observation files.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.1 beta
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
% Portions of code contributed by Damiano Triglione (2012)
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

% Check the input arguments
if (nargin < 3)
    wait_dlg_PresenceFlag = false;
else
    wait_dlg_PresenceFlag = true;
end

if (nargin < 2 || isempty(constellations)) %then use only GPS as default
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

%number of satellite slots for enabled constellations
nSatTot = constellations.nEnabledSat;

%number of RINEX files to be read
if (iscell(filename))
    nFiles = size(filename,1);
else
    nFiles = 1;
end

%variable initialization
nEpochs = 86400;
time = zeros(nEpochs,1,nFiles);
week = zeros(nEpochs,1,nFiles);
pr1 = zeros(nSatTot,nEpochs,nFiles);
pr2 = zeros(nSatTot,nEpochs,nFiles);
ph1 = zeros(nSatTot,nEpochs,nFiles);
ph2 = zeros(nSatTot,nEpochs,nFiles);
dop1 = zeros(nSatTot,nEpochs,nFiles);
dop2 = zeros(nSatTot,nEpochs,nFiles);
snr1 = zeros(nSatTot,nEpochs,nFiles);
snr2 = zeros(nSatTot,nEpochs,nFiles);
date = zeros(nEpochs,6,nFiles);
pos = zeros(3,1,nFiles);
interval = zeros(1,1,nFiles);
antoff = zeros(3,1,nFiles);

for f = 1 : nFiles

    if (iscell(filename))
        current_file = filename{f,1};
    else
        current_file = filename;
    end
    
    fprintf(['Reading RINEX file ' current_file ': ... ']);
    
    %open RINEX observation file
    fid = fopen(current_file,'r');
    
    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg,['RINEX file ' current_file ': parsing header...'])
    end
    
    %parse RINEX header
    [obs_type, pos(:,1,f), basic_info, interval(1,1,f), sysId, antoff(:,1,f)] = RINEX_parse_hdr(fid);
    
    %check the availability of basic data to parse the RINEX file
    if (basic_info == 0)
        error(['RINEX file ' current_file ': basic data is missing in the file header'])
    end
    
    %find observation type columns
    [obsColumns, nObsTypes] = obs_type_find(obs_type, sysId);
    
    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end
    
    %-------------------------------------------------------------------------------
    
    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg,['RINEX file ' current_file ': reading observations...'])
    end
    
    k = 1;
    while (~feof(fid))
        
        %read data for the current epoch (ROVER)
        [time(k,1,f), date(k,:,f), num_sat, sat, sat_types] = RINEX_get_epoch(fid);
        
        if (k > nEpochs)
            %variable initialization (GPS)
            pr1(:,k,f) = zeros(nSatTot,1);
            pr2(:,k,f) = zeros(nSatTot,1);
            ph1(:,k,f) = zeros(nSatTot,1);
            ph2(:,k,f) = zeros(nSatTot,1);
            dop1(:,k,f) = zeros(nSatTot,1);
            dop2(:,k,f) = zeros(nSatTot,1);
            snr1(:,k,f) = zeros(nSatTot,1);
            snr2(:,k,f) = zeros(nSatTot,1);
            
            nEpochs = nEpochs  + 1;
        end
        
        %read ROVER observations
        obs = RINEX_get_obs(fid, num_sat, sat, sat_types, obsColumns, nObsTypes, constellations);
        
        %read ROVER observations
        if (sum(obs.P1 ~= 0) == sum(obs.C1 ~= 0))
            pr1(:,k,f) = obs.P1;
        else
            pr1(:,k,f) = obs.C1;
        end
        pr2(:,k,f) = obs.P2;
        ph1(:,k,f) = obs.L1;
        ph2(:,k,f) = obs.L2;
        dop1(:,k,f) = obs.D1;
        dop2(:,k,f) = obs.D2;
        snr1(:,k,f) = obs.S1;
        snr2(:,k,f) = obs.S2;
        
        k = k + 1;
    end
    
    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end
    
    %GPS week number
    week(:,1,f) = date2gps(date(:,:,f));
    
    %observation rate
    if (interval(:,1,f) == 0)
        interval(:,1,f) = median(time(2:k-1,1,f) - time(1:k-2,1,f));
    end
    
    %-------------------------------------------------------------------------------
    
    %close RINEX files
    fclose(fid);
    
    fprintf('done\n');
end

%sync observations
[time_ref, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, interval] = ...
sync_obs(time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, interval);
