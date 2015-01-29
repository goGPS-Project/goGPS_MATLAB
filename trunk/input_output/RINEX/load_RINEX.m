function [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
          dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
          snr2_R, snr2_M, time, time_R, time_M, week_R, week_M, ...
          date_R, date_M, pos_R, pos_M, Eph, iono, interval, antoff_R, antoff_M] = ...
          load_RINEX(filename_nav, filename_R_obs, filename_M_obs, constellations, flag_SP3, wait_dlg)

% SYNTAX:
%   [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%    dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
%    snr2_R, snr2_M, time, time_R, time_M, week_R, week_M, ...
%    date_R, date_M, pos_R, pos_M, Eph, iono, interval, antoff_R, antoff_R] = ...
%    load_RINEX(filename_nav, filename_R_obs, filename_M_obs, constellations, flag_SP3, wait_dlg);
%
% INPUT:
%   filename_nav = RINEX navigation file
%   filename_R_obs = RINEX observation file (ROVER)
%   filename_M_obs = RINEX observation file (MASTER) (empty if not available)
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%   flag_SP3 = boolean flag to indicate SP3 availability
%   wait_dlg = optional handler to waitbar figure (optional)
%
% OUTPUT:
%   pr1_R = code observation (L1 carrier, ROVER)
%   pr1_M = code observation (L1 carrier, MASTER)
%   ph1_R = phase observation (L1 carrier, ROVER)
%   ph1_M = phase observation (L1 carrier, MASTER)
%   pr2_R = code observation (L2 carrier, ROVER)
%   pr2_M = code observation (L2 carrier, MASTER)
%   ph2_R = phase observation (L2 carrier, ROVER)
%   ph2_M = phase observation (L2 carrier, MASTER)
%   dop1_R = Doppler observation (L1 carrier, ROVER)
%   dop1_M = Doppler observation (L1 carrier, MASTER)
%   dop2_R = Doppler observation (L2 carrier, ROVER)
%   dop2_M = Doppler observation (L2 carrier, MASTER)
%   snr1_R = signal-to-noise ratio (L1 carrier, ROVER)
%   snr1_M = signal-to-noise ratio (L1 carrier, MASTER)
%   snr2_R = signal-to-noise ratio (L2 carrier, ROVER)
%   snr2_M = signal-to-noise ratio (L2 carrier, MASTER)
%   time = reference time
%   time_R = rover time
%   time_M = master time
%   date = date (year,month,day,hour,minute,second)
%   pos_R = rover approximate position
%   pos_M = master station position
%   Eph = matrix containing 33 navigation parameters for each satellite
%   iono = vector containing ionosphere parameters
%   interval = observation rate [s]
%   antoff_R = antenna offset (ROVER)
%   antoff_M = antenna offset (MASTER)
%
% DESCRIPTION:
%   Parses RINEX files (observation and navigation) for both the ROVER
%   and the MASTER. Selects epochs they have in common.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni,Eugenio Realini
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
if (nargin < 6)
    wait_dlg_PresenceFlag = false;
else
    wait_dlg_PresenceFlag = true;
end
if (isempty(filename_M_obs))
    filename_M_obs_PresenceFlag = false;
else
    filename_M_obs_PresenceFlag = true;
end
if (isempty(constellations)) %then use only GPS as default
    [constellations] = multi_constellation_settings(1, 0, 0, 0, 0, 0);
end
tic;

%number of satellite slots for enabled constellations
nSatTot = constellations.nEnabledSat;

%fraction of INTERVAL (epoch-to-epoch timespan, as specified in the header)
%that is allowed as maximum difference between rover and master timings
%during synchronization
max_desync_frac = 0.1;

%read navigation files
if (~flag_SP3)
    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg,'Reading navigation files...')
    end
    
    Eph_G = []; iono_G = zeros(8,1);
    Eph_R = []; iono_R = zeros(8,1);
    Eph_E = []; iono_E = zeros(8,1);
    Eph_C = []; iono_C = zeros(8,1);
    Eph_J = []; iono_J = zeros(8,1);
    
    if (strcmpi(filename_nav(end),'p'))
        flag_mixed = 1;
    else
        flag_mixed = 0;
    end
    
    if (constellations.GPS.enabled || flag_mixed)
        if (exist(filename_nav,'file'))
            %parse RINEX navigation file (GPS) NOTE: filename expected to
            %end with 'n' or 'N' (GPS) or with 'p' or 'P' (mixed GNSS)
            [Eph_G, iono_G] = RINEX_get_nav(filename_nav, constellations);
        else
            fprintf('... WARNING: GPS navigation file not found. Disabling GPS positioning. \n');
            constellations.GPS.enabled = 0;
        end
    end
    
    if (constellations.GLONASS.enabled)
        if (exist([filename_nav(1:end-1) 'g'],'file'))
            %parse RINEX navigation file (GLONASS)
            [Eph_R, iono_R] = RINEX_get_nav([filename_nav(1:end-1) 'g'], constellations);
        elseif (~flag_mixed)
            fprintf('... WARNING: GLONASS navigation file not found. Disabling GLONASS positioning. \n');
            constellations.GLONASS.enabled = 0;
        end
    end
    
    if (constellations.Galileo.enabled)
        if (exist([filename_nav(1:end-1) 'l'],'file'))
            %parse RINEX navigation file (Galileo)
            [Eph_E, iono_E] = RINEX_get_nav([filename_nav(1:end-1) 'l'], constellations);
        elseif (~flag_mixed)
            fprintf('... WARNING: Galileo navigation file not found. Disabling Galileo positioning. \n');
            constellations.Galileo.enabled = 0;
        end
    end
    
    if (constellations.BeiDou.enabled)
        if (exist([filename_nav(1:end-1) 'b'],'file'))
            parse RINEX navigation file (BeiDou)
            [Eph_C, iono_C] = RINEX_get_nav([filename_nav(1:end-1) 'b'], constellations);
        elseif (~flag_mixed)
            fprintf('... WARNING: BeiDou navigation file not found. Disabling BeiDou positioning. \n');
            constellations.BeiDou.enabled = 0;
        end
    end
    
    if (constellations.QZSS.enabled)
        if (exist([filename_nav(1:end-1) 'q'],'file'))
            %parse RINEX navigation file (QZSS)
            [Eph_J, iono_J] = RINEX_get_nav([filename_nav(1:end-1) 'q'], constellations);
        elseif (~flag_mixed)
            fprintf('... WARNING: QZSS navigation file not found. Disabling QZSS positioning. \n');
            constellations.QZSS.enabled = 0;
        end
    end

    Eph = [Eph_G Eph_R Eph_E Eph_C Eph_J];
    
    if (any(iono_G))
        iono = iono_G;
    elseif (any(iono_R))
        iono = iono_R;
    elseif (any(iono_E))
        iono = iono_E;
    elseif (any(iono_C))
        iono = iono_C;
    elseif (any(iono_J))
        iono = iono_J;
    else
        iono = zeros(8,1);
        fprintf('... WARNING: ionosphere parameters not found in navigation file(s).\n');
    end
    
    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end
else
    Eph = zeros(33,nSatTot);
    iono = zeros(8,1);
end

%-------------------------------------------------------------------------------

%open RINEX observation file (ROVER)
FR_oss = fopen(filename_R_obs,'r');

if (filename_M_obs_PresenceFlag)
    %open RINEX observation file (MASTER)
    FM_oss = fopen(filename_M_obs,'r');
end

%-------------------------------------------------------------------------------

if (wait_dlg_PresenceFlag)
    waitbar(0.5,wait_dlg,'Parsing RINEX headers...')
end

%parse RINEX header
[obs_typ_R, pos_R, info_base_R, interval_R, sysId, antoff_R] = RINEX_parse_hdr(FR_oss);

%check RINEX version (rover)
if (isempty(sysId))
    RINEX_version_R = 2;
else
    RINEX_version_R = 3;
end

%check the availability of basic data to parse the RINEX file (ROVER)
if (info_base_R == 0)
    error('Basic data is missing in the ROVER RINEX header')
end

%find observation type columns
[obs_col_R, nObsTypes_R] = obs_type_find(obs_typ_R, sysId);

%number of lines to be read for each epoch (only for RINEX v2.xx)
if (RINEX_version_R == 2)
    nLinesToRead_R = ceil(nObsTypes_R/5);  %maximum of 5 obs per line
end

if (filename_M_obs_PresenceFlag)
    [obs_typ_M, pos_M, info_base_M, interval_M, sysId, antoff_M] = RINEX_parse_hdr(FM_oss);
    
    %check RINEX version (rover)
    if (isempty(sysId))
        RINEX_version_M = 2;
    else
        RINEX_version_M = 3;
    end
    
    %check the availability of basic data to parse the RINEX file (MASTER)
    if (info_base_M == 0)
        error('Basic data is missing in the ROVER RINEX header')
    end
    
    %find observation type columns
    [obs_col_M, nObsTypes_M] = obs_type_find(obs_typ_M, sysId);
    
    %number of lines to be read for each epoch (only for RINEX v2.xx)
    if (~isstruct(nObsTypes_M))
        nLinesToRead_M = ceil(nObsTypes_M/5);  %maximum of 5 obs per line
    end
else
    pos_M = zeros(3,1);
    interval_M = [];
end

if (wait_dlg_PresenceFlag)
    waitbar(1,wait_dlg)
end

interval = min([interval_R, interval_M]);

%-------------------------------------------------------------------------------

nEpochs = 86400;

%variable initialization (GPS)
time_R = zeros(nEpochs,1);
time_M = zeros(nEpochs,1);
pr1_R = zeros(nSatTot,nEpochs);
pr2_R = zeros(nSatTot,nEpochs);
ph1_R = zeros(nSatTot,nEpochs);
ph2_R = zeros(nSatTot,nEpochs);
dop1_R = zeros(nSatTot,nEpochs);
dop2_R = zeros(nSatTot,nEpochs);
snr1_R = zeros(nSatTot,nEpochs);
snr2_R = zeros(nSatTot,nEpochs);
pr1_M = zeros(nSatTot,nEpochs);
pr2_M = zeros(nSatTot,nEpochs);
ph1_M = zeros(nSatTot,nEpochs);
ph2_M = zeros(nSatTot,nEpochs);
snr1_M = zeros(nSatTot,nEpochs);
snr2_M = zeros(nSatTot,nEpochs);
dop1_M = zeros(nSatTot,nEpochs);
dop2_M = zeros(nSatTot,nEpochs);
date_R = zeros(nEpochs,6);
date_M = zeros(nEpochs,6);

%read data for the first epoch (ROVER)
[time_R(1), date_R(1,:), num_sat_R, sat_R, sat_types_R] = RINEX_get_epoch(FR_oss);

%-------------------------------------------------------------------------------

if (filename_M_obs_PresenceFlag)
    %read data for the first epoch (MASTER)
    [time_M(1), date_M(1,:), num_sat_M, sat_M, sat_types_M] = RINEX_get_epoch(FM_oss);
end
%-------------------------------------------------------------------------------

if (wait_dlg_PresenceFlag)
    waitbar(0.5,wait_dlg,'Parsing RINEX headers...')
end

if (filename_M_obs_PresenceFlag)
    while ((time_M(1) - time_R(1)) < 0 && abs(time_M(1) - time_R(1)) >= max_desync_frac*interval)
        
        %number of lines to be skipped
        if (RINEX_version_M == 2)
            nSkipLines = num_sat_M*nLinesToRead_M;
        else
            nSkipLines = num_sat_M;
        end
        
        %skip observations
        for s = 1 : nSkipLines
            fgetl(FM_oss);
        end
        
        %read data for the current epoch (MASTER)
        [time_M(1), date_M(1,:), num_sat_M, sat_M, sat_types_M] = RINEX_get_epoch(FM_oss);
    end
    
    while ((time_R(1) - time_M(1)) < 0 && abs(time_R(1) - time_M(1)) >= max_desync_frac*interval)

        %number of lines to be skipped
        if (RINEX_version_R == 2)
            nSkipLines = num_sat_R*nLinesToRead_R;
        else
            nSkipLines = num_sat_R;
        end
        
        %skip observations
        for s = 1 : nSkipLines
            fgetl(FR_oss);
        end
        
        %read data for the current epoch (ROVER)
        [time_R(1), date_R(1,:), num_sat_R, sat_R, sat_types_R] = RINEX_get_epoch(FR_oss);
    end
end

%read first batch of observations
%ROVER
obs_R = RINEX_get_obs(FR_oss, num_sat_R, sat_R, sat_types_R, obs_col_R, nObsTypes_R, constellations);

%read ROVER observations
if (sum(obs_R.P1 ~= 0) == sum(obs_R.C1 ~= 0))
    fprintf('Rover: using P1 code observations.\n');
    pr1_R(:,1) = obs_R.P1;
else
    fprintf('Rover: using C1 code observations.\n');
    pr1_R(:,1) = obs_R.C1;
end
pr2_R(:,1) = obs_R.P2;
ph1_R(:,1) = obs_R.L1;
ph2_R(:,1) = obs_R.L2;
dop1_R(:,1) = obs_R.D1;
dop2_R(:,1) = obs_R.D2;
snr1_R(:,1) = obs_R.S1;
snr2_R(:,1) = obs_R.S2;

if (filename_M_obs_PresenceFlag)
    %MASTER
    obs_M = RINEX_get_obs(FM_oss, num_sat_M, sat_M, sat_types_M, obs_col_M, nObsTypes_M, constellations);
    
    %read MASTER observations
    if (sum(obs_M.P1 ~= 0) == sum(obs_M.C1 ~= 0))
        fprintf('Master: using P1 code observations.\n');
        pr1_M(:,1) = obs_M.P1;
    else
        fprintf('Master: using C1 code observations.\n');
        pr1_M(:,1) = obs_M.C1;
    end
    pr2_M(:,1) = obs_M.P2;
    ph1_M(:,1) = obs_M.L1;
    ph2_M(:,1) = obs_M.L2;
    dop1_M(:,1) = obs_M.D1;
    dop2_M(:,1) = obs_M.D2;
    snr1_M(:,1) = obs_M.S1;
    snr2_M(:,1) = obs_M.S2;
end

if (wait_dlg_PresenceFlag)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

%define the reference time
time(1,1) = roundmod(time_R(1),interval);

if (wait_dlg_PresenceFlag)
    waitbar(0.5,wait_dlg,'Reading RINEX observations...')
end

k = 2;
n_mis_epo_R = 0;
n_mis_epo_M = 0;
while (~feof(FR_oss))

    if (abs((time_R(k-1) - time(k-1))) < max_desync_frac*interval)
        %display previously missing epochs, if any
        if (n_mis_epo_R > 0)
            for n = n_mis_epo_R : -1 : 1
                mis_epo = datevec(datenum(date_R(k-1,:))-datenum([0 0 0 0 0 n*interval]));
                fprintf('Missing epoch %d/%02d/%02d %02d:%02d:%02.3f (ROVER)\n', mis_epo);
            end
        end

        n_mis_epo_R = 0;
        %read data for the current epoch (ROVER)
        [time_R(k), date_R(k,:), num_sat_R, sat_R, sat_types_R] = RINEX_get_epoch(FR_oss);
    else
        n_mis_epo_R = n_mis_epo_R + 1;
        time_R(k) = time_R(k-1);
        date_R(k,:) = date_R(k-1,:);
        time_R(k-1) = 0;
        date_R(k-1,:) = [0 0 0 0 0 0];
    end

    if (filename_M_obs_PresenceFlag)
        if (abs((time_M(k-1) - time(k-1))) < max_desync_frac*interval)
            %display previously missing epochs, if any
            if (n_mis_epo_M > 0)
                for n = n_mis_epo_M : -1 : 1
                    mis_epo = datevec(datenum(date_M(k-1,:))-datenum([0 0 0 0 0 n*interval]));
                    fprintf('Missing epoch %d/%02d/%02d %02d:%02d:%02.3f (MASTER)\n', mis_epo);
                end
            end
            
            n_mis_epo_M = 0;
            %read data for the current epoch (MASTER)
            [time_M(k), date_M(k,:), num_sat_M, sat_M, sat_types_M] = RINEX_get_epoch(FM_oss);
        else
            n_mis_epo_M = n_mis_epo_M + 1;
            time_M(k) = time_M(k-1);
            date_M(k,:) = date_M(k-1,:);
            time_M(k-1) = 0;
            date_M(k-1,:) = [0 0 0 0 0 0];
        end
    end

    if (k > nEpochs)
        %variable initialization (GPS)
        pr1_R(:,k) = zeros(nSatTot,1);
        pr2_R(:,k) = zeros(nSatTot,1);
        ph1_R(:,k) = zeros(nSatTot,1);
        ph2_R(:,k) = zeros(nSatTot,1);
        dop1_R(:,k) = zeros(nSatTot,1);
        dop2_R(:,k) = zeros(nSatTot,1);
        snr1_R(:,k) = zeros(nSatTot,1);
        snr2_R(:,k) = zeros(nSatTot,1);
        pr1_M(:,k) = zeros(nSatTot,1);
        pr2_M(:,k) = zeros(nSatTot,1);
        ph1_M(:,k) = zeros(nSatTot,1);
        ph2_M(:,k) = zeros(nSatTot,1);
        snr1_M(:,k) = zeros(nSatTot,1);
        snr2_M(:,k) = zeros(nSatTot,1);
        dop1_M(:,k) = zeros(nSatTot,1);
        dop2_M(:,k) = zeros(nSatTot,1);

        nEpochs = nEpochs  + 1;
    end

    time(k,1) = time(k-1,1) + interval;

    if (abs(time_R(k)-time(k)) < max_desync_frac*interval)

        %read ROVER observations
        obs_R = RINEX_get_obs(FR_oss, num_sat_R, sat_R, sat_types_R, obs_col_R, nObsTypes_R, constellations);

        %read ROVER observations
        if (sum(obs_R.P1 ~= 0) == sum(obs_R.C1 ~= 0))
            pr1_R(:,k) = obs_R.P1;
        else
            pr1_R(:,k) = obs_R.C1;
        end
        pr2_R(:,k) = obs_R.P2;
        ph1_R(:,k) = obs_R.L1;
        ph2_R(:,k) = obs_R.L2;
        dop1_R(:,k) = obs_R.D1;
        dop2_R(:,k) = obs_R.D2;
        snr1_R(:,k) = obs_R.S1;
        snr2_R(:,k) = obs_R.S2;
%     else
%         %number of lines to be skipped
%         if (RINEX_version_R == 2)
%             nSkipLines = num_sat_R*nLinesToRead_R;
%         else
%             nSkipLines = num_sat_R;
%         end
%         
%         %skip observations
%         for s = 1 : nSkipLines
%             fgetl(FR_oss);
%         end
    end

    if (filename_M_obs_PresenceFlag)

        if (abs(time_M(k) - time(k)) < max_desync_frac*interval)
            
            %read MASTER observations
            obs_M = RINEX_get_obs(FM_oss, num_sat_M, sat_M, sat_types_M, obs_col_M, nObsTypes_M, constellations);
            
            %read MASTER observations
            if (sum(obs_M.P1 ~= 0) == sum(obs_M.C1 ~= 0))
                pr1_M(:,k) = obs_M.P1;
            else
                pr1_M(:,k) = obs_M.C1;
            end
            pr2_M(:,k) = obs_M.P2;
            ph1_M(:,k) = obs_M.L1;
            ph2_M(:,k) = obs_M.L2;
            dop1_M(:,k) = obs_M.D1;
            dop2_M(:,k) = obs_M.D2;
            snr1_M(:,k) = obs_M.S1;
            snr2_M(:,k) = obs_M.S2;
%         else
%             %number of lines to be skipped
%             if (RINEX_version_M == 2)
%                 nSkipLines = num_sat_M*nLinesToRead_M;
%             else
%                 nSkipLines = num_sat_M;
%             end
%             
%             %skip observations
%             for s = 1 : nSkipLines
%                 fgetl(FM_oss);
%             end
        end
    end
    
    k = k + 1;
end

%remove empty slots
time_R(k:nEpochs) = [];
time_M(k:nEpochs) = [];
pr1_R(:,k:nEpochs) = [];
pr2_R(:,k:nEpochs) = [];
ph1_R(:,k:nEpochs) = [];
ph2_R(:,k:nEpochs) = [];
dop1_R(:,k:nEpochs) = [];
dop2_R(:,k:nEpochs) = [];
snr1_R(:,k:nEpochs) = [];
snr2_R(:,k:nEpochs) = [];
pr1_M(:,k:nEpochs) = [];
pr2_M(:,k:nEpochs) = [];
ph1_M(:,k:nEpochs) = [];
ph2_M(:,k:nEpochs) = [];
snr1_M(:,k:nEpochs) = [];
snr2_M(:,k:nEpochs) = [];
dop1_M(:,k:nEpochs) = [];
dop2_M(:,k:nEpochs) = [];
date_R(k:nEpochs,:) = [];
date_M(k:nEpochs,:) = [];

%remove rover tail
if (filename_M_obs_PresenceFlag)
    flag_tail = 1;
    while (flag_tail)
        if (time_M(end) == 0)
            date_R(end,:) = [];
            date_M(end,:) = [];
            time(end) = [];
            time_R(end) = [];
            time_M(end) = [];
            pr1_R(:,end) = [];
            pr2_R(:,end) = [];
            ph1_R(:,end) = [];
            ph2_R(:,end) = [];
            dop1_R(:,end) = [];
            dop2_R(:,end) = [];
            snr1_R(:,end) = [];
            snr2_R(:,end) = [];
            pr1_M(:,end) = [];
            pr2_M(:,end) = [];
            ph1_M(:,end) = [];
            ph2_M(:,end) = [];
            snr1_M(:,end) = [];
            snr2_M(:,end) = [];
            dop1_M(:,end) = [];
            dop2_M(:,end) = [];
        else
            flag_tail = 0;
        end
    end
end

if (wait_dlg_PresenceFlag)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

%close RINEX files
fclose(FR_oss);
if (filename_M_obs_PresenceFlag)
    fclose(FM_oss);
end

%GPS week number
week_R = date2gps(date_R);
week_M = date2gps(date_M);

fprintf('The RINEX file has been read.\n');
toc
