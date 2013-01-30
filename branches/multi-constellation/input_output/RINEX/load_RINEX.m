function [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
          dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
          snr2_R, snr2_M, time_GPS, time_GPS_R, time_GPS_M, ...
          date, pos_R, pos_M, Eph, iono, interval] = ...
          load_RINEX(filename_nav, filename_R_obs, filename_M_obs, constellations, flag_SP3, wait_dlg)

% SYNTAX:
%   [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%    dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
%    snr2_R, snr2_M, time_GPS, time_GPS_R, time_GPS_M, ...
%    date, pos_R, pos_M, Eph, iono, interval] = ...
%    load_RINEX(filename_nav, filename_R_obs, filename_M_obs, constellations, flag_SP3, wait_dlg);
%
% INPUT:
%   filename_nav = RINEX navigation file
%   filename_R_obs = RINEX observation file (ROVER)
%   filename_M_obs = RINEX observation file (MASTER) (empty if not available)
%   constellations = struct with multi-constellation settings (empty if not available)
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
%   time_GPS = reference GPS time
%   time_GPS_R = rover GPS time
%   time_GPS_M = master GPS time
%   date = date (year,month,day,hour,minute,second)
%   pos_R = rover approximate position
%   pos_M = master station position
%   Eph = matrix containing 29 ephemerides for each satellite
%   iono = vector containing ionosphere parameters
%
% DESCRIPTION:
%   Parses RINEX files (observation and navigation) for both the ROVER
%   and the MASTER. Selects epochs they have in common.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni,Eugenio Realini
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
    constellations.GPS = struct('numSat', 32, 'enabled', 1, 'indexes', [1:32]);
    constellations.nEnabledSat = 32;
end

%number of satellite slots for enabled constellations
nSatTot = constellations.nEnabledSat;

%fraction of INTERVAL (epoch-to-epoch timespan, as specified in the header)
%that is allowed as maximum difference between rover and master timings
%during synchronization
max_desync_frac = 0.1;

if (~flag_SP3)
    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg,'Reading navigation files...')
    end
    
    %parse RINEX navigation file (ROVER)
    [Eph, iono] = RINEX_get_nav(filename_nav);
    %[Eph, iono] = RINEX_get_nav_ORIGINALE(filename_nav);
    
    %parse RINEX navigation file (ROVER)
    % [Eph_R] = RINEX_get_nav_GLO(filename_nav_GLO);
    
    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end
else
    Eph = zeros(29,nSatTot);
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
[obs_typ_R,  pos_R, info_base_R, interval_R] = RINEX_parse_hdr(FR_oss);

%check the availability of basic data to parse the RINEX file (ROVER)
if (info_base_R == 0)
    error('Basic data is missing in the ROVER RINEX header')
end

if (filename_M_obs_PresenceFlag)
    [obs_typ_M, pos_M, info_base_M, interval_M] = RINEX_parse_hdr(FM_oss);
    
    %check the availability of basic data to parse the RINEX file (MASTER)
    if (info_base_M == 0)
        error('Basic data is missing in the ROVER RINEX header')
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
time_GPS_R = zeros(nEpochs,1);
time_GPS_M = zeros(nEpochs,1);
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
date = zeros(nEpochs,6);

%read data for the first epoch (ROVER)
[time_GPS_R(1), sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);

%read ROVER observations
obs_GPS_R = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R);

%read ROVER observations
if (obs_GPS_R.P1 == 0)
    pr1_R(:,1) = obs_GPS_R.C1;
else
    pr1_R(:,1) = obs_GPS_R.P1;
end
pr2_R(:,1) = obs_GPS_R.P2;
ph1_R(:,1) = obs_GPS_R.L1;
ph2_R(:,1) = obs_GPS_R.L2;
dop1_R(:,1) = obs_GPS_R.D1;
dop2_R(:,1) = obs_GPS_R.D2;
snr1_R(:,1) = obs_GPS_R.S1;
snr2_R(:,1) = obs_GPS_R.S2;

%-------------------------------------------------------------------------------

if (filename_M_obs_PresenceFlag)
    %read data for the first epoch (MASTER)
    [time_GPS_M(1), sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
    
    %read MASTER observations
    obs_GPS_M = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M);
    
    %read MASTER observations
    if (obs_GPS_M.P1 == 0)
        pr1_M(:,1) = obs_GPS_M.C1;
    else
        pr1_M(:,1) = obs_GPS_M.P1;
    end
    pr2_M(:,1) = obs_GPS_M.P2;
    ph1_M(:,1) = obs_GPS_M.L1;
    ph2_M(:,1) = obs_GPS_M.L2;
    dop1_M(:,1) = obs_GPS_M.D1;
    dop2_M(:,1) = obs_GPS_M.D2;
    snr1_M(:,1) = obs_GPS_M.S1;
    snr2_M(:,1) = obs_GPS_M.S2;

end
%-------------------------------------------------------------------------------

if (wait_dlg_PresenceFlag)
    waitbar(0.5,wait_dlg,'Parsing RINEX headers...')
end

if (filename_M_obs_PresenceFlag)
    while ((time_GPS_M(1) - time_GPS_R(1)) < 0 && abs(time_GPS_M(1) - time_GPS_R(1)) >= max_desync_frac*interval)
        
        %read data for the current epoch (MASTER)
        [time_GPS_M(1), sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
        
        %read MASTER observations
        obs_GPS_M = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M);
        
        %read MASTER observations
        if (obs_GPS_M.P1 == 0)
            pr1_M(:,1) = obs_GPS_M.C1;
        else
            pr1_M(:,1) = obs_GPS_M.P1;
        end
        pr2_M(:,1) = obs_GPS_M.P2;
        ph1_M(:,1) = obs_GPS_M.L1;
        ph2_M(:,1) = obs_GPS_M.L2;
        dop1_M(:,1) = obs_GPS_M.D1;
        dop2_M(:,1) = obs_GPS_M.D2;
        snr1_M(:,1) = obs_GPS_M.S1;
        snr2_M(:,1) = obs_GPS_M.S2;
    end
    
    while ((time_GPS_R(1) - time_GPS_M(1)) < 0 && abs(time_GPS_R(1) - time_GPS_M(1)) >= max_desync_frac*interval)
        
        %read data for the current epoch (ROVER)
        [time_GPS_R(1), sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);
        
        %read ROVER observations
        obs_GPS_R = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R);
        
        %read ROVER observations
        if (obs_GPS_R.P1 == 0)
            pr1_R(:,1) = obs_GPS_R.C1;
        else
            pr1_R(:,1) = obs_GPS_R.P1;
        end
        pr2_R(:,1) = obs_GPS_R.P2;
        ph1_R(:,1) = obs_GPS_R.L1;
        ph2_R(:,1) = obs_GPS_R.L2;
        dop1_R(:,1) = obs_GPS_R.D1;
        dop2_R(:,1) = obs_GPS_R.D2;
        snr1_R(:,1) = obs_GPS_R.S1;
        snr2_R(:,1) = obs_GPS_R.S2;
    end
end

if (wait_dlg_PresenceFlag)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------


time_GPS(1,1) = roundmod(time_GPS_R(1),interval);
date(1,:) = date_R(1,:);

if (wait_dlg_PresenceFlag)
    waitbar(0.5,wait_dlg,'Reading RINEX observations...')
end

k = 2;
while (~feof(FR_oss))

    if (abs((time_GPS_R(k-1) - time_GPS(k-1))) < max_desync_frac*interval)
        %read data for the current epoch (ROVER)
        [time_GPS_R(k), sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);
    else
        time_GPS_R(k) = time_GPS_R(k-1);
        if (time_GPS_R(k-1) ~= 0)
            fprintf('Missing epoch %f (ROVER)\n', time_GPS(k-1));
        end
        time_GPS_R(k-1) = 0;
    end

    if (filename_M_obs_PresenceFlag)
        if (abs((time_GPS_M(k-1) - time_GPS(k-1))) < max_desync_frac*interval)
            %read data for the current epoch (MASTER)
            [time_GPS_M(k), sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
        else
            time_GPS_M(k) = time_GPS_M(k-1);
            if (time_GPS_M(k-1) ~= 0)
                fprintf('Missing epoch %f (MASTER)\n', time_GPS(k-1));
            end
            time_GPS_M(k-1) = 0;
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
    
    date(k,:) = date_R(1,:);

    time_GPS(k,1) = time_GPS(k-1,1) + interval;
    
    if (abs(time_GPS_R(k)-time_GPS(k)) < max_desync_frac*interval)

        %read ROVER observations
        obs_GPS_R = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R);

        %read ROVER observations
        if (obs_GPS_R.P1 == 0)
            pr1_R(:,k) = obs_GPS_R.C1;
        else
            pr1_R(:,k) = obs_GPS_R.P1;
        end
        pr2_R(:,k) = obs_GPS_R.P2;
        ph1_R(:,k) = obs_GPS_R.L1;
        ph2_R(:,k) = obs_GPS_R.L2;
        dop1_R(:,k) = obs_GPS_R.D1;
        dop2_R(:,k) = obs_GPS_R.D2;
        snr1_R(:,k) = obs_GPS_R.S1;
        snr2_R(:,k) = obs_GPS_R.S2;
    end

    if (filename_M_obs_PresenceFlag)

        if (abs(time_GPS_M(k) - time_GPS(k)) < max_desync_frac*interval)
            
            %read MASTER observations
            obs_GPS_M = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M);
            
            %read MASTER observations
            if (obs_GPS_M.P1 == 0)
                pr1_M(:,k) = obs_GPS_M.C1;
            else
                pr1_M(:,k) = obs_GPS_M.P1;
            end
            pr2_M(:,k) = obs_GPS_M.P2;
            ph1_M(:,k) = obs_GPS_M.L1;
            ph2_M(:,k) = obs_GPS_M.L2;
            dop1_M(:,k) = obs_GPS_M.D1;
            dop2_M(:,k) = obs_GPS_M.D2;
            snr1_M(:,k) = obs_GPS_M.S1;
            snr2_M(:,k) = obs_GPS_M.S2;
        end
    end
    
    k = k+1;
end

%remove empty slots
time_GPS_R(k:nEpochs) = [];
time_GPS_M(k:nEpochs) = [];
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
date(k:nEpochs,:) = [];

%remove rover tail
if (filename_M_obs_PresenceFlag)
    flag_tail = 1;
    while (flag_tail)
        if (time_GPS_M(end) == 0)
            date(end,:) = [];
            time_GPS(end) = [];
            time_GPS_R(end) = [];
            time_GPS_M(end) = [];
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
