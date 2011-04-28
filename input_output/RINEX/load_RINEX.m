function [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
          dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
          snr2_R, snr2_M, pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
          dop1_RR, dop1_MR, dop2_RR, dop2_MR, snr_RR, snr_MR, ...
          time_GPS, time_GPS_R, time_GPS_M, date, pos, Eph, iono, Eph_R] = ...
          load_RINEX(filename_R_obs, filename_nav, filename_M_obs, wait_dlg)

% SYNTAX:
%   [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%    dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
%    snr1_R, snr1_M, pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
%    Eph_R, Eph_MR, snr_RR, snr_MR, ...
%    time_GPS, date, pos, Eph, iono, Eph_R] = ...
%   load_RINEX(filename_R_obs, filename_nav, filename_M_obs, wait_dlg);
%
% INPUT:
%   filename_R_obs = RINEX observation file (ROVER)
%   filename_M_obs = RINEX observation file (MASTER)
%   filename_nav = RINEX navigation file
%   wait_dlg = optional handler to waitbar figure
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
%   pr1_RR = code observation (GLONASS, L1 carrier, ROVER)
%   pr1_MR = code observation (GLONASS, L1 carrier, MASTER)
%   ph1_RR = phase observation (GLONASS, L1 carrier, ROVER)
%   ph1_MR = phase observation (GLONASS, L1 carrier, MASTER)
%   pr2_RR = code observation (GLONASS, L2 carrier, ROVER)
%   pr2_MR = code observation (GLONASS, L2 carrier, MASTER)
%   ph2_RR = phase observation (GLONASS, L2 carrier, ROVER)
%   ph2_MR = phase observation (GLONASS, L2 carrier, MASTER)
%   dop1_RR = Doppler observation (GLONASS, L1 carrier, ROVER)
%   dop1_MR = Doppler observation (GLONASS, L1 carrier, MASTER)
%   dop2_RR = Doppler observation (GLONASS, L2 carrier, ROVER)
%   dop2_MR = Doppler observation (GLONASS, L2 carrier, MASTER)
%   snr_RR = signal-to-noise ratio (GLONASS, ROVER)
%   snr_MR = signal-to-noise ratio (GLONASS, MASTER)
%   time_GPS = reference GPS time
%   time_GPS_R = rover GPS time
%   time_GPS_M = master GPS time
%   date = date (year,month,day,hour,minute,second)
%   pos = master station position (rover if stand-alone)
%   Eph = matrix containing 29 ephemerides for each satellite
%   iono = vector containing ionosphere parameters
%   Eph_R = matrix containing 29 ephemerides for each satellite (GLONASS)
%
% DESCRIPTION:
%   Parses RINEX files (observation and navigation) for both the ROVER
%   and the MASTER. Selects epochs they have in common.

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

Eph_R = zeros(17,32);

if (nargin == 5)
    waitbar(0.5,wait_dlg,'Reading navigation files...')
end

%parse RINEX navigation file (ROVER)
[Eph, iono] = RINEX_get_nav(filename_nav);

%parse RINEX navigation file (ROVER)
% [Eph_R] = RINEX_get_nav_GLO(filename_nav_GLO);

if (nargin == 5)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

%open RINEX observation file (ROVER)
FR_oss = fopen(filename_R_obs,'r');

if (nargin > 2)
    %open RINEX observation file (MASTER)
    FM_oss = fopen(filename_M_obs,'r');
end

%-------------------------------------------------------------------------------

if (nargin == 5)
    waitbar(0.5,wait_dlg,'Parsing RINEX headers...')
end

%parse RINEX header
[obs_typ_R,  pos, info_base_R] = RINEX_parse_hdr(FR_oss);

%check the availability of basic data to parse the RINEX file (ROVER)
if (info_base_R == 0)
    error('Basic data is missing in the ROVER RINEX header')
end

if (nargin > 2)
    [obs_typ_M, pos_M, info_base_M] = RINEX_parse_hdr(FM_oss);
    
    %check the availability of basic data to parse the RINEX file (MASTER)
    if (info_base_M == 0)
        error('Basic data is missing in the ROVER RINEX header')
    end
    
    pos = pos_M;
end

if (nargin == 5)
    waitbar(1,wait_dlg)
end

interval = 1;

%-------------------------------------------------------------------------------

nEpochs = 30000;

%variable initialization (GPS)
time_GPS_R = zeros(nEpochs,1);
time_GPS_M = zeros(nEpochs,1);
pr1_R = zeros(32,nEpochs);
pr2_R = zeros(32,nEpochs);
ph1_R = zeros(32,nEpochs);
ph2_R = zeros(32,nEpochs);
dop1_R = zeros(32,nEpochs);
dop2_R = zeros(32,nEpochs);
snr1_R = zeros(32,nEpochs);
snr2_R = zeros(32,nEpochs);
pr1_M = zeros(32,nEpochs);
pr2_M = zeros(32,nEpochs);
ph1_M = zeros(32,nEpochs);
ph2_M = zeros(32,nEpochs);
snr1_M = zeros(32,nEpochs);
snr2_M = zeros(32,nEpochs);
dop1_M = zeros(32,nEpochs);
dop2_M = zeros(32,nEpochs);

%variable initialization (GLONASS)
pr1_RR = zeros(32,1);
pr2_RR = zeros(32,1);
ph1_RR = zeros(32,1);
ph2_RR = zeros(32,1);
dop1_RR = zeros(32,1);
dop2_RR = zeros(32,1);
snr_RR = zeros(32,1);
pr1_MR = zeros(32,1);
pr2_MR = zeros(32,1);
ph1_MR = zeros(32,1);
ph2_MR = zeros(32,1);
dop1_MR = zeros(32,1);
dop2_MR = zeros(32,1);
snr_MR = zeros(32,1);

%read data for the first epoch (ROVER)
[time_GPS_R(1), sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);

%read ROVER observations
[obs_GPS_R, obs_GLO_R, obs_SBS_R] = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R); %#ok<NASGU>

%read ROVER observations (GPS)
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

%read ROVER observations (GLONASS)
% pr1_RR(:,1) = obs_GLO_R.C1;
% %pr2_RR(:,1) = obs_GLO_R.P2;
% ph1_RR(:,1) = obs_GLO_R.L1;
% %ph2_RR(:,1) = obs_GLO_R.L2;
% dop1_RR(:,1) = obs_GLO_R.D1;
% %dop2_RR(:,1) = obs_GLO_R.D2;
% snr_RR(:,1) = obs_GLO_R.S1;

%-------------------------------------------------------------------------------

if (nargin > 2)
    %read data for the first epoch (MASTER)
    [time_GPS_M(1), sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
    
    %read MASTER observations
    [obs_GPS_M, obs_GLO_M, obs_SBS_M] = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M); %#ok<NASGU>
    
    %read MASTER observations (GPS)
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
    
    %read MASTER observations (GLONASS)
    % pr1_MR(:,1) = obs_GLO_M.C1;
    % %pr2_MR(:,1) = obs_GLO_M.P2;
    % ph1_MR(:,1) = obs_GLO_M.L1;
    % %ph2_MR(:,1) = obs_GLO_M.L2;
    % dop1_MR(:,1) = obs_GLO_M.D1;
    % %dop2_MR(:,1) = obs_GLO_M.D2;
    % snr_MR(:,1) = obs_GLO_M.S1;
end
%-------------------------------------------------------------------------------

if (nargin == 5)
    waitbar(0.5,wait_dlg,'Parsing RINEX headers...')
end

if (nargin > 2)
    while (round(time_GPS_M(1)) < round(time_GPS_R(1)))
        
        %read data for the current epoch (MASTER)
        [time_GPS_M(1), sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
        
        %read MASTER observations
        [obs_GPS_M, obs_GLO_M, obs_SBS_M] = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M); %#ok<NASGU>
        
        %read MASTER observations (GPS)
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
        
        %read MASTER observations (GLONASS)
        % pr1_MR(:,1) = obs_GLO_M.C1;
        % %pr2_MR(:,1) = obs_GLO_M.P2;
        % ph1_MR(:,1) = obs_GLO_M.L1;
        % %ph2_MR(:,1) = obs_GLO_M.L2;
        % dop1_MR(:,1) = obs_GLO_M.D1;
        % %dop2_MR(:,1) = obs_GLO_M.D2;
        % snr_MR(:,1) = obs_GLO_M.S1;
    end
    
    while (round(time_GPS_R(1)) < round(time_GPS_M(1)))
        
        %read data for the current epoch (ROVER)
        [time_GPS_R(1), sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);
        
        %read ROVER observations
        [obs_GPS_R, obs_GLO_R, obs_SBS_R] = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R); %#ok<NASGU>
        
        %read ROVER observations (GPS)
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
        
        %read ROVER observations (GLONASS)
        % pr1_RR(:,1) = obs_GLO_R.C1;
        % %pr2_RR(:,1) = obs_GLO_R.P2;
        % ph1_RR(:,1) = obs_GLO_R.L1;
        % %ph2_RR(:,1) = obs_GLO_R.L2;
        % dop1_RR(:,1) = obs_GLO_R.D1;
        % %dop2_RR(:,1) = obs_GLO_R.D2;
        % snr_RR(:,1) = obs_GLO_R.S1;
    end
end

if (nargin == 5)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

k = 2;
time_GPS(1,1) = round(time_GPS_R(1));
date(1,:) = date_R(1,:);

if (nargin == 5)
    waitbar(0.5,wait_dlg,'Reading RINEX observations...')
end

while (~feof(FR_oss))

    if (round(time_GPS_R(k-1)) == time_GPS(k-1))
        %read data for the current epoch (ROVER)
        [time_GPS_R(k), sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);
    else
        time_GPS_R(k) = time_GPS_R(k-1);
        if (time_GPS_R(k-1) ~= 0)
            fprintf('Missing epoch %f (ROVER)\n', time_GPS(k-1));
        end
        time_GPS_R(k-1) = 0;
    end

    if (nargin > 2)
        if (round(time_GPS_M(k-1)) == time_GPS(k-1))
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
        pr1_R(:,k) = zeros(32,1);
        pr2_R(:,k) = zeros(32,1);
        ph1_R(:,k) = zeros(32,1);
        ph2_R(:,k) = zeros(32,1);
        dop1_R(:,k) = zeros(32,1);
        dop2_R(:,k) = zeros(32,1);
        snr1_R(:,k) = zeros(32,1);
        snr2_R(:,k) = zeros(32,1);
        pr1_M(:,k) = zeros(32,1);
        pr2_M(:,k) = zeros(32,1);
        ph1_M(:,k) = zeros(32,1);
        ph2_M(:,k) = zeros(32,1);
        snr1_M(:,k) = zeros(32,1);
        snr2_M(:,k) = zeros(32,1);
        dop1_M(:,k) = zeros(32,1);
        dop2_M(:,k) = zeros(32,1);
        
%         %variable initialization (GLONASS)
%         pr1_RR(:,k) = zeros(32,1);
%         pr2_RR(:,k) = zeros(32,1);
%         ph1_RR(:,k) = zeros(32,1);
%         ph2_RR(:,k) = zeros(32,1);
%         dop1_RR(:,k) = zeros(32,1);
%         dop2_RR(:,k) = zeros(32,1);
%         snr_RR(:,k) = zeros(32,1);
%         pr1_MR(:,k) = zeros(32,1);
%         pr2_MR(:,k) = zeros(32,1);
%         ph1_MR(:,k) = zeros(32,1);
%         ph2_MR(:,k) = zeros(32,1);
%         dop1_MR(:,k) = zeros(32,1);
%         dop2_MR(:,k) = zeros(32,1);
%         snr_MR(:,k) = zeros(32,1);

        nEpochs = nEpochs  + 1;
    end
    
    date(k,:) = date_R(1,:);

%     if((date_R(1,5) == 0 & round(date_R(1,6)) == 0) | (date_R(1,5) == 59 & round(date_R(1,6)) == 60))
%         fprintf('\n%dy %dm %dd %dh:%dm:%fs', date_R(1,1), date_R(1,2), date_R(1,3), date_R(1,4), date_R(1,5), date_R(1,6));
%     end
%     if(mod(round(date_R(1,5)),10) == 0 & (round(date_R(1,6)) == 0 | round(date_R(1,6)) == 60))
%         fprintf('.');
%     end

    time_GPS(k,1) = time_GPS(k-1,1) + interval;
    
    if (round(time_GPS_R(k)) == time_GPS(k))

        %read ROVER observations
        [obs_GPS_R, obs_GLO_R, obs_SBS_R] = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R); %#ok<NASGU>

        %read ROVER observations (GPS)
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

        %read ROVER observations (GLONASS)
        % pr1_RR(:,k) = obs_GLO_R.C1;
        % %pr2_RR(:,k) = obs_GLO_R.P2;
        % ph1_RR(:,k) = obs_GLO_R.L1;
        % %ph2_RR(:,k) = obs_GLO_R.L2;
        % dop1_RR(:,k) = obs_GLO_R.D1;
        % %dop2_RR(:,k) = obs_GLO_R.D2;
        % snr_RR(:,k) = obs_GLO_R.S1;
    end

    if (nargin > 2)

        if (round(time_GPS_M(k)) == time_GPS(k))
            
            %read MASTER observations
            [obs_GPS_M, obs_GLO_M, obs_SBS_M] = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M); %#ok<NASGU>
            
            %read MASTER observations (GPS)
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
            
            %read MASTER observations (GLONASS)
            % pr1_MR(:,k) = obs_GLO_M.C1;
            % %pr2_MR(:,k) = obs_GLO_M.P2;
            % ph1_MR(:,k) = obs_GLO_M.L1;
            % %ph2_MR(:,k) = obs_GLO_M.L2;
            % dop1_MR(:,k) = obs_GLO_M.D1;
            % %dop2_MR(:,k) = obs_GLO_M.D2;
            % snr_MR(:,k) = obs_GLO_M.S1;
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
% pr1_RR(:,k:nEpochs) = [];
% pr2_RR(:,k:nEpochs) = [];
% ph1_RR(:,k:nEpochs) = [];
% ph2_RR(:,k:nEpochs) = [];
% dop1_RR(:,k:nEpochs) = [];
% dop2_RR(:,k:nEpochs) = [];
% snr_RR(:,k:nEpochs) = [];
% pr1_MR(:,k:nEpochs) = [];
% pr2_MR(:,k:nEpochs) = [];
% ph1_MR(:,k:nEpochs) = [];
% ph2_MR(:,k:nEpochs) = [];
% dop1_MR(:,k:nEpochs) = [];
% dop2_MR(:,k:nEpochs) = [];
% snr_MR(:,k:nEpochs) = [];

%remove rover tail
if (nargin > 2)
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

if (nargin == 5)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

%close RINEX files
fclose(FR_oss);
if (nargin > 2)
    fclose(FM_oss);
end
