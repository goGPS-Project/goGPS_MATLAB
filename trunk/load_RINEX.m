function [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
          dop1_R, dop1_M, dop2_R, dop2_M, Eph_R, Eph_M, iono_R, iono_M, snr1_R, snr1_M, ...
          snr2_R, snr2_M, pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
          dop1_RR, dop1_MR, dop2_RR, dop2_MR, Eph_RR, Eph_MR, snr_RR, snr_MR, ...
          time_GPS, time_GPS_R, time_GPS_M, date, pos] = ...
          load_RINEX(nome_FR_oss, nome_FR_nav, nome_FM_oss, nome_FM_nav, wait_dlg)

% SYNTAX:
%   [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%    dop1_R, dop1_M, dop2_R, dop2_M, Eph_R, Eph_M, iono_R, iono_M, snr1_R, snr1_M, ...
%    snr1_R, snr1_M, pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
%    Eph_RR, Eph_MR, snr_RR, snr_MR, ...
%    time_GPS, date, pos] = ...
%   load_RINEX(nome_FR_oss, nome_FR_nav, nome_FM_oss, nome_FM_nav, wait_dlg);
%
% INPUT:
%   nome_FR_oss = RINEX observation file (ROVER)
%   nome_FR_nav = RINEX navigation file (ROVER)
%   nome_FM_oss = RINEX observation file (MASTER)
%   nome_FM_nav = RINEX navigation file (MASTER)
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
%   Eph_R = matrix containing 29 ephemerides for each satellite (ROVER)
%   Eph_M = matrix containing 29 ephemerides for each satellite (MASTER)
%   iono_R = matrix containing ionosphere parameters (ROVER)
%   iono_M = matrix containing ionosphere parameters (MASTER)
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
%   Eph_RR = matrix containing 29 ephemerides for each satellite (GLONASS, ROVER)
%   Eph_MR = matrix containing 29 ephemerides for each satellite (GLONASS, MASTER)
%   snr_RR = signal-to-noise ratio (GLONASS, ROVER)
%   snr_MR = signal-to-noise ratio (GLONASS, MASTER)
%   time_GPS = reference GPS time
%   time_GPS_R = rover GPS time
%   time_GPS_M = master GPS time
%   date = date (year,month,day,hour,minute,second)
%   pos = master station position (rover if stand-alone)
%
% DESCRIPTION:
%   Parse RINEX files (both observation and navigation) for both the ROVER
%   and the MASTER. Select epochs they have in common.

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

Eph_RR = zeros(17,32);
Eph_MR = zeros(17,32);

time_GPS_M = 0;
Eph_M = zeros(29,32);
iono_M = zeros(8,1);

if (nargin == 5)
    waitbar(0.33,wait_dlg,'Reading navigation files...')
end

%parse RINEX navigation file (ROVER)
[Eph_R, iono_R] = RINEX_get_nav(nome_FR_nav);

%parse RINEX navigation file (ROVER)
% [Eph_RR] = RINEX_get_nav_GLO(nome_FR_glo);

if (nargin == 5)
    waitbar(0.66,wait_dlg)
end

if (nargin > 2)
    %parse RINEX navigation file (MASTER)
    [Eph_M, iono_M] = RINEX_get_nav(nome_FM_nav);
    
    %parse RINEX navigation file (MASTER)
    % [Eph_MR] = RINEX_get_nav_GLO(nome_FM_glo);
end

if (nargin == 5)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

%open RINEX observation file (ROVER)
FR_oss = fopen(nome_FR_oss,'r');

if (nargin > 2)
    %open RINEX observation file (MASTER)
    FM_oss = fopen(nome_FM_oss,'r');
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

%variable initialization (GPS)
pr1_R(:,1) = zeros(32,1);
pr2_R(:,1) = zeros(32,1);
ph1_R(:,1) = zeros(32,1);
ph2_R(:,1) = zeros(32,1);
dop1_R(:,1) = zeros(32,1);
dop2_R(:,1) = zeros(32,1);
snr1_R(:,1) = zeros(32,1);
snr2_R(:,1) = zeros(32,1);
pr1_M(:,1) = zeros(32,1);
pr2_M(:,1) = zeros(32,1);
ph1_M(:,1) = zeros(32,1);
ph2_M(:,1) = zeros(32,1);
snr1_M(:,1) = zeros(32,1);
snr2_M(:,1) = zeros(32,1);
dop1_M(:,1) = zeros(32,1);
dop2_M(:,1) = zeros(32,1);

%variable initialization (GLONASS)
pr1_RR(:,1) = zeros(32,1);
pr2_RR(:,1) = zeros(32,1);
ph1_RR(:,1) = zeros(32,1);
ph2_RR(:,1) = zeros(32,1);
dop1_RR(:,1) = zeros(32,1);
dop2_RR(:,1) = zeros(32,1);
snr_RR(:,1) = zeros(32,1);
pr1_MR(:,1) = zeros(32,1);
pr2_MR(:,1) = zeros(32,1);
ph1_MR(:,1) = zeros(32,1);
ph2_MR(:,1) = zeros(32,1);
dop1_MR(:,1) = zeros(32,1);
dop2_MR(:,1) = zeros(32,1);
snr_MR(:,1) = zeros(32,1);

%read data for the first epoch (ROVER)
[time_GPS_R, sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);

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
    [time_GPS_M, sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
    
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
    while (round(time_GPS_M) < round(time_GPS_R))
        
        %read data for the current epoch (MASTER)
        [time_GPS_M, sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
        
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
    
    while (round(time_GPS_R) < round(time_GPS_M))
        
        %read data for the current epoch (ROVER)
        [time_GPS_R, sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);
        
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
time_GPS(1,1) = round(time_GPS_R);
date(1,:) = date_R(1,:);

if (nargin == 5)
    waitbar(0.5,wait_dlg,'Reading RINEX observations...')
end

while (~feof(FR_oss))

    %read data for the current epoch (ROVER)
    [time_GPS_R(k), sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);

    if (nargin > 2)
        %read data for the current epoch (MASTER)
        [time_GPS_M(k), sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>

        if (round(time_GPS_R(k)) > round(time_GPS_M(k)))
            %ignore rover tail
            time_GPS_R(k) = [];
            time_GPS_M(k) = [];
            break
        end
    end
    
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

    %variable initialization (GLONASS)
    pr1_RR(:,k) = zeros(32,1);
    pr2_RR(:,k) = zeros(32,1);
    ph1_RR(:,k) = zeros(32,1);
    ph2_RR(:,k) = zeros(32,1);
    dop1_RR(:,k) = zeros(32,1);
    dop2_RR(:,k) = zeros(32,1);
    snr_RR(:,k) = zeros(32,1);
    pr1_MR(:,k) = zeros(32,1);
    pr2_MR(:,k) = zeros(32,1);
    ph1_MR(:,k) = zeros(32,1);
    ph2_MR(:,k) = zeros(32,1);
    dop1_MR(:,k) = zeros(32,1);
    dop2_MR(:,k) = zeros(32,1);
    snr_MR(:,k) = zeros(32,1);
    
    date(k,:) = date_R(1,:);
    
    %read ROVER observations
    [obs_GPS_R, obs_GLO_R, obs_SBS_R] = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R); %#ok<NASGU>
    
    time_GPS(k,1) = time_GPS(k-1,1) + interval;

    if (round(time_GPS_R(k)) == time_GPS(k))

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
        
        %read MASTER observations
        [obs_GPS_M, obs_GLO_M, obs_SBS_M] = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M); %#ok<NASGU>

        if (round(time_GPS_M(k)) == time_GPS(k))
            
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

if (nargin == 5)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

%close RINEX files
fclose(FR_oss);
if (nargin > 2)
    fclose(FM_oss);
end
