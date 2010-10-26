function [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
          Eph_R, Eph_M, iono_R, iono_M, snr_R, snr_M, ...
          pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
          Eph_RR, Eph_MR, snr_RR, snr_MR, ...
          time_GPS, date, pos_M] = ...
          load_RINEX(nome_FR_oss, nome_FR_nav, nome_FM_oss, nome_FM_nav, wait_dlg)

% SYNTAX:
%   [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%    Eph_R, Eph_M, iono_R, iono_M, snr_R, snr_M, ...
%    pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
%    Eph_RR, Eph_MR, snr_RR, snr_MR, ...
%    time_GPS, date, pos_M] = ...
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
%   Eph_R = matrix containing 29 ephemerides for each satellite (ROVER)
%   Eph_M = matrix containing 29 ephemerides for each satellite (MASTER)
%   iono_R = matrix containing ionosphere parameters (ROVER)
%   iono_M = matrix containing ionosphere parameters (MASTER)
%   time_GPS = GPS time of ROVER observations
%   date = date (year,month,day,hour,minute,second)
%   pos_M = master station approximate position
%
% DESCRIPTION:
%   Parse RINEX files (both observation and navigation) for both the ROVER
%   and the MASTER. Select epochs they have in common.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
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

Eph_RR = zeros(17,32);
Eph_MR = zeros(17,32);

Eph_M = zeros(29,32);
iono_M = zeros(8,1);
pos_M = zeros(3,1);

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
[obs_typ_R,  null, info_base_R] = RINEX_parse_hdr(FR_oss); %#ok<ASGLU>

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
end

if (nargin == 5)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

%read data for the first epoch (ROVER)
[time_GPS_R, sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);

%read ROVER observations
[obs_GPS_R, obs_GLO_R, obs_SBS_R] = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R); %#ok<NASGU>

%-------------------------------------------------------------------------------

if (nargin > 2)
    %read data for the first epoch (MASTER)
    [time_GPS_M, sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
    
    %read MASTER observations
    [obs_GPS_M, obs_GLO_M, obs_SBS_M] = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M); %#ok<NASGU>
end
%-------------------------------------------------------------------------------

if (nargin == 5)
    waitbar(0.5,wait_dlg,'Parsing RINEX headers...')
end

if (nargin > 2)
    while (time_GPS_M < time_GPS_R)
        
        %read data for the current epoch (MASTER)
        [time_GPS_M, sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
        
        %read MASTER observations
        [obs_GPS_M, obs_GLO_M, obs_SBS_M] = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M); %#ok<NASGU>
    end
    
    while (time_GPS_R < time_GPS_M)
        
        %read data for the current epoch (ROVER)
        [time_GPS_R, sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);
        
        %read ROVER observations
        [obs_GPS_R, obs_GLO_R, obs_SBS_R] = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R); %#ok<NASGU>
    end
end

if (nargin == 5)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

k = 1;
time_GPS(1,1) = time_GPS_R;
date(1,:) = date_R(1,:);

if (nargin == 5)
    waitbar(0.5,wait_dlg,'Reading RINEX observations...')
end

while (~feof(FR_oss))

    %variable initialization (GPS)
    pr1_R(:,k) = zeros(32,1);
    pr2_R(:,k) = zeros(32,1);
    ph1_R(:,k) = zeros(32,1);
    ph2_R(:,k) = zeros(32,1);
    snr_R(:,k) = zeros(32,1);
    pr1_M(:,k) = zeros(32,1);
    pr2_M(:,k) = zeros(32,1);
    ph1_M(:,k) = zeros(32,1);
    ph2_M(:,k) = zeros(32,1);
    snr_M(:,k) = zeros(32,1);

    %variable initialization (GLONASS)
    pr1_RR(:,k) = zeros(32,1);
    pr2_RR(:,k) = zeros(32,1);
    ph1_RR(:,k) = zeros(32,1);
    ph2_RR(:,k) = zeros(32,1);
    snr_RR(:,k) = zeros(32,1);
    pr1_MR(:,k) = zeros(32,1);
    pr2_MR(:,k) = zeros(32,1);
    ph1_MR(:,k) = zeros(32,1);
    ph2_MR(:,k) = zeros(32,1);
    snr_MR(:,k) = zeros(32,1);

    if (time_GPS_R == time_GPS(k))

        %read ROVER observations (GPS)
        pr1_R(:,k) = obs_GPS_R.C1;
        pr2_R(:,k) = obs_GPS_R.P2;
        ph1_R(:,k) = obs_GPS_R.L1;
        ph2_R(:,k) = obs_GPS_R.L2;
        snr_R(:,k) = obs_GPS_R.S1;
        
        %read ROVER observations (GLONASS)
        % pr1_RR(:,k) = obs_GLO_R.C1;
        % %pr2_RR(:,k) = obs_GLO_R.P2;
        % ph1_RR(:,k) = obs_GLO_R.L1;
        % %ph2_RR(:,k) = obs_GLO_R.L2;
        % snr_RR(:,k) = obs_GLO_R.S1;

        %read data for the current epoch (ROVER)
        [time_GPS_R, sat_R, sat_types_R, date_R] = RINEX_get_epoch(FR_oss);

        %read ROVER observations
        [obs_GPS_R, obs_GLO_R, obs_SBS_R] = RINEX_get_obs(FR_oss, sat_R, sat_types_R, obs_typ_R); %#ok<NASGU>

    end

    if (nargin > 2)
        if (time_GPS_M == time_GPS(k))
            
            %read MASTER observations (GPS)
            pr1_M(:,k) = obs_GPS_M.C1;
            pr2_M(:,k) = obs_GPS_M.P2;
            ph1_M(:,k) = obs_GPS_M.L1;
            ph2_M(:,k) = obs_GPS_M.L2;
            snr_M(:,k) = obs_GPS_M.S1;
            
            %read MASTER observations (GLONASS)
            % pr1_MR(:,k) = obs_GLO_M.C1;
            % %pr2_MR(:,k) = obs_GLO_M.P2;
            % ph1_MR(:,k) = obs_GLO_M.L1;
            % %ph2_MR(:,k) = obs_GLO_M.L2;
            % snr_MR(:,k) = obs_GLO_M.S1;
            
            %read data for the current epoch (MASTER)
            [time_GPS_M, sat_M, sat_types_M, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>
            
            %read MASTER observations
            [obs_GPS_M, obs_GLO_M, obs_SBS_M] = RINEX_get_obs(FM_oss, sat_M, sat_types_M, obs_typ_M); %#ok<NASGU>
            
        end
        
        %ignore rover tail
        if (time_GPS_R > time_GPS_M)
            break
        end
    end

    k = k+1;
    time_GPS(k,1) = time_GPS(k-1,1) + 1;
    date(k,:) = date_R(1,:);

end

if (nargin == 5)
    waitbar(1,wait_dlg)
end

time_GPS(end) = [];

%-------------------------------------------------------------------------------

%close RINEX files
fclose(FR_oss);
if (nargin > 2)
    fclose(FM_oss);
end
