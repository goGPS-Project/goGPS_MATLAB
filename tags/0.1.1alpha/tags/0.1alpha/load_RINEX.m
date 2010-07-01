function [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
          Eph_R, Eph_M, iono_R, iono_M, snr_R, snr_M, ...
          pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
          Eph_RR, Eph_MR, snr_RR, snr_MR, ...
          time_GPS, date] = ...
          load_RINEX(nome_FR_oss, nome_FR_nav, nome_FM_oss, nome_FM_nav)

% SYNTAX:
%   [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%   Eph_R, Eph_M, iono_R, iono_M, snr_R, snr_M, time_GPS, date] = ...
%   load_RINEX(nome_FR_oss, nome_FR_nav, nome_FM_oss, nome_FM_nav);
%
% INPUT:
%   nome_FR_oss = RINEX observation file (ROVER)
%   nome_FR_nav = RINEX navigation file (ROVER)
%   nome_FM_oss = RINEX observation file (MASTER)
%   nome_FM_nav = RINEX navigation file (MASTER)
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
%   Eph_R = matrix containing 21 ephemerides for each satellite (ROVER)
%   Eph_M = matrix containing 21 ephemerides for each satellite (MASTER)
%   iono_R = matrix containing ionosphere parameters (ROVER)
%   iono_M = matrix containing ionosphere parameters (MASTER)
%   time_GPS = GPS time of ROVER observations
%   date = date (year,month,day,hour,minute,second)
%
% DESCRIPTION:
%   Parse RINEX files (both observation and navigation) for both the ROVER
%   and the MASTER. Select epochs they have in common.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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

Eph_RR = zeros(21,32);
Eph_MR = zeros(21,32);

%parse RINEX navigation file (ROVER)
[Eph_R, iono_R] = RINEX_get_nav(nome_FR_nav);

%parse RINEX navigation file (ROVER)
% [Eph_RR] = RINEX_get_nav_GLO(nome_FR_glo);

%parse RINEX navigation file (MASTER)
[Eph_M, iono_M] = RINEX_get_nav(nome_FM_nav);

%parse RINEX navigation file (MASTER)
% [Eph_MR] = RINEX_get_nav_GLO(nome_FM_glo);

%-------------------------------------------------------------------------------

%open RINEX observation file (ROVER)
FR_oss = fopen(nome_FR_oss,'r');

%open RINEX observation file (MASTER)
FM_oss = fopen(nome_FM_oss,'r');

%-------------------------------------------------------------------------------

%parse RINEX header
[obs_typ_R,info_base_R,eof_R] = obs_type_list(nome_FR_oss);
[obs_typ_M,info_base_M,eof_M] = obs_type_list(nome_FM_oss);

%check the availability of basic data to parse the RINEX file (ROVER)
if ((info_base_R == 0) | (eof_R == 1))
    error('Basic data is missing in the ROVER RINEX header')
end;

%check the availability of basic data to parse the RINEX file (MASTER)
if ((info_base_M == 0) | (eof_M == 1))
    error('Basic data is missing in the ROVER RINEX header')
end;

%computation of the number of observation types (ROVER) e.g. L1 C1 P1...
num_type_obs_R = size(obs_typ_R,2)/2;

%computation of the number of observation types (MASTER) e.g. L1 C1 P1...
num_type_obs_M = size(obs_typ_M,2)/2;

%detect the observation column index (ROVER) P-C1 L1 P2 L2
[col_ph1_R, col_ph2_R, col_cod1_R, col_cod2_R] = obs_type_find(obs_typ_R); %#ok<NASGU,ASGLU>

%detect the observation column index (ROVER) P-C1 L1 P2 L2
[col_ph1_M, col_ph2_M, col_cod1_M, col_cod2_M] = obs_type_find(obs_typ_M); %#ok<NASGU,ASGLU>

%jump the RINEX header
RINEX_jump_hdr(FR_oss);
RINEX_jump_hdr(FM_oss);

%-------------------------------------------------------------------------------

%read data for the first epoch (ROVER)
[time_GPS_R, sat_R, sat_RS, sat_RR, date_R] = RINEX_get_epoch(FR_oss);

%number of satellites (ROVER)
num_sat_R = size(sat_R,1);
num_sat_RS = size(sat_RS,1);
num_sat_RR = size(sat_RR,1);

%read ROVER observations
[mat_oss_R, snr_R_single] = RINEX_get_obs(FR_oss, num_sat_R, num_type_obs_R, col_ph1_R);
%[mat_oss_R, snr_R_single] = RINEX_get_obs(FR_oss, num_sat_R, num_type_obs_R, col_cod1_R);
[mat_oss_RS, snr_RS_single] = RINEX_get_obs(FR_oss, num_sat_RS, num_type_obs_R, col_ph1_R); %#ok<NASGU>
[mat_oss_RR, snr_RR_single] = RINEX_get_obs(FR_oss, num_sat_RR, num_type_obs_R, col_ph1_R); %#ok<NASGU>

%-------------------------------------------------------------------------------

%read data for the first epoch (MASTER)
[time_GPS_M, sat_M, sat_MS, sat_MR, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>

%number of satellites (MASTER)
num_sat_M = size(sat_M,1);
num_sat_MS = size(sat_MS,1);
num_sat_MR = size(sat_MR,1);

%read MASTER observations
[mat_oss_M, snr_M_single] = RINEX_get_obs(FM_oss, num_sat_M, num_type_obs_M, col_ph1_M);
[mat_oss_MS, snr_MS_single] = RINEX_get_obs(FM_oss, num_sat_MS, num_type_obs_M, col_ph1_M); %#ok<NASGU>
[mat_oss_MR, snr_MR_single] = RINEX_get_obs(FM_oss, num_sat_MR, num_type_obs_M, col_ph1_M); %#ok<NASGU>

%-------------------------------------------------------------------------------

while (time_GPS_M < time_GPS_R)

    %read data for the current epoch (MASTER)
    [time_GPS_M, sat_M, sat_MS, sat_MR, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>

    %number of satellites (MASTER)
    num_sat_M = size(sat_M,1);
    num_sat_MS = size(sat_MS,1);
    num_sat_MR = size(sat_MR,1);

    %read MASTER observations
    [mat_oss_M, snr_M_single] = RINEX_get_obs(FM_oss, num_sat_M, num_type_obs_M, col_ph1_M);
    [mat_oss_MS, snr_MS_single] = RINEX_get_obs(FM_oss, num_sat_MS, num_type_obs_M, col_ph1_M); %#ok<NASGU>
    [mat_oss_MR, snr_MR_single] = RINEX_get_obs(FM_oss, num_sat_MR, num_type_obs_M, col_ph1_M); %#ok<NASGU>
end

%-------------------------------------------------------------------------------

k = 1;
time_GPS(1,1) = time_GPS_R;
date(1,:) = date_R(1,:);

while (~feof(FR_oss))

    %variable initialization (GPS)
    pr1_R(:,k) = zeros(32,1);
    pr2_R(:,k) = zeros(32,1);
    ph1_R(:,k) = zeros(32,1);
    ph2_R(:,k) = zeros(32,1);
    pr1_M(:,k) = zeros(32,1);
    pr2_M(:,k) = zeros(32,1);
    ph1_M(:,k) = zeros(32,1);
    ph2_M(:,k) = zeros(32,1);
    snr_R(:,k) = zeros(32,1);
    snr_M(:,k) = zeros(32,1);

    %variable initialization (GLONASS)
    pr1_RR(:,k) = zeros(32,1);
    pr2_RR(:,k) = zeros(32,1);
    ph1_RR(:,k) = zeros(32,1);
    ph2_RR(:,k) = zeros(32,1);
    pr1_MR(:,k) = zeros(32,1);
    pr2_MR(:,k) = zeros(32,1);
    ph1_MR(:,k) = zeros(32,1);
    ph2_MR(:,k) = zeros(32,1);
    snr_RR(:,k) = zeros(32,1);
    snr_MR(:,k) = zeros(32,1);

    if (time_GPS_R == time_GPS(k))

        %read ROVER observations (GPS)
        pr1_R(sat_R,k) = mat_oss_R(:,col_cod1_R);
        %pr2_R(sat_R,k) = mat_oss_R(:,col_cod2_R);
        ph1_R(sat_R,k) = mat_oss_R(:,col_ph1_R);
        %ph2_R(sat_R,k) = mat_oss_R(:,col_ph2_R);
        snr_R(sat_R,k) = snr_R_single(:,1);

        %read ROVER observations (GLONASS)
%         pr1_RR(sat_R,k) = mat_oss_RR(:,col_cod1_R);
%         %pr2_RR(sat_R,k) = mat_oss_RR(:,col_cod2_R);
%         ph1_RR(sat_R,k) = mat_oss_RR(:,col_ph1_R);
%         %ph2_RR(sat_R,k) = mat_oss_RR(:,col_ph2_R);
%         snr_RR(sat_R,k) = snr_RR_single(:,1);

        %read data for the current epoch (ROVER)
        [time_GPS_R, sat_R, sat_RS, sat_RR, date_R] = RINEX_get_epoch(FR_oss);

        %number of satellites (ROVER)
        num_sat_R = size(sat_R,1);
        num_sat_RS = size(sat_RS,1);
        num_sat_RR = size(sat_RR,1);

        %read ROVER observations
        [mat_oss_R, snr_R_single] = RINEX_get_obs(FR_oss, num_sat_R, num_type_obs_R, col_ph1_R);
        %[mat_oss_R, snr_R_single] = RINEX_get_obs(FR_oss, num_sat_R, num_type_obs_R, col_cod1_R);
        [mat_oss_RS, snr_RS_single] = RINEX_get_obs(FR_oss, num_sat_RS, num_type_obs_R, col_ph1_R); %#ok<NASGU>
        [mat_oss_RR, snr_RR_single] = RINEX_get_obs(FR_oss, num_sat_RR, num_type_obs_R, col_ph1_R); %#ok<NASGU>

    end

    if (time_GPS_M == time_GPS(k))

        %read MASTER observations (GPS)
        pr1_M(sat_M,k) = mat_oss_M(:,col_cod1_M);
        %pr2_M(sat_M,k) = mat_oss_M(:,col_cod2_M);
        ph1_M(sat_M,k) = mat_oss_M(:,col_ph1_M);
        %ph2_M(sat_M,k) = mat_oss_M(:,col_ph2_M);
        snr_M(sat_M,k) = snr_M_single(:,1);

        %read MASTER observations (GLONASS)
%         pr1_MR(sat_M,k) = mat_oss_MR(:,col_cod1_M);
%         %pr2_MR(sat_M,k) = mat_oss_MR(:,col_cod2_M);
%         ph1_MR(sat_M,k) = mat_oss_MR(:,col_ph1_M);
%         %ph2_MR(sat_M,k) = mat_oss_MR(:,col_ph2_M);
%         snr_MR(sat_M,k) = snr_MR_single(:,1);

        %read data for the current epoch (MASTER)
        [time_GPS_M, sat_M, sat_MS, sat_MR, date_M] = RINEX_get_epoch(FM_oss); %#ok<NASGU>

        %number of satellites (MASTER)
        num_sat_M = size(sat_M,1);
        num_sat_MS = size(sat_MS,1);
        num_sat_MR = size(sat_MR,1);

        %read ROVER observations
        [mat_oss_M, snr_M_single] = RINEX_get_obs(FM_oss, num_sat_M, num_type_obs_M, col_ph1_M);
        [mat_oss_MS, snr_MS_single] = RINEX_get_obs(FM_oss, num_sat_MS, num_type_obs_M, col_ph1_M); %#ok<NASGU>
        [mat_oss_MR, snr_MR_single] = RINEX_get_obs(FM_oss, num_sat_MR, num_type_obs_M, col_ph1_M); %#ok<NASGU>

    end

    k = k+1;
    time_GPS(k,1) = time_GPS(k-1,1) + 1;
    date(k,:) = date_R(1,:);

end

time_GPS(end) = [];

%-------------------------------------------------------------------------------

%close RINEX files
fclose(FR_oss);
fclose(FM_oss);
