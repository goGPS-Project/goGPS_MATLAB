function [sbas] = load_ems(data_dir_ems, gps_week, time_R)

% SYNTAX:
%   [sbas] = load_ems(data_dir_ems, gps_week, time_R);
%
% INPUT:
%   data_dir_ems = path to the directory containing EMS files [string]
%   gps_week = reference vector of GPS week numbers
%   time_R = reference vector of GPS time
%
% OUTPUT:
%   sbas = struct containing SBAS data
%
% DESCRIPTION:
%   Tool for loading EMS files and reading SBAS data.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
% Portions of code contributed by Antonio Herrera Olmo, 2012
% Adapted by Eugenio Realini, 2013
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

%EMS data initialization
sv_e   = [];
year   = []; 
month  = []; 
day    = [];
hour   = []; 
minute = []; 
second = []; 
MT     = []; 
msg    = [];

%output initialization
sbas = [];

%convert GPS time to time-of-week
gps_tow = weektime2tow(gps_week, time_R);

%directory containing EMS files
data_dir = dir(data_dir_ems);

%check the number of files contained in the directory
nmax = size(data_dir,1);

%file counter
n = 0;

%find files with ".ems" extension
for j = 1 : nmax
    
    %read the name of the j-th file
    ems_file_name = getfield(data_dir,{j,1},'name');
    
    %get the number of characters in the filename
    ems_fn_length = size(ems_file_name,2);
    
    %check if the filename corresponds to that expected from a standard EMS file (i.e. "h??.ems", with 0 <= '??' <= 23),
    % or if it corresponds to that expected from .ems files downloaded by goGPS (e.g. 2012_d347_h11.ems)
    if (((ems_fn_length == 7  && strcmpi(ems_file_name(1), 'h') && (str2num(ems_file_name(2:3))   >= 0 && str2num(ems_file_name(2:3))   <= 23))  || ...
         (ems_fn_length == 17 && strcmpi(ems_file_name(11),'h') && (str2num(ems_file_name(12:13)) >= 0 && str2num(ems_file_name(12:13)) <= 23))) && ...
          strcmpi(ems_file_name(ems_fn_length - 3 : ems_fn_length), '.ems'))
        
        n = n + 1;
        
        %full path to the target file
        ems_file_target  = strcat(data_dir_ems, '/', ems_file_name);
        
        %open .ems file
        fid_fd = fopen(ems_file_target,'r');
        
        %warnings
        if (fid_fd ~= -1)
            %fprintf(['Reading EMS file ', ems_file_name, '\n']);
            if (n == 1)
                fprintf(['Reading EMS files...\n']);
            end
        else
            fprintf(['WARNING: impossible to open EMS file ', ems_file_name, '\n']);
            break
        end

        %rewind pointer
        %frewind(fid_fd)
        
        %read the whole file. Each line contains:
        %PRN SV EGNOS, year, month, day, hour, minutes, seconds, MT, message
        data = textscan(fid_fd, '%n%n%n%n%n%n%n%n%64c');
        
        fclose(fid_fd);
        
        sv_e   = [sv_e;   data{1}];
        year   = [year;   data{2}];
        month  = [month;  data{3}];
        day    = [day;    data{4}];
        hour   = [hour;   data{5}];
        minute = [minute; data{6}];
        second = [second; data{7}];
        MT     = [MT;     data{8}];
        msg    = [msg;    data{9}];
    end
end

%if no .ems files are available, return
if (isempty(sv_e))
    fprintf(['No EMS files found in ' data_dir_ems ' directory.\n'])
    return
end

%computation of the GPS time in weeks and seconds of week
year = four_digit_year(year);
[GPS_wk_e, GPS_sec_wk_e] = date2gps([year, month, day, hour, minute, second]);

time_E(:,1) = GPS_wk_e;     % GPS week
time_E(:,2) = GPS_sec_wk_e; % GPS seconds

time_E_sec = time_E(:,1)*7*86400 + time_E(:,2); 

index = find(time_E_sec >= (time_R(1)-60*60) & time_E_sec <= (time_R(end)+60*60));

if isempty(index)
    fprintf(['No matching epochs found in ' data_dir_ems ' directory.\n'])
    return
else
    sv_e   = sv_e(index);
    year   = year(index);
    month  = month(index);
    day    = day(index);
    hour   = hour(index);
    minute = minute(index);
    second = second(index);
    MT     = MT(index);
    msg    = msg(index,:);
    time_E = time_E(index,:);
end


%%count lost messages
% dt = diff(time_E(:,2)); 
% tt = dt - 1;
% lost_mess = sum(tt);

clear GPS_wk_e GPS_sec_wk_e sv_e year month day hour minute second time_E_sec

%CRC check
check_mt = zeros(length(MT), 1);
for i = 1 : length(MT)
     
     [crc, parity] = ems_parity(msg(i,:));
     
     if (crc == parity)
         check_mt(i) = 1;
     end
end
%check_mt = ones(length(MT), 1);

%keep only the messages with valid CRC
mt_ok = find(check_mt);
MT = MT(mt_ok);
msg = msg(mt_ok,:);
time_E = time_E(mt_ok,:);

%sort messages based on time
t_abs = time_E(:,1)*604800 + time_E(:,2);
[~,ii] = sort(t_abs, 1, 'ascend');
time_E(:,1) = time_E(ii,1);
time_E(:,2) = time_E(ii,2);
MT = MT(ii);
msg(:,:) = msg(ii,:);

%%number of non-corrupted messages
%nn = length(mt_ok);

clear check_mt crc parity mt_ok t_abs yy ii

%load PRN masks
[iodp_mask, prn_mask] = load_prnmask(MT, msg); 

%load IGP masks
[iodi_mask, band_mask, igp_mask] = load_igpmask(MT, msg); 

%load fast corrections (or pseudorange corrections)
[prc_E, GPS_time_fc] = load_fc(iodp_mask, prn_mask, MT, msg, time_E); 

%load long term corrections (delta position SV and delta clock SV)
[dx_E, dy_E, dz_E, doffset_E, iode_E, GPS_time_ltc] = load_ltc(iodp_mask, prn_mask, MT, msg, time_E); 

%load the coordinates of IGP nodes and the corresponding vertical ionospheric correction
[igp, ivd_E, lat_igp, lon_igp, GPS_time_ic] = load_ic(iodi_mask, band_mask, igp_mask, MT, msg, time_E); 

%check that there are sufficient EMS data matching the survey timespan
[ems_data_available] = check_ems_availability(GPS_time_fc, GPS_time_ltc, GPS_time_ic, gps_week, gps_tow);

%if EMS data are not sufficient, return for standard (i.e. not SBAS-corrected) positioning
if (~ems_data_available)
    return
end

%matrix synchronization
[E_prc] = sync_ER(prc_E, GPS_time_fc, gps_week, gps_tow);
[E_dx]  = sync_ER(dx_E, GPS_time_ltc, gps_week, gps_tow);
[E_dy]  = sync_ER(dy_E, GPS_time_ltc, gps_week, gps_tow);
[E_dz]  = sync_ER(dz_E, GPS_time_ltc, gps_week, gps_tow);
[E_doffset] = sync_ER(doffset_E, GPS_time_ltc, gps_week, gps_tow);
[E_iode] = sync_ER(iode_E, GPS_time_ltc, gps_week, gps_tow);
[E_ivd]  = sync_ER(ivd_E, GPS_time_ic, gps_week, gps_tow);

%collect all the SBAS data in a structure
sbas.prc = E_prc;
sbas.dx  = E_dx;
sbas.dy  = E_dy;
sbas.dz  = E_dz;
sbas.doffset = E_doffset;
sbas.iode = E_iode;
sbas.ivd  = E_ivd;
sbas.igp  = igp;
sbas.lat_igp = lat_igp;
sbas.lon_igp = lon_igp;
