function [file_ems] = download_ems(prn, gps_week, gps_time)

% SYNTAX:
%   [file_ems] = download_ems(prn, gps_week, gps_time);
%
% INPUT:
%   prn      = EGNOS satellite PRN (120,124 or 126)
%   gps_week = starting and ending GPS week [vector]
%   gps_time = starting and ending GPS time [vector]
%
% OUTPUT:
%   file_ems = donwloaded .ems file names 
%
% DESCRIPTION:
%   Download of .ems files from the EGNOS Message Server, via Internet.
%
%   WARNING: .ems files are available with about 1 hour latency

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
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

file_ems = {};

if (prn ~= 120 & prn ~= 124 & prn ~= 126)
    %fprintf('ERROR: Valid EGNOS PRNs are only 120, 124 and 126.\n')
    return
end

%EGNOS Message Server IP address
ems_ip = '131.176.49.48'; % ems.estec.esa.int

%download directory
down_dir = '../data/EMS';

%buffer in minutes
buf_min = 13;

%convert GPS time to time-of-week
gps_tow = weektime2tow(gps_week, gps_time);

% starting time
[date_f, day_of_year_f] = gps2date(gps_week(1), gps_tow(1));

% ending time
[date_l, day_of_year_l] = gps2date(gps_week(end), gps_tow(end));

%store the initial values in case the year changes
date_f_temp = date_f;
date_l_temp = date_l;

%check leap years
if ((mod(date_f(1),4) == 0) && ((mod(date_f(1),100) ~= 0) || (mod(date_f(1),400) == 0)))
    doy_max_f = 366;
else
    doy_max_f = 365;
end
if ((mod(date_l(1),4) == 0) && ((mod(date_l(1),100) ~= 0) || (mod(date_l(1),400) == 0)))
    doy_max_l = 366;
else
    doy_max_l = 365;
end
%check if the previous year was a leap year
if ((mod(date_f(1)-1,4) == 0) && ((mod(date_f(1)-1,100) ~= 0) || (mod(date_f(1)-1,400) == 0)))
    doy_max_p = 366;
else
    doy_max_p = 365;
end

flag_year = 0;

%if the requested timespan starts within the first 13 minutes of the hour,
% download also the file of the previous hour
if (date_f(5) < buf_min)
    
    %hour update
    date_f(4) = date_f(4) - 1;
    
    %check if the day has changed
    if (date_f(4) < 0)
        
        date_f(4) = 23;
        
        %update the day
        day_of_year_f = day_of_year_f - 1;
        
        %check if the year has changed
        if (day_of_year_f < 1)
            
            %restore initial values
            date_f = date_f_temp;
            
            %update year, month, day, hour, minutes
            date_f(1) = date_f(1) - 1;
            date_f(2) = 12;
            date_f(3) = 31;
            date_f(4) = 23;
            date_f(5) = 60 - (buf_min - date_f(5));
            
            %update day of year
            day_of_year_f = doy_max_p;
            
            flag_year = 1;
        end
    end
end

%if the requested timespan ends within the last 13 minutes of the hour,
% download also the file of the next hour
if (date_l(5) > 60-buf_min)
    
    %hour update
    date_l(4) = date_l(4) + 1;
    
    %check if the day has changed
    if (date_l(4) > 23)
        
        date_l(4) = 0;
        
        %update the day
        day_of_year_l = day_of_year_l + 1;
        
        %check if the year has changed
        if (day_of_year_l > doy_max_l)
            
            %restore initial values
            date_l = date_l_temp;
            
            %update year, month, day, hour, minutes
            date_l(1) = date_l(1) + 1;
            date_l(2) = 1;
            date_l(3) = 1;
            date_l(4) = 0;
            date_l(5) = buf_min - (60 - date_l(5));
            
            %update day of year
            day_of_year_l = 1;
            
            flag_year = 2;
        end
    end
end

fprintf(['FTP connection to the EGNOS Message Server (http://' ems_ip '). Please wait...'])

if (day_of_year_f <= day_of_year_l) %within the same year
    
    doy   = [day_of_year_f : 1 : day_of_year_l];
    years = ones(1,length(doy))*date_f(1);
    
elseif (flag_year == 1) %including one day in the previous year
    
    doy = [1 : day_of_year_l];
    doy = [doy_max_p doy];
    years = [date_f(1), ones(1,length(doy)-1)*date_l(1)];
    
elseif (flag_year == 2) %including one day in the following year
    
    doy   = [day_of_year_f : doy_max_f];
    doy   = [doy 1];
    years = [ones(1,length(doy)-1)*date_f(1), date_l(1)];
end

if (length(doy) < 1)
    %fprintf('ERROR: Data range not valid.\n')
    return
end

%connect to the EMS server
try
    ftp_server = ftp(ems_ip);
catch
    fprintf(' connection failed.\n');
    return
end

fprintf('\n');

for i = 1 : length(doy)
    
    %target directory
    s = ['pub/PRN', num2str(prn), '/y', num2str(years(i)), '/d', num2str(doy(i),'%03d')];
    
    if (length(doy) == 1) %only one day
        hours = [date_f(4) : 1 : date_l(4)];
    else %multiple days
        if (i == 1) %start day
            hours = [date_f(4) : 1 : 23];
        elseif (i == length(doy)) %end day
            hours = [0 : 1 : date_l(4)];
        else %days within start and end days
            hours = [0 : 1 : 23];
        end
    end

    cd(ftp_server, '/');
    cd(ftp_server, s);
    
    for j = 1 : length(hours)
        
        %target file
        s2 = ['h', num2str(hours(j),'%02d'), '.ems'];
        mget(ftp_server,s2,down_dir);
        
        %add year and day-of-year to the downloaded file
        filename = [num2str(years(i)), '_d', num2str(doy(i),'%03d'), '_', s2];
        movefile([down_dir, '/', s2], [down_dir, '/', filename]);

        %cell array with the paths to the downloaded files
        entry = {[down_dir, '/', filename]};
        file_ems = [file_ems; entry];
        
        %if the downloaded file is empty, stop here and let the caller switch to the next PRN
        s3 = dir([down_dir, '/', filename]);
        if (s3.bytes == 0)
            file_ems = {};
            return
        end
        
        fprintf(['Downloaded EGNOS .ems file: PRN ', num2str(prn), ' Year ', num2str(years(i)), ' DOY ', num2str(doy(i),'%03d'), ' Hour ', num2str(hours(j),'%02d'), '\n'])
    end
end

close(ftp_server);

fprintf('Download complete.\n')
