function [file_dcb, compressed] = download_dcb(gps_week, gps_time)

% SYNTAX:
%   [file_dcb, compressed] = download_dcb(gps_week, gps_time);
%
% INPUT:
%   gps_week = starting and ending GPS week [vector]
%   gps_time = starting and ending GPS time [vector]
%
% OUTPUT:
%   file_dcb = donwloaded .DCB file names 
%   compressed = flag to let the calling function know whether the
%                downloaded files are still compressed
%
% DESCRIPTION:
%   Download of .DCB files from the AIUB FTP server.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3 beta
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

file_dcb = {};
compressed = 0;

%AIUB FTP server IP address
aiub_ip = '130.92.9.78'; % ftp.unibe.ch

%download directory
down_dir = '../data/DCB';

%convert GPS time to time-of-week
gps_tow = weektime2tow(gps_week, gps_time);

% starting time
date_f = gps2date(gps_week(1), gps_tow(1));

% ending time
date_l = gps2date(gps_week(end), gps_tow(end));

fprintf(['FTP connection to the AIUB server (http://' aiub_ip '). Please wait...'])

year_orig  = date_f(1) : 1 : date_l(1);
if (length(year_orig) < 1)
    %fprintf('ERROR: Data range not valid.\n')
    return
elseif (length(year_orig) == 1)
    month = date_f(2) : 1 : date_l(2);
    year = year_orig;
else
    month = date_f(2) : 1 : 12;
    year  = date_f(1).*ones(size(month));
    for y = 2 : length(year_orig)-1
        month = [month 1 : 1 : 12]; 
        year = [year (year+y-1).*ones(1,12)];
    end
    month = [month 1 : 1 : date_l(2)];
    year  = [year date_l(1).*ones(1,date_l(2))];
end

%connect to the DCB server
try
    ftp_server = ftp(aiub_ip);
catch
    fprintf(' connection failed.\n');
    return
end

fprintf('\n');

m = 0;

for y = 1 : length(year_orig)
    
    %target directory
    s = ['/aiub/CODE/', num2str(year_orig(y))];
    
    cd(ftp_server, '/');
    cd(ftp_server, s);
    
    while(m <= length(month)-1)
        
        m = m + 1;
        
        ff = {'P1C1','P1P2'};
        
        for p = 1 : length(ff)
            %target file
            s2 = [ff{p} num2str(two_digit_year(year(m)),'%02d') num2str(month(m),'%02d') '.DCB.Z'];
            mget(ftp_server,s2,down_dir);
            if (isunix())
                system(['uncompress -f ' down_dir '/' s2]);
            end
            
            %cell array with the paths to the downloaded files
            entry = {[down_dir, '/', s2]};
            file_dcb = [file_dcb; entry];
            
            fprintf(['Downloaded DCB file: ' s2 '\n']);
            
            if (~isunix())
                fprintf(['Please decompress the ' s2 ' file before trying to use it in goGPS.\n']);
                compressed = compressed + 1;
            end
        end
        
        if (month(m) == 12)
            break
        end
    end
end

close(ftp_server);

fprintf('Download complete.\n')
