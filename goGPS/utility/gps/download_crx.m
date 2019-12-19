function [file_crx] = download_crx(gps_week, gps_time)

% SYNTAX:
%   [file_crx] = download_crx(gps_week, gps_time);
%
% INPUT:
%   gps_week = starting and ending GPS week [vector]
%   gps_time = starting and ending GPS time [vector]
%
% OUTPUT:
%   file_crx = donwloaded .CRX file names
%
% DESCRIPTION:
%   Download of .CRX files from the AIUB FTP server.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

% Pointer to the global settings:
state = Core.getCurrentSettings();

%AIUB FTP server IP address
% aiub_ip = '130.92.9.78'; % ftp.aiub.unibe.ch
aiub_ip = 'ftp.aiub.unibe.ch';

%download directory
down_dir = state.crx_dir;

%convert GPS time to time-of-week
gps_tow = weektime2tow(gps_week, gps_time);

% starting time
date_f = gps2date(gps_week(1), gps_tow(1));

% ending time
date_l = gps2date(gps_week(end), gps_tow(end));

% Check / create output folder
if not(exist(down_dir, 'dir'))
    mkdir(down_dir);
end

fprintf(['FTP connection to the AIUB server (ftp://' aiub_ip '). Please wait...'])

year  = date_f(1) : 1 : date_l(1);
file_crx = cell(length(year),1);

%connect to the CRX server
try
    ftp_server = ftp(aiub_ip);
catch
    fprintf(' connection failed.\n');
    return
end

fprintf('\n');

for y = 1 : length(year)

    %target directory
    s = '/BSWUSER52/GEN';

    cd(ftp_server, '/');
    cd(ftp_server, s);

    %target file
    s2 = ['SAT_' num2str(year(y),'%04d') '.CRX'];
    mget(ftp_server,s2,down_dir);

    %cell array with the paths to the downloaded files
    file_crx{y} = [down_dir, '/', s2];

    fprintf(['Downloaded CRX file: ' s2 '\n']);
end

close(ftp_server);

fprintf('Download complete.\n')
