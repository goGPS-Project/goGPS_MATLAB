function [file_tro] = download_tro(week_start, week_end, down_dir)
% SYNTAX:
%   [file_tro] = download_tro(week_start, week_end, down_dir);
%
% INPUT:
%   week_start = starting GPS week
%   week_end   = ending GPS week
%   down_dir = download directory
%
% OUTPUT:
%   file_tro = donwloaded .tro file names
%
% DESCRIPTION:
%   Download of .tro files from the EUREF FTP server.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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


file_tro = {};

%EUREF FTP server IP address
euref_ip = '141.74.33.23'; % igs.bkg.bund.de


fprintf(['FTP connection to the EUREF server (ftp://' euref_ip '). Please wait...\n'])

%connect to the EUREF server
try
    ftp_server = ftp(euref_ip);
catch
    fprintf(' connection failed.\n');
    return
end

s = '/EUREF/products/';

for week = week_start : week_end
    for day = 0 : 6
        cd(ftp_server, '/');
        cd(ftp_server, s);

        cd(ftp_server, num2str(week,'%04d'));
        s2 = ['asi' num2str(week,'%04d') num2str(day) '.tro.Z'];
        try
            mget(ftp_server,s2,down_dir);
            fprintf(['Downloaded TRO file: ' s2 '\n']);
        catch
            fprintf(['Problem downloading TRO file: ' s2 '\n']);
        end
    end
end

close(ftp_server);
