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
