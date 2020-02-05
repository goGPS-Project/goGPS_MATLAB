function [file_met] = download_met(doy_start, doy_end, down_dir)

file_met = {};

%EUREF FTP server IP address
igs_ip = 'igs.bkg.bund.de'; %  www.epncb.oma.be


fprintf(['FTP connection to the EPN server (ftp://' igs_ip '). Please wait...\n'])

%connect to the EPN server
try
    ftp_server = ftp(igs_ip);
catch
    fprintf(' connection failed.\n');
    return
end

s = 'pub/obs/2011';

for doy = doy_start : doy_end
    
        cd(ftp_server, '/');
        cd(ftp_server, s);
        
        cd(ftp_server, num2str(doy,'%03d'));
        s2 = ['M0SE' num2str(doy,'%03d') num2str(0) '.11D.Z'];
        try
            mget(ftp_server,s2,down_dir);
            fprintf(['Downloaded OBS file: ' s2 '\n']);
        catch
            fprintf(['Problem downloading OBS file: ' s2 '\n']);
        end
   
end