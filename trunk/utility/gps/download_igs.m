function download_igs(gpsweek, dow, items, path_download)

% SYNTAX:
%   download_igs(gpsweek, dow, items, path_download);
%
% INPUT:
%   gpsweek = 4digit gpsweek
%   dow = day of week
%   items = vector specifying IGS products to be downloaded 
%   path_download = full path to the download directory
%
% OUTPUT:
%
% DESCRIPTION:
%   Download and decompress IGS products as specified in items for day
%   corresponding to gpsweek and dow.
%	items is a vector filled with 1 (download) - 0 (don't download) as
%	follows:
%   items =  [igswwwwd.sp3              IGS final orbits           
%             igrwwwwd.sp3              IGS rapid orbtis           
%             iguwwwwd_00.sp3           IGS ultrarapid orbits 00   
%             iguwwwwd_06.sp3           IGS ultrarapid orbits 06   
%             iguwwwwd_12.sp3           IGS ultrarapid orbits 12   
%             iguwwwwd_18.sp3]          IGS ultrarapid orbits 18
%   e.g. items=[1 0 0 0 0 0] to download final orbits only

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
% Routine by Stefano Caldera, 2012-05-10
%----------------------------------------------------------------------------------------------

addpath(genpath(pwd));

fprintf('\n');
fprintf('--> IGSDownload <--\n');
fprintf('    -----------\n\n');
    

% checking if the download is needed
if sum(items)==0
    % nothing do download
    fprintf('    --> Warning: no file requested.\n\n');
else
    % ftp address and credential
    username='anonymous';
    password='stefano@geomatica.como.polimi.it';
    
    
    items_filename={};
    items_filename{1}=sprintf('igs%04d%d.sp3.Z',gpsweek,dow);
    items_filename{2}=sprintf('igr%04d%d.sp3.Z',gpsweek,dow);
    items_filename{3}=sprintf('igu%04d%d_00.sp3.Z',gpsweek,dow);
    items_filename{4}=sprintf('igu%04d%d_06.sp3.Z',gpsweek,dow);
    items_filename{5}=sprintf('igu%04d%d_12.sp3.Z',gpsweek,dow);
    items_filename{6}=sprintf('igu%04d%d_18.sp3.Z',gpsweek,dow);    

    items_status=items-1;   % -1: not to be donwloaded, 0: not found, 1: downloaded
    
    
    %              ================================
    %  SERVER # 1  International GNSS Service (IGS)
    %              ================================
    
    fprintf('    --> Donwload from International GNSS Service (IGS)\n');
    fprintf('        --> Connecting to remote host ... ');
    
    host='igscb.jpl.nasa.gov';
    
    try
        f = ftp(host,username,password);
        fprintf('OK\n');
        
        fprintf('        --> Changing directory ... ');
        cd(f,'pub');
        cd(f,'gps');
        cd(f,sprintf('%04d',gpsweek));
        fprintf('OK\n');
        fprintf('        --> Download:\n');
        
        for i=1:length(items)
            if items(i)==1
                fprintf('              - %-17s ... ', char(items_filename{i}));
                mget(f,char(items_filename{i}),path_download);
                
                % check download result
                if exist(sprintf('%s/%s',path_download, char(items_filename{i})))
                    items(i)=0;
                    items_status(i)=1;
                    fprintf('OK\n');
                else
                    items_status(i)=0;
                    fprintf('NOT FOUND\n');
                end
            end
        end
        close(f);
        
    catch connection_status
        fprintf('%s\n',connection_status.message);        
    end

    if sum(items)>0
        %              ===================================
        %  SERVER # 2  Goddard Space Flight Center (CCDIS)
        %              ===================================
        fprintf('\n    --> Donwload from Goddard Space Flight Center (CCDIS)\n');
        fprintf('        --> Connecting to remote host ... ');
        
        host='cddis.gsfc.nasa.gov';
        
        try
            f = ftp(host,username,password);
            fprintf('OK\n');
            
            fprintf('        --> Changing directory ... ');
            cd(f,'pub');
            cd(f,'gps');
            cd(f,'products');
            cd(f,sprintf('%04d',gpsweek));
            fprintf('OK\n');
            fprintf('        --> Download:\n');
            for i=1:length(items)
                if items(i)==1
                    fprintf('              - %-17s ... ', char(items_filename{i}));
                    mget(f,char(items_filename{i}),path_download);
                    
                    % check download result
                    if exist(sprintf('%s/%s',path_download, char(items_filename{i})))
                        items(i)=0;
                        items_status(i)=1;
                        fprintf('OK\n');
                    else
                        items_status(i)=0;
                        fprintf('NOT FOUND\n');
                    end
                end
            end
            close(f);

        catch connection_status
            fprintf('%s\n',connection_status.message);
        end
    end

    if sum(items)>0
        %              ================================================
        %  SERVER # 2  Scripps Orbit and Permanent Array Center (SOPAC)
        %              ================================================
        fprintf('\n    --> Donwload from Scripps Orbit and Permanent Array Center (SOPAC)\n');
        fprintf('        --> Connecting to remote host ... ');
        
        host='garner.ucsd.edu';
        
        try
            f = ftp(host,username,password);
            fprintf('OK\n');
            
            fprintf('        --> Changing directory ... ');
            cd(f,'pub');
            cd(f,'products');
            cd(f,sprintf('%04d',gpsweek));
            fprintf('OK\n');
            fprintf('        --> Download:\n');
            
            for i=1:length(items)
                if items(i)==1
                    fprintf('              - %-17s ... ', char(items_filename{i}));
                    mget(f,char(items_filename{i}),path_download);
                    
                    % check download result
                    if exist(sprintf('%s/%s',path_download, char(items_filename{i})))
                        items(i)=0;
                        items_status(i)=1;
                        fprintf('OK\n');
                    else
                        items_status(i)=0;
                        fprintf('NOT FOUND\n');
                    end
                end
            end
        catch connection_status
            fprintf('%s\n',connection_status.message);
        end
        close(f);
    end

%      % decompression
%      fprintf('\n    --> Decompression\n');
%      
%      if(isunix)
%          % decompression using unix
%          
%          
%      else
%          % decompression in windows
%          for i=1:length(items)
%              if items_status(i)==1
%                  name_char=char(items_filename{i});
%                  fprintf('        - %17s --> %17s ... ', name_char, name_char(1:end-2));
%                  command=sprintf('gzip -d %s/%s',path_download, name_char);
%                  dos(command);
%                  fprintf('OK\n');
%              end
%          end
%      end
end
