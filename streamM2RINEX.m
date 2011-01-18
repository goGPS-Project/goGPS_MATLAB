function streamM2RINEX(fileroot, filename, week, wait_dlg)

% SYNTAX:
%   streamM2RINEX(fileroot, filename, week, wait_dlg);
%
% INPUT:
%   fileroot = input file root (master data, binary stream)
%   filename = output file name (master data, RINEX format)
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%
% DESCRIPTION:
%   File conversion from master stream (RTCM 3.x) to RINEX format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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

if (nargin == 4)
    waitbar(0.5,wait_dlg,'Reading master stream files...')
end

%MASTER stream reading
data_master_all = [];                                                %overall stream
hour = 0;                                                            %hour index (integer)
hour_str = num2str(hour,'%02d');                                     %hour index (string)
d = dir([fileroot '_master_' hour_str '.bin']);                      %file to be read
while ~isempty(d)
    if (nargin == 3)
        fprintf(['Reading: ' fileroot '_master_' hour_str '.bin\n']);
    end
    num_bytes = d.bytes;                                             %file size (number of bytes)
    fid_master = fopen([fileroot '_master_' hour_str '.bin']);       %file opening
    data_master = fread(fid_master,num_bytes,'uint8');               %file reading
    data_master = dec2bin(data_master,8);                            %conversion into binary number (N x 8bits matrix)
    data_master = data_master';                                      %transposed (8bits x N matrix)
    data_master = data_master(:)';                                   %conversion into a string (8N bits vector)
    fclose(fid_master);                                              %file closing
    data_master_all = [data_master_all data_master];                 %stream concatenation
    hour = hour+1;                                                   %hour increase
    hour_str = num2str(hour,'%02d');
    d = dir([fileroot '_master_' hour_str '.bin']);                  %file to be read
end
clear hour hour_str d
clear data_master fid_master

if (nargin == 4)
    waitbar(1,wait_dlg)
end

%----------------------------------------------------------------------------------------------

if (~isempty(data_master_all))

    %displaying
    if (nargin == 3)
        fprintf('Decoding master data \n');
    end

    pos = 1;
    sixofeight = [];
    is_rtcm2 = 1;

    while (pos + 7 <= length(data_master_all))
        if (~strcmp(data_master_all(pos:pos+1),'01'))
            is_rtcm2 = 0;
            break
        end
        sixofeight = [sixofeight fliplr(data_master_all(pos+2:pos+7))];
        pos = pos + 8;
    end
    
    %stream decodification
    if(is_rtcm2)
        error('RTCM2.x conversion not supported yet!');
        %     [cell_master] = decode_rtcm2(sixofeight);
    else
        if (nargin == 4)
            [cell_master] = decode_rtcm3(data_master_all, wait_dlg);
        else
            [cell_master] = decode_rtcm3(data_master_all);
        end
    end
    clear data_master_all
    
    %initialization (to make the writing faster)
    Ncell   = size(cell_master,2);                        %number of read RTCM packets
    time_M  = zeros(Ncell,1);                             %GPS time
    pr1_M   = zeros(32,Ncell);                            %code observations
    ph1_M   = zeros(32,Ncell);                            %phase observations
    snr1_M  = zeros(32,Ncell);                            %signal-to-noise ratio
    pr2_M   = zeros(32,Ncell);                            %code observations
    ph2_M   = zeros(32,Ncell);                            %phase observations
    snr2_M  = zeros(32,Ncell);                            %signal-to-noise ratio
    Eph_M = zeros(29,32,Ncell);                           %ephemerides
    pos_M = zeros(3,1);                                   %master station position
    
    if (nargin == 4)
        waitbar(0,wait_dlg,'Reading master data...')
    end
    
    flag_L2 = 0;
    
    i = 1;
    for j = 1 : Ncell
        if (nargin == 4)
            waitbar(j/Ncell,wait_dlg)
        end
        
        if (cell_master{1,j} == 1002)                     %RTCM 1002 message
            
            time_M(i)   = cell_master{2,j}(2);            %GPS time logging
            pr1_M(:,i)  = cell_master{3,j}(:,2);          %code observations logging
            ph1_M(:,i)  = cell_master{3,j}(:,3);          %phase observations logging
            snr1_M(:,i) = cell_master{3,j}(:,5);          %signal-to-noise ratio logging
            
            i = i+1;                                      %epoch counter increase
            
            if (cell_master{3,j}(:,1) == 0)
                code_type = 'C1';
            else
                code_type = 'P1';
            end
            
        elseif (cell_master{1,j} == 1004)                 %RTCM 1004 message
            
            time_M(i)   = cell_master{2,j}(2);            %GPS time logging
            pr1_M(:,i)  = cell_master{3,j}(:,2);          %code observations logging (L1)
            ph1_M(:,i)  = cell_master{3,j}(:,3);          %phase observations logging (L1)
            snr1_M(:,i) = cell_master{3,j}(:,5);          %signal-to-noise ratio logging (L1)
            pr2_M(:,i)  = cell_master{3,j}(:,7);          %code observations logging (L2)
            ph2_M(:,i)  = cell_master{3,j}(:,8);          %phase observations logging (L2)
            snr2_M(:,i) = cell_master{3,j}(:,10);         %signal-to-noise ratio logging (L2)
            
            flag_L2 = 1;
            
            i = i+1;
            
            if (cell_master{3,j}(:,1) == 0)
                code_type = 'C1';
            else
                code_type = 'P1';
            end

        elseif ((cell_master{1,j} == 1005) | (cell_master{1,j} == 1006)) & (pos_M == 0)
                
                coordX_M = cell_master{2,j}(8);
                coordY_M = cell_master{2,j}(9);
                coordZ_M = cell_master{2,j}(10);
                
                pos_M(:,1) = [coordX_M; coordY_M; coordZ_M];
            
        elseif (cell_master{1,j} == 1019)                 %RTCM 1019 message
            
            %satellite number
            sat = cell_master{2,j}(1);                    %satellite number
            tom = cell_master{2,j}(21);                   %time of measurement
            
            %if the ephemerides are not already available
            if (isempty(find(Eph_M(21,sat,:) ==  tom, 1)))
                Eph_M(:,sat,i) = cell_master{2,j}(:);     %single satellite ephemerides logging
            end
        end
    end
    clear Ncell pos sat tom
    
    %residual data erase (after initialization)
    time_M(i:end)   = [];
    pr1_M(:,i:end)  = [];
    ph1_M(:,i:end)  = [];
    snr1_M(:,i:end) = [];
    pr2_M(:,i:end)  = [];
    ph2_M(:,i:end)  = [];
    snr2_M(:,i:end) = [];
    Eph_M(:,:,i:end) = [];
    
    %manage "nearly null" data
    ph1_M(ph1_M < 1e-100) = 0;
    ph2_M(ph2_M < 1e-100) = 0;

    %date decoding
    date = datevec(time_M/(3600*24) + 7*week + datenum([1980,1,6,0,0,0]));

    %----------------------------------------------------------------------------------------------
    % RINEX OBSERVATION FILE
    %----------------------------------------------------------------------------------------------

    %displaying
    if (nargin == 3)
        fprintf(['Writing: ' filename '.obs\n']);
    end

    %create RINEX observation file
    fid_obs = fopen([filename '.obs'],'wt');
    
    %write header
    fprintf(fid_obs,'     2.10           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE\n');
    fprintf(fid_obs,'goGPS                                                       PGM / RUN BY / DATE \n');
    fprintf(fid_obs,'                                                            MARKER NAME         \n');
    fprintf(fid_obs,'                                                            OBSERVER / AGENCY   \n');
    fprintf(fid_obs,'                                                            REC # / TYPE / VERS \n');
    fprintf(fid_obs,'                                                            ANT # / TYPE        \n');
    fprintf(fid_obs,'%14.4f%14.4f%14.4f                  APPROX POSITION XYZ \n', pos_M(1,1), pos_M(2,1), pos_M(3,1));
    fprintf(fid_obs,'        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N\n');
    fprintf(fid_obs,'     1     1                                                WAVELENGTH FACT L1/2\n');
    if (flag_L2)
        fprintf(fid_obs,['     6    ' code_type '    P2    L1    L2    S1    S2                  # / TYPES OF OBSERV \n']);
    else
        fprintf(fid_obs,['     3    ' code_type '    L1    S1                                    # / TYPES OF OBSERV \n']);
    end
    fprintf(fid_obs,'     1                                                      INTERVAL            \n');
    fprintf(fid_obs,'%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS   \n', ...
        date(1,1), date(1,2), date(1,3), date(1,4), date(1,5), date(1,6));
    fprintf(fid_obs,'                                                            END OF HEADER       \n');
    
    clear cell_master
    
    %-------------------------------------------------------------------------------
    
    %number of records
    N = length(time_M);
    
    if (nargin == 4)
        waitbar(0,wait_dlg,'Writing master observation file...')
    end
    
    %write data
    for i = 1 : N
        if (nargin == 4)
            waitbar(i/N,wait_dlg)
        end
        
        sat = find(pr1_M(:,i) ~= 0);
        n = length(sat);
        
        %if no observations are available, do not write anything
        if (n > 0)
            fprintf(fid_obs,' %02d %2d %2d %2d %2d %10.7f  0 %2d', ...
                date(i,1)-2000, date(i,2), date(i,3), date(i,4), date(i,5), date(i,6), n);
            if (n>12)
                for j = 1 : 12
                    fprintf(fid_obs,'G%02d',sat(j));
                end
                fprintf(fid_obs,'\n');
                fprintf(fid_obs,'%32s','');
                for j = 13 : n
                    fprintf(fid_obs,'G%02d',sat(j));
                end
            else
                for j = 1 : n
                    fprintf(fid_obs,'G%02d',sat(j));
                end
            end
            fprintf(fid_obs,'\n');
            for j = 1 : n
                fprintf(fid_obs,'%14.3f %1d',pr1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                if (flag_L2)
                    fprintf(fid_obs,'%14.3f %1d',pr2_M(sat(j),i),floor(snr2_M(sat(j),i)/6));
                end
                if (abs(ph1_M(sat(j),i)) > 1e-100)
                    fprintf(fid_obs,'%14.3f %1d',ph1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                else
                    fprintf(fid_obs,'                ');
                end
                if (flag_L2)
                    if (abs(ph2_M(sat(j),i)) > 1e-100)
                        fprintf(fid_obs,'%14.3f %1d',ph2_M(sat(j),i),floor(snr2_M(sat(j),i)/6));
                    else
                        fprintf(fid_obs,'                ');
                    end
                end
                fprintf(fid_obs,'%14.3f %1d',snr1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                if (flag_L2)
                    fprintf(fid_obs,'\n');
                end
                if (flag_L2)
                    fprintf(fid_obs,'%14.3f %1d',snr2_M(sat(j),i),floor(snr2_M(sat(j),i)/6));
                end
                fprintf(fid_obs,'\n');
            end
        end
    end
    
    %close RINEX observation file
    fclose(fid_obs);
    
    %----------------------------------------------------------------------------------------------
    % RINEX NAVIGATION FILE
    %----------------------------------------------------------------------------------------------
    
    %if ephemerides are available
    if (~isempty(find(Eph_M(1,:,:) ~= 0, 1)))
        
        %displaying
        if (nargin == 3)
            fprintf(['Writing: ' filename '.nav\n']);
        end
        
        %create RINEX observation file
        fid_nav = fopen([filename '.nav'],'wt');
        
        %write header
        fprintf(fid_nav,'     2.10           NAVIGATION DATA                         RINEX VERSION / TYPE\n');
        fprintf(fid_nav,'goGPS                                                       PGM / RUN BY / DATE \n');
        fprintf(fid_nav,'                                                            END OF HEADER       \n');
        
        if (nargin == 4)
            waitbar(0,wait_dlg,'Writing master navigation file...')
        end
        
        for i = 1 : N
            if (nargin == 4)
                waitbar(i/N,wait_dlg)
            end
            
            satEph = find(Eph_M(1,:,i) ~= 0);
            for j = 1 : length(satEph)
                af2      = Eph_M(2,satEph(j),i);
                M0       = Eph_M(3,satEph(j),i);
                roota    = Eph_M(4,satEph(j),i);
                deltan   = Eph_M(5,satEph(j),i);
                ecc      = Eph_M(6,satEph(j),i);
                omega    = Eph_M(7,satEph(j),i);
                cuc      = Eph_M(8,satEph(j),i);
                cus      = Eph_M(9,satEph(j),i);
                crc      = Eph_M(10,satEph(j),i);
                crs      = Eph_M(11,satEph(j),i);
                i0       = Eph_M(12,satEph(j),i);
                idot     = Eph_M(13,satEph(j),i);
                cic      = Eph_M(14,satEph(j),i);
                cis      = Eph_M(15,satEph(j),i);
                Omega0   = Eph_M(16,satEph(j),i);
                Omegadot = Eph_M(17,satEph(j),i);
                toe      = Eph_M(18,satEph(j),i);
                af0      = Eph_M(19,satEph(j),i);
                af1      = Eph_M(20,satEph(j),i);
                tom      = Eph_M(21,satEph(j),i);
                IODE     = Eph_M(22,satEph(j),i);
                codes    = Eph_M(23,satEph(j),i);
                weekno   = Eph_M(24,satEph(j),i); %#ok<NASGU>
                L2flag   = Eph_M(25,satEph(j),i);
                svaccur  = Eph_M(26,satEph(j),i);
                svhealth = Eph_M(27,satEph(j),i);
                tgd      = Eph_M(28,satEph(j),i);
                fit_int  = Eph_M(29,satEph(j),i);
                
                %time of measurement decoding
                date = datevec(tom/(3600*24) + 7*week + datenum([1980,1,6,0,0,0]));
                
                lineE(1,:) = sprintf('%2d %02d %2d %2d %2d %2d%5.1f% 18.12E% 18.12E% 18.12E\n', ...
                    satEph(j),date(1)-2000, date(2), date(3), date(4), date(5), date(6), ...
                    af0, af1, af2);
                linesE(1,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', IODE , crs, deltan, M0);
                linesE(2,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', cuc, ecc, cus, roota);
                linesE(3,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', toe, cic, Omega0, cis);
                linesE(4,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', i0, crc, omega, Omegadot);
                linesE(5,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', idot, codes, week, L2flag);
                linesE(6,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', svaccur, svhealth, tgd, IODE);
                linesE(7,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', tom, fit_int, 0, 0);
                
                %if running on Windows, convert three-digits exponential notation
                %to two-digits; in any case, replace 'E' with 'D' and print the string
                if (~isunix)
                    lineD = strrep(lineE(1,:),'E+0','D+');
                    lineD = strrep(lineD,'E-0','D-');
                    fprintf(fid_nav,'%s',lineD);
                    for k = 1 : 7
                        lineD = strrep(linesE(k,:),'E+0','D+');
                        lineD = strrep(lineD,'E-0','D-');
                        fprintf(fid_nav,'%s',lineD);
                    end
                else
                    lineD = strrep(lineE(1,:),'E+','D+');
                    lineD = strrep(lineD,'E-','D-');
                    fprintf(fid_nav,'%s',lineD);
                    for k = 1 : 7
                        lineD = strrep(linesE(k,:),'E+','D+');
                        lineD = strrep(lineD,'E-','D-');
                        fprintf(fid_nav,'%s',lineD);
                    end
                end
            end
        end
        
        %close RINEX navigation file
        fclose(fid_nav);
    end
    
else
    %displaying
    if (nargin == 4)
        msgbox('No master data acquired.');
    else
        fprintf('No master data acquired! \n');
    end
end
