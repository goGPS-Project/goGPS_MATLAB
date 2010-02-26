function streamM2RINEX(fileroot, filename, week)

% SYNTAX:
%   streamM2RINEX(fileroot, filename);
%
% INPUT:
%   fileroot = input file root (master data, binary stream)
%   filename = output file name (master data, RINEX format)
%
% OUTPUT:
%
% DESCRIPTION:
%   File conversion from binary to RINEX format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
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

% global XM YM ZM

%----------------------------------------------------------------------------------------------

%MASTER stream reading
data_master_all = [];                                                %overall stream
hour = 0;                                                            %hour index (integer)
hour_str = num2str(hour,'%02d');                                     %hour index (string)
d = dir([fileroot '_master_' hour_str '.bin']);                      %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_master_' hour_str '.bin\n']);
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

%----------------------------------------------------------------------------------------------

%displaying
fprintf('Decoding master data \n');

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
    error('RTCM2.x conversion not supported!');
%     [cell_master] = decode_rtcm2(sixofeight);
else
    [cell_master] = decode_rtcm3(data_master_all);
end

%initialization (to make the writing faster)
Ncell   = size(cell_master,2);                        %number of read RTCM packets
time_M  = zeros(Ncell,1);                             %GPS time
pr1_M   = zeros(32,Ncell);                            %code observations
ph1_M   = zeros(32,Ncell);                            %phase observations
snr1_M  = zeros(32,Ncell);                            %signal-to-noise ratio
pr2_M   = zeros(32,Ncell);                            %code observations
ph2_M   = zeros(32,Ncell);                            %phase observations
snr2_M  = zeros(32,Ncell);                            %signal-to-noise ratio

i = 1;
for j = 1 : Ncell
    if (cell_master{1,j} == 1002)                 %RTCM 1002 message
        
        time_M(i)   = cell_master{2,j}(2);            %GPS time logging
        pr1_M(:,i)  = cell_master{3,j}(:,2);          %code observations logging
        ph1_M(:,i)  = cell_master{3,j}(:,3);          %phase observations logging
        snr1_M(:,i) = cell_master{3,j}(:,5);          %signal-to-noise ratio logging
        
        i = i+1;                                      %epoch counter increase

    elseif (cell_master{1,j} == 1004)                 %RTCM 1004 message
        
        time_M(i)   = cell_master{2,j}(2);            %GPS time logging
        pr1_M(:,i)  = cell_master{3,j}(:,2);          %code observations logging (L1)
        ph1_M(:,i)  = cell_master{3,j}(:,3);          %phase observations logging (L1)
        snr1_M(:,i) = cell_master{3,j}(:,5);          %signal-to-noise ratio logging (L1)
        pr2_M(:,i)  = cell_master{3,j}(:,7);          %code observations logging (L2)
        ph2_M(:,i)  = cell_master{3,j}(:,8);          %phase observations logging (L2)
        snr2_M(:,i) = cell_master{3,j}(:,10);         %signal-to-noise ratio logging (L2)
        
        i = i+1;

    end
end

%residual data erase (after initialization)
time_M(i:end)   = [];
pr1_M(:,i:end)  = [];
ph1_M(:,i:end)  = [];
snr1_M(:,i:end) = [];
pr2_M(:,i:end)  = [];
ph2_M(:,i:end)  = [];
snr2_M(:,i:end) = [];

%manage "nearly null" data
ph1_M(ph1_M < 1e-100) = 0;
ph2_M(ph2_M < 1e-100) = 0;

%date decoding
date = datevec(time_M/(3600*24) + 7*week + datenum([1980,1,6,0,0,0]));

%----------------------------------------------------------------------------------------------

%displaying
fprintf(['Writing: ' filename '\n\n']);

%create RINEX observation file
fid = fopen(filename,'wt');	

%write header
fprintf(fid,'     2.10           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE\n');
fprintf(fid,'goGPS               Geomatics Lab.                          PGM / RUN BY / DATE \n');
fprintf(fid,'Antenna marker                                              MARKER NAME         \n');
fprintf(fid,'Geomatics Lab.      Politecnico Milano                      OBSERVER / AGENCY   \n');
fprintf(fid,'                    ublox                                   REC # / TYPE / VERS \n');
fprintf(fid,'                    ANN-MS                                  ANT # / TYPE        \n');
% fprintf(fid,'%14.4f%14.4f%14.4f                  APPROX POSITION XYZ \n', XM, YM, ZM);
fprintf(fid,'        0.0000        0.0000        0.0000                  APPROX POSITION XYZ \n');
fprintf(fid,'        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N\n');
fprintf(fid,'     2     0                                                WAVELENGTH FACT L1/2\n');
fprintf(fid,'     3    C1    L1    S1                                    # / TYPES OF OBSERV \n');
fprintf(fid,'     1                                                      INTERVAL            \n');
fprintf(fid,'%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS   \n', ...
        date(1,1), date(1,2), date(1,3), date(1,4), date(1,5), date(1,6));
fprintf(fid,'                                                            END OF HEADER       \n');

%-------------------------------------------------------------------------------

%number of records
N = length(time_M);

%write data
for i = 1 : N
    sat = find(pr1_M(:,i) ~= 0);
    n = length(sat);
    fprintf(fid,' %02d %2d %2d %2d %2d %10.7f  0 %2d', ...
            date(i,1)-2000, date(i,2), date(i,3), date(i,4), date(i,5), round(date(i,6)), n);
    for j = 1 : n
        fprintf(fid,'G%02d',sat(j));
    end
    fprintf(fid,'\n');
    for j = 1 : n
        fprintf(fid,'%14.3f %1d',pr1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
        if (ph1_M(sat(j),i) > 1e-100)
            fprintf(fid,'%14.3f %1d',ph1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
        else
            fprintf(fid,'                ');
        end
        fprintf(fid,'%14.3f %1d',snr1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
        fprintf(fid,'\n');
    end
end

%-------------------------------------------------------------------------------

%close RINEX file
fclose(fid);
