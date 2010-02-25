function streamR2RINEX(fileroot, filename)

% SYNTAX:
%   streamR2RINEX(fileroot, filename);
%
% INPUT:
%   fileroot = input file root (rover data, binary stream)
%   filename = output file name (rover data, RINEX format)
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
global lambda1

%----------------------------------------------------------------------------------------------

%ROVER stream reading
data_rover_all = [];                                                 %overall stream
hour = 0;                                                            %hour index (integer)
hour_str = num2str(hour,'%02d');                                     %hour index (string)
d = dir([fileroot '_rover_' hour_str '.bin']);                       %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_rover_' hour_str '.bin\n']);
    num_bytes = d.bytes;                                             %file size (number of bytes)
    fid_rover = fopen([fileroot '_rover_' hour_str '.bin']);         %file opening
    data_rover = fread(fid_rover,num_bytes,'uint8');                 %file reading
    data_rover = dec2bin(data_rover,8);                              %conversion in binary number (N x 8bits matrix)
    data_rover = data_rover';                                        %transposed (8bits x N matrix)
    data_rover = data_rover(:)';                                     %conversion into a string (8N bits vector)
    fclose(fid_rover);                                               %file closing
    data_rover_all = [data_rover_all data_rover];                    %stream concatenation
    hour = hour+1;                                                   %hour increase
    hour_str = num2str(hour,'%02d');
    d = dir([fileroot '_rover_' hour_str '.bin']);                   %file to be read
end
clear hour hour_str d
clear data_rover fid_rover

%----------------------------------------------------------------------------------------------

%displaying
fprintf('Decoding rover data \n');

%message decoding
[cell_rover] = decode_ublox(data_rover_all);
clear data_rover_all

%initialization (to make writing faster)
Ncell  = size(cell_rover,2);                          %number of read RTCM packets
time_R = zeros(Ncell,1);                              %GPS time of week
week_R = zeros(Ncell,1);                              %GPS week
ph1_R  = zeros(32,Ncell);                             %phase observations
pr1_R  = zeros(32,Ncell);                             %code observations
dop_R  = zeros(32,Ncell);                             %doppler measurements
snr_R  = zeros(32,Ncell);                             %signal-to-noise ratio
lock_R = zeros(32,Ncell);                             %loss of lock indicator

i = 1;
for j = 1 : Ncell
    if (strcmp(cell_rover{1,j},'RXM-RAW'))            %RXM-RAW message data

        time_R(i)   = round(cell_rover{2,j}(1));
        week_R(i)   = cell_rover{2,j}(2);
        ph1_R(:,i)  = cell_rover{3,j}(:,1);
        pr1_R(:,i)  = cell_rover{3,j}(:,2);
        dop_R(:,i)  = cell_rover{3,j}(:,3);
        snr_R(:,i)  = cell_rover{3,j}(:,6);
        lock_R(:,i) = cell_rover{3,j}(:,7);

        %manage "nearly null" data
        pos = abs(ph1_R(:,i)) < 1e-100;
        ph1_R(pos,i) = 0;

        %phase adjustement
        pos = abs(ph1_R(:,i)) > 0 & abs(ph1_R(:,i)) < 1e7;
        if(sum(pos) ~= 0)
            ambig = 2^23;
            n = floor( (pr1_R(pos,i)/lambda1-ph1_R(pos,i)) / ambig + 0.5 );
            ph1_R(pos,i) = ph1_R(pos,i) + n*ambig;
        end

        i = i + 1;
    end
end
clear cell_rover
clear Ncell pos

%residual data erase (after initialization)
time_R(i:end)  = [];
week_R(i:end)  = [];
ph1_R(:,i:end) = [];
pr1_R(:,i:end) = [];
dop_R(:,i:end) = [];
snr_R(:,i:end) = [];
lock_R(:,i:end) = [];

%date decoding
date = datevec(time_R/(3600*24) + 7*week_R + datenum([1980,1,6,0,0,0]));

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
fprintf(fid,'     4    C1    L1    S1    D1                              # / TYPES OF OBSERV \n');
fprintf(fid,'     1                                                      INTERVAL            \n');
fprintf(fid,'%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS   \n', ...
        date(1,1), date(1,2), date(1,3), date(1,4), date(1,5), date(1,6));
fprintf(fid,'                                                            END OF HEADER       \n');

%-------------------------------------------------------------------------------

%number of records
N = length(time_R);

%write data
for i = 1 : N
    sat = find(pr1_R(:,i) ~= 0);
    n = length(sat);
    fprintf(fid,' %02d %2d %2d %2d %2d %10.7f  0 %2d', ...
            date(i,1)-2000, date(i,2), date(i,3), date(i,4), date(i,5), round(date(i,6)), n);
    for j = 1 : n
        fprintf(fid,'G%02d',sat(j));
    end
    fprintf(fid,'\n');
    for j = 1 : n
        fprintf(fid,'%14.3f %1d',pr1_R(sat(j),i),floor(snr_R(sat(j),i)/6));
        if (ph1_R(sat(j),i) > 1e-100)
            fprintf(fid,'%14.3f%1d%1d',ph1_R(sat(j),i),lock_R(sat(j),i),floor(snr_R(sat(j),i)/6));
        else
            fprintf(fid,'                ');
        end
        fprintf(fid,'%14.3f %1d',snr_R(sat(j),i),floor(snr_R(sat(j),i)/6));
        fprintf(fid,'%14.3f %1d',dop_R(sat(j),i),floor(snr_R(sat(j),i)/6));
        fprintf(fid,'\n');
    end
end

%-------------------------------------------------------------------------------

%close RINEX file
fclose(fid);
