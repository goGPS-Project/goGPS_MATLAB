function [time_GPS, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, snr_R, snr_M, Eph_M, ...
          loss_R, loss_M, data_rover_all, data_master_all] = load_stream (fileroot)

% SYNTAX:
%   [time_GPS, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, snr_R, snr_M, Eph_M, ...
%    loss_R, loss_M, data_rover_all, data_master_all] = load_stream (fileroot);
%
% INPUT:
%   fileroot = name of the file to be read
%
% OUTPUT:
%   time_GPS = reference GPS time
%   time_R   = GPS time for the ROVER observations
%   time_M   = GPS time for the MASTER observations
%   pr1_R    = ROVER-SATELLITE code-pseudorange (carrier L1)
%   pr1_M    = MASTER-SATELLITE code-pseudorange (carrier L1)
%   ph1_R    = ROVER-SATELLITE phase observations (carrier L1)
%   ph1_M    = MASTER-SATELLITE phase observations (carrier L1)
%   Eph_M    = matrix of 21 parameters each satellite (MASTER)
%   loss_R   = flag for the ROVER loss of signal
%   loss_M   = flag for the MASTER loss of signal
%   data_rover_all  = ROVER overall stream 
%   data_master_all = MASTER overall stream 
%
% DESCRIPTION:
%   Reading of the streams trasmitted by the u-blox receiver (ROVER) and
%   by the permanent station (MASTER) in RTCM format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
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

%-------------------------------------------------------------------------------

if ~isempty(data_rover_all)

    %displaying
    fprintf(['Decoding rover data \n']);

    %message decoding
    [cell_rover] = decode_ublox(data_rover_all);

    %initialization (to make the writing faster)
    Ncell  = size(cell_rover,2);                          %number of read RTCM packets
    time_R = zeros(Ncell,1);                              %GPS time
    pr1_R  = zeros(32,Ncell);                             %code observations
    ph1_R  = zeros(32,Ncell);                             %phase observations
    snr_R  = zeros(32,Ncell);                             %signal-to-noise ratio

    i = 1;
    for j = 1 : Ncell
        if (strcmp(cell_rover{1,j},'RXM-RAW'))            %RXM-RAW message data

            time_R(i) = round(cell_rover{2,j}(1));
            pr1_R(:,i) = cell_rover{3,j}(:,2);
            ph1_R(:,i) = cell_rover{3,j}(:,1);
            snr_R(:,i) = cell_rover{3,j}(:,6);

            i = i + 1;
        end
    end
    
    %residual data erase (after initialization)
    time_R(i:end)  = [];
    pr1_R(:,i:end) = [];
    ph1_R(:,i:end) = [];

    %manage "nearly null" data
    ph1_R(find(ph1_R < 1e-100)) = 0;

else
    %displaying
    fprintf('No rover data acquired! \n');

    time_R = [];                         %GPS time
    pr1_R  = [];                         %code observations
    ph1_R  = [];                         %phase observations
    snr_R  = [];                         %signal-to-noise ratio

end

%-------------------------------------------------------------------------------

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

%-------------------------------------------------------------------------------

if ~isempty(data_master_all)

    %displaying
    fprintf(['Decoding master data \n']);

    %stream decodification
    [cell_master] = decode_rtcm(data_master_all);

    %initialization (to make the writing faster)
    Ncell  = size(cell_master,2);                         %number of read RTCM packets
    time_M = zeros(Ncell,1);                              %GPS time
    pr1_M  = zeros(32,Ncell);                             %code observations
    ph1_M  = zeros(32,Ncell);                             %phase observations
    snr_M  = zeros(32,Ncell);                             %signal-to-noise ratio
    Eph_M  = zeros(21,32,Ncell);                          %ephemerides

    i = 1;
    for j = 1 : Ncell
        if (cell_master{1,j} == 1002) | (cell_master{1,j} == 1004) %RTCM 1002 or 1004 message

            time_M(i)  = cell_master{2,j}(2);             %GPS time logging
            pr1_M(:,i) = cell_master{3,j}(:,2);           %code observations logging
            ph1_M(:,i) = cell_master{3,j}(:,3);           %phase observations logging
            snr_M(:,i) = cell_master{3,j}(:,5);           %signal-to-noise ratio logging

            i = i+1;                                      %epoch counter increase

            Eph_M(:,:,i)   = Eph_M(:,:,i-1);              %previous epoch ephemerides copying

        elseif (cell_master{1,j} == 1019)                 %RTCM 1019 message

            sat = cell_master{2,j}(1);                    %satellite number
            Eph_M(:,sat,i) = cell_master{2,j}(:);         %single satellite ephemerides logging
        end
    end

    %residual data erase (after initialization)
    time_M(i:end)  = [];
    pr1_M(:,i:end) = [];
    ph1_M(:,i:end) = [];
    snr_M(:,i:end) = [];
    Eph_M(:,:,i:end) = [];

    ph1_M(find(ph1_M < 1e-100)) = 0;

else
    %displaying
    fprintf('No master data acquired! \n');

    time_M = [];                         %GPS time
    pr1_M  = [];                         %code observations
    ph1_M  = [];                         %phase observations
    snr_M  = [];                         %signal-to-noise ratio
    Eph_M  = [];                         %ephemerides
end

%-------------------------------------------------------------------------------

if ~isempty(time_R) & ~isempty(time_M)

    %initial synchronization
    while (time_R(1) < time_M(1))
        time_R(1)  = [];                         %GPS time
        pr1_R(:,1) = [];                         %code observations
        ph1_R(:,1) = [];                         %phase observations
        snr_R(:,1) = [];                         %signal-to-noise ratio
    end

    while (time_M(1) < time_R(1))
        time_M(1)  = [];                         %GPS time
        pr1_M(:,1) = [];                         %code observations
        ph1_M(:,1) = [];                         %phase observations
        snr_M(:,1) = [];                         %signal-to-noise ratio
        Eph_M(:,1) = [];                         %ephemerides
    end

end

%-------------------------------------------------------------------------------

%signal losses
time_GPS = union(time_R,time_M);                     %overall reference time
if ~isempty(time_GPS)

    time_GPS = (time_GPS(1) : 1 : time_GPS(end))';   %GPS time without interruptions

    loss_R = 1 - ismember(time_GPS,time_R);          %losses of signal (ROVER)
    loss_M = 1 - ismember(time_GPS,time_M);          %losses of signal (MASTER)

    if ~isempty(time_R)

        newtime_R = setdiff(time_GPS, time_R);       %ROVER missing epochs
        for i = 1 : length(newtime_R)

            pos = find(time_R == newtime_R(i) - 1);  %position before the "holes"

            time_R = [time_R(1:pos);  newtime_R(i);  time_R(pos+1:end)];
            pr1_R  = [pr1_R(:,1:pos)  zeros(32,1)    pr1_R(:,pos+1:end)];
            ph1_R  = [ph1_R(:,1:pos)  zeros(32,1)    ph1_R(:,pos+1:end)];
            snr_R  = [snr_R(:,1:pos)  zeros(32,1)    snr_R(:,pos+1:end)];
        end
    else
        time_R = time_GPS;
        pr1_R  = zeros(32,length(time_GPS));
        ph1_R  = zeros(32,length(time_GPS));
        snr_R  = zeros(32,length(time_GPS));
    end

    if ~isempty(time_M)

        newtime_M = setdiff(time_GPS, time_M);       %MASTER missing epochs
        for i = 1 : length(newtime_M)

            pos = find(time_M == newtime_M(i) - 1);  %position before the "holes"

            time_M = [time_M(1:pos);  newtime_M(i);  time_M(pos+1:end)];
            pr1_M  = [pr1_M(:,1:pos)  zeros(32,1)    pr1_M(:,pos+1:end)];
            ph1_M  = [ph1_M(:,1:pos)  zeros(32,1)    ph1_M(:,pos+1:end)];
            snr_M  = [snr_M(:,1:pos)  zeros(32,1)    snr_M(:,pos+1:end)];

            Eph_M  = cat(3, Eph_M(:,:,1:pos), zeros(21,32,1), Eph_M(:,:,pos+1:end));
        end
    else
        time_M = time_GPS;
        pr1_M  = zeros(32,length(time_GPS));
        ph1_M  = zeros(32,length(time_GPS));
        snr_M  = zeros(32,length(time_GPS));
        Eph_M  = zeros(21,32,length(time_GPS));
    end

else
    loss_R = [];          %losses of signal (ROVER)
    loss_M = [];          %losses of signal (MASTER)    
end

%-------------------------------------------------------------------------------