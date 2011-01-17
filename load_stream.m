function [time_GPS, week_R, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, snr_R, snr_M, pos_M, Eph, ...
          iono, loss_R, loss_M, data_rover_all, data_master_all, nmea_sentences] = load_stream (fileroot, wait_dlg)

% SYNTAX:
%   [time_GPS, week_R, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, snr_R, snr_M, pos_M, Eph, ...
%    iono, loss_R, loss_M, data_rover_all, data_master_all, nmea_sentences] = load_stream (fileroot, wait_dlg);
%
% INPUT:
%   fileroot = name of the file to be read
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   time_GPS = reference GPS time
%   week_R   = GPS week
%   time_R   = GPS time for the ROVER observations
%   time_M   = GPS time for the MASTER observations
%   pr1_R    = ROVER-SATELLITE code-pseudorange (carrier L1)
%   pr1_M    = MASTER-SATELLITE code-pseudorange (carrier L1)
%   ph1_R    = ROVER-SATELLITE phase observations (carrier L1)
%   ph1_M    = MASTER-SATELLITE phase observations (carrier L1)
%   Eph      = matrix of 29 ephemerides for each satellite
%   iono     = ionosphere parameters
%   loss_R   = flag for the ROVER loss of signal
%   loss_M   = flag for the MASTER loss of signal
%   data_rover_all  = ROVER overall stream
%   data_master_all = MASTER overall stream
%   nmea_sentences = cell-array containing the receiver NMEA sentences
%
% DESCRIPTION:
%   Reading of the streams trasmitted by the u-blox receiver (ROVER) and
%   by the permanent station (MASTER) in RTCM format.

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

global lambda1

nmea_sentences = cell(0);

if (nargin == 2)
    waitbar(0.5,wait_dlg,'Reading rover stream files...')
end

%ROVER stream reading
data_rover_all = [];                                                 %overall stream
hour = 0;                                                            %hour index (integer)
hour_str = num2str(hour,'%02d');                                     %hour index (string)
d = dir([fileroot '_rover_' hour_str '.bin']);                       %file to be read
while ~isempty(d)
    if (nargin == 1)
        fprintf(['Reading: ' fileroot '_rover_' hour_str '.bin\n']);
    end
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

%ROVER stream reading (.ubx file)
d = dir(fileroot);                                                   %file to be read
if ~isempty(d)
    if (nargin == 1)
        fprintf(['Reading: ' fileroot '\n']);
    end
    num_bytes = d.bytes;                                             %file size (number of bytes)
    fid_rover = fopen(fileroot);                                     %file opening
    data_rover_all = fread(fid_rover,num_bytes,'uint8');             %file reading
    data_rover_all = dec2bin(data_rover_all,8);                      %conversion in binary number (N x 8bits matrix)
    data_rover_all = data_rover_all';                                %transposed (8bits x N matrix)
    data_rover_all = data_rover_all(:)';                             %conversion into a string (8N bits vector)
    fclose(fid_rover);                                               %file closing
end

if (nargin == 2)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

if ~isempty(data_rover_all)
    
    %display
    if (nargin == 1)
        fprintf('Decoding rover stream\n');
    end
    
    %detect binary format
    header1 = 'B5';      % UBX header (hexadecimal value)
    header2 = '62';      % UBX header (hexadecimal value)
    codeHEX = [header1 header2];                % initial hexadecimal stream
    codeBIN = dec2bin(hex2dec(codeHEX),16);     % initial binary stream
    pos_UBX = findstr(data_rover_all, codeBIN); % message initial index

    header1 = 'A0';      % SkyTraq header (hexadecimal value)
    header2 = 'A1';      % SkyTraq header (hexadecimal value)
    codeHEX = [header1 header2];                % initial hexadecimal stream
    codeBIN = dec2bin(hex2dec(codeHEX),16);     % initial binary stream
    pos_STQ = findstr(data_rover_all, codeBIN); % message initial index

    if (length(pos_UBX) > length(pos_STQ))
        %UBX format decoding
        if (nargin == 2)
            [cell_rover, nmea_sentences] = decode_ublox(data_rover_all, wait_dlg);
        else
            [cell_rover, nmea_sentences] = decode_ublox(data_rover_all);
        end
    else
        %SkyTraq format decoding
        if (nargin == 2)
            [cell_rover, nmea_sentences] = decode_skytraq(data_rover_all, wait_dlg);
        else
            [cell_rover, nmea_sentences] = decode_skytraq(data_rover_all);
        end
    end

    %initialization (to make the writing faster)
    Ncell  = size(cell_rover,2);                          %number of read UBX messages
    time_R = zeros(Ncell,1);                              %GPS time of week
    week_R = zeros(Ncell,1);                              %GPS week
    pr1_R  = zeros(32,Ncell);                             %code observations
    ph1_R  = zeros(32,Ncell);                             %phase observations
    snr_R  = zeros(32,Ncell);                             %signal-to-noise ratio
    Eph_R  = zeros(29,32,Ncell);                          %ephemerides
    iono   = zeros(8,Ncell);                              %ionosphere parameters

    if (nargin == 2)
        waitbar(0,wait_dlg,'Reading rover data...')
    end

    i = 1;
    for j = 1 : Ncell
        if (nargin == 2)
            waitbar(j/Ncell,wait_dlg)
        end
    
        %%%%%%%%%%%%%%%%%%%%%% UBX messages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (strcmp(cell_rover{1,j},'RXM-RAW'))            %RXM-RAW message data
            time_R(i) = round(cell_rover{2,j}(1));
            week_R(i) = cell_rover{2,j}(2);
            pr1_R(:,i) = cell_rover{3,j}(:,2);
            ph1_R(:,i) = cell_rover{3,j}(:,1);
            snr_R(:,i) = cell_rover{3,j}(:,6);

            %manage "nearly null" data
            pos = abs(ph1_R(:,i)) < 1e-100;
            ph1_R(pos,i) = 0;

            %phase adjustement
            pos = abs(ph1_R(:,i)) > 0 & abs(ph1_R(:,i)) < 1e7;
            if(sum(pos) ~= 0)
                ambig = 2^23;
                n = floor((pr1_R(pos,i)/lambda1-ph1_R(pos,i)) / ambig + 0.5 );
                ph1_R(pos,i) = ph1_R(pos,i) + n*ambig;
            end

            i = i + 1;

            Eph_R(:,:,i) = Eph_R(:,:,i-1);          %previous epoch ephemerides copying
            iono(:, i) = iono(:, i-1);              %previous epoch iono parameters copying
            
        %RXM-SFRB message data save
        elseif (strcmp(cell_rover{1,j},'RXM-SFRB'))
            
            if (sum(cell_rover{2,j}(1:8)) ~= 0)
                %ionosphere parameters
                iono(:, i) = cell_rover{2,j}(1:8);
            end

        %RXM-EPH message data save
        elseif (strcmp(cell_rover{1,j},'RXM-EPH'))

            %satellite number
            sat = cell_rover{2,j}(1);

            Eph_R(:,sat,i) = cell_rover{2,j}(:);
            
        %AID-EPH message data save
        elseif (strcmp(cell_rover{1,j},'AID-EPH'))

            %satellite number
            sat = cell_rover{2,j}(1);

            Eph_R(:,sat,i) = cell_rover{2,j}(:);
            
        %AID-HUI message data save
        elseif (strcmp(cell_rover{1,j},'AID-HUI'))

            %ionosphere parameters
            iono(:, i) = cell_rover{3,j}(9:16);
            
        %%%%%%%%%%%%%%%%%% SkyTraq messages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %MEAS_TIME message data save
        elseif (strcmp(cell_rover{1,j},'MEAS_TIME'))

            time_R(i) = round(cell_rover{2,j}(3));
            week_R(i) = cell_rover{2,j}(2);
            
        %RAW_MEAS message data save
        elseif (strcmp(cell_rover{1,j},'RAW_MEAS'))

            pr1_R(:,i) = cell_rover{3,j}(:,3);
            ph1_R(:,i) = cell_rover{3,j}(:,4);
            snr_R(:,i) = cell_rover{3,j}(:,2);

            %manage "nearly null" data
            pos = abs(ph1_R(:,i)) < 1e-100;
            ph1_R(pos,i) = 0;

            %phase adjustement
            pos = abs(ph1_R(:,i)) > 0 & abs(ph1_R(:,i)) < 1e8;
            if(sum(pos) ~= 0)
                ambig = 2^23;
                n = floor((pr1_R(pos,i)/lambda1-ph1_R(pos,i)) / ambig + 0.5 );
                ph1_R(pos,i) = ph1_R(pos,i) + n*ambig;
            end

            i = i + 1;

            Eph_R(:,:,i) = Eph_R(:,:,i-1);          %previous epoch ephemerides copying
            iono(:, i) = iono(:, i-1);              %previous epoch iono parameters copying

        end
    end

    %residual data erase (after initialization)
    time_R(i:end)  = [];
    week_R(i:end)  = [];
    pr1_R(:,i:end) = [];
    ph1_R(:,i:end) = [];
    Eph_R(:,:,i:end) = [];
    iono(:,i:end) = [];

else
    %displaying
    if (nargin == 2)
        msgbox('No rover data acquired.');
    else
        fprintf('No rover data acquired! \n');
    end

    time_R = [];                         %GPS time
    week_R = [];                         %GPS week
    pr1_R  = [];                         %code observations
    ph1_R  = [];                         %phase observations
    snr_R  = [];                         %signal-to-noise ratio
    Eph_R  = [];                         %ephemerides
    iono   = [];                         %ionosphere parameters

end

%-------------------------------------------------------------------------------

if (nargin == 2)
    waitbar(0.5,wait_dlg,'Reading master stream files...')
end

%MASTER stream reading
data_master_all = [];                                                %overall stream
hour = 0;                                                            %hour index (integer)
hour_str = num2str(hour,'%02d');                                     %hour index (string)
d = dir([fileroot '_master_' hour_str '.bin']);                      %file to be read
while ~isempty(d)
    if (nargin == 1)
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

if (nargin == 2)
    waitbar(1,wait_dlg)
end

%-------------------------------------------------------------------------------

if ~isempty(data_master_all)

    %display
    if (nargin == 1)
        fprintf('Decoding master stream\n');
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
        [cell_master] = decode_rtcm2(sixofeight);
    else
        if (nargin == 2)
            [cell_master] = decode_rtcm3(data_master_all, wait_dlg);
        else
            [cell_master] = decode_rtcm3(data_master_all);
        end
    end

    %initialization (to make the writing faster)
    Ncell  = size(cell_master,2);                         %number of read RTCM packets
    time_M = zeros(Ncell,1);                              %GPS time
    pr1_M  = zeros(32,Ncell);                             %code observations
    ph1_M  = zeros(32,Ncell);                             %phase observations
    snr_M  = zeros(32,Ncell);                             %signal-to-noise ratio
    pos_M  = zeros(3,Ncell);                              %master station position
    Eph_M  = zeros(29,32,Ncell);                          %ephemerides

    if (nargin == 2)
        waitbar(0,wait_dlg,'Reading master data...')
    end

    i = 1;
    for j = 1 : Ncell
        if (nargin == 2)
            waitbar(j/Ncell,wait_dlg)
        end
        if (cell_master{1,j} == 3)
        elseif (cell_master{1,j} == 18)
        elseif (cell_master{1,j} == 19)
        elseif (cell_master{1,j} == 1002) | (cell_master{1,j} == 1004) %RTCM 1002 or 1004 message

            time_M(i)  = cell_master{2,j}(2);             %GPS time logging
            pr1_M(:,i) = cell_master{3,j}(:,2);           %code observations logging
            ph1_M(:,i) = cell_master{3,j}(:,3);           %phase observations logging
            snr_M(:,i) = cell_master{3,j}(:,5);           %signal-to-noise ratio logging

            i = i+1;                                      %epoch counter increase

            pos_M(:,i) = pos_M(:,i-1);                    %previous epoch master coordinates copying
            Eph_M(:,:,i) = Eph_M(:,:,i-1);                %previous epoch ephemerides copying

        elseif (cell_master{1,j} == 1005)                 %RTCM 1005 message

            coordX_M = cell_master{2,j}(8);
            coordY_M = cell_master{2,j}(9);
            coordZ_M = cell_master{2,j}(10);

            pos_M(:,i) = [coordX_M; coordY_M; coordZ_M];
            
            if (i == 2) & (pos_M(:,1) == 0)
                pos_M(:,i-1) = pos_M(:,i);
            end

        elseif (cell_master{1,j} == 1006)                 %RTCM 1006 message

            coordX_M = cell_master{2,j}(8);
            coordY_M = cell_master{2,j}(9);
            coordZ_M = cell_master{2,j}(10);

            pos_M(:,i) = [coordX_M; coordY_M; coordZ_M];
            
            if (i == 2) & (pos_M(:,1) == 0)
                pos_M(:,i-1) = pos_M(:,i);
            end

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
    pos_M(:,i:end) = [];
    Eph_M(:,:,i:end) = [];

    ph1_M(ph1_M < 1e-100) = 0;

else
    %displaying
    if (nargin == 2)
        msgbox('No master data acquired.');
    else
        fprintf('No master data acquired! \n');
    end

    time_M = [];                         %GPS time
    pr1_M  = [];                         %code observations
    ph1_M  = [];                         %phase observations
    snr_M  = [];                         %signal-to-noise ratio
    pos_M  = [];                         %master station position
    Eph_M  = [];                         %ephemerides
end

%-------------------------------------------------------------------------------

if (nargin == 2)
    waitbar(0.33,wait_dlg,'Synchronizing data...')
end

if ~isempty(time_R) & ~isempty(time_M)

    %head synchronization
    if (time_R(1) < time_M(1))
        pos = find(time_R < time_M(1));
        time_R(pos)    = [];                       %GPS time
        week_R(pos)    = [];                       %GPS week
        pr1_R(:,pos)   = [];                       %code observations
        ph1_R(:,pos)   = [];                       %phase observations
        snr_R(:,pos)   = [];                       %signal-to-noise ratio
        Eph_R(:,:,pos) = [];                       %ephemerides
        iono(:,pos) = [];                          %ionosphere parameters
    end

    if (time_M(1) < time_R(1))
        pos = find(time_M < time_R(1));
        time_M(pos)    = [];                       %GPS time
        pr1_M(:,pos)   = [];                       %code observations
        ph1_M(:,pos)   = [];                       %phase observations
        snr_M(:,pos)   = [];                       %signal-to-noise ratio
        pos_M(:,pos)   = [];                       %master station position
        Eph_M(:,:,pos) = [];                       %ephemerides
    end

    %tail synchronization
    if (time_R(end) > time_M(end))
        pos = find(time_R > time_M(end));
        time_R(pos)    = [];                       %GPS time
        week_R(pos)    = [];                       %GPS week
        pr1_R(:,pos)   = [];                       %code observations
        ph1_R(:,pos)   = [];                       %phase observations
        snr_R(:,pos)   = [];                       %signal-to-noise ratio
        Eph_R(:,:,pos) = [];                       %ephemerides
        iono(:,pos) = [];                          %ionosphere parameters
    end

    if (time_M(end) > time_R(end))
        pos = find(time_M > time_R(end));
        time_M(pos)    = [];                       %GPS time
        pr1_M(:,pos)   = [];                       %code observations
        ph1_M(:,pos)   = [];                       %phase observations
        snr_M(:,pos)   = [];                       %signal-to-noise ratio
        pos_M(:,pos)   = [];                       %master station position
        Eph_M(:,:,pos) = [];                       %ephemerides
    end

end

%-------------------------------------------------------------------------------

if (nargin == 2)
    waitbar(0.66,wait_dlg)
end

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
            week_R = [week_R(1:pos);  0;             week_R(pos+1:end)];
            pr1_R  = [pr1_R(:,1:pos)  zeros(32,1)    pr1_R(:,pos+1:end)];
            ph1_R  = [ph1_R(:,1:pos)  zeros(32,1)    ph1_R(:,pos+1:end)];
            snr_R  = [snr_R(:,1:pos)  zeros(32,1)    snr_R(:,pos+1:end)];
            iono   = [iono(:,1:pos)   zeros(8,1)     iono(:,pos+1:end)];

            Eph_R  = cat(3, Eph_R(:,:,1:pos), zeros(29,32,1), Eph_R(:,:,pos+1:end));
        end
    else
        time_R = time_GPS;
        week_R = zeros(1,length(time_GPS));
        pr1_R  = zeros(32,length(time_GPS));
        ph1_R  = zeros(32,length(time_GPS));
        snr_R  = zeros(32,length(time_GPS));
        Eph_R  = zeros(29,32,length(time_GPS));
        iono   = zeros(8,length(time_GPS));
    end

    if ~isempty(time_M)

        newtime_M = setdiff(time_GPS, time_M);       %MASTER missing epochs
        for i = 1 : length(newtime_M)

            pos = find(time_M == newtime_M(i) - 1);  %position before the "holes"

            time_M = [time_M(1:pos);  newtime_M(i);  time_M(pos+1:end)];
            pr1_M  = [pr1_M(:,1:pos)  zeros(32,1)    pr1_M(:,pos+1:end)];
            ph1_M  = [ph1_M(:,1:pos)  zeros(32,1)    ph1_M(:,pos+1:end)];
            snr_M  = [snr_M(:,1:pos)  zeros(32,1)    snr_M(:,pos+1:end)];
            pos_M  = [pos_M(:,1:pos)  zeros(3,1)     pos_M(:,pos+1:end)];

            Eph_M  = cat(3, Eph_M(:,:,1:pos), zeros(29,32,1), Eph_M(:,:,pos+1:end));
        end
    else
        time_M = time_GPS;
        pr1_M  = zeros(32,length(time_GPS));
        ph1_M  = zeros(32,length(time_GPS));
        snr_M  = zeros(32,length(time_GPS));
        pos_M  = zeros(3,length(time_GPS));
        Eph_M  = zeros(29,32,length(time_GPS));
    end

else
    loss_R = [];          %losses of signal (ROVER)
    loss_M = [];          %losses of signal (MASTER)
end

if (nargin == 2)
    waitbar(1,wait_dlg)
end

%if ephemerides coming from RTCM stream are not available, use rover ones
if (~Eph_M)
    Eph = Eph_R;
else
    Eph = Eph_M;
end
