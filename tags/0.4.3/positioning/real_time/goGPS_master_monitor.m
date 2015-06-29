function goGPS_master_monitor(filerootOUT, flag_NTRIP)

% SYNTAX:
%   goGPS_master_monitor(filerootOUT)
%
% INPUT:
%   filerootOUT = output file prefix
%   flag_NTRIP = use/don't use NTRIP flag
%
% DESCRIPTION:
%   Master station monitor: stream reading, output data saving.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

global master_ip master_port
global nmea_init nmea_update_rate
global master
global server_delay

%------------------------------------------------------
% initialization
%------------------------------------------------------

num_sat = 32;
Eph = zeros(33,num_sat);

%counter for creating hourly files
hour = 0;

%do not overwrite existing files
i = 1;
j = length(filerootOUT);
while (~isempty(dir([filerootOUT '_obs*.bin'])) || ...
        ~isempty(dir([filerootOUT '_eph*.bin'])) )
    
    filerootOUT(j+1:j+4) = ['_' num2str(i,'%03d')];
    i = i + 1;
end

%------------------------------------------------------
% data file creation
%------------------------------------------------------

%master binary stream (uint8)
fid_master = fopen([filerootOUT '_master_000.bin'],'w+');

%input observations
%  time_GPS --> double, [1,1]   --> zeros(1,1)
%  time_M   --> double, [1,1]
%  time_R   --> double, [1,1]   --> zeros(1,1)
%  pr_M     --> double, [num_sat,1]
%  pr_R     --> double, [num_sat,1]  --> zeros(num_sat,1)
%  ph_M     --> double, [num_sat,1]
%  ph_R     --> double, [num_sat,1]  --> zeros(num_sat,1)
%  snr_M    --> double, [num_sat,1]
%  snr_R    --> double, [num_sat,1]  --> zeros(num_sat,1)
fid_obs = fopen([filerootOUT '_obs_000.bin'],'w+');

%input ephemerides
%  timeGPS  --> double, [1,1]   --> zeros(1,1)
%  Eph      --> double, [33,num_sat]
fid_eph = fopen([filerootOUT '_eph_000.bin'],'w+');

%write number of satellites
fwrite(fid_obs, num_sat, 'int8');
fwrite(fid_eph, num_sat, 'int8');

%------------------------------------------------------
% creation of the connection to the master
%------------------------------------------------------

master = tcpip(master_ip,master_port);
set(master,'InputBufferSize', 16384);
fopen(master);

if (flag_NTRIP)
    ntripstring = NTRIP_string_generator(nmea_init);
    %fprintf('NTRIP request [%s]',ntripstring);
    fwrite(master,ntripstring);
end

%------------------------------------------------------
% absolute time startup
%------------------------------------------------------

tic

%file writing initialization
delete('../data/MASTER_monitor.txt');
diary('../data/MASTER_monitor.txt');
diary on

%--------------------------------------------------------
% read 1st message (used only for synchronization)
%--------------------------------------------------------

%visualization
fprintf('\n');
fprintf('LOCK-PHASE (FIRST DATA PACKAGE)\n');

%initialization
master_1 = 0;
master_2 = 0;

%starting epoch determination
while (master_1 ~= master_2) || (master_1 == 0)

    %starting time
    current_time = toc;

    %TCP/IP port checking
    master_1 = get(master,'BytesAvailable');
    pause(0.05);
    master_2 = get(master,'BytesAvailable');

    %visualization
    fprintf('master: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time, master_1, master_2);

end

%clear TCP/IP port (data not decoded)
data_master = fread(master,master_1,'uint8'); %#ok<NASGU>

%set starting time
safety_lag = 0.1;                       %safety lag on MASTER data reading
start_time = current_time-safety_lag;   %starting time

%--------------------------------------------------------
% data reading and saving
%--------------------------------------------------------

%visualization
fprintf('\n');
fprintf('ACQUISITION-PHASE\n');

%counter initialization
t = 0;

%master position initialization
pos_M = zeros(3,1);

%loop control initialization
f1 = figure;
s1 = get(0,'ScreenSize');
set(f1, 'position', [s1(3)-240-20 s1(4)-80-40 240 80], 'menubar', 'none', 'name', 'MASTER monitor');
h1 = uicontrol(gcf, 'style', 'pushbutton', 'position', [80 20 80 40], 'string', 'STOP', ...
    'callback', 'setappdata(gcf, ''run'', 0)'); %#ok<NASGU>
flag = 1;
setappdata(gcf, 'run', flag);

%nmea flag
nmea_sent = 0;

%waiting time
waiting_time_start = toc;

%approximate message rate
approx_msg_rate = 1;

%infinite loop
while flag
    
    %time reading
    current_time = toc;

    %-------------------------------------
    % hourly files
    %-------------------------------------
    if (floor(current_time/3600) > hour)
        
        hour = floor(current_time/3600);
        hour_str = num2str(hour,'%03d');
        
        fclose(fid_master);
        fclose(fid_obs);
        fclose(fid_eph);
        
        fid_master = fopen([filerootOUT '_master_'  hour_str '.bin'],'w+');
        fid_obs    = fopen([filerootOUT '_obs_'    hour_str '.bin'],'w+');
        fid_eph    = fopen([filerootOUT '_eph_'    hour_str '.bin'],'w+');
        
        %write number of satellites
        fwrite(fid_obs, num_sat, 'int8');
        fwrite(fid_eph, num_sat, 'int8');
    end

    %TCP/IP port checking
    master_1 = get(master,'BytesAvailable');
    pause(server_delay);
    master_2 = get(master,'BytesAvailable');

    %test if the package writing is finished
    if (master_1 == master_2) && (master_1 ~= 0)
        
        %approximate message rate
        approx_msg_rate = max([1 round(current_time - waiting_time_start)]);
        
        %reset waiting time start
        waiting_time_start = current_time;

        data_master = fread(master,master_1,'uint8');     %TCP/IP port reading
        fwrite(fid_master,data_master,'uint8');           %transmitted stream saving
        data_master = dec2bin(data_master,8)';            %conversion to binary (N x 8bit matrix)
        data_master = data_master (:);                    %transpose (8bit x N matrix)
        data_master = data_master';                       %conversion to string (8N bit vector)

        pos = 1;
        sixofeight = [];
        is_rtcm2 = 1;

        while (pos + 7 <= length(data_master))
            if (~strcmp(data_master(pos:pos+1),'01'))
                is_rtcm2 = 0;
                break
            end
            sixofeight = [sixofeight fliplr(data_master(pos+2:pos+7))];
            pos = pos + 8;
        end

        if(is_rtcm2)
            [cell_master] = decode_rtcm2(sixofeight);

            %read data type
            type = '';

            pr1_M = zeros(num_sat,1);
            pr2_M = zeros(num_sat,1);
            ph1_M = zeros(num_sat,1);
            ph2_M = zeros(num_sat,1);

            time_18L1 = 0;
            time_18L2 = 0;
            time_19CA = 0;
            time_19P  = 0;

            %data type counters
            n3 = 0;
            n18L1 = 0; n18L2 = 0;
            n19CA = 0; n19P = 0;

            for i = 1 : size(cell_master,2)

                %message 3 data save
                if (cell_master{1,i} == 3)

                    coordX_M = cell_master{2,i}(1);
                    coordY_M = cell_master{2,i}(2);
                    coordZ_M = cell_master{2,i}(3);

                    type = [type '3 '];
                    n3 = n3 + 1;
                end

                %message 18 data save
                if (cell_master{1,i} == 18)

                    time_M = cell_master{2,i}(2);
                    if (cell_master{2,i}(1) == 0)
                        ph1_M   = cell_master{3,i}(:,7);

                        %manage "nearly null" data
                        ph1_M(abs(ph1_M) < 1e-100) = 0;

                        time_18L1 = time_M;

                        n18L1 = n18L1 + 1;
                    elseif (cell_master{2,i}(1) == 2)
                        ph2_M   = cell_master{3,i}(:,7);

                        %manage "nearly null" data
                        ph2_M(abs(ph2_M) < 1e-100) = 0;

                        time_18L2 = time_M;

                        n18L2 = n18L2 + 1;
                    end

                    type = [type '18 '];

                end

                %message 19 data save
                if (cell_master{1,i} == 19)

                    time_M  = cell_master{2,i}(2);
                    if (cell_master{2,i}(1) == 0)
                        pr1_M   = cell_master{3,i}(:,7);

                        time_19CA = time_M;

                        %counter increment
                        t = t+1;

                        n19CA = n19CA + 1;
                    elseif (cell_master{2,i}(1) == 2)
                        pr2_M   = cell_master{3,i}(:,7);

                        time_19P = time_M;

                        n19P = n19P + 1;
                    end

                    type = [type '19 '];

                end
            end

            %keep just last packet for visualization
            time_M = max([time_18L1, time_18L2, time_19CA, time_19P]);
            if (time_18L1 ~= time_M)
                ph1_M = zeros(num_sat,1);
            end
            if (time_18L2 ~= time_M)
                ph2_M = zeros(num_sat,1);
            end
            if (time_19CA ~= time_M)
                pr1_M = zeros(num_sat,1);
                ph1_M = zeros(num_sat,1);
                ph2_M = zeros(num_sat,1);
            end
            if (time_19P ~= time_M)
                pr2_M = zeros(num_sat,1);
            end

            %visualization
            fprintf('\n');
            fprintf('master: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time-start_time, master_1, master_2);
            fprintf('MSG types: %s\n', type);

            if ((n18L1 > 0 || n18L2 > 0) && n19CA > 0)
                sat_pr = find(pr1_M ~= 0);
                sat_ph = find(ph1_M ~= 0);
                sat1 = union(sat_pr,sat_ph);

                sat_pr = find(pr2_M ~= 0);
                sat_ph = find(ph2_M ~= 0);
                sat2 = union(sat_pr,sat_ph);

                fprintf('Epoch %3d: time=%.3f\n', t, time_M);
                fprintf('GPS L1 (%d satellites)\n', length(sat1));
                for i = 1 : length(sat1)
                    fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f\n', sat1(i), pr1_M(sat1(i)), ph1_M(sat1(i)));
                end
                fprintf('GPS L2 (%d satellites)\n', length(sat2));
                for i = 1 : length(sat2)
                    fprintf('   SAT %02d:  P2=%11.2f  L2=%12.2f\n', sat2(i), pr2_M(sat2(i)), ph2_M(sat2(i)));
                end
            end

            %visualization (msg 3 information)
            if (n3 > 0)
                fprintf('Coord:  X=%12.4f  Y=%12.4f  Z=%12.4f\n', coordX_M, coordY_M, coordZ_M);
            end

        else
            [cell_master] = decode_rtcm3(data_master);

            %read data type
            type = '';

            %data type counters
            n02 = 0; n03 = 0; n04 = 0;
            n05 = 0; n06 = 0;
            n07 = 0; n08 = 0;
            n10 = 0; n11 = 0; n12 = 0;
            n19 = 0;

            for i = 1 : size(cell_master,2)

                %1002 message data save
                if (cell_master{1,i} == 1002)

                    time_M = cell_master{2,i}(2);
                    pr_M   = cell_master{3,i}(:,2);
                    ph_M   = cell_master{3,i}(:,3);
                    lock_M = cell_master{3,i}(:,4);
                    snr_M  = cell_master{3,i}(:,5);

                    %manage "nearly null" data
                    ph_M(ph_M < 1e-100) = 0;

                    %counter increment
                    t = t+1;

                    type = [type '1002 '];
                    n02 = n02 + 1;

                %1003 message data save
                elseif (cell_master{1,i} == 1003)

                    time_M  = cell_master{2,i}(2);
                    pr1_M   = cell_master{3,i}(:,2);
                    ph1_M   = cell_master{3,i}(:,3);
                    lock1_M = cell_master{3,i}(:,4);
                    pr2_M   = cell_master{3,i}(:,6);
                    ph2_M   = cell_master{3,i}(:,7);
                    lock2_M = cell_master{3,i}(:,8);
                    snr1_M = zeros(num_sat,1);
                    snr2_M = zeros(num_sat,1);

                    %manage "nearly null" data
                    ph1_M(abs(ph1_M) < 1e-100) = 0;
                    ph2_M(abs(ph2_M) < 1e-100) = 0;

                    %counter increment
                    t = t+1;

                    type = [type '1003 '];
                    n03 = n03 + 1;

                %1004 message data save
                elseif (cell_master{1,i} == 1004)

                    time_M  = cell_master{2,i}(2);
                    pr1_M   = cell_master{3,i}(:,2);
                    ph1_M   = cell_master{3,i}(:,3);
                    lock1_M = cell_master{3,i}(:,4);
                    snr1_M  = cell_master{3,i}(:,5);
                    pr2_M   = cell_master{3,i}(:,7);
                    ph2_M   = cell_master{3,i}(:,8);
                    lock2_M = cell_master{3,i}(:,9);
                    snr2_M  = cell_master{3,i}(:,10);

                    %manage "nearly null" data
                    ph1_M(abs(ph1_M) < 1e-100) = 0;
                    ph2_M(abs(ph2_M) < 1e-100) = 0;

                    %counter increment
                    t = t+1;

                    type = [type '1004 '];
                    n04 = n04 + 1;

                %1005 message data save
                elseif (cell_master{1,i} == 1005)

                    coordX_M = cell_master{2,i}(8);
                    coordY_M = cell_master{2,i}(9);
                    coordZ_M = cell_master{2,i}(10);
                    
                    pos_M(:,1) = [coordX_M; coordY_M; coordZ_M];

                    type = [type '1005 '];
                    n05 = n05 + 1;

                %1006 message data save
                elseif (cell_master{1,i} == 1006)

                    coordX_M = cell_master{2,i}(8);
                    coordY_M = cell_master{2,i}(9);
                    coordZ_M = cell_master{2,i}(10);
                    height_M = cell_master{2,i}(11);
                    
                    pos_M(:,1) = [coordX_M; coordY_M; coordZ_M];

                    type = [type '1006 '];
                    n06 = n06 + 1;

                %1007 message data save
                elseif (cell_master{1,i} == 1007)

                    setup_ant_M = cell_master{2,i}(2);
                    descr_ant_M = cell_master{3,i};

                    type = [type '1007 '];
                    n07 = n07 + 1;

                %1008 message data save
                elseif (cell_master{1,i} == 1008)

                    setup_ant_M  = cell_master{2,i}(2);
                    descr_ant_M  = cell_master{3,i}(1,:);
                    serial_ant_M = cell_master{3,i}(2,:);

                    type = [type '1008 '];
                    n08 = n08 + 1;

                %1010 message data save
                elseif (cell_master{1,i} == 1010)

                    time_M_GLO  = cell_master{2,i}(2);
                    pr_M_GLO   = cell_master{3,i}(:,2);
                    ph_M_GLO   = cell_master{3,i}(:,3);
                    lock_M_GLO = cell_master{3,i}(:,4);
                    snr_M_GLO  = cell_master{3,i}(:,5);

                    %manage "nearly null" data
                    ph_M_GLO(ph_M_GLO < 1e-100) = 0;

                    type = [type '1010 '];
                    n10 = n10 + 1;

                %1011 message data save
                elseif (cell_master{1,i} == 1011)

                    time_M_GLO  = cell_master{2,i}(2);
                    pr1_M_GLO   = cell_master{3,i}(:,2);
                    ph1_M_GLO   = cell_master{3,i}(:,3);
                    lock1_M_GLO = cell_master{3,i}(:,4);
                    pr2_M_GLO   = cell_master{3,i}(:,6);
                    ph2_M_GLO   = cell_master{3,i}(:,7);
                    lock2_M_GLO = cell_master{3,i}(:,8);
                    snr1_M_GLO = zeros(num_sat,1);
                    snr2_M_GLO = zeros(num_sat,1);

                    %manage "nearly null" data
                    ph1_M_GLO(ph1_M_GLO < 1e-100) = 0;
                    ph2_M_GLO(ph2_M_GLO < 1e-100) = 0;

                    type = [type '1011 '];
                    n11 = n11 + 1;

                %1012 message data save
                elseif (cell_master{1,i} == 1012)

                    time_M_GLO  = cell_master{2,i}(2);
                    pr1_M_GLO   = cell_master{3,i}(:,2);
                    ph1_M_GLO   = cell_master{3,i}(:,3);
                    lock1_M_GLO = cell_master{3,i}(:,4);
                    snr1_M_GLO  = cell_master{3,i}(:,5);
                    pr2_M_GLO   = cell_master{3,i}(:,7);
                    ph2_M_GLO   = cell_master{3,i}(:,8);
                    lock2_M_GLO = cell_master{3,i}(:,9);
                    snr2_M_GLO  = cell_master{3,i}(:,10);

                    %manage "nearly null" data
                    ph1_M_GLO(ph1_M_GLO < 1e-100) = 0;
                    ph2_M_GLO(ph2_M_GLO < 1e-100) = 0;

                    type = [type '1012 '];
                    n12 = n12 + 1;

                %1019 message data save
                elseif (cell_master{1,i} == 1019)

                    %satellite number
                    sat = cell_master{2,i}(1);

                    Eph(:,sat) = cell_master{2,i}(:);

                    type = [type '1019 '];
                    n19 = n19 + 1;

                end
            end
            
            if (t > 0) && (any(pos_M))
                %data save
                fwrite(fid_obs, [0; time_M; 0; 0; pr1_M; zeros(num_sat,1); ph1_M; zeros(num_sat,1); zeros(num_sat,1); snr1_M; zeros(num_sat,1); pos_M(:,1); zeros(8,1)], 'double');
                fwrite(fid_eph, [0; Eph(:)], 'double');
            end

            %----------------------------------

            %visualization
            fprintf('\n');
            fprintf('master: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time-start_time, master_1, master_2);
            fprintf('MSG types: %s\n', type);

            %visualization (1002 msg information)
            if (n02 > 0)
                sat_pr = find(pr_M ~= 0);
                sat_ph = find(ph_M ~= 0);
                sat = union(sat_pr,sat_ph);

                fprintf('Epoch %3d: GPStime=%.3f (%d satellites)\n', t, time_M, length(sat));
                for i = 1 : length(sat)
                    fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat(i), pr_M(sat(i)), ph_M(sat(i)), snr_M(sat(i)), lock_M(sat(i)));
                end
            end

            %visualization (1003 msg information)
            if (n03 > 0)
                sat_pr = find(pr1_M ~= 0);
                sat_ph = find(ph1_M ~= 0);
                sat1 = union(sat_pr,sat_ph);

                sat_pr = find(pr2_M ~= 0);
                sat_ph = find(ph2_M ~= 0);
                sat2 = union(sat_pr,sat_ph);

                fprintf('Epoch %3d: GPStime=%.3f\n', t, time_M);
                fprintf('GPS L1 (%d satellites)\n', length(sat1));
                for i = 1 : length(sat1)
                    fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat1(i), pr1_M(sat1(i)), ph1_M(sat1(i)), snr1_M(sat1(i)), lock1_M(sat1(i)));
                end
                fprintf('GPS L2 (%d satellites)\n', length(sat2));
                for i = 1 : length(sat2)
                    fprintf('   SAT %02d:  P2=%11.2f  L2=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat2(i), pr2_M(sat2(i)), ph2_M(sat2(i)), snr2_M(sat2(i)), lock2_M(sat2(i)));
                end
            end

            %visualization (1004 msg information)
            if (n04 > 0)
                sat_pr = find(pr1_M ~= 0);
                sat_ph = find(ph1_M ~= 0);
                sat1 = union(sat_pr,sat_ph);

                sat_pr = find(pr2_M ~= 0);
                sat_ph = find(ph2_M ~= 0);
                sat2 = union(sat_pr,sat_ph);

                fprintf('Epoch %3d: GPStime=%.3f\n', t, time_M);
                fprintf('GPS L1 (%d satellites)\n', length(sat1));
                for i = 1 : length(sat1)
                    fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat1(i), pr1_M(sat1(i)), ph1_M(sat1(i)), snr1_M(sat1(i)), lock1_M(sat1(i)));
                end
                fprintf('GPS L2 (%d satellites)\n', length(sat2));
                for i = 1 : length(sat2)
                    fprintf('   SAT %02d:  P2=%11.2f  L2=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat2(i), pr2_M(sat2(i)), ph2_M(sat2(i)), snr2_M(sat2(i)), lock2_M(sat2(i)));
                end
            end

            %visualization (1005 msg information)
            if (n05 > 0)
                fprintf('Coord:  X=%12.4f  Y=%12.4f  Z=%12.4f\n', coordX_M, coordY_M, coordZ_M);
            end

            %visualization (1006 msg information)
            if (n06 > 0)
                fprintf('Coord:  X=%12.4f  Y=%12.4f  Z=%12.4f  h=%6.4f\n', coordX_M, coordY_M, coordZ_M, height_M);
            end

            %visualization (1007 msg information)
            if (n07 > 0)
                fprintf('Antenna: %s (setup=%d)\n', descr_ant_M, setup_ant_M);
            end

            %visualization (1008 msg information)
            if (n08 > 0)
                fprintf('Antenna: %s (setup=%d), serial: %s\n', descr_ant_M, setup_ant_M, serial_ant_M);
            end

            %visualization (1010 msg information)
            if (n10 > 0)
                sat_pr = find(pr_M_GLO ~= 0);
                sat_ph = find(ph_M_GLO ~= 0);
                sat = union(sat_pr,sat_ph);

                fprintf('Epoch %3d: GPStime=%.3f (%d satellites)\n', t, time_M, length(sat));
                for i = 1 : length(sat)
                    fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat(i), pr_M_GLO(sat(i)), ph_M_GLO(sat(i)), snr_M_GLO(sat(i)), lock_M_GLO(sat(i)));
                end
            end

            %visualization (1011 msg information)
            if (n11 > 0)
                sat_pr = find(pr1_M_GLO ~= 0);
                sat_ph = find(ph1_M_GLO ~= 0);
                sat1 = union(sat_pr,sat_ph);

                sat_pr = find(pr2_M_GLO ~= 0);
                sat_ph = find(ph2_M_GLO ~= 0);
                sat2 = union(sat_pr,sat_ph);

                fprintf('Epoch %3d: GLONASStime=%.3f\n', t, time_M_GLO);
                fprintf('GLONASS L1 (%d satellites)\n', length(sat1));
                for i = 1 : length(sat1)
                    fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat1(i), pr1_M_GLO(sat1(i)), ph1_M_GLO(sat1(i)), snr1_M_GLO(sat1(i)), lock1_M_GLO(sat1(i)));
                end
                fprintf('GLONASS L2 (%d satellites)\n', length(sat2));
                for i = 1 : length(sat2)
                    fprintf('   SAT %02d:  P2=%11.2f  L2=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat2(i), pr2_M_GLO(sat2(i)), ph2_M_GLO(sat2(i)), snr2_M_GLO(sat2(i)), lock2_M_GLO(sat2(i)));
                end
            end

            %visualization (1012 msg information)
            if (n12 > 0)
                sat_pr = find(pr1_M_GLO ~= 0);
                sat_ph = find(ph1_M_GLO ~= 0);
                sat1 = union(sat_pr,sat_ph);

                sat_pr = find(pr2_M_GLO ~= 0);
                sat_ph = find(ph2_M_GLO ~= 0);
                sat2 = union(sat_pr,sat_ph);

                fprintf('Epoch %3d: GLONASStime=%.3f\n', t, time_M_GLO);
                fprintf('GLONASS L1 (%d satellites)\n', length(sat1));
                for i = 1 : length(sat1)
                    fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat1(i), pr1_M_GLO(sat1(i)), ph1_M_GLO(sat1(i)), snr1_M_GLO(sat1(i)), lock1_M_GLO(sat1(i)));
                end
                fprintf('GLONASS L2 (%d satellites)\n', length(sat2));
                for i = 1 : length(sat2)
                    fprintf('   SAT %02d:  P2=%11.2f  L2=%12.2f  SNR=%5.2f  LOCK=%3d\n', sat2(i), pr2_M_GLO(sat2(i)), ph2_M_GLO(sat2(i)), snr2_M_GLO(sat2(i)), lock2_M_GLO(sat2(i)));
                end
            end

            %visualization (1019 msg information)
            if (n19 > 0)
                sat = find(sum(abs(Eph))>0);
                fprintf('Eph: ');
                for i = 1 : length(sat)
                    fprintf('%d ', sat(i));
                end
                fprintf('\n');
            end

        end

        %send a new NMEA string
        if (flag_NTRIP) && (mod(current_time-start_time,nmea_update_rate) < 1)
            if (nmea_sent == 0)
                nmea_update = sprintf('%s\r\n',nmea_init);
                fwrite(master,nmea_update);
                nmea_sent = 1;
            end
        else
            nmea_sent = 0;
        end
    else
        %check waiting time
        waiting_time = current_time - waiting_time_start;
        
        if (waiting_time > 10*approx_msg_rate)
            
            %display message
            fprintf('Not receiving data. Trying to reconnect... ');

            %close master connection
            fclose(master);
            
            reconnected = 0;

            while (~reconnected)
                try
                    master = tcpip(master_ip,master_port);
                    set(master,'InputBufferSize', 16384);
                    fopen(master);
                    reconnected = 1;
                catch
                    pause(5);
                end
            end
            
            if (flag_NTRIP)
                ntripstring = NTRIP_string_generator(nmea_init);
                %fprintf('NTRIP request [%s]',ntripstring);
                fwrite(master,ntripstring);
            end
            
            %wait until the buffer writing is started before continuing
            while get(master,'BytesAvailable') == 0, end;
            
            %reset waiting time start
            waiting_time_start = current_time;
            
            %display message
            fprintf('done.\n');
        end
    end

    %----------------------------------

    %test if the cycle execution has ended
    flag = getappdata(gcf, 'run');
    drawnow

end

%------------------------------------------------------
% tasks at the end of the cycle
%------------------------------------------------------

%close connection
fclose(master);
delete(master);

%close data file
fclose(fid_master);
fclose(fid_obs);
fclose(fid_eph);

%close log file
diary off

%close figure
close(f1);