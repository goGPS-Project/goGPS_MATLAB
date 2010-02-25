function goGPS_ublox_monitor(filerootOUT)

% SYNTAX:
%   goGPS_ublox_monitor(filerootOUT)
%
% INPUT:
%   filerootOUT = output file prefix
%
% DESCRIPTION:
%   u-blox receiver monitor: stream reading, output data saving.

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

global COMportR
global rover

%------------------------------------------------------
% initialization
%------------------------------------------------------

Eph = zeros(21,32);

%------------------------------------------------------
% initialization to save data
%------------------------------------------------------

%dep_rover = [];      % rover binary stream save

%dep_time_R = [];     % time variable save
%dep_pr_R   = [];     % code variable save
%dep_ph_R   = [];     % phase variable save
%dep_snr_R  = [];     % s/n ratio variable save

%------------------------------------------------------
% data file creation
%------------------------------------------------------

%rover binary stream (uint8)
fid_rover = fopen([filerootOUT '_rover_00.bin'],'w+');

%input observations
%  time_GPS --> double, [1,1]  --> zeros(1,1)
%  time_M   --> double, [1,1]  --> zeros(1,1)
%  time_R   --> double, [1,1]
%  pr_M     --> double, [32,1] --> zeros(32,1)
%  pr_R     --> double, [32,1]
%  ph_M     --> double, [32,1] --> zeros(32,1)
%  ph_R     --> double, [32,1]
%  snr_M    --> double, [32,1] --> zeros(32,1)
%  snr_R    --> double, [32,1]
fid_obs = fopen([filerootOUT '_obs_00.bin'],'w+');

%input ephemerides
%  timeGPS  --> double, [1,1]   --> zeros(1,1)
%  time_M   --> double, [21,32]
fid_eph = fopen([filerootOUT '_eph_00.bin'],'w+');

%nmea sentences
fid_nmea = fopen([filerootOUT '_ublox_NMEA.txt'],'wt');

%------------------------------------------------------
% creation of the rover connection (u-blox)
%------------------------------------------------------

% find a serial port object.
obj1 = instrfind('Type', 'serial', 'Port', COMportR, 'Tag', '');

% if a serial object already exists, delete it before creating a new one
if ~isempty(obj1)
    delete(obj1);
end

% serial object creation
rover = serial (COMportR,'BaudRate',57600);
set(rover,'InputBufferSize',16384);
set(rover,'FlowControl','hardware');
set(rover,'RequestToSend','on');
fopen(rover);

%------------------------------------------------------
% u-blox rover configuration
%------------------------------------------------------

% save receiver configuration
fprintf('Saving u-blox receiver configuration...\n');

reply_save = ublox_CFG_CFG(rover, 'save');
tries = 0;

while (~reply_save)
    tries = tries + 1;
    if (tries > 3)
        disp('It was not possible to save the receiver configuration.');
        break
    end
    %close and delete old serial object
    try
        fclose(rover);
        delete(rover);
    catch
        stopasync(rover);
        fclose(rover);
        delete(rover);
    end
    % create new serial object
    rover = serial (COMportR,'BaudRate',57600);
    set(rover,'InputBufferSize',16384);
    set(rover,'FlowControl','hardware');
    set(rover,'RequestToSend','on');
    fopen(rover);
    reply_save = ublox_CFG_CFG(rover, 'save');
end

% set output rate to 1Hz
fprintf('Setting measurement rate to 1Hz...\n');

% ublox_poll_message(rover, '06', '08', 0);

reply_RATE = ublox_CFG_RATE(rover, 1000, 1, 1);
tries = 0;

while (~reply_RATE)
    tries = tries + 1;
    if (tries > 3)
        disp('It was not possible to set the receiver output rate to 1Hz.');
        break
    end
    %close and delete old serial object
    try
        fclose(rover);
        delete(rover);
    catch
        stopasync(rover);
        fclose(rover);
        delete(rover);
    end
    % create new serial object
    rover = serial (COMportR,'BaudRate',57600);
    set(rover,'InputBufferSize',16384);
    set(rover,'FlowControl','hardware');
    set(rover,'RequestToSend','on');
    fopen(rover);
    reply_RATE = ublox_CFG_RATE(rover, 1000, 1, 1);
end

% enable raw measurements output
fprintf('Enabling u-blox receiver RAW measurements...\n');

reply_RAW = ublox_CFG_MSG(rover, 'RXM', 'RAW', 1);
tries = 0;

while (~reply_RAW)
    tries = tries + 1;
    if (tries > 3)
        disp('It was not possible to configure the receiver to provide RAW data.');
        break
    end
    %close and delete old serial object
    try
        fclose(rover);
        delete(rover);
    catch
        stopasync(rover);
        fclose(rover);
        delete(rover);
    end
    % create new serial object
    rover = serial (COMportR,'BaudRate',57600);
    set(rover,'InputBufferSize',16384);
    set(rover,'FlowControl','hardware');
    set(rover,'RequestToSend','on');
    fopen(rover);
    reply_RAW = ublox_CFG_MSG(rover, 'RXM', 'RAW', 1);
end

% disable all NMEA messages
fprintf('Disabling u-blox receiver NMEA messages:\n');

% ublox_CFG_MSG(rover, 'NMEA', 'GGA', 0); fprintf('Disabling GGA... ');
ublox_CFG_MSG(rover, 'NMEA', 'GLL', 0); fprintf('GLL ');
ublox_CFG_MSG(rover, 'NMEA', 'GSA', 0); fprintf('GSA ');
ublox_CFG_MSG(rover, 'NMEA', 'GSV', 0); fprintf('GSV ');
ublox_CFG_MSG(rover, 'NMEA', 'RMC', 0); fprintf('RMC ');
ublox_CFG_MSG(rover, 'NMEA', 'VTG', 0); fprintf('VTG ');
ublox_CFG_MSG(rover, 'NMEA', 'GRS', 0); fprintf('GRS ');
ublox_CFG_MSG(rover, 'NMEA', 'GST', 0); fprintf('GST ');
ublox_CFG_MSG(rover, 'NMEA', 'ZDA', 0); fprintf('ZDA ');
ublox_CFG_MSG(rover, 'NMEA', 'GBS', 0); fprintf('GBS ');
ublox_CFG_MSG(rover, 'NMEA', 'DTM', 0); fprintf('DTM ');
ublox_CFG_MSG(rover, 'PUBX', '00', 0); fprintf('PUBX00 ');
ublox_CFG_MSG(rover, 'PUBX', '01', 0); fprintf('PUBX01 ');
ublox_CFG_MSG(rover, 'PUBX', '03', 0); fprintf('PUBX03 ');
ublox_CFG_MSG(rover, 'PUBX', '04', 0); fprintf('PUBX04\n');

%------------------------------------------------------
% absolute time startup
%------------------------------------------------------

tic

%initialize log file
delete('../data/UBLOX_monitor.txt');
diary('../data/UBLOX_monitor.txt');
diary on

%------------------------------------------------------
% read header package (transmission start)
%------------------------------------------------------

%visualization
fprintf('\n');
fprintf('LOCK-PHASE (HEADER PACKAGE)\n');

%initialization
rover_1 = 0;
rover_2 = 0;

%starting epoch determination
while (rover_1 ~= rover_2) | (rover_1 == 0)

    %starting time
    current_time = toc;

    %serial port checking
    rover_1 = get(rover,'BytesAvailable');
    pause(0.05);
    rover_2 = get(rover,'BytesAvailable');

    %visualization
    fprintf('u-blox: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time, rover_1, rover_2);

end

%clear the serial port (data not decoded)
data_rover = fread(rover,rover_1,'uint8'); %#ok<NASGU>

%--------------------------------------------------------
% read 1st message (used only for synchronization)
%--------------------------------------------------------

%visualization
fprintf('\n');
fprintf('LOCK-PHASE (FIRST DATA PACKAGE)\n');

%initialization
rover_1 = 0;
rover_2 = 0;

%starting epoch determination
while (rover_1 ~= rover_2) | (rover_1 == 0)

    %starting time
    current_time = toc;

    %serial port checking
    rover_1 = get(rover,'BytesAvailable');
    pause(0.05);
    rover_2 = get(rover,'BytesAvailable');

    %visualization
    fprintf('u-blox: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time, rover_1, rover_2);

end

%clear the serial port (data not decoded)
data_rover = fread(rover,rover_1,'uint8'); %#ok<NASGU>

%set the starting time
safety_lag = 0.1;                       %safety lag on ROVER data reading
start_time = current_time-safety_lag;   %starting time

%--------------------------------------------------------
% data reading and saving
%--------------------------------------------------------

%visualization
fprintf('\n');
fprintf('ACQUISITION-PHASE\n');

%counter initialization
t = 0;

%loop control initialization
f1 = figure;
s1 = get(0,'ScreenSize');
set(f1, 'position', [s1(3)-240-20 s1(4)-80-40 240 80], 'menubar', 'none', 'name', 'UBLOX monitor');
h1 = uicontrol(gcf, 'style', 'pushbutton', 'position', [80 20 80 40], 'string', 'STOP', ...
    'callback', 'setappdata(gcf, ''run'', 0)'); %#ok<NASGU>
flag = 1;
setappdata(gcf, 'run', flag);

%ephemerides poll flag
eph_polled = 0;

%infinite loop
while flag

    %time reading (relative to start_time)
    current_time = toc;

    %serial port checking
    rover_1 = get(rover,'BytesAvailable');
    pause(0.05);
    rover_2 = get(rover,'BytesAvailable');

    %test if the package writing is finished
    if (rover_1 == rover_2) & (rover_1 ~= 0)

        data_rover = fread(rover,rover_1,'uint8');     %serial port reading
        fwrite(fid_rover,data_rover,'uint8');          %transmitted stream save
        data_rover = dec2bin(data_rover,8);            %conversion to binary (N x 8bit matrix)
        data_rover = data_rover';                      %transpose (8bit x N matrix)
        data_rover = data_rover(:)';                   %conversion to string (8N bit vector)
        %dep_rover = strcat(dep_rover,data_rover);

        [cell_rover, nmea_string] = decode_ublox(data_rover);

        %read data type
        type = '';

        %data type counters
        nRAW = 0;
        nEPH = 0;

        for i = 1 : size(cell_rover,2)

            %RXM-RAW message
            if (strcmp(cell_rover{1,i},'RXM-RAW'))

                time_R = cell_rover{2,i}(1);
                week_R = cell_rover{2,i}(2);
                ph_R   = cell_rover{3,i}(:,1);
                pr_R   = cell_rover{3,i}(:,2);
                dop_R  = cell_rover{3,i}(:,3);
                qual_R = cell_rover{3,i}(:,5);
                snr_R  = cell_rover{3,i}(:,6);
                lock_R = cell_rover{3,i}(:,7);

                %manage "nearly null" data
                ph_R(abs(ph_R) < 1e-100) = 0;

                %counter increment
                t = t+1;

                %data save
                fwrite(fid_obs, [0; 0; time_R; zeros(32,1); pr_R; zeros(32,1); ph_R; zeros(32,1); snr_R; zeros(3,1)], 'double');
                fwrite(fid_eph, [0; Eph(:)], 'double');
                %dep_time_R (t) = time_R;
                %dep_pr_R(:,t)  = pr_R;
                %dep_ph_R(:,t)  = ph_R;
                %dep_snr_R(:,t) = snr_R;

                type = [type 'RXM-RAW '];
                nRAW = nRAW + 1;

            %RXM-EPH message
            elseif (strcmp(cell_rover{1,i},'RXM-EPH'))

                %satellite number
                sat = cell_rover{2,i}(1);

                Eph(:, sat) = cell_rover{2,i}(:);

                if (nEPH == 0)
                    type = [type 'RXM-EPH '];
                end
                nEPH = nEPH + 1;

            end

        end

        if (~isempty(nmea_string))
            fprintf(fid_nmea, '%s', nmea_string);

            type = [type 'NMEA '];
        end

        %----------------------------------

        %visualization
        fprintf('\n');
        fprintf('u-blox: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time-start_time, rover_1, rover_2);
        fprintf('MSG types: %s\n', type);

        %visualization (RXM-RAW information)
        if (nRAW > 0)
            sat_pr = find(pr_R ~= 0);       %satellites with code available
            sat_ph = find(ph_R ~= 0);       %satellites with phase available
            sat = union(sat_pr,sat_ph);     %satellites with code or phase available

            if (i < length(time_R)), fprintf(' DELAYED\n'); else fprintf('\n'); end
            fprintf('Epoch %3d:  GPStime=%d:%.3f (%d satellites)\n', t, week_R, time_R, length(sat));
            for j = 1 : length(sat)
                fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f  D1=%7.1f  QI=%1d  SNR=%2d  LOCK=%1d\n', ...
                    sat(j), pr_R(sat(j)), ph_R(sat(j)), dop_R(sat(j)), qual_R(sat(j)), snr_R(sat(j)), lock_R(sat(j)));
            end
        end

        %visualization (RXM-EPH information)
        if (nEPH > 0)
            sat = find(sum(abs(Eph))>0);
            fprintf('Eph: ');
            for i = 1 : length(sat)
                fprintf('%d ', sat(i));
            end
            fprintf('\n');
        end

        %poll a new RXM-EPH message every 10 epochs
        if (mod(current_time-start_time,10) < 1)
            if (eph_polled == 0)
                ublox_poll_message(rover, 'RXM', 'EPH', 0);
                eph_polled = 1;
            end
        else
            eph_polled = 0;
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

%load u-blox saved configuration
if (reply_save)
    fprintf('Restoring saved u-blox receiver configuration...\n');

    reply_load = ublox_CFG_CFG(rover, 'load');
    tries = 0;

    while (~reply_load)
        tries = tries + 1;
        if (tries > 3)
            disp('It was not possible to reload the receiver previous configuration.');
            break
        end
        reply_load = ublox_CFG_CFG(rover, 'load');
    end
end

%connection closing
fclose(rover);
delete(rover);

%data files closing
fclose(fid_rover);
fclose(fid_obs);
fclose(fid_eph);
fclose(fid_nmea);

%log file closing
diary off

%figure closing
close(f1);