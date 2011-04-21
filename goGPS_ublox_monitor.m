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

global COMportR
global rover

%------------------------------------------------------
% initialization
%------------------------------------------------------

Eph = zeros(29,32);

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
%  Eph      --> double, [29,32]
fid_eph = fopen([filerootOUT '_eph_00.bin'],'w+');

%nmea sentences
fid_nmea = fopen([filerootOUT '_rover_NMEA.txt'],'wt');

%------------------------------------------------------
% initialization
%------------------------------------------------------

%ionosphere parameters
iono = zeros(8,1);

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

% disable subframe buffer output
fprintf('Disabling u-blox receiver subframe buffer (SFRB) messages...\n');

reply_SFRB = ublox_CFG_MSG(rover, 'RXM', 'SFRB', 0);
tries = 0;

while (~reply_SFRB)
    tries = tries + 1;
    if (tries > 3)
        disp('It was not possible to disable SFRB messages.');
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
    reply_SFRB = ublox_CFG_MSG(rover, 'RXM', 'SFRB', 0);
end

% enable GGA messages, disable all other NMEA messages
fprintf('Configuring u-blox receiver NMEA messages:\n');

ublox_CFG_MSG(rover, 'NMEA', 'GGA', 1); fprintf('Enabling GGA...\n');
ublox_CFG_MSG(rover, 'NMEA', 'GLL', 0); fprintf('Disabling GLL ');
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

%first poll
ublox_poll_message(rover, 'AID', 'EPH', 0);
pause(0.1);
ublox_poll_message(rover, 'AID', 'HUI', 0);
%poll flags
eph_polled = 1;
hui_polled = 1;

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

        [cell_rover, nmea_sentences] = decode_ublox(data_rover);

        %read data type
        type = '';

        %data type counters
        nRAW = 0;
        nEPH = 0;
        nHUI = 0;

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
                
                %satellites with ephemerides available
                satEph = find(sum(abs(Eph))~=0);
                
                %delete data if ephemerides are not available
                delsat = setdiff(1:32,satEph);
                pr_R(delsat,1)  = 0;
                ph_R(delsat,1)  = 0;
                dop_R(delsat,1) = 0;
                snr_R(delsat,1) = 0;
                
                %satellites with observations available
                satObs = find(pr_R(:,1) ~= 0);
                
                %if all the visible satellites ephemerides have been transmitted
                %and the total number of satellites is >= 4
                if (ismember(satObs,satEph)) & (length(satObs) >= 4)
                    
                    %data save
                    fwrite(fid_obs, [0; 0; time_R; week_R; zeros(32,1); pr_R; zeros(32,1); ph_R; dop_R; zeros(32,1); snr_R; zeros(3,1); iono(:,1)], 'double');
                    fwrite(fid_eph, [0; Eph(:)], 'double');
                end

                type = [type 'RXM-RAW '];
                nRAW = nRAW + 1;

            %AID-HUI message data save
            elseif (strcmp(cell_rover{1,i},'AID-HUI'))
                
                %ionosphere parameters
                iono(:, 1) = cell_rover{3,i}(9:16);
                
                if (nHUI == 0)
                    type = [type 'AID-HUI '];
                end
                nHUI = nHUI + 1;
                
            %AID-EPH message
            elseif (strcmp(cell_rover{1,i},'AID-EPH'))

                %satellite number
                sat = cell_rover{2,i}(1);

                if (~isempty(sat) & sat > 0)
                    Eph(:, sat) = cell_rover{2,i}(:);
                end

                if (nEPH == 0)
                    type = [type 'AID-EPH '];
                end
                nEPH = nEPH + 1;

            end

        end

        if (~isempty(nmea_sentences))
            n = size(nmea_sentences,1);
            for i = 1 : n
                fprintf(fid_nmea, '%s', char(nmea_sentences(i,1)));
            end

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
        
        %visualization (AID-HUI information)
        if (nHUI > 0)
            fprintf('Ionosphere parameters: ');
            if (sum(iono) ~= 0)
                fprintf('\n');
                fprintf('    alpha0: %12.4E\n', iono(1));
                fprintf('    alpha1: %12.4E\n', iono(2));
                fprintf('    alpha2: %12.4E\n', iono(3));
                fprintf('    alpha3: %12.4E\n', iono(4));
                fprintf('    beta0 : %12.4E\n', iono(5));
                fprintf('    beta1 : %12.4E\n', iono(6));
                fprintf('    beta2 : %12.4E\n', iono(7));
                fprintf('    beta3 : %12.4E\n', iono(8));
            else
                fprintf('not sent\n');
            end
        end

        %visualization (AID-EPH information)
        if (nEPH > 0)
            sat = find(sum(abs(Eph))>0);
            fprintf('Eph: ');
            for i = 1 : length(sat)
                fprintf('%d ', sat(i));
            end
            fprintf('\n');
        end

        %poll a new AID-EPH message every 10 epochs
        if (mod(current_time-start_time,10) < 1)
            if (eph_polled == 0)
                ublox_poll_message(rover, 'AID', 'EPH', 0);
                eph_polled = 1;
            end
        else
            eph_polled = 0;
        end
        
        %wait for asynchronous write to finish
        pause(0.1);
        
        %poll a new AID-HUI message every 60 epochs
        if (mod(current_time-start_time,60) < 1)
            if (hui_polled == 0)
                ublox_poll_message(rover, 'AID', 'HUI', 0);
                hui_polled = 1;
            end
        else
            hui_polled = 0;
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