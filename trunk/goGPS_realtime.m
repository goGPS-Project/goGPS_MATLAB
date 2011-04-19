function goGPS_realtime(filerootOUT, protocol, mode_vinc, flag_ms, flag_ge, flag_cov, flag_NTRIP, flag_ms_pos, flag_skyplot, flag_plotproc, ref_path, mat_path, pos_M, dop1_M, pr2_M, pr2_R, ph2_M, ph2_R, dop2_M, dop2_R)

% SYNTAX:
%   goGPS_realtime(filerootOUT, protocol, mode_vinc, flag_ms, flag_ge, 
%   flag_cov, flag_NTRIP, flag_ms_pos, flag_skyplot, ref_path, mat_path,
%   pos_M, dop1_M, pr2_M, pr2_R, ph2_M, ph2_R, dop2_M, dop2_R);
%
% INPUT:
%   filerootOUT = output file prefix
%   protocol    = protocol (e.g. 0-Ublox, 1-Fastrax, 2-SkyTraq)
%   mode_vinc   = constraint flag
%   flag_ms     = plot master station flag
%   flag_ge     =  google earth flag
%   flag_cov    = plot error ellipse flag
%   flag_NTRIP  = use/don't use NTRIP flag
%   flag_ms_pos = use/don't use RTCM master position
%   flag_skyplot  = use/don't use CPU saving mode (no skyplot, no SNR graph)
%   flag_plotproc = display/don't display plot figure
%   ref_path    = reference path
%   mat_path    = reference path adjacency matrix
%   pos_M  = master station position (X,Y,Z)
%   dop1_M = Doppler observation MASTER-SATELLITE (carrier L1)
%   pr2_M  = code pseudorange MASTER-SATELLITE (carrier L2)
%   pr2_R  = code pseudorange ROVER-SATELLITE (carrier L2)
%   ph2_M  = phase observation MASTER-SATELLITE (carrier L2)
%   ph2_R  = phase observation ROVER-SATELLITE (carrier L2)
%   dop2_M = Doppler observation MASTER-SATELLITE (carrier L2)
%   dop2_R = Doppler observation ROVER-SATELLITE (carrier L2)
%
% DESCRIPTION:
%   goGPS real-time algorithm: stream reading and synchronization,
%   positioning, plotting, output data saving.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Ivan Reguzzoni
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

%global lambda1
global o1 o2 nN
global COMportR master_ip master_port server_delay
global nmea_init nmea_update_rate hui_poll_rate
global azR elR distR azM elM distM
global PDOP HDOP VDOP KPDOP KHDOP KVDOP
global Xhat_t_t Cee conf_sat conf_cs pivot Yhat_t_t
global master rover

%------------------------------------------------------
% read protocol parameters
%------------------------------------------------------
% u-blox
if (protocol == 0)
    prot_par = param_ublox;
elseif (protocol == 1)
    prot_par = param_fastrax;
elseif (protocol == 2)
    prot_par = param_skytraq;
end

%------------------------------------------------------
% data file creation
%------------------------------------------------------

%master binary stream (uint8)
fid_master = fopen([filerootOUT '_master_00.bin'],'w+');

%rover binary stream (uint8)
fid_rover = fopen([filerootOUT '_rover_00.bin'],'w+');

%input observations (master & rover)
%  time_GPS --> double, [1,1]
%  time_M   --> double, [1,1]
%  time_R   --> double, [1,1]
%  pr_M     --> double, [32,1]
%  pr_R     --> double, [32,1]
%  ph_M     --> double, [32,1]
%  ph_R     --> double, [32,1]
%  dop_R    --> double, [32,1]
%  snr_M    --> double, [32,1]
%  snr_R    --> double, [32,1]
%  XM       --> double, [1,1]
%  YM       --> double, [1,1]
%  ZM       --> double, [1,1]
fid_obs = fopen([filerootOUT '_obs_00.bin'],'w+');

%input ephemerides
%  time_GPS --> double, [1,1]
%  Eph      --> double, [29,32]
fid_eph = fopen([filerootOUT '_eph_00.bin'],'w+');

%Kalman filter output data
%free motion
%  Xhat_t_t --> double, [o3+nN,1]
%  Cee      --> double, [o3+nN,o3+nN]
%constrained motion
%  Xhat_t_t --> double, [o1+nN,1]
%  Yhat_t_t --> double, [3,1]
%  Cee      --> double, [o1+nN,o1+nN]
fid_kal = fopen([filerootOUT '_kal_00.bin'],'w+');

%satellite azimuth, elevation and distance
%  azM   --> double, [32,1]
%  azR   --> double, [32,1]
%  elM   --> double, [32,1]
%  elR   --> double, [32,1]
%  distM --> double, [32,1]
%  distR --> double, [32,1]
fid_sat = fopen([filerootOUT '_sat_00.bin'],'w+');

%dilution of precision
%  PDOP     --> double, [1,1]
%  HDOP     --> double, [1,1]
%  VDOP     --> double, [1,1]
%  KPDOP    --> double, [1,1]
%  KHDOP    --> double, [1,1]
%  KVDOP    --> double, [1,1]
fid_dop = fopen([filerootOUT '_dop_00.bin'],'w+');

%satellite configuration
%  conf_sat --> int8, [32,1]
%  conf_cs  --> int8, [32,1]
%  pivot    --> int8, [1,1]
fid_conf = fopen([filerootOUT '_conf_00.bin'],'w+');

%nmea sentences 
fid_nmea = fopen([filerootOUT '_', prot_par{1,1},'_NMEA.txt'],'wt');

%"file hour" variable
hour = 0;

%------------------------------------------------------
% initialization
%------------------------------------------------------

%number of unknown phase ambiguities
nN = 32;

%ionosphere parameters
iono = zeros(8,1);

%------------------------------------------------------
% creation of the connection to the ROVER (u-blox)
%------------------------------------------------------

% find a serial port object.
obj1 = instrfind('Type', 'serial', 'Port', COMportR, 'Tag', '');

% if a serial object already exists, delete it before creating a new one
if ~isempty(obj1)
    delete(obj1);
end

% serial object creation
rover = serial (COMportR,'BaudRate',prot_par{2,1});
set(rover,'InputBufferSize',prot_par{3,1});
set(rover,'FlowControl','hardware');
set(rover,'RequestToSend','on');
fopen(rover);

if (protocol == 0)
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
        rover = serial (COMportR,'BaudRate',prot_par{2,1});
        set(rover,'InputBufferSize',prot_par{3,1});
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
        rover = serial (COMportR,'BaudRate',prot_par{2,1});
        set(rover,'InputBufferSize',prot_par{3,1});
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
        rover = serial (COMportR,'BaudRate',prot_par{2,1});
        set(rover,'InputBufferSize',prot_par{3,1});
        set(rover,'FlowControl','hardware');
        set(rover,'RequestToSend','on');
        fopen(rover);
        reply_RAW = ublox_CFG_MSG(rover, 'RXM', 'RAW', 1);
    end

    % disable NMEA messages
    fprintf('Disabling u-blox receiver NMEA messages...\n');

    % ublox_CFG_MSG(rover, 'NMEA', 'GGA', 0); fprintf('GGA ');
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
end

%------------------------------------------------------
% absolute time start
%------------------------------------------------------

tic

%log file initialization
delete([filerootOUT '_log.txt']);
diary([filerootOUT '_log.txt']);
diary on

%------------------------------------------------------
% rover header package acquisition
%------------------------------------------------------

%visualization
fprintf('\n');
fprintf('ROVER LOCK-PHASE (HEADER PACKAGE)\n');

%initialization
rover_1 = 0;
rover_2 = 0;

%starting epoch determination
while (rover_1 ~= rover_2) || (rover_1 == 0) || (rover_1 < prot_par{4,1})

    %starting time
    current_time = toc;

    %serial port check
    rover_1 = get(rover,'BytesAvailable');
    pause(0.05);
    rover_2 = get(rover,'BytesAvailable');
    
    %visualization
    fprintf([prot_par{1,1},': %7.4f sec (%4d bytes --> %4d bytes)\n'], current_time, rover_1, rover_2);

end

%empty serial port (data not decoded)
data_rover = fread(rover,rover_1,'uint8'); %#ok<NASGU>

%-----------------------------------------------------------
% rover initial positioning (stand-alone)
%-----------------------------------------------------------

%visualization
fprintf('\n');
fprintf('ROVER POSITIONING (STAND-ALONE)...\n');
fprintf('note: it might take some time to acquire signal from 4 satellites\n');

%pseudoranges
pr_R = zeros(32,1);
%ephemerides
Eph = zeros(29,32);
%satellites with observations available
satObs = [];
nsatObs_old = [];
%satellites with ephemerides available
satEph = [];

while ((length(satObs) < 4) | (~ismember(satObs,satEph)))

    if (protocol == 0)
        %poll available ephemerides
        ublox_poll_message(rover, 'AID', 'EPH', 0);
        
        %wait for asynchronous write to finish
        pause(0.1);

        %poll AID-HUI message (sat. Health / UTC / Ionosphere)
        ublox_poll_message(rover, 'AID', 'HUI', 0);

        %wait for asynchronous write to finish
        pause(0.1);
    end

    %initialization
    rover_1 = 0;
    rover_2 = 0;

    %starting epoch determination
    while (rover_1 ~= rover_2) || (rover_1 == 0)

        %starting time
        current_time = toc;

        %serial port check
        rover_1 = get(rover,'BytesAvailable');
        pause(0.1);
        rover_2 = get(rover,'BytesAvailable');
    end

    data_rover = fread(rover,rover_1,'uint8');     %serial port reading
    fwrite(fid_rover,data_rover,'uint8');          %transmitted stream save
    data_rover = dec2bin(data_rover,8);            %conversion to binary (N x 8bit matrix)
    data_rover = data_rover';                      %transpose (8bit x N matrix)
    data_rover = data_rover(:)';                   %conversion to string (8N bit vector)

    % -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (protocol == 0)
        [cell_rover] = decode_ublox(data_rover);
    elseif (protocol == 1)
        [cell_rover] = decode_fastrax_it03(data_rover);
    elseif (protocol == 2)
        [cell_rover] = decode_skytraq(data_rover);
    end
	
	%for SkyTraq
	IOD_time = -1;
    
    %cell_rover
    for i = 1 : size(cell_rover,2)

        %Timing/raw message data save (RXM-RAW | PSEUDO)
        if (strcmp(cell_rover{1,i},prot_par{1,2}))

            %just information needed for basic positioning is saved
            time_GPS  = round(cell_rover{2,i}(1));
            pr_R(:,1) = cell_rover{3,i}(:,2);
        
        %Timing message data save (MEAS_TIME)
        elseif (strcmp(cell_rover{1,i},prot_par{4,2}))
            
            IOD_time = cell_rover{2,i}(1);
            time_GPS = cell_rover{2,i}(3);
            
        %Raw message data save (RAW_MEAS)
        elseif (strcmp(cell_rover{1,i},prot_par{1,2}))
            
            IOD_raw = cell_rover{2,i}(1);
            if (IOD_raw == IOD_time)
                pr_R = cell_rover{3,i}(:,3);
            end
            
        %Eph message data save (AID-EPH | FTX-EPH)
        elseif (strcmp(cell_rover{1,i},prot_par{2,2}))

            %satellite number
            sat = cell_rover{2,i}(1);

            if (~isempty(sat) & sat > 0)
                Eph(:, sat) = cell_rover{2,i}(:);
            end

            
        %Hui message data save (AID-HUI)
        elseif (strcmp(cell_rover{1,i},prot_par{3,2}))

            %ionosphere parameters
            iono(:, 1) = cell_rover{3,i}(9:16);

        end
    end

    %satellites with ephemerides available
    satEph = find(sum(abs(Eph))~=0);

    %delete data if ephemerides are not available
    delsat = setdiff(1:32,satEph);
    pr_R(delsat)  = 0;

    %satellites with observations available
    satObs = find(pr_R ~= 0);
    
    if (isempty(nsatObs_old) | length(satObs) ~= nsatObs_old)
        %display current number of satellites
        fprintf('Number of visible satellites with ephemerides: %d\n', length(satObs));
        nsatObs_old = length(satObs);
    end
end

%positioning by Bancroft algorithm
[pos_R] = input_bancroft(pr_R(satObs,1), satObs, time_GPS, Eph);

fprintf('ROVER approximate position computed using %d satellites\n', sum(pr_R ~= 0));

%NMEA sentence with initial approximate position
nmea_init = NMEA_GGA_gen([pos_R(1) pos_R(2) pos_R(3)],10);

%------------------------------------------------------------
% acquisition of the next rover message (for synchronization)
%------------------------------------------------------------

%visualization
fprintf('\n');
fprintf('ROVER SYNCHRONIZATION...\n');

%initialization
rover_1 = 0;
rover_2 = 0;
sync_rover = 0;

while (~sync_rover)
    
    %starting epoch determination
    while (rover_1 ~= rover_2) || (rover_1 == 0) || (rover_1 < prot_par{4,1})
        
        %starting time
        current_time = toc;
        
        %serial port check
        rover_1 = get(rover,'BytesAvailable');
        pause(0.1);
        rover_2 = get(rover,'BytesAvailable');
        
        %visualization
        fprintf([prot_par{1,1},': %7.4f sec (%4d bytes --> %4d bytes)\n'], current_time, rover_1, rover_2);
        
    end
    
    data_rover = fread(rover,rover_1,'uint8');     %serial port reading
    data_rover = dec2bin(data_rover,8);            %conversion to binary (N x 8bit matrix)
    data_rover = data_rover';                      %transpose (8bit x N matrix)
    data_rover = data_rover(:)';                   %conversion to string (8N bit vector)
    
    %message decoding
    if (protocol == 0)
        [cell_rover] = decode_ublox(data_rover);
    elseif (protocol == 1)
        [cell_rover] = decode_fastrax_it03(data_rover);
    elseif (protocol == 2)
        [cell_rover] = decode_skytraq(data_rover);
    end
    
    for i = 1 : size(cell_rover,2)
        
        %RXM-RAW message data save
        if (strcmp(cell_rover{1,i},prot_par{1,2}))
            
            %just information about the epoch is saved
            time_GPS = round(cell_rover{2,i}(1));
            week_GPS = cell_rover{2,i}(2);
            
            sync_rover = 1;
            
        %Timing message data save (MEAS_TIME)
        elseif (strcmp(cell_rover{1,i},prot_par{4,2}))

            time_GPS = cell_rover{2,i}(3);
            
            sync_rover = 1;
        end
    end
end

%starting time is set
safety_lag = 0.1;                       %safety lag for reading ROVER data
start_time = current_time-safety_lag;   %starting time

%------------------------------------------------------
% creation of the connection to the MASTER
%------------------------------------------------------

ntripstring = NTRIP_string_generator(nmea_init);

master = tcpip(master_ip,master_port);
set(master,'InputBufferSize', 16384);
fopen(master);
fwrite(master,ntripstring);

%wait until the buffer is written before continuing
while get(master,'BytesAvailable') == 0, end

%--------------------------------------------------------
% acquisition of the 1st master message (dropped)
%--------------------------------------------------------

%visualization
fprintf('\n');
fprintf('MASTER INITIALIZATION\n');

%wait until ALL the master packages have arrived
while (current_time-start_time < 0.9)
    current_time = toc;
end

%check of the TCP-IP port
master_1 = get(master,'BytesAvailable');
pause(0.05);
master_2 = get(master,'BytesAvailable');

%visualization
fprintf('master: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time-start_time, master_1, master_2);

%empty TCP-IP port (data not decoded)
if (master_1 == master_2) & (master_1 ~= 0)
    data_master = fread(master,master_1,'uint8'); %#ok<NASGU>
end

%the master stopped!
if (master_1 == master_2) & (master_1 == 0)
    fclose(master);
    fopen(master);
end

%--------------------------------------------------------
% buffer settings
%--------------------------------------------------------

%buffer current position
b = 1;

%buffer dimension
B = 20;

%buffer initialization
tick_M = zeros(B,1);      % empty/full master buffer
tick_R = zeros(B,1);      % empty/full rover buffer
time_M = zeros(B,1);      % master time buffer
time_R = zeros(B,1);      % rover time buffer
week_R = zeros(B,1);      % rover week buffer
pr_M   = zeros(32,B);     % master code buffer
pr_R   = zeros(32,B);     % rover code buffer
ph_M   = zeros(32,B);     % master phase buffer
ph_R   = zeros(32,B);     % rover phase buffer
dop_R  = zeros(32,B);     % rover Doppler buffer
snr_M  = zeros(32,B);     % master SNR buffer
snr_R  = zeros(32,B);     % rover SNR buffer
if (flag_ms_pos)
    pos_M  = zeros(3, B);        % master station coordinates read from RTCM
else
    for i = 2 : B
        pos_M(:,i) = pos_M(:,1); % master station coordinates set manually
    end
end

%--------------------------------------------------------
% master position update (VRS) management
%--------------------------------------------------------

%master position update expected (i.e. a NMEA string was sent)
master_update = 1;

%master position awaiting indexing (i.e. time tag from an observation)
master_waiting = 0;

%--------------------------------------------------------
% figure management
%--------------------------------------------------------

%loop control initialization
if (flag_plotproc)
    h1 = uicontrol(gcf, 'style', 'pushbutton', 'position', [10 10 40 20], 'string', 'STOP', ...
        'callback', 'setappdata(gcf, ''run'', 0)');
    set(gcf, 'name', 'goGPS', 'toolbar', 'figure');
    flag = 1;
    setappdata(gcf, 'run', flag);
else
    f1 = figure;
    s1 = get(0,'ScreenSize');
    set(f1, 'position', [s1(3)-240-20 s1(4)-80-40 240 80], 'menubar', 'none', 'name', 'Navigation');
    h1 = uicontrol(gcf, 'style', 'pushbutton', 'position', [80 20 80 40], 'string', 'STOP', ...
        'callback', 'setappdata(gcf, ''run'', 0)');
    flag = 1;
    setappdata(gcf, 'run', flag);
end

%--------------------------------------------------------
% start time synchronization
%--------------------------------------------------------

%go to the subsequent epoch(s)
dtime = ceil(current_time-start_time);
while (current_time-start_time < dtime)
    current_time = toc;
end

%DEBUG tick(0) bug
if (dtime - 1) > 1
    fprintf('WARNING! Master connection delay=%d sec\n', dtime - 1);
end

%GPS epoch increment
time_GPS = time_GPS + dtime;

%starting time re-initialization
start_time = start_time + dtime - 1;

%--------------------------------------------------------
% master/rover data acquisition and position computation
%--------------------------------------------------------

%counter initialization
t = 1;

%time increment initialization (default 1 sec)
dtime = 1;

%infinite loop
while flag

    %visualization
    fprintf('\n');
    fprintf('-----------------------------------------------\n');
    fprintf('TIMING\n');

    %visualization
    fprintf('epoch %d: GPStime=%d:%d\n', t, week_GPS, time_GPS);

    %-------------------------------------
    % file management
    %-------------------------------------

    if (floor(t/3600) > hour)

        hour = floor(t/3600);
        hour_str = num2str(hour,'%02d');

        fclose(fid_master);
        fclose(fid_rover);
        fclose(fid_obs);
        fclose(fid_eph);
        fclose(fid_kal);
        fclose(fid_sat);
        fclose(fid_dop);
        fclose(fid_conf);
        
        fid_master = fopen([filerootOUT '_master_' hour_str '.bin'],'w+');
        fid_rover  = fopen([filerootOUT '_rover_'  hour_str '.bin'],'w+');
        fid_obs    = fopen([filerootOUT '_obs_'    hour_str '.bin'],'w+');
        fid_eph    = fopen([filerootOUT '_eph_'    hour_str '.bin'],'w+');
        fid_kal    = fopen([filerootOUT '_kal_'    hour_str '.bin'],'w+');
        fid_sat    = fopen([filerootOUT '_sat_'    hour_str '.bin'],'w+');
        fid_dop    = fopen([filerootOUT '_dop_'    hour_str '.bin'],'w+');
        fid_conf   = fopen([filerootOUT '_conf_'   hour_str '.bin'],'w+');

    end

    %-------------------------------------
    % rover data
    %-------------------------------------

    %visualization
    fprintf('\n');
    fprintf('ROVER DATA\n');

    %time acquisition
    current_time = toc;
    step_time = round(current_time-start_time);

    %initialization
    %rover_1 = 0;
    %rover_2 = 0;

    %initialization
    rover_init = get(rover,'BytesAvailable');   % if previous data are present
    rover_1 = rover_init;
    rover_2 = rover_init;

    % maximum waiting time for the rover
    dtMax_rover = 0.2;

    %multiple condition: while (package not available) AND (waiting time not expired)
    %while ((rover_1 ~= rover_2) | (rover_1 == 0)) & (current_time-start_time-step_time < dtMax_rover)
    while ((rover_1 ~= rover_2) | (rover_1 == rover_init)) & (current_time-start_time-step_time < dtMax_rover)

        %time acquisition
        current_time = toc;

        %serial port check
        rover_1 = get(rover,'BytesAvailable');
        pause(0.05);
        rover_2 = get(rover,'BytesAvailable');

    end

    %visualization
    fprintf([prot_par{1,1},': %7.4f sec (%4d bytes --> %4d bytes)\n'], current_time-start_time, rover_1, rover_2);

    %-------------------------------------

    if (dtime < B)

        %shift of the rover buffers
        tick_R(1+dtime:end)  = tick_R(1:end-dtime);
        time_R(1+dtime:end)  = time_R(1:end-dtime);
        week_R(1+dtime:end)  = week_R(1:end-dtime);
        pr_R(:,1+dtime:end)  = pr_R(:,1:end-dtime);
        ph_R(:,1+dtime:end)  = ph_R(:,1:end-dtime);
        dop_R(:,1+dtime:end) = dop_R(:,1:end-dtime);
        snr_R(:,1+dtime:end) = snr_R(:,1:end-dtime);

        %current cell to zero
        tick_R(1:dtime)  = zeros(dtime,1);
        time_R(1:dtime)  = zeros(dtime,1);
        week_R(1:dtime)  = zeros(dtime,1);
        pr_R(:,1:dtime)  = zeros(32,dtime);
        ph_R(:,1:dtime)  = zeros(32,dtime);
        dop_R(:,1:dtime) = zeros(32,dtime);
        snr_R(:,1:dtime) = zeros(32,dtime);

    else

        %buffer to zero
        tick_R = zeros(B,1);
        time_R = zeros(B,1);
        week_R = zeros(B,1);
        pr_R   = zeros(32,B);
        ph_R   = zeros(32,B);
        dop_R  = zeros(32,B);
        snr_R  = zeros(32,B);

    end

    %-------------------------------------

    %read message type
    type = '';

    %check if the writing of the package has ended
    if (rover_1 == rover_2) & (rover_1 ~= 0)

        data_rover = fread(rover,rover_1,'uint8');     %serial port reading
        fwrite(fid_rover,data_rover,'uint8');          %transmitted stream save
        data_rover = dec2bin(data_rover,8);            %conversion to binary (N x 8bit matrix)
        data_rover = data_rover';                      %transpose (8bit x N matrix)
        data_rover = data_rover(:)';                   %conversion to string (8N bit vector)
        %dep_rover = strcat(dep_rover,data_rover);

        %message decoding
        if (protocol == 0)
            [cell_rover, nmea_sentences] = decode_ublox(data_rover);
        elseif (protocol == 1)
            %-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            [cell_rover] = decode_fastrax_it03(data_rover);
            nmea_sentences = [];
        end

        %data type counters
        nEPH = 0;
        nHUI = 0;
        
        %for Fastrax
        tick_TRACK = 0;
        
        %for SkyTraq
        IOD_time = -1;

        for i = 1 : size(cell_rover,2)
            
            % TRACK message data save (only for Fastrax)
            if (strcmp(cell_rover{1,i},prot_par{6,2}))

                tick_TRACK    = cell_rover{2,i}(1);
                phase_TRACK   = cell_rover{3,i}(:,6);
                
                type = [type prot_par{6,2} ' '];

            %RXM-RAW message data save
            elseif (strcmp(cell_rover{1,i},prot_par{1,2}))

                %buffer index computation
                index = time_GPS - round(cell_rover{2,i}(1)) + 1;
                
                if (index <= B)
                    while (index < 1)
                        time_GPS = time_GPS + 1;
                        index = time_GPS - round(cell_rover{2,i}(1)) + 1;
                    end

                    %buffer writing
                    tick_R(index)  = 1;
                    time_R(index)  = round(cell_rover{2,i}(1));
                    week_R(index)  = cell_rover{2,i}(2);
                    pr_R(:,index)  = cell_rover{3,i}(:,2);
                    ph_R(:,index)  = cell_rover{3,i}(:,1);
                    dop_R(:,index) = cell_rover{3,i}(:,3);
                    snr_R(:,index) = cell_rover{3,i}(:,6);

                    %phase computation (only for Fastrax)
                    tick_PSEUDO = cell_rover{2,i}(4);
                    if (protocol == 1 & tick_TRACK == tick_PSEUDO)
                        ph_R(abs(pr_R(:,index)) > 0) = phase_TRACK(abs(pr_R(:,index)) > 0);
                    end

                    %manage "nearly null" data
                    pos = abs(ph_R(:,index)) < 1e-100;
                    ph_R(pos,index) = 0;

%                     %phase adjustement
%                     pos = abs(ph_R(:,index)) > 0 & abs(ph_R(:,index)) < 1e7;
%                     ambig = 2^23;
%                     n = floor( (pr_R(pos,index)/lambda1-ph_R(pos,index)) / ambig + 0.5 );
%                     ph_R(pos,index) = ph_R(pos,index) + n*ambig;

                    type = [type prot_par{1,2} ' '];
                end
                
            %Timing message data save (MEAS_TIME)
            elseif (strcmp(cell_rover{1,i},prot_par{4,2}))
                
                IOD_time = cell_rover{2,i}(1);
                time_stq = cell_rover{2,i}(3);
                week_stq = cell_rover{2,i}(2);
                
                type = [type prot_par{4,2} ' '];
                
            %Raw message data save (RAW_MEAS)
            elseif (strcmp(cell_rover{1,i},prot_par{1,2}))
                
                IOD_raw = cell_rover{2,i}(1);
                if (IOD_raw == IOD_time)

                    %buffer index computation
                    index = time_GPS - time_stq + 1;
                    
                    if (index <= B)
                        while (index < 1)
                            time_GPS = time_GPS + 1;
                            index = time_GPS - time_stq + 1;
                        end
                        
                        %buffer writing
                        tick_R(index)  = 1;
                        time_R(index)  = round(time_stq);
                        week_R(index)  = week_stq;
                        pr_R(:,index)  = cell_rover{3,i}(:,3);
                        ph_R(:,index)  = cell_rover{3,i}(:,4);
                        snr_R(:,index) = cell_rover{3,i}(:,2);
                        dop_R(:,index) = cell_rover{3,i}(:,5);
                        
                        %manage "nearly null" data
                        pos = abs(ph_R(:,index)) < 1e-100;
                        ph_R(pos,index) = 0;
                        
                        type = [type prot_par{1,2} ' '];
                    end
                end

            %Eph message data save
            elseif (strcmp(cell_rover{1,i},prot_par{2,2}))

                %satellite number
                sat = cell_rover{2,i}(1);

                if (~isempty(sat) & sat > 0)
                    Eph(:, sat) = cell_rover{2,i}(:);
                end

                if (nEPH == 0)
                    type = [type prot_par{2,2} ' '];
                end

                nEPH = nEPH + 1;
                
            %Hui message data save
            elseif (strcmp(cell_rover{1,i},prot_par{3,2}))
                
                %ionosphere parameters
                iono(:, 1) = cell_rover{3,i}(9:16);
                
                if (nHUI == 0)
                    type = [type prot_par{3,2} ' '];
                end
                nHUI = nHUI + 1;

            end
        end

        %NMEA data save
        if (~isempty(nmea_sentences))
            n = size(nmea_sentences,1);
            for i = 1 : n
                fprintf(fid_nmea, '%s', char(nmea_sentences(i,1)));
            end

            type = [type 'NMEA '];
        end
    end

    %time acquisition (at the end of the rover decoding)
    current_time = toc;

    %-------------------------------------

    %visualization
    i = min(b,B);                        %pointer to the buffer or last buffer cell
    sat_pr = find(pr_R(:,i) ~= 0);       %satellites with code available
    sat_ph = find(ph_R(:,i) ~= 0);       %satellites with phase available
    sat = union(sat_pr,sat_ph);          %satellites with code or phase available

    fprintf('decoding: %7.4f sec (%smessages)\n', current_time-start_time, type);
    fprintf('GPStime=%7.4f (%d satellites)\n', time_R(i), length(sat));

    fprintf('P1 SAT:');
    for j = 1 : length(sat_pr)
        fprintf(' %02d', sat_pr(j));
    end
    fprintf('\n');

    fprintf('L1 SAT:');
    k = 1;
    for j = 1 : length(sat_ph)
        while (sat_ph(j) ~= sat_pr(k))
            fprintf('   ');
            k = k + 1;
        end
        fprintf(' %02d', sat_ph(j));
        k = k + 1;
    end
    fprintf('\n');

    %--------------------------------------------------------------
    %ephemerides request
    %--------------------------------------------------------------

    if (~isempty(sat) & index > 0)
        %satellites with observations available for ephemerides polling
        conf_sat_eph = zeros(32,1);
        conf_sat_eph(sat_pr) = 1;

        %ephemerides update cycle
        conf_eph = (sum(abs(Eph),1) == 0);

        [null, sat_index] = sort(snr_R(:, index),1,'descend'); %#ok<ASGLU>

        conf_sat_eph = conf_sat_eph(sat_index);
        conf_eph = conf_eph(sat_index);

        check = 0;
        i = 1;

        while ((check == 0) & (i<=32))

            s = sat_index(i);

            %if satellite i is available
            if (abs(conf_sat_eph(i)) == 1)

                %time from the ephemerides reference epoch
                if (conf_eph(i) == 0)
                    toe = Eph(18,s);
                    tk = check_t(time_GPS-toe);
                end

                %if ephemeris i is not present OR ephemeris i is too old
                if (conf_eph(i) == 1) | (tk > 3600)
                    if (protocol == 0)
                        ublox_poll_message(rover, 'AID', 'EPH', 1, dec2hex(s,2));
                        fprintf('Satellite %d ephemeris polled\n', s);
                    end
                    % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    check = 1;
                end
            end
            i = i + 1;
        end
    end

    %-------------------------------------
    % master data
    %-------------------------------------

    %visualization
    fprintf('\n');
    fprintf('MASTER DATA\n');

    %time acquisition
    current_time = toc;

    %output data initialization
    cell_master = cell(0);

    %test condition initialization
    test_master = 0;

    %maximum master waiting time
    dtMax_master = 0.8;

    %multiple condition: while (I have not received the 19/1002/1004 message for the time_GPS epoch) AND (time is not expired)
    while ((test_master == 0) && (current_time-start_time-step_time < dtMax_master))

        %time acquisition
        current_time = toc;

        %TCP-IP port check
        master_1 = get(master,'BytesAvailable');
        pause(server_delay);
        master_2 = get(master,'BytesAvailable');
                
        %check if package writing is finished
        if (master_1 == master_2) & (master_1 ~= 0)

            data_master = fread(master,master_1,'uint8');     %TCP-IP port reading
            fwrite(fid_master,data_master,'uint8');           %transmitted stream save
            data_master = dec2bin(data_master,8)';            %conversion to binary (N x 8bit matrix)
            data_master = data_master(:);                     %transpose (8bit x N matrix)
            data_master = data_master';                       %conversion to string (8N bit vector)
            %dep_master = strcat(dep_master,data_master);

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
                cell_master = [cell_master decode_rtcm2(sixofeight,time_GPS)]; %RTCM 2 decoding
            else
                cell_master = [cell_master decode_rtcm3(data_master)];         %RTCM 3 decoding and appending
            end
        end

        %detect the last read 19/1002/1004 message
        i = size(cell_master,2);
        while (i > 0) & (isempty(cell_master{1,i}) | ((cell_master{1,i} ~= 19) & (cell_master{1,i} ~= 1002) & (cell_master{1,i} ~= 1004)))
            i = i - 1;
        end

        %check the exit condition
        if (i > 0) && (round(cell_master{2,i}(2)) == time_GPS)
            test_master = 1;
        end

    end

    %visualization 
    fprintf('master: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time-start_time, master_1, master_2);

    %-------------------------------------

    if (dtime < B)

        %shift of the master buffers
        tick_M(1+dtime:end)  = tick_M(1:end-dtime);
        time_M(1+dtime:end)  = time_M(1:end-dtime);
        pr_M(:,1+dtime:end)  = pr_M(:,1:end-dtime);
        ph_M(:,1+dtime:end)  = ph_M(:,1:end-dtime);
        snr_M(:,1+dtime:end) = snr_M(:,1:end-dtime);
        pos_M(:,1+dtime:end) = pos_M(:,1:end-dtime);
        
        %current cell to zero
        tick_M(1:dtime)  = zeros(dtime,1);
        time_M(1:dtime)  = zeros(dtime,1);
        pr_M(:,1:dtime)  = zeros(32,dtime);
        ph_M(:,1:dtime)  = zeros(32,dtime);
        snr_M(:,1:dtime) = zeros(32,dtime);
        %pos_M current cell keeps the latest value(s), until it is updated
        % by a new RTCM message (3, 1005 or 1006)
        pos = find(sum(pos_M(:,1+dtime:end)) ~= 0);
        if (~isempty(pos))
           % pos = pos(3) / 3;
            pos = pos(1);
            pos_M(1,1:dtime) = pos_M(1,pos+dtime);
            pos_M(2,1:dtime) = pos_M(2,pos+dtime);
            pos_M(3,1:dtime) = pos_M(3,pos+dtime);
        end
    else

        %buffer to zero
        tick_M = zeros(B,1);
        time_M = zeros(B,1);
        pr_M   = zeros(32,B);
        ph_M   = zeros(32,B);
        snr_M  = zeros(32,B);
        if (flag_ms_pos)
            % master station coordinates read from RTCM
            pos_M  = zeros(3, B);
        else
            % master station coordinates set manually
            pos_M  = [pos_M(:,1) zeros(3, B-1)];
        end

    end

    %-------------------------------------

    %read message type
    type = '';

    index_ph = [];

    for i = 1 : size(cell_master,2)

        if (~isempty(cell_master{1,i}))
            switch cell_master{1,i}

                %message 18 (RTCM2)
                case 18

                    %buffer index computation
                    index = time_GPS - round(cell_master{2,i}(2)) + 1;

                    index_ph = [index_ph index];

                    if (index <= B)

                        %%buffer writing
                        %tick_M(index) = 1;
                        %time_M(index) = cell_master{2,i}(2);
                        %
                        %if L1
                        if (cell_master{2,i}(1) == 0)

                            ph_M(:,index) = cell_master{3,i}(:,7);

                            %manage "nearly null" data
                            pos = abs(ph_M(:,index)) < 1e-100;
                            ph_M(pos,index) = 0;
                        end

                        type = [type '18 '];
                    end

                %message 19 (RTCM2)
                case 19

                    %buffer index computation
                    index = time_GPS - round(cell_master{2,i}(2)) + 1;

                    if (index <= B)

                        %buffer writing
                        tick_M(index) = 1;
                        time_M(index) = round(cell_master{2,i}(2));

                        %if L1
                        if (cell_master{2,i}(1) == 0)
                            pr_M(:,index) = cell_master{3,i}(:,7);
                        end

                        type = [type '19 '];
                    end

                %message 3 (RTCM2)
                case 3

                    coordX_M = cell_master{2,i}(1);
                    coordY_M = cell_master{2,i}(2);
                    coordZ_M = cell_master{2,i}(3);

                    if (flag_ms_pos & master_update)

                        if(index ~= 0)
                            pos_M(:,index) = [coordX_M; coordY_M; coordZ_M];
                            master_update = 0;
                            master_waiting = 0;
                        else
                            master_waiting = 1;
                        end
                    end

                    type = [type '3 '];

                %message 1002/1004 (RTCM3)
                case {1002, 1004}

                    %buffer index computation
                    index = time_GPS - round(cell_master{2,i}(2)) + 1;

                    if (index <= B)
                        while (index < 1)
                            time_GPS = time_GPS + 1;
                            index = time_GPS - round(cell_master{2,i}(2)) + 1;
                        end

                        %buffer writing
                        tick_M(index)  = 1;
                        time_M(index)  = cell_master{2,i}(2);
                        pr_M(:,index)  = cell_master{3,i}(:,2);
                        ph_M(:,index)  = cell_master{3,i}(:,3);
                        snr_M(:,index) = cell_master{3,i}(:,5);

                        %manage "nearly null" data
                        pos = abs(ph_M(:,index)) < 1e-100;
                        ph_M(pos,index) = 0;

                        type = [type num2str(cell_master{1,i}) ' '];

                    end

                %message 1005 (RTCM3)
                case 1005

                    coordX_M = cell_master{2,i}(8);
                    coordY_M = cell_master{2,i}(9);
                    coordZ_M = cell_master{2,i}(10);

                    if (flag_ms_pos & master_update)

                        if(index ~= 0)
                            pos_M(:,index) = [coordX_M; coordY_M; coordZ_M];
                            master_update = 0;
                            master_waiting = 0;
                        else
                            master_waiting = 1;
                        end
                    end

                    type = [type '1005 '];

                %message 1006 (RTCM3)
                case 1006

                    coordX_M = cell_master{2,i}(8);
                    coordY_M = cell_master{2,i}(9);
                    coordZ_M = cell_master{2,i}(10);
                    height_M = cell_master{2,i}(11); %#ok<NASGU>

                    if (flag_ms_pos & master_update)

                        if(index ~= 0)
                            pos_M(:,index) = [coordX_M; coordY_M; coordZ_M];
                            master_update = 0;
                            master_waiting = 0;
                        else
                            master_waiting = 1;
                        end
                    end

                    type = [type '1006 '];

                %message 1019 (RTCM3)
                case 1019

                    %satellite number
                    sat = cell_master{2,i}(1);

                    Eph(:,sat) = cell_master{2,i}(:);

                    type = [type '1019 '];

            end
            %if no master position is awaiting indexing
            if(~master_waiting)
                index = 0;
            end

            %if a master position is awaiting indexing
            if(index ~= 0 & master_waiting)
                pos_M(1,1:index) = coordX_M;
                pos_M(2,1:index) = coordY_M;
                pos_M(3,1:index) = coordZ_M;
                master_update = 0;
                master_waiting = 0;
            end
        end
    end

%     %Resolution of 2^23 cy carrier phase ambiguity
%     %caused by 32-bit data field restrictions (RTCM2)
%     if(test_master & is_rtcm2)
%         for i = 1 : length(index_ph)
%             pos = find(ph_M(:,index_ph(i)) & pr_M(:,index_ph(i)));
%             ambig = 2^23;
%             n = floor( (pr_M(pos,index_ph(i))/lambda1-ph_M(pos,index_ph(i))) / ambig + 0.5 );
%             ph_M(pos,index_ph(i)) = ph_M(pos,index_ph(i)) + n*ambig;
%         end
%     end

    %time reading (end of master decoding)
    current_time = toc;

    %----------------------------------

    %visualization
    i = min(b,B);                            %pointer to the last buffer cell
    if ~isempty(type)
        sat_pr = find(pr_M(:,i) ~= 0);       %satellites with code available
        sat_ph = find(ph_M(:,i) ~= 0);       %satellites with phase available
        sat = union(sat_pr,sat_ph);          %satellites with code or phase available

        %fprintf('MSG types: %s\n', type);
        fprintf('decoding: %7.4f sec (%smessages)\n', current_time-start_time, type);
        fprintf('GPStime=%d (%d satellites)\n', time_M(i), length(sat));

        fprintf('P1 SAT:');
        for p = 1 : length(sat_pr)
            fprintf(' %02d', sat_pr(p));
        end
        fprintf('\n');

        fprintf('L1 SAT:');
        r = 1;
        for p = 1 : length(sat_ph)
            while (sat_ph(p) ~= sat_pr(r))
                fprintf('   ');
                r = r + 1;
            end
            fprintf(' %02d', sat_ph(p));
            r = r + 1;
        end
        fprintf('\n');

    else
        fprintf('no messages\n');
    end

    fprintf('Station position:');
    if (sum(sum(abs(pos_M(:,i:end)))) ~= 0)
        fprintf(' X=%.4f, Y=%.4f, Z=%.4f km\n', pos_M(1,i)/1000, pos_M(2,i)/1000, pos_M(3,i)/1000);
    else
        fprintf(' not available\n');
    end

    %-------------------------------------
    % buffer situation
    %-------------------------------------

    fprintf('\n');
    fprintf('BUFFER (ROVER):  ');
    for i = B : -1 : 1
        if (tick_R(i) == 1)
            fprintf('x');
        else
            fprintf('o');
        end
    end
    fprintf('  --- time --->\n');
    fprintf('BUFFER (MASTER): ');
    for i = B : -1 : 1
        if (tick_M(i) == 1)
            fprintf('x');
        else
            fprintf('o');
        end
    end
    fprintf('  --- time --->\n');
    fprintf('                 ');
    for i = B : -1 : 1
        if (i == min(b,B))
            fprintf('^');
        else
            fprintf(' ');
        end
    end
    %fprintf('\n');

    %-------------------------------------
    % positioning
    %-------------------------------------

    %enable/disable positioning
    if 0==0

        %visualization
        fprintf('\n');
        fprintf('POSITIONING\n');

        %if the Kalman filter is not yet initialized
        if (t == 1)

            %current date reading for Google Earth visualization
            date = clock;

            %satellites with ephemerides available
            satEph = find(sum(abs(Eph))~=0);

            %delete data if ephemerides are not available
            %the buffer is activated only after the Kalman filter initialization
            delsat = setdiff(1:32,satEph);
            pr_R(delsat,1)  = 0;
            pr_M(delsat,1)  = 0;
            ph_R(delsat,1)  = 0;
            ph_M(delsat,1)  = 0;
            dop_R(delsat,1) = 0;
            snr_R(delsat,1) = 0;
            snr_M(delsat,1) = 0;

            %satellites with observations available
            %satObs = find( (pr_R(:,1) ~= 0) & (ph_R(:,1) ~= 0) & (pr_M(:,1) ~= 0) & (ph_M(:,1) ~= 0));
            satObs = find( (pr_R(:,1) ~= 0) & (pr_M(:,1) ~= 0));

            %if all the visible satellites ephemerides have been transmitted
            %and the total number of satellites is >= 4 and the master
            %station position is available
            if (ismember(satObs,satEph)) & (length(satObs) >= 4) & (sum(abs(pos_M(:,1))) ~= 0)
            %if (length(satObs_M) == length(satEph)) & (length(satObs) >= 4)

                %input data save
                fwrite(fid_obs, [time_GPS; time_M(1); time_R(1); week_R(1); pr_M(:,1); pr_R(:,1); ph_M(:,1); ph_R(:,1); dop_R(:,1); snr_M(:,1); snr_R(:,1); pos_M(:,1); iono(:,1)], 'double');
                fwrite(fid_eph, [time_GPS; Eph(:)], 'double');

                %WARNING: with just 4 satellites the least squares problem
                %for double differences is not solvable and the covariance matrix
                %of the estimation error cannot be computed

                %Kalman filter
                if (mode_vinc == 0)
                    kalman_goGPS_init (pos_M(:,1), time_M(1), Eph, iono, pr_R(:,1), pr_M(:,1), ph_R(:,1), ph_M(:,1), dop_R(:,1), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,1), snr_M(:,1), 1);
                else
                    kalman_goGPS_vinc_init (pos_M(:,1), time_M(1), Eph, iono, pr_R(:,1), pr_M(:,1), ph_R(:,1), ph_M(:,1), dop_R(:,1), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,1), snr_M(:,1), 1, ref_path);
                end
                
                %output data save
                if (mode_vinc == 0)
                    fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                    fwrite(fid_dop, [PDOP, HDOP, VDOP, KPDOP, KHDOP, KVDOP], 'double');
                else
                    fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
                    fwrite(fid_dop, [PDOP, HDOP, VDOP, 0, 0, 0], 'double');
                end
                fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');

                %estimated position and velocity
                if (mode_vinc == 0)
                    pos_t(1,1) = Xhat_t_t(1);
                    pos_t(2,1) = Xhat_t_t(o1+1);
                    pos_t(3,1) = Xhat_t_t(o2+1);
                    vel_t = sqrt(Xhat_t_t(2)^2+Xhat_t_t(o1+2)^2+Xhat_t_t(o2+2)^2);
                else
                    pos_t(1,1) = Yhat_t_t(1);
                    pos_t(2,1) = Yhat_t_t(2);
                    pos_t(3,1) = Yhat_t_t(3);
                    vel_t = Xhat_t_t(2);
                end

                %graphical representations
                if (flag_plotproc)
                    if (flag_cov == 0)
                        rtplot_matlab (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
                        if (flag_ge == 1), rtplot_googleearth (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), date), end
                    else
                        rtplot_matlab_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
                        if (flag_ge == 1), rtplot_googleearth_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date), end
                    end
                    if (flag_skyplot == 1)
                        rtplot_skyplot (t, azR, elR, conf_sat, pivot);
                        rtplot_snr (snr_R(:,1));
                    else
                        rttext_sat (t, azR, elR, snr_R(:,1), conf_sat, pivot);
                    end
                end

                %time reading
                current_time = toc;

                %visualization
                fprintf('kalman(%d): %7.4f sec (initialization)\n', t, current_time-start_time);
                fprintf('pos = (X=%.4f, Y=%.4f, Z=%.4f) km\n', pos_t(1)/1000, pos_t(2)/1000, pos_t(3)/1000);
                if (o1 ~= 1)
                    fprintf('vel = %.4f km/h\n', vel_t*3.6);
                end

                %send a new NMEA string
                if (flag_NTRIP)
                    nmea_update = sprintf('%s\r\n',NMEA_GGA_gen([pos_t(1); pos_t(2); pos_t(3)], sum(abs(conf_sat)), time_M(1), HDOP));
                    fwrite(master,nmea_update);
                    fprintf('NMEA sent\n');
                    master_update = 1;
                end

                %counter increment
                t = t + 1;

            else

                %visualization
                fprintf('no position/velocity are computed\n');

                %check Internet connection
                connected = 0;
                try
                    java.net.InetAddress.getByName(master_ip);
                catch
                    %close master connection
                    fclose(master);

                    %visualization
                    fprintf('wait for reconnection...\n');

                    %wait for connection
                    while ~connected
                        try java.net.InetAddress.getByName(master_ip)
                            connected = 1;
                        catch
                        end
                    end

                    %start a new connection
                    master = tcpip(master_ip,master_port);
                    set(master,'InputBufferSize', 16384);
                    fopen(master);
					if (flag_NTRIP) && (~isempty(nmea_sentences))
					    ntripstring = NTRIP_string_generator(nmea_init);
                        fwrite(master,ntripstring);
					end

                    %wait until the buffer writing is started before continuing
                    while get(master,'BytesAvailable') == 0, end;
                end
            end

            %buffer pointer to zero
            b = 0;

            %---------------------------------------------------------------------------

        %if the Kalman filter was already initialized
        else %if (t > 1)

            %signal loss because the buffer is too small
            while (b > B)

                %current date for visualization on Google Earth
                date = clock;

                %input data save
                fwrite(fid_obs, [time_GPS; 0; 0; 0; zeros(32,1); zeros(32,1); zeros(32,1); zeros(32,1); zeros(32,1); zeros(32,1); zeros(32,1); zeros(3,1); zeros(8,1)], 'double');
                fwrite(fid_eph, [time_GPS; Eph(:)], 'double');

                %Kalman filter
                if (mode_vinc == 0)
                    [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop (zeros(3,1), 0, Eph, iono, zeros(32,1), zeros(32,1), zeros(32,1), zeros(32,1), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, zeros(32,1), zeros(32,1), 1); %#ok<NASGU>
                else
                    [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop (zeros(3,1), 0, Eph, iono, zeros(32,1), zeros(32,1), zeros(32,1), zeros(32,1), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, zeros(32,1), zeros(32,1), 1, ref_path); %#ok<NASGU>
                end

                %output data save
                if (mode_vinc == 0)
                    fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                    fwrite(fid_dop, [PDOP, HDOP, VDOP, KPDOP, KHDOP, KVDOP], 'double');
                else
                    fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
                    fwrite(fid_dop, [PDOP, HDOP, VDOP, 0, 0, 0], 'double');
                end
                fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');

                %estimated position and velocity
                if (mode_vinc == 0)
                    pos_t(1,1) = Xhat_t_t(1);
                    pos_t(2,1) = Xhat_t_t(o1+1);
                    pos_t(3,1) = Xhat_t_t(o2+1);
                    vel_t = sqrt(Xhat_t_t(2)^2+Xhat_t_t(o1+2)^2+Xhat_t_t(o2+2)^2);
                else
                    pos_t(1,1) = Yhat_t_t(1);
                    pos_t(2,1) = Yhat_t_t(2);
                    pos_t(3,1) = Yhat_t_t(3);
                    vel_t = Xhat_t_t(2);
                end

                %graphical representations
                if (flag_plotproc)
                    if (flag_cov == 0)
                        rtplot_matlab (t, [pos_t(1); pos_t(2); pos_t(3)], zeros(3,1), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
                        if (flag_ge == 1), rtplot_googleearth (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), date), end
                    else
                        rtplot_matlab_cov (t, [pos_t(1); pos_t(2); pos_t(3)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
                        if (flag_ge == 1), rtplot_googleearth_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date), end
                    end
                    if (flag_skyplot == 1)
                        rtplot_skyplot (t, azR, elR, conf_sat, pivot);
                        rtplot_snr (zeros(32,1));
                    else
                        rttext_sat (t, azR, elR, zeros(32,1), conf_sat, pivot);
                    end
                end

                %time reading
                current_time = toc;

                %visualization
                fprintf('kalman(%d): %7.4f sec (over-buffering)\n', t, current_time-start_time);
                fprintf('pos = (X=%.4f, Y=%.4f, Z=%.4f) km\n', pos_t(1)/1000, pos_t(2)/1000, pos_t(3)/1000);
                if (o1 ~= 1)
                    fprintf('vel = %.4f km/h\n', vel_t*3.6);
                end

                %send a new NMEA string
                if ((flag_NTRIP) && (mod(t,nmea_update_rate) == 0))
                    nmea_update = sprintf('%s\r\n',NMEA_GGA_gen([pos_t(1); pos_t(2); pos_t(3)], sum(abs(conf_sat))));
                    fwrite(master,nmea_update);
                    fprintf('NMEA sent\n');
                    master_update = 1;
                end
                
                %poll a new AID-HUI message
                if (mod(t,hui_poll_rate) == 0)
                    if (protocol == 0)
                        ublox_poll_message(rover, 'AID', 'HUI', 0);
                        fprintf('AID-HUI msg polled\n');
                    end
                end

                %counter increment
                t = t + 1;

                %buffer pointer decrement
                b = b - 1;

            end %after this point the pointer is inside the buffer

            %-----------------------------------------------------------------------

            %loss of master signal (data is not available "now" but it is available "in the future")
            if (tick_M(b) == 0) & (sum(tick_M(1:b)) > 0)

                %data loss management
                while (tick_M(b) == 0)

                    %current date for visualization on Google Earth
                    date = clock;

                    %satellites with available ephemerides
                    satEph = find(sum(abs(Eph))~=0);

                    %delete data if ephemerides are not available
                    delsat = setdiff(1:32,satEph);
                    pr_R(delsat,b)  = 0;
                    pr_M(delsat,b)  = 0;
                    ph_R(delsat,b)  = 0;
                    ph_M(delsat,b)  = 0;
                    dop_R(delsat,b) = 0;
                    snr_R(delsat,b) = 0;
                    snr_M(delsat,b) = 0;

                    %satellites with available observations
                    satObs = find( (pr_R(:,b) ~= 0) & (pr_M(:,b) ~= 0));

                    %input data save
                    fwrite(fid_obs, [time_GPS; time_M(b); time_R(b); week_R(b); pr_M(:,b); pr_R(:,b); ph_M(:,b); ph_R(:,b); dop_R(:,b); snr_M(:,b); snr_R(:,b); pos_M(:,b); iono(:,1)], 'double');
                    fwrite(fid_eph, [time_GPS; Eph(:)], 'double');

                    %Kalman filter
                    if (mode_vinc == 0)
                        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop (pos_M(:,b), time_M(b), Eph, iono, pr_R(:,b), pr_M(:,b), ph_R(:,b), ph_M(:,b), dop_R(:,b), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,b), snr_M(:,b), 1);
                    else
                        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop (pos_M(:,b), time_M(b), Eph, iono, pr_R(:,b), pr_M(:,b), ph_R(:,b), ph_M(:,b), dop_R(:,b), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,b), snr_M(:,b), 1, ref_path);
                    end

                    %output data save
                    if (mode_vinc == 0)
                        fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                        fwrite(fid_dop, [PDOP, HDOP, VDOP, KPDOP, KHDOP, KVDOP], 'double');
                    else
                        fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
                        fwrite(fid_dop, [PDOP, HDOP, VDOP, 0, 0, 0], 'double');
                    end
                    fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                    fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');

                    %estimated position and velocity
                    if (mode_vinc == 0)
                        pos_t(1,1) = Xhat_t_t(1);
                        pos_t(2,1) = Xhat_t_t(o1+1);
                        pos_t(3,1) = Xhat_t_t(o2+1);
                        vel_t = sqrt(Xhat_t_t(2)^2+Xhat_t_t(o1+2)^2+Xhat_t_t(o2+2)^2);
                    else
                        pos_t(1,1) = Yhat_t_t(1);
                        pos_t(2,1) = Yhat_t_t(2);
                        pos_t(3,1) = Yhat_t_t(3);
                        vel_t = Xhat_t_t(2);
                    end

                    %graphical representations
                    if (flag_plotproc)
                        if (flag_cov == 0)
                            rtplot_matlab (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,b), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
                            if (flag_ge == 1), rtplot_googleearth (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), date), end
                        else
                            rtplot_matlab_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,b), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
                            if (flag_ge == 1), rtplot_googleearth_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date), end
                        end
                        if (flag_skyplot == 1)
                            rtplot_skyplot (t, azR, elR, conf_sat, pivot);
                            rtplot_snr (snr_R(:,b));
                        else
                            rttext_sat (t, azR, elR, snr_R(:,b), conf_sat, pivot);
                        end
                    end

                    %time reading
                    current_time = toc;

                    %visualization
                    fprintf('kalman(%d): %7.4f sec (master loss)\n', t, current_time-start_time);
                    fprintf('pos = (X=%.4f, Y=%.4f, Z=%.4f) km\n', pos_t(1)/1000, pos_t(2)/1000, pos_t(3)/1000);
                    if (o1 ~= 1)
                        fprintf('vel = %.4f km/h\n', vel_t*3.6);
                    end

                    %send a new NMEA string
                    if ((flag_NTRIP) && (mod(t,nmea_update_rate) == 0))
                        nmea_update = sprintf('%s\r\n',NMEA_GGA_gen([pos_t(1); pos_t(2); pos_t(3)], sum(abs(conf_sat)), time_R(b), HDOP));
                        fwrite(master,nmea_update);
                        fprintf('NMEA sent\n');
                        master_update = 1;
                    end
                    
                    %poll a new AID-HUI message
                    if (mod(t,hui_poll_rate) == 0)
                        if (protocol == 0)
                            ublox_poll_message(rover, 'AID', 'HUI', 0);
                        end
                    end

                    %counter increment
                    t = t + 1;

                    %buffer pointer decrement
                    b = b - 1;
                end
            end %after this point there are no more master losses

            %-----------------------------------------------------------------------

            %loss of rover signal (data is not available "now" but it is available "in the future")
            if (tick_R(b) == 0) & (sum(tick_R(1:b)) > 0)

                %data loss management
                while (tick_R(b) == 0)

                    %current date for visualization on Google Earth
                    date = clock;

                    %satellites with available ephemerides
                    satEph = find(sum(abs(Eph))~=0);

                    %delete data if ephemerides are not available
                    delsat = setdiff(1:32,satEph);
                    pr_R(delsat,b)  = 0;
                    pr_M(delsat,b)  = 0;
                    ph_R(delsat,b)  = 0;
                    ph_M(delsat,b)  = 0;
                    dop_R(delsat,b) = 0;
                    snr_R(delsat,b) = 0;
                    snr_M(delsat,b) = 0;

                    %satellites with available observations
                    satObs = find( (pr_R(:,b) ~= 0) & (pr_M(:,b) ~= 0));

                    %input data save
                    fwrite(fid_obs, [time_GPS; time_M(b); time_R(b); week_R(b); pr_M(:,b); pr_R(:,b); ph_M(:,b); ph_R(:,b); dop_R(:,b); snr_M(:,b); snr_R(:,b); pos_M(:,b); iono(:,1)], 'double');
                    fwrite(fid_eph, [time_GPS; Eph(:)], 'double');

                    %Kalman filter
                    if (mode_vinc == 0)
                        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop (pos_M(:,b), time_M(b), Eph, iono, pr_R(:,b), pr_M(:,b), ph_R(:,b), ph_M(:,b), dop_R(:,b), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,b), snr_M(:,b), 1);
                    else
                        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop (pos_M(:,b), time_M(b), Eph, iono, pr_R(:,b), pr_M(:,b), ph_R(:,b), ph_M(:,b), dop_R(:,b), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,b), snr_M(:,b), 1, ref_path);
                    end

                    %output data save
                    if (mode_vinc == 0)
                        fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                        fwrite(fid_dop, [PDOP, HDOP, VDOP, KPDOP, KHDOP, KVDOP], 'double');
                    else
                        fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
                        fwrite(fid_dop, [PDOP, HDOP, VDOP, 0, 0, 0], 'double');
                    end
                    fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                    fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');

                    %estimated position and velocity
                    if (mode_vinc == 0)
                        pos_t(1,1) = Xhat_t_t(1);
                        pos_t(2,1) = Xhat_t_t(o1+1);
                        pos_t(3,1) = Xhat_t_t(o2+1);
                        vel_t = sqrt(Xhat_t_t(2)^2+Xhat_t_t(o1+2)^2+Xhat_t_t(o2+2)^2);
                    else
                        pos_t(1,1) = Yhat_t_t(1);
                        pos_t(2,1) = Yhat_t_t(2);
                        pos_t(3,1) = Yhat_t_t(3);
                        vel_t = Xhat_t_t(2);
                    end

                    %graphical representations
                    if (flag_plotproc)
                        if (flag_cov == 0)
                            rtplot_matlab (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,b),check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
                            if (flag_ge == 1), rtplot_googleearth (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), date), end
                        else
                            rtplot_matlab_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,b), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
                            if (flag_ge == 1), rtplot_googleearth_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date), end
                        end
                        if (flag_skyplot == 1)
                            rtplot_skyplot (t, azR, elR, conf_sat, pivot);
                            rtplot_snr (snr_R(:,b));
                        else
                            rttext_sat (t, azR, elR, snr_R(:,b), conf_sat, pivot);
                        end
                    end

                    %time reading
                    current_time = toc;

                    %visualization
                    fprintf('kalman(%d): %7.4f sec (rover loss)\n', t, current_time-start_time);
                    fprintf('pos = (X=%.4f, Y=%.4f, Z=%.4f) km\n', pos_t(1)/1000, pos_t(2)/1000, pos_t(3)/1000);
                    if (o1 ~= 1)
                        fprintf('vel = %.4f km/h\n', vel_t*3.6);
                    end

                    %send a new NMEA string
                    if ((flag_NTRIP) && (mod(t,nmea_update_rate) == 0))
                        nmea_update = sprintf('%s\r\n',NMEA_GGA_gen([pos_t(1); pos_t(2); pos_t(3)], sum(abs(conf_sat)), time_R(b), HDOP));
                        fwrite(master,nmea_update);
                        fprintf('NMEA sent\n');
                        master_update = 1;
                    end
                    
                    %poll a new AID-HUI message
                    if (mod(t,hui_poll_rate) == 0)
                        if (protocol == 0)
                            ublox_poll_message(rover, 'AID', 'HUI', 0);
                        end
                    end

                    %counter increment
                    t = t + 1;

                    %buffer pointer decrement
                    b = b - 1;
                end
            end %after this point there are no more master nor rover losses

            %-----------------------------------------------------------------------

            %delay in master or rover signal
            if (tick_M(b) == 0) | (tick_R(b) == 0)

                %safety threshold for the buffer (probably no more useful)
                %safety_B = 2;
                %safety_B = min(B,safety_B);

                %if there is still space in the buffer, wait!
                if (b < B)
                %if (b < B - safety_B)
                    
                    %visualization
                    fprintf('wait for data (delay=%d sec)\n',b);
                    
                %otherwise make one or more step using just dynamics
                else
                    
                    %if the master is interrupted
                    if (tick_M(b) == 0)
                        
                        %check Internet connection
                        try
                            java.net.InetAddress.getByName(master_ip);

                            %clear the whole buffer
                            lastB = 1;
                            
                            %close master connection
                            fclose(master);
                            
                            %start a new connection
                            master = tcpip(master_ip,master_port);
                            set(master,'InputBufferSize', 16384);
                            fopen(master);
                            if ((flag_NTRIP) && (~isempty(nmea_sentences)))
                                nmea_update = sprintf('%s\r\n',NMEA_GGA_gen([pos_t(1); pos_t(2); pos_t(3)], sum(abs(conf_sat))));
                                ntripstring = NTRIP_string_generator(nmea_update);
                                fwrite(master,ntripstring);
                            end
                        catch
                            %clear just the last buffer cell
                            lastB = B;
                        end
                        
                    %if the rover is interrupted
                    else
                        
                        %clear just the last buffer cell
                        lastB = B;
                        %lastB = B - safety_B;
                        
                    end

                    %clear the buffer up to the desired position
                    while (b >= lastB)

                        %current date for visualization on Google Earth
                        date = clock;

                        %satellites for which there are avilable ephemerides
                        satEph = find(sum(abs(Eph))~=0);

                        %delete data if ephemerides are not available (b=B)
                        delsat = setdiff(1:32,satEph);
                        pr_R(delsat,b)  = 0;
                        pr_M(delsat,b)  = 0;
                        ph_R(delsat,b)  = 0;
                        ph_M(delsat,b)  = 0;
                        dop_R(delsat,b) = 0;
                        snr_R(delsat,b) = 0;
                        snr_M(delsat,b) = 0;

                        %satellites with available observations
                        satObs = find( (pr_R(:,b) ~= 0) & (pr_M(:,b) ~= 0));

                        %output data save
                        fwrite(fid_obs, [time_GPS; time_M(b); time_R(b); week_R(b); pr_M(:,b); pr_R(:,b); ph_M(:,b); ph_R(:,b); dop_R(:,b); snr_M(:,b); snr_R(:,b); pos_M(:,b); iono(:,1)], 'double');
                        fwrite(fid_eph, [time_GPS; Eph(:)], 'double');

                        %Kalman filter
                        if (mode_vinc == 0)
                            [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop (pos_M(:,b), time_M(b), Eph, iono, pr_R(:,b), pr_M(:,b), ph_R(:,b), ph_M(:,b), dop_R(:,b), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,b), snr_M(:,b), 1);
                        else
                            [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop (pos_M(:,b), time_M(b), Eph, iono, pr_R(:,b), pr_M(:,b), ph_R(:,b), ph_M(:,b), dop_R(:,b), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,b), snr_M(:,b), 1, ref_path);
                        end

                        %output data save
                        if (mode_vinc == 0)
                            fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                            fwrite(fid_dop, [PDOP, HDOP, VDOP, KPDOP, KHDOP, KVDOP], 'double');
                        else
                            fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
                            fwrite(fid_dop, [PDOP, HDOP, VDOP, 0, 0, 0], 'double');
                        end
                        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                        fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');

                        %estimated position and velocity
                        if (mode_vinc == 0)
                            pos_t(1,1) = Xhat_t_t(1);
                            pos_t(2,1) = Xhat_t_t(o1+1);
                            pos_t(3,1) = Xhat_t_t(o2+1);
                            vel_t = sqrt(Xhat_t_t(2)^2+Xhat_t_t(o1+2)^2+Xhat_t_t(o2+2)^2);
                        else
                            pos_t(1,1) = Yhat_t_t(1);
                            pos_t(2,1) = Yhat_t_t(2);
                            pos_t(3,1) = Yhat_t_t(3);
                            vel_t = Xhat_t_t(2);
                        end

                        %graphical representations
                        if (flag_plotproc)
                            if (flag_cov == 0)
                                rtplot_matlab (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,b), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
                                if (flag_ge == 1), rtplot_googleearth (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), date), end
                            else
                                rtplot_matlab_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,b), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
                                if (flag_ge == 1), rtplot_googleearth_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date), end
                            end
                            if (flag_skyplot == 1)
                                rtplot_skyplot (t, azR, elR, conf_sat, pivot);
                                rtplot_snr (snr_R(:,b));
                            else
                                rttext_sat (t, azR, elR, snr_R(:,b), conf_sat, pivot);
                            end
                        end

                        %time reading
                        current_time = toc;

                        %visualization
                        fprintf('kalman(%d): %7.4f sec (dyn-only, delay=%d sec)\n', t, current_time-start_time, b-1);
                        fprintf('pos = (X=%.4f, Y=%.4f, Z=%.4f) km\n', pos_t(1)/1000, pos_t(2)/1000, pos_t(3)/1000);
                        if (o1 ~= 1)
                            fprintf('vel = %.4f km/h\n', vel_t*3.6);
                        end

                        %send a new NMEA string
                        if ((flag_NTRIP) && (mod(t,nmea_update_rate) == 0))
                            nmea_update = sprintf('%s\r\n',NMEA_GGA_gen([pos_t(1); pos_t(2); pos_t(3)], sum(abs(conf_sat)), time_R(b), HDOP));
                            fwrite(master,nmea_update);
                            fprintf('NMEA sent\n');
                            master_update = 1;
                        end
                        
                        %poll a new AID-HUI message
                        if (mod(t,hui_poll_rate) == 0)
                            if (protocol == 0)
                                ublox_poll_message(rover, 'AID', 'HUI', 0);
                            end
                        end

                        %counter increment
                        t = t + 1;

                        %buffer pointer decrement
                        b = b - 1;
                    end
                end %after this point there are no more losses nor delays

                %-----------------------------------------------------------------------

            %data available for both the master and the rover
            else

                %Kalman filter
                while (b > 0) & (tick_M(b) == tick_R(b)) & (tick_M(b) == 1)

                    %current date for visualization on Google Earth
                    date = clock;

                    %satellites for which there are available ephemerides
                    satEph = find(sum(abs(Eph))~=0);

                    %delete data if ephemerides are not available
                    delsat = setdiff(1:32,satEph);
                    pr_R(delsat,b)  = 0;
                    pr_M(delsat,b)  = 0;
                    ph_R(delsat,b)  = 0;
                    ph_M(delsat,b)  = 0;
                    dop_R(delsat,b) = 0;
                    snr_R(delsat,b) = 0;
                    snr_M(delsat,b) = 0;

                    %satellites with available observations
                    satObs = find( (pr_R(:,b) ~= 0) & (pr_M(:,b) ~= 0));

                    %input data save
                    fwrite(fid_obs, [time_GPS; time_M(b); time_R(b); week_R(b); pr_M(:,b); pr_R(:,b); ph_M(:,b); ph_R(:,b); dop_R(:,b); snr_M(:,b); snr_R(:,b); pos_M(:,b); iono(:,1)], 'double');
                    fwrite(fid_eph, [time_GPS; Eph(:)], 'double');

                    %Kalman filter
                    if (mode_vinc == 0)
                        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop (pos_M(:,b), time_M(b), Eph, iono, pr_R(:,b), pr_M(:,b), ph_R(:,b), ph_M(:,b), dop_R(:,b), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,b), snr_M(:,b), 1);
                    else
                        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop (pos_M(:,b), time_M(b), Eph, iono, pr_R(:,b), pr_M(:,b), ph_R(:,b), ph_M(:,b), dop_R(:,b), dop1_M, pr2_R, pr2_M, ph2_R, ph2_M, dop2_R, dop2_M, snr_R(:,b), snr_M(:,b), 1, ref_path);
                    end

                    %output data save
                    if (mode_vinc == 0)
                        fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                        fwrite(fid_dop, [PDOP, HDOP, VDOP, KPDOP, KHDOP, KVDOP], 'double');
                    else
                        fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
                        fwrite(fid_dop, [PDOP, HDOP, VDOP, 0, 0, 0], 'double');
                    end
                    fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                    fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');

                    %estimated position and velocity
                    if (mode_vinc == 0)
                        pos_t(1,1) = Xhat_t_t(1);
                        pos_t(2,1) = Xhat_t_t(o1+1);
                        pos_t(3,1) = Xhat_t_t(o2+1);
                        vel_t = sqrt(Xhat_t_t(2)^2+Xhat_t_t(o1+2)^2+Xhat_t_t(o2+2)^2);
                    else
                        pos_t(1,1) = Yhat_t_t(1);
                        pos_t(2,1) = Yhat_t_t(2);
                        pos_t(3,1) = Yhat_t_t(3);
                        vel_t = Xhat_t_t(2);
                    end

                    %graphical representations
                    if (flag_plotproc)
                        if (flag_cov == 0)
                            rtplot_matlab (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,b), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
                            if (flag_ge == 1), rtplot_googleearth (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), date), end
                        else
                            rtplot_matlab_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,b), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
                            if (flag_ge == 1), rtplot_googleearth_cov (t, [pos_t(1); pos_t(2); pos_t(3)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date), end
                        end
                        if (flag_skyplot == 1)
                            rtplot_skyplot (t, azR, elR, conf_sat, pivot);
                            rtplot_snr (snr_R(:,b));
                        else
                            rttext_sat (t, azR, elR, snr_R(:,b), conf_sat, pivot);
                        end
                    end

                    %time reading
                    current_time = toc;

                    %visualization
                    if (b > 1)
                        fprintf('kalman(%d): %7.4f sec (delay=%d sec)\n', t, current_time-start_time, b-1);
                    else
                        fprintf('kalman(%d): %7.4f sec\n', t, current_time-start_time);
                    end
                    fprintf('pos = (X=%.4f, Y=%.4f, Z=%.4f) km\n', pos_t(1)/1000, pos_t(2)/1000, pos_t(3)/1000);
                    if (o1 ~= 1)
                        fprintf('vel = %.4f km/h\n', vel_t*3.6);
                    end
                    fprintf('pivot = %02d\n', pivot)
                    if (check_cs == 1)
                        fprintf('cycle-slip =');
                        fprintf(' %02d', find(conf_cs==1));
                        fprintf('\n');
                    end

                    %send a new NMEA string
                    if ((flag_NTRIP) && (mod(t,nmea_update_rate) == 0))
                        nmea_update = sprintf('%s\r\n',NMEA_GGA_gen([pos_t(1); pos_t(2); pos_t(3)], sum(abs(conf_sat)), time_R(b), HDOP));
                        fwrite(master,nmea_update);
                        fprintf('NMEA sent\n');
                        master_update = 1;
                    end
                    
                    %poll a new AID-HUI message
                    if (mod(t,hui_poll_rate) == 0)
                        if (protocol == 0)
                            ublox_poll_message(rover, 'AID', 'HUI', 0);
                        end
                    end

                    %counter increment
                    t = t + 1;

                    %buffer pointer decrement
                    b = b - 1;
                end
            end %end of data processing
        end

        %pos_KAL = [Xhat_t_t(1,:); Xhat_t_t(o1+1,:); Xhat_t_t(o2+1,:)];

        %-------------------------------------

        %visualization
        fprintf('EPH SAT:');
        for i = 1 : length(satEph)
            fprintf(' %02d', satEph(i));
        end
        fprintf('\n');

        fprintf('OBS SAT:');
        j = 1;
        for i = 1 : length(satObs)
            while (satObs(i) ~= satEph(j))
                fprintf('   ');
                j = j + 1;
            end
            fprintf(' %02d', satObs(i));
            j = j + 1;
        end
        fprintf('\n');

        %-------------------------------------

    %disable positioning (debug)
    end

    %----------------------------------

    %test if the cycle execution has ended
    flag = getappdata(gcf, 'run');
    drawnow

    %-------------------------------------

    %computation of the delays due to data processing (by default dtime=1)
    dtime1 = ceil(current_time-start_time-step_time);

    %computation of the delays due to external causes after data processing
    current_time = toc;
    dtime2 = ceil(current_time-start_time-step_time);

    if (dtime2 > dtime1)
        fprintf('WARNING! System slowdown: %7.4f sec (delay=%d sec)\n',current_time-start_time, dtime2-dtime1);
    end
    dtime = dtime2;

    %go to next epoch
    while (current_time-start_time-step_time < dtime)
        current_time = toc;
    end

    %GPS epoch increment
    %time_GPS = time_GPS + 1;
    time_GPS = time_GPS + dtime;

    %buffer pointer increment
    b = b + dtime;

    %starting time re-initialization
    if (t == 1)
        start_time = start_time + dtime;
    end

    %clear screen
    clc

end

%------------------------------------------------------
% tasks at the end of the cycle
%------------------------------------------------------

if (protocol == 0)
    %load u-blox saved configuration
    if (reply_save)
        fprintf('Restoring saved u-blox receiver configuration...\n');

        reply_load = ublox_CFG_CFG(rover, 'load');
        tries = 0;

        while (reply_save & ~reply_load)
            tries = tries + 1;
            if (tries > 3)
                disp('It was not possible to reload the receiver previous configuration.');
                break
            end
            reply_load = ublox_CFG_CFG(rover, 'load');
        end
    end
end

%connection closing
fclose(master);
fclose(rover);
delete(master);
delete(rover);

%data files closing
fclose(fid_master);
fclose(fid_rover);
fclose(fid_obs);
fclose(fid_eph);
fclose(fid_kal);
fclose(fid_sat);
fclose(fid_dop);
fclose(fid_conf);
fclose(fid_nmea);

%log file closing
diary off

if (flag_plotproc)
    %hide STOP button
    set(h1,'visible','off');
else
    %close figure
    close(f1);
end
