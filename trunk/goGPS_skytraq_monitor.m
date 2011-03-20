function goGPS_skytraq_monitor(filerootOUT)

% SYNTAX:
%   goGPS_skytraq_monitor(filerootOUT)
%
% INPUT:
%   filerootOUT = output file prefix
%
% DESCRIPTION:
%   SkyTraq receiver monitor: stream reading, output data saving.

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
rover = serial (COMportR,'BaudRate',115200);
set(rover,'InputBufferSize',16384);
set(rover,'FlowControl','hardware');
set(rover,'RequestToSend','on');
fopen(rover);

%------------------------------------------------------
% absolute time startup
%------------------------------------------------------

tic

%initialize log file
delete('../data/STQ_monitor.txt');
diary('../data/STQ_monitor.txt');
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
    fprintf('SkyTraq: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time, rover_1, rover_2);

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
    fprintf('SkyTraq: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time, rover_1, rover_2);

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

%poll flags
eph_polled = 0;
hui_polled = 0;

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

        [cell_rover] = decode_skytraq(data_rover);

        %read data type
        type = '';

        %data type counters
        nRAW = 0;
        nTIM = 0;
        
        for i = 1 : size(cell_rover,2)
            
            %MEAS_TIME message data save
            if (strcmp(cell_rover{1,i},'MEAS_TIME'))
                
                time_R = cell_rover{2,i}(3);
                week_R = cell_rover{2,i}(2);
                
                type = [type 'MEAS_TIME '];
                nTIM = nTIM + 1;
                
            %RAW_MEAS message data save
            elseif (strcmp(cell_rover{1,i},'RAW_MEAS'))
                
                pr_R = cell_rover{3,i}(:,3);
                ph_R = cell_rover{3,i}(:,4);
                snr_R = cell_rover{3,i}(:,2);
                dop_R = cell_rover{3,i}(:,5);
                
                %manage "nearly null" data
                ph_R(abs(ph_R) < 1e-100) = 0;
                
                type = [type 'RAW_MEAS '];
                nRAW = nRAW + 1;
                
                %counter increment
                t = t+1;
            end
            
            if (nRAW > 0 & nTIM > 0)
                
                %data save
                fwrite(fid_obs, [0; 0; time_R; week_R; zeros(32,1); pr_R; zeros(32,1); ph_R; dop_R; zeros(32,1); snr_R; zeros(3,1); zeros(8,1)], 'double');
            end
        end
        
        %----------------------------------

        %visualization
        fprintf('\n');
        fprintf('SkyTraq: %7.4f sec (%4d bytes --> %4d bytes)\n', current_time-start_time, rover_1, rover_2);
        fprintf('MSG types: %s\n', type);

        %visualization (RXM-RAW information)
        if (nRAW > 0 & nTIM > 0)
            sat_pr = find(pr_R ~= 0);       %satellites with code available
            sat_ph = find(ph_R ~= 0);       %satellites with phase available
            sat = union(sat_pr,sat_ph);     %satellites with code or phase available

            if (i < length(time_R)), fprintf(' DELAYED\n'); else fprintf('\n'); end
            fprintf('Epoch %3d:  GPStime=%d:%.3f (%d satellites)\n', t, week_R, time_R, length(sat));
            for j = 1 : length(sat)
                fprintf('   SAT %02d:  P1=%11.2f  L1=%12.2f  D1=%7.1f  QI=%1d  SNR=%2d  LOCK=%1d\n', ...
                    sat(j), pr_R(sat(j)), ph_R(sat(j)), dop_R(sat(j)), 0, snr_R(sat(j)), 0);
            end
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