function [time_GPS, week_R, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, ...
          dop1_R, snr_R, snr_M, pos_M, Eph, iono, delay, loss_R, loss_M] = load_goGPSinput (fileroot)

% SYNTAX:
%   [time_GPS, week_R, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, ...
%    dop1_R, snr_R, snr_M, pos_M, Eph, iono, delay, loss_R, loss_M] = load_goGPSinput (fileroot);
%
% INPUT:
%   fileroot = name of the file to be read
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
%   dop1_R    = ROVER-SATELLITE Doppler observations (carrier L1)
%   snr_R    = ROVER signal-to-noise ratio
%   snr_M    = MASTER signal-to-fnoise ratio
%   pos_M    = MASTER station coordinates
%   Eph      = matrix of 33 parameters each satellite (MASTER)
%   iono     = ionosphere parameters
%   delay    = delay in observations processing
%   loss_R   = flag for the ROVER loss of signal
%   loss_M   = flag for the MASTER loss of signal
%
% DESCRIPTION:
%   Kalman filter input data reading.

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

%initialization
time_GPS = [];                         %reference GPS time
time_M = [];                           %MASTER GPS time
time_R = [];                           %ROVER GPS time
week_R = [];                           %GPS week number
pr1_M = [];                            %MASTER code observations
pr1_R = [];                            %ROVER code observations
ph1_M = [];                            %MASTER phase observations
ph1_R = [];                            %ROVER phase observations
dop1_R = [];                           %ROVER Doppler observations
snr_M = [];                            %MASTER signal-to-noise ratio
snr_R = [];                            %ROVER signal-to-noise ratio
pos_M = [];                            %MASTER station coordinates
iono = [];                             %ionosphere parameters

%default number of satellites (only GPS)
num_sat = 32;

%observations reading
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%03d');                                    %hour index (string)
d = dir([fileroot '_obs_' hour_str '.bin']);                        %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_obs_' hour_str '.bin\n']);
    fid_obs = fopen([fileroot '_obs_' hour_str '.bin']);            %file opening
    num_sat = fread(fid_obs,1,'int8');                              %read number of satellites
    num_bytes = d.bytes-1;                                          %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / (4+num_sat*7+11);                       %file size (number of packets)
    buf_obs = fread(fid_obs,num_words,'double');                    %file reading
    fclose(fid_obs);                                                %file closing
    time_GPS = [time_GPS; zeros(num_packs,1)];                      %observations concatenation
    time_M   = [time_M;   zeros(num_packs,1)];
    time_R   = [time_R;   zeros(num_packs,1)];
    week_R   = [week_R;   zeros(num_packs,1)];
    pr1_M    = [pr1_M     zeros(num_sat,num_packs)];
    pr1_R    = [pr1_R     zeros(num_sat,num_packs)];
    ph1_M    = [ph1_M     zeros(num_sat,num_packs)];
    ph1_R    = [ph1_R     zeros(num_sat,num_packs)];
    dop1_R   = [dop1_R    zeros(num_sat,num_packs)];
    snr_M    = [snr_M     zeros(num_sat,num_packs)];
    snr_R    = [snr_R     zeros(num_sat,num_packs)];
    pos_M    = [pos_M     zeros(3,num_packs)];
    iono     = [iono      zeros(8,num_packs)];
    for j = 0 : (4+num_sat*7+11) : num_words-1
        i = i+1;                                                    %epoch counter increase
        time_GPS(i,1) = buf_obs(j + 1);                             %observations logging
        time_M(i,1)   = buf_obs(j + 2);
        time_R(i,1)   = buf_obs(j + 3);
        week_R(i,1)   = buf_obs(j + 4);
        pr1_M(:,i)    = buf_obs(j + [4+1:4+1*num_sat]);
        pr1_R(:,i)    = buf_obs(j + [4+1*num_sat+1:4+2*num_sat]);
        ph1_M(:,i)    = buf_obs(j + [4+2*num_sat+1:4+3*num_sat]);
        ph1_R(:,i)    = buf_obs(j + [4+3*num_sat+1:4+4*num_sat]);
        dop1_R(:,i)   = buf_obs(j + [4+4*num_sat+1:4+5*num_sat]);
        snr_M(:,i)    = buf_obs(j + [4+5*num_sat+1:4+6*num_sat]);
        snr_R(:,i)    = buf_obs(j + [4+6*num_sat+1:4+7*num_sat]);
        pos_M(:,i)    = buf_obs(j + [4+7*num_sat+1:4+7*num_sat+3]);
        iono(:,i)     = buf_obs(j + [4+7*num_sat+4:4+7*num_sat+11]);
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%03d');                                %conversion into a string
    d = dir([fileroot '_obs_' hour_str '.bin']);                    %file to be read
end

%-------------------------------------------------------------------------------

%initialization
Eph = [];
num_eph = 33;

%read ephemerides
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%03d');                                    %hour index (string)
d = dir([fileroot '_eph_' hour_str '.bin']);                        %file to be read
if isempty(d)
    Eph = zeros(num_eph,num_sat,length(time_GPS(:,1)));
end
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_eph_' hour_str '.bin\n']);
    fid_eph = fopen([fileroot '_eph_' hour_str '.bin']);            %file opening
    num_sat = fread(fid_eph,1,'int8');                              %read number of satellites
    num_bytes = d.bytes-1;                                          %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / (1+num_eph*num_sat);                    %file size (number of packets)
    buf_eph = fread(fid_eph,num_words,'double');                    %file reading
    fclose(fid_eph);                                                %file closing
    Eph = cat(3,Eph,zeros(num_eph,num_sat,num_packs));              %ephemerides concatenation
    for j = 0 : (1+num_eph*num_sat) : num_words-1
        i = i+1;                                                    %epoch counter increase
        %time_GPS(i,1) = buf_eph(j + 1);                            %GPS time logging
        Eph(:,:,i) = reshape(buf_eph(j + [2:1+num_eph*num_sat]), ...
                             [num_eph,num_sat]);                    %ephemerides concatenation
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%03d');                                %conversion into a string
    d = dir([fileroot '_eph_' hour_str '.bin']);                    %file to be read
end

%-------------------------------------------------------------------------------

interval = median(time_GPS(2:end) - time_GPS(1:end-1));

delay = time_GPS - roundmod(time_R,interval); %processing delays

loss_R = zeros(length(time_R),1);             %ROVER losses initialization
loss_M = zeros(length(time_M),1);             %MASTER losses initialization

pos = find(time_R == 0);                      %ROVER epochs with losses
loss_R(pos) = 1;                              %flag for ROVER losses
delay(pos) = -1;                              %delay corrections in case of losses

pos = find(time_M == 0);                      %MASTER epochs with losses
loss_M(pos) = 1;              %#ok<FNDSB>     %flag for MASTER losses
