function [time_GPS, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, ...
          snr_R, snr_M, pos_M, Eph, delay, loss_R, loss_M] = load_goGPSinput (fileroot)

% SYNTAX:
%   [time_GPS, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, ...
%    snr_R, snr_M, pos_M, Eph, delay, loss_R, loss_M] = load_goGPSinput (fileroot);
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
%   snr_R    = ROVER signal-to-noise ratio
%   snr_M    = MASTER signal-to-noise ratio
%   pos_M    = MASTER station coordinates
%   Eph      = matrix of 21 parameters each satellite (MASTER)
%   delay    = delay in observations processing
%   loss_R   = flag for the ROVER loss of signal
%   loss_M   = flag for the MASTER loss of signal
%
% DESCRIPTION:
%   Kalman filter input data reading.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
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

%initialization
time_GPS = [];                         %reference GPS time
time_M = [];                           %MASTER GPS time
time_R = [];                           %ROVER GPS time
pr1_M = [];                            %MASTER code observations
pr1_R = [];                            %ROVER code observations
ph1_M = [];                            %MASTER phase observations
ph1_R = [];                            %ROVER phase observations
snr_M = [];                            %MASTER signal-to-noise ratio
snr_R = [];                            %ROVER signal-to-noise ratio
pos_M = [];                            %MASTER station coordinates

%observations reading
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%02d');                                    %hour index (string)
d = dir([fileroot '_obs_' hour_str '.bin']);                        %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_obs_' hour_str '.bin\n']);
    num_bytes = d.bytes;                                            %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / (3+32*6);                               %file size (number of packets)
    fid_obs = fopen([fileroot '_obs_' hour_str '.bin']);            %file opening
    buf_obs = fread(fid_obs,num_words,'double');                    %file reading
    fclose(fid_obs);                                                %file closing
    time_GPS = [time_GPS; zeros(num_packs,1)];                      %observations concatenation
    time_M   = [time_M;   zeros(num_packs,1)];
    time_R   = [time_R;   zeros(num_packs,1)];
    pr1_M    = [pr1_M     zeros(32,num_packs)];
    pr1_R    = [pr1_R     zeros(32,num_packs)];
    ph1_M    = [ph1_M     zeros(32,num_packs)];
    ph1_R    = [ph1_R     zeros(32,num_packs)];
    snr_M    = [snr_M     zeros(32,num_packs)];
    snr_R    = [snr_R     zeros(32,num_packs)];
    for j = 0 : (3+32*6) : num_words-1
        i = i+1;                                                    %epoch counter increase
        time_GPS(i,1) = buf_obs(j + 1);                             %observations logging
        time_M(i,1)   = buf_obs(j + 2);
        time_R(i,1)   = buf_obs(j + 3);
        pr1_M(:,i)    = buf_obs(j + [4:35]);
        pr1_R(:,i)    = buf_obs(j + [36:67]);
        ph1_M(:,i)    = buf_obs(j + [68:99]);
        ph1_R(:,i)    = buf_obs(j + [100:131]);
        snr_M(:,i)    = buf_obs(j + [132:163]);
        snr_R(:,i)    = buf_obs(j + [164:195]);
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%02d');                                %conversion into a string
    d = dir([fileroot '_obs_' hour_str '.bin']);                    %file to be read
end

%-------------------------------------------------------------------------------

%inizializzazione
Eph = [];

%lettura effemeridi
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%02d');                                    %hour index (string)
d = dir([fileroot '_eph_' hour_str '.bin']);                        %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_eph_' hour_str '.bin\n']);
    num_bytes = d.bytes;                                            %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / (1+672);                                %file size (number of packets)
    fid_eph = fopen([fileroot '_eph_' hour_str '.bin']);            %file opening
    buf_eph = fread(fid_eph,num_words,'double');                    %file reading
    fclose(fid_eph);                                                %file closing
    Eph = cat(3,Eph,zeros(21,32,num_packs));                        %ephemerides concatenation
    for j = 0 : (1+672) : num_words-1
        i = i+1;                                                    %epoch counter increase
        %time_GPS(i,1) = buf_eph(j + 1);                            %GPS time logging
        Eph(:,:,i) = reshape(buf_eph(j + [2:673]), [21,32]);        %ephemerides concatenation
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%02d');                                %conversion into a string
    d = dir([fileroot '_eph_' hour_str '.bin']);                    %file to be read
end

%-------------------------------------------------------------------------------

delay = time_GPS - time_R;                %processing delays

loss_R = zeros(length(time_R),1);         %ROVER losses initialization
loss_M = zeros(length(time_M),1);         %MASTER losses initialization

pos = find(time_R == 0);                  %ROVER epochs with losses
loss_R(pos) = 1;                          %flag for ROVER losses
delay(pos) = -1;                          %delay corrections in case of losses

pos = find(time_M == 0);                  %MASTER epochs with losses
loss_M(pos) = 1;                          %flag for MASTER losses
