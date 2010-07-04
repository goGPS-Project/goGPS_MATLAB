function [Xhat_t_t, Yhat_t_t, Cee, azM, azR, elM, elR, distM, distR, ...
          conf_sat, conf_cs, pivot, PDOP, HDOP, VDOP, KPDOP, ...
          KHDOP, KVDOP]= load_goGPSoutput (fileroot, mode, mode_vinc)

% SYNTAX:
%   [Xhat_t_t, Yhat_t_t, Cee, azM, azR, elM, elR, distM, distR, ...
%    conf_sat, conf_cs, pivot, PDOP, HDOP, VDOP, KPDOP, ...
%    KHDOP, KVDOP]= load_goGPSoutput (fileroot, mode, mode_vinc);
%
% INPUT:
%   fileroot  = name of the file to be read
%   mode      = functioning mode
%   mode_vinc = navigation mode (free=0, constrained=1)
%
% OUTPUT:
%   Xhat_t_t = state variables estimate
%   Yhat_t_t = receiver positions estimate (only constrained path)
%   Cee      = estimation error covariance matrix
%   azM      = satellite - MASTER azimuth
%   azR      = satellite - ROVER azimuth
%   elM      = satellite - MASTER elevation
%   elR      = satellite - ROVER elevation
%   distM    = satellite - MASTER distance
%   distR    = satellite - ROVER distance
%   conf_sat = satellites-in-view configuration
%   conf_cs  = cycle-slips configuration
%   pivot    = pivot satellite
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%   KPDOP = Kalman filter PDOP
%   KHDOP = Kalman filter HDOP
%   KVDOP = Kalman filter VDOP
%
% DESCRIPTION:
%   Kalman filter output data reading.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
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

%global variables loading
global o1 o3 nN

%-------------------------------------------------------------------------------

%initialization
Xhat_t_t = [];                     %state variables estimate
Yhat_t_t = [];                     %receiver positions estimate
Cee = [];                          %estimation error covariance matrix

%observations reading
if (mode == 1 & mode_vinc == 1)
    i = 0;                                                              %epoch counter
    hour = 0;                                                           %hour index (integer)
    hour_str = num2str(hour,'%02d');                                    %hour index (string)
    d = dir([fileroot '_kal_' hour_str '.bin']);                        %file to be read
    while ~isempty(d)
        fprintf(['Reading: ' fileroot '_kal_' hour_str '.bin\n']);
        num_bytes = d.bytes;                                            %file size (number of bytes)
        num_words = num_bytes / 8;                                      %file size (number of words)
        dim_packs = (o1+nN)+3+(o1+nN)^2;                                %packets size
        num_packs = num_words / dim_packs;                              %file size (number of packets)
        fid_kal = fopen([fileroot '_kal_' hour_str '.bin'],'r+');       %file opening
        buf_kal = fread(fid_kal,num_words,'double');                    %file reading
        fclose(fid_kal);                                                %file closing
        Xhat_t_t = [Xhat_t_t  zeros(o1+nN,num_packs)];                  %observations concatenation
        Yhat_t_t = [Yhat_t_t  zeros(3,num_packs)];
        Cee = cat(3,Cee,zeros(o1+nN,o1+nN,num_packs));
        for j = 0 : dim_packs : num_words-1
            i = i+1;                                                    %epoch counter increase
            Xhat_t_t(:,i) = buf_kal(j + [1:o1+nN]);                     %observations logging
            Yhat_t_t(:,i) = buf_kal(j + [o1+nN+1:o1+nN+3]);
            Cee(:,:,i) = reshape(buf_kal(j + [o1+nN+4:dim_packs]), o1+nN, o1+nN);
        end
        hour = hour+1;                                                  %hour increase
        hour_str = num2str(hour,'%02d');                                %conversion into a string
        d = dir([fileroot '_kal_' hour_str '.bin']);                    %file to be read
    end
else
    i = 0;                                                              %epoch counter
    hour = 0;                                                           %hour index (integer)
    hour_str = num2str(hour,'%02d');                                    %hour index (string)
    d = dir([fileroot '_kal_' hour_str '.bin']);                        %file to be read
    while ~isempty(d)
        fprintf(['Reading: ' fileroot '_kal_' hour_str '.bin\n']);
        num_bytes = d.bytes;                                            %file size (number of bytes)
        num_words = num_bytes / 8;                                      %file size (number of words)
        dim_packs = (o3+nN)+(o3+nN)^2;                                  %packets size
        num_packs = num_words / dim_packs;                              %file size (number of packets)
        fid_kal = fopen([fileroot '_kal_' hour_str '.bin'],'r+');       %file opening
        buf_kal = fread(fid_kal,num_words,'double');                    %file reading
        fclose(fid_kal);                                                %file closing
        Xhat_t_t = [Xhat_t_t  zeros(o3+nN,num_packs)];                  %observations concatenation
        Cee = cat(3,Cee,zeros(o3+nN,o3+nN,num_packs));
        for j = 0 : dim_packs : num_words-1
            i = i+1;                                                    %epoch counter increase
            Xhat_t_t(:,i) = buf_kal(j + [1:o3+nN]);                     %observations logging
            Cee(:,:,i) = reshape(buf_kal(j + [o3+nN+1:dim_packs]), o3+nN, o3+nN);
        end
        hour = hour+1;                                                  %hour increase
        hour_str = num2str(hour,'%02d');                                %conversion into a string
        d = dir([fileroot '_kal_' hour_str '.bin']);                    %file to be read
    end
end

%-------------------------------------------------------------------------------

%initialization
azM = [];                            %satellite - MASTER azimuth
azR = [];                            %satellite - ROVER azimuth
elM = [];                            %satellite - MASTER elevation
elR = [];                            %satellite - ROVER elevation
distM = [];                          %satellite - MASTER distance
distR = [];                          %satellite - ROVER distance

%observations reading
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%02d');                                    %hour index (string)
d = dir([fileroot '_sat_' hour_str '.bin']);                        %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_sat_' hour_str '.bin\n']);
    num_bytes = d.bytes;                                            %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / (32*6);                                 %file size (number of packets)
    fid_sat = fopen([fileroot '_sat_' hour_str '.bin'],'r+');       %file opening
    buf_sat = fread(fid_sat,num_words,'double');                    %file reading
    fclose(fid_sat);                                                %file closing
    azM   = [azM    zeros(32,num_packs)];                           %observations concatenation
    azR   = [azR    zeros(32,num_packs)];
    elM   = [elM    zeros(32,num_packs)];
    elR   = [elR    zeros(32,num_packs)];
    distM = [distM  zeros(32,num_packs)];
    distR = [distR  zeros(32,num_packs)];
    for j = 0 : (32*6) : num_words-1
        i = i+1;                                                    %epoch counter increase
        azM(:,i)   = buf_sat(j + [1:32]);                           %observations logging
        azR(:,i)   = buf_sat(j + [33:64]);
        elM(:,i)   = buf_sat(j + [65:96]);
        elR(:,i)   = buf_sat(j + [97:128]);
        distM(:,i) = buf_sat(j + [129:160]);
        distR(:,i) = buf_sat(j + [161:192]);
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%02d');                                %conversion into a string
    d = dir([fileroot '_sat_' hour_str '.bin']);                    %file to be read
end

%-------------------------------------------------------------------------------

%initialization
PDOP = [];                           %position dilution of precision
HDOP = [];                           %horizontal dilution of precision
VDOP = [];                           %vertical dilution of precision
KPDOP = [];                          %Kalman filter PDOP
KHDOP = [];                          %Kalman filter HDOP
KVDOP = [];                          %Kalman filter VDOP

%observations reading
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%02d');                                    %hour index (string)
d = dir([fileroot '_dop_' hour_str '.bin']);                        %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_dop_' hour_str '.bin\n']);
    num_bytes = d.bytes;                                            %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / 6;                                      %file size (number of packets)
    fid_dop = fopen([fileroot '_dop_' hour_str '.bin'],'r+');       %file opening
    buf_dop = fread(fid_dop,num_words,'double');                    %file reading
    fclose(fid_dop);                                                %file closing
    PDOP  = [PDOP;  zeros(num_packs,1)];                            %observations concatenation
    HDOP  = [HDOP;  zeros(num_packs,1)];                            %observations concatenation
    VDOP  = [VDOP;  zeros(num_packs,1)];                            %observations concatenation
    KPDOP  = [KPDOP;  zeros(num_packs,1)];                          %observations concatenation
    KHDOP  = [KHDOP;  zeros(num_packs,1)];                          %observations concatenation
    KVDOP  = [KVDOP;  zeros(num_packs,1)];                          %observations concatenation
    for j = 0 : 6 : num_words-1
        i = i+1;                                                    %epoch counter increase
        PDOP(i,1)  = buf_dop(j + 1);                                %observations logging
        HDOP(i,1)  = buf_dop(j + 2);                                %observations logging
        VDOP(i,1)  = buf_dop(j + 3);                                %observations logging
        KPDOP(i,1)  = buf_dop(j + 4);                               %observations logging
        KHDOP(i,1)  = buf_dop(j + 5);                               %observations logging
        KVDOP(i,1)  = buf_dop(j + 6);                               %observations logging
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%02d');                                %conversion into a string
    d = dir([fileroot '_dop_' hour_str '.bin']);                    %file to be read
end

%-------------------------------------------------------------------------------

%initialization
conf_sat = [];                          %satellites-in-view configuration
conf_cs = [];                           %cycle-slips configuration
pivot = [];                             %pivot satellite

%observations reading
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%02d');                                    %hour index (string)
d = dir([fileroot '_conf_' hour_str '.bin']);                       %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_conf_' hour_str '.bin\n']);
    num_bytes = d.bytes;                                            %file size (number of bytes)
    num_packs = num_bytes / (32*2+1);                               %file size (number of packets)
    fid_conf = fopen([fileroot '_conf_' hour_str '.bin'],'r+');     %file opening
    buf_conf = fread(fid_conf,num_bytes,'int8');                    %file reading
    fclose(fid_conf);                                               %file closing
    conf_sat = [conf_sat  zeros(32,num_packs)];                     %observations concatenation
    conf_cs  = [conf_cs   zeros(32,num_packs)];
    pivot    = [pivot;    zeros(num_packs,1)];
    for j = 0 : (32*2+1) : num_bytes-1
        i = i+1;                                                    %epoch counter increase
        conf_sat(:,i) = buf_conf(j + [1:32]);                       %observations logging
        conf_cs(:,i)  = buf_conf(j + [33:64]);
        pivot(i,1)    = buf_conf(j + 65);
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%02d');                                %conversion into a string
    d = dir([fileroot '_conf_' hour_str '.bin']);                   %file to be read
end

%-------------------------------------------------------------------------------