function [Xhat_t_t, Yhat_t_t, Cee, azM, azR, elM, elR, distM, distR, ...
          conf_sat, conf_cs, pivot, PDOP, HDOP, VDOP, KPDOP, ...
          KHDOP, KVDOP, RES_CODE_FIX, RES_PHASE_FIX, RES_CODE_FLOAT, RES_PHASE_FLOAT, outliers_CODE, outliers_PHASE]= load_goGPSoutput (fileroot, mode, mode_vinc)

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

%global variables loading
global o1 o3 nN

%-------------------------------------------------------------------------------

%initialization
Xhat_t_t = [];                     %state variables estimate
Yhat_t_t = [];                     %receiver positions estimate
Cee = [];                          %estimation error covariance matrix

%observations reading
if (mode == 14 & mode_vinc == 1)
    i = 0;                                                              %epoch counter
    hour = 0;                                                           %hour index (integer)
    hour_str = num2str(hour,'%03d');                                    %hour index (string)
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
        hour_str = num2str(hour,'%03d');                                %conversion into a string
        d = dir([fileroot '_kal_' hour_str '.bin']);                    %file to be read
    end
else
    i = 0;                                                              %epoch counter
    hour = 0;                                                           %hour index (integer)
    hour_str = num2str(hour,'%03d');                                    %hour index (string)
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
        hour_str = num2str(hour,'%03d');                                %conversion into a string
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
hour_str = num2str(hour,'%03d');                                    %hour index (string)
d = dir([fileroot '_sat_' hour_str '.bin']);                        %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_sat_' hour_str '.bin\n']);
    fid_sat = fopen([fileroot '_sat_' hour_str '.bin'],'r+');       %file opening
    num_sat = fread(fid_sat,1,'int8');                              %read number of satellites
    num_bytes = d.bytes-1;                                          %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / (num_sat*6);                            %file size (number of packets)
    buf_sat = fread(fid_sat,num_words,'double');                    %file reading
    fclose(fid_sat);                                                %file closing
    azM   = [azM    zeros(num_sat,num_packs)];                      %observations concatenation
    azR   = [azR    zeros(num_sat,num_packs)];
    elM   = [elM    zeros(num_sat,num_packs)];
    elR   = [elR    zeros(num_sat,num_packs)];
    distM = [distM  zeros(num_sat,num_packs)];
    distR = [distR  zeros(num_sat,num_packs)];
    for j = 0 : (num_sat*6) : num_words-1
        i = i+1;                                                    %epoch counter increase
        azM(:,i)   = buf_sat(j + [1:num_sat]);                      %observations logging
        azR(:,i)   = buf_sat(j + [1*num_sat+1:2*num_sat]);
        elM(:,i)   = buf_sat(j + [2*num_sat+1:3*num_sat]);
        elR(:,i)   = buf_sat(j + [3*num_sat+1:4*num_sat]);
        distM(:,i) = buf_sat(j + [4*num_sat+1:5*num_sat]);
        distR(:,i) = buf_sat(j + [5*num_sat+1:6*num_sat]);
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%03d');                                %conversion into a string
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
hour_str = num2str(hour,'%03d');                                    %hour index (string)
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
    hour_str = num2str(hour,'%03d');                                %conversion into a string
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
hour_str = num2str(hour,'%03d');                                    %hour index (string)
d = dir([fileroot '_conf_' hour_str '.bin']);                       %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_conf_' hour_str '.bin\n']);
    fid_conf = fopen([fileroot '_conf_' hour_str '.bin'],'r+');     %file opening
    num_sat = fread(fid_conf,1,'int8');                             %read number of satellites
    num_bytes = d.bytes-1;                                          %file size (number of bytes)
    num_packs = num_bytes / (num_sat*2+1);                          %file size (number of packets)
    buf_conf = fread(fid_conf,num_bytes,'int8');                    %file reading
    fclose(fid_conf);                                               %file closing
    conf_sat = [conf_sat  zeros(num_sat,num_packs)];                %observations concatenation
    conf_cs  = [conf_cs   zeros(num_sat,num_packs)];
    pivot    = [pivot;    zeros(num_packs,1)];
    for j = 0 : (num_sat*2+1) : num_bytes-1
        i = i+1;                                                    %epoch counter increase
        conf_sat(:,i) = buf_conf(j + [1:num_sat]);                  %observations logging
        conf_cs(:,i)  = buf_conf(j + [1*num_sat+1:2*num_sat]);
        pivot(i,1)    = buf_conf(j + 2*num_sat+1);
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%03d');                                %conversion into a string
    d = dir([fileroot '_conf_' hour_str '.bin']);                   %file to be read
end

%-------------------------------------------------------------------------------
%initialization
RES_CODE_FIX  = [];                      %double differences code residuals (fixed solution)
RES_PHASE_FIX = [];                      %phase differences phase residuals (fixed solution)
RES_CODE_FLOAT  = [];                    %double differences code residuals (float solution)
RES_PHASE_FLOAT = [];                    %phase differences phase residuals (float solution)
outliers_CODE = [];                      %code double difference outlier? (fixed solution)
outliers_PHASE = [];                     %phase double difference outlier? (fixed solution)
%observations reading
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%03d');                                    %hour index (string)
d = dir([fileroot '_res_' hour_str '.bin']);                        %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_res_' hour_str '.bin\n']);
    fid_sat = fopen([fileroot '_res_' hour_str '.bin'],'r+');       %file opening
    num_sat = fread(fid_sat,1,'int8');                              %read number of satellites
    num_bytes = d.bytes-1;                                          %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / (num_sat*6);                            %file size (number of packets)
    buf_sat = fread(fid_sat,num_words,'double');                    %file reading
    fclose(fid_sat);                                                %file closing
    RES_CODE_FIX    = [RES_CODE_FIX    zeros(num_sat,num_packs)];   %observations concatenation
    RES_PHASE_FIX   = [RES_PHASE_FIX   zeros(num_sat,num_packs)];
    RES_CODE_FLOAT  = [RES_CODE_FLOAT  zeros(num_sat,num_packs)];
    RES_PHASE_FLOAT = [RES_PHASE_FLOAT zeros(num_sat,num_packs)];
    outliers_CODE   = [outliers_CODE   zeros(num_sat,num_packs)];
    outliers_PHASE  = [outliers_PHASE  zeros(num_sat,num_packs)];
    for j = 0 : (num_sat*6) : num_words-1
        i = i+1;                                                    %epoch counter increase
        RES_CODE_FIX(:,i)    = buf_sat(j + [1:num_sat]);            %observations logging
        RES_PHASE_FIX(:,i)   = buf_sat(j + [1*num_sat+1:2*num_sat]);
        RES_CODE_FLOAT(:,i)  = buf_sat(j + [2*num_sat+1:3*num_sat]);
        RES_PHASE_FLOAT(:,i) = buf_sat(j + [3*num_sat+1:4*num_sat]);
        outliers_CODE(:,i)   = buf_sat(j + [4*num_sat+1:5*num_sat]);
        outliers_PHASE(:,i)  = buf_sat(j + [5*num_sat+1:6*num_sat]);
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%03d');                                %conversion into a string
    d = dir([fileroot '_res_' hour_str '.bin']);                    %file to be read
end
