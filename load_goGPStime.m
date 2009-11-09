function [dt_acqR, dt_decR, dt_acqM, dt_decM, dt_saveI, dt_kal, dt_saveO, dt_plot, dt_ge, dt_sky, dt_snr] = ...
          load_goGPStime (fileroot)

% SYNTAX:
%   [dt_acqR, dt_decR, dt_acqM, dt_decM, dt_saveI, dt_kal, dt_saveO, dt_plot, dt_ge, dt_sky, dt_snr] = ...
%    load_goGPStime (fileroot);
%
% INPUT:
%   fileroot  = name of the file to be read
%
% OUTPUT:
%   dt_acqR  = ROVER data capture duration
%   dt_decR  = ROVER data decodification duration
%   dt_acqM  = MASTER data capture duration
%   dt_decM  = MASTER data decodification duration
%   dt_saveI = input data saving duration
%   dt_kal   = Kalman filter execution duration
%   dt_saveO = output data saving duration
%   dt_plot  = path plotting in MATLAB duration
%   dt_ge    = path plotting in GoogleEarth duration
%   dt_sky   = sky-plotting duration 
%   dt_snr   = signal-to-noise plotting duration
%
% DESCRIPTION:
%   Estimate of processing time.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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
dt_acqR  = [];                       %processing and plotting duration
dt_decR  = [];
dt_acqM  = [];
dt_decM  = [];
dt_saveI = [];
dt_kal   = [];
dt_saveO = [];
dt_plot  = [];
dt_ge    = [];
dt_sky   = [];
dt_snr   = [];

%observations reading
i = 0;                                                              %epoch counter
hour = 0;                                                           %hour index (integer)
hour_str = num2str(hour,'%02d');                                    %hour index (string)
d = dir([fileroot '_dt_' hour_str '.bin']);                         %file to be read
while ~isempty(d)
    fprintf(['Reading: ' fileroot '_dt_' hour_str '.bin\n']);
    num_bytes = d.bytes;                                            %file size (number of bytes)
    num_words = num_bytes / 8;                                      %file size (number of words)
    num_packs = num_words / 11;                                     %file size (number of packets)
    fid_dt = fopen([fileroot '_dt_' hour_str '.bin'],'r+');         %file opening
    buf_dt = fread(fid_dt,num_words,'double');                      %file reading
    fclose(fid_dt);                                                 %file closing
    dt_acqR  = [dt_acqR;   zeros(num_packs,1)];                     %observations concatenation
    dt_decR  = [dt_decR;   zeros(num_packs,1)];
    dt_acqM  = [dt_acqM;   zeros(num_packs,1)];
    dt_decM  = [dt_decM;   zeros(num_packs,1)];
    dt_saveI = [dt_saveI;  zeros(num_packs,1)];
    dt_kal   = [dt_kal;    zeros(num_packs,1)];
    dt_saveO = [dt_saveO;  zeros(num_packs,1)];
    dt_plot  = [dt_plot;   zeros(num_packs,1)];
    dt_ge    = [dt_ge;     zeros(num_packs,1)];
    dt_sky   = [dt_sky;    zeros(num_packs,1)];
    dt_snr   = [dt_snr;    zeros(num_packs,1)];
    for j = 0 : 11 : num_words-1
        i = i+1;                                                    %epoch counter increase
        dt_acqR(i,1)  = buf_dt(j + 1);                              %observations logging
        dt_decR(i,1)  = buf_dt(j + 2);
        dt_acqM(i,1)  = buf_dt(j + 3);
        dt_decM(i,1)  = buf_dt(j + 4);
        dt_saveI(i,1) = buf_dt(j + 5);
        dt_kal(i,1)   = buf_dt(j + 6);
        dt_saveO(i,1) = buf_dt(j + 7);
        dt_plot(i,1)  = buf_dt(j + 8);
        dt_ge(i,1)    = buf_dt(j + 9);
        dt_sky(i,1)   = buf_dt(j + 10);
        dt_snr(i,1)   = buf_dt(j + 11);
    end
    hour = hour+1;                                                  %hour increase
    hour_str = num2str(hour,'%02d');                                %conversion into a string
    d = dir([fileroot '_dt_' hour_str '.bin']);                     %file to be read
end

%-------------------------------------------------------------------------------