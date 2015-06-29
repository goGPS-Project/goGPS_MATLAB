function [SP3] = load_SP3(filename_SP3, time, week, constellations, wait_dlg)

% SYNTAX:
%   [SP3] = load_SP3(filename_SP3, time, week, constellations, wait_dlg);
%
% INPUT:
%   filename_SP3 = SP3 file
%   time = time window (GPS time)
%   week = GPS week
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   SP3 = structure with the following fields:
%      .time  = precise ephemeris timestamps (GPS time)
%      .coord = satellite coordinates  [m]
%      .clock = satellite clock errors [s]
%
% DESCRIPTION:
%   SP3 (precise ephemeris) file parser.
%   NOTE: at the moment the parser reads only time, coordinates and clock;
%         it does not take into account all the other flags and parameters
%         available according to the SP3c format specified in the document
%         http://www.ngs.noaa.gov/orbits/sp3c.txt

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

if (isempty(constellations)) %then use only GPS as default
    [constellations] = multi_constellation_settings(1, 0, 0, 0, 0, 0);
end

%starting index in the total array for the various constellations
idGPS = constellations.GPS.indexes(1);
idGLONASS = constellations.GLONASS.indexes(1);
idGalileo = constellations.Galileo.indexes(1);
idBeiDou = constellations.BeiDou.indexes(1);
idQZSS = constellations.QZSS.indexes(1);

%degree of interpolation polynomial (Lagrange)
n = 10;

%number of seconds in a quarter of an hour
quarter_sec = 900;

if (nargin > 4)
    waitbar(0.5,wait_dlg,'Reading SP3 (precise ephemeris) file...')
end

%extract containing folder
% if (isunix)
%     slash = '/';
% else
%     slash = '\';
% end
pos = strfind(filename_SP3, '/');
if (isempty(pos))
    pos = strfind(filename_SP3, '\');
end
filename_SP3 = filename_SP3(1:pos(end)+3);

%define time window
time_start = time(1);
time_end   = time(end);
week_start = week(1);
week_end   = week(end);

%day-of-week
[~, ~, dow_start] = gps2date(week_start, weektime2tow(week_start, time_start));
[~, ~, dow_end] = gps2date(week_end, weektime2tow(week_end, time_end));

%add a buffer before and after
if (time_start - weektow2time(week_start, dow_start*86400, 'G') < n/2*quarter_sec)
    if (dow_start == 0)
        week_start = week_start - 1;
        dow_start = 6;
    else
        dow_start = dow_start - 1;
    end
else
end

if (time_end - weektow2time(week_end, dow_end*86400, 'G') > 86400-n/2*quarter_sec)
    if (dow_end == 6)
        week_end = week_end + 1;
        dow_end = 0;
    else
        dow_end = dow_end + 1;
    end
else
end

week_dow  = [];
week_curr = week_start;
dow_curr  = dow_start;
while (week_curr <= week_end)
    while ((week_curr < week_end & dow_curr <= 6) | (week_curr == week_end & dow_curr <= dow_end))
        week_dow = [week_dow; week_curr dow_curr];
        dow_curr = dow_curr + 1;
    end
    week_curr = week_curr + 1;
    dow_curr = 0;
end

nEpochs  = 96*size(week_dow,1);
SP3.time = zeros(nEpochs,1);
SP3.coord = zeros(3,constellations.nEnabledSat,nEpochs);
SP3.clock = zeros(constellations.nEnabledSat,nEpochs);
SP3.avail = zeros(constellations.nEnabledSat,1);
SP3.prn   = zeros(constellations.nEnabledSat,1);
SP3.sys   = zeros(constellations.nEnabledSat,1);

k = 0;
flag_unavail = 0;

for p = 1 : size(week_dow,1)
    
    %open SP3 file
    f_sp3 = fopen([filename_SP3 num2str(week_dow(p,1)) num2str(week_dow(p,2)) '.sp3'],'r');
    
    if (f_sp3 ~= -1)
        %skip the SP3 header (first 22 lines)
        for i = 1 : 22
            fgetl(f_sp3);
        end
        
        while (~feof(f_sp3))
            
            %get the next line
            lin = fgetl(f_sp3);
            
            if (strcmp(lin(1),'*'))
                
                k = k + 1;
                
                %read the epoch header
                %example 1: "*  1994 12 17  0  0  0.00000000"
                data   = textscan(lin(2:31),'%f%f%f%f%f%f');
                year   = data{1};
                month  = data{2};
                day    = data{3};
                hour   = data{4};
                minute = data{5};
                second = data{6};

                %computation of the GPS time in weeks and seconds of week
                [week, time] = date2gps([year, month, day, hour, minute, second]);
                SP3.time(k,1) = weektow2time(week, time, 'G');
                
            elseif (strcmp(lin(1),'P'))
                %read position and clock
                %example 1: "P  1  16258.524750  -3529.015750 -20611.427050    -62.540600"
                %example 2: "PG01   8953.350886  12240.218129 -21918.986611 999999.999999"
                %example 3: "PG02 -13550.970765 -16758.347434 -15825.576525    274.198680  7  8  8 138"
                sys_id = lin(2);
                if (strcmp(sys_id,' ') | strcmp(sys_id,'G') | strcmp(sys_id,'R') | strcmp(sys_id,'E') | ...
                    strcmp(sys_id,'C') | strcmp(sys_id,'J'))

                    PRN = sscanf(lin(3:4),'%f');
                    X   = sscanf(lin(5:18),'%f');
                    Y   = sscanf(lin(19:32),'%f');
                    Z   = sscanf(lin(33:46),'%f');
                    clk = sscanf(lin(47:60),'%f');
                    
                    switch (sys_id)
                        case 'G'
                            index = idGPS;
                        case 'R'
                            index = idGLONASS;
                        case 'E'
                            index = idGalileo;
                        case 'C'
                            index = idBeiDou;
                        case 'J'
                            index = idQZSS;
                    end
                    
                    index = index + PRN - 1;
                    
                    SP3.coord(1, index, k) = X*1e3;
                    SP3.coord(2, index, k) = Y*1e3;
                    SP3.coord(3, index, k) = Z*1e3;
                    
                    SP3.clock(index,k) = clk/1e6; %NOTE: clk >= 999999 stands for bad or absent clock values
                    
                    SP3.prn(index) = PRN;
                    SP3.sys(index) = sys_id;
                    
                    if (SP3.clock(index,k) < 0.9)
                        SP3.avail(index) = index;
                    end
                end
            end
        end
        
        %close SP3 file
        fclose(f_sp3);
    else
        fprintf('Missing SP3 file: %s\n', [filename_SP3 num2str(week_dow(p,1)) num2str(week_dow(p,2)) '.sp3']);
        flag_unavail = 1;
    end
end

%if the required SP3 files are not available, stop the execution
if (flag_unavail)
    error('Error: required SP3 files not available.');
end

%remove empty slots
SP3.time(k+1:nEpochs) = [];
SP3.coord(:,:,k+1:nEpochs) = [];
SP3.clock(:,k+1:nEpochs) = [];

if (nargin > 4)
    waitbar(1,wait_dlg)
end
