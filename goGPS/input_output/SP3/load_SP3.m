function [SP3] = load_SP3(filename_SP3, filename_CLK, time, week, constellations, date_start, date_stop, wait_dlg)

% SYNTAX:
%   [SP3] = load_SP3(filename_SP3, filename_CLK, time, week, constellations, wait_dlg);
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

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

log = Core.getLogger();

log.addMarkedMessage('Reading SP3s (precise ephemeris) and clocks');

if (isempty(constellations)) %then use only GPS as default
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

% starting index in the total array for the various constellations
idGPS = constellations.GPS.indexes(1);
idGLONASS = constellations.GLONASS.indexes(1);
idGalileo = constellations.Galileo.indexes(1);
idBeiDou = constellations.BeiDou.indexes(1);
idQZSS = constellations.QZSS.indexes(1);

% degree of interpolation polynomial (Lagrange)
n = 10;

% number of seconds in a quarter of an hour
quarter_sec = 900;

if (nargin > 7)
    waitbar(0.5,wait_dlg,'Reading SP3 (precise ephemeris) file...')
end

% extract containing folder
% [data_dir, file_name, file_ext] = fileparts(filename_SP3{1});
% filename_SP3 = File_Name_Processor.checkPath(strcat(data_dir, filesep, file_name(1:3)));

% define time window
[week_start, time_start] = time2weektow(time(1));
[week_end, time_end] = time2weektow(time(end));

% day-of-week
[~, ~, dow_start] = gps2date(week_start, time_start);
[~, ~, dow_end] = gps2date(week_end, time_end);

% add a buffer before and after
% if the first obs is close to the beginning of the week
if (time(1) - weektow2time(week_start, dow_start * 86400, 'G') <= n / 2 * quarter_sec)
    if (dow_start == 0)
        week_start = week_start - 1;
        dow_start = 6;
    else
        dow_start = dow_start - 1;
    end
end

% if the last obs is close to the end of the week
if (time(end) - weektow2time(week_end, dow_end * 86400, 'G') >= 86400 - n / 2 * quarter_sec)
    if (dow_end == 6)
        week_end = week_end + 1;
        dow_end = 0;
    else
        dow_end = dow_end + 1;
    end
else
end

% find the SP3 files dates needed for the processing
% an SP3 file contains data for the entire day, but to interpolate it at
% the beginning and and I also need the previous and the following file
week_dow  = [];
week_curr = week_start;
dow_curr  = dow_start;
while (week_curr <= week_end)
    while ((week_curr < week_end && dow_curr <= 6) || (week_curr == week_end && dow_curr <= dow_end))
        week_dow = [week_dow; week_curr dow_curr];
        dow_curr = dow_curr + 1;
    end
    week_curr = week_curr + 1;
    dow_curr = 0;
end

% init SP3 variables
nEpochs  = 96*size(week_dow,1);
SP3.time = zeros(nEpochs,1);
SP3.coord = zeros(3,constellations.nEnabledSat,nEpochs);
SP3.clock = zeros(constellations.nEnabledSat,nEpochs);
SP3.avail = zeros(constellations.nEnabledSat,1);
SP3.prn   = zeros(constellations.nEnabledSat,1);
SP3.sys   = zeros(constellations.nEnabledSat,1);
SP3.time_hr = [];
SP3.clock_hr = [];
SP3.coord_rate = 900;
SP3.clock_rate = 900;

k = 0; % current epoch
new_file = false;
flag_unavail = 0;

% for each part (SP3 file)
for p = 1 : numel(filename_SP3)

    %SP3 file
    f_sp3 = fopen(filename_SP3{p},'r');

    if p > 1
        new_file = true;
    end
    if (f_sp3 ~= -1)

        % Read the entire sp3 file in memory
        sp3_file = textscan(f_sp3,'%s','Delimiter', '\n');
        if (length(sp3_file) == 1)
            sp3_file = sp3_file{1};
        end
        sp3_cur_line = 1;
        fclose(f_sp3);

        % while there are lines to process
        while (sp3_cur_line <= length(sp3_file))

            % get the next line
            lin = sp3_file{sp3_cur_line};  sp3_cur_line = sp3_cur_line + 1;

            if (strcmp(lin(1:2),'##'))
                rate = str2num(lin(25:38));
                SP3.coord_rate = rate;
                SP3.clock_rate = rate;
            end

            if (lin(1) == '*')

                k = k + 1;

                % read the epoch header
                % example 1: "*  1994 12 17  0  0  0.00000000"
                data   = sscanf(lin(2:31),'%f');
                year   = data(1);
                month  = data(2);
                day    = data(3);
                hour   = data(4);
                minute = data(5);
                second = data(6);

                %computation of the GPS time in weeks and seconds of week
                [w, t] = date2gps([year, month, day, hour, minute, second]);
                %convert GPS time-of-week to continuous time
                if new_file
                    new_k = find(SP3.time(1:k-1,1) >= weektow2time(w, t, 'G'), 1, 'first');
                    if ~isempty(new_k)
                        k = new_k;
                    end
                    new_file = false;
                end
                SP3.time(k,1) = weektow2time(w, t, 'G');
                clear w t;

            elseif (strcmp(lin(1),'P'))
                %read position and clock
                %example 1: "P  1  16258.524750  -3529.015750 -20611.427050    -62.540600"
                %example 2: "PG01   8953.350886  12240.218129 -21918.986611 999999.999999"
                %example 3: "PG02 -13550.970765 -16758.347434 -15825.576525    274.198680  7  8  8 138"
                sys_id = lin(2);
                if (strcmp(sys_id,' ') || strcmp(sys_id,'G') || strcmp(sys_id,'R') || ...
                    strcmp(sys_id,'E') || strcmp(sys_id,'C') || strcmp(sys_id,'J'))

                    index = -1;
                    switch (sys_id)
                        case {'G', ' '}
                            if (constellations.GPS.enabled)
                                index = idGPS;
                            end
                        case 'R'
                            if (constellations.GLONASS.enabled)
                                index = idGLONASS;
                            end
                        case 'E'
                            if (constellations.Galileo.enabled)
                                index = idGalileo;
                            end
                        case 'C'
                            if (constellations.BeiDou.enabled)
                                index = idBeiDou;
                            end
                        case 'J'
                            if (constellations.QZSS.enabled)
                                index = idQZSS;
                            end
                    end

                    % If the considered line is referred to an active constellation
                    if (index >= 0)
                        PRN = sscanf(lin(3:4),'%f');
                        X   = sscanf(lin(5:18),'%f');
                        Y   = sscanf(lin(19:32),'%f');
                        Z   = sscanf(lin(33:46),'%f');
                        clk = sscanf(lin(47:60),'%f');

                        index = index + PRN - 1;

                        SP3.coord(1, index, k) = X*1e3;
                        SP3.coord(2, index, k) = Y*1e3;
                        SP3.coord(3, index, k) = Z*1e3;

                        SP3.clock(index,k) = clk/1e6; % NOTE: clk >= 999999 stands for bad or absent clock values

                        SP3.prn(index) = PRN;
                        SP3.sys(index) = sys_id;

                        if (SP3.clock(index,k) < 0.9)
                            SP3.avail(index) = index;
                        end
                    end
                end
            end
        end
        clear sp3_file;
    end
end

if SP3.time(1) > date_start.getGpsTime()
    log.addWarning(sprintf('Some SP3 files for the processing of the beginning of the interval are missing!!!\n'));
    flag_unavail = 1;
end
if SP3.time(end) < date_stop.getGpsTime()
    log.addWarning(sprintf('Some SP3 files for the processing of the end of the interval are missing!!!\n'));
    flag_unavail = 1;
end

if (~flag_unavail)

    w = zeros(constellations.nEnabledSat,1);
    t = zeros(constellations.nEnabledSat,1);
    clk = zeros(constellations.nEnabledSat,1);
    q = zeros(constellations.nEnabledSat,1);

    % for each part (SP3 file)
    for p = 1 : numel(filename_CLK)

        if strcmp(filename_CLK(p), filename_SP3(p))
            f_clk = -1;
        else
            f_clk = fopen(filename_CLK{p},'r');
        end

        if (f_clk == -1)
            log.addWarning(sprintf('No clk files have been found at %s', filename_CLK{p}));
        else
            [~, fn, fn_ext] = fileparts(filename_CLK{p});
            log.addMessage(sprintf('         Using as clock file: %s%s', fn, fn_ext));
            % read the entire clk file in memory
            clk_file = textscan(f_clk,'%s','Delimiter', '\n');
            if (length(clk_file) == 1)
                clk_file = clk_file{1};
            end
            clk_cur_line = 1;
            fclose(f_clk);

            % while there are lines to process
            while (clk_cur_line <= length(clk_file))

                % get the next line
                lin = clk_file{clk_cur_line};  clk_cur_line = clk_cur_line + 1;

                if (strcmp(lin(1:3),'AS '))

                    sys_id = lin(4);
                    if (strcmp(sys_id,' ') || strcmp(sys_id,'G') || strcmp(sys_id,'R') || ...
                        strcmp(sys_id,'E') || strcmp(sys_id,'C') || strcmp(sys_id,'J'))

                        data = sscanf(lin([5:35 41:59]),'%f'); % sscanf can read multiple data in one step

                        % read PRN
                        PRN = data(1);

                        index = -1;
                        switch (sys_id)
                            case {'G', ' '}
                            if (constellations.GPS.enabled && PRN <= constellations.GPS.numSat)
                                index = idGPS;
                            end
                        case 'R'
                            if (constellations.GLONASS.enabled && PRN <= constellations.GLONASS.numSat)
                                index = idGLONASS;
                            end
                        case 'E'
                            if (constellations.Galileo.enabled && PRN <= constellations.Galileo.numSat)
                                index = idGalileo;
                            end
                        case 'C'
                            if (constellations.BeiDou.enabled && PRN <= constellations.BeiDou.numSat)
                                index = idBeiDou;
                            end
                        case 'J'
                            if (constellations.QZSS.enabled && PRN <= constellations.QZSS.numSat)
                                index = idQZSS;
                            end
                        end

                        % If the considered line is referred to an active constellation
                        if (index >= 0)
                            %read epoch
                            year   = data(2);
                            month  = data(3);
                            day    = data(4);
                            hour   = data(5);
                            minute = data(6);
                            second = data(7);
                            index = index + PRN - 1;
                            q(index) = q(index) + 1;

                            % computation of the GPS time in weeks and seconds of week
                            [w(index,q(index)), t(index,q(index))] = date2gps([year, month, day, hour, minute, second]);
                            clk(index,q(index)) = data(8);
                        end
                    end
                end
            end

            clear clk_file;

            SP3.clock_rate = median(serialize(diff(t(sum(t,2)~=0,:),1,2)));
            % rmndr tells how many observations are needed to reach the end of the last day (criptic)
            rmndr = 86400 / SP3.clock_rate - mod((SP3.time(k,1) - SP3.time(1,1)) / SP3.clock_rate, 86400 / SP3.clock_rate) - 1;
            % SP3.time_hr is the reference clock time from the first value till the end of the day of the last observation
            SP3.time_hr = (SP3.time(1,1) : SP3.clock_rate : (SP3.time(k,1) + rmndr * SP3.clock_rate))';
            SP3.clock_hr = zeros(constellations.nEnabledSat,length(SP3.time_hr));

            s = 1 : constellations.nEnabledSat;
            idx = round((weektow2time(w(s,:), t(s,:), 'G')-SP3.time(1,1)) / SP3.clock_rate) + 1;
            for s = 1 : constellations.nEnabledSat
                SP3.clock_hr(s,idx(s,w(s,:) > 0)) = clk(s,w(s,:) > 0);
            end
        end
    end

    if (SP3.clock_rate >= 60)
        log.addMarkedMessage(sprintf('Satellite clock rate: %f minutes', SP3.clock_rate/60));
    else
        log.addMarkedMessage(sprintf('Satellite clock rate: %f seconds', SP3.clock_rate));
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

if (nargin > 7)
    waitbar(1,wait_dlg)
end
log.newLine();
