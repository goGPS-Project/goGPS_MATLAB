function [CRX, found] = load_crx(data_dir_crx, time, cc)

% SYNTAX:
%   [CRX, found] = load_crx(data_dir_crx, gps_week, time_R, nSatTot, constellations);
%
% INPUT:
%   data_dir_crx = path to the directory containing CRX files [string]
%   gps_week = reference vector of GPS week numbers
%   time_R = reference vector of GPS time
%   cc = Constellation_Collector object, contains the satus of the satellite systems in use
%
% OUTPUT:
%   CRX = matrix containing CRX data
%   found = flag to check if the required file was found
%
% DESCRIPTION:
%   Tool for loading .CRX files: information on satellite problems.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:
%  Contributors:     Andrea Gatti ...
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

%output initialization
%CRX = sparse(false(cc.getNumSat, length(time_R)));
CRX = sparse(false(cc.getNumSat, time.length));
log = Core.getLogger();

%detect starting and ending year/month
date6_start = time.first.get6ColDate();
date6_stop = time.last.get6ColDate();
year_start = date6_start(1);
year_end   = date6_stop(1);
date_start = time.first.getMatlabTime();
date_stop  = time.last.getMatlabTime();

% directory containing CRX files
data_dir = dir(data_dir_crx);

% check the number of files contained in the directory
n_max = size(data_dir,1);

% file counter
n = 0;

% CRX file found
found = 0;

log.addMarkedMessage('Reading CRX data');

% find files with ".CRX" extension
for j = 1 : n_max

    % read the name of the j-th file
    crx_file_name = getfield(data_dir,{j,1},'name');

    % get the number of characters in the filename
    crx_fn_length = size(crx_file_name,2);

    if (crx_fn_length < 12)
        continue
    end

    year = str2num(crx_file_name(5:8)); %#ok<ST2NM>

    % check if the filename corresponds to that expected from a standard CRX file required by goGPS
    % (e.g. "SAT_yyyy.CRX", with 'yyyy' = four-digit year)
    if (crx_fn_length == 12  && (strcmpi(crx_file_name(1:4), 'SAT_') && ...
       (year >=  year_start && year  <=  year_end)  && ...
        strcmpi(crx_file_name(crx_fn_length - 3 : crx_fn_length), '.CRX')))

        n = n + 1;
        log.addMessage(sprintf('       - %s', crx_file_name));

        % full path to the target file
        crx_file_target  = strcat(data_dir_crx, '/', crx_file_name);

        % open .crx file
        fid_fd = fopen(crx_file_target,'rt');

        % warnings
        if (fid_fd ~= -1)
            log.addStatusOk();
        else
            log.addWarning(['impossible to open CRX file - ', crx_file_name]);
            break
        end

        line = fgetl(fid_fd);
        while(~feof(fid_fd) && (isempty(line) || ~strcmp(line(1:5), '  ***')))
            line = fgetl(fid_fd);
        end
        fgetl(fid_fd);

        % for each data line
        while(~feof(fid_fd))
            line = fgetl(fid_fd);

            if (~isempty(line))

                prn = abs(str2num(line(3:5))); %#ok<ST2NM>
                e = [];

                if (~isempty(prn))

                    offset = -1;
                    if (prn < 100) % GPS
                        if cc.getGPS().isActive() && (prn < cc.getGPS().N_SAT)
                            offset = cc.getGPS().getFirstId(); % starting index in the total array for the various constellations
                        end
                    elseif (prn < 200) % GLONASS
                        prn = prn - 100;
                        if cc.getGLONASS().isActive() && (prn < cc.getGLONASS().N_SAT)
                            offset = cc.getGLONASS().getFirstId(); % starting index in the total array for the various constellations
                        end
                    elseif (prn < 300) % Galileo
                        prn = prn - 200;
                        if cc.getGalileo().isActive() && (prn < cc.getGalileo().N_SAT)
                            offset = cc.getGalileo().getFirstId(); % starting index in the total array for the various constellations
                        end
                    elseif (prn < 400) % SBAS
                        prn = prn - 300;
                        if cc.getSBAS().isActive() && (prn < cc.getSBAS().N_SAT)
                            offset = cc.getSBAS().getFirstId(); % starting index in the total array for the various constellations
                        end
                    elseif (prn < 500) % BeiDou
                        prn = prn - 400;
                        if cc.getBeiDou().isActive() && (prn < cc.getBeiDou().N_SAT)
                            offset = cc.getBeiDou().getFirstId(); % starting index in the total array for the various constellations
                        end
                    elseif (prn < 600) % QZSS
                        prn = prn - 500;
                        if cc.getQZSS().isActive() && (prn < cc.getQZSS().N_SAT)
                            offset = cc.getQZSS().getFirstId(); % starting index in the total array for the various constellations
                        end
                    end

                    if (offset >= 0)
                        index = offset + prn - 1;

                        p = str2num(line(12:18)); %problem
                        s = datenum(str2num(line(33:51))); %start date
                        if (length(line) >= 72)
                            e = datenum(str2num(line(54:72))); %end date
                        end
                        if (isempty(e))
                            if (p == 0) %satellite maneuver: remove 15 minutes before and after
                                s = s - datenum([0 0 0 0 15 0]);
                                e = s + datenum([0 0 0 0 15 0]);
                            else
                                e = floor(s) + 1; %arc split: exclude the satellite for the rest of the processing
                                %e = date_stop; %arc split: exclude the satellite for the rest of the processing
                            end
                        end
                        if ((p == 0 &&           (s <= date_stop && e >= date_start)) || ... %satellite maneuver
                                (p >= 1 && p <= 3 && (s <= date_stop && e >= date_start)))   %bad code and/or phase data

                                %(p == 4 &&           (s <= date_stop && e >= date_start))) % arc split <--> not needed by goGPS

                            [~, idx_start] = min(abs(s - time.getMatlabTime()));
                            [~, idx_end]   = min(abs(e - time.getMatlabTime()));
                            CRX(index, idx_start:idx_end) = 1;
                        end
                    end
                end
            end

        end

        fclose(fid_fd);
    end
end

% if no .CRX files are available, return
if (n == 0)
    log.addWarning(['The required (updated) CRX files were not found in ' data_dir_crx ' directory.\n']);
else
    % CRX file found
    found = 1;
    bad_sat = find(any(CRX,2));
    if ~isempty(bad_sat)
        n_bad_epochs = sum(CRX(bad_sat,:)~=0,2);
        for b = 1 : numel(bad_sat)
            log.addWarning(sprintf('%5d bad epochs (sat %c%02d) have been discovered into CRX', nonzeros(n_bad_epochs(b)), cc.system(bad_sat(b)), cc.prn(bad_sat(b)) ));
        end
    end
end

log.newLine();
