function [CRX, found] = load_crx(data_dir_crx, gps_week, time_R, nSatTot, constellations)

% SYNTAX:
%   [CRX, found] = load_crx(data_dir_crx, gps_week, time_R, nSatTot, constellations);
%
% INPUT:
%   data_dir_crx = path to the directory containing CRX files [string]
%   gps_week = reference vector of GPS week numbers
%   time_R = reference vector of GPS time
%   nSatTot = total number of satellites
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   CRX = matrix containing CRX data
%   found = flag to check if the required file was found
%
% DESCRIPTION:
%   Tool for loading .CRX files: information on satellite problems.

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
idSBAS = constellations.SBAS.indexes(1);

%output initialization
CRX = zeros(nSatTot, length(time_R));

%convert GPS time to time-of-week
gps_tow = weektime2tow(gps_week, time_R);

%detect starting and ending year/month
date = gps2date(gps_week, gps_tow);
year_start = date(1,1);
year_end   = date(end,1);
dnum = datenum(date);
date_start = dnum(1);
date_end   = dnum(end);

%directory containing CRX files
data_dir = dir(data_dir_crx);

%check the number of files contained in the directory
nmax = size(data_dir,1);

%file counter
n = 0;

%CRX file found
found = 0;

%find files with ".CRX" extension
for j = 1 : nmax
    
    %read the name of the j-th file
    crx_file_name = getfield(data_dir,{j,1},'name');
    
    %get the number of characters in the filename
    crx_fn_length = size(crx_file_name,2);

    if (crx_fn_length < 12)
        continue
    end
    
    year = str2num(crx_file_name(5:8)); %#ok<ST2NM>
    
    %check if the filename corresponds to that expected from a standard CRX file required by goGPS
    % (e.g. "SAT_yyyy.CRX", with 'yyyy' = four-digit year)
    if (crx_fn_length == 12  && (strcmpi(crx_file_name(1:4), 'SAT_') && ...
       (year >=  year_start && year  <=  year_end)  && ...
        strcmpi(crx_file_name(crx_fn_length - 3 : crx_fn_length), '.CRX')))
        
        n = n + 1;
        
        %full path to the target file
        crx_file_target  = strcat(data_dir_crx, '/', crx_file_name);
        
        %open .crx file
        fid_fd = fopen(crx_file_target,'r');
        
        %warnings
        if (fid_fd ~= -1)
            if (n == 1)
                fprintf('Reading CRX files...\n');
            end
        else
            fprintf(['WARNING: impossible to open CRX file ', crx_file_name, '\n']);
            break
        end
        
        line = fgetl(fid_fd);
        while(~feof(fid_fd) && (isempty(line) || ~strcmp(line(1:5), '  ***')))
            line = fgetl(fid_fd);
        end
        fgetl(fid_fd);
        
        while(~feof(fid_fd))
            line = fgetl(fid_fd);
            
            if (isempty(line))
                continue
            end
            
            PRN = abs(str2num(line(3:5))); %#ok<ST2NM>
            e = [];
            
            if (isempty(PRN))
                continue
            end
            
            if (PRN <= 0+constellations.GPS.numSat && constellations.GPS.enabled) %GPS
                index = idGPS;
            elseif (PRN > 100 && PRN <= 100+constellations.GLONASS.numSat && constellations.GLONASS.enabled) %GLONASS
                index = idGLONASS;
            elseif (PRN > 200 && PRN <= 200+constellations.Galileo.numSat && constellations.Galileo.enabled) %Galileo
                index = idGalileo;
            elseif (PRN > 300 && PRN <= 300+constellations.SBAS.numSat && constellations.SBAS.enabled)    %SBAS
                index = idSBAS;
            elseif (PRN > 400 && PRN <= 400+constellations.BeiDou.numSat && constellations.BeiDou.enabled)  %BeiDou
                index = idBeiDou;
            elseif (PRN > 500 && PRN <= 500+constellations.QZSS.numSat && constellations.QZSS.enabled)    %QZSS
                index = idQZSS;
            else
                continue
            end
            
            index = index + PRN - 1;
            
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
                    e = date_end; %arc split: exclude the satellite for the rest of the processing
                end
            end
            if ((p == 0 && ((s >= date_start && s <= date_end) || (e >= date_start && e <= date_end))) || ... %satellite maneuver
                    (p >= 1 && p <= 3 && (s <= date_end && e >= date_start)) || ... %bad code and/or phase data
                    (p == 4 && s >= date_start && e <= date_end)) % arc split
                [~, idx_start] = min(abs(s - dnum));
                [~, idx_end]   = min(abs(e - dnum));
                CRX(index, idx_start:idx_end) = 1;
            end
        end
        
        fclose(fid_fd);
    end
end

%if no .CRX files are available, return
if (n == 0)
    fprintf(['The required (updated) CRX files were not found in ' data_dir_crx ' directory.\n'])
    return
else
    %CRX file found
    found = 1;
end
