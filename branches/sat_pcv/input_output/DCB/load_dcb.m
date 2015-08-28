function [DCB] = load_dcb(data_dir_dcb, gps_week, time_R, constellations)

% SYNTAX:
%   [DCB] = load_dcb(data_dir_dcb, gps_week, time_R, constellations);
%
% INPUT:
%   data_dir_dcb = path to the directory containing DCB files [string]
%   gps_week = reference vector of GPS week numbers
%   time_R = reference vector of GPS time
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   DCB = matrix containing DCB data
%
% DESCRIPTION:
%   Tool for loading .DCB files DCB data.

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

%output initialization
DCB = [];

%convert GPS time to time-of-week
gps_tow = weektime2tow(gps_week, time_R);

%detect starting and ending year/month
date = gps2date(gps_week, gps_tow);
year_start  = two_digit_year(date(1,1));
year_end    = two_digit_year(date(end,1));
month_start = date(1,2);
month_end   = date(end,2);

%directory containing DCB files
data_dir = dir(data_dir_dcb);

%check the number of files contained in the directory
nmax = size(data_dir,1);

%file counter
n = 0;

%find files with ".DCB" extension
for j = 1 : nmax
    
    %read the name of the j-th file
    dcb_file_name = getfield(data_dir,{j,1},'name');
    
    %get the number of characters in the filename
    dcb_fn_length = size(dcb_file_name,2);

    if (dcb_fn_length < 12)
        continue
    end
    
    year = str2num(dcb_file_name(5:6));
    month = str2num(dcb_file_name(7:8));
    
    %check if the filename corresponds to that expected from a standard DCB file required by goGPS (e.g. "P1C1xxyy.DCB",
    % with 'xx' = two-digit year and 'yy' = two-digit month)
    if (dcb_fn_length == 12  && strcmpi(dcb_file_name(1:4), 'P1P2') && ...
       ((year >  year_start && year  <  year_end)    || ...
        (year == year_start && month >= month_start) || ...
        (year == year_end   && month <= month_end))  && ...
         strcmpi(dcb_file_name(dcb_fn_length - 3 : dcb_fn_length), '.DCB')) %#ok<*ST2NM>
        
        n = n + 1;
        
        %full path to the target file
        dcb_file_target  = strcat(data_dir_dcb, '/', dcb_file_name);
        
        %open .dcb file
        fid_fd = fopen(dcb_file_target,'r');
        
        %warnings
        if (fid_fd ~= -1)
            %fprintf(['Reading DCB file ', dcb_file_name, '\n']);
            if (n == 1)
                fprintf(['Reading DCB files...\n']);
            end
        else
            fprintf(['WARNING: impossible to open DCB file ', dcb_file_name, '\n']);
            break
        end
        
        line = '';
        while(~feof(fid_fd) && ~strcmp(line, '***   ****************    *****.***   *****.***'))
            line = fgetl(fid_fd);
        end
        
        while(~feof(fid_fd))
            line = fgetl(fid_fd);
            
            if (isempty(line))
                continue
            end
            
            sys_id = line(1);
            if (strcmp(sys_id,'G') && constellations.GPS.enabled || ...
                strcmp(sys_id,'R') && constellations.GLONASS.enabled || ...
                strcmp(sys_id,'E') && constellations.Galileo.enabled || ...
                strcmp(sys_id,'C') && constellations.BeiDou.enabled || ...
                strcmp(sys_id,'J') && constellations.QZSS.enabled)
                
                PRN   = sscanf(line(2:3),'%f');
                value = sscanf(line(30:35),'%f');
                rms   = sscanf(line(43:47),'%f');
                
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
                
                [w, s] = date2gps([four_digit_year(year) month 15 0 0 0]);

                DCB.time(n, 1)      = weektow2time(w, s, 'G');
                DCB.value(index, n) = value;
                DCB.rms(index, n)   = rms;
                DCB.prn(index)      = PRN;
                DCB.sys(index)      = sys_id;
            end
        end
        
        fclose(fid_fd);
    end
end

%if no .DCB files are available, return
if (isempty(DCB))
    fprintf(['No DCB files found in ' data_dir_dcb ' directory.\n'])
    return
end
