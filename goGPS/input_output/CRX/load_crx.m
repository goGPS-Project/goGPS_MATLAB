function [CRX] = load_crx(data_dir_crx, gps_week, time_R)

% SYNTAX:
%   [CRX] = load_crx(data_dir_crx, gps_week, time_R);
%
% INPUT:
%   data_dir_crx = path to the directory containing CRX files [string]
%   gps_week = reference vector of GPS week numbers
%   time_R = reference vector of GPS time
%
% OUTPUT:
%   CRX = matrix containing CRX data
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

%output initialization
CRX = zeros(nSatTot, length(time_R));

%convert GPS time to time-of-week
gps_tow = weektime2tow(gps_week, time_R);

%detect starting and ending year/month
date = gps2date(gps_week, gps_tow);
year_start  = two_digit_year(date(1,1));
year_end    = two_digit_year(date(end,1));

%directory containing CRX files
data_dir = dir(data_dir_crx);

%check the number of files contained in the directory
nmax = size(data_dir,1);

%file counter
n = 0;

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
                fprintf(['Reading CRX files...\n']);
            end
        else
            fprintf(['WARNING: impossible to open CRX file ', crx_file_name, '\n']);
            break
        end
        
        line = '';
        while(~feof(fid_fd) && ~strcmp(line(1:5), '  ***'))
            line = fgetl(fid_fd);
        end
        fgetl(fid_fd);
        
        while(~feof(fid_fd))
            line = fgetl(fid_fd);
            
            if (isempty(line))
                continue
            end
            
            PRN = str2num(line(3:5)); %#ok<ST2NM>
            if (PRN <= 32)
                
            end
        end
        
        fclose(fid_fd);
    end
end

%if no .CRX files are available, return
if (n == 0)
    fprintf(['The required (updated) CRX files were not found in ' data_dir_crx ' directory.\n'])
    return
end
