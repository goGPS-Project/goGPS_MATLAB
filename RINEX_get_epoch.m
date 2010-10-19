function [time, sat, sat_types, datee] = RINEX_get_epoch(fid)

% SYNTAX:
%   [time, sat, sat_types, datee] = RINEX_get_epoch(fid);
%
% INPUT:
%   fid = pointer to the observation RINEX file
%
% OUTPUT:
%   time = observation GPS time
%   sat  = list of all visible satellites
%   sat_types = ordered list of satellite types ('G' = GPS, 'R' = GLONASS, 'S' = SBAS)
%   datee = date (year,month,day,hour,minute,second)
%
% DESCRIPTION:
%   Scan the first line of each epoch (RINEX) and return
%   the information it contains.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.1 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
%
% Partially based on FEPOCH_0.M (EASY suite) by Kai Borre
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

%variable initialization
time = 0;
sat = [];
sat_types = [];
datee=[0 0 0 0 0 0];
eof = 0;

%search data
while (eof==0)
    %read the string
    lin = fgets(fid);
    answer = findstr(lin,'COMMENT');
    %if it is a comment line read the following one
    if ~isempty(answer)
        lin = fgetl(fid);
    end
    %check if the end of file is reached
    if (feof(fid) == 1);
        return
    end
    %check if it is a string that should be analyzed
    if (strcmp(lin(29),'0') == 1) | (strcmp(lin(29),'1') == 1)
    % if lin(2) ~= ' '

        %save line information
        data = textscan(lin(1:26),'%f%f%f%f%f%f');
        year = data{1};
        month = data{2};
        day = data{3};
        hour = data{4};
        minute = data{5};
        second = data{6};
        %computation of the corresponding julian day
        jd = julday(year+2000, month, day, 0);
        %computation of the GPS time in weeks and seconds of week
        [week, sec_of_week] = gps_time(jd); %#ok<ASGLU>
        time = sec_of_week + hour*3600+minute*60+second;

        %number of visible satellites
        [num_sat] = sscanf(lin(30:32),'%d');

        %keep just the satellite data
        lin = lin(33:end);
        
        %remove 'blank spaces' and unwanted characters at the end of the string
        while (lin(end) == ' ') | (lin(end) == double(10)) | (lin(end) == double(13))
            lin = lin(1:end-1);
        end
        
        %add the second line in case there are more than 12 satellites
        if (num_sat > 12)
            lin_add = fgetl(fid);
            lin_add = lin_add(33:end);
            lin = [lin lin_add];
        end

        pos = 1;
        for i = 1 : num_sat
            %check if GPS satellites are labeled 'G' or not labeled
            if (strcmp(lin(pos),' '))
                type = 'G';
            else
                type = lin(pos);
            end
            sat_types = [sat_types; type];
            sat = [sat; sscanf(lin(pos+1:pos+2),'%d')];
            pos = pos + 3;
        end

        eof = 1;
    end
end

datee = [year month day hour minute second];