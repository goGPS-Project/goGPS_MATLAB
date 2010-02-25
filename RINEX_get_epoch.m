function [time, satG, satS, satR, datee] = RINEX_get_epoch(fid)

% SYNTAX:
%   [time, satG, satS, satR, datee] = RINEX_get_epoch(fid);
%
% INPUT:
%   fid = pointer to the observation RINEX file
%
% OUTPUT:
%   time = observation GPS time
%   satG = list of visible GPS satellites
%   satS = list of visible SBAS satellites
%   satR = list of visible GLONASS satellites
%   datee = date (year,month,day,hour,minute,second)
%
% DESCRIPTION:
%   Scan the first line of each epoch (RINEX) and return
%   the information it contains.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
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
satG = [];
satS = [];
satR = [];
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
        lin = lin(30:end);
        
        %add the second line in case there are more than 12 satellites
        if (num_sat > 12)
            lin_doppler = fgetl(fid);
            lin = [lin, lin_doppler];
        end

        %remove 'blank spaces' and unwanted characters from beginning and end of the string
        while (lin(1) == ' ')
            lin = lin(2:end);
        end
        while (lin(end) == ' ') | (lin(end) == double(10)) | (lin(end) == double(13))
            lin = lin(1:end-1);
        end
        
        %check if GPS satellites are labeled 'G' or not labeled
        %ex. with G:       6G 1G 4G10G12G28G29
        %    not labeled:  6  1  4 10 12 28 29

        %search for 'G'
        [trovatoG posG] = find(lin=='G');
        if (~isempty(trovatoG))
            for k = 1 : length(posG)
                satG = [satG; sscanf(lin(posG(k)+1:posG(k)+2),'%d')];
            end
        else
            %search for 'blank'
            [trovato pos] = find(lin==' ');
            if (~isempty(trovato))
                for k = 1 : length(pos)
                    satG = [satG; sscanf(lin(pos(k)+1:pos(k)+2),'%d')];
                end
            end
        end

        %search for 'S'
        [trovatoS posS] = find(lin=='S');
        if (~isempty(trovatoS))
            for k = 1 : length(posS)
                satS = [satS; sscanf(lin(posS(k)+1:posS(k)+2),'%d')];
            end
        end

        %search for 'R'
        [trovatoR posR] = find(lin=='R');
        if (~isempty(trovatoR))
            for k = 1 : length(posR)
                satR = [satR; sscanf(lin(posR(k)+1:posR(k)+2),'%d')];
            end
        end

        eof = 1;
    end
end

datee = [year month day hour minute second];