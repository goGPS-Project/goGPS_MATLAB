function [time, datee, num_sat, sat, sat_types, tow] = RINEX_get_epoch(fid)

% SYNTAX:
%   [time, datee, num_sat, sat, sat_types, tow] = RINEX_get_epoch(fid);
%
% INPUT:
%   fid = pointer to the observation RINEX file
%
% OUTPUT:
%   time = observation GPS time (continuous)
%   datee = date (year,month,day,hour,minute,second)
%   num_sat = number of available satellites (NOTE: RINEX v3.xx does not output 'sat' and 'sat_types')
%   sat = list of all visible satellites
%   sat_types = ordered list of satellite types ('G' = GPS, 'R' = GLONASS, 'S' = SBAS)
%   tow = observation GPS time (seconds-of-week)
%
% DESCRIPTION:
%   Scan the first line of each epoch (RINEX) and return
%   the information it contains.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini.
%
% Portions of code contributed by Damiano Triglione (2012).
% Portions of code contributed by Andrea Gatti (2013).
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
time = NaN;
sat = [];
sat_types = [];
num_sat = 0;
datee=[NaN NaN NaN NaN NaN NaN]; %Preallocation not useful (see last line of code)
eof = 0;
tow = NaN;
if (nargout > 3)
    datee_RequestedInOutputFlag = true;
else
    datee_RequestedInOutputFlag = false;
end% if

%search data
while (eof==0)
    %read the string
    lin = fgets(fid);
    %answer = strfind(lin,'COMMENT');
    keywords = {'COMMENT', 'MARKER NAME', 'MARKER NUMBER', 'APPROX POSITION XYZ', 'ANTENNA: DELTA H/E/N'};
    answer = [];
    s = 1;
    while (s <= length(keywords) && isempty(answer))
        answer = strfind(lin,keywords{s});
        s = s + 1;
    end
    %if it is a line that should be skipped read the following one
    while (~isempty(answer) && ~feof(fid))
        lin = fgetl(fid);
        %check again
        answer = [];
        s = 1;
        while (s <= length(keywords) && isempty(answer))
            answer = strfind(lin,keywords{s});
            s = s + 1;
        end
    end
    %check if the end of file is reached
    if (feof(fid) == 1);
        return
    end

    %check RINEX version
    if (~strcmp(lin(1),'>')) %RINEX v2.xx
        
        %check if it is a string that should be analyzed
        if (strcmp(lin(29),'0') || strcmp(lin(29),'1') || strcmp(lin(29),'2'))
            
            %save time information
            data   = textscan(lin(1:26),'%f%f%f%f%f%f');
            year   = data{1};
            month  = data{2};
            day    = data{3};
            hour   = data{4};
            minute = data{5};
            second = data{6};
            
            %computation of the GPS time in weeks and seconds of week
            year = four_digit_year(year);
            [week, tow] = date2gps([year, month, day, hour, minute, second]);
            [time] = weektow2time(week, tow, 'G');
            
            %number of visible satellites
            [num_sat] = sscanf(lin(30:32),'%d');
            
            %keep just the satellite data
            lin = ExtractSubstring(lin, 33, 68);
            
            %remove 'blank spaces' and unwanted characters at the end of the string
            lin = RemoveUnwantedTrailingSpaces(lin);
            
            %read additional lines, depending on the number of satellites
            nlines = ceil(num_sat/12);
            for n = 1 : nlines - 1
                lin = [lin ExtractSubstring(fgetl(fid), 33, 68)];
                lin = RemoveUnwantedTrailingSpaces(lin);
            end
            
            pos = 1;
            sat = zeros(num_sat,1);
            sat_types = char(32*uint8(ones(num_sat,1))');
            for i = 1 : num_sat
                %check if GPS satellites are labeled 'G' or not labeled
                if (strcmp(lin(pos),' '))
                    type = 'G';
                else
                    type = lin(pos);
                end
                % sat_types = [sat_types; type];
                sat_types(i) = type;
                % sat(i) = sscanf(lin(pos+1:pos+2),'%d');
                sat(i) = mod((lin(pos+1)-48)*10+(lin(pos+2)-48),160);
                pos = pos + 3;
            end
            
            eof = 1;
        end
        
    else %RINEX v3.xx
        
        %check if it is a string that should be analyzed
        if (strcmp(lin(32),'0') || strcmp(lin(32),'1') || strcmp(lin(32),'2'))
            
            %save time information
            data   = textscan(lin(2:29),'%f%f%f%f%f%f');
            year   = data{1};
            month  = data{2};
            day    = data{3};
            hour   = data{4};
            minute = data{5};
            second = data{6};
            
            %computation of the GPS time in weeks and seconds of week
            [week, tow] = date2gps([year, month, day, hour, minute, second]);
            [time] = weektow2time(week, tow, 'G');
            
            %number of visible satellites
            [num_sat] = sscanf(lin(33:35),'%d');
            
            eof = 1;
        end
    end
end

if datee_RequestedInOutputFlag
    datee = [year month day hour minute second];
end %if
