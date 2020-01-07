function [datee, num_sat, sat, sat_types] = RINEX_get_epoch_old(fid)

% SYNTAX:
%   [datee, num_sat, sat, sat_types] = RINEX_get_epoch(fid);
%
% INPUT:
%   fid = pointer to the observation RINEX file
%
% OUTPUT:
%   datee = date (year,month,day,hour,minute,second)
%   num_sat = number of available satellites (NOTE: RINEX v3.xx does not output 'sat' and 'sat_types')
%   sat = list of all visible satellites
%   sat_types = ordered list of satellite types ('G' = GPS, 'R' = GLONASS, 'S' = SBAS)
%
% DESCRIPTION:
%   Scan the first line of each epoch (RINEX) and return
%   the information it contains.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     Andrea Gatti, ...
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

% variable initialization
sat = [];
sat_types = [];
num_sat = 0;
datee=[NaN NaN NaN NaN NaN NaN];
eof = 0;


% search data
while (eof==0)
    % read the string
    lin = fgets(fid);
    % answer = strfind(lin,'COMMENT');
    keywords = {'COMMENT', 'MARKER NAME', 'MARKER NUMBER', 'APPROX POSITION XYZ', 'ANTENNA: DELTA H/E/N'};
    answer = [];
    s = 1;
    while (s <= length(keywords) && isempty(answer))
        answer = strfind(lin,keywords{s});
        s = s + 1;
    end
    % if it is a line that should be skipped read the following one
    while (~isempty(answer) && ~feof(fid))
        lin = fgetl(fid);
        % check again
        answer = [];
        s = 1;
        while (s <= length(keywords) && isempty(answer))
            answer = strfind(lin,keywords{s});
            s = s + 1;
        end
    end
    % check if the end of file is reached
    if (feof(fid) == 1)
        return
    end

    % check RINEX version
    if strcmp(lin(1),' ') %RINEX v2.xx

        % check if it is a string that should be analyzed
        if (strcmp(lin(28:30),' 0 ') || strcmp(lin(28:30),' 1 ') || strcmp(lin(28:30),' 2 '))

            % save time information
            datee = sscanf(lin(1:min(26,end)),'%f%f%f%f%f%f')';

            % year format from 2 to 4 digits
            [datee(1)] = four_digit_year(datee(1));

            % number of visible satellites
            [num_sat] = sscanf(lin(30:32),'%d');

            % keep just the satellite data
            lin = ExtractSubstring(lin, 33, 68);

            % remove 'blank spaces' and unwanted characters at the end of the string
            lin = deblank(lin);

            % read additional lines, depending on the number of satellites
            nlines = ceil(num_sat/12);
            for n = 1 : nlines - 1
                lin = deblank([lin ExtractSubstring(fgetl(fid), 33, 68)]);
            end

            pos = 1;
            sat = zeros(num_sat,1);
            sat_types = char(32*uint8(ones(num_sat,1))');
            for i = 1 : num_sat
                % check if GPS satellites are labeled 'G' or not labeled
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

    elseif strcmp(lin(1),'>') %RINEX v3.xx

        % check if it is a string that should be analyzed
        if (strcmp(lin(32),'0') || strcmp(lin(32),'1') || strcmp(lin(32),'2'))

            % save time information
            datee = sscanf(lin(3:min(30,end)),'%f%f%f%f%f%f');

            % number of visible satellites
            [num_sat] = sscanf(lin(33:35),'%d');

            eof = 1;
        end
    end
end
