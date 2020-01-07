function [cur_line, Obs_types, pos_M, ifound_types, interval, sysId, antoff, antmod, marker] = RINEX_parse_hdr(buf, cur_line)

% SYNTAX:
%   [Obs_types, pos_M, ifound_types, interval, sysId, antoff, antmod, marker] = RINEX_parse_hdr(file);
%
% INPUT:
%   file = pointer to RINEX observation file
%
% OUTPUT:
%   Obs_types = cell of strings containing observation types
%               RINEX v2.xx --> e.g. L1C1P1...
%               RINEX v3.xx --> e.g. C1CL1CD1C...
%   pos_M = master station approximate position
%   ifound_types = boolean variable to check the correct acquisition of basic information
%   interval = observation interval in seconds
%   sysId = cell-array containing one-letter identifiers for constellations
%   antoff = antenna offset [m]
%   antmod = antenna model [string]
%   marker = marker name [string]
%
% DESCRIPTION:
%   RINEX observation file header analysis.

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
%  Contributors:     Damiano Triglione 2012, ...
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
ifound_types = 0;
Obs_types = cell(0,0);
sysId = cell(0,0);
pos_M = [];
interval = 0;
antmod = '';

%parse first line
cur_line = cur_line + 1; line = buf{cur_line};

%constellation counter for RINEX v3.xx
c = 1;

%check if the end of the header or the end of the file has been reached
while isempty(strfind(line,'END OF HEADER')) && ischar(line)
    %NOTE: ischar is better than checking if line is the number -1.

    answer = strfind(line,'# / TYPES OF OBSERV'); %RINEX v2.xx
    if ~isempty(answer)
        Obs_types{1} = [];
        nObs = sscanf(line(1:6),'%d');
        nLinObs = ceil(nObs/9);
        for i = 1 : nLinObs
            if (i > 1)
                cur_line = cur_line + 1; line = buf{cur_line};
            end
            n = min(nObs,9);
            for k = 1 : n
                ot = sscanf(line(k*6+1:k*6+6),'%s');
                Obs_types{1} = [Obs_types{1} ot];
            end
            nObs = nObs - 9;
        end

        ifound_types = 1;
    end

    answer = strfind(line,'SYS / # / OBS TYPES'); %RINEX v3.xx
    if ~isempty(answer)
        sysId{c} = sscanf(line(1),'%s');
        nObs = sscanf(line(2:6),'%d');
        Obs_types.(sysId{c}) = [];
        nLinObs = ceil(nObs/13);
        for i = 1 : nLinObs
            if (i > 1)
                cur_line = cur_line + 1; line = buf{cur_line};
            end
            n = min(nObs,13);
            for k = 0 : n-1
                ot = sscanf(line(6+k*4+1:6+k*4+4),'%s');
                Obs_types.(sysId{c}) = [Obs_types.(sysId{c}) ot];
            end
            nObs = nObs - 13;
        end

        c = c + 1;
        ifound_types = 1;
    end

    answer = strfind(line,'APPROX POSITION XYZ');
    if ~isempty(answer)
        X = sscanf(line(1:14),'%f');
        Y = sscanf(line(15:28),'%f');
        Z = sscanf(line(29:42),'%f');
        pos_M = [X; Y; Z];
    end

    answer = strfind(line,'ANTENNA: DELTA H/E/N');
    if ~isempty(answer)
        dU = sscanf(line(1:14),'%f');
        dE = sscanf(line(15:28),'%f');
        dN = sscanf(line(29:42),'%f');
        antoff = [dE; dN; dU];
    end

    answer = strfind(line,'ANT # / TYPE');
    if ~isempty(answer)
        antmod = sscanf(line(21:35),'%c');
        radtyp = sscanf(line(36:40),'%c');
        if (isempty(find(radtyp ~= ' ', 1)))
            radtyp = ' NONE';
        end
        antmod = [antmod radtyp];
    end

    answer = strfind(line,'INTERVAL');
    if ~isempty(answer)
        interval = sscanf(line(1:10),'%f');
    end

    answer = strfind(line,'MARKER NAME');
    if ~isempty(answer)
        marker = sscanf(line(1:4),'%s');
    end

    %parse next line
    cur_line = cur_line + 1; line = buf{cur_line};
end

% %apply the antenna offset from the marker (if available)
% if (any(pos_M) && any(antoff))
%     pos_M = local2globalPos(antoff, pos_M);
% end
