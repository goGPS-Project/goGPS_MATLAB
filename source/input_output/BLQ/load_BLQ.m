function [ocean_load_disp, status] = load_BLQ(filename, marker)

% SYNTAX:
%   [ocean_load_disp, status] = load_BLQ(filename, marker);
%
% INPUT:
%   filename = ocean loading displacement file (.BLQ)
%   marker   = marker name(s)
%
% OUTPUT:
%   ocean_load_disp = ocean loading displacement values read from .BLQ file
%   status = data found status flag (<0 ocean loding not found =0 marker not found >0 found)
%
% DESCRIPTION:
%   Reads ocean loading displacement values from a file in BLQ format.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
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

% ocean_load_disp struct definition
% ocean_load_disp.marker : marker name
% ocean_load_disp.matrix : ocean loading displacement matrix
log = Core.getLogger();
ocean_load_disp=[];
for m = 1 : length(marker)
    ocean_load_disp(m).name = marker{m}; %#ok<*AGROW>
    ocean_load_disp(m).matrix = zeros(6,11);
    ocean_load_disp(m).available = 0;
end

status = 0;
for file_blq = 1 : size(filename,1)
    if (~isempty(filename))
        fid = fopen(char(filename(file_blq,:)),'rt');
        if (fid ~= -1)
            if (file_blq == 1)
                log.addMessage(log.indent(['Reading ocean loading file ', File_Name_Processor.getFileName(char(filename(file_blq,:))), '...']));
            end
            while (~feof(fid) && status < length(marker))
                line = fgetl(fid);
                for m = 1 : length(marker)
                    if ~isempty(strtrim(line))  && (~strcmp(line(1:2),'$$'))
                        if length(line) >= 6 && ~isempty((strfind(upper(line(1:6)), ['  ' upper(marker{m})])))
                            line = fgetl(fid);
                            while(strcmp(line(1:2),'$$'))
                                line = fgetl(fid);
                            end
                            for l = 1 : 6
                                ocean_load_disp(m).matrix(l,:) = cell2mat(textscan(line, '  %f %f %f %f %f %f %f %f %f %f %f'));
                                line = fgetl(fid);
                            end
                            ocean_load_disp(m).available = 1;
                            status = status + 1;
                            break
                        end
                    end
                end
            end
            fclose(fid);
        else
            ocean_load_disp = [];
            log.addWarning(['Ocean loading file ', char(filename(file_blq,:)), ' could not be read.']);
            status = -1;
        end
    else
        ocean_load_disp = [];
        log.addWarning('Ocean loading file not provided');
        status = -2;
    end
end

if (status == 0)
    ocean_load_disp = [];
end
