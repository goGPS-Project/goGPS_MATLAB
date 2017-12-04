function [ocean_load_disp] = load_BLQ(filename, marker)

% SYNTAX:
%   [ocean_load_disp] = load_BLQ(filename, marker);
%
% INPUT:
%   filename = ocean loading displacement file (.BLQ)
%   marker   = marker name(s)
%
% OUTPUT:
%   ocean_load_disp = ocean loading displacement values read from .BLQ file
%
% DESCRIPTION:
%   Reads ocean loading displacement values from a file in BLQ format.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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

ocean_load_disp=[];
for m = 1 : length(marker)
    ocean_load_disp(m).name = marker{m}; %#ok<*AGROW>
    ocean_load_disp(m).matrix = zeros(6,11);
    ocean_load_disp(m).available = 0;
end

found = 0;
for file_blq=1:size(filename,1)
    if (~isempty(filename))
        fid = fopen(char(filename(file_blq,:)),'r');
        if (fid ~= -1)
            if (file_blq == 1)
                fprintf('%s', ['Reading ocean loading file ', char(filename(file_blq,:)), '...']); fprintf('\n');
            end
            while (~feof(fid) && found < length(marker))
                line = fgetl(fid);
                for m = 1 : length(marker)
                    if (~strcmp(line(1:2),'$$'))
                        if (strfind(line,marker{m}))
                            line = fgetl(fid);
                            while(strcmp(line(1:2),'$$'))
                                line = fgetl(fid);
                            end
                            for l = 1 : 6
                                ocean_load_disp(m).matrix(l,:) = cell2mat(textscan(line, '  %f %f %f %f %f %f %f %f %f %f %f'));
                                line = fgetl(fid);
                            end
                            ocean_load_disp(m).available = 1;
                            found = found + 1;
                            break
                        end
                    end
                end
            end
            fclose(fid);
        else
            ocean_load_disp=[];
            fprintf('%s', ['... WARNING: Ocean loading file ', char(filename(file_blq,:)), ' could not be read.']); fprintf('\n');
        end
    else
        ocean_load_disp=[];
        fprintf('%s', ['... WARNING: Ocean loading file ', char(filename(file_blq,:)), ' could not be read.']); fprintf('\n');
    end
end

if (found == 0)
    ocean_load_disp=[];
    fprintf('%s', '... WARNING: Ocean loading parameters not found.'); fprintf('\n');
end
