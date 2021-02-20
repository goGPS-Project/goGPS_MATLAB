function versionChanger(new_version_str, base_dir)
% SYNTAX:
%    versionChanger(new_version_str, base_dir);
% EXAMPLE:
%    versionChanger('0.5.0');
%
% DESCRIPTION:
%    Change the version number in all the goGPS source files with standard
% header - it requires a unix system
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
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

% find all the m files in goGPS directory
if (nargin <= 1)
    base_dir = '.';
end

% Cleaning spaces
fprintf('Cleaning spaces\n');
sourceCleaner(base_dir);

fprintf('\nChanging version\n');
if (nargin == 0)
    new_version_str = Core.GO_GPS_VERSION;
end
[~, list] = dos(['find ' base_dir ' -name \*.m']);
list = textscan(list,'%s','Delimiter','\n','whitespace','');
list = list{1};

tic
for i = 1 : length(list)
    file_name = list{i};
    fprintf('Opening file %3d/%3d: %s', i, length(list), file_name);
    fid = fopen(file_name, 'r');
    txt = char(fread(fid))';
    fclose(fid);
    te = regexp(txt, '(?<=(\n%\s*\|___\/\s*v\s*))(.[^\n]*)','once','tokenExtents');
    if not(isempty(te))
        txt = [txt(1:te(1)) new_version_str txt(te(2)+1:end)];
        fid = fopen(file_name, 'w');
        fwrite(fid, txt);
        fclose(fid);
        fprintf(' -> changed\n');
    else
        fprintf('\n');
    end
end
toc;
