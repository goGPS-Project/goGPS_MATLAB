function fileRep(base_dir, expression, replace, filter)
% SYNTAX:
%    fileRep( base_dir, expression, replace, <filter = '*.m'>)
% EXAMPLE:
%    fileRep('.','expression', 'replace phrase', '*.m');
%
% DESCRIPTION:
%    Replace an exprassion in all the '.m' files under base_dir - it requires a unix system
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

if nargin < 4 || isempty(filter)
    filter = '*.m';
end

% find all the m files in goGPS directory
if isempty(base_dir)
    base_dir = '.';
end

if ~exist(base_dir, 'dir')
    if iscell(base_dir)
        list = base_dir;
    else
        list = {base_dir};
    end
else
    [~, list] = dos(['find ' base_dir ' -name \' filter]);
    if ~isempty(list)
        list = textscan(list,'%s','Delimiter','\n','whitespace','');
        list = list{1};
    end
end
tic
for i = 1 : length(list)
    file_name = list{i};
    fid = fopen(file_name, 'r');
    txt = fread(fid,'*char')';
    fclose(fid);
    %occurencies = regexp(txt, expression, 'once');
    %clean_txt = regexprep(txt, expression, replace);
    occurencies = strfind(txt, expression);
    clean_txt = strrep(txt, expression, replace);
    
    if not(isempty(clean_txt)) && ~isempty(occurencies)
        fprintf('Opening file %3d/%3d: %s', i, length(list), file_name);
        fid = fopen(file_name, 'w');
        fwrite(fid, clean_txt);
        fclose(fid);
        fprintf(' -> changed\n');
    end
end
toc;
