function [universal_path, is_valid] = checkPath(path)

% SYNTAX:
%   universal_path = checkPath(path)
%
% INPUT:
%   path
%
% OUTPUT:
%   universal_path
%   < is_valid >            optional, contains the status of existence
%
% DESCRIPTION:
%   Conversion of path OS specific to a universal one: "/" or "\" are converted to filesep 
%   if the second parameter is present is_valid contains the status of existence 
%       2 => is a file
%       7 => is a folder

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti
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

if not(isempty(path))
    if (iscell(path))
        % for each line of the cell+
        universal_path = cell(size(path));
        for c = 1 : length(path)
            universal_path{c} = regexprep(path{c}, '(\\(?![ ]))|(\/)', filesep);
            universal_path{c} = regexprep(universal_path{c}, ['\' filesep '\' filesep], filesep);
        end
        if (nargout == 2)
            is_valid = zeros(size(path));
            for c = 1 : length(path)
                is_valid(c) = exist(universal_path{c}, 'file'); % if it is a file is_valid contains 2, if it is a dir it contains 7
            end
        end
    else
        universal_path = regexprep(path, '(\\(?![ ]))|(\/)', filesep);
        universal_path = regexprep(universal_path, ['\' filesep '\' filesep], filesep);
        if (nargout == 2)
            is_valid = exist(path, 'file'); % if it is a file is_valid contains 2, if it is a dir it contains 7
        end
    end    
else
    universal_path = [];
    is_valid = 0;
end
