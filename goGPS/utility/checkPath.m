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

%----------------------------------------------------------------------------------------------
%                           goGPS v0.5.9
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%   File contributed by Andrea Gatti
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

if not(isempty(path))
    if (iscell(path))
        % for each line of the cell+
        universal_path = cell(size(path));
        for c = 1 : length(path)
            universal_path{c} = regexprep(path{c}, '(\\(?=[a-zA-Z0-9]))|(\/)', filesep);
        end
        if (nargout == 2)
            is_valid = zeros(size(path));
            for c = 1 : length(path)
                is_valid(c) = exist(universal_path{c}, 'file'); % if it is a file is_valid contains 2, if it is a dir it contains 7
            end
        end
    else
        universal_path = regexprep(path, '(\\(?=[a-zA-Z0-9]))|(\/)', filesep);
        if (nargout == 2)
            is_valid = exist(path, 'file'); % if it is a file is_valid contains 2, if it is a dir it contains 7
        end
    end    
else
    universal_path = [];
    is_valid = 0;
end
