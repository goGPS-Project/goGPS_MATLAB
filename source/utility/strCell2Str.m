function str = strCell2Str(str_cell, separator)
% SYNTAX:
%   [str] = strCell2Str(str_cell, <separator == ' '>);
%
% INPUT:
%   str_cell  = cell array of strings
%   separator = <optional> contains the separator string (white space as default)
%
% OUTPUT:
%   str = single string of the element list separated by the "separator"
%
% DESCRIPTION:
%   Convert to a character array the input cell of strings

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
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
    if nargin == 1
        separator = ' ';
    end
    if ischar(str_cell)
        str_cell = {str_cell};
    end
    if ~isempty(str_cell)
        str = str_cell{1};
        for i = 2 : numel(str_cell)
            str = sprintf('%s%s%s', str, separator, str_cell{i});
        end
    else
        str = '';
    end
end
