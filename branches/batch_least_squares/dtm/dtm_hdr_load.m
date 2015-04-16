function [ncols, nrows, cellsize, xll, yll, nodata] = dtm_hdr_load(hdr_target)

% SYNTAX:
%   [ncols, nrows, cellsize, xll, yll, nodata] = dtm_hdr_load(hdr_target);
%
% INPUT:
%   hdr_target = path to the .hdr file
%
% OUTPUT:
%   ncols = number of grid columns
%   nrows = number of grid rows
%   cellsize = ground size of a cell
%   xll = X coordinate of the center of the lower left cell
%   yll = Y coordinate of the center of the lower left cell
%   nodata = value used for cells not containing data
%
% DESCRIPTION:
%   Function that reads a .hdr file containing and ASCII-GRID header.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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


fn_length = size(hdr_target,2);

fprintf('\n');
fprintf('Loading %s...', hdr_target(1,1:fn_length - 4));
fprintf('\n');

%checking file existence
if (fopen(hdr_target,'r') == -1)
    error('... ERROR: header file (.hdr) not found.');
end

%header reading
[hdr_1, hdr_2] = textread(hdr_target,'%s %f');

%saving the number of lines that compose each header column
[hdr_rows_1, null] = size(hdr_1); %#ok<NASGU>
[hdr_rows_2, null] = size(hdr_2); %#ok<NASGU>

%if the header follows the ASCII-GRID format
if  (isnumeric(hdr_2) & (hdr_rows_1 == 5 | hdr_rows_1 == 6) & ...
    (hdr_rows_2 == 5 | hdr_rows_2 == 6) & hdr_rows_1 == hdr_rows_2)

    %save values
    flag_corner = 0;
    for i = 1 : hdr_rows_1
        variable = cell2mat(hdr_1(i));
        switch variable
            case {'ncols', 'NCOLS'}
                %number of columns
                ncols = hdr_2(i);
            case {'nrows', 'NROWS'}
                %number of lines
                nrows = hdr_2(i);
            case {'cellsize', 'CELLSIZE'}
                %cell dimension
                cellsize = hdr_2(i);
            case {'xllcenter', 'XLLCENTER'}
                %X coordinate of the center of the lower left cell
                xll = hdr_2(i);
            case {'yllcenter', 'YLLCENTER'}
                %Y coordinate of the center of the lower left cell
                yll = hdr_2(i);
            case {'xllcorner', 'XLLCORNER'}
                %X coordinate of the lower left corner of the lower left cell
                xll = hdr_2(i);
                flag_corner = 1;
            case {'yllcorner', 'YLLCORNER'}
                %Y coordinate of the lower left corner of the lower left cell
                yll = hdr_2(i);
                flag_corner = 1;
            case {'nodata_value', 'NODATA_VALUE', 'NODATA_value'}
                %value used for cells not containing data
                nodata = hdr_2(i);
            otherwise
                error('... ERROR: header is not standard ASCII-GRID.');
        end
    end
    %if cell corners are used
    if(flag_corner == 1)
        xll = xll + cellsize / 2;
        yll = yll + cellsize / 2;
    end
else
    %the header does not follow the standard ASCII-GRID format
    error('... ERROR: header is not standard ASCII-GRID.');
end
