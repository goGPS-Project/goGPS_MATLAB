function dtm_bulk_load(dtm_dir_path, tile_size)

% SYNTAX:
%   dtm_bulk_load(dtm_dir_path, tile_size);
%
% INPUT:
%   dtm_dir_path = path to ASCII-GRID files
%   tile_size = size of the tile side (number of cells, default=100)
%
% DESCRIPTION:
%   This function loads all the ASCII-GRID dtm files from the directory
%   passed in input, saves them in .MAT binary format and splits them in tiles.

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

%computation time evaluation
tic

%----------------------------------------------------------------------------------------------
% LOAD ASCII-GRID FILES AND SAVE THEM IN .MAT BINARY FORMAT
%----------------------------------------------------------------------------------------------

%directory containing dtm files
dtm_dir = dir(dtm_dir_path);

%check the number of files contained in the directory
nmax = size(dtm_dir,1);

j = 0;
for i = 1 : nmax

    %read the name of the i-th file
    got = getfield(dtm_dir,{i,1},'name');

    %get the number of characters in the filename
    fn_length = size(got,2);

    %check that the filename has the ".hdr" extension
    if (fn_length >= 4 & strcmp(got(fn_length - 3 : fn_length), '.hdr'))

        j = j+1;

        dtm_target = strcat(dtm_dir_path, '/', got);

        %read the header information
        [ncols(j), nrows(j), cellsize(j), xll(j), yll(j), nodata(j)] = dtm_hdr_load(dtm_target);

        %check consistency between cellsize and nodata values
        if (j > 1 & cellsize(j) ~= cellsize(j - 1))
            error('... ERROR: cellsize value not consistent among dtm files');
        end

        if (j > 1 & nodata(j) ~= nodata(j - 1))
            error('... ERROR: nodata value not consistent among dtm files');
        end

        target_length = size(dtm_target,2);

        %load the file containing height values
        dtm_cell{j} = load(strcat(dtm_target(1,1:target_length - 4),'.dtm'));

        %computation of the border coordinates for the loaded file
        W_border(j) = xll(j);
        N_border(j) = yll(j) + (nrows(j)-1)*cellsize(j);
        S_border(j) = yll(j);
        E_border(j) = xll(j) + (ncols(j)-1)*cellsize(j);
    end
end

if (j == 0)
    error(['... ERROR: no proper DTM file found in folder ' dtm_dir_path '.']);
end

%resolution and nodata value must be the same for all the dtm files, thus
% they are saved in simple variables
cellsize = cellsize(1);
nodata = nodata(1);

%----------------------------------------------------------------------------------------------
% TILING
%----------------------------------------------------------------------------------------------

%if the tile size is not specified, use default value
if (nargin < 2)
    tile_size = 100;
end

%if 'tiles' directory does not exist, create it
if (~isdir([dtm_dir_path '/tiles']))
    mkdir(dtm_dir_path,'tiles');
end

%if 'tiles' directory is not empty, delete its content
if (~isempty(dir([dtm_dir_path '/tiles'])))
    delete([dtm_dir_path '/tiles/*.*']);
end

%computation of the border coordinates of the whole dataset
extents = [min(W_border), max(N_border), min(S_border), max(E_border)];

%computation of the number of tiles (number of rows and columns)
tile_col_num = ceil((extents(4) - extents(1))/(tile_size*cellsize));
tile_row_num = ceil((extents(2) - extents(3))/(tile_size*cellsize));

tile_header = struct('ncols',tile_size,'nrows',tile_size,'cellsize',cellsize,'nodata',nodata); %#ok<NASGU>
tile_georef = [];

%bounding box initialization (W N S E ground coordinates) for the first tile (upper left)
tile_extents = [extents(1), extents(2), extents(2)-(tile_size-1)*cellsize, extents(1)+(tile_size-1)*cellsize];

%cycle to process all the tiles
for i = 1 : tile_row_num
    for j = 1 : tile_col_num
        fprintf('\n');
        fprintf('Tile generation (row:%d/%d column:%d/%d)...',i,tile_row_num,j,tile_col_num);
        fprintf('\n');

        %computation of the lower left corner coordinates
        tile_georef(i,j,1) = tile_extents(1) - cellsize/2;
        tile_georef(i,j,2) = tile_extents(2) + cellsize/2;
        tile_georef(i,j,3) = tile_extents(3) - cellsize/2;
        tile_georef(i,j,4) = tile_extents(4) + cellsize/2;

        %cycle to process all the cells inside a tile
        for tile_row = 1 : tile_size
            for tile_col = 1 : tile_size

                %coordinates of the central point of the cell
                cell_E = tile_extents(1) + (tile_col-1)*cellsize;
                cell_N = tile_extents(2) - (tile_row-1)*cellsize;

                %if the cell is contained within the bounding box of a portion of the DTM...
                found = find ( (W_border <= cell_E) & (N_border > cell_N) & (S_border <= cell_N) & (E_border > cell_E) );
                if (~isempty(found))
                    dtm_row = round((N_border(found) - cell_N)/cellsize + 1);
                    dtm_col = round((cell_E - W_border(found))/cellsize + 1);
                    % ...save the value...
                    tile(tile_row,tile_col) = dtm_cell{found}(dtm_row, dtm_col);
                else
                    % ...otherwise set the nodata value
                    tile(tile_row,tile_col) = nodata;
                end
            end
        end

        %save the tile in a .mat file only if it contains at least one non-null cell
        if (~isempty(find(tile ~= nodata, 1)))
            save(strcat(dtm_dir_path, '/tiles/tile_', num2str(i), '_', num2str(j)), 'tile', '-V6');
        end

        %compute the western and eastern borders of the bounding box of the next tile (new column)
        tile_extents(1) = tile_extents(4)+cellsize;
        tile_extents(4) = tile_extents(1)+(tile_size-1)*cellsize;
    end

    %compute the northern and southern borders of the bounding box of the next tile (new row)
    tile_extents(2) = tile_extents(3)-cellsize;
    tile_extents(3) = tile_extents(2)-(tile_size-1)*cellsize;

    %reset the western and eastern borders of the bounding box of the next tile to the initial values (first column)
    tile_extents(1) = extents(1);
    tile_extents(4) = extents(1)+(tile_size-1)*cellsize;
end

%save the files used to manage tiles during the execution of goGPS
save(strcat(dtm_dir_path, '/tiles/tile_header'), 'tile_header', '-V6');
save(strcat(dtm_dir_path, '/tiles/tile_georef'), 'tile_georef', '-V6');

%----------------------------------------------------------------------------------------------

%computation time evaluation
toc
