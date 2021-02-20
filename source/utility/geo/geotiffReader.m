function [img, georef, lat, lon, info] = geotiffReader(file_name)
% Custom simplyfied implementatation of geotiffread
% useful when the mapping toolbox is not available
%
% INPUT
%   file_name   full path of the geotiff
%
% OUTPUT
%   img         matrix containing the tiff image
%   georef      GeographicPostingsReference as struct
%                  - LatitudeLimits
%                  - LongitudeLimits
%                  - RasterSize
%                  - RasterInterpretation
%                  - ColumnsStartFrom
%                  - RowsStartFrom
%                  - SampleSpacingInLatitude
%                  - SampleSpacingInLongitude
%                  - RasterExtentInLatitude
%                  - RasterExtentInLongitude
%                  - XIntrinsicLimits
%                  - YIntrinsicLimits
%                  - CoordinateSystemType
%                  - AngleUnit
%   lat         array of latitude coordinates of the image [deg]
%   lon         array of longitude coordinates of the image [deg]
%   info        as output of infinfo
%
%
% SYNTAX
%   [img, georef, lat, lon, info] = geotiffReader(file_name)
    
%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

    img = imread(file_name);
    info = imfinfo(file_name);

    lat = flipud(info.ModelTiepointTag(5) - ((0 : info.Height - 1)' .* info.ModelPixelScaleTag(2)));
    lon = info.ModelTiepointTag(4) + ((0 : info.Width - 1)' .* info.ModelPixelScaleTag(1));

    georef = struct('LatitudeLimits',           [lat(1) lat(end)], ...
        'LongitudeLimits',          [lon(1) lon(end)], ...
        'RasterSize',               [info.Height info.Width], ...
        'RasterInterpretation',     'postings', ...
        'ColumnsStartFrom',         'north', ...
        'RowsStartFrom',            'west', ...
        'SampleSpacingInLatitude',  info.ModelPixelScaleTag(2), ...
        'SampleSpacingInLongitude', info.ModelPixelScaleTag(1), ...
        'RasterExtentInLatitude',   diff([lat(1) lat(end)]), ...
        'RasterExtentInLongitude',  diff([lon(1) lon(end)]), ...
        'XIntrinsicLimits',         [1 info.Width], ...
        'YIntrinsicLimits',         [1 info.Height], ...
        'CoordinateSystemType',     'geographic', ...
        'AngleUnit',                'degree'); % in geotiffread this is an object "GeographicPostingsReference"

    if isfield(info, 'GDAL_NODATA')
        img(img == str2num(info.GDAL_NODATA)) = nan;
    end
end
