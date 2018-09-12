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
    

    img = imread(file_name);
    info = imfinfo(file_name);

    lat = info.ModelTiepointTag(5) - ((0 : info.Height - 1)' .* info.ModelPixelScaleTag(2));
    lon = info.ModelTiepointTag(4) + ((0 : info.Width - 1)' .* info.ModelPixelScaleTag(1));
    
    georef = struct('LatitudeLimits',           [lat(end) lat(1)], ...
                    'LongitudeLimits',          [lon(1) lon(end)], ...
                    'RasterSize',               [info.Height info.Width], ...
                    'RasterInterpretation',     'postings', ...
                    'ColumnsStartFrom',         'north', ...
                    'RowsStartFrom',            'west', ...                    
                    'SampleSpacingInLatitude',  info.ModelPixelScaleTag(2), ...
                    'SampleSpacingInLongitude', info.ModelPixelScaleTag(1), ...
                    'RasterExtentInLatitude',   diff([lat(1) lat(end)]), ...
                    'RasterExtentInLongitude',  diff([lon(end) lon(1)]), ...
                    'XIntrinsicLimits',         [1 info.Width], ...
                    'YIntrinsicLimits',         [1 info.Height], ...
                    'CoordinateSystemType',     'geographic', ...
                    'AngleUnit',                'degree'); % in geotiffread this is an object "GeographicPostingsReference"
