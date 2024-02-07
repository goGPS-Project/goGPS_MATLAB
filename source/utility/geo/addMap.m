% Adds a background map to a given axis using various map providers.
%
% DESCRIPTION:
% This function fetches map tiles from the specified provider and displays them as a background on the given axis. It offers flexibility in specifying various parameters such as latitude and longitude limits, zoom level, and map provider.
%
% INPUT:
%   'demo'      If this parameter is passed with the value 'demo', the function enters a demonstration mode.
%
%   'ax'        Handle to the target axes. Default is the current axes (gca).
%
%   'lat_lim'   Vector specifying the latitude limits [min_lat, max_lat]. Default is ylim of the current axes.
%   'lon_lim'   Vector specifying the longitude limits [min_lon, max_lon]. Default is xlim of the current axes.
%
%   'zoom_lev'  Scalar specifying the zoom level for the map. If not provided, it's calculated automatically.
%
%   'alpha'     Alpha value of the map [0..1]
%
%   'img_only'  to get only the image without plotting it (an ax should be present)
%
%   'provider': String specifying the name of the map provider. Options include:
%       - Tested:
%         - OpenStreetMap: Free editable map of the world.
%         - ArcGIS: Platform for maps, applications, and data.
%         - OpenTopoMap: Open-source topographic map.
%         - GoogleRoad: Google's road map.
%         - GoogleSatellite: Google's satellite imagery.
%         - GoogleTerrain: Google's terrain map with physical features.
%         - GoogleHybrid: Combination of Google's satellite and road map.
%         
% OUTPUT:
% - 'im_h': Handle to the background map image.
% - 'img_bg': map_image
% - 'latitudes': latitudes
% - 'longitudes': longitudes
%
% SYNTAX:
%   [im_h, longitudes, latitudes, img_bg] = addMap(Name, Value, ...);
%   Where Name-Value pairs can be any of the INPUT parameters specified above.
%
% EXAMPLE:
%   im_h = addMap('ax', gca, 'provider', 'OpenStreetMap', 'lat_lim', [45.5, 46], 'lon_lim', [8.8, 9.3]);
%
% NOTE:
% Providers:
% - Tested:
%   - OpenStreetMap: https://www.openstreetmap.org/copyright
%   - ArcGIS: https://www.esri.com/en-us/legal/copyright-trademarks
%   - OpenTopoMap: https://opentopomap.org/about#verwendung
%   - GoogleRoad: https://cloud.google.com/maps-platform/terms
%   - GoogleSatellite: https://cloud.google.com/maps-platform/terms
%   - GoogleTerrain: https://cloud.google.com/maps-platform/terms
%   - GoogleHybrid: https://cloud.google.com/maps-platform/terms
    
%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
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

function [im_h, longitudes, latitudes, img_bg] = addMap(varargin)    

    % Check for 'demo' call
    if nargin == 1 && strcmpi(varargin{1}, 'demo')
        demoMapTiles();
        return;
    end

    % Default values
    ax = [];
    provider = 'ArcGIS';
    zoom_lev = [];
    img_only = false;
    alpha = [];
    lat_lim = [];
    lon_lim = [];
    flag_m_map = false;

    % Parse varargin
    k = 1;
    while k < length(varargin)
        switch lower(varargin{k})
            case 'img_only'
                img_only = varargin{k + 1};
            case 'ax'
                ax = varargin{k + 1};
            case {'provider', 'map', 'map_type'}
                provider = varargin{k + 1};
            case 'lat_lim'
                lat_lim = varargin{k + 1};
            case 'lon_lim'
                lon_lim = varargin{k + 1};
            case 'zoom_lev'
                zoom_lev = varargin{k + 1};
            case 'alpha'
                alpha = varargin{k + 1};
            case 'm_map'
                flag_m_map = varargin{k + 1};
            otherwise
                warning('Unknown parameter: %s', varargin{k});
        end
        k = k + 2;        
    end

    % Detect if m_map is in use
    if img_only
        using_m_map = false;
    else 
        if isempty(ax)
            ax = gca;
        end
        if flag_m_map
            if isempty(ax.Tag)
                ax.Tag = 'm_map';
            end
        end
        using_m_map = isUsingMMap(ax);        
        if isempty(lat_lim) || isempty(lon_lim)
            lat_lim = ylim(ax);
            if isempty(ax)
                ax = gca;
            end
            lon_lim = xlim(ax);
            if using_m_map
                [x_lim, y_lim] = m_xy2ll(lon_lim, lat_lim);
                [lon_lim, lat_lim] = m_xy2ll([x_lim(1), x_lim(2)], [y_lim(1), y_lim(2)]);
            end
        end
    end

    lon_lim = sort(lon_lim);
    lat_lim = sort(lat_lim);

    lon_diff = lon_lim(2) - lon_lim(1);
    lat_diff = lat_lim(2) - lat_lim(1);
    if ~img_only
        if isempty(ax)
            ax = gca;
        end
        if isstruct(ax.UserData)
            ax.UserData.provider =  provider;
            ax.UserData.lat_lim =  lat_lim;
            ax.UserData.lon_lim =  lon_lim;
        else
            ax.UserData = struct('provider', provider, ...
                'lat_lim', lat_lim, ...
                'lon_lim', lon_lim);
        end
    end
    
    buffer = 0.1 * max(lon_diff, lat_diff); % Buffer value in degrees
        
    % Calculate the appropriate zoom level if not provided
    if nargin < 5 || isempty(zoom_lev)
        zoom_lon = -log2(lon_diff/360);
        zoom_lat = -log2(lat_diff/180);
        zoom_lev = ceil(min(min(zoom_lon+2, zoom_lat+2), max(zoom_lon, zoom_lat)));
    end    

    lat_lim = sort(lat_lim);
    lon_lim = sort(lon_lim);

    % Extend the lat/lon limits with the buffer
    lat_lim_buffered = [lat_lim(1) - buffer, lat_lim(2) + buffer];
    lon_lim_buffered = [lon_lim(1) - buffer, lon_lim(2) + buffer];
    
    % Ensure the buffered limits are within the valid range
    lat_lim_buffered(1) = max(lat_lim_buffered(1), -85);
    lat_lim_buffered(2) = min(lat_lim_buffered(2), 85);
    lon_lim_buffered(1) = max(lon_lim_buffered(1), -180);
    lon_lim_buffered(2) = min(lon_lim_buffered(2), 180);

    % Download the tiles with the buffered lat/lon limits
    try
        [img_bg, latitudes, longitudes] = downloadMapTiles(lat_lim_buffered, lon_lim_buffered, buffer, zoom_lev, provider);
    catch
        Core.getLogger.addError('No map downloaded, check internet connection!');
        return
    end
    
    if ~img_only
        if isempty(ax)
            ax = gca;
        end

        % Display the image using latitudes and longitudes
        cur_child = get(ax,'children');
        if using_m_map
            bgmap = findobj(cur_child,'Tag','m_image');
            for i = numel(bgmap):-1:1
                if ~isfield(bgmap(i).UserData, 'bgmap')
                    bgmap(i) = [];
                end
            end
        else
            bgmap = findobj(cur_child,'tag','bgmap');
        end
        if ~isempty(bgmap)
            if using_m_map
                [x,~] = m_ll2xy(longitudes, mean(latitudes)*ones(size(longitudes)));
                [~,y] = m_ll2xy(mean(longitudes)*ones(size(latitudes)), latitudes);
                bgmap(1).XData = x(~isnan(x));
                bgmap(1).YData = y(~isnan(y));
                bgmap(1).CData = img_bg(~isnan(y),~isnan(x),:);
                if ~isempty(alpha)
                    bgmap(1).AlphaData = alpha;
                else
                    bgmap(1).AlphaData = mean(bgmap(1).AlphaData(:), 'omitnan');
                end
            else
                bgmap(1).XData = longitudes;
                bgmap(1).YData = (latitudes);
                bgmap(1).CData = img_bg;
                if ~isempty(alpha)
                    bgmap(1).AlphaData = alpha;
                end
            end
        else
            if using_m_map
                if ~isempty(alpha)
                    im_h = m_image(longitudes, latitudes, img_bg);
                    set(im_h, 'AlphaData', alpha, 'UserData', struct('bgmap', true));
                else
                    im_h = m_image(longitudes, latitudes, img_bg, 'UserData', struct('bgmap', true));
                end                
                hold on; im_hpan = image([-pi pi], [-pi pi], ones(2,2), 'AlphaData', 0); % add fake layer to allow pannig
            else
                if ~isempty(alpha)
                    im_h = image(ax, 'XData', longitudes, 'YData', (latitudes(:)), 'CData', img_bg, 'Parent', ax, 'tag', 'bgmap', 'AlphaData', alpha);
                else
                    im_h = image(ax, 'XData', longitudes, 'YData', (latitudes(:)), 'CData', img_bg, 'Parent', ax, 'tag', 'bgmap');
                end
                hold on; im_hpan = image([-180 180], [-90 90], ones(2,2), 'AlphaData', 0); % add fake layer to allow pannig
            end
            im_hpan.UserData.fakelayer = 'panning';
            uistack(im_h,'bottom')
            uistack(im_hpan,'bottom')
            if using_m_map
                uistack(im_h)
            end
        end
        if ~using_m_map
            axis(ax, 'on');
            axis equal;
            set(ax, 'xlim', lon_lim, 'ylim', lat_lim, 'YDir', 'normal');
        end

        h_zoom = zoom(ax.Parent); % Get the zoom object of the figure
        h_pan = pan(ax.Parent);   % Get the pan object of the figure

        % Set the callbacks for zoom and pan actions
        set(h_zoom, 'ActionPostCallback', @(src, event) handleAxisChanges(ax));
        set(h_pan, 'ActionPostCallback', @(src, event) handleAxisChanges(ax));
    end
    if img_only
        im_h = img_bg;
    end
end

function handleAxisChanges(ax)
    % Handles changes in the axis limits or zoom level.
    %
    % Syntax:
    %   handleAxisChanges(ax, lon_diff, lat_diff)
    %
    % Input:
    %   ax - Handle to the target axes.
    %
    % Output:
    %   None.
    %
    % Example:
    %   This function is generally used as a callback and not called directly.

    persistent prev_provider prev_zoom_level prev_lat_lim prev_lon_lim

    provider = ax.UserData.provider;
    same_provider = ~isempty(prev_provider) && strcmp(prev_provider, provider);

    lat_lim = get(ax, 'YLim');
    lon_lim = get(ax, 'XLim');
    
    % Calculate the current zoom level
    lon_diff = lon_lim(2) - lon_lim(1);
    lat_diff = lat_lim(2) - lat_lim(1);
    zoom_lon = -log2(lon_diff/360);
    zoom_lat = -log2(lat_diff/180);
    current_zoom_level = ceil(max(zoom_lon, zoom_lat));

    using_m_map = isUsingMMap(ax);

    % Check if the change in zoom level is significant or if panning occurred
    if same_provider && ~isempty(prev_zoom_level) && ~isempty(prev_lat_lim) && ~isempty(prev_lon_lim)
        % Check for panning:
        % Compute buffer
        if using_m_map
            % Using m_map, the image background have basically no buffer
            buffer_x = 1e-8 * m_xy2ll(lon_diff, 0);
            buffer_y = 1e-8 * m_xy2ll(0,lat_diff);
            d_lon = prev_lon_lim - m_xy2ll(lon_lim, 0);
            d_lat = prev_lat_lim - m_xy2ll(0, lat_lim);
        else
            buffer_x = 0.01 * lon_diff;
            buffer_y = 0.01 * lat_diff;
            d_lon = prev_lon_lim - lon_lim;
            d_lat = prev_lat_lim - lat_lim;            
        end

        if current_zoom_level == prev_zoom_level && max(abs(d_lon)) < buffer_x && max(abs(d_lat)) < buffer_y
            if using_m_map
                %[lon_lim, lat_lim] = resetMMapProj(ax);
            end
            return;
        end
    end

    prev_provider = provider;
    prev_zoom_level = current_zoom_level;
    if using_m_map
        prev_lon_lim = m_xy2ll(lon_lim, 0);
        prev_lat_lim = m_xy2ll(0, lat_lim);
    else
        prev_lon_lim = lon_lim;
        prev_lat_lim = lat_lim;
    end
    
    % Fetch the provider from the UserData field of the axis
    
    % Reload the tiles
    try
        using_m_map = isUsingMMap(ax);    
        if using_m_map
            % Since we are here let's fix auto
            [lon_lim, lat_lim] = resetMMapProj(ax);
        end
        addMap('ax', ax, 'provider', provider, 'lat_lim', lat_lim, 'lon_lim', lon_lim);
    catch ex
        Core_Utils.printEx(ex);  
    end
end

function [lon_lim, lat_lim] = resetMMapProj(ax, lon_lim, lat_lim)
    global MAP_PROJECTION %#ok<GVMIS>
    if nargin < 3
        x = xlim();
        y = ylim();
        [lon_lim, lat_lim] = m_xy2ll(x, y);
    end
    lon_lim = sort(lon_lim);
    lat_lim = sort(lat_lim);
    axis normal
    % proj_info = m_proj('get');
    current_proj_type = MAP_PROJECTION.name;
    m_proj(current_proj_type, 'lon', lon_lim, 'lat', lat_lim);
    % delete previous grids
    is_fancy = false;
    for element = ax.Children(:)'
        if strfind(element.Tag, 'm_grid') == 1
            if contains(element.Tag, 'fancybox')
                is_fancy = true;
            end
            delete(element);
        end
    end
    m_grid('box', iif(is_fancy, 'fancy', 'on')); % Add grid lines
end


function [img, latitudes, longitudes] = downloadMapTiles(lat_lim, lon_lim, buffer, zoom_lev, provider)
    % Limit zoom level
    zoom_lev = min(zoom_lev, 19);
    
    % Limit latitude and longitude
    lat_lim(1) = max(lat_lim(1), -85);
    lat_lim(2) = min(lat_lim(2), 85);
    lon_lim(1) = max(lon_lim(1), -179.9999);
    lon_lim(2) = min(lon_lim(2), 179.9999);

    % Compute tile limits
    tile_x_min = floor((lon_lim(1) + 180) / 360 * (2^zoom_lev));
    tile_x_max = floor((lon_lim(2) + 180) / 360 * (2^zoom_lev));
    tile_y_min = floor((1 - log(tan(lat_lim(2) * pi / 180) + 1 / cos(lat_lim(2) * pi / 180)) / pi) / 2 * (2^zoom_lev));
    tile_y_max = floor((1 - log(tan(lat_lim(1) * pi / 180) + 1 / cos(lat_lim(1) * pi / 180)) / pi) / 2 * (2^zoom_lev));

    % PRecache an empty image to contain the tiles
    single_tile = imread(sprintf('https://tile.openstreetmap.org/%d/%d/%d.png', zoom_lev, tile_x_min, tile_y_min));
    img = zeros(size(single_tile, 1) * (tile_y_max - tile_y_min + 1), size(single_tile, 2) * (tile_x_max - tile_x_min + 1), 3, 'uint8');

    % Download the data
    options = weboptions('CertificateFilename', '');
    for tx = tile_x_min:tile_x_max
        for ty = tile_y_min:tile_y_max
            switch provider
                case 'OpenStreetMap'
                    url = sprintf('https://tile.openstreetmap.org/%d/%d/%d.png', zoom_lev, tx, ty);
                case 'ArcGIS'
                    url = sprintf('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/%d/%d/%d', zoom_lev, ty, tx);
                case 'OpenTopoMap'
                    url = sprintf('https://tile.opentopomap.org/%d/%d/%d.png', zoom_lev, tx, ty);
                % need keys integrated with the software:
                case {'GoogleRoad', 'road'}
                    map_type = 'roadmap';
                case {'GoogleSatellite', 'satellite'}
                    map_type = 'satellite';
                case {'GoogleTerrain', 'terrain'}
                    map_type = 'terrain';
                case {'GoogleHybrid', 'hybrid'}
                    map_type = 'hybrid';
                otherwise
                    error('Provider not supported');
            end

            if ismember(provider, {'GoogleRoad', 'GoogleSatellite', 'GoogleTerrain', 'GoogleHybrid', 'road', 'satellite', 'terrain', 'hybrid'})
                [img, latitudes, longitudes] = getGoogleMapImage(lat_lim, lon_lim, buffer, map_type, options);
                if ~any(img(:))
                    % Fallback on ArcGIS in case the google image fails to be download (missing API key?)
                    Core.getLogger.addWarning('Google map failed to be download, check your API key, switching to ArcGIS map');
                    [img, latitudes, longitudes] = downloadMapTiles(lat_lim, lon_lim, buffer, zoom_lev, 'ArcGIS');
                    return;
                else
                    return;
                end
            end

            try
                [tile_img, map] = webread(url, options);
            catch
                map = [];
                try
                    [tile_img] = webread(url, options);                
                catch
                    tile_img  = zeros(1,1,3);
                end
            end

            if ~isempty(map)
                % If the image has a colormap, convert it to RGB
                tile_img = uint8(ind2rgb(tile_img, map) * 255);
            end

            
            % merge the tile
            x_offset = (tx - tile_x_min) * size(single_tile, 2);
            y_offset = (ty - tile_y_min) * size(single_tile, 1);
            img(y_offset+1:y_offset+size(tile_img, 1), x_offset+1:x_offset+size(tile_img, 2), :) = tile_img;
        end
    end

    % Get tile latitudes
    mercator_latitudes = zeros(size(img,1),1);
    for ty = tile_y_min:tile_y_max
        % For latitudes:
        % merge the tile
        y_offset = (ty - tile_y_min) * size(single_tile, 1);
        mercator_latitudes(y_offset+1:y_offset+size(tile_img, 1)) = getMercatoreTileLatitudes(ty, zoom_lev, size(tile_img, 1));
    end
            
    % Generate new equally spaced latitudes for stretching
    latitudes = flipud(linspace(min(mercator_latitudes), max(mercator_latitudes), 2*size(img, 1))');
    % Initialize the output image
    img_stretched = zeros([2*size(img,1) size(img,2)], 'uint8');
    % Stretch each channel of the image
    if any(img(:))
        for ch = 1:size(img, 3)
            img_stretched(:,:,ch) = interp1(mercator_latitudes, double(img(:,:,ch)), latitudes, 'linear');
        end
    end
    img = img_stretched;
    
    longitudes = tile2lon([tile_x_min tile_x_max+1], zoom_lev);
    longitudes = linspace(longitudes(1), longitudes(end), size(img_stretched,2));
end

function lon = tile2lon(tx, zoom_lev)
    % convert tx to latitude
    lon = tx / 2^zoom_lev * 360.0 - 180.0;
end

function latitudes = getMercatoreTileLatitudes(ty, zoom_lev, tile_size)
    % Given a Mercatore tile number and the zoom level
    % Return the latitude of each line of the tile
    if nargin < 3
        tile_size = 256; % default tile size
    end

    y = 0:(tile_size-1);
    lat_rad = atan(sinh(pi * (1 - 2 * (ty * tile_size + y) / (2^zoom_lev * tile_size))));
    latitudes = rad2deg(lat_rad);
end

function latitudes = getTileLatitudes(top_lat, bottom_lat, numRows)
    % Compute the Mercator Y values corresponding to the given latitudes
    top_y = log(tan((pi / 4) + (top_lat * (pi / 180) / 2)));
    bottom_y = log(tan((pi / 4) + (bottom_lat * (pi / 180) / 2)));

    % Generate the Mercator Y values for each row within the tile
    y_values = linspace(top_y, bottom_y, numRows);

    % Convert the Mercator Y values back to latitudes
    latitudes = (2 * atan(exp(y_values)) - (pi / 2)) * (180 / pi);
end

function [img, latitudes, longitudes] = getGoogleMapImage(lat_lim, lon_lim, buffer, map_type, options)
    
    zoom_lon = -log2((diff(lon_lim)-2*buffer)/360);
    zoom_lat = -log2((diff(lat_lim)-2*buffer)/180);
    zoom_lev = ceil(min(zoom_lon, zoom_lat));
    zoom_lev = min(zoom_lev, 19);

    % Calculate center coordinate
    lat = (lat_lim(1) + lat_lim(2)) / 2;
    lon = (lon_lim(1) + lon_lim(2)) / 2;

    % Default Google settings
    width = 640;  % maximum width for free tier
    height = 640; % maximum height for free tier
    scale = 2;    % default scale
    rrm = Remote_Resource_Manager.getInstance;
    api_key = rrm.getGoogleMapsAPI;  % fetch the API key
    if isempty(api_key)
        img = nan(2,2,3);
        latitudes = [];
        longitudes = [];
        return
    end

    % Compute the latitudes and longitudes for each pixel
    tile_size = 256; % default tile size
    initial_resolution = 2 * pi * 6378137 / tile_size; % for a tile size of 256x256
    %zoom = ceil(log2(initialResolution/minRes));
    %zoom = max(1,min(zoom, 19));

    % Construct query URL for Google Maps
    url = sprintf('http://maps.googleapis.com/maps/api/staticmap?center=%f,%f&zoom=%d&scale=%d&size=%dx%d&maptype=%s&format=png&key=%s', ...
        lat, lon, zoom_lev, scale, width, height, map_type, api_key);
    try
        [img , map] = webread(url, options);
        if ~isempty(map)
            % If the image has a colormap, convert it to RGB
            img = uint8(ind2rgb(img, map) * 255);
        end
    catch
        img = nan(2, 2, 3);  % return a blank image on error
    end

     
    % Calculate a meshgrid of pixel coordinates in EPSG:900913
    width = size(img,2);
    height = size(img,1);
    center_pixel_y = round(height/2);
    center_pixel_x = round(width/2);
    [center_x, center_y] = latLonToMeters(lat, lon); % center coordinates in EPSG:900913
    cur_resolution = initial_resolution / 2^zoom_lev / scale; % meters/pixel (EPSG:900913)
    x_vec = center_x + ((1:width)-center_pixel_x) * cur_resolution; % x vector
    y_vec = center_y + ((height:-1:1)-center_pixel_y) * cur_resolution; % y vector
    
    % convert meshgrid to WGS1984
    [lon_vec,lat_vec] = metersToLatLon(x_vec,y_vec);

    new_lat_vec = linspace(min(lat_vec), max(lat_vec), size(img,1));
    img_stretched = zeros(length(new_lat_vec), size(img,2), size(img,3), 'uint8');
    for ch = 1:size(img,3)
        img_stretched(:,:,ch) = interp1(lat_vec, double(img(:,:,ch)), new_lat_vec, 'linear');
    end
    img = img_stretched;

    latitudes = new_lat_vec; % Since it's uniformly spaced in latitude for each column
    longitudes = lon_vec(:)'; % Since it's uniformly spaced in longitude for each row
end

% Coordinate transformation functions for google maps ()
% References:
%  http://www.mathworks.com/matlabcentral/fileexchange/24113
%  http://www.maptiler.org/google-maps-coordinates-tile-bounds-projection/
%  http://developers.google.com/maps/documentation/staticmaps/
%  https://github.com/zoharby/plot_google_map
%
% Acknowledgements:
%  Val Schmidt for his submission of get_google_map.m
%
% Author:
%  Zohar Bar-Yehuda

function [lon,lat] = metersToLatLon(x,y)
    % Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum
    origin_shift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
    lon = (x ./ origin_shift) * 180;
    lat = (y ./ origin_shift) * 180;
    lat = 180 / pi * (2 * atan( exp( lat * pi / 180)) - pi / 2);
end

function using_m_map = isUsingMMap(ax)
    % Search for m_map handles to understand if m_map is in use 
    try
        using_m_map = strfind(ax.Tag, 'm_') == 1;
        if isempty(using_m_map)
            using_m_map = false;
        end
    catch
        using_m_map = false;
    end
end

function [x,y] = latLonToMeters(lat, lon )
    % Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"
    origin_shift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
    x = lon * origin_shift / 180;
    y = log(tan((90 + lat) * pi / 360 )) / (pi / 180);
    y = y * origin_shift / 180;
end


% ====================================================================================
%  DEMO
% ====================================================================================

function demoMapTiles()
    % Initial map settings
    lat_lim = [45.6, 46.3]; % e.g., latitude of lake Como
    lon_lim = [8.8, 9.5];   % e.g., longitude of lake Como

    flag_m_map = true;

    % Create a figure
    f = figure('Name', 'Map Tiles Viewer', 'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);

    % UI controls for latitude and longitude limits
    uicontrol('Style', 'text', 'String', 'Lat Min:', 'Position', [10, 570, 60, 20]);
    lat_min_edit = uicontrol('Style', 'edit', 'String', num2str(lat_lim(1)), 'Position', [70, 570, 60, 20]);
    
    uicontrol('Style', 'text', 'String', 'Lat Max:', 'Position', [10, 540, 60, 20]);
    lat_max_edit = uicontrol('Style', 'edit', 'String', num2str(lat_lim(2)), 'Position', [70, 540, 60, 20]);

    uicontrol('Style', 'text', 'String', 'Lon Min:', 'Position', [10, 510, 60, 20]);
    lon_min_edit = uicontrol('Style', 'edit', 'String', num2str(lon_lim(1)), 'Position', [70, 510, 60, 20]);
    
    uicontrol('Style', 'text', 'String', 'Lon Max:', 'Position', [10, 480, 60, 20]);
    lon_max_edit = uicontrol('Style', 'edit', 'String', num2str(lon_lim(2)), 'Position', [70, 480, 60, 20]);

    % List of map providers that don't require an API key
    providers = {'GoogleSatellite', 'GoogleRoad', 'GoogleTerrain', ...
                 'OpenStreetMap', 'ArcGIS', 'OpenTopoMap', 'CartoDB', 'CartoDBDark', ...
                 'ESRITopo', 'ESRIShadedRelief', 'OSMDE', ...
                 'USGSTopo', 'USGSImagery', 'USGSImageryTopo', 'USGSShadedReliefOnly', 'USGSHydro', 'USGSHistoricalTopo'};
    
    % Dropdown menu for map providers
    uicontrol('Style', 'text', 'String', 'Map Provider:', 'Position', [10, 400, 80, 20]);
    mapProviderMenu = uicontrol('Style', 'popupmenu', 'String', providers, 'Position', [100, 400, 150, 25], ...
                                'Callback', @(src, event) updateTiles(flag_m_map));
    
    % Refresh button to update the tiles
    uicontrol('Style', 'pushbutton', 'String', 'Refresh', 'Position', [10, 360, 120, 30], 'Callback', @(src, event) updateTiles(flag_m_map));

    ax = axes('Parent', f, 'Position', [0.3, 0.2, 0.65, 0.65]);

    if flag_m_map
        ylim(sort(lat_lim));
        xlim(sort(lon_lim));
        m_proj('mercator', 'lon', lon_lim, 'lat', lat_lim); % Initialize m_map with mercator projection
        m_grid('box', 'fancy'); % Add grid lines
    end

    % Function to update the tiles
    function updateTiles(flag_m_map)
        lon_lim = [str2double(get(lon_min_edit, 'String')), str2double(get(lon_max_edit, 'String'))];
        lat_lim = [str2double(get(lat_min_edit, 'String')), str2double(get(lat_max_edit, 'String'))];
        
        if flag_m_map
            [lon_lim, lat_lim] = resetMMapProj(ax, lon_lim, lat_lim);
        end
        selectedProvider = providers{get(mapProviderMenu, 'Value')};
        addMap('ax', ax, 'provider', selectedProvider, 'lat_lim', lat_lim, 'lon_lim', lon_lim, 'alpha', 0.8);
    end

    % Initialize with the default selected map provider
    updateTiles(flag_m_map);
end
