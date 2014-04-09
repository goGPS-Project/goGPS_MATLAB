function varargout = plot_google_map(varargin)
% function h = plot_google_map(varargin)
% Plots a google map on the current axes using the Google Static Maps API
%
% USAGE:
% h = plot_google_map(Property, Value,...)
% Plots the map on the given axes. Used also if no output is specified
%
% Or:
% [lonVec latVec imag] = plot_google_map(Property, Value,...)
% Returns the map without plotting it
%
% PROPERTIES:
%    Height (640)   - Height of the image in pixels (max 640)
%    Width  (640)    - Width of the image in pixels (max 640)
%    Scale (2)          - (1/2) Resolution scale factor . using Scale=2 will
%                                  double the resulotion of the downloaded image (up
%                                  to 1280x1280) and will result in finer rendering,
%                                  but processing time will be longer.
%    MapType         - ('roadmap')  Type of map to return. Any of [roadmap, 
%                                  satellite, terrain, hybrid) See the Google Maps API for
%                                  more information. 
%    Alpha (1)         - (0-1) Transparency level of the map (0 is fully
%                                  transparent). While the map is always
%                                  moved to the bottom of the plot (i.e. will
%                                  not hide previously drawn items), this can
%                                  be useful in order to increase readability
%                                  if many colors are ploted (using SCATTER
%                                  for example).
%    Marker            - The marker argument is a text string with fields
%                                  conforming to the Google Maps API. The
%                                  following are valid examples:
%                                  '43.0738740,-70.713993' (dflt midsize orange marker)
%                                  '43.0738740,-70.713993,blue' (midsize blue marker)
%                                  '43.0738740,-70.713993,yellowa' (midsize yellow
%                                  marker with label "A")
%                                  '43.0738740,-70.713993,tinyredb' (tiny red marker
%                                  with label "B")
%    Refresh (1)    - (0/1) defines whether to automatically refresh the
%                                  map upon zoom action on the figure.
%    AutoAxis (1) - (0/1) defines whether to automatically adjust the axis
%                                  of the plot to avoid the map being stretched.
%                                  This will adjust the span to be correct
%                                  only if the figure window is square. If
%                                  the figure window is rectangular, it will
%                                  still appear somewhat stretched.
%    APIKey            - (string) set your own API key which you obtained from Google: 
%                                   http://developers.google.com/maps/documentation/staticmaps/#api_key
%                                   This will enable up to 25,000 map
%                                   requests per day, compared to 400
%                                   requests without a key. 
%                                   To set the key, use:
%                                   plot_google_map('APIKey','SomeLongStringObtaindFromGoogle')
%                                   To disable the use of a key, use:
%                                   plot_google_map('APIKey','')
%
% OUTPUT:
%    h                         -  Handle to the plotted map
%
%    lonVect         -  Vector of Longidute coordinates (WGS84) of the image 
%    latVect         -  Vector of Latidute coordinates (WGS84) of the image 
%    imag                 -  Image matrix (height,width,3) of the map
%
% EXAMPLE - plot a map showing some capitals in Europe:
%    lat = [48.8708   51.5188   41.9260   40.4312   52.523   37.982];
%    lon = [2.4131    -0.1300    12.4951   -3.6788    13.415   23.715];
%    plot(lon,lat,'.r','MarkerSize',20)
%    plot_google_map
%
% References:
% http://www.mathworks.com/matlabcentral/fileexchange/24113
% http://www.maptiler.org/google-maps-coordinates-tile-bounds-projection/
% http://developers.google.com/maps/documentation/staticmaps/
%
%  Acknowledgement to Val Schmidt for his submission of get_google_map.m
%
%  Author:
%  Zohar Bar-Yehuda
% Version 1.2 - 16/06/2012
%       - Support use of the "scale=2" parameter by default for finer rendering (set scale=1 if too slow).
%       - Auto-adjust axis extent so the map isn't stretched.
%       - Set and use an API key which enables a much higher usage volume per day.
%  Version 1.1 - 25/08/2011
%
% Adapted for goGPS by Eugenio Realini (2013)

% store parameters in global variable (used for auto-refresh)
global inputParams
persistent apiKey
if isnumeric(apiKey)
    % first run, check if API key file exists
    if exist('api_key.mat','file')
        load api_key
    else
        apiKey = '';
    end
end
axHandle = gca;
inputParams.(['ax' num2str(axHandle*1e6,'%.0f')]) = varargin;

% Handle input arguments

height = 640;
width = 640;
scale = 2;
maptype = 'roadmap';
alphaData = 1;
autoRferesh = 1;
autoAxis = 1;
hold on

markeridx = 1;
markerlist = {};
if nargin >= 2
    for idx = 1:2:length(varargin)
        switch lower(varargin{idx})
            case 'height'
                height = varargin{idx+1};
            case 'width'
                width = varargin{idx+1};
            case 'maptype'
                maptype = varargin{idx+1};
            case 'alpha'
                alphaData = varargin{idx+1};
            case 'refresh'
                autoRferesh = varargin{idx+1};
            case 'marker'
                markerlist{markeridx} = varargin{idx+1};
                markeridx = markeridx + 1;
            case 'autoaxis'
                autoAxis = varargin{idx+1};
            case 'apikey'
                apiKey = varargin{idx+1}; % set new key
                % save key to file
                funcFile = which('plot_google_map.m');
                pth = fileparts(funcFile);
                keyFile = fullfile(pth,'api_key.mat');
                save(keyFile,'apiKey')
            otherwise
                error(['Unrecognized variable: ' varargin{idx}])
        end
    end
end
if height > 640
    height = 640;
end
if width > 640
    width = 640;
end

curAxis = axis;
if isequal(curAxis,[0 1 0 1]) % probably an empty figure
    % display world map
    curAxis = [-180 180 -85 85];
    axis(curAxis)
end

if autoAxis
    % adjust current axis limit to avoid strectched maps
    [xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
    xExtent = diff(xExtent); % just the size of the span
    yExtent = diff(yExtent); 
    if xExtent > yExtent
        % enlarge the X extent
        centerY = mean(curAxis(3:4));
        spanY = (curAxis(4)-curAxis(3))/2;
        curAxis(3) = centerY-spanY*xExtent/yExtent;
        curAxis(4) = centerY+spanY*xExtent/yExtent;
    elseif yExtent > xExtent
        % enlarge the Y extent
        centerX = mean(curAxis(1:2));
        spanX = (curAxis(2)-curAxis(1))/2;
        curAxis(1) = centerX-spanX*yExtent/xExtent;
        curAxis(2) = centerX+spanX*yExtent/xExtent;
    end
    axis(curAxis) % update axis as quickly as possible, before downloading new image
    drawnow
end

% Delete previous map from plot (if exists)
if nargout <= 1 % only if in plotting mode
    curChildren = get(axHandle,'children');
    delete(findobj(curChildren,'tag','gmap'))
end

% Enforce Latitude constraints of EPSG:900913 
if curAxis(3) < -85
    curAxis(3) = -85;
end
if curAxis(4) > 85
    curAxis(4) = 85;
end

% Calculate zoom level for current axis limits
[xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
minResX = diff(xExtent) / width;
minResY = diff(yExtent) / height;
minRes = max([minResX minResY]);
tileSize = 256;
initialResolution = 2 * pi * 6378137 / tileSize; % 156543.03392804062 for tileSize 256 pixels
zoomlevel = floor(log2(initialResolution/minRes));

% Enforce valid zoom levels
if zoomlevel < 0 
    zoomlevel = 0;
end
if zoomlevel > 19 
    zoomlevel = 19;
end

% Calculate center coordinate in WGS1984
lat = (curAxis(3)+curAxis(4))/2;
lon = (curAxis(1)+curAxis(2))/2;

% CONSTRUCT QUERY URL
preamble = 'http://maps.googleapis.com/maps/api/staticmap';
location = ['?center=' num2str(lat,10) ',' num2str(lon,10)];
zoomStr = ['&zoom=' num2str(zoomlevel)];
sizeStr = ['&scale=' num2str(scale) '&size=' num2str(width) 'x' num2str(height)];
maptypeStr = ['&maptype=' maptype ];
if ~isempty(apiKey)
    keyStr = ['&key=' apiKey];
else
    keyStr = '';
end
markers = '&markers=';
for idx = 1:length(markerlist)
    if idx < length(markerlist)
        markers = [markers markerlist{idx} '%7C'];
    else
        markers = [markers markerlist{idx}];
    end
end
if ismember(maptype,{'satellite','hybrid'})
    filename = 'tmp.jpg';
    format = '&format=jpg';
    convertNeeded = 0;
else
    filename = 'tmp.png';
    format = '&format=png';
    convertNeeded = 1;
end
sensor = '&sensor=false';
url = [preamble location zoomStr sizeStr maptypeStr format markers sensor keyStr];

% Get the image
try
    java.net.InetAddress.getByName('maps.googleapis.com');
    try
        urlwrite(url,filename);
    catch % error downloading map
        fprintf('Unable to download map form Google Servers.\nPossible reasons: no network connection, or quota exceeded (25000 map requests per day).\n')
        return
    end
catch
    fprintf('Unable to connect to the Internet for downloading Google Maps background.\n')
end
[M Mcolor] = imread(filename);
M = cast(M,'double');
delete(filename); % delete temp file
width = size(M,2);
height = size(M,1);

% Calculate a meshgrid of pixel coordinates in EPSG:900913
centerPixelY = round(height/2);
centerPixelX = round(width/2);
[centerX,centerY] = latLonToMeters(lat, lon ); % center coordinates in EPSG:900913
curResolution = initialResolution / 2^zoomlevel/scale; % meters/pixel (EPSG:900913)
xVec = centerX + ((1:width)-centerPixelX) * curResolution; % x vector
yVec = centerY + ((height:-1:1)-centerPixelY) * curResolution; % y vector
[xMesh,yMesh] = meshgrid(xVec,yVec); % construct meshgrid 

% convert meshgrid to WGS1984
[lonMesh,latMesh] = metersToLatLon(xMesh,yMesh);

% We now want to convert the image from a colormap image with an uneven
% mesh grid, into an RGB truecolor image with a uniform grid.
% This would enable displaying it with IMAGE, instead of PCOLOR.
% Advantages are:
% 1) faster rendering
% 2) makes it possible to display together with other colormap annotations (PCOLOR, SCATTER etc.)

% Convert image from colormap type to RGB truecolor (if PNG is used)
if convertNeeded
    imag = zeros(height,width,3);
    for idx = 1:3
        imag(:,:,idx) = reshape(Mcolor(M(:)+1+(idx-1)*size(Mcolor,1)),height,width);
    end
else
    imag = M/255;
end

% Next, project the data into a uniform WGS1984 grid
sizeFactor = 1; % factoring of new image
uniHeight = round(height*sizeFactor);
uniWidth = round(width*sizeFactor);
latVect = linspace(latMesh(1,1),latMesh(end,1),uniHeight);
lonVect = linspace(lonMesh(1,1),lonMesh(1,end),uniWidth);
[uniLonMesh,uniLatMesh] = meshgrid(lonVect,latVect);
uniImag = zeros(uniHeight,uniWidth,3);

% old version (projection using INTERP2)
% for idx = 1:3
%      % 'nearest' method is the fastest. difference from other methods is neglible
%          uniImag(:,:,idx) =  interp2(lonMesh,latMesh,imag(:,:,idx),uniLonMesh,uniLatMesh,'nearest');
% end
uniImag =  myTurboInterp2(lonMesh,latMesh,imag,uniLonMesh,uniLatMesh);

if nargout <= 1 % plot map
    % display image
    h = image(lonVect,latVect,uniImag);
    set(gca,'YDir','Normal')
    set(h,'tag','gmap')
    set(h,'AlphaData',alphaData)
    
    % older version (display without conversion to uniform grid)
    % h =pcolor(lonMesh,latMesh,(M));
    % colormap(Mcolor)
    % caxis([0 255])
    % warning off % to avoid strange rendering warnings
    % shading flat
   
    uistack(h,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
    axis(curAxis) % restore original zoom
    if nargout == 1
        varargout{1} = h;
    end
    % if auto-refresh mode - override zoom callback to allow autumatic 
    % refresh of map upon zoom actions.
    zoomHandle = zoom; 
    panHandle = pan;
    if autoRferesh        
        set(zoomHandle,'ActionPostCallback',@update_google_map);        
        set(panHandle, 'ActionPostCallback',@update_google_map);
    else % disable zoom override
        set(zoomHandle,'ActionPostCallback',[]);
        set(panHandle, 'ActionPostCallback',[]);
    end
else % don't plot, only return map
    varargout{1} = lonVect;
    varargout{2} = latVect;
    varargout{3} = uniImag;
end


% Coordinate transformation functions

function [lon,lat] = metersToLatLon(x,y)
% Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
lon = (x ./ originShift) * 180;
lat = (y ./ originShift) * 180;
lat = 180 / pi * (2 * atan( exp( lat * pi / 180)) - pi / 2);

function [x,y] = latLonToMeters(lat, lon )
% Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
x = lon * originShift / 180;
y = log(tan((90 + lat) * pi / 360 )) / (pi / 180);
y = y * originShift / 180;


function ZI = myTurboInterp2(X,Y,Z,XI,YI)
% An extremely fast nearest neighbour 2D interpolation, assuming both input
% and output grids consist only of squares, meaning:
% - uniform X for each column
% - uniform Y for each row
XI = XI(1,:);
X = X(1,:);
YI = YI(:,1);
Y = Y(:,1);

xiPos = nan*ones(size(XI));
xLen = length(X);
yiPos = nan*ones(size(YI));
yLen = length(Y);
% find x conversion
xPos = 1;
for idx = 1:length(xiPos)
    if XI(idx) >= X(1) && XI(idx) <= X(end)
        while xPos < xLen && X(xPos+1)<XI(idx)
            xPos = xPos + 1;
        end
        diffs = abs(X(xPos:xPos+1)-XI(idx));
        if diffs(1) < diffs(2)
            xiPos(idx) = xPos;
        else
            xiPos(idx) = xPos + 1;
        end
    end
end
% find y conversion
yPos = 1;
for idx = 1:length(yiPos)
    if YI(idx) <= Y(1) && YI(idx) >= Y(end)
        while yPos < yLen && Y(yPos+1)>YI(idx)
            yPos = yPos + 1;
        end
        diffs = abs(Y(yPos:yPos+1)-YI(idx));
        if diffs(1) < diffs(2)
            yiPos(idx) = yPos;
        else
            yiPos(idx) = yPos + 1;
        end
    end
end
ZI = Z(yiPos,xiPos,:);


function update_google_map(obj,evd)
% callback function for auto-refresh
drawnow;
global inputParams
params = inputParams.(['ax' num2str(gca*1e6,'%.0f')]);
plot_google_map(params{:});
    