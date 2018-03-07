% prettyQuiver
%--------------------------------------------------------------------------
%
% plot with projection an area of the world
%
% POSSIBLE SINTAXES:
%   quiverHandler = prettyQuiver(dataQuiver, phi, lambda);
%
%   prettyQuiver(dataQuiver, phi, lambda, shape);
%   prettyQuiver(dataQuiver, phi, lambda, projection);
%   prettyQuiver(dataQuiver, phi, lambda, lineCol);
%
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax);
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, shape);
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection);
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, lineCol);
%
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, shape, projection);
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, shape, lineCol);
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection, shape);
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection, lineCol);
%
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection, shape, lineCol);
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, shape, projection, lineCol);
%
% EXAMPLE:
%   prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, 'Miller Cylindrical');
%
% INPUT:
%   dataQuiver      is a two column array [u v]
%   phi             array [degree]
%   lambda          array [degree]
%   phiMin          minimum latitude    [degree]
%   phiMax          maximum latitude    [degree]
%   lambdaMin       minimum longitude   [degree]
%   lambdaMax       maximum longitude   [degree]
%   projection      type of projection to be used "standard" is the default
%   shape           shapefile to load as coast (or country) contour
%                       - coast         only coasts coarse
%                       - 50m           1:50000000 scale country contours
%                       - 30m           1:30000000 scale country contours
%                       - 10m           1:10000000 scale country contours
%   lineCol         [1 1 1] array of RGB component to draw the contour lines
%
% DEFAULT VALUES:
%    projection = 'Lambert'
%
% AVAILABLE PROJECTION:
%    * Lambert
%      Stereographic
%      Orthographic
%      Azimuthal Equal-area
%      Azimuthal Equidistant
%      Gnomonic
%      Satellite
%      Albers Equal-Area Conic
%      Lambert Conformal Conic
%      Mercator
%    * Miller Cylindrical
%    * Equidistant Cylindrical (world dataQuiver)
%      Oblique Mercator
%      Transverse Mercator
%      Sinusoidal
%      Gall-Peters
%      Hammer-Aitoff
%      Mollweide
%      Robinson
%    * UTM
%
% SEE ALSO:
%   mapPlot, mapPlot3D, quiver
%
% REQUIREMENTS:
%   M_Map: http://www.eos.ubc.ca/~rich/map.html
%   shape files with contours
%
% VERSION: 2.1
%
% CREDITS:
%   http://www.eos.ubc.ca/~rich/map.html
%
%   Andrea Gatti
%   DIIAR - Politecnico di Milano
%   2013-12-19
%

function quiverHandler = prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection, shape, lineCol)

%shape = 'coast';
%shape = '10m';
%shape = '30m';
%shape = '50m';

% lineCol = [0 0 0];

if size(dataQuiver,1) == 2
    dataQuiver = dataQuiver';
end
phi = phi(:);
lambda = lambda(:);

limitsOk = false;

% Manage opening a new figure;
tohold = false;
if length(findall(0,'Type','figure'))>=1
    if ishold
        tohold = true;
    else
        figure;
    end
end

switch (nargin)
    case 3
        shape = 'coast';
        lineCol = [0 0 0];
        if (ischar(phi))
            projection = 'Miller Cylindrical';
            if (sum(strcmp(phi,[{'coast'},{'10m'},{'30m'},{'50m'}])))
                shape = phi;
                if (ischar(lambda))
                    projection = lambda;                              % prettyQuiver(dataQuiver, shape, projection);
                else
                    lineCol = lambda;                                 % prettyQuiver(dataQuiver, shape, lineCol);
                end
            else
                projection = phi;
                if (ischar(lambda))
                    shape = lambda;                                   % prettyQuiver(dataQuiver, projection, shape);
                else
                    lineCol = lambda;                                 % prettyQuiver(dataQuiver, projection, lineCol);
                end
            end
            phiMin = 90;
            phiMax = -90;
            lambdaMin = -180;
            lambdaMax = 180;
            
            deltaPhi = (phiMax-phiMin)/size(dataQuiver,1);
            deltaLambda = (lambdaMax-lambdaMin)/size(dataQuiver,2);
            
            phi    = (phiMin + deltaPhi/2 : deltaPhi : phiMax - deltaPhi/2)';
            lambda = (lambdaMin + deltaLambda/2 :  deltaLambda :  lambdaMax  - deltaLambda/2)';
        else                                                              % prettyQuiver(dataQuiver, phi, lambda);
            projection = 'Miller Cylindrical';
            phiMin = max(phi);
            phiMax = min(phi);
            lambdaMin = min(lambda);
            lambdaMax = max(lambda);
        end
    case 4
        shape = 'coast';
        lineCol = [0 0 0];
        projection = 'Miller Cylindrical';
        if (ischar(phiMin))
            if (sum(strcmp(phiMin,[{'coast'},{'10m'},{'30m'},{'50m'}])))  % prettyQuiver(dataQuiver, phi, lambda, shape);
                shape = phiMin;
            else                                                          % prettyQuiver(dataQuiver, phi, lambda, projection);
                projection = phiMin;
            end
        elseif (length(phiMin) == 3)                                      % prettyQuiver(dataQuiver, phi, lambda, lineCol);
            lineCol = phiMin;
        end
        
        phiMin = max(phi);
        phiMax = min(phi);
        lambdaMin = min(lambda);
        lambdaMax = max(lambda);
    case 5
        shape = 'coast';
        lineCol = [0 0 0];
        projection = 'Miller Cylindrical';
        if (ischar(phiMin))
            if (sum(strcmp(phiMin,[{'coast'},{'10m'},{'30m'},{'50m'}])))
                shape = phiMin;
                if (ischar(phiMax))
                    projection = phiMax;                                  % prettyQuiver(dataQuiver, phiMin, phiMax, shape, projection);
                else
                    lineCol = phiMax;                                     % prettyQuiver(dataQuiver, phiMin, phiMax, shape, lineCol);
                end
            else
                projection = phiMin;
                if (ischar(phiMax))
                    shape = phiMax;                                       % prettyQuiver(dataQuiver, phiMin, phiMax, projection, shape);
                else
                    lineCol = phiMax;                                     % prettyQuiver(dataQuiver, phiMin, phiMax, projection, lineCol);
                end
            end
            phiMin = max(phi);
            phiMax = min(phi);
            lambdaMin = min(lambda);
            lambdaMax = max(lambda);
        else                                                             %  prettyQuiver(dataQuiver, phiMin, phiMax, lambdaMin, lambdaMax);
            limitsOk = true;
            lambdaMin = phiMin;
            lambdaMax = phiMax;
            phiMin = phi;
            phiMax = lambda;
            
            if (phiMin < phiMax)
                tmp = phiMin;
                phiMin = phiMax;
                phiMax = tmp;
            end
            
            deltaPhi = (phiMax-phiMin)/size(dataQuiver,1);
            deltaLambda = (lambdaMax-lambdaMin)/size(dataQuiver,2);
            
            phi    = (phiMin + deltaPhi/2 : deltaPhi : phiMax - deltaPhi/2)';
            lambda = (lambdaMin + deltaLambda/2 :  deltaLambda :  lambdaMax  - deltaLambda/2)';
        end
    case 6
        shape = 'coast';
        lineCol = [0 0 0];
        limitsOk = true;
        projection = 'Lambert';
        if (ischar(lambdaMin))
            if (sum(strcmp(lambdaMin,[{'coast'},{'10m'},{'30m'},{'50m'}])))
                shape = lambdaMin;                                        % prettyQuiver(dataQuiver, phiMin, phiMax, lambdaMin, lambdaMax, shape);
            else
                projection = lambdaMin;
            end
        elseif (length(lambdaMin) == 3)
            lineCol = lambdaMin;                                          % prettyQuiver(dataQuiver, phiMin, phiMax, lambdaMin, lambdaMax, lineCol);
        end
        
        lambdaMax = phiMax;
        lambdaMin = phiMin;
        phiMin = phi;
        phiMax = lambda;
        
        if (phiMin < phiMax)
            tmp = phiMin;
            phiMin = phiMax;
            phiMax = tmp;
        end
        
        deltaPhi = (phiMax-phiMin)/size(dataQuiver,1);
        deltaLambda = (lambdaMax-lambdaMin)/size(dataQuiver,2);
        
        phi    = (phiMin + deltaPhi/2 : deltaPhi : phiMax - deltaPhi/2)';
        lambda = (lambdaMin + deltaLambda/2 :  deltaLambda :  lambdaMax  - deltaLambda/2)';
    case 7
        shape = 'coast';
        lineCol = [0 0 0];
        limitsOk = true;
        projection = 'Lambert';
        if (ischar(lambdaMin))
            if (sum(strcmp(lambdaMin,[{'coast'},{'10m'},{'30m'},{'50m'}])))
                shape = lambdaMin;
                if (ischar(lambdaMax))
                    projection = lambdaMax;                               % prettyQuiver(dataQuiver, phiMin, phiMax, lambdaMin, lambdaMax, shape, projection);
                else
                    lineCol = lambdaMax;                                  % prettyQuiver(dataQuiver, phiMin, phiMax, lambdaMin, lambdaMax, shape, lineCol);
                end
            else
                projection = lambdaMin;
                if (ischar(lambdaMax))
                    shape = lambdaMax;                                    % prettyQuiver(dataQuiver, phiMin, phiMax, lambdaMin, lambdaMax, projection, shape);
                else
                    lineCol = lambdaMax;                                  % prettyQuiver(dataQuiver, phiMin, phiMax, lambdaMin, lambdaMax, projection, lineCol);
                end
            end
            
            lambdaMin = phiMin;
            lambdaMax = phiMax;
            phiMin = phi;
            phiMax = lambda;
            
            if (phiMin < phiMax)
                tmp = phiMin;
                phiMin = phiMax;
                phiMax = tmp;
            end
            
            deltaPhi = (phiMax-phiMin)/size(dataQuiver,1);
            deltaLambda = (lambdaMax-lambdaMin)/size(dataQuiver,2);
            
            phi    = (phiMin + deltaPhi/2 : deltaPhi : phiMax - deltaPhi/2)';
            lambda = (lambdaMin + deltaLambda/2 :  deltaLambda :  lambdaMax  - deltaLambda/2)';
        else
            projection = 'lambert';                                       % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax);
        end
    case 8
        shape = 'coast';
        lineCol = [0 0 0];
        limitsOk = true;
        if (ischar(projection))
            if (sum(strcmp(projection,[{'coast'},{'10m'},{'30m'},{'50m'}])))
                shape = projection;                                       % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, shape);
                projection = 'Lambert';
            else
                % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection);
            end
        elseif (length(projection) == 3)
            lineCol = projection;                                         % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, lineCol);
            projection = 'Lambert';
        end
    case 9
        lineCol = [0 0 0];
        limitsOk = true;
        if (ischar(projection))
            if (sum(strcmp(projection,[{'coast'},{'10m'},{'30m'},{'50m'}])))
                tmp = shape;
                shape = projection;
                if (ischar(tmp))
                    projection = tmp;                                     % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, shape, projection);
                else
                    lineCol = tmp;                                        % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, shape, lineCol);
                    projection = 'UTM';
                end
            else
                if (ischar(shape))
                    % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection, shape);
                else
                    lineCol = shape;                                      % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection, lineCol);
                    shape = 'coast';
                end
            end
        end
    case 10                                                               % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, projection, shape, lineCol)
        limitsOk = true;
        if (sum(strcmp(projection,[{'coast'},{'10m'},{'30m'},{'50m'}])))   % prettyQuiver(dataQuiver, phi, lambda, phiMin, phiMax, lambdaMin, lambdaMax, shape, projection, lineCol)
            tmp = shape;
            shape = projection;
            projection = tmp;
        end
end

if (phiMin < phiMax)
    tmp = phiMin;
    phiMin = phiMax;
    phiMax = tmp;
end

lambdaTmp = sort(lambda);
[val, idMax] = max(diff(lambdaTmp));
if (sum(diff(lambdaTmp) == val) == 1) && val > 10
    lambdaTmp(1:idMax) = lambdaTmp(1:idMax)+360;
    if ~limitsOk
        lambdaMax = lambdaTmp(idMax);
        lambdaMin = lambdaTmp(idMax+1);
    end
end

if(lambdaMax<lambdaMin)
    lambdaMax = lambdaMax+360;
end

if (~limitsOk)
    lambdaMin = lambdaMin-3;
    lambdaMax = lambdaMax+3;
    phiMin = min(phiMin+3,90);
    phiMax = max(phiMax-3,-90);
end

% setup the projection
if ~tohold
    if (strcmpi(projection,{'lambert'}) && abs(phiMax==-phiMin))
        projection='Miller Cylindrical';
    end
    
    if (sum(strcmpi(projection,[{'lambert'},{'UTM'},{'Sinusoidal'},{'Transverse Mercator'},{'Oblique Mercator'},{'Miller Cylindrical'}])))
        m_proj(projection,'long',[lambdaMin lambdaMax],'lat',[phiMax phiMin]);
    else
        m_proj(projection);
    end
    
    % Printing projection
    fprintf('Using projection: %s\n', projection);
    
    if sum(diff(lambda)<-200)
        lambda(lambda<0)=lambda(lambda<0)+360;
    end
    % plot the dataQuiver
    m_pcolor([lambdaMin lambdaMax],[phiMin phiMax], nan(2));
    % set the light
    shading flat;
    
end

ids = (lambda<lambdaMin) & (lambda < 0);
lambda(ids) = lambda(ids)+360;
hold on;
[xlocal,ylocal] = m_ll2xy(lambda(:),phi(:));
if isfield(dataQuiver,'rgb') % I also have colors!!!
    s=1;
    for c = 1:length(dataQuiver)
        quiverHandler = quiver(xlocal(s:(s+length(dataQuiver(c).u)-1)),ylocal(s:(s+length(dataQuiver(c).u)-1)),dataQuiver(c).u,dataQuiver(c).v,0, 'Color', [dataQuiver(c).rgb(1) dataQuiver(c).rgb(2) dataQuiver(c).rgb(3)]); % <========================= QUIVER function is here
        s = s+length(dataQuiver(c).u);
    end
else
    quiverHandler = quiver(xlocal,ylocal,dataQuiver.u(:),dataQuiver.v(:),0); % <========================= QUIVER function is here
end

% read shapefile
if (~strcmp(shape,'coast'))
    if (strcmp(shape,'10m'))
        M=m_shaperead('countries_10m');
    elseif (strcmp(shape,'30m'))
        M=m_shaperead('countries_30m');
    else
        M=m_shaperead('countries_50m');
    end
    [xMin,yMin] = m_ll2xy(lambdaMin,phiMin);
    [xMax,yMax] = m_ll2xy(lambdaMax,phiMax);
    for k=1:length(M.ncst)
        lamC = M.ncst{k}(:,1);
        ids = lamC < lambdaMin;
        lamC(ids) = lamC(ids) + 360;
        phiC = M.ncst{k}(:,2);
        [x,y] = m_ll2xy(lamC,phiC);
        if sum(~isnan(x))>1
            x(find(abs(diff(x))>=abs(xMax-xMin)*0.90)+1) = nan; % Romove lines that occupy more than th 90% of the plot
            line(x,y,'color', lineCol);
        end
    end;
else
    m_coast('line','color', lineCol);
end

m_grid('box','fancy','tickdir','in');
colorbar;

if tohold
    hold on;
else
    hold off;
end
