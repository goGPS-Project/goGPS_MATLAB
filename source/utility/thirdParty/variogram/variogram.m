function S = variogram(x,y,varargin)

% isotropic and anisotropic experimental (semi-)variogram
%
% Syntax:
%   d = variogram(x,y)
%   d = variogram(x,y,'propertyname','propertyvalue',...)
%
% Description:
%   variogram calculates the experimental variogram in various 
%   dimensions. 
%
% Input:
%   x - array with coordinates. Each row is a location in a 
%       size(x,2)-dimensional space (e.g. [x y elevation])
%   y - column vector with values of the locations in x. 
%
% Propertyname/-value pairs:
%   nrbins - number bins the distance should be grouped into
%            (default = 20)
%   maxdist - maximum distance for variogram calculation
%            (default = maximum distance in the dataset / 2)
%   type -   'gamma' returns the variogram value (default)
%            'cloud1' returns the binned variogram cloud
%            'cloud2' returns the variogram cloud
%   plotit - true -> plot variogram
%            false -> don't plot (default)
%   subsample - number of randomly drawn points if large datasets are used.
%               scalar (positive integer, e.g. 3000)
%               inf (default) = no subsampling
%   anisotropy - false (default), true (works only in two dimensions)
%   thetastep - if anisotropy is set to true, specifying thetastep 
%            allows you the angle width (default 30°)
%   
%   
% Output:
%   d - structure array with distance and gamma - vector
%   
% Example: Generate a random field with periodic variation in x direction
% 
%     x = rand(1000,1)*4-2;  
%     y = rand(1000,1)*4-2;
%     z = 3*sin(x*15)+ randn(size(x));
%
%     subplot(2,2,1)
%     scatter(x,y,4,z,'filled'); box on;
%     ylabel('y'); xlabel('x')
%     title('data (coloring according to z-value)')
%     subplot(2,2,2)
%     hist(z,20)
%     ylabel('frequency'); xlabel('z')
%     title('histogram of z-values')
%     subplot(2,2,3)
%     d = variogram([x y],z,'plotit',true,'nrbins',50);
%     title('Isotropic variogram')
%     subplot(2,2,4)
%     d2 = variogram([x y],z,'plotit',true,'nrbins',50,'anisotropy',true);
%     title('Anisotropic variogram')
%
% Requirements:
%   The function uses parseargs (objectId=10670) 
%   by Malcolm wood as subfunction.
%
% See also: KRIGING, VARIOGRAMFIT
%
% Date: 21. July, 2011
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)


% error checking
if size(y,1) ~= size(x,1);
    error('x and y must have the same number of rows')
end

% check for nans
II = any(isnan(x),2) | isnan(y);
x(II,:) = [];
y(II)   = [];

% extent of dataset
minx = min(x,[],1);
maxx = max(x,[],1);
maxd = sqrt(sum((maxx-minx).^2));
nrdims = size(x,2);

% check input using PARSEARGS
params.nrbins      = 20;
params.maxdist     = maxd/2;
params.type        = {'default','gamma','cloud1','cloud2'};
params.plotit      = false;
params.anisotropy  = false;
params.thetastep   = 30;
params.subsample   = inf;
params = parseargs(params,varargin{:});

if params.maxdist > maxd;
    warning('Matlab:Variogram',...
            ['Maximum distance exceeds maximum distance \n' ... 
             'in the dataset. maxdist was decreased to ' num2str(maxd) ]);
    params.maxdist  = maxd;
end

if params.anisotropy && nrdims ~= 2 
    params.anisotropy = false;
    warning('Matlab:Variogram',...
            'Anistropy is only supported for 2D data');
end

% take only a subset of the data;
if ~isinf(params.subsample) && numel(y)>params.subsample;
    IX = randperm(numel(y));
    x  = x(IX(1:params.subsample),:);
    y  = y(IX(1:params.subsample),:);
end

% calculate bin tolerance
tol = params.maxdist/params.nrbins;

% calculate distance matrix
iid = distmat(x,params.maxdist);

% calculate squared difference between values of coordinate pairs
lam      = (y(iid(:,1))-y(iid(:,2))).^2;

% anisotropy
if params.anisotropy 
    nrthetaedges = floor(180/params.thetastep);
  
    % calculate with radians, not degrees
    params.thetastep = params.thetastep/180*pi;

    % calculate angles, note that angle is calculated clockwise from top
    theta    = atan2(x(iid(:,2),1)-x(iid(:,1),1),...
                     x(iid(:,2),2)-x(iid(:,1),2));
    
    % only the semicircle is necessary for the directions
    I        = theta < 0;
    theta(I) = theta(I)+pi;
    I        = theta >= pi-params.thetastep/2;
    theta(I) = 0;
        
    % create a vector with edges for binning of theta
    % directions go from 0 to 180 degrees;
    thetaedges = linspace(-params.thetastep/2,pi-params.thetastep/2,nrthetaedges);
    
    % bin theta
    [ntheta,ixtheta] = histc(theta,thetaedges);
    
    % bin centers
    thetacents = thetaedges(1:end)+params.thetastep/2;
    thetacents(end) = pi; %[];
end

% calculate variogram
switch params.type
    case {'default','gamma'}
        % variogram anonymous function
        fvar     = @(x) 1./(2*numel(x)) * sum(x);
        
        % distance bins
        edges      = linspace(0,params.maxdist,params.nrbins+1);
        edges(end) = inf;

        [nedge,ixedge] = histc(iid(:,3),edges);
        
        if params.anisotropy
            S.val      = accumarray([ixedge ixtheta],lam,...
                                 [numel(edges) numel(thetaedges)],fvar,nan);
            S.val(:,end)=S.val(:,1); 
            S.theta    = thetacents;
            S.num      = accumarray([ixedge ixtheta],ones(size(lam)),...
                                 [numel(edges) numel(thetaedges)],@sum,nan);
            S.num(:,end)=S.num(:,1);                 
        else
            S.val      = accumarray(ixedge,lam,[numel(edges) 1],fvar,nan);     
            S.num      = accumarray(ixedge,ones(size(lam)),[numel(edges) 1],@sum,nan);
        end
        S.distance = (edges(1:end-1)+tol/2)';
        S.val(end,:) = [];
        S.num(end,:) = [];

    case 'cloud1'
        edges      = linspace(0,params.maxdist,params.nrbins+1);
        edges(end) = inf;
        
        [nedge,ixedge] = histc(iid(:,3),edges);
        
        S.distance = edges(ixedge) + tol/2;
        S.distance = S.distance(:);
        S.val      = lam;  
        if params.anisotropy            
            S.theta   = thetacents(ixtheta);
        end
    case 'cloud2'
        S.distance = iid(:,3);
        S.val      = lam;
        if params.anisotropy            
            S.theta   = thetacents(ixtheta);
        end
end


% create plot if desired
if params.plotit
    switch params.type
        case {'default','gamma'}
            marker = 'o--';
        otherwise
            marker = '.';
    end
    
    if ~params.anisotropy
        plot(S.distance,S.val,marker);
        axis([0 params.maxdist 0 max(S.val)*1.1]);
        xlabel('h');
        ylabel('\gamma (h)');
        title('(Semi-)Variogram');
    else
        [Xi,Yi] = pol2cart(repmat(S.theta,numel(S.distance),1),repmat(S.distance,1,numel(S.theta)));
        surf(Xi,Yi,S.val)
        xlabel('h y-direction')
        ylabel('h x-direction')
        zlabel('\gamma (h)')
        title('directional variogram')
%         set(gca,'DataAspectRatio',[1 1 1/30])
    end
end
        
end


% subfunction distmat

function iid = distmat(X,dmax)

% constrained distance function
%
% iid -> [rows, columns, distance]
 

n     = size(X,1);
nrdim = size(X,2);
if size(X,1) < 1000;
    [i,j] = find(triu(true(n)));
    if nrdim == 1;
        d = abs(X(i)-X(j));
    elseif nrdim == 2;
        d = hypot(X(i,1)-X(j,1),X(i,2)-X(j,2));
    else
        d = sqrt(sum((X(i,:)-X(j,:)).^2));
    end
    I = d<=dmax;
    iid = [i(I) j(I) d(I)];
else
    ix = (1:n)';
    if nrdim == 1;
        iid = arrayfun(@distmatsub1d,(1:n)','UniformOutput',false);
    elseif nrdim == 2;
        % if needed change distmatsub to distmatsub2d which is numerically
        % better but slower
        iid = arrayfun(@distmatsub,(1:n)','UniformOutput',false);
    else
        iid = arrayfun(@distmatsub,(1:n)','UniformOutput',false);
    end
    nn  = cellfun(@(x) size(x,1),iid,'UniformOutput',true);  
    I   = nn>0;
    ix  = ix(I);
    nn  = nn(I);
    nncum = cumsum(nn);
    c     = zeros(nncum(end),1);
    c([1;nncum(1:end-1)+1]) = 1;
    i = ix(cumsum(c));
    iid = [i cell2mat(iid)];
    
end

function iid = distmatsub1d(i) 
    j  = (i+1:n)'; 
    d  = abs(X(i)-X(j));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end

function iid = distmatsub2d(i)  %#ok<DEFNU>
    j  = (i+1:n)'; 
    d = hypot(X(i,1) - X(j,1),X(i,2) - X(j,2));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end
    
function iid = distmatsub(i)
    j  = (i+1:n)'; 
    d = sqrt(sum(bsxfun(@minus,X(i,:),X(j,:)).^2,2));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end
end



% subfunction parseargs

function X = parseargs(X,varargin)

%PARSEARGS - Parses name-value pairs
%
% Behaves like setfield, but accepts multiple name-value pairs and provides
% some additional features:
% 1) If any field of X is an cell-array of strings, it can only be set to
%    one of those strings.  If no value is specified for that field, the
%    first string is selected.
% 2) Where the field is not empty, its data type cannot be changed
% 3) Where the field contains a scalar, its size cannot be changed.
%
% X = parseargs(X,name1,value1,name2,value2,...) 
%
% Intended for use as an argument parser for functions which multiple options.
% Example usage:
%
% function my_function(varargin)
%   X.StartValue = 0;
%   X.StopOnError = false;
%   X.SolverType = {'fixedstep','variablestep'};
%   X.OutputFile = 'out.txt';
%   X = parseargs(X,varargin{:});
%
% Then call (e.g.):
%
% my_function('OutputFile','out2.txt','SolverType','variablestep');

% The various #ok comments below are to stop MLint complaining about
% inefficient usage.  In all cases, the inefficient usage (of error, getfield, 
% setfield and find) is used to ensure compatibility with earlier versions
% of MATLAB.

remaining = nargin-1; % number of arguments other than X
count = 1;
fields = fieldnames(X);
modified = zeros(size(fields));
% Take input arguments two at a time until we run out.
while remaining>=2
    fieldname = varargin{count};
    fieldind = find(strcmp(fieldname,fields));
    if ~isempty(fieldind)
        oldvalue = getfield(X,fieldname); %#ok
        newvalue = varargin{count+1};
        if iscell(oldvalue)
            % Cell arrays must contain strings, and the new value must be
            % a string which appears in the list.
            if ~iscellstr(oldvalue)
                error(sprintf('All allowed values for "%s" must be strings',fieldname));  %#ok
            end
            if ~ischar(newvalue)
                error(sprintf('New value for "%s" must be a string',fieldname));  %#ok
            end
            if isempty(find(strcmp(oldvalue,newvalue))) %#ok
                error(sprintf('"%s" is not allowed for field "%s"',newvalue,fieldname));  %#ok
            end
        elseif ~isempty(oldvalue)
            % The caller isn't allowed to change the data type of a non-empty property,
            % and scalars must remain as scalars.
            if ~strcmp(class(oldvalue),class(newvalue))
                error(sprintf('Cannot change class of field "%s" from "%s" to "%s"',...
                    fieldname,class(oldvalue),class(newvalue))); %#ok
            elseif numel(oldvalue)==1 & numel(newvalue)~=1 %#ok
                error(sprintf('New value for "%s" must be a scalar',fieldname));  %#ok
            end
        end
        X = setfield(X,fieldname,newvalue); %#ok
        modified(fieldind) = 1;
    else
        error(['Not a valid field name: ' fieldname]);
    end
    remaining = remaining - 2;
    count = count + 2;
end
% Check that we had a value for every name.
if remaining~=0
    error('Odd number of arguments supplied.  Name-value pairs required');
end

% Now find cell arrays which were not modified by the above process, and select
% the first string.
notmodified = find(~modified);
for i=1:length(notmodified)
    fieldname = fields{notmodified(i)};
    oldvalue = getfield(X,fieldname); %#ok
    if iscell(oldvalue)
        if ~iscellstr(oldvalue)
            error(sprintf('All allowed values for "%s" must be strings',fieldname)); %#ok
        elseif isempty(oldvalue)
            error(sprintf('Empty cell array not allowed for field "%s"',fieldname)); %#ok
        end
        X = setfield(X,fieldname,oldvalue{1}); %#ok
    end
end
end







