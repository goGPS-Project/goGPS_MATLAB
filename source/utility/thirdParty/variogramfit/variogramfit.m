function [a,c,n,S] = variogramfit(h,gammaexp,a0,c0,numobs,varargin)

% fit a theoretical variogram to an experimental variogram
%
% Syntax:
%
%     [a,c,n] = variogramfit(h,gammaexp,a0,c0)
%     [a,c,n] = variogramfit(h,gammaexp,a0,c0,numobs)
%     [a,c,n] = variogramfit(h,gammaexp,a0,c0,numobs,'pn','pv',...)
%     [a,c,n,S] = variogramfit(...)
%
% Description:
%
%     variogramfit performs a least squares fit of various theoretical 
%     variograms to an experimental, isotropic variogram. The user can
%     choose between various bounded (e.g. spherical) and unbounded (e.g.
%     exponential) models. A nugget variance can be modelled as well, but
%     higher nested models are not supported.
%
%     The function works best with the function fminsearchbnd available on
%     the FEX. You should download it from the File Exchange (File ID:
%     #8277). If you don't have fminsearchbnd, variogramfit uses
%     fminsearch. The problem with fminsearch is, that it might return 
%     negative variances or ranges.
%
%     The variogram fitting algorithm is in particular sensitive to initial
%     values below the optimal solution. In case you have no idea of
%     initial values variogramfit calculates initial values for you
%     (c0 = max(gammaexp); a0 = max(h)*2/3;). If this is a reasonable
%     guess remains to be answered. Hence, visually inspecting your data
%     and estimating a theoretical variogram by hand should always be
%     your first choice.
%
%     Note that for unbounded models, the supplied parameter a0 (range) is
%     the distance where gamma equals 95% of the sill variance. The
%     returned parameter a0, however, is the parameter r in the model. The
%     range at 95% of the sill variance is then approximately 3*r.
%
% Input arguments:
%
%     h         lag distance of the experimental variogram
%     gammaexp  experimental variogram values (gamma)
%     a0        initial value (scalar) for range
%     c0        initial value (scalar) for sill variance
%     numobs    number of observations per lag distance (used for weight
%               function)
%
% Output arguments:
%
%     a         range
%     c         sill
%     n         nugget (empty if nugget variance is not applied)
%     S         structure array with additional information
%               .range
%               .sill
%               .nugget
%               .model - theoretical variogram
%               .func - anonymous function of variogram model (only the
%               function within range for bounded models)
%               .h  - distance
%               .gamma  - experimental variogram values
%               .gammahat - estimated variogram values
%               .residuals - residuals
%               .Rs - R-square of goodness of fit
%               .weights - weights
%               .exitflag - see fminsearch
%               .algorithm - see fminsearch
%               .funcCount - see fminsearch
%               .iterations - see fminsearch
%               .message - see fminsearch
%
% Property name/property values:
% 
%     'model'   a string that defines the function that can be fitted 
%               to the experimental variogram. 
% 
%               Supported bounded functions are:
%               'blinear' (bounded linear) 
%               'circular' (circular model)
%               'spherical' (spherical model, =default)
%               'pentaspherical' (pentaspherical model)
% 
%               Supported unbounded functions are:
%               'exponential' (exponential model)
%               'gaussian' (gaussian variogram)
%               'whittle' Whittle's elementary correlation (involves a
%                         modified Bessel function of the second kind.
%               'stable' (stable models sensu Wackernagel 1995). Same as
%                         gaussian, but with different exponents. Supply 
%                         the exponent alpha (<2) in an additional pn,pv 
%                         pair: 
%                        'stablealpha',alpha (default = 1.5).
%               'matern' Matern model. Requires an additional pn,pv pair. 
%                        'nu',nu (shape parameter > 0, default = 1)
%                        Note that for particular values of nu the matern 
%                        model reduces to other authorized variogram models.
%                        nu = 0.5 : exponential model
%                        nu = 1 : Whittles model
%                        nu -> inf : Gaussian model
%               
%               See Webster and Oliver (2001) for an overview on variogram 
%               models. See Minasny & McBratney (2005) for an introduction
%               to the Matern variogram.
%           
%     'nugget'  initial value (scalar) for nugget variance. The default
%               value is []. In this case variogramfit doesn't fit a nugget
%               variance. 
% 
%     'plotit'  true (default), false: plot experimental and theoretical 
%               variogram together.
% 
%     'solver'  'fminsearchbnd' (default) same as fminsearch , but with  
%               bound constraints by transformation (function by John 
%               D'Errico, File ID: #8277 on the FEX). The advantage in 
%               applying fminsearchbnd is that upper and lower bound 
%               constraints can be applied. That prevents that nugget 
%               variance or range may become negative.           
%               'fminsearch'
%
%     'weightfun' 'none' (default). 'cressie85' and 'mcbratney86' require
%               you to include the number of observations per experimental
%               gamma value (as returned by VARIOGRAM). 
%               'cressie85' uses m(hi)/gammahat(hi)^2 as weights
%               'mcbratney86' uses m(hi)*gammaexp(hi)/gammahat(hi)^3
%               
%
% Example: fit a variogram to experimental data
%
%     load variogramexample
%     a0 = 15; % initial value: range 
%     c0 = 0.1; % initial value: sill 
%     [a,c,n] = variogramfit(h,gammaexp,a0,c0,[],...
%                            'solver','fminsearchbnd',...
%                            'nugget',0,...
%                            'plotit',true);
%
%           
% See also: VARIOGRAM, FMINSEARCH, FMINSEARCHBND
%           
%
% References: Wackernagel, H. (1995): Multivariate Geostatistics, Springer.
%             Webster, R., Oliver, M. (2001): Geostatistics for
%             Environmental Scientists. Wiley & Sons.
%             Minsasny, B., McBratney, A. B. (2005): The Matérn function as
%             general model for soil variograms. Geoderma, 3-4, 192-207.
% 
% Date: 7. October, 2010
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)




% check input arguments

if nargin == 0
    help variogramfit
    return
elseif nargin>0 && nargin < 2;
    error('Variogramfit:inputargs',...
          'wrong number of input arguments');
end
if ~exist('a0','var') || isempty(a0)
    a0 = max(h)*2/3;
end
if ~exist('c0','var') || isempty(c0)
    c0 = max(gammaexp);
end
if ~exist('numobs','var') || isempty(a0)
    numobs = [];
end
      

% check input parameters
params.model       = 'spherical';
params.nugget      = [];
params.plotit      = true;
params.solver      = {'fminsearchbnd','fminsearch'};
params.stablealpha = 1.5;
params.weightfun   = {'none','cressie85','mcbratney86'};
params.nu          = 1;
params = parseargs(params,varargin{:});

% check if fminsearchbnd is in the search path
switch lower(params.solver)
    case 'fminsearchbnd'
        if ~exist('fminsearchbnd.m','file')==2
            params.solver = 'fminsearch';
            warning('Variogramfit:fminsearchbnd',...
            'fminsearchbnd was not found. fminsearch is used instead')
        end
end

% check if h and gammaexp are vectors and have the same size
if ~isvector(h) || ~isvector(gammaexp)
    error('Variogramfit:inputargs',...
          'h and gammaexp must be vectors');
end

% force column vectors
h = h(:);
gammaexp = gammaexp(:);

% check size of supplied vectors 
if numel(h) ~= numel(gammaexp)
    error('Variogramfit:inputargs',...
          'h and gammaexp must have same size');
end

% remove nans;
nans = isnan(h) | isnan(gammaexp);
if any(nans);
    h(nans) = [];
    gammaexp(nans) = [];
    if ~isempty(numobs)
        numobs(nans) = [];
    end
end

% check weight inputs
if isempty(numobs);
    params.weightfun = 'none';
end
    

% create options for fminsearch
options = optimset('MaxFunEvals',1000000);

% create vector with initial values
% b(1) range
% b(2) sill
% b(3) nugget if supplied
b0 = [a0 c0 params.nugget];

% variogram function definitions
switch lower(params.model)    
    case 'spherical'
        type = 'bounded';
        func = @(b,h)b(2)*((3*h./(2*b(1)))-1/2*(h./b(1)).^3);
    case 'pentaspherical'
        type = 'bounded';
        func = @(b,h)b(2)*(15*h./(8*b(1))-5/4*(h./b(1)).^3+3/8*(h./b(1)).^5);
    case 'blinear'
        type = 'bounded';
        func = @(b,h)b(2)*(h./b(1));
    case 'circular'
        type = 'bounded';
        func = @(b,h)b(2)*(1-(2./pi)*acos(h./b(1))+2*h/(pi*b(1)).*sqrt(1-(h.^2)/(b(1)^2)));
    case 'exponential'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-exp(-h./b(1)));
    case 'gaussian'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-exp(-(h.^2)/(b(1)^2)));
    case 'stable'
        type = 'unbounded';
        stablealpha = params.stablealpha;
        func = @(b,h)b(2)*(1-exp(-(h.^stablealpha)/(b(1)^stablealpha)));  
    case 'whittle'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-h/b(1).*besselk(1,h/b(1)));
    case 'matern'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-(1/((2^(params.nu-1))*gamma(params.nu))) * (h/b(1)).^params.nu .* besselk(params.nu,h/b(1)));
    otherwise
        error('unknown model')
end


% check if there are zero distances 
% if yes, remove them, since the besselk function returns nan for
% zero
switch lower(params.model) 
    case {'whittle','matern'}
        izero = h==0;
        if any(izero)
            flagzerodistances = true;
        else
            flagzerodistances = false;
        end
    otherwise
        flagzerodistances = false;
end
        

% if model type is unbounded, then the parameter b(1) is r, which is
% approximately range/3. 
switch type
    case 'unbounded'
        b0(1) = b0(1)/3;
end


% nugget variance
if isempty(params.nugget)
    nugget = false;
    funnugget = @(b) 0;
else
    nugget = true;
    funnugget = @(b) b(3);
end

% generate upper and lower bounds when fminsearchbnd is used
switch lower(params.solver)
    case {'fminsearchbnd'};
        % lower bounds
        lb = zeros(size(b0));
        % upper bounds
        if nugget;
            ub = [inf max(gammaexp) max(gammaexp)]; %
        else
            ub = [inf max(gammaexp)];
        end
end

% create weights (see Webster and Oliver)
switch params.weightfun
    case 'cressie85'
        weights = @(b,h) (numobs./variofun(b,h).^2)./sum(numobs./variofun(b,h).^2);
    case 'mcbratney86'
        weights = @(b,h) (numobs.*gammaexp./variofun(b,h).^3)/sum(numobs.*gammaexp./variofun(b,h).^3);
    otherwise
        weights = @(b,h) 1;
end


% create objective function: weighted least square
objectfun = @(b)sum(((variofun(b,h)-gammaexp).^2).*weights(b,h));

% call solver
switch lower(params.solver)
    case 'fminsearch'                
        % call fminsearch
        [b,fval,exitflag,output] = fminsearch(objectfun,b0,options);
    case 'fminsearchbnd'
        % call fminsearchbnd
        [b,fval,exitflag,output] = fminsearchbnd(objectfun,b0,lb,ub,options);
    otherwise
        error('Variogramfit:Solver','unknown or unsupported solver')
end


% prepare output
a = b(1); %range
c = b(2); %sill
if nugget;
    n = b(3);%nugget
else
    n = [];
end


% Create structure array with results 
if nargout == 4;    
    S.model     = lower(params.model); % model
    S.func      = func;
    S.type      = type;
    switch S.model 
        case 'matern';
            S.nu = params.nu;
        case 'stable';
            S.stablealpha = params.stablealpha;
    end
        
    
    S.range     = a;
    S.sill      = c;
    S.nugget    = n;
    S.h         = h; % distance
    S.gamma     = gammaexp; % experimental values
    S.gammahat  = variofun(b,h); % estimated values
    S.residuals = gammaexp-S.gammahat; % residuals
    COVyhaty    = cov(S.gammahat,gammaexp);
    S.Rs        = (COVyhaty(2).^2) ./...
                  (var(S.gammahat).*var(gammaexp)); % Rsquare
    S.weights   = weights(b,h); %weights
    S.weightfun = params.weightfun;
    S.exitflag  = exitflag; % exitflag (see doc fminsearch)
    S.algorithm = output.algorithm;
    S.funcCount = output.funcCount;
    S.iterations= output.iterations;
    S.message   = output.message;
end



% if you want to plot the results...
if params.plotit
    switch lower(type)
        case 'bounded'
            plot(h,gammaexp,'rs','MarkerSize',10);
            hold on
            fplot(@(h) funnugget(b) + func(b,h),[0 b(1)])
            fplot(@(h) funnugget(b) + b(2),[b(1) max(h)])
            
        case 'unbounded'
            plot(h,gammaexp,'rs','MarkerSize',10);
            hold on
            fplot(@(h) funnugget(b) + func(b,h),[0 max(h)])
    end
    axis([0 max(h) 0 max(gammaexp)])
    xlabel('lag distance h')
    ylabel('\gamma(h)')
    hold off
end


% fitting functions for  fminsearch/bnd
function gammahat = variofun(b,h)
    
    switch type
        % bounded model
        case 'bounded'
            I = h<=b(1);
            gammahat     = zeros(size(I));
            gammahat(I)  = funnugget(b) + func(b,h(I));
            gammahat(~I) = funnugget(b) + b(2);        
        % unbounded model
        case 'unbounded'
            gammahat = funnugget(b) + func(b,h);
            if flagzerodistances
                gammahat(izero) = funnugget(b);
            end    
    end
end

end

% and that's it...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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