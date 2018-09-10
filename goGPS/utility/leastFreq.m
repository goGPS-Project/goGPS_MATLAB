function [M,F,C] = leastFreq(x,dim)
%MODE   Mode, or most frequent value in a sample.
%   M=MODE(X) for vector X computes M as the sample mode, or most frequently
%   occurring value in X.  For a matrix X, M is a row vector containing
%   the mode of each column.  For N-D arrays, MODE(X) is the mode of the
%   elements along the first non-singleton dimension of X.
%
%   When there are multiple values occurring equally frequently, MODE
%   returns the smallest of those values.  For complex inputs, this is taken
%   to be the first value in a sorted list of values.
%
%   [M,F]=MODE(X) also returns an array F, of the same size as M.
%   Each element of F is the number of occurrences of the corresponding
%   element of M.
%
%   [M,F,C]=MODE(X) also returns a cell array C, of the same size
%   as M.  Each element of C is a sorted vector of all the values having
%   the same frequency as the corresponding element of M.
%
%   [...]=MODE(X,DIM) takes the mode along the dimension DIM of X.
%
%   This function is most useful with discrete or coarsely rounded data.
%   The mode for a continuous probability distribution is defined as
%   the peak of its density function.  Applying the MODE function to a
%   sample from that distribution is unlikely to provide a good estimate
%   of the peak; it would be better to compute a histogram or density
%   estimate and calculate the peak of that estimate.  Also, the MODE
%   function is not suitable for finding peaks in distributions having
%   multiple modes.
%
%   Example: If X = [3 3 1 4
%                    0 0 1 1
%                    0 1 2 4]
%
%   then mode(X) is [0 0 1 4] and mode(X,2) is [3
%                                               0
%                                               0]
%
%   To find the mode of a continuous variable grouped into bins:
%      y = randn(1000,1);
%      edges = -6:.25:6;
%      bin = discretize(y,edges);
%      m = mode(bin);
%      edges([m, m+1])
%      histogram(y,edges)
%
%   Class support for input X:
%      float:  double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also MEAN, MEDIAN, HISTOGRAM, HISTCOUNTS.

%   Copyright 2005-2015 The MathWorks, Inc.

dofreq = nargout>=2;
docell = nargout>=3;

if nargin<2
    % Special case to make mode and mean behave similarly
    if isequal(x, [])
        if isinteger(x)
            M = zeros('like',x);
        else
            M = NaN('like',x);
        end
        if dofreq
            F = zeros(1,1,'like',full(double(x([]))));
        end
        if docell
            C = {zeros(0,1,'like',x)};
        end
        warning(message('MATLAB:mode:EmptyInput'))
        return
    end
    
    % Determine the first non-singleton dimension
    dim = find(size(x)~=1, 1);
    if isempty(dim)
        dim = 1;
    end
else
    if ~isscalar(dim)   || ~isnumeric(dim)  || dim~=floor(dim)  || ...
            dim<1       || ~isreal(dim)     || ~isfinite(dim)
        error(message('MATLAB:mode:BadDim'));
    end
end

sizex = size(x);
if dim>length(sizex)
    sizex(length(sizex)+1:dim) = 1;
end

sizem = sizex;
sizem(dim) = 1;

% Dispose of empty arrays right away
if isempty(x)
    M = zeros(sizem,'like',x);
    if prod(sizem)>0
        M(:) = NaN;
    end
    if dofreq
        F = zeros(sizem,'like',full(double(x([]))));
    end
    if docell
        C = cell(sizem);
        C(:) = {M(1:0)};  % fill C with empties of the proper type
    end
    return
end

if isvector(x) && dim <=2
    % Treat vectors separately
    if (iscolumn(x) && dim == 2) || (~iscolumn(x) && dim == 1)
        % No computation needed for mode(col,2) and mode(row,1)
        M = x;
        
        if dofreq
            F = ones(sizex,'like',full(double(x([]))));
        end
        if docell
            C = num2cell(x);
        end
    else
        % Sort the vector and compute the mode
        x = sort(x(:));
        % start of run of equal values
        start = find([1; x(1:end-1)~=x(2:end)]);
        % frequencies for each run (force to double datatype)
        freq = zeros(numel(x),1,'like',full(double(x([]))));
        freq(start) = [diff(start); numel(x)+1-start(end)];
        [maxfreq,firstloc] = min(freq);
        
        M = x(firstloc);                % Mode
        
        if dofreq
            F = maxfreq;                % Frequency
        end
        if docell
            C = {x(freq == maxfreq)};   % Cell array with modes
        end
    end
else
    % Permute data and reshape into a 2-D array
    perm = [dim, (1:dim-1), (dim+1:length(sizex))];
    sizem = sizem(perm);
    x = permute(x, perm);
    x = reshape(x,[sizex(dim),prod(sizem)]);
    [nrows,ncols] = size(x);
    
    % Compute the modes for each column of the 2-D array
    x = sort(x,1);
    % start of run of equal values
    start = [ones(1,ncols); x(1:end-1,:)~=x(2:end,:)];
    start = find(start(:));
    % frequencies for each run (force to double datatype)
    freq = zeros([nrows,ncols],'like',full(double(x([]))));
    freq(start) = [start(2:end); numel(x)+1]-start;
    [maxfreq,firstloc] = max(freq,[],1);
    
    M = x((0:nrows:numel(x)-1)+firstloc);           % Modes for each column
    M = ipermute(reshape(M,sizem), perm);           % Reshape and permute back
    
    if dofreq
        F = ipermute(reshape(maxfreq,sizem), perm); % Frequencies
    end
    if docell
        C = cell(size(M));                          % Cell array with modes
        selection = bsxfun(@eq, freq, maxfreq);
        for j = 1:numel(M)
            C{j} = x(selection(:,j),j);
        end
    end
end
