function d = ipdm(data1,varargin)
% ipdm: Inter-Point Distance Matrix
% usage: d = ipdm(data1)
% usage: d = ipdm(data1,data2)
% usage: d = ipdm(data1,prop,value)
% usage: d = ipdm(data1,data2,prop,value)
%
% Arguments: (input)
%  data1 - array of data points, each point is one row. p dimensional
%          data will be represented by matrix with p columns.
%          If only data1 is provided, then the distance matrix
%          is computed between all pairs of rows of data1.
%
%          If your data is one dimensional, it MUST form a column
%          vector. A row vector of length n will be interpreted as
%          an n-dimensional data set.
%
%  data2 - second array, supplied only if distances are to be computed
%          between two sets of points.
%
%
% Class support: data1 and data2 are assumed to be either
% single or double precision. I have not tested this code to
% verify its success on integer data of any class.
%
%
% Additional parameters are expected to be property/value pairs.
% Property/value pairs are pairs of arguments, the first of which
% (properties) must always be a character string. These strings
% may be shortened as long as the shortening is unambiguous.
% Capitalization is ignored. Valid properties for ipdm are:
%
%  'Metric', 'Subset', 'Limit', 'Result'
%
%  'Metric' - numeric flag - defines the distance metric used
%          metric = 2 --> (DEFAULT) Euclidean distance = 2-norm
%                         The standard distance metric.
%
%          metric = 1 --> 1-norm = sum of absolute differences
%                         Also sometimes known as the "city block
%                         metric", since this is the sum of the
%                         differences in each dimension.
%
%          metric = inf --> infinity-norm = maximum difference
%                         over all dimensions. The name refers
%                         to the limit of the p-norm, as p
%                         approaches infinity.
%
%          metric = 0 --> minimum difference over all dimensions.
%                         This is not really a useful norm in
%                         practice.
%
%          metric = 5 --> use spherical coordinates.
%
%          Note: while other distance metrics exist, IMHO, these
%          seemed to be the common ones.
%
%
%  'Result' - A string variable that denotes the style of returned
%          result. Valid result types are 'Array', 'Structure'.
%          Capitalization is ignored, and the string may be
%          shortened if you wish.
%
%          result = 'Array' --> (DEFAULT) A matrix of all
%                         interpoint distances will be generated.
%                         This array may be large. If this option
%                         is specified along with a minimum or
%                         maximum value, then those elements above
%                         or below the limiting values will be
%                         set as -inf or +inf, as appropriate.
%
%                         When any of 'LargestFew', 'SmallestFew',
%                         or 'NearestNeighbor' are set, then the
%                         resulting array will be a sparse matrix
%                         if 'array' is specified as the result.
%
%          result = 'Structure' --> A list of all computed distances,
%                         defined as a structure. This structure
%                         will have fields named 'rowindex',
%                         'columnindex', and 'distance'.
%
%                         This option will be useful when a subset
%                         criterion for the distances has been
%                         specified, since then the distance matrix
%                         may be very sparsely populated. Distances
%                         for pairs outside of the criterion will
%                         not be returned.
%
%
%  'Subset' - Character string, any of:
%
%          'All', 'Maximum', 'Minimum', 'LargestFew', 'SmallestFew',
%          'NearestNeighbor', 'FarthestNeighbor', or empty
%
%          Like properties, capitalization is ignored here, and
%          any unambiguous shortening of the word is acceptable.
%
%          DEFAULT = 'All'
%
%          Some interpoint distance matrices can be huge. Often
%          these matrices are too large to be fully retained in
%          memory, yet only the pair of points with the largest
%          or smallest distance may be needed. When only some
%          subset of the complete set of distances is of interest,
%          these options allow you to specify which distances will
%          be returned.
%
%          If 'result' is defined to be an array, then a sparse
%          matrix will be returned for the 'LargestFew', 'SmallestFew',
%          'NearestNeighbor', and 'FarthestNeighbor' subset classes.
%          'Minimum' and 'Maximum' will yield full matrices by
%          default. If a structure is specified, then only those
%          elements which have been identified will be returned.
%
%          Where a subset is specified, its limiting value is
%          specified by the 'Limit' property. Call that value k.
%
%
%          'All' -->     (DEFAULT) Return all interpoint distances
%
%          'Minimum' --> Only look for those distances above
%                        the cutoff k. All other distances will
%                        be returned as -inf.
%
%          'Maximum' --> Only look for those distances below
%                        the cutoff k. All other distances will
%                        be returned as +inf.
%
%          'SmallestFew' --> Only return the subset of the k
%                        smallest distances. Where only one data
%                        set is provided, only the upper triangle
%                        of the inter-point distance matrix will
%                        be generated since that matrix is symmetric.
%
%          'LargestFew' --> Only return the subset of the k
%                        largest distances. Where only one data
%                        set is provided, only the upper triangle
%                        of the inter-point distance matrix will
%                        be generated since that matrix is symmetric.
%
%          'NearestNeighbor' --> Only return the single nearest
%                        neighbor in data2 to each point in data1.
%                        No limiting value is required for this
%                        option. If multiple points have the same
%                        nearest distance, then return the first
%                        such point found. With only one input set,
%                        a point will not be its own nearest
%                        neighbor.
%
%                        Note that exact replicates in a single set
%                        will cause problems, since a sparse matrix
%                        is returned by default. Since they will have
%                        a zero distance, they will not show up in
%                        the sparse matrix. A structure return will
%                        show those points as having a zero distance
%                        though.
%
%          'FarthestNeighbor' --> Only return the single farthest
%                        neighbor to each point. No limiting value
%                        is required for this option. If multiple
%                        points have the same farthest distance,
%                        then return the first such point found.
%
%
%  'Limit' - scalar numeric value or []. Used only when some
%           Subset is specified.
%
%          DEFAULT = []
%
%
%  'ChunkSize' - allows a user with lower RAM limits
%          to force the code to only grab smaller chunks of RAM
%          at a time (where possible). This parameter is specified
%          in bytes of RAM. The default is 32 megabytes, or 2^22
%          elements in any piece of the distance matrix. Only some
%          options will break the problem into chunks, thus as long
%          as a full matrix is expected to be returned, there seems
%          no reason to break the problem up into pieces.
%
%          DEFAULT = 2^25
%
%
% Arguments: (output)
%  d     - array of interpoint distances, or a struct wth the
%          fields {'rowindex', 'columnindex', 'distance'}.
%
%          d(i,j) represents the distance between point i
%          (from data1) and point j (from data2).
%
%          If only one (n1 x p) array is supplied, then d will
%          be an array of size == [n1,n1].
%
%          If two arrays (of sizes n1 x p and n2 x p) then d
%          will be an array of size == [n1,n2].
%
%
% Efficiency considerations:
%  Where possible, this code will use bsxfun to compute its
%  distances.
%
%
% Example:
%  Compute the interpoint distances between all pairs of points
%  in a list of 5 points, in 2 dimensions and using Euclidean
%  distance as the distance metric.
%
%  A = randn(5,2);
%  d = ipdm(A,'metric',2)
%  d =
%            0       2.3295       3.2263       2.0263       2.8244
%       2.3295            0       1.1485      0.31798       1.0086
%       3.2263       1.1485            0       1.4318       1.8479
%       2.0263      0.31798       1.4318            0       1.0716
%       2.8244       1.0086       1.8479       1.0716            0
%
% (see the demo file for many other examples)
%
% See also: pdist
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/26/08

% Default property values
params.Metric = 2;
params.Result = 'array';
params.Subset = 'all';
params.Limit = [];
params.ChunkSize = 2^25;

% untangle the arguments
if nargin<1
  % if called with no arguments, then the user probably
  % needs help. Give it to them.
  help ipdm
  return
end

% were two sets of data provided?
pvpairs = {};
if nargin==1
  % only 1 set of data provided
  dataflag = 1;
  data2 = [];
else
  if ischar(varargin{1})
    dataflag = 1;
    data2 = [];
    pvpairs = varargin;
  else
    dataflag = 2;
    data2 = varargin{1};
    if nargin>2
      pvpairs = varargin(2:end);
    end
  end
end

% get data sizes for later
[n1,dim] = size(data1);
if dataflag == 2
  n2 = size(data2,1);
end

% Test the class of the input variables
if ~(isa(data1,'double') || isa(data1,'single')) || ...
    ((dataflag == 2) && ~(isa(data2,'double') || isa(data2,'single')))
  error('data points must be either single or double precision variables.')
end

% do we need to process any property/value pairs?
if nargin>2
  params = parse_pv_pairs(params,pvpairs);
  
  % check for problems in the properties
  
  % was a legal Subset provided?
  if ~isempty(params.Subset) && ~ischar(params.Subset)
    error('If provided, ''Subset'' must be character')
  elseif isempty(params.Subset)
    params.Subset = 'all';
  end
  valid = {'all','maximum','minimum','largestfew','smallestfew', ...
    'nearestneighbor','farthestneighbor'};
  ind = find(strncmpi(params.Subset,valid,length(params.Subset)));
  if (length(ind)==1)
    params.Subset = valid{ind};
  else
    error(['Invalid Subset: ',params.Subset])
  end
  
  % was a limit provided?
  if ~ismember(params.Subset,{'all','nearestneighbor','farthestneighbor'}) && ...
      isempty(params.Limit)
    error('No limit provided, but a Subset that requires a limit value was specified')
  end
  % check the limit values for validity
  if length(params.Limit)>1
    error('Limit must be scalar or empty')
  end
  
  switch params.Subset
    case {'largestfew', 'smallestfew'}
      % must be at least 1, and an integer
      if (params.Limit<1) || (round(params.Limit)~=params.Limit)
        error('Limit must be a positive integer for LargestFew or NearestFew')
      end
  end
  
  % was a legal Result provided?
  if isempty(params.Result)
    params.result = 'Array';
  elseif ~ischar(params.Result)
    error('If provided, ''Result'' must be character or empty')
  end
  valid = {'array','structure'};
  ind = find(strncmpi(params.Result,valid,length(params.Result)));
  if (length(ind)==1)
    params.Result = valid{ind};
  else
    error(['Invalid Result: ',params.Subset])
  end

  % check for the metric
  if isempty(params.Metric)
    params.Metric = 2;
  elseif (length(params.Metric)~=1) || ~ismember(params.Metric,[0 1 2 inf 5])
    error('If supplied, ''Metric'' must be a scalar, and one of [0 1 2 inf]')
  end
end % if nargin>2
  
% If Metric was given as 2, but the dimension is only 1, then it will
% be slightly faster (and equivalent) to use the 1-norm Metric.
if (dim == 1) && (params.Metric == 2)
  params.Metric = 1;
end

% Can we use bsxfun to compute the interpoint distances?
% Older Matlab releases will not have bsxfun, but if it is
% around, it will ne both faster and less of a memory hog.
params.usebsxfun = (5==exist('bsxfun','builtin'));

% check for dimension mismatch if 2 sets
if (dataflag==2) && (size(data2,2)~=dim)
  error('If 2 point sets provided, then both must have the same number of columns')
end

% Total number of distances to compute, in case I must do it in batches
if dataflag==1
  n2 = n1;
end
ntotal = n1*n2;

% FINALLY!!! Compute inter-point distances
switch params.Subset
  case 'all'
    % The complete set of interpoint distances. There is no need
    % to break this into chunks, since we must return all distances.
    % If that is too much to compute in memory, then it will fail
    % anyway when we try to store the result. bsxfun will at least
    % do the computation efficiently.
    
    % One set or two?
    if dataflag == 1
      d = distcomp(data1,data1,params);
    else
        d = distcomp(data1,data2,params);
    end
    
    % Must we return it as a struct?
    if params.Result(1) == 's'
      [rind,cind] = ndgrid(1:size(d,1),1:size(d,2));
      ds.rowindex = rind(:);
      ds.columnindex = cind(:);
      ds.distance = d(:);
      d = ds;
    end
    
  case {'minimum' 'maximum'}
    % There is no reason to break this into pieces if the result
    % sill be filled in the end with +/- inf. Only break it up
    % if the final result is a struct.
    if ((ntotal*8)<=params.ChunkSize) || (params.Result(1) == 'a')
      % its small enough to do it all at once
      
      % One set or two?
      if dataflag == 1
        d = distcomp(data1,data1,params);
      else
          d = distcomp(data1,data2,params);
      end
      
      % Must we return it as a struct?
      if params.Result(1) == 'a'
        % its an array, fill the unwanted distances with +/- inf
        if params.Subset(2) == 'i'
          % minimum
          d(d<=params.Limit) = -inf;
        else
          % maximum
          d(d>=params.Limit) = +inf;
        end
      else
        % a struct will be returned
        if params.Subset(2) == 'i'
          % minimum
          [dist.rowindex,dist.columnindex] = find(d>=params.Limit);
        else
          % maximum
          [dist.rowindex,dist.columnindex] = find(d<=params.Limit);
        end
        dist.distance = d(dist.rowindex + n1*(dist.columnindex-1));
        d = dist;
      end
      
    else
      % we need to break this into chunks. This branch
      % will always return a struct.
      
      % this is the number of rows of data1 that we will
      % process at a time.
      bs = floor(params.ChunkSize/(8*n2));
      bs = min(n1,max(1,bs));
      
      % Accumulate the result into a cell array. Do it this
      % way because we don't know in advance how many elements
      % that we will find satisfying the minimum or maximum
      % limit specified.
      accum = cell(0,1);
      
      % now loop over the chunks
      batch = 1:bs;
      while ~isempty(batch)
        
        % One set or two?
        if dataflag == 1
          dist = distcomp(data1(batch,:),data1,params);
        else
          dist = distcomp(data1(batch,:),data2,params);
        end
        
        % big or small as requested
        if ('i'==params.Subset(2))
          % minimum value specified
          [I,J,V] = find(dist>=params.Limit);
        else
          % maximum limit
          [I,J] = find(dist<=params.Limit);
          I = I(:);
          J = J(:);
          V = dist(I + (J-1)*length(batch));
          I = I + (batch(1)-1);
        end

        % and stuff them into the cell structure
        if ~isempty(V)
          accum{end+1,1} = [I,J,V(:)]; %#ok
        end

        % increment the batch
        batch = batch + bs;
        if batch(end)>n1
          batch(batch>n1) = [];
        end

      end

      % convert the cells into one flat array
      accum = cell2mat(accum);

      if isempty(accum)
        d.rowindex = [];
        d.columnindex = [];
        d.distance = [];
      else
        % we found something

        % sort on the second column, to put them in a reasonable order
        accum = sortrows(accum,[2 1]);

        d.rowindex = accum(:,1);
        d.columnindex = accum(:,2);
        d.distance = accum(:,3);
      end

    end

  case {'smallestfew' 'largestfew'}
    % find the k smallest/largest distances. k is
    % given by params.Limit

    % if only 1 set, params.Limit must be less than n*(n-1)/2
    if dataflag == 1
      params.Limit = min(params.Limit,n1*(n1-1)/2);
    end

    % is this a large problem?
    if ((ntotal*8) <= params.ChunkSize)
      % small potatoes

      % One set or two?
      if dataflag == 1
        dist = distcomp(data1,data1,params);
        % if only one data set, set the diagonal and
        % below that to +/- inf so we don't find it.
        temp = find(tril(ones(n1,n1),0));
        if params.Subset(1) == 's'
          dist(temp) = inf;
        else
          dist(temp) = -inf;
        end
      else
        dist = distcomp(data1,data2,params);
      end

      % sort the distances to find those we need
      if ('s'==params.Subset(1))
        % smallestfew
        [val,tags] = sort(dist(:),'ascend');
      else
        % largestfew
        [val,tags] = sort(dist(:),'descend');
      end
      val = val(1:params.Limit);
      tags = tags(1:params.Limit);

      % recover the row and column index from the linear
      % index returned by sort in tags.
      [d.rowindex,d.columnindex] = ind2sub([n1,size(dist,2)],tags);

      % create the matrix as a sparse one or a struct?
      if params.Result(1)=='a'
        % its an array, so make the array sparse.
        d = sparse(d.rowindex,d.columnindex,val,n1,size(dist,2));
      else
        % a structure
        d.distance = val;
      end

    else
      % chunks

      % this is the number of rows of data1 that we will
      % process at a time.
      bs = floor(params.ChunkSize/(8*n2));
      bs = min(n1,max(1,bs));

      % We need to find the extreme cases. There are two possible
      % algorithms, depending on how many total elements we will
      % search for.
      % 1. Only a very few total elements.
      % 2. A relatively large number of total elements, forming
      %    a significant fraction of the total set.
      %
      % Case #1 would suggest to retain params.Limit numberr of
      % elements from each batch, then at the end, sort them all
      % to find the best few. Case #2 will result in too many
      % elements to retain, so we must distinguish between these
      % alternatives.
      if (8*params.Limit*n1/bs) <= params.ChunkSize
        % params.Limit is small enough to fall into case #1.

        % Accumulate the result into a cell array. Do it this
        % way because we don't know in advance how many elements
        % that we will find satisfying the minimum or maximum
        % limit specified.
        accum = cell(0,1);

        % now loop over the chunks
        batch = (1:bs)';
        while ~isempty(batch)
          % One set or two?
          if dataflag == 1
            dist = distcomp(data1(batch,:),data1,params);
            k = find(tril(ones(length(batch),n2),batch(1)-1));
            if ('s'==params.Subset(1))
              dist(k) = inf;
            else
              dist(k) = -inf;
            end
          else
            dist = distcomp(data1(batch,:),data2,params);
          end

          % big or small as requested, keeping only the best
          % params.Limit number of elements
          if ('s'==params.Subset(1))
            % minimum value specified
            [tags,tags] = sort(dist(:),1,'ascend');
            tags = tags(1:bs);
            [I,J] = ndgrid(batch,1:n2);
            ijv = [I(tags),J(tags),dist(tags)];
          else
            % maximum limit
            [tags,tags] = sort(dist(:),1,'descend');
            tags = tags(1:bs);
            [I,J] = ndgrid(batch,1:n2);
            ijv = [I(tags),J(tags),dist(tags)];
          end
          % and stuff them into the cell structure
          accum{end+1,1} = ijv; %#ok

          % increment the batch
          batch = batch + bs;
          if batch(end)>n1
            batch(batch>n1) = [];
          end
        end

        % convert the cells into one flat array
        accum = cell2mat(accum);

        % keep only the params.Limit best of those singled out
        accum = sortrows(accum,3);
        if ('s'==params.Subset(1))
          % minimum value specified
          accum = accum(1:params.Limit,:);
        else
          % minimum value specified
          accum = accum(end + 1 - (1:params.Limit),:);
        end
        d.rowindex = accum(:,1);
        d.columnindex = accum(:,2);
        d.distance = accum(:,3);

        % create the matrix as a sparse one or a struct?
        if params.Result(1)=='a'
          % its an array, so make the array sparse.
          d = sparse(d.rowindex,d.columnindex,d.distance,n1,size(dist,2));
        end

      else
        % params.Limit forces us into the domain of case #2.
        % Here we cannot retain params.Limit elements from each chunk.
        % so we will grab each chunk and append it to the best elements
        % found so far, then filter out the best after each chunk is
        % done. This may be slower than we want, but its the only way.
        ijv = zeros(0,3);

        % loop over the chunks
        batch = (1:bs)';
        while ~isempty(batch)
          % One set or two?
          if dataflag == 1
            dist = distcomp(data1(batch,:),data1,params);
            k = find(tril(ones(length(batch),n2),batch(1)-1));
            if ('s'==params.Subset(1))
              dist(k) = inf;
            else
              dist(k) = -inf;
            end
          else
            dist = distcomp(data1(batch,:),data2,params);
          end

          [I,J] = ndgrid(batch,1:n2);
          ijv = [ijv;[I(:),J(:),dist(:)]]; %#ok

          % big or small as requested, keeping only the best
          % params.Limit number of elements
          if size(ijv,1) > params.Limit
            if ('s'==params.Subset(1))
              % minimum value specified
              [tags,tags] = sort(ijv(:,3),1,'ascend');
            else
              [tags,tags] = sort(ijv(:,3),1,'ascend');
            end
            ijv = ijv(tags(1:params.Limit),:);
          end

          % increment the batch
          batch = batch + bs;
          if batch(end)>n1
            batch(batch>n1) = [];
          end
        end

        % They are fully trimmed down. stuff a structure
        d.rowindex = ijv(:,1);
        d.columnindex = ijv(:,2);
        d.distance = ijv(:,3);

        % create the matrix as a sparse one or a struct?
        if params.Result(1)=='a'
          % its an array, so make the array sparse.
          d = sparse(d.rowindex,d.columnindex,d.distance,n1,size(dist,2));
        end

      end

    end
    
  case {'nearestneighbor' 'farthestneighbor'}
    % find the closest/farthest neighbor for every point

    % is this a large problem? Or a 1-d problem?
    if dim == 1
      % its a 1-d nearest/farthest neighbor problem. we can
      % special case these easily enough, and all the distance
      % metric options are the same in 1-d.

      % first split it into the farthest versus nearest cases.
      if params.Subset(1) == 'f'
        % farthest away

        % One set or two?
        if dataflag == 1
          [d2min,minind] = min(data1);
          [d2max,maxind] = max(data1);
        else
          [d2min,minind] = min(data2);
          [d2max,maxind] = max(data2);
        end

        d.rowindex = (1:n1)';
        d.columnindex = repmat(maxind,n1,1);
        d.distance = repmat(d2max,n1,1);

        % which endpoint was further away?
        k = abs((data1 - d2min)) >= abs((data1 - d2max));
        if any(k)
          d.columnindex(k) = minind;
          d.distance(k) = d2min;
        end

      else
        % nearest. this is mainly a sort and some fussing around.
        d.rowindex = (1:n1)';
        d.columnindex = ones(n1,1);
        d.distance = zeros(n1,1);

        % One set or two?
        if dataflag == 1
          % if only one data point, then we are done
          if n1 == 2
            % if exactly two data points, its trivial
            d.columnindex = [2 1];
            d.distance = repmat(abs(diff(data1)),2,1);
          elseif n1>2
            % at least three points. do a sort.
            [sorted_data,tags] = sort(data1);

            % handle the first and last points separately
            d.columnindex(tags(1)) = tags(2);
            d.distance(tags(1)) = sorted_data(2) - sorted_data(1);
            d.columnindex(tags(end)) = tags(end-1);
            d.distance(tags(end)) = sorted_data(end) - sorted_data(end-1);

            ind = (2:(n1-1))';

            d1 = sorted_data(ind) - sorted_data(ind-1);
            d2 = sorted_data(ind+1) - sorted_data(ind);

            k = d1 < d2;
            d.distance(tags(ind(k))) = d1(k);
            d.columnindex(tags(ind(k))) = tags(ind(k)-1);
            k = ~k;
            d.distance(tags(ind(k))) = d2(k);
            d.columnindex(tags(ind(k))) = tags(ind(k)+1);
          end % if n1 == 2
        else
          % Two sets of data. still really a sort and some fuss.
          if n2 == 1
            % there is only one point in data2
            d.distance = abs(data1 - data2);
            % d.columnindex is already set correctly
          else
            % At least two points in data2
            % We need to sort all the data points together, but also
            % know which points from each set went where. ind12 and
            % bool12 will help keep track.
            ind12 = [1:n1,1:n2]';
            bool12 = [zeros(n1,1);ones(n2,1)];
            [sorted_data,tags] = sort([data1;data2]);

            ind12 = ind12(tags);
            bool12 = bool12(tags);

            % where did each point end up after the sort?
            loc1 = find(~bool12);
            loc2 = find(bool12);

            % for each point in data1, what is the (sorted) data2
            % element which appears most nearly to the left of it?
            cs = cumsum(bool12);
            leftelement = cs(loc1);

            % any points which fell below the minimum element in data2
            % will have a zero for the index of the element on their
            % left. fix this.
            leftelement = max(1,leftelement);

            % likewise, any point greater than the max in data2 will
            % have an n2 in left element. this too will be a problem
            % later, so fix it.
            leftelement = min(n2-1,leftelement);

            % distance to the left hand element
            dleft = abs(sorted_data(loc1) - sorted_data(loc2(leftelement)));
            dright = abs(sorted_data(loc1) - sorted_data(loc2(leftelement+1)));

            % find the points which are closer to the left element in data2
            k = (dleft < dright);
            d.distance(ind12(loc1(k))) = dleft(k);
            d.columnindex(ind12(loc1(k))) = ind12(loc2(leftelement(k)));
            k = ~k;
            d.distance(ind12(loc1(k))) = dright(k);
            d.columnindex(ind12(loc1(k))) = ind12(loc2(leftelement(k)+1));

          end % if n2 == 1
        end % if dataflag == 1
      end % if params.Subset(1) == 'f'

      % create the matrix as a sparse one or a struct?
      if params.Result(1)=='a'
        % its an array, so make the array sparse.
        d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
      end

    elseif (ntotal>1000) && (((params.Metric == 0) && (params.Subset(1) == 'n')) || ...
        ((params.Metric == inf) && (params.Subset(1) == 'f')))
      % nearest/farthest neighbour in n>1 dimensions, but for an
      % infinity norm metric. Reduce this to a sequence of
      % 1-d problems, each of which will be faster in general.
      % do this only if the problem is moderately large, since
      % we must overcome the extra overhead of the recursive
      % calls to ipdm.

      % do the first dimension
      if dataflag == 1
        d = ipdm(data1(:,1),data1(:,1),'subset',params.Subset,'metric',params.Metric,'result','struct');
      else
        d = ipdm(data1(:,1),data2(:,1),'subset',params.Subset,'metric',params.Metric,'result','struct');
      end

      % its slightly different for nearest versus farthest here
      % now, loop over dimensions
      for i = 2:dim
        if dataflag == 1
          di = ipdm(data1(:,i),data1(:,i),'subset',params.Subset,'metric',params.Metric,'result','struct');
        else
          di = ipdm(data1(:,i),data2(:,i),'subset',params.Subset,'metric',params.Metric,'result','struct');
        end

        % did any of the distances change?
        if params.Metric == 0
          % the 0 norm, with nearest neighbour, so take the
          % smallest distance in any dimension.
          k = d.distance > di.distance;
        else
          % inf norm. so take the largest distance across dimensions
          k = d.distance < di.distance;
        end

        if any(k)
          d.distance(k) = di.distance(k);
          d.columnindex(k) = di.columnindex(k);
        end
      end

      % create the matrix as a sparse one or a struct?
      if params.Result(1)=='a'
        % its an array, so make the array sparse.
        d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
      end

    elseif ((ntotal*8) <= params.ChunkSize)
      % None of the other special cases apply, so do it using brute
      % force for the small potatoes problem.

      % One set or two?
      if dataflag == 1
        dist = distcomp(data1,data1,params);
      else
        dist = distcomp(data1,data2,params);
      end

      % if only one data set and if a nearest neighbor
      % problem, set the diagonal to +inf so we don't find it.
      if (dataflag==1) && (n1>1) && ('n'==params.Subset(1))
        diagind = (1:n1) + (0:n1:(n1^2-1));
        dist(diagind) = +inf;
      end

      if ('n'==params.Subset(1))
        % nearest
        [val,j] = min(dist,[],2);
      else
        % farthest
        [val,j] = max(dist,[],2);
      end

      % create the matrix as a sparse one or a struct?
      if params.Result(1)=='a'
        % its an array, so make the array sparse.
        d = sparse((1:n1)',j,val,n1,size(dist,2));
      else
        % a structure
        d.rowindex = (1:n1)';
        d.columnindex = j;
        d.distance = val;
      end

    else

      % break it into chunks
      bs = floor(params.ChunkSize/(8*n2));
      bs = min(n1,max(1,bs));

      % pre-allocate the result
      d.rowindex = (1:n1)';
      d.columnindex = zeros(n1,1);
      d.distance = zeros(n1,1);

      % now loop over the chunks
      batch = 1:bs;
      while ~isempty(batch)

        % One set or two?
        if dataflag == 1
          dist = distcomp(data1(batch,:),data1,params);
        else
          dist = distcomp(data1(batch,:),data2,params);
        end

        % if only one data set and if a nearest neighbor
        % problem, set the diagonal to +inf so we don't find it.
        if (dataflag==1) && (n1>1) && ('n'==params.Subset(1))
          diagind = 1:length(batch);
          diagind = diagind + (diagind-2+batch(1))*length(batch);
          dist(diagind) = +inf;
        end

        % big or small as requested
        if ('n'==params.Subset(1))
          % nearest
          [val,j] = min(dist,[],2);
        else
          % farthest
          [val,j] = max(dist,[],2);
        end

        % and stuff them into the result structure
        d.columnindex(batch) = j;
        d.distance(batch) = val;

        % increment the batch
        batch = batch + bs;
        if batch(end)>n1
          batch(batch>n1) = [];
        end

      end

      % did we need to return a struct or an array?
      if params.Result(1) == 'a'
        % an array. make it a sparse one
        d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
      end

    end % if dim == 1

end  % switch params.Subset

% End of mainline

% ======================================================
% begin subfunctions
% ======================================================
function d = distcomp(set1,set2,params)
% Subfunction to compute all distances between two sets of points
dim = size(set1,2);
% can we take advantage of bsxfun?
% Note: in theory, there is no need to loop over the dimensions. We
% could Just let bsxfun do ALL the work, then wrap a sum around the
% outside. In practice, this tends to create large intermediate
% arrays, especially in higher numbers of dimensions. Its also when
% we might gain here by use of a vectorized code. This will only be
% a serious gain when the number of points is relatively small and
% the dimension is large.
if params.usebsxfun
  % its a recent enough version of matlab that we can
  % use bsxfun at all.
  n1 = size(set1,1);
  n2 = size(set2,1);
  if (dim>1) && ((n1*n2*dim)<=params.ChunkSize)
    % its a small enough problem that we might gain by full
    % use of bsxfun
    switch params.Metric
      case 2
        d = sum(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim])).^2,3);
      case 1
        d = sum(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),3);
      case inf
        d = max(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),[],3);
      case 0
        d = min(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),[],3);
      case 5
          d = acosd(cosd(repmat(set1(:,2),1,n2)) .* cosd(repmat(set2(:,2),1,n1))' .* cosd(repmat(set1(:,1),1,n2) - repmat(set2(:,1),1,n1)') + sind(repmat(set1(:,2),1,n2)) .* sind(repmat(set2(:,2),1,n1)'));
    end
  else
    % too big, so that the ChunkSize will have been exceeded, or just 1-d
    if params.Metric == 2
      d = bsxfun(@minus,set1(:,1),set2(:,1)').^2;
    else
      d = abs(bsxfun(@minus,set1(:,1),set2(:,1)'));
    end
    for i=2:dim
      switch params.Metric
        case 2
          d = d + bsxfun(@minus,set1(:,i),set2(:,i)').^2;
        case 1
          d = d + abs(bsxfun(@minus,set1(:,i),set2(:,i)'));
        case inf
          d = max(d,abs(bsxfun(@minus,set1(:,i),set2(:,i)')));
        case 0
          d = min(d,abs(bsxfun(@minus,set1(:,i),set2(:,i)')));
      end
    end
  end
else
  % Cannot use bsxfun. Sigh. Do things the hard (and slower) way.
  n1 = size(set1,1);
  n2 = size(set2,1);
  if params.Metric == 2
    % Note: While some people might use a different Euclidean
    % norm computation based on expanding the square of the
    % difference of two numbers, that computation is inherantly
    % inaccurate when implemented in floating point arithmetic.
    % While it might be faster, I won't use it here. Sorry.
    d = (repmat(set1(:,1),1,n2) - repmat(set2(:,1)',n1,1)).^2;
  else
    d = abs(repmat(set1(:,1),1,n2) - repmat(set2(:,1)',n1,1));
  end
  for i=2:dim
    switch params.Metric
      case 2
        d = d + (repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)).^2;
      case 1
        d = d + abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1));
      case inf
        d = max(d,abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)));
      case 0
        d = min(d,abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)));
    end
  end
end
% if 2 norm, then we must sqrt at the end
if params.Metric==2
  d = sqrt(d);
end

% ==============================================================
%    end main ipdm
%    begin included function - parse_pv_pairs
% ==============================================================
function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we
% have four parameters that we wish to use optionally in
% the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be
% supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params =
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity'
% was truncated as supplied. Also note that the order the pairs were
% supplied was arbitrary.

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
  error 'Property/value pairs must come in PAIRS.'
end
if n<=0
  % just return the defaults
  return
end

if ~isstruct(params)
  error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
  p_i = lower(pv_pairs{2*i-1});
  v_i = pv_pairs{2*i};

  ind = strmatch(p_i,lpropnames,'exact');
  if isempty(ind)
    ind = find(strncmp(p_i,lpropnames,length(p_i)));
    if isempty(ind)
      error(['No matching property found for: ',pv_pairs{2*i-1}])
    elseif length(ind)>1
      error(['Ambiguous property name: ',pv_pairs{2*i-1}])
    end
  end
  p_i = propnames{ind};

  % override the corresponding default in params.
  % Use setfield for comptability issues with older releases.
  params = setfield(params,p_i,v_i); %#ok

end

