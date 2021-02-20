function set( obj, varargin )
%uix.set  Set property values
%
%  uix.set(o,p1,v1,p2,v2,...) sets property p1 of the object o to value v1,
%  property p2 to value v2, etc.
%
%  In contrast to builtin set, querying possible values is not supported.

%  Copyright 2017 The MathWorks, Inc.
%  $Revision: 1435 $ $Date: 2016-11-17 17:50:34 +0000 (Thu, 17 Nov 2016) $

if nargin == 1, return, end
assert( rem( nargin, 2 ) == 1, 'uix:InvalidArgument', ...
    'Parameters and values must be provided in pairs.' )
set( obj, varargin{:} )

end % set