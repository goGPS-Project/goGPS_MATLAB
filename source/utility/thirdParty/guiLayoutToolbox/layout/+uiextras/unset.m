function unset( ~, ~, ~ )
%uiextras.unset  Clear a default property value from a parent object
%
%  This functionality has been removed.

%  Copyright 2009-2014 The MathWorks, Inc.
%  $Revision: 979 $ $Date: 2014-09-28 14:26:12 -0400 (Sun, 28 Sep 2014) $

% Check inputs
narginchk( 2, 2 )

% Warn
warning( 'uiextras:Deprecated', 'uiextras.unset has been removed.' )

end % uiextras.unset