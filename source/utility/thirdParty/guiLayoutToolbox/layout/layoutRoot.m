function folder = layoutRoot()
%layoutRoot  Folder containing the GUI Layout Toolbox
%
%   folder = layoutRoot() returns the full path to the folder containing
%   the GUI Layout Toolbox.
%
%   Examples:
%   >> folder = layoutRoot()
%   folder = 'C:\tools\layouts2\layout'
%
%   See also: layoutVersion

%  Copyright 2009-2014 The MathWorks, Inc.
%  $Revision: 980 $ $Date: 2014-09-28 14:27:26 -0400 (Sun, 28 Sep 2014) $

folder = fileparts( mfilename( 'fullpath' ) );

end % layoutRoot