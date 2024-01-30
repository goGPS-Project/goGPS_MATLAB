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

%  Copyright 2009-2020 The MathWorks, Inc.

folder = fileparts( mfilename( 'fullpath' ) );

end % layoutRoot