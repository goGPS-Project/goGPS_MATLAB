function obj = Empty( varargin )
%uiextras.Empty  Create an empty space
%
%   obj = uiextras.Empty() creates an empty space that can be used to add
%   gaps between elements in layouts.
%
%   obj = uiextras.Empty(param,value,...) also sets one or more property
%   values.
%
%   See the <a href="matlab:doc uiextras.Empty">documentation</a> for more detail and the list of properties.
%
%   Examples:
%   >> f = figure();
%   >> box = uiextras.HBox( 'Parent', f );
%   >> uicontrol( 'Parent', box, 'Background', 'r' )
%   >> uiextras.Empty( 'Parent', box )
%   >> uicontrol( 'Parent', box, 'Background', 'b' )

%  Copyright 2009-2020 The MathWorks, Inc.

% Call uix construction function
obj = uix.Empty( varargin{:} );

% Auto-parent
if ~ismember( 'Parent', varargin(1:2:end) )
    obj.Parent = gcf();
end

end % uiextras.Empty