function obj = Empty( varargin )
%uix.Empty  Create an empty space
%
%   obj = uix.Empty() creates an empty space that can be used to add gaps
%   between elements in layouts.
%
%   obj = uix.Empty(param,value,...) also sets one or more property
%   values.
%
%   See the <a href="matlab:doc uix.Empty">documentation</a> for more detail and the list of properties.
%
%   Examples:
%   >> f = figure();
%   >> box = uix.HBox( 'Parent', f );
%   >> uicontrol( 'Parent', box, 'Background', 'r' )
%   >> uix.Empty( 'Parent', box )
%   >> uicontrol( 'Parent', box, 'Background', 'b' )

%   Copyright 2009-2016 The MathWorks, Inc.
%   $Revision: 1436 $ $Date: 2016-11-17 17:53:29 +0000 (Thu, 17 Nov 2016) $

% Create uicontainer
obj = matlab.ui.container.internal.UIContainer( 'Tag', 'empty', varargin{:} );

% Create property for Parent listener
p = addprop( obj, 'ParentListener' );
p.Hidden = true;

% Create Parent listener
obj.ParentListener = event.proplistener( obj, ...
    findprop( obj, 'Parent' ), 'PostSet', @(~,~)onParentChanged(obj) );

% Create property for Parent color listener
p = addprop( obj, 'ParentColorListener' );
p.Hidden = true;

% Initialize color and listener
updateColor( obj )
updateListener( obj )

end % uix.Empty

function onParentChanged( obj )
%onParentColorChanged  Event handler

% Update color and listener
updateColor( obj )
updateListener( obj )

end % onParentChanged

function onParentColorChanged( obj )
%onParentColorChanged  Event handler

% Update color
updateColor( obj )

end % onParentColorChanged

function name = getColorProperty( obj )
%getColorProperty  Get color property

names = {'Color','BackgroundColor'}; % possible names
for ii = 1:numel( names ) % loop over possible names
    name = names{ii};
    if isprop( obj, name )
        return
    end
end
error( 'Cannot find color property for %s.', class( obj ) )

end % getColorProperty

function updateColor( obj )
%updateColor  Set uicontainer BackgroundColor to match Parent

parent = obj.Parent;
if isempty( parent ), return, end
property = getColorProperty( parent );
color = parent.( property );
try
    obj.BackgroundColor = color;
catch e
    warning( e.identifier, e.message ) % rethrow as warning
end

end % updateColor

function updateListener( obj )
%updateListener  Create listener to parent color property

parent = obj.Parent;
if isempty( parent )
    obj.ParentColorListener = [];
else
    property = getColorProperty( parent );
    obj.ParentColorListener = event.proplistener( parent, ...
        findprop( parent, property ), 'PostSet', ...
        @(~,~)onParentColorChanged(obj) );
end

end % updateListener