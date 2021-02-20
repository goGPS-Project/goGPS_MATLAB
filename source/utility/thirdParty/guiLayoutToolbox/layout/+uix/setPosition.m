function setPosition( o, p, u )
%setPosition  Set position of graphics object
%
%  setPosition(o,p,u) sets the position of a graphics object o to value p
%  with units u.
%
%  In contrast to setting the Position property directly, this function
%  honors the ActivePositionProperty of axes.

%  Copyright 2009-2016 The MathWorks, Inc.
%  $Revision: 1765 $ $Date: 2018-10-25 17:11:26 +0100 (Thu, 25 Oct 2018) $

% Set units
if ~isempty( findprop( o, 'Units' ) )
    o.Units = u;
else
    assert( strcmp( u, 'pixels' ), 'uix:InvalidOperation', ...
        'Objects without property ''Units'' have units of ''pixels''.' )
end

% Set position
if ~isempty( findprop( o, 'ActivePositionProperty' ) )
    switch o.ActivePositionProperty
        case 'position'
            o.Position = p;
        case 'outerposition'
            o.OuterPosition = p;
        otherwise
            error( 'uix:InvalidState', ...
                'Unknown value ''%s'' for property ''ActivePositionProperty'' of %s.', ...
                o.ActivePositionProperty, class( o ) )
    end
else
    o.Position = p;
end

end % setPosition