function p = getPosition( o, u )
%getPosition  Get position of graphics object
%
%  p = getPosition(o,u) sets the position of a graphics object o with units
%  u.
%
%  In contrast to getting the Position property directly, this function
%  honors the ActivePositionProperty of axes.

%  Copyright 2009-2016 The MathWorks, Inc.
%  $Revision: 1435 $ $Date: 2016-11-17 17:50:34 +0000 (Thu, 17 Nov 2016) $

% Get position
if ~isempty( findprop( o, 'ActivePositionProperty' ) )
    switch o.ActivePositionProperty
        case 'position'
            q = o.Position;
        case 'outerposition'
            q = o.OuterPosition;
        otherwise
            error( 'uix:InvalidState', ...
                'Unknown value ''%s'' for property ''ActivePositionProperty'' of %s.', ...
                o.ActivePositionProperty, class( o ) )
    end
else
    q = o.Position;
end

% Get units
if ~isempty( findprop( o, 'Units' ) )
    v = o.Units;
else
    v = 'pixels';
end

% Convert
if strcmp( u, v ) % trivial
    p = q;
else % conversion required
    p = hgconvertunits( ancestor( o, 'figure' ), q, v, u, o.Parent );
end

end % getPosition