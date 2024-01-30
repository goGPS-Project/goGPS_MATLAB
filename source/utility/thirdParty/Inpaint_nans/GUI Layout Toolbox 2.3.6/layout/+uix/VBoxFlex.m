classdef VBoxFlex < uix.VBox & uix.mixin.Flex
    %uix.VBoxFlex  Flexible vertical box
    %
    %  b = uix.VBoxFlex(p1,v1,p2,v2,...) constructs a flexible vertical box
    %  and sets parameter p1 to value v1, etc.
    %
    %  A vertical box lays out contents from top to bottom.  Users can
    %  resize contents by dragging the dividers.
    %
    %  See also: uix.HBoxFlex, uix.GridFlex, uix.VBox, uix.VButtonBox
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    properties( Access = public, Dependent, AbortSet )
        DividerMarkings % divider markings [on|off]
    end
    
    properties( Access = private )
        RowDividers = uix.Divider.empty( [0 1] ) % row dividers
        FrontDivider % front divider
        DividerMarkings_ = 'on' % backing for DividerMarkings
        MousePressListener = event.listener.empty( [0 0] ) % mouse press listener
        MouseReleaseListener = event.listener.empty( [0 0] ) % mouse release listener
        MouseMotionListener = event.listener.empty( [0 0] ) % mouse motion listener
        ActiveDivider = 0 % active divider index
        ActiveDividerPosition = [NaN NaN NaN NaN] % active divider position
        MousePressLocation = [NaN NaN] % mouse press location
        BackgroundColorListener % background color listener
    end
    
    methods
        
        function obj = VBoxFlex( varargin )
            %uix.VBoxFlex  Flexible vertical box constructor
            %
            %  b = uix.VBoxFlex() constructs a flexible vertical box.
            %
            %  b = uix.VBoxFlex(p1,v1,p2,v2,...) sets parameter p1 to value
            %  v1, etc.
            
            % Create front divider
            frontDivider = uix.Divider( 'Parent', obj, ...
                'Orientation', 'horizontal', ...
                'BackgroundColor', obj.BackgroundColor * 0.75, ...
                'Visible', 'off' );
            
            % Create listeners
            backgroundColorListener = event.proplistener( obj, ...
                findprop( obj, 'BackgroundColor' ), 'PostSet', ...
                @obj.onBackgroundColorChanged );
            
            % Store properties
            obj.FrontDivider = frontDivider;
            obj.BackgroundColorListener = backgroundColorListener;
            
            % Set Spacing property to 5 (may be overwritten by uix.set)
            obj.Spacing = 5;
            
            % Set properties
            try
                uix.set( obj, varargin{:} )
            catch e
                delete( obj )
                e.throwAsCaller()
            end
            
        end % constructor
        
    end % structors
    
    methods
        
        function value = get.DividerMarkings( obj )
            
            value = obj.DividerMarkings_;
            
        end % get.DividerMarkings
        
        function set.DividerMarkings( obj, value )
            
            % Check
            value = uix.validateScalarStringOrCharacterArray( value, ...
                'DividerMarkings' );
            assert( any( strcmp( value, {'on','off'} ) ), ...
                'uix:InvalidArgument', ...
                'Property ''DividerMarkings'' must be ''on'' or ''off'.' )
            
            % Set
            obj.DividerMarkings_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.DividerMarkings
        
    end % accessors
    
    methods( Access = protected )
        
        function onMousePress( obj, source, eventData )
            %onMousePress  Handler for WindowMousePress events
            
            % Check whether mouse is over a divider
            loc = find( obj.RowDividers.isMouseOver( eventData ) );
            if isempty( loc ), return, end
            
            % Capture state at button down
            divider = obj.RowDividers(loc);
            obj.ActiveDivider = loc;
            obj.ActiveDividerPosition = divider.Position;
            root = groot();
            obj.MousePressLocation = root.PointerLocation;
            
            % Make sure the pointer is appropriate
            obj.updateMousePointer( source, eventData );
            
            % Activate divider
            frontDivider = obj.FrontDivider;
            frontDivider.Position = divider.Position;
            divider.Visible = 'off';
            frontDivider.Parent = [];
            frontDivider.Parent = obj;
            frontDivider.Visible = 'on';
            
        end % onMousePress
        
        function onMouseRelease( obj, ~, ~ )
            %onMousePress  Handler for WindowMouseRelease events
            
            % Compute new positions
            loc = obj.ActiveDivider;
            if loc > 0
                root = groot();
                delta = root.PointerLocation(2) - obj.MousePressLocation(2);
                ih = loc;
                jh = loc + 1;
                ic = loc;
                jc = loc + 1;
                divider = obj.RowDividers(loc);
                contents = obj.Contents_;
                ip = uix.getPosition( contents(ic), 'pixels' );
                jp = uix.getPosition( contents(jc), 'pixels' );
                oldPixelHeights = [ip(4); jp(4)];
                minimumHeights = obj.MinimumHeights_(ih:jh,:);
                if delta < 0 % limit to minimum distance from lower neighbor
                    delta = max( delta, minimumHeights(2) - oldPixelHeights(2) );
                else % limit to minimum distance from upper neighbor
                    delta = min( delta, oldPixelHeights(1) - minimumHeights(1) );
                end
                oldHeights = obj.Heights_(loc:loc+1);
                newPixelHeights = oldPixelHeights - delta * [1;-1];
                if oldHeights(1) < 0 && oldHeights(2) < 0 % weight, weight
                    newHeights = oldHeights .* newPixelHeights ./ oldPixelHeights;
                elseif oldHeights(1) < 0 && oldHeights(2) >= 0 % weight, pixels
                    newHeights = [oldHeights(1) * newPixelHeights(1) / ...
                        oldPixelHeights(1); newPixelHeights(2)];
                elseif oldHeights(1) >= 0 && oldHeights(2) < 0 % pixels, weight
                    newHeights = [newPixelHeights(1); oldHeights(2) * ...
                        newPixelHeights(2) / oldPixelHeights(2)];
                else % sizes(1) >= 0 && sizes(2) >= 0 % pixels, pixels
                    newHeights = newPixelHeights;
                end
                obj.Heights_(loc:loc+1) = newHeights;
            else
                return
            end
            
            % Deactivate divider
            obj.FrontDivider.Visible = 'off';
            divider.Visible = 'on';
            
            % Reset state at button down
            obj.ActiveDivider = 0;
            obj.ActiveDividerPosition = [NaN NaN NaN NaN];
            obj.MousePressLocation = [NaN NaN];
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % onMouseRelease
        
        function onMouseMotion( obj, source, eventData )
            %onMouseMotion  Handler for WindowMouseMotion events
            
            loc = obj.ActiveDivider;
            if loc == 0 % hovering, update pointer
                obj.updateMousePointer( source, eventData );
            else % dragging row divider
                root = groot();
                delta = root.PointerLocation(2) - obj.MousePressLocation(2);
                ih = loc;
                jh = loc + 1;
                ic = loc;
                jc = loc + 1;
                contents = obj.Contents_;
                ip = uix.getPosition( contents(ic), 'pixels' );
                jp = uix.getPosition( contents(jc), 'pixels' );
                oldPixelHeights = [ip(4); jp(4)];
                minimumHeights = obj.MinimumHeights_(ih:jh,:);
                if delta < 0 % limit to minimum distance from lower neighbor
                    delta = max( delta, minimumHeights(2) - oldPixelHeights(2) );
                else % limit to minimum distance from upper neighbor
                    delta = min( delta, oldPixelHeights(1) - minimumHeights(1) );
                end
                obj.FrontDivider.Position = ...
                    obj.ActiveDividerPosition + [0 delta 0 0];
            end
            
        end % onMouseMotion
        
        function onBackgroundColorChanged( obj, ~, ~ )
            %onBackgroundColorChanged  Handler for BackgroundColor changes
            
            backgroundColor = obj.BackgroundColor;
            highlightColor = min( [backgroundColor / 0.75; 1 1 1] );
            shadowColor = max( [backgroundColor * 0.75; 0 0 0] );
            rowDividers = obj.RowDividers;
            for ii = 1:numel( rowDividers )
                rowDivider = rowDividers(ii);
                rowDivider.BackgroundColor = backgroundColor;
                rowDivider.HighlightColor = highlightColor;
                rowDivider.ShadowColor = shadowColor;
            end
            frontDivider = obj.FrontDivider;
            frontDivider.BackgroundColor = shadowColor;
            
        end % onBackgroundColorChanged
        
    end % event handlers
    
    methods( Access = protected )
        
        function redraw( obj )
            %redraw  Redraw contents
            %
            %  c.redraw() redraws the container c.
            
            % Call superclass method
            redraw@uix.VBox( obj )
            
            % Create or destroy row dividers
            q = numel( obj.RowDividers ); % current number of dividers
            r = max( [numel( obj.Heights_ )-1 0] ); % required number of dividers
            if q < r % create
                for ii = q+1:r
                    divider = uix.Divider( 'Parent', obj, ...
                        'Orientation', 'horizontal', ...
                        'BackgroundColor', obj.BackgroundColor );
                    obj.RowDividers(ii,:) = divider;
                end
            elseif q > r % destroy
                % Destroy dividers
                delete( obj.RowDividers(r+1:q,:) )
                obj.RowDividers(r+1:q,:) = [];
                % Update pointer
                if r == 0 && strcmp( obj.Pointer, 'top' )
                    obj.unsetPointer()
                end
            end
            
            % Compute container bounds
            bounds = hgconvertunits( ancestor( obj, 'figure' ), ...
                [0 0 1 1], 'normalized', 'pixels', obj );
            
            % Retrieve size properties
            heights = obj.Heights_;
            minimumHeights = obj.MinimumHeights_;
            padding = obj.Padding_;
            spacing = obj.Spacing_;
            
            % Compute row divider positions
            xRowPositions = [padding + 1, max( bounds(3) - 2 * padding, 1 )];
            xRowPositions = repmat( xRowPositions, [r 1] );
            yRowSizes = uix.calcPixelSizes( bounds(4), heights, ...
                minimumHeights, padding, spacing );
            yRowPositions = [bounds(4) - cumsum( yRowSizes(1:r,:) ) - padding - ...
                spacing * transpose( 1:r ) + 1, repmat( spacing, [r 1] )];
            rowPositions = [xRowPositions(:,1), yRowPositions(:,1), ...
                xRowPositions(:,2), yRowPositions(:,2)];
            
            % Position row dividers
            for ii = 1:r
                rowDivider = obj.RowDividers(ii);
                rowDivider.Position = rowPositions(ii,:);
                switch obj.DividerMarkings_
                    case 'on'
                        rowDivider.Markings = rowPositions(ii,3)/2;
                    case 'off'
                        rowDivider.Markings = zeros( [0 1] );
                end
            end
            
        end % redraw
        
        function reparent( obj, oldFigure, newFigure )
            %reparent  Reparent container
            %
            %  c.reparent(a,b) reparents the container c from the figure a
            %  to the figure b.
            
            % Update listeners
            if isempty( newFigure )
                mousePressListener = event.listener.empty( [0 0] );
                mouseReleaseListener = event.listener.empty( [0 0] );
                mouseMotionListener = event.listener.empty( [0 0] );
            else
                mousePressListener = event.listener( newFigure, ...
                    'WindowMousePress', @obj.onMousePress );
                mouseReleaseListener = event.listener( newFigure, ...
                    'WindowMouseRelease', @obj.onMouseRelease );
                mouseMotionListener = event.listener( newFigure, ...
                    'WindowMouseMotion', @obj.onMouseMotion );
            end
            obj.MousePressListener = mousePressListener;
            obj.MouseReleaseListener = mouseReleaseListener;
            obj.MouseMotionListener = mouseMotionListener;
            
            % Call superclass method
            reparent@uix.VBox( obj, oldFigure, newFigure )
            
            % Update pointer
            if ~isempty( oldFigure ) && ~strcmp( obj.Pointer, 'unset' )
                obj.unsetPointer()
            end
            
        end % reparent
        
    end % template methods
    
    methods( Access = protected )
        
        function updateMousePointer ( obj, source, eventData  )
            
            oldPointer = obj.Pointer;
            if any( obj.RowDividers.isMouseOver( eventData ) )
                newPointer = 'top';
            else
                newPointer = 'unset';
            end
            switch newPointer
                case oldPointer % no change
                    % do nothing
                case 'unset' % change, unset
                    obj.unsetPointer()
                otherwise % change, set
                    obj.setPointer( source, newPointer )
            end
            
        end % updateMousePointer
        
    end % helpers methods
    
end % classdef