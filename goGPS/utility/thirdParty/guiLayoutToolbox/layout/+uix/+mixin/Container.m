classdef Container < handle
    %uix.mixin.Container  Container mixin
    %
    %  uix.mixin.Container is a mixin class used by containers to provide
    %  various properties and template methods.
    
    %  Copyright 2009-2016 The MathWorks, Inc.
    %  $Revision: 1778 $ $Date: 2019-01-23 23:08:35 +0000 (Wed, 23 Jan 2019) $
    
    properties( Dependent, Access = public )
        Contents % contents in layout order
    end
    
    properties( Access = public, Dependent, AbortSet )
        Padding % space around contents, in pixels
    end
    
    properties( Access = protected )
        Contents_ = gobjects( [0 1] ) % backing for Contents
        Padding_ = 0 % backing for Padding
    end
    
    properties( Dependent, Access = protected )
        Dirty % needs redraw
    end
    
    properties( Access = private )
        Dirty_ = false % backing for Dirty
        FigureObserver % observer
        FigureListener % listener
        ChildObserver % observer
        ChildAddedListener % listener
        ChildRemovedListener % listener
        SizeChangedListener % listener
        ActivePositionPropertyListeners = cell( [0 1] ) % listeners
    end
    
    methods
        
        function obj = Container()
            %uix.mixin.Container  Initialize
            %
            %  c@uix.mixin.Container() initializes the container c.
            
            % Create observers and listeners
            figureObserver = uix.FigureObserver( obj );
            figureListener = event.listener( figureObserver, ...
                'FigureChanged', @obj.onFigureChanged );
            childObserver = uix.ChildObserver( obj );
            childAddedListener = event.listener( ...
                childObserver, 'ChildAdded', @obj.onChildAdded );
            childRemovedListener = event.listener( ...
                childObserver, 'ChildRemoved', @obj.onChildRemoved );
            sizeChangedListener = event.listener( ...
                obj, 'SizeChanged', @obj.onSizeChanged );
            
            % Store observers and listeners
            obj.FigureObserver = figureObserver;
            obj.FigureListener = figureListener;
            obj.ChildObserver = childObserver;
            obj.ChildAddedListener = childAddedListener;
            obj.ChildRemovedListener = childRemovedListener;
            obj.SizeChangedListener = sizeChangedListener;
            
            % Track usage
            obj.track()
            
        end % constructor
        
    end % structors
    
    methods
        
        function value = get.Contents( obj )
            
            value = obj.Contents_;
            
        end % get.Contents
        
        function set.Contents( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Check
            [tf, indices] = ismember( value, obj.Contents_ );
            assert( isequal( size( obj.Contents_ ), size( value ) ) && ...
                numel( value ) == numel( unique( value ) ) && all( tf ), ...
                'uix:InvalidOperation', ...
                'Property ''Contents'' may only be set to a permutation of itself.' )
            
            % Call reorder
            obj.reorder( indices )
            
        end % set.Contents
        
        function value = get.Padding( obj )
            
            value = obj.Padding_;
            
        end % get.Padding
        
        function set.Padding( obj, value )
            
            % Check
            assert( isa( value, 'double' ) && isscalar( value ) && ...
                isreal( value ) && ~isinf( value ) && ...
                ~isnan( value ) && value >= 0, ...
                'uix:InvalidPropertyValue', ...
                'Property ''Padding'' must be a non-negative scalar.' )
            
            % Set
            obj.Padding_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.Padding
        
        function value = get.Dirty( obj )
            
            value = obj.Dirty_;
            
        end % get.Dirty
        
        function set.Dirty( obj, value )
            
            if value
                if obj.isDrawable() % drawable
                    obj.redraw() % redraw now
                else % not drawable
                    obj.Dirty_ = true; % flag for future redraw
                end
            end
            
        end % set.Dirty
        
    end % accessors
    
    methods( Access = private, Sealed )
        
        function onFigureChanged( obj, ~, eventData )
            %onFigureChanged  Event handler
            
            % Call template method
            obj.reparent( eventData.OldFigure, eventData.NewFigure )
            
            % Redraw if possible and if dirty
            if obj.Dirty_ && obj.isDrawable()
                obj.redraw()
                obj.Dirty_ = false;
            end
            
        end % onFigureChanged
        
        function onChildAdded( obj, ~, eventData )
            %onChildAdded  Event handler
            
            % Do nothing if add is internal tree surgery
            if isTreeSurgery( eventData.Child ), return, end
            
            % Call template method
            obj.addChild( eventData.Child )
            
        end % onChildAdded
        
        function onChildRemoved( obj, ~, eventData )
            %onChildRemoved  Event handler
            
            % Do nothing if container is being deleted
            if strcmp( obj.BeingDeleted, 'on' ), return, end
            
            % Do nothing if remove is internal tree surgery
            if isTreeSurgery( eventData.Child ), return, end
            
            % Call template method
            obj.removeChild( eventData.Child )
            
        end % onChildRemoved
        
        function onSizeChanged( obj, ~, ~ )
            %onSizeChanged  Event handler
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % onSizeChanged
        
        function onActivePositionPropertyChanged( obj, ~, ~ )
            %onActivePositionPropertyChanged  Event handler
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % onActivePositionPropertyChanged
        
    end % event handlers
    
    methods( Abstract, Access = protected )
        
        redraw( obj )
        
    end % abstract template methods
    
    methods( Access = protected )
        
        function addChild( obj, child )
            %addChild  Add child
            %
            %  c.addChild(d) adds the child d to the container c.
            
            % Add to contents
            obj.Contents_(end+1,:) = child;
            
            % Add listeners
            if isa( child, 'matlab.graphics.axis.Axes' )
                obj.ActivePositionPropertyListeners{end+1,:} = ...
                    event.proplistener( child, ...
                    findprop( child, 'ActivePositionProperty' ), ...
                    'PostSet', @obj.onActivePositionPropertyChanged );
            else
                obj.ActivePositionPropertyListeners{end+1,:} = [];
            end
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % addChild
        
        function removeChild( obj, child )
            %removeChild  Remove child
            %
            %  c.removeChild(d) removes the child d from the container c.
            
            % Remove from contents
            contents = obj.Contents_;
            tf = contents == child;
            obj.Contents_(tf,:) = [];
            
            % Remove listeners
            obj.ActivePositionPropertyListeners(tf,:) = [];
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % removeChild
        
        function reparent( obj, oldFigure, newFigure ) %#ok<INUSD>
            %reparent  Reparent container
            %
            %  c.reparent(a,b) reparents the container c from the figure a
            %  to the figure b.
            
        end % reparent
        
        function reorder( obj, indices )
            %reorder  Reorder contents
            %
            %  c.reorder(i) reorders the container contents using indices
            %  i, c.Contents = c.Contents(i).
            
            % Reorder contents
            obj.Contents_ = obj.Contents_(indices,:);
            
            % Reorder listeners
            obj.ActivePositionPropertyListeners = ...
                obj.ActivePositionPropertyListeners(indices,:);
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % reorder
        
        function tf = isDrawable( obj )
            %isDrawable  Test for drawability
            %
            %  c.isDrawable() is true if the container c is drawable, and
            %  false otherwise.  To be drawable, a container must be
            %  rooted.
            
            tf = ~isempty( obj.FigureObserver.Figure );
            
        end % isDrawable
        
        function track( obj )
            %track  Track usage
            
            persistent TRACKED % single shot
            if isempty( TRACKED )
                v = ver( 'layout' );
                try %#ok<TRYNC>
                    uix.tracking( 'UA-82270656-2', v(1).Version, class( obj ) )
                end
                TRACKED = true;
            end
            
        end % track
        
    end % template methods
    
end % classdef

function tf = isTreeSurgery( child )
%isTreeSurgery  Test for internal tree surgery
%
%  Certain internal operations perform tree surgery, removing a child
%  temporarily and adding it back under an additional node.

if ~isa( child, 'matlab.graphics.axis.Axes' ) % only certain types
    tf = false;
elseif verLessThan( 'MATLAB', '9.5' ) % only certain versions
    tf = false;
else % check stack
    s = dbstack();
    tf = any( strcmp( 'ScribeStackManager.createLayer', {s.name} ) );
end

end % isTreeSurgery