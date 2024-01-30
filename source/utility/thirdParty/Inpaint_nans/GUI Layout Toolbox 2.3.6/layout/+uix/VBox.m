classdef VBox < uix.Box
    %uix.VBox  Vertical box
    %
    %  b = uix.VBox(p1,v1,p2,v2,...) constructs a vertical box and sets
    %  parameter p1 to value v1, etc.
    %
    %  A vertical box lays out contents from top to bottom.
    %
    %  See also: uix.HBox, uix.Grid, uix.VButtonBox, uix.VBoxFlex
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    properties( Access = public, Dependent, AbortSet )
        Heights % heights of contents, in pixels and/or weights
        MinimumHeights % minimum heights of contents, in pixels
    end
    
    properties( Access = protected )
        Heights_ = zeros( [0 1] ) % backing for Heights
        MinimumHeights_ = zeros( [0 1] ) % backing for MinimumHeights
    end
    
    methods
        
        function obj = VBox( varargin )
            %uix.VBox  Vertical box constructor
            %
            %  b = uix.VBox() constructs a horizontal box.
            %
            %  b = uix.VBox(p1,v1,p2,v2,...) sets parameter p1 to value v1,
            %  etc.
            
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
        
        function value = get.Heights( obj )
            
            value = obj.Heights_;
            
        end % get.Heights
        
        function set.Heights( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''Heights'' must be of type double.' )
            assert( all( isreal( value ) ) && ~any( isinf( value ) ) && ...
                ~any( isnan( value ) ), 'uix:InvalidPropertyValue', ...
                'Elements of property ''Heights'' must be real and finite.' )
            assert( isequal( size( value ), size( obj.Contents_ ) ), ...
                'uix:InvalidPropertyValue', ...
                'Size of property ''Heights'' must match size of contents.' )
            
            % Set
            obj.Heights_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.Heights
        
        function value = get.MinimumHeights( obj )
            
            value = obj.MinimumHeights_;
            
        end % get.MinimumHeights
        
        function set.MinimumHeights( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''MinimumHeights'' must be of type double.' )
            assert( all( isreal( value ) ) && ~any( isinf( value ) ) && ...
                all( value >= 0 ), 'uix:InvalidPropertyValue', ...
                'Elements of property ''MinimumHeights'' must be non-negative.' )
            assert( isequal( size( value ), size( obj.Heights_ ) ), ...
                'uix:InvalidPropertyValue', ...
                'Size of property ''MinimumHeights'' must match size of contents.' )
            
            % Set
            obj.MinimumHeights_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.MinimumHeights
        
    end % accessors
    
    methods( Access = protected )
        
        function redraw( obj )
            %redraw  Redraw
            %
            %  c.redraw() redraws the container c.
            
            % Compute positions
            bounds = hgconvertunits( ancestor( obj, 'figure' ), ...
                [0 0 1 1], 'normalized', 'pixels', obj );
            heights = obj.Heights_;
            minimumHeights = obj.MinimumHeights_;
            padding = obj.Padding_;
            spacing = obj.Spacing_;
            r = numel( heights );
            xPositions = [padding + 1, max( bounds(3) - 2 * padding, 1 )];
            xPositions = repmat( xPositions, [r 1] );
            ySizes = uix.calcPixelSizes( bounds(4), heights, ...
                minimumHeights, padding, spacing );
            yPositions = [bounds(4) - cumsum( ySizes ) - padding - ...
                spacing * transpose( 0:r-1 ) + 1, ySizes];
            positions = [xPositions(:,1), yPositions(:,1), ...
                xPositions(:,2), yPositions(:,2)];
            
            % Set positions
            children = obj.Contents_;
            for ii = 1:numel( children )
                uix.setPosition( children(ii), positions(ii,:), 'pixels' )
            end
            
        end % redraw
        
        function addChild( obj, child )
            %addChild  Add child
            %
            %  c.addChild(d) adds the child d to the container c.
            
            % Add to sizes
            obj.Heights_(end+1,:) = -1;
            obj.MinimumHeights_(end+1,:) = 1;
            
            % Call superclass method
            addChild@uix.Box( obj, child )
            
        end % addChild
        
        function removeChild( obj, child )
            %removeChild  Remove child
            %
            %  c.removeChild(d) removes the child d from the container c.
            
            % Remove from sizes
            tf = obj.Contents_ == child;
            obj.Heights_(tf,:) = [];
            obj.MinimumHeights_(tf,:) = [];
            
            % Call superclass method
            removeChild@uix.Box( obj, child )
            
        end % removeChild
        
        function reorder( obj, indices )
            %reorder  Reorder contents
            %
            %  c.reorder(i) reorders the container contents using indices
            %  i, c.Contents = c.Contents(i).
            
            % Reorder
            obj.Heights_ = obj.Heights_(indices,:);
            obj.MinimumHeights_ = obj.MinimumHeights_(indices,:);
            
            % Call superclass method
            reorder@uix.Box( obj, indices )
            
        end % reorder
        
    end % template methods
    
end % classdef