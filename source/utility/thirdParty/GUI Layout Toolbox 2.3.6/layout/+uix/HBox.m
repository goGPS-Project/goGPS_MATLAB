classdef HBox < uix.Box
    %uix.HBox  Horizontal box
    %
    %  b = uix.HBox(p1,v1,p2,v2,...) constructs a horizontal box and sets
    %  parameter p1 to value v1, etc.
    %
    %  A horizontal box lays out contents from left to right.
    %
    %  See also: uix.VBox, uix.Grid, uix.HButtonBox, uix.HBoxFlex
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    properties( Access = public, Dependent, AbortSet )
        Widths % widths of contents, in pixels and/or weights
        MinimumWidths % minimum widths of contents, in pixels
    end
    
    properties( Access = protected )
        Widths_ = zeros( [0 1] ) % backing for Widths
        MinimumWidths_ = zeros( [0 1] ) % backing for MinimumWidths
    end
    
    methods
        
        function obj = HBox( varargin )
            %uix.HBox  Horizontal box constructor
            %
            %  b = uix.HBox() constructs a horizontal box.
            %
            %  b = uix.HBox(p1,v1,p2,v2,...) sets parameter p1 to value v1,
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
        
        function value = get.Widths( obj )
            
            value = obj.Widths_;
            
        end % get.Widths
        
        function set.Widths( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''Widths'' must be of type double.' )
            assert( all( isreal( value ) ) && ~any( isinf( value ) ) && ...
                ~any( isnan( value ) ), 'uix:InvalidPropertyValue', ...
                'Elements of property ''Widths'' must be real and finite.' )
            assert( isequal( size( value ), size( obj.Contents_ ) ), ...
                'uix:InvalidPropertyValue', ...
                'Size of property ''Widths'' must match size of contents.' )
            
            % Set
            obj.Widths_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.Widths
        
        function value = get.MinimumWidths( obj )
            
            value = obj.MinimumWidths_;
            
        end % get.MinimumWidths
        
        function set.MinimumWidths( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''MinimumWidths'' must be of type double.' )
            assert( all( isreal( value ) ) && ~any( isinf( value ) ) && ...
                all( value >= 0 ), 'uix:InvalidPropertyValue', ...
                'Elements of property ''MinimumWidths'' must be non-negative.' )
            assert( isequal( size( value ), size( obj.Widths_ ) ), ...
                'uix:InvalidPropertyValue', ...
                'Size of property ''MinimumWidths'' must match size of contents.' )
            
            % Set
            obj.MinimumWidths_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.MinimumWidths
        
    end % accessors
    
    methods( Access = protected )
        
        function redraw( obj )
            %redraw  Redraw
            %
            %  c.redraw() redraws the container c.
            
            % Compute positions
            bounds = hgconvertunits( ancestor( obj, 'figure' ), ...
                [0 0 1 1], 'normalized', 'pixels', obj );
            widths = obj.Widths_;
            minimumWidths = obj.MinimumWidths_;
            padding = obj.Padding_;
            spacing = obj.Spacing_;
            c = numel( widths );
            xSizes = uix.calcPixelSizes( bounds(3), widths, ...
                minimumWidths, padding, spacing );
            xPositions = [cumsum( [0; xSizes(1:c-1,:)] ) + padding + ...
                spacing * transpose( 0:c-1 ) + 1, xSizes];
            yPositions = [padding + 1, max( bounds(4) - 2 * padding, 1 )];
            yPositions = repmat( yPositions, [c 1] );
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
            obj.Widths_(end+1,:) = -1;
            obj.MinimumWidths_(end+1,:) = 1;
            
            % Call superclass method
            addChild@uix.Box( obj, child )
            
        end % addChild
        
        function removeChild( obj, child )
            %removeChild  Remove child
            %
            %  c.removeChild(d) removes the child d from the container c.
            
            % Remove from sizes
            tf = obj.Contents_ == child;
            obj.Widths_(tf,:) = [];
            obj.MinimumWidths_(tf,:) = [];
            
            % Call superclass method
            removeChild@uix.Box( obj, child )
            
        end % removeChild
        
        function reorder( obj, indices )
            %reorder  Reorder contents
            %
            %  c.reorder(i) reorders the container contents using indices
            %  i, c.Contents = c.Contents(i).
            
            % Reorder
            obj.Widths_ = obj.Widths_(indices,:);
            obj.MinimumWidths_ = obj.MinimumWidths_(indices,:);
            
            % Call superclass method
            reorder@uix.Box( obj, indices )
            
        end % reorder
        
    end % template methods
    
end % classdef