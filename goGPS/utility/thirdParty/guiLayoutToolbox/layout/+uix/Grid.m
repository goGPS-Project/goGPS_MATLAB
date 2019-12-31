classdef Grid < uix.Box
    %uix.Grid  Grid
    %
    %  b = uix.Grid(p1,v1,p2,v2,...) constructs a grid and sets parameter
    %  p1 to value v1, etc.
    %
    %  A grid lays out contents from top to bottom and left to right.
    %
    %  See also: uix.HBox, uix.VBox, uix.GridFlex
    
    %  Copyright 2009-2016 The MathWorks, Inc.
    %  $Revision: 1594 $ $Date: 2018-03-28 02:27:52 +1100 (Wed, 28 Mar 2018) $
    
    properties( Access = public, Dependent, AbortSet )
        Widths % widths of contents, in pixels and/or weights
        MinimumWidths % minimum widths of contents, in pixels
        Heights % heights of contents, in pixels and/or weights
        MinimumHeights % minimum heights of contents, in pixels
    end
    
    properties( Access = protected )
        Widths_ = zeros( [0 1] ) % backing for Widths
        MinimumWidths_ = zeros( [0 1] ) % backing for MinimumWidths
        Heights_ = zeros( [0 1] ) % backing for Heights
        MinimumHeights_ = zeros( [0 1] ) % backing for MinimumHeights
    end
    
    methods
        
        function obj = Grid( varargin )
            %uix.Grid  Grid constructor
            %
            %  b = uix.Grid() constructs a grid.
            %
            %  b = uix.Grid(p1,v1,p2,v2,...) sets parameter p1 to value v1,
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
            n = numel( obj.Contents_ );
            b = numel( obj.Widths_ );
            q = numel( obj.Heights_ );
            c = numel( value );
            r = ceil( n / c );
            if c < min( [1 n] )
                error( 'uix:InvalidPropertyValue' , ...
                    'Property ''Widths'' must be non-empty for non-empty contents.' )
            elseif ceil( n / r ) < c
                error( 'uix:InvalidPropertyValue' , ...
                    'Size of property ''Widths'' must not lead to empty columns.' )
            elseif c > n
                error( 'uix:InvalidPropertyValue' , ...
                    'Size of property ''Widths'' must be no larger than size of contents.' )
            end
            
            % Set
            obj.Widths_ = value;
            if c < b % number of columns decreasing
                obj.MinimumWidths_(c+1:end,:) = [];
                if r > q % number of rows increasing
                    obj.Heights_(end+1:r,:) = -1;
                    obj.MinimumHeights_(end+1:r,:) = 1;
                end
            elseif c > b % number of columns increasing
                obj.MinimumWidths_(end+1:c,:) = -1;
                if r < q % number of rows decreasing
                    obj.Heights_(r+1:end,:) = [];
                    obj.MinimumHeights_(r+1:end,:) = [];
                end
            end
            
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
            n = numel( obj.Contents_ );
            b = numel( obj.Widths_ );
            q = numel( obj.Heights_ );
            r = numel( value );
            c = ceil( n / r );
            if r < min( [1 n] )
                error( 'uix:InvalidPropertyValue' , ...
                    'Property ''Heights'' must be non-empty for non-empty contents.' )
            elseif r > n
                error( 'uix:InvalidPropertyValue' , ...
                    'Size of property ''Heights'' must be no larger than size of contents.' )
            end
            
            % Set
            obj.Heights_ = value;
            if r < q % number of rows decreasing
                obj.MinimumHeights_(r+1:end,:) = [];
                if c > b % number of columns increasing
                    obj.Widths_(end+1:c,:) = -1;
                    obj.MinimumWidths_(end+1:c,:) = 1;
                end
            elseif r > q % number of rows increasing
                obj.MinimumHeights_(end+1:r,:) = 1;
                if c < b % number of columns decreasing
                    obj.Widths_(c+1:end,:) = [];
                    obj.MinimumWidths_(c+1:end,:) = [];
                end
            end
            
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
            widths = obj.Widths_;
            minimumWidths = obj.MinimumWidths_;
            heights = obj.Heights_;
            minimumHeights = obj.MinimumHeights_;
            padding = obj.Padding_;
            spacing = obj.Spacing_;
            c = numel( widths );
            r = numel( heights );
            n = numel( obj.Contents_ );
            xSizes = uix.calcPixelSizes( bounds(3), widths, ...
                minimumWidths, padding, spacing );
            xPositions = [cumsum( [0; xSizes(1:end-1,:)] ) + padding + ...
                spacing * transpose( 0:c-1 ) + 1, xSizes];
            ySizes = uix.calcPixelSizes( bounds(4), heights, ...
                minimumHeights, padding, spacing );
            yPositions = [bounds(4) - cumsum( ySizes ) - padding - ...
                spacing * transpose( 0:r-1 ) + 1, ySizes];
            [iy, ix] = ind2sub( [r c], transpose( 1:n ) );
            positions = [xPositions(ix,1), yPositions(iy,1), ...
                xPositions(ix,2), yPositions(iy,2)];
            
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
            
            % Add column and even a row if necessary
            n = numel( obj.Contents_ );
            c = numel( obj.Widths_ );
            r = numel( obj.Heights_ );
            if n == 0
                obj.Widths_(end+1,:) = -1;
                obj.MinimumWidths_(end+1,:) = 1;
                obj.Heights_(end+1,:) = -1;
                obj.MinimumHeights_(end+1,:) = 1;
            elseif ceil( (n+1)/r ) > c
                obj.Widths_(end+1,:) = -1;
                obj.MinimumWidths_(end+1,:) = 1;
            end
            
            % Call superclass method
            addChild@uix.Box( obj, child )
            
        end % addChild
        
        function removeChild( obj, child )
            %removeChild  Remove child
            %
            %  c.removeChild(d) removes the child d from the container c.
            
            % Remove column and even row if necessary
            n = numel( obj.Contents_ );
            c = numel( obj.Widths_ );
            r = numel( obj.Heights_ );
            if n == 1
                obj.Widths_(end,:) = [];
                obj.MinimumWidths_(end,:) = [];
                obj.Heights_(end,:) = [];
                obj.MinimumHeights_(end,:) = [];
            elseif c == 1
                obj.Heights_(end,:) = [];
                obj.MinimumHeights_(end,:) = [];
            elseif ceil( (n-1)/r ) < c
                obj.Widths_(end,:) = [];
                obj.MinimumWidths_(end,:) = [];
            end
            
            % Call superclass method
            removeChild@uix.Box( obj, child )
            
        end % removeChild
        
    end % template methods
    
end % classdef