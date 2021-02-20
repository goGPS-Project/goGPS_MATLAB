classdef Grid < uix.Grid
    %uiextras.Grid  Container with contents arranged in a grid
    %
    %   obj = uiextras.Grid() creates a new new grid layout with all
    %   properties set to defaults. The number of rows and columns to use
    %   is determined from the number of elements in the RowSizes and
    %   ColumnSizes properties respectively. Child elements are arranged
    %   down column one first, then column two etc. If there are
    %   insufficient columns then a new one is added. The output is a new
    %   layout object that can be used as the parent for other
    %   user-interface components. The output is a new layout object that
    %   can be used as the parent for other user-interface components.
    %
    %   obj = uiextras.Grid(param,value,...) also sets one or more
    %   parameter values.
    %
    %   See the <a href="matlab:doc uiextras.Grid">documentation</a> for more detail and the list of properties.
    %
    %   Examples:
    %   >> f = figure();
    %   >> g = uiextras.Grid( 'Parent', f, 'Spacing', 5 );
    %   >> uicontrol( 'Style', 'frame', 'Parent', g, 'Background', 'r' )
    %   >> uicontrol( 'Style', 'frame', 'Parent', g, 'Background', 'b' )
    %   >> uicontrol( 'Style', 'frame', 'Parent', g, 'Background', 'g' )
    %   >> uiextras.Empty( 'Parent', g )
    %   >> uicontrol( 'Style', 'frame', 'Parent', g, 'Background', 'c' )
    %   >> uicontrol( 'Style', 'frame', 'Parent', g, 'Background', 'y' )
    %   >> set( g, 'ColumnSizes', [-1 100 -2], 'RowSizes', [-1 100] );
    %
    %   See also: uiextras.GridFlex
    
    %  Copyright 2009-2014 The MathWorks, Inc.
    %  $Revision: 1077 $ $Date: 2015-03-19 16:44:14 +0000 (Thu, 19 Mar 2015) $
    
    properties( Hidden, Access = public, Dependent )
        Enable % deprecated
        RowSizes % heights of contents, in pixels and/or weights
        MinimumRowSizes % minimum heights of contents, in pixels
        ColumnSizes % widths of contents, in pixels and/or weights
        MinimumColumnSizes % minimum widths of contents, in pixels
    end
    
    methods
        
        function obj = Grid( varargin )
            
            % Call uix constructor
            obj@uix.Grid( varargin{:} )
            
            % Auto-parent
            if ~ismember( 'Parent', varargin(1:2:end) )
                obj.Parent = gcf();
            end
            
        end % constructor
        
    end % structors
    
    methods
        
        function value = get.Enable( ~ )
            
            % Warn
            % warning( 'uiextras:Deprecated', ...
            %     'Property ''Enable'' will be removed in a future release.' )
            
            % Return
            value = 'on';
            
        end % get.Enable
        
        function set.Enable( ~, value )
            
            % Check
            assert( ischar( value ) && any( strcmp( value, {'on','off'} ) ), ...
                'uiextras:InvalidPropertyValue', ...
                'Property ''Enable'' must be ''on'' or ''off''.' )
            
            % Warn
            % warning( 'uiextras:Deprecated', ...
            %     'Property ''Enable'' will be removed in a future release.' )
            
        end % set.Enable
        
        function value = get.RowSizes( obj )
            
            % Get
            value = obj.Heights;
            
        end % get.RowSizes
        
        function set.RowSizes( obj, value )
            
            % Set
            obj.Heights = value;
            
        end % set.RowSizes
        
        function value = get.MinimumRowSizes( obj )
            
            % Get
            value = obj.MinimumHeights;
            
        end % get.MinimumRowSizes
        
        function set.MinimumRowSizes( obj, value )
            
            % Set
            obj.MinimumHeights = value;
            
        end % set.MinimumRowSizes
        
        function value = get.ColumnSizes( obj )
            
            % Get
            value = obj.Widths;
            
        end % get.ColumnSizes
        
        function set.ColumnSizes( obj, value )
            
            % Get
            obj.Widths = value;
            
        end % set.ColumnSizes
        
        function value = get.MinimumColumnSizes( obj )
            
            % Get
            value = obj.MinimumWidths;
            
        end % get.MinimumColumnSizes
        
        function set.MinimumColumnSizes( obj, value )
            
            % Get
            obj.MinimumWidths = value;
            
        end % set.MinimumColumnSizes
        
    end % accessors
    
end % classdef