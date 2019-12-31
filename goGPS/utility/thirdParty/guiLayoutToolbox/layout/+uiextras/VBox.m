classdef VBox < uix.VBox
    %uiextras.VBox  Arrange elements vertically in a single column
    %
    %   obj = uiextras.VBox() creates a new vertical box layout with all
    %   parameters set to defaults. The output is a new layout object that
    %   can be used as the parent for other user-interface components.
    %
    %   obj = uiextras.VBox(param,value,...) also sets one or more
    %   parameter values.
    %
    %   See the <a href="matlab:doc uiextras.VBox">documentation</a> for more detail and the list of properties.
    %
    %   Examples:
    %   >> f = figure();
    %   >> b = uiextras.VBox( 'Parent', f );
    %   >> uicontrol( 'Parent', b, 'Background', 'r' )
    %   >> uicontrol( 'Parent', b, 'Background', 'b' )
    %   >> uicontrol( 'Parent', b, 'Background', 'g' )
    %   >> set( b, 'Sizes', [-1 100 -2], 'Spacing', 5 );
    %
    %   >> f = figure();
    %   >> b1 = uiextras.VBox( 'Parent', f );
    %   >> b2 = uiextras.HBox( 'Parent', b1, 'Padding', 5, 'Spacing', 5 );
    %   >> uicontrol( 'Style', 'frame', 'Parent', b1, 'Background', 'r' )
    %   >> uicontrol( 'Parent', b2, 'String', 'Button1' )
    %   >> uicontrol( 'Parent', b2, 'String', 'Button2' )
    %   >> set( b1, 'Sizes', [30 -1] );
    %
    %   See also: uiextras.HBox
    %             uiextras.VBoxFlex
    %             uiextras.Grid
    
    %  Copyright 2009-2014 The MathWorks, Inc.
    %  $Revision: 1077 $ $Date: 2015-03-19 16:44:14 +0000 (Thu, 19 Mar 2015) $
    
    properties( Hidden, Access = public, Dependent )
        Enable % deprecated
        Sizes
        MinimumSizes
    end
    
    methods
        
        function obj = VBox( varargin )
            
            % Call uix constructor
            obj@uix.VBox( varargin{:} )
            
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
        
        function value = get.Sizes( obj )
            
            % Get
            value = transpose( obj.Heights );
            
        end % get.Sizes
        
        function set.Sizes( obj, value )
            
            % Set
            obj.Heights = value;
            
        end % set.Sizes
        
        function value = get.MinimumSizes( obj )
            
            % Get
            value = transpose( obj.MinimumHeights );
            
        end % get.MinimumSizes
        
        function set.MinimumSizes( obj, value )
            
            % Get
            obj.MinimumHeights = value;
            
        end % set.MinimumSizes
        
    end % accessors
    
end % classdef