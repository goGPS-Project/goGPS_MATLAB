classdef HBoxFlex < uix.HBoxFlex
    %uiextras.HBoxFlex  A dynamically resizable horizontal layout
    %
    %   obj = uiextras.HBoxFlex() creates a new dynamically resizable
    %   horizontal box layout with all parameters set to defaults. The
    %   output is a new layout object that can be used as the parent for
    %   other user-interface components.
    %
    %   obj = uiextras.HBoxFlex(param,value,...) also sets one or more
    %   parameter values.
    %
    %   See the <a href="matlab:doc uiextras.HBoxFlex">documentation</a> for more detail and the list of properties.
    %
    %   Examples:
    %   >> f = figure( 'Name', 'uiextras.HBoxFlex example' );
    %   >> b = uiextras.HBoxFlex( 'Parent', f );
    %   >> uicontrol( 'Parent', b, 'Background', 'r' )
    %   >> uicontrol( 'Parent', b, 'Background', 'b' )
    %   >> uicontrol( 'Parent', b, 'Background', 'g' )
    %   >> uicontrol( 'Parent', b, 'Background', 'y' )
    %   >> set( b, 'Sizes', [-1 100 -2 -1], 'Spacing', 5 );
    %
    %   See also: uiextras.VBoxFlex
    %             uiextras.HBox
    %             uiextras.Grid
    
    %  Copyright 2009-2014 The MathWorks, Inc.
    %  $Revision: 979 $ $Date: 2014-09-28 14:26:12 -0400 (Sun, 28 Sep 2014) $
    
    properties( Hidden, Access = public, Dependent )
        Enable % deprecated
        Sizes
        MinimumSizes
        ShowMarkings
    end
    
    methods
        
        function obj = HBoxFlex( varargin )
            
            % Call uix constructor
            obj@uix.HBoxFlex( varargin{:} )
            
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
            value = transpose( obj.Widths );
            
        end % get.Sizes
        
        function set.Sizes( obj, value )
            
            % Set
            obj.Widths = value;
            
        end % set.Sizes
        
        function value = get.MinimumSizes( obj )
            
            % Get
            value = transpose( obj.MinimumWidths );
            
        end % get.MinimumSizes
        
        function set.MinimumSizes( obj, value )
            
            % Get
            obj.MinimumWidths = value;
            
        end % set.MinimumSizes
        
        function value = get.ShowMarkings( obj )
            
            % Get
            value = obj.DividerMarkings;
            
        end % get.ShowMarkings
        
        function set.ShowMarkings( obj, value )
            
            % Set
            obj.DividerMarkings = value;
            
        end % set.ShowMarkings
        
    end % accessors
    
end % classdef