classdef CardPanel < uix.CardPanel
    %uiextras.CardPanel  Show one element (card) from a list
    %
    %   obj = uiextras.CardPanel() creates a new card panel which allows
    %   selection between the different child objects contained, making the
    %   selected child fill the space available and all other children
    %   invisible. This is commonly used for creating wizards or quick
    %   switching between different views of a single data-set.
    %
    %   obj = uiextras.CardPanel(param,value,...) also sets one or more
    %   property values.
    %
    %   See the <a href="matlab:doc uiextras.CardPanel">documentation</a> for more detail and the list of properties.
    %
    %   Examples:
    %   >> f = figure();
    %   >> p = uiextras.CardPanel( 'Parent', f, 'Padding', 5 );
    %   >> uicontrol( 'Style', 'frame', 'Parent', p, 'Background', 'r' );
    %   >> uicontrol( 'Style', 'frame', 'Parent', p, 'Background', 'b' );
    %   >> uicontrol( 'Style', 'frame', 'Parent', p, 'Background', 'g' );
    %   >> p.SelectedChild = 2;
    %
    %   See also: uiextras.Panel
    %             uiextras.BoxPanel
    %             uiextras.TabPanel
    
    %  Copyright 2009-2014 The MathWorks, Inc.
    %  $Revision: 979 $ $Date: 2014-09-28 14:26:12 -0400 (Sun, 28 Sep 2014) $
    
    properties( Hidden, Access = public, Dependent )
        Enable % deprecated
        SelectedChild
    end
    
    methods
        
        function obj = CardPanel( varargin )
            
            % Call uix constructor
            obj@uix.CardPanel( varargin{:} )
            
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
        
        function value = get.SelectedChild( obj )
            
            % Get
            value = obj.Selection;
            
        end % get.SelectedChild
        
        function set.SelectedChild( obj, value )
            
            % Set
            obj.Selection = value;
            
        end % set.SelectedChild
        
    end % accessors
    
end % classdef