classdef Panel < uix.Panel
    %uiextras.Panel  Show one element inside a panel
    %
    %   obj = uiextras.Panel() creates a standard UIPANEL object but with
    %   automatic management of the contained widget or layout. The
    %   properties available are largely the same as the builtin UIPANEL
    %   object. Where more than one child is added, the currently visible
    %   child is determined using the SelectedChild property.
    %
    %   obj = uiextras.Panel(param,value,...) also sets one or more
    %   property values.
    %
    %   See the <a href="matlab:doc uiextras.Panel">documentation</a> for more detail and the list of properties.
    %
    %   Examples:
    %   >> f = figure();
    %   >> p = uiextras.Panel( 'Parent', f, 'Title', 'A Panel', 'Padding', 5 );
    %   >> uicontrol( 'Parent', p, 'Background', 'r' )
    %
    %   >> f = figure();
    %   >> p = uiextras.Panel( 'Parent', f, 'Title', 'A Panel', 'Padding', 5 );
    %   >> b = uiextras.HBox( 'Parent', p, 'Spacing', 5 );
    %   >> uicontrol( 'Style', 'listbox', 'Parent', b, 'String', {'Item 1','Item 2'} );
    %   >> uicontrol( 'Parent', b, 'Background', 'b' );
    %   >> set( b, 'Sizes', [100 -1] );
    %
    %   See also: uipanel
    %             uiextras.BoxPanel
    %             uiextras.HBox
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    properties( Hidden, Access = public, Dependent )
        SelectedChild
    end
    
    methods
        
        function obj = Panel( varargin )
            
            % Call uix constructor
            obj@uix.Panel( varargin{:} )
            
            % Add Enable property
            if ~isprop( obj, 'Enable' )
                p = addprop( obj, 'Enable' );
                p.GetMethod = @getEnable;
                p.SetMethod = @setEnable;
            end
            
            % Auto-parent
            if ~ismember( 'Parent', varargin(1:2:end) )
                obj.Parent = gcf();
            end
            
        end % constructor
        
    end % structors
    
    methods
        
        function value = getEnable( ~ )
            
            % Warn
            % warning( 'uiextras:Deprecated', ...
            %     'Property ''Enable'' will be removed in a future release.' )
            
            % Return
            value = 'on';
            
        end % getEnable
        
        function setEnable( ~, value )
            
            % Check
            assert( ischar( value ) && any( strcmp( value, {'on','off'} ) ), ...
                'uiextras:InvalidPropertyValue', ...
                'Property ''Enable'' must be ''on'' or ''off''.' )
            
            % Warn
            % warning( 'uiextras:Deprecated', ...
            %     'Property ''Enable'' will be removed in a future release.' )
            
        end % setEnable
        
        function value = get.SelectedChild( obj )
            
            % Warn
            % warning( 'uiextras:Deprecated', ...
            %     'Property ''SelectedChild'' will be removed in a future release.' )
            
            % Get
            if isempty( obj.Contents_ )
                value = [];
            else
                value = 1;
            end
            
        end % get.SelectedChild
        
        function set.SelectedChild( ~, ~ )
            
            % Warn
            % warning( 'uiextras:Deprecated', ...
            %     'Property ''SelectedChild'' will be removed in a future release.' )
            
        end % set.SelectedChild
        
    end % accessors
    
end % classdef