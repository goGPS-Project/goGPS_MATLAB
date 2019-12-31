classdef TabPanel < uix.TabPanel
    %TabPanel  Show one element inside a tabbed panel
    %
    %   obj = uiextras.TabPanel() creates a panel with tabs along one edge
    %   to allow selection between the different child objects contained.
    %
    %   obj = uiextras.TabPanel(param,value,...) also sets one or more
    %   property values.
    %
    %   See the <a href="matlab:doc uiextras.TabPanel">documentation</a> for more detail and the list of properties.
    %
    %   Examples:
    %   >> f = figure();
    %   >> p = uiextras.TabPanel( 'Parent', f, 'Padding', 5 );
    %   >> uicontrol( 'Style', 'frame', 'Parent', p, 'Background', 'r' );
    %   >> uicontrol( 'Style', 'frame', 'Parent', p, 'Background', 'b' );
    %   >> uicontrol( 'Style', 'frame', 'Parent', p, 'Background', 'g' );
    %   >> p.TabNames = {'Red', 'Blue', 'Green'};
    %   >> p.SelectedChild = 2;
    %
    %   See also: uiextras.Panel
    %             uiextras.BoxPanel
    
    %  Copyright 2009-2014 The MathWorks, Inc.
    %  $Revision: 979 $ $Date: 2014-09-28 14:26:12 -0400 (Sun, 28 Sep 2014) $
    
    properties( Hidden, Access = public, Dependent )
        Callback
    end
    
    properties( Access = private )
        Callback_ = '' % backing for Callback
    end
    
    properties( Hidden, Access = public, Dependent )
        Enable % deprecated
        SelectedChild
        TabEnable
        TabNames
        TabPosition
        TabSize
    end
    
    properties( Access = private )
        SelectionChangedListener % listener
    end
    
    methods
        
        function obj = TabPanel( varargin )
            
            % Call uix constructor
            obj@uix.TabPanel( varargin{:} )
            
            % Auto-parent
            if ~ismember( 'Parent', varargin(1:2:end) )
                obj.Parent = gcf();
            end
            
            % Create listeners
            selectionChangedListener = event.listener( obj, ...
                'SelectionChanged', @obj.onSelectionChanged );
            
            % Store properties
            obj.SelectionChangedListener = selectionChangedListener;
            
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
        
        function value = get.Callback( obj )
            
            % Get
            value = obj.Callback_;
            
        end % get.Callback
        
        function set.Callback( obj, value )
            
            % Check
            if ischar( value ) % string
                % OK
            elseif isa( value, 'function_handle' ) && ...
                    isequal( size( value ), [1 1] ) % function handle
                % OK
            elseif iscell( value ) && ndims( value ) == 2 && ...
                    size( value, 1 ) == 1 && size( value, 2 ) > 0 && ...
                    isa( value{1}, 'function_handle' ) && ...
                    isequal( size( value{1} ), [1 1] ) %#ok<ISMAT> % cell callback
                % OK
            else
                error( 'uiextras:InvalidPropertyValue', ...
                    'Property ''Callback'' must be a valid callback.' )
            end
            
            % Set
            obj.Callback_ = value;
            
        end % set.Callback
        
        function value = get.SelectedChild( obj )
            
            % Get
            value = obj.Selection;
            
        end % get.SelectedChild
        
        function set.SelectedChild( obj, value )
            
            % Set
            obj.Selection = value;
            
        end % set.SelectedChild
        
        function value = get.TabEnable( obj )
            
            % Get
            value = transpose( obj.TabEnables );
            
        end % get.TabEnable
        
        function set.TabEnable( obj, value )
            
            % Set
            obj.TabEnables = value;
            
        end % set.TabEnable
        
        function value = get.TabNames( obj )
            
            % Get
            value = transpose( obj.TabTitles );
            
        end % get.TabNames
        
        function set.TabNames( obj, value )
            
            % Set
            obj.TabTitles = value;
            
        end % set.TabNames
        
        function value = get.TabPosition( obj )
            
            % Get
            value = obj.TabLocation;
            
        end % get.TabPosition
        
        function set.TabPosition( obj, value )
            
            % Set
            obj.TabLocation = value;
            
        end % set.TabPosition
        
        function value = get.TabSize( obj )
            
            % Get
            value = obj.TabWidth;
            
        end % get.TabSize
        
        function set.TabSize( obj, value )
            
            % Set
            obj.TabWidth = value;
            
        end % set.TabSize
        
    end % accessors
    
    methods( Access = private )
        
        function onSelectionChanged( obj, source, eventData )
            
            % Create legacy event data structure
            oldEventData = struct( 'Source', eventData.Source, ...
                'PreviousChild', eventData.OldValue, ...
                'SelectedChild', eventData.NewValue );
            
            % Call callback
            callback = obj.Callback_;
            if ischar( callback ) && isequal( callback, '' )
                % do nothing
            elseif ischar( callback )
                feval( callback, source, oldEventData )
            elseif isa( callback, 'function_handle' )
                callback( source, oldEventData )
            elseif iscell( callback )
                feval( callback{1}, source, oldEventData, callback{2:end} )
            end
            
        end % onSelectionChanged
        
    end % event handlers
    
end % classdef