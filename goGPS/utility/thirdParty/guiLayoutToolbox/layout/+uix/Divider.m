classdef ( Hidden ) Divider < matlab.mixin.SetGet
    %uix.Divider  Draggable divider
    %
    %  d = uix.Divider() creates a divider.
    %
    %  d = uix.Divider(p1,v1,p2,v2,...) creates a divider and sets
    %  specified property p1 to value v1, etc.
    
    %  Copyright 2009-2016 The MathWorks, Inc.
    %  $Revision: 1601 $ $Date: 2018-05-01 10:22:53 +0100 (Tue, 01 May 2018) $
    
    properties( Dependent )
        Parent % parent
        Units % units [inches|centimeters|characters|normalized|points|pixels]
        Position % position
        Visible % visible [on|off]
        BackgroundColor % background color [RGB]
        HighlightColor % border highlight color [RGB]
        ShadowColor % border shadow color [RGB]
        Orientation % orientation [vertical|horizontal]
        Markings % markings [pixels]
    end
    
    properties( Access = private )
        Control % uicontrol
        BackgroundColor_ = get( 0, 'DefaultUicontrolBackgroundColor' ) % backing for BackgroundColor
        HighlightColor_ = [1 1 1] % backing for HighlightColor
        ShadowColor_ = [0.7 0.7 0.7] % backing for ShadowColor
        Orientation_ = 'vertical' % backing for Orientation
        Markings_ = zeros( [0 1] ) % backing for Markings
        SizeChangedListener % listener
    end
    
    methods
        
        function obj = Divider( varargin )
            %uix.Divider  Draggable divider
            %
            %  d = uix.Divider() creates a divider.
            %
            %  d = uix.Divider(p1,v1,p2,v2,...) creates a dividerand sets
            %  specified property p1 to value v1, etc.
            
            % Create control
            control = matlab.ui.control.UIControl( ...
                'Style', 'checkbox', 'Internal', true, ...
                'Enable', 'inactive', 'DeleteFcn', @obj.onDeleted,...
                'Tag', 'uix.Divider' );
            
            % Store control
            obj.Control = control;
            
            % Set properties
            try
                uix.set( obj, varargin{:} )
            catch e
                delete( obj )
                e.throwAsCaller()
            end
            
            % Force update
            obj.update()
            
            % Create listener
            sizeChangedListener = event.listener( control, 'SizeChanged', ...
                @obj.onSizeChanged );
            
            % Store listener
            obj.SizeChangedListener = sizeChangedListener;
            
        end % constructor
        
        function delete( obj )
            %delete  Destructor
            
            control = obj.Control;
            if isgraphics( control ) && strcmp( control.BeingDeleted, 'off' )
                delete( control )
            end
            
        end % destructor
        
    end % structors
    
    methods
        
        function value = get.Parent( obj )
            
            value = obj.Control.Parent;
            
        end % get.Parent
        
        function set.Parent( obj, value )
            
            obj.Control.Parent = value;
            
        end % set.Parent
        
        function value = get.Units( obj )
            
            value = obj.Control.Units;
            
        end % get.Units
        
        function set.Units( obj, value )
            
            obj.Control.Units = value;
            
        end % set.Units
        
        function value = get.Position( obj )
            
            value = obj.Control.Position;
            
        end % get.Position
        
        function set.Position( obj, value )
            
            obj.Control.Position = value;
            
        end % set.Position
        
        function value = get.Visible( obj )
            
            value = obj.Control.Visible;
            
        end % get.Visible
        
        function set.Visible( obj, value )
            
            obj.Control.Visible = value;
            
        end % set.Visible
        
        function value = get.BackgroundColor( obj )
            
            value = obj.BackgroundColor_;
            
        end % get.BackgroundColor
        
        function set.BackgroundColor( obj, value )
            
            % Check
            assert( isa( value, 'double' ) && ...
                isequal( size( value ), [1 3] ) && ...
                all( value >= 0 ) && all( value <= 1 ), ...
                'uix:InvalidArgument', ...
                'Property ''BackgroundColor'' must be a valid colorspec.' )
            
            % Set
            obj.BackgroundColor_ = value;
            
            % Update
            obj.update()
            
        end % set.BackgroundColor
        
        function value = get.HighlightColor( obj )
            
            value = obj.HighlightColor_;
            
        end % get.HighlightColor
        
        function set.HighlightColor( obj, value )
            
            % Check
            assert( isnumeric( value ) && isequal( size( value ), [1 3] ) && ...
                all( isreal( value ) ) && all( value >= 0 ) && all( value <= 1 ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''HighlightColor'' must be an RGB triple.' )
            
            % Set
            obj.HighlightColor_ = value;
            
            % Update
            obj.update()
            
        end % set.HighlightColor
        
        function value = get.ShadowColor( obj )
            
            value = obj.ShadowColor_;
            
        end % get.ShadowColor
        
        function set.ShadowColor( obj, value )
            
            % Check
            assert( isnumeric( value ) && isequal( size( value ), [1 3] ) && ...
                all( isreal( value ) ) && all( value >= 0 ) && all( value <= 1 ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''ShadowColor'' must be an RGB triple.' )
            
            % Set
            obj.ShadowColor_ = value;
            
            % Update
            obj.update()
            
        end % set.ShadowColor
        
        function value = get.Orientation( obj )
            
            value = obj.Orientation_;
            
        end % get.Orientation
        
        function set.Orientation( obj, value )
            
            % Check
            assert( ischar( value ) && ismember( value, ...
                {'horizontal','vertical'} ) )
            
            % Set
            obj.Orientation_ = value;
            
            % Update
            obj.update()
            
        end % set.Orientation
        
        function value = get.Markings( obj )
            
            value = obj.Markings_;
            
        end % get.Markings
        
        function set.Markings( obj, value )
            
            % Check
            assert( isa( value, 'double' ) && ndims( value ) == 2 && ...
                size( value, 2 ) == 1 && all( isreal( value ) ) && ...
                all( ~isinf( value ) ) && all( ~isnan( value ) ) && ...
                all( value > 0 ), 'uix:InvalidPropertyValue', ...
                'Property ''Markings'' must be a vector of positive values.' ) %#ok<ISMAT>
            
            % Set
            obj.Markings_ = value;
            
            % Update
            obj.update()
            
        end % set.Markings
        
    end % accessors
    
    methods
        
        function tf = isMouseOver( obj, eventData )
            %isMouseOver  Test for mouse over
            %
            %  tf = d.isMouseOver(wmd) tests whether the WindowMouseData
            %  wmd is consistent with the mouse pointer being over the
            %  divider d.
            %
            %  This method returns false for dividers that are being
            %  deleted.
            
            tf = isvalid( obj ); % initialize
            for ii = 1:numel( obj )
                tf(ii) = tf(ii) && obj(ii).Control == eventData.HitObject;
            end
            
        end % isMouseOver
        
    end % methods
    
    methods( Access = private )
        
        function onDeleted( obj, ~, ~ )
            %onDeleted  Event handler
            
            % Call destructor
            obj.delete()
            
        end % onDeleted
        
        function onSizeChanged( obj, ~, ~ )
            %onSizeChanged  Event handler
            
            % Update
            obj.update()
            
        end % onSizeChanged
        
    end % event handlers
    
    methods( Access = private )
        
        function update( obj )
            %update  Update divider
            %
            %  d.update() updates the divider markings.
            
            % Get properties
            control = obj.Control;
            position = control.Position;
            backgroundColor = obj.BackgroundColor;
            highlightColor = obj.HighlightColor;
            shadowColor = obj.ShadowColor;
            orientation = obj.Orientation;
            markings = obj.Markings;
            
            % Assemble mask
            mask = zeros( floor( position([4 3]) ) - [1 1] ); % initialize
            switch orientation
                case 'vertical'
                    markings(markings < 4) = [];
                    markings(markings > position(4)-6) = [];
                    for ii = 1:numel( markings )
                        marking = markings(ii);
                        mask(floor( marking ) + [-3 0 3],1:end-1) = 1;
                        mask(floor( marking ) + [-2 1 4],1:end-1) = 2;
                    end
                case 'horizontal'
                    markings(markings < 4) = [];
                    markings(markings > position(3)-6) = [];
                    for ii = 1:numel( markings )
                        marking = markings(ii);
                        mask(2:end,floor( marking ) + [-3 0 3]) = 1;
                        mask(2:end,floor( marking ) + [-2 1 4]) = 2;
                    end
            end
            
            % Assemble color data
            cData1 = repmat( backgroundColor(1), size( mask ) );
            cData1(mask==1) = highlightColor(1);
            cData1(mask==2) = shadowColor(1);
            cData2 = repmat( backgroundColor(2), size( mask ) );
            cData2(mask==1) = highlightColor(2);
            cData2(mask==2) = shadowColor(2);
            cData3 = repmat( backgroundColor(3), size( mask ) );
            cData3(mask==1) = highlightColor(3);
            cData3(mask==2) = shadowColor(3);
            cData = cat( 3, cData1, cData2, cData3 );
            
            % Set properties
            control.ForegroundColor = backgroundColor;
            control.BackgroundColor = backgroundColor;
            control.CData = cData;
            
        end % update
        
    end % methods
    
end % classdef