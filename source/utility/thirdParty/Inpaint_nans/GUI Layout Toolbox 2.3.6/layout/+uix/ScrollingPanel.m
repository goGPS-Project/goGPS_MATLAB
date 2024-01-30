classdef ScrollingPanel < uix.Container & uix.mixin.Panel
    %uix.ScrollingPanel  Scrolling panel
    %
    %  p = uix.ScrollingPanel(p1,v1,p2,v2,...) constructs a scrolling panel
    %  and sets parameter p1 to value v1, etc.
    %
    %  A scrolling panel is a standard container (uicontainer) that shows
    %  one its contents and hides the others.
    %
    %  See also: uix.Panel, uix.BoxPanel, uix.TabPanel, uicontainer
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    properties( Dependent )
        Heights % heights of contents, in pixels and/or weights
        MinimumHeights % minimum heights of contents, in pixels
        VerticalOffsets % vertical offsets of contents, in pixels
        VerticalSteps % vertical slider steps, in pixels
        Widths % widths of contents, in pixels and/or weights
        MinimumWidths % minimum widths of contents, in pixels
        HorizontalOffsets % horizontal offsets of contents, in pixels
        HorizontalSteps % horizontal slider steps, in pixels
        MouseWheelEnabled % mouse wheel scrolling enabled [on|off]
    end
    
    properties( Access = protected )
        Heights_ = zeros( [0 1] ) % backing for Heights
        MinimumHeights_ = zeros( [0 1] ) % backing for MinimumHeights
        Widths_ = zeros( [0 1] ) % backing for Widths
        MinimumWidths_ = zeros( [0 1] ) % backing for MinimumWidths
        HorizontalSliders = matlab.ui.control.UIControl.empty( [0 1] ) % sliders
        VerticalSliders = matlab.ui.control.UIControl.empty( [0 1] ) % sliders
        BlankingPlates = matlab.ui.control.UIControl.empty( [0 1] ) % blanking plates
        HorizontalSteps_ = zeros( [0 1] ) % steps
        VerticalSteps_ = zeros( [0 1] ) % steps
    end
    
    properties( Access = private )
        MouseWheelListener = [] % mouse listener
        MouseWheelEnabled_ = 'on' % backing for MouseWheelEnabled
        ScrollingListener = [] % slider listener
        ScrolledListener = [] % slider listener
        BackgroundColorListener % property listener
    end
    
    properties( Constant, Access = protected )
        SliderSize = 20 % slider size, in pixels
        SliderStep = 10 % slider step, in pixels
    end
    
    events( NotifyAccess = private )
        Scrolling
        Scrolled
    end
    
    methods
        
        function obj = ScrollingPanel( varargin )
            %uix.ScrollingPanel  Scrolling panel constructor
            %
            %  p = uix.ScrollingPanel() constructs a scrolling panel.
            %
            %  p = uix.ScrollingPanel(p1,v1,p2,v2,...) sets parameter p1 to
            %  value v1, etc.
            
            % Create listeners
            backgroundColorListener = event.proplistener( obj, ...
                findprop( obj, 'BackgroundColor' ), 'PostSet', ...
                @obj.onBackgroundColorChanged );
            
            % Store properties
            obj.BackgroundColorListener = backgroundColorListener;
            
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
        
        function value = get.VerticalOffsets( obj )
            
            sliders = obj.VerticalSliders;
            if isempty( sliders )
                value = zeros( size( sliders ) );
            else
                value = -vertcat( sliders.Value ) - 1;
                value(value<0) = 0;
            end
            
        end % get.VerticalOffsets
        
        function set.VerticalOffsets( obj, value )
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''VerticalOffsets'' must be of type double.' )
            assert( all( isreal( value ) ) && ~any( isinf( value ) ) && ...
                ~any( isnan( value ) ), 'uix:InvalidPropertyValue', ...
                'Elements of property ''VerticalOffsets'' must be real and finite.' )
            assert( isequal( size( value ), size( obj.Contents_ ) ), ...
                'uix:InvalidPropertyValue', ...
                'Size of property ''VerticalOffsets'' must match size of contents.' )
            
            % Set
            sliders = obj.VerticalSliders;
            for ii = 1:numel( sliders )
                sliders(ii).Value = -value(ii) - 1;
            end
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.VerticalOffsets
        
        function value = get.VerticalSteps( obj )
            
            value = obj.VerticalSteps_;
            
        end % get.VerticalSteps
        
        function set.VerticalSteps( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''VerticalSteps'' must be of type double.' )
            assert( all( isreal( value ) ) && ~any( isinf( value ) ) && ...
                ~any( isnan( value ) ) && all( value > 0 ), ...
                'uix:InvalidPropertyValue', ...
                'Elements of property ''VerticalSteps'' must be real, finite and positive.' )
            assert( isequal( size( value ), size( obj.Contents_ ) ), ...
                'uix:InvalidPropertyValue', ...
                'Size of property ''VerticalSteps'' must match size of contents.' )
            
            % Set
            obj.VerticalSteps_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.VerticalSteps
        
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
        
        function value = get.HorizontalOffsets( obj )
            
            sliders = obj.HorizontalSliders;
            if isempty( sliders )
                value = zeros( size( sliders ) );
            else
                value = vertcat( sliders.Value );
                value(value<0) = 0;
            end
            
        end % get.HorizontalOffsets
        
        function set.HorizontalOffsets( obj, value )
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''HorizontalOffsets'' must be of type double.' )
            assert( all( isreal( value ) ) && ~any( isinf( value ) ) && ...
                ~any( isnan( value ) ), 'uix:InvalidPropertyValue', ...
                'Elements of property ''HorizontalOffsets'' must be real and finite.' )
            assert( isequal( size( value ), size( obj.Contents_ ) ), ...
                'uix:InvalidPropertyValue', ...
                'Size of property ''HorizontalOffsets'' must match size of contents.' )
            
            % Set
            sliders = obj.HorizontalSliders;
            for ii = 1:numel( sliders )
                sliders(ii).Value = value(ii);
            end
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.HorizontalOffsets
        
        function value = get.HorizontalSteps( obj )
            
            value = obj.HorizontalSteps_;
            
        end % get.HorizontalSteps
        
        function set.HorizontalSteps( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''HorizontalSteps'' must be of type double.' )
            assert( all( isreal( value ) ) && ~any( isinf( value ) ) && ...
                ~any( isnan( value ) ) && all( value > 0 ), ...
                'uix:InvalidPropertyValue', ...
                'Elements of property ''HorizontalSteps'' must be real, finite and positive.' )
            assert( isequal( size( value ), size( obj.Contents_ ) ), ...
                'uix:InvalidPropertyValue', ...
                'Size of property ''HorizontalSteps'' must match size of contents.' )
            
            % Set
            obj.HorizontalSteps_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.HorizontalSteps
        
        function value = get.MouseWheelEnabled( obj )
            
            value = obj.MouseWheelEnabled_;
            
        end % get.MouseWheelEnabled
        
        function set.MouseWheelEnabled( obj, value )
            
            value = uix.validateScalarStringOrCharacterArray( value, ...
                'MouseWheelEnabled' );
            assert( any( strcmp( value, {'on','off'} ) ), ...
                'uix:InvalidArgument', ...
                'Property ''MouseWheelEnabled'' must ''on'' or ''off''.' )
            listener = obj.MouseWheelListener;
            if ~isempty( listener )
                listener.Enabled = strcmp( value, 'on' );
            end
            obj.MouseWheelEnabled_ = value;
            
        end % set.MouseWheelEnabled
        
    end % accessors
    
    methods( Access = protected )
        
        function redraw( obj )
            %redraw  Redraw
            
            % Return if no contents
            selection = obj.Selection_;
            if selection == 0, return, end
            
            % Retrieve width and height of selected contents
            contentsWidth = obj.Widths_(selection);
            minimumWidth = obj.MinimumWidths_(selection);
            contentsHeight = obj.Heights_(selection);
            minimumHeight = obj.MinimumHeights_(selection);
            
            % Retrieve selected contents and corresponding decorations
            child = obj.Contents_(selection);
            vSlider = obj.VerticalSliders(selection);
            hSlider = obj.HorizontalSliders(selection);
            plate = obj.BlankingPlates(selection);
            
            % Compute dimensions
            bounds = hgconvertunits( ancestor( obj, 'figure' ), ...
                [0 0 1 1], 'normalized', 'pixels', obj );
            width = bounds(3);
            height = bounds(4);
            sliderSize = obj.SliderSize; % slider size
            vSliderWidth = sliderSize * ...
                (contentsHeight > height | ...
                minimumHeight > height); % first pass
            hSliderHeight = sliderSize * ...
                (contentsWidth > width - vSliderWidth | ...
                minimumWidth > width - vSliderWidth);
            vSliderWidth = sliderSize * ...
                (contentsHeight > height - hSliderHeight | ...
                minimumHeight > height - hSliderHeight); % second pass
            vSliderWidth = min( vSliderWidth, width ); % limit
            hSliderHeight = min( hSliderHeight, height ); % limit
            vSliderHeight = height - hSliderHeight;
            hSliderWidth = width - vSliderWidth;
            widths = uix.calcPixelSizes( width, ...
                [contentsWidth;vSliderWidth], ...
                [minimumWidth;vSliderWidth], 0, 0 );
            contentsWidth = widths(1); % to be offset
            heights = uix.calcPixelSizes( height, ...
                [contentsHeight;hSliderHeight], ...
                [minimumHeight;hSliderHeight], 0, 0 );
            contentsHeight = heights(1); % to be offset
            
            % Compute positions
            contentsPosition = [1 1+hSliderHeight+vSliderHeight-contentsHeight contentsWidth contentsHeight];
            vSliderPosition = [1+hSliderWidth 1+hSliderHeight vSliderWidth vSliderHeight];
            hSliderPosition = [1 1 hSliderWidth hSliderHeight];
            platePosition = [1+hSliderWidth 1 vSliderWidth hSliderHeight];
            
            % Compute and set vertical slider properties
            if vSliderWidth == 0 || vSliderHeight == 0 || vSliderHeight <= vSliderWidth
                % Slider is invisible or incorrectly oriented
                set( vSlider, 'Style', 'text', 'Enable', 'inactive', ...
                    'Position', vSliderPosition, ...
                    'Min', 0, 'Max', 1, 'Value', 1 )
            else
                % Compute properties
                vSliderMin = 0;
                vSliderMax = contentsHeight - vSliderHeight;
                vSliderValue = -vSlider.Value; % negative sign convention
                vSliderValue = max( vSliderValue, vSliderMin ); % limit
                vSliderValue = min( vSliderValue, vSliderMax ); % limit
                vStep = obj.VerticalSteps_(selection);
                vSliderStep(1) = min( vStep / vSliderMax, 1 );
                vSliderStep(2) = max( vSliderHeight / vSliderMax, vSliderStep(1) );
                contentsPosition(2) = contentsPosition(2) + vSliderValue;
                % Set properties
                set( vSlider, 'Style', 'slider', 'Enable', 'on', ...
                    'Position', vSliderPosition, ...
                    'Min', -vSliderMax, 'Max', -vSliderMin, ...
                    'Value', -vSliderValue, 'SliderStep', vSliderStep )
            end
            
            % Compute and set horizontal slider properties
            if hSliderHeight == 0 || hSliderWidth == 0 || hSliderWidth <= hSliderHeight
                % Slider is invisible or incorrectly oriented
                set( hSlider, 'Style', 'text', 'Enable', 'inactive', ...
                    'Position', hSliderPosition, ...
                    'Min', -1, 'Max', 0, 'Value', -1 )
            else
                % Compute properties
                hSliderMin = 0;
                hSliderMax = contentsWidth - hSliderWidth;
                hSliderValue = hSlider.Value; % positive sign convention
                hSliderValue = max( hSliderValue, hSliderMin ); % limit
                hSliderValue = min( hSliderValue, hSliderMax ); % limit
                hStep = obj.HorizontalSteps_(selection);
                hSliderStep(1) = min( hStep / hSliderMax, 1 );
                hSliderStep(2) = max( hSliderWidth / hSliderMax, hSliderStep(1) );
                contentsPosition(1) = contentsPosition(1) - hSliderValue;
                % Set properties
                set( hSlider, 'Style', 'slider', 'Enable', 'on', ...
                    'Position', hSliderPosition, ...
                    'Min', hSliderMin, 'Max', hSliderMax, ...
                    'Value', hSliderValue, 'SliderStep', hSliderStep )
            end
            
            % Set contents and blanking plate positions
            uix.setPosition( child, contentsPosition, 'pixels' )
            set( plate, 'Position', platePosition )
            
        end % redraw
        
        function addChild( obj, child )
            %addChild  Add child
            %
            %  c.addChild(d) adds the child d to the container c.
            
            % Create decorations
            verticalSlider = matlab.ui.control.UIControl( ...
                'Internal', true, 'Parent', obj, ...
                'Units', 'pixels', 'Style', 'slider', ...
                'BackgroundColor', obj.BackgroundColor );
            horizontalSlider = matlab.ui.control.UIControl( ...
                'Internal', true, 'Parent', obj, ...
                'Units', 'pixels', 'Style', 'slider', ...
                'BackgroundColor', obj.BackgroundColor );
            blankingPlate = matlab.ui.control.UIControl( ...
                'Internal', true, 'Parent', obj, ...
                'Units', 'pixels', 'Style', 'text', 'Enable', 'inactive', ...
                'BackgroundColor', obj.BackgroundColor );
            
            % Add to sizes
            obj.Widths_(end+1,:) = -1;
            obj.MinimumWidths_(end+1,:) = 1;
            obj.Heights_(end+1,:) = -1;
            obj.MinimumHeights_(end+1,:) = 1;
            obj.VerticalSliders(end+1,:) = verticalSlider;
            obj.HorizontalSliders(end+1,:) = horizontalSlider;
            obj.BlankingPlates(end+1,:) = blankingPlate;
            obj.VerticalSteps_(end+1,:) = obj.SliderStep;
            obj.HorizontalSteps_(end+1,:) = obj.SliderStep;
            obj.updateSliderListeners()
            
            % Call superclass method
            addChild@uix.mixin.Panel( obj, child )
            
        end % addChild
        
        function removeChild( obj, child )
            %removeChild  Remove child
            %
            %  c.removeChild(d) removes the child d from the container c.
            
            % Identify child
            tf = obj.Contents_ == child;
            
            % Destroy decorations
            delete( obj.VerticalSliders(tf,:) )
            delete( obj.HorizontalSliders(tf,:) )
            delete( obj.BlankingPlates(tf,:) )
            
            % Remove from sizes
            obj.Widths_(tf,:) = [];
            obj.MinimumWidths_(tf,:) = [];
            obj.Heights_(tf,:) = [];
            obj.MinimumHeights_(tf,:) = [];
            obj.VerticalSliders(tf,:) = [];
            obj.HorizontalSliders(tf,:) = [];
            obj.BlankingPlates(tf,:) = [];
            obj.VerticalSteps_(tf,:) = [];
            obj.HorizontalSteps_(tf,:) = [];
            obj.updateSliderListeners()
            
            % Call superclass method
            removeChild@uix.mixin.Panel( obj, child )
            
        end % removeChild
        
        function reparent( obj, ~, newFigure )
            %reparent  Reparent container
            %
            %  c.reparent(a,b) reparents the container c from the figure a
            %  to the figure b.
            
            if isempty( newFigure )
                obj.MouseWheelListener = [];
            else
                listener = event.listener( newFigure, ...
                    'WindowScrollWheel', @obj.onMouseScrolled );
                listener.Enabled = strcmp( obj.MouseWheelEnabled_, 'on' );
                obj.MouseWheelListener = listener;
            end
            
        end % reparent
        
        function reorder( obj, indices )
            %reorder  Reorder contents
            %
            %  c.reorder(i) reorders the container contents using indices
            %  i, c.Contents = c.Contents(i).
            
            % Reorder
            obj.Widths_ = obj.Widths_(indices,:);
            obj.MinimumWidths_ = obj.MinimumWidths_(indices,:);
            obj.Heights_ = obj.Heights_(indices,:);
            obj.MinimumHeights_ = obj.MinimumWidths_(indices,:);
            obj.VerticalSliders = obj.VerticalSliders(indices,:);
            obj.HorizontalSliders = obj.HorizontalSliders(indices,:);
            obj.BlankingPlates = obj.BlankingPlates(indices,:);
            obj.VerticalSteps_ = obj.VerticalSteps_(indices,:);
            obj.HorizontalSteps_ = obj.HorizontalSteps_(indices,:);
            
            % Call superclass method
            reorder@uix.mixin.Panel( obj, indices )
            
        end % reorder
        
        function showSelection( obj )
            %showSelection  Show selected child, hide the others
            %
            %  c.showSelection() shows the selected child of the container
            %  c, and hides the others.
            
            % Call superclass method
            showSelection@uix.mixin.Panel( obj )
            
            % Show and hide sliders based on selection
            selection = obj.Selection_;
            for ii = 1:numel( obj.Contents_ )
                if ii == selection
                    obj.VerticalSliders(ii).Visible = 'on';
                    obj.HorizontalSliders(ii).Visible = 'on';
                    obj.BlankingPlates(ii).Visible = 'on';
                else
                    obj.VerticalSliders(ii).Visible = 'off';
                    obj.HorizontalSliders(ii).Visible = 'off';
                    obj.BlankingPlates(ii).Visible = 'off';
                end
            end
            
        end % showSelection
        
    end % template methods
    
    methods( Access = ?matlab.unittest.TestCase )
        
        function onSliderScrolling( obj, ~, ~ )
            %onSliderScrolling  Event handler
            
            % Mark as dirty
            obj.Dirty = true;
            
            % Raise event
            notify( obj, 'Scrolling' )
            
        end % onSliderScrolling
        
        function onSliderScrolled( obj, ~, ~ )
            %onSliderScrolled  Event handler
            
            % Mark as dirty
            obj.Dirty = true;
            
            % Raise event
            notify( obj, 'Scrolled' )
            
        end % onSliderScrolled
        
        function onMouseScrolled( obj, ~, eventData )
            %onMouseScrolled  Event handler
            
            sel = obj.Selection_;
            if sel == 0
                return
            else
                % Get pointer position and panel bounds
                pp = getpixelposition( obj, true );
                f = ancestor( obj, 'figure' );
                cpu = f.CurrentPoint; % figure Units
                cpwhu = [cpu 0 0]; % [x y] to [x y w h] for hgconvertunits
                cpwh = hgconvertunits( f, cpwhu, f.Units, 'pixels', obj ); % pixels
                cp = cpwh(1:2); % [x y w h] to [x y]
                
                % Check that pointer is over panel
                if cp(1) < pp(1) || cp(1) > pp(1) + pp(3) || ...
                        cp(2) < pp(2) || cp(2) > pp(2) + pp(4), return, end
                % Scroll
                if strcmp( obj.VerticalSliders(sel).Enable, 'on' ) % scroll vertically
                    delta = eventData.VerticalScrollCount * ...
                        eventData.VerticalScrollAmount * obj.VerticalSteps(sel);
                    obj.VerticalOffsets(sel) = obj.VerticalOffsets(sel) + delta;
                elseif strcmp( obj.HorizontalSliders(sel).Enable, 'on' ) % scroll horizontally
                    delta = eventData.VerticalScrollCount * ...
                        eventData.VerticalScrollAmount * obj.HorizontalSteps(sel);
                    obj.HorizontalOffsets(sel) = obj.HorizontalOffsets(sel) + delta;
                end
                % Raise event
                notify( obj, 'Scrolled' )
            end
            
        end % onMouseScrolled
        
        function onBackgroundColorChanged( obj, ~, ~ )
            %onBackgroundColorChanged  Handler for BackgroundColor changes
            
            set( obj.HorizontalSliders, 'BackgroundColor', obj.BackgroundColor )
            set( obj.VerticalSliders, 'BackgroundColor', obj.BackgroundColor )
            set( obj.BlankingPlates, 'BackgroundColor', obj.BackgroundColor )
            
        end % onBackgroundColorChanged
        
    end % event handlers
    
    methods( Access = private )
        
        function updateSliderListeners( obj )
            %updateSliderListeners  Update listeners to slider events
            
            if isempty( obj.VerticalSliders )
                obj.ScrollingListener = [];
                obj.ScrolledListener = [];
            else
                obj.ScrollingListener = event.listener( ...
                    [obj.VerticalSliders; obj.HorizontalSliders], ...
                    'ContinuousValueChange', @obj.onSliderScrolling );
                obj.ScrolledListener = event.listener( ...
                    [obj.VerticalSliders; obj.HorizontalSliders], ...
                    'Action', @obj.onSliderScrolled );
            end
            
        end % updateSliderListeners
        
    end % helpers
    
end % classdef