classdef TabPanel < uix.Container & uix.mixin.Panel
    %uix.TabPanel  Tab panel
    %
    %  p = uix.TabPanel(p1,v1,p2,v2,...) constructs a tab panel and sets
    %  parameter p1 to value v1, etc.
    %
    %  A tab panel shows one of its contents and hides the others according
    %  to which tab is selected.
    %
    %  From R2014b, MATLAB provides uitabgroup and uitab as standard
    %  components.  Consider using uitabgroup and uitab for new code if
    %  these meet your requirements.
    %
    %  See also: uitabgroup, uitab, uix.CardPanel
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    properties( Access = public, Dependent, AbortSet )
        FontAngle % font angle
        FontName % font name
        FontSize % font size
        FontWeight % font weight
        FontUnits % font weight
        ForegroundColor % tab text color [RGB]
        HighlightColor % border highlight color [RGB]
        ShadowColor % border shadow color [RGB]
    end
    
    properties
        SelectionChangedFcn = '' % selection change callback
    end
    
    properties( Access = public, Dependent, AbortSet )
        TabEnables % tab enable states
        TabLocation % tab location [top|bottom]
        TabTitles % tab titles
        TabContextMenus % tab context menus
        TabWidth % tab width
    end
    
    properties( Access = private )
        FontAngle_ = get( 0, 'DefaultUicontrolFontAngle' ) % backing for FontAngle
        FontName_ = get( 0, 'DefaultUicontrolFontName' ) % backing for FontName
        FontSize_ = get( 0, 'DefaultUicontrolFontSize' ) % backing for FontSize
        FontWeight_ = get( 0, 'DefaultUicontrolFontWeight' ) % backing for FontWeight
        FontUnits_ = get( 0, 'DefaultUicontrolFontUnits' ) % backing for FontUnits
        ForegroundColor_ = get( 0, 'DefaultUicontrolForegroundColor' ) % backing for ForegroundColor
        HighlightColor_ = [1 1 1] % backing for HighlightColor
        ShadowColor_ = [0.7 0.7 0.7] % backing for ShadowColor
        ParentBackgroundColor = get( 0, 'DefaultUicontrolForegroundColor' ) % default parent background color
        Tabs = gobjects( [0 1] ) % tabs
        TabListeners = event.listener.empty( [0 1] ) % tab listeners
        TabLocation_ = 'top' % backing for TabPosition
        TabHeight = -1 % cache of tab height (-1 denotes stale cache)
        TabWidth_ = 50 % backing for TabWidth
        Dividers % tab dividers
        BackgroundColorListener % listener
        SelectionChangedListener % listener
        ParentListener % listener
        ParentBackgroundColorListener % listener
    end
    
    properties( Access = private, Constant )
        FontNames = listfonts() % all available font names
        DividerMask = uix.TabPanel.getDividerMask() % divider image data
        DividerWidth = 8 % divider width
        TabMinimumHeight = 9 % tab minimum height
        Tint = 0.85 % tint factor for unselected tabs
    end
    
    methods
        
        function obj = TabPanel( varargin )
            %uix.TabPanel  Tab panel constructor
            %
            %  p = uix.TabPanel() constructs a tab panel.
            %
            %  p = uix.TabPanel(p1,v1,p2,v2,...) sets parameter p1 to value
            %  v1, etc.
            
            % Create dividers
            dividers = matlab.ui.control.UIControl( 'Internal', true, ...
                'Parent', obj, 'Units', 'pixels', 'Style', 'checkbox',...
                'Tag', 'TabPanelDividers' );
            
            % Create listeners
            backgroundColorListener = event.proplistener( obj, ...
                findprop( obj, 'BackgroundColor' ), 'PostSet', ...
                @obj.onBackgroundColorChanged );
            selectionChangedListener = event.listener( obj, ...
                'SelectionChanged', @obj.onSelectionChanged );
            parentListener = event.proplistener( obj, ...
                findprop( obj, 'Parent' ), 'PostSet', ...
                @obj.onParentChanged );
            
            % Store properties
            obj.Dividers = dividers;
            obj.BackgroundColorListener = backgroundColorListener;
            obj.SelectionChangedListener = selectionChangedListener;
            obj.ParentListener = parentListener;
            
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
        
        function value = get.FontAngle( obj )
            
            value = obj.FontAngle_;
            
        end % get.FontAngle
        
        function set.FontAngle( obj, value )
            
            % Check
            assert( ischar( value ) && any( strcmp( value, {'normal','italic','oblique'} ) ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''FontAngle'' must be ''normal'', ''italic'' or ''oblique''.' )
            
            % Set
            obj.FontAngle_ = value;
            
            % Update existing tabs
            tabs = obj.Tabs;
            n = numel( tabs );
            for ii = 1:n
                tab = tabs(ii);
                tab.FontAngle = value;
            end
            
            % Mark as dirty
            obj.TabHeight = -1;
            obj.Dirty = true;
            
        end % set.FontAngle
        
        function value = get.FontName( obj )
            
            value = obj.FontName_;
            
        end % get.FontName
        
        function set.FontName( obj, value )
            
            % Check
            assert( ischar( value ) && any( strcmp( value, obj.FontNames ) ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''FontName'' must be a valid font name.' )
            
            % Set
            obj.FontName_ = value;
            
            % Update existing tabs
            tabs = obj.Tabs;
            n = numel( tabs );
            for ii = 1:n
                tab = tabs(ii);
                tab.FontName = value;
            end
            
            % Mark as dirty
            obj.TabHeight = -1;
            obj.Dirty = true;
            
        end % set.FontName
        
        function value = get.FontSize( obj )
            
            value = obj.FontSize_;
            
        end % get.FontSize
        
        function set.FontSize( obj, value )
            
            % Check
            assert( isa( value, 'double' ) && isscalar( value ) && ...
                isreal( value ) && ~isinf( value ) && ...
                ~isnan( value ) && value > 0, ...
                'uix:InvalidPropertyValue', ...
                'Property ''FontSize'' must be a positive scalar.' )
            
            % Set
            obj.FontSize_ = value;
            
            % Update existing tabs
            tabs = obj.Tabs;
            n = numel( tabs );
            for ii = 1:n
                tab = tabs(ii);
                tab.FontSize = value;
            end
            
            % Mark as dirty
            obj.TabHeight = -1;
            obj.Dirty = true;
            
        end % set.FontSize
        
        function value = get.FontWeight( obj )
            
            value = obj.FontWeight_;
            
        end % get.FontWeight
        
        function set.FontWeight( obj, value )
            
            % Check
            assert( ischar( value ) && any( strcmp( value, {'normal','bold'} ) ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''FontWeight'' must be ''normal'' or ''bold''.' )
            
            % Set
            obj.FontWeight_ = value;
            
            % Update existing tabs
            tabs = obj.Tabs;
            n = numel( tabs );
            for ii = 1:n
                tab = tabs(ii);
                tab.FontWeight = value;
            end
            
            % Mark as dirty
            obj.TabHeight = -1;
            obj.Dirty = true;
            
        end % set.FontWeight
        
        function value = get.FontUnits( obj )
            
            value = obj.FontUnits_;
            
        end % get.FontUnits
        
        function set.FontUnits( obj, value )
            
            % Check
            assert( ischar( value ) && ...
                any( strcmp( value, {'inches','centimeters','points','pixels'} ) ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''FontUnits'' must be ''inches'', ''centimeters'', ''points'' or ''pixels''.' )
            
            % Compute size in new units
            oldUnits = obj.FontUnits_;
            oldSize = obj.FontSize_;
            newUnits = value;
            newSize = oldSize * convert( oldUnits ) / convert( newUnits );
            
            % Set size and units
            obj.FontSize_ = newSize;
            obj.FontUnits_ = newUnits;
            
            % Update existing tabs
            tabs = obj.Tabs;
            n = numel( tabs );
            for ii = 1:n
                tab = tabs(ii);
                tab.FontUnits = newUnits;
            end
            
            % Mark as dirty
            obj.TabHeight = -1;
            obj.Dirty = true;
            
            function factor = convert( units )
                %convert  Compute conversion factor to points
                %
                %  f = convert(u) computes the conversion factor from units
                %  u to points.  For example, convert('inches') since 1
                %  inch equals 72 points.
                
                persistent SCREEN_PIXELS_PER_INCH
                if isequal( SCREEN_PIXELS_PER_INCH, [] ) % uninitialized
                    SCREEN_PIXELS_PER_INCH = get( 0, 'ScreenPixelsPerInch' );
                end
                
                switch units
                    case 'inches'
                        factor = 72;
                    case 'centimeters'
                        factor = 72 / 2.54;
                    case 'points'
                        factor = 1;
                    case 'pixels'
                        factor = 72 / SCREEN_PIXELS_PER_INCH;
                end
                
            end % convert
            
        end % set.FontUnits
        
        function value = get.ForegroundColor( obj )
            
            value = obj.ForegroundColor_;
            
        end % get.ForegroundColor
        
        function set.ForegroundColor( obj, value )
            
            % Check
            assert( isnumeric( value ) && isequal( size( value ), [1 3] ) && ...
                all( isreal( value ) ) && all( value >= 0 ) && all( value <= 1 ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''ForegroundColor'' must be an RGB triple.' )
            
            % Set
            obj.ForegroundColor_ = value;
            
            % Update existing tabs
            tabs = obj.Tabs;
            n = numel( tabs );
            for ii = 1:n
                tab = tabs(ii);
                tab.ForegroundColor = value;
            end
            
        end % set.ForegroundColor
        
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
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.HighlightColor
        
        function set.SelectionChangedFcn( obj, value )
            
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
                error( 'uix:InvalidPropertyValue', ...
                    'Property ''SelectionChangedFcn'' must be a valid callback.' )
            end
            
            % Set
            obj.SelectionChangedFcn = value;
            
        end % set.SelectionChangedFcn
        
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
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.ShadowColor
        
        function value = get.TabEnables( obj )
            
            value = get( obj.Tabs, {'Enable'} );
            value(strcmp( value, 'inactive' )) = {'on'};
            
        end % get.TabEnables
        
        function set.TabEnables( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Retrieve tabs
            tabs = obj.Tabs;
            tabListeners = obj.TabListeners;
            
            % Check
            assert( iscellstr( value ) && ...
                isequal( size( value ), size( tabs ) ) && ...
                all( strcmp( value, 'on' ) | strcmp( value, 'off' ) ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''TabEnables'' should be a cell array of strings ''on'' or ''off'', one per tab.' )
            
            % Set
            tf = strcmp( value, 'on' );
            value(tf) = {'inactive'};
            for ii = 1:numel( tabs )
                tabs(ii).Enable = value{ii};
                tabListeners(ii).Enabled = tf(ii);
            end
            
            % Show selected child
            obj.showSelection()
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.TabEnables
        
        function value = get.TabLocation( obj )
            
            value = obj.TabLocation_;
            
        end % get.TabLocation
        
        function set.TabLocation( obj, value )
            
            % Check
            assert( ischar( value ) && ...
                any( strcmp( value, {'top','bottom'} ) ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''TabLocation'' should be ''top'' or ''bottom''.' )
            
            % Set
            obj.TabLocation_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.TabLocation
        
        function value = get.TabTitles( obj )
            
            value = get( obj.Tabs, {'String'} );
            
        end % get.TabTitles
        
        function set.TabTitles( obj, value )
            
            % For those who can't tell a column from a row...
            if isrow( value )
                value = transpose( value );
            end
            
            % Retrieve tabs
            tabs = obj.Tabs;
            
            % Check
            assert( iscellstr( value ) && ...
                isequal( size( value ), size( tabs ) ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''TabTitles'' should be a cell array of strings, one per tab.' )
            
            % Set
            n = numel( tabs );
            for ii = 1:n
                tabs(ii).String = value{ii};
            end
            
            % Mark as dirty
            obj.TabHeight = -1;
            obj.Dirty = true;
            
        end % set.TabTitles
        
        function value = get.TabContextMenus( obj )
            
            tabs = obj.Tabs;
            n = numel( tabs );
            value = cell( [n 1] );
            for ii = 1:n
                value{ii} = tabs(ii).UIContextMenu;
            end
            
        end % get.TabContextMenus
        
        function set.TabContextMenus( obj, value )
            
            tabs = obj.Tabs;
            n = numel( tabs );
            for ii = 1:n
                tabs(ii).UIContextMenu = value{ii};
            end
            
        end % set.TabContextMenus
        
        function value = get.TabWidth( obj )
            
            value = obj.TabWidth_;
            
        end % get.TabWidth
        
        function set.TabWidth( obj, value )
            
            % Check
            assert( isa( value, 'double' ) && isscalar( value ) && ...
                isreal( value ) && ~isinf( value ) && ...
                ~isnan( value ) && value ~= 0, ...
                'uix:InvalidPropertyValue', ...
                'Property ''TabWidth'' must be a non-zero scalar.' )
            
            % Set
            obj.TabWidth_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.TabWidth
        
    end % accessors
    
    methods( Access = protected )
        
        function redraw( obj )
            
            % Compute positions
            bounds = hgconvertunits( ancestor( obj, 'figure' ), ...
                [0 0 1 1], 'normalized', 'pixels', obj );
            w = ceil( bounds(1) + bounds(3) ) - floor( bounds(1) ); % width
            h = ceil( bounds(2) + bounds(4) ) - floor( bounds(2) ); % height
            p = obj.Padding_; % padding
            tabs = obj.Tabs;
            n = numel( tabs ); % number of tabs
            tH = obj.TabHeight; % tab height
            if tH == -1 % cache stale, refresh
                if n > 0
                    cTabExtents = get( tabs, {'Extent'} );
                    tabExtents = vertcat( cTabExtents{:} );
                    tH = max( tabExtents(:,4) );
                end
                tH = max( tH, obj.TabMinimumHeight ); % apply minimum
                tH = ceil( tH ); % round up
                obj.TabHeight = tH; % store
            end
            cH = max( [h - 2 * p - tH, 1] ); % contents height
            switch obj.TabLocation_
                case 'top'
                    cY = 1 + p; % contents y
                    tY = cY + cH + p; % tab y
                case 'bottom'
                    tY = 1; % tab y
                    cY = tY + tH + p; % contents y
            end
            cX = 1 + p; % contents x
            cW = max( [w - 2 * p, 1] ); % contents width
            tW = obj.TabWidth_; % tab width
            dW = obj.DividerWidth; % tab divider width
            if tW < 0 && n > 0 % relative
                tW = max( ( w - (n+1) * dW ) / n, 1 );
            end
            tW = ceil( tW ); % round up
            for ii = 1:n
                tabs(ii).Position = [1 + (ii-1) * tW + ii * dW, tY, tW, tH];
            end
            obj.Dividers.Position = [0 tY w+1 tH];
            contentsPosition = [cX cY cW cH];
            
            % Redraw tabs
            obj.redrawTabs()
            
            % Redraw contents
            selection = obj.Selection_;
            if selection ~= 0 && strcmp( obj.TabEnables{selection}, 'on' )
                uix.setPosition( obj.Contents_(selection), contentsPosition, 'pixels' )
            end
            
        end % redraw
        
        function addChild( obj, child )
            %addChild  Add child
            %
            %  c.addChild(d) adds the child d to the container c.
            
            % Create new tab
            n = numel( obj.Tabs );
            tab = matlab.ui.control.UIControl( 'Internal', true, ...
                'Parent', obj, 'Style', 'text', 'Enable', 'inactive', ...
                'Units', 'pixels', 'FontUnits', obj.FontUnits_, ...
                'FontSize', obj.FontSize_, 'FontName', obj.FontName_, ...
                'FontAngle', obj.FontAngle_, 'FontWeight', obj.FontWeight_, ...
                'ForegroundColor', obj.ForegroundColor_, ...
                'String', sprintf( 'Page %d', n + 1 ) );
            tabListener = event.listener( tab, 'ButtonDown', @obj.onTabClicked );
            obj.Tabs(n+1,:) = tab;
            obj.TabListeners(n+1,:) = tabListener;
            
            % Mark as dirty
            obj.TabHeight = -1;
            
            % Check for bug
            if verLessThan( 'MATLAB', '8.5' ) && strcmp( child.Visible, 'off' )
                obj.G1218142 = true;
            end
            
            % Select new content
            oldSelection = obj.Selection_;
            if numel( obj.Contents_ ) == 0
                newSelection = 1;
                obj.Selection_ = newSelection;
            else
                newSelection = oldSelection;
            end
            
            % Call superclass method
            addChild@uix.mixin.Container( obj, child )
            
            % Show selected child
            obj.showSelection()
            
            % Notify selection change
            if oldSelection ~= newSelection
                obj.notify( 'SelectionChanged', ...
                    uix.SelectionData( oldSelection, newSelection ) )
            end
            
        end % addChild
        
        function removeChild( obj, child )
            %removeChild  Remove child
            %
            %  c.removeChild(d) removes the child d from the container c.
            
            % Find index of removed child
            contents = obj.Contents_;
            index = find( contents == child );
            
            % Remove tab
            delete( obj.Tabs(index) )
            obj.Tabs(index,:) = [];
            obj.TabListeners(index,:) = [];
            
            % Call superclass method
            removeChild@uix.mixin.Panel( obj, child )
            
        end % removeChild
        
        function reorder( obj, indices )
            %reorder  Reorder contents
            %
            %  c.reorder(i) reorders the container contents using indices
            %  i, c.Contents = c.Contents(i).
            
            % Reorder
            obj.Tabs = obj.Tabs(indices,:);
            obj.TabListeners = obj.TabListeners(indices,:);
            
            % Call superclass method
            reorder@uix.mixin.Panel( obj, indices )
            
        end % reorder
        
        function reparent( obj, oldFigure, newFigure )
            %reparent  Reparent container
            %
            %  c.reparent(a,b) reparents the container c from the figure a
            %  to the figure b.
            
            if ~isequal( oldFigure, newFigure )
                contextMenus = obj.TabContextMenus;
                for ii = 1:numel( contextMenus )
                    contextMenu = contextMenus{ii};
                    if ~isempty( contextMenu )
                        contextMenu.Parent = newFigure;
                    end
                end
            end
            
            % Call superclass method
            reparent@uix.mixin.Panel( obj, oldFigure, newFigure )
            
        end % reparent
        
        function showSelection( obj )
            %showSelection  Show selected child, hide the others
            %
            %  c.showSelection() shows the selected child of the container
            %  c, and hides the others.
            
            % Call superclass method
            showSelection@uix.mixin.Panel( obj )
            
            % If not enabled, hide selected contents too
            selection = obj.Selection_;
            if selection ~= 0 && strcmp( obj.TabEnables{selection}, 'off' )
                child = obj.Contents_(selection);
                child.Visible = 'off';
                if isa( child, 'matlab.graphics.axis.Axes' )
                    child.ContentsVisible = 'off';
                end
                % As a remedy for g1100294, move off-screen too
                margin = 1000;
                if isa( child, 'matlab.graphics.axis.Axes' ) ...
                        && strcmp(child.ActivePositionProperty, 'outerposition' )
                    child.OuterPosition(1) = -child.OuterPosition(3)-margin;
                else
                    child.Position(1) = -child.Position(3)-margin;
                end
            end
            
        end % showSelection
        
    end % template methods
    
    methods( Access = private )
        
        function redrawTabs( obj )
            %redrawTabs  Redraw tabs
            %
            %  p.redrawTabs() redraws the tabs.
            
            % Get relevant properties
            selection = obj.Selection_;
            tabs = obj.Tabs;
            t = numel( tabs );
            dividers = obj.Dividers;
            
            % Handle no tabs as a special case
            if t == 0
                dividers.Visible = 'off'; % hide
                return
            end
            
            % Repaint tabs
            backgroundColor = obj.BackgroundColor;
            for ii = 1:t
                tab = tabs(ii);
                if ii == selection
                    tab.BackgroundColor = backgroundColor;
                else
                    tab.BackgroundColor = obj.Tint * backgroundColor;
                end
            end
            
            % Repaint dividers
            d = t + 1;
            dividerNames = repmat( 'F', [d 2] ); % initialize
            dividerNames(1,1) = 'E'; % end
            dividerNames(end,2) = 'E'; % end
            if selection ~= 0
                dividerNames(selection,2) = 'T'; % selected
                dividerNames(selection+1,1) = 'T'; % selected
            end
            tH = obj.TabHeight;
            assert( tH >= obj.TabMinimumHeight, 'uix:InvalidState', ...
                'Cannot redraw tabs with invalid TabHeight.' )
            tW = obj.Tabs(1).Position(3);
            dW = obj.DividerWidth;
            allCData = zeros( [tH 0 3] ); % initialize
            map = [obj.ShadowColor; obj.BackgroundColor; ...
                obj.Tint * obj.BackgroundColor; obj.HighlightColor;...
                obj.ParentBackgroundColor];
            for ii = 1:d
                % Select mask
                iMask = obj.DividerMask.( dividerNames(ii,:) );
                % Resize
                iData = repmat( iMask(5,:), [tH 1] );
                iData(1:4,:) = iMask(1:4,:);
                iData(end-3:end,:) = iMask(end-3:end,:);
                % Convert to RGB
                cData = ind2rgb( iData+1, map );
                % Orient
                switch obj.TabLocation_
                    case 'bottom'
                        cData = flipud( cData );
                end
                % Insert
                allCData(1:tH,(ii-1)*(dW+tW)+(1:dW),:) = cData; % center
                if ii > 1 % extend left under transparent uicontrol edge
                    allCData(1:tH,(ii-1)*(dW+tW),:) = cData(:,1,:);
                end
                if ii < d % extend right under transparent uicontrol edge
                    allCData(1:tH,(ii-1)*(dW+tW)+dW+1,:) = cData(:,end,:);
                end
            end
            dividers.CData = allCData; % paint
            dividers.BackgroundColor = obj.ParentBackgroundColor;
            dividers.Visible = 'on'; % show
            
        end % redrawTabs
        
    end % helper methods
    
    methods( Access = private )
        
        function onTabClicked( obj, source, ~ )
            
            % Update selection
            oldSelection = obj.Selection_;
            newSelection = find( source == obj.Tabs );
            if oldSelection == newSelection, return, end % abort set
            obj.Selection_ = newSelection;
            
            % Show selected child
            obj.showSelection()
            
            % Mark as dirty
            obj.Dirty = true;
            
            % Notify selection change
            obj.notify( 'SelectionChanged', ...
                uix.SelectionData( oldSelection, newSelection ) )
            
        end % onTabClicked
        
        function onBackgroundColorChanged( obj, ~, ~ )
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % onBackgroundColorChanged
        
        function onSelectionChanged( obj, source, eventData )
            
            % Call callback
            callback = obj.SelectionChangedFcn;
            if ischar( callback ) && isequal( callback, '' )
                % do nothing
            elseif ischar( callback )
                feval( callback, source, eventData )
            elseif isa( callback, 'function_handle' )
                callback( source, eventData )
            elseif iscell( callback )
                feval( callback{1}, source, eventData, callback{2:end} )
            end
            
        end % onSelectionChanged
        
        function onParentChanged( obj, ~, ~ )
            
            % Update ParentBackgroundColor and ParentBackgroundColor
            if isprop( obj.Parent, 'BackgroundColor' )
                prop = 'BackgroundColor';
            elseif isprop( obj.Parent, 'Color' )
                prop = 'Color';
            else
                prop = [];
            end
            
            if ~isempty( prop )
                obj.ParentBackgroundColorListener = event.proplistener( obj.Parent, ...
                    findprop( obj.Parent, prop ), 'PostSet', ...
                    @( src, evt ) obj.updateParentBackgroundColor( prop ) );
            else
                obj.ParentBackgroundColorListener = [];
            end
            
            obj.updateParentBackgroundColor( prop );
            
        end % onParentChanged
        
        function updateParentBackgroundColor( obj, prop )
            
            if isempty( prop )
                obj.ParentBackgroundColor = obj.BackgroundColor;
            else
                obj.ParentBackgroundColor = obj.Parent.(prop);
            end
            
            % Mark as dirty
            obj.Dirty = true;
            
        end
        
    end % event handlers
    
    methods( Access = private, Static )
        
        function mask = getDividerMask()
            %getDividerMask  Get divider image data
            %
            %  m = uix.TabPanel.getDividerMask() returns the image masks
            %  for tab panel dividers.  Mask entries are 0 (shadow), 1
            %  (background), 2 (tint) and 3 (highlight).
            
            mask.EF = indexColor( uix.loadIcon( 'tab_NoEdge_NotSelected.png' ) );
            mask.ET = indexColor( uix.loadIcon( 'tab_NoEdge_Selected.png' ) );
            mask.FE = indexColor( uix.loadIcon( 'tab_NotSelected_NoEdge.png' ) );
            mask.FF = indexColor( uix.loadIcon( 'tab_NotSelected_NotSelected.png' ) );
            mask.FT = indexColor( uix.loadIcon( 'tab_NotSelected_Selected.png' ) );
            mask.TE = indexColor( uix.loadIcon( 'tab_Selected_NoEdge.png' ) );
            mask.TF = indexColor( uix.loadIcon( 'tab_Selected_NotSelected.png' ) );
            
            function mask = indexColor( rgbMap )
                %indexColor  Returns a map of index given an RGB map
                %
                %  mask = indexColor( rgbMap ) returns a mask of color
                %  index based on the supplied rgbMap.
                %  black  : 0
                %  red    : 1
                %  yellow : 2
                %  white  : 3
                %  blue   : 4
                mask = nan( size( rgbMap, 1 ),size( rgbMap, 2 ) );
                % Black
                colorIndex = isColor( rgbMap, [0 0 0] );
                mask(colorIndex) = 0;
                % Red
                colorIndex = isColor( rgbMap, [1 0 0] );
                mask(colorIndex) = 1;
                % Yellow
                colorIndex = isColor( rgbMap, [1 1 0] );
                mask(colorIndex) = 2;
                % White
                colorIndex = isColor( rgbMap, [1 1 1] );
                mask(colorIndex) = 3;
                % Blue
                colorIndex = isColor( rgbMap, [0 0 1] );
                mask(colorIndex) = 4;
                % Nested
                function boolMap = isColor( map, color )
                    %isColor  Return a map of boolean where map is equal to color
                    boolMap = all( bsxfun( @eq, map, permute( color, [1 3 2] ) ), 3 );
                end
            end
            
        end % getDividerMask
        
    end % static helper methods
    
end % classdef