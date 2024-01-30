classdef (Sealed) Toolbar < matlab.ui.componentcontainer.ComponentContainer & ...
        wt.mixin.BackgroundColorable & ...
        wt.mixin.TitleColorable & ...
        wt.mixin.FontStyled & ...
        wt.mixin.PropertyViewable

    % A configurable toolbar

    % Copyright 2020-2023 The MathWorks Inc.

    %% Events
    events (HasCallbackProperty, NotifyAccess = protected)

        % Triggered when a button is pushed
        ButtonPushed

    end %events



    %% Public properties
    properties (AbortSet)

        % Sections that are part of the toolbar
        Section (:,1) wt.toolbar.HorizontalSection

    end %properties


    properties (Dependent, AbortSet, UsedInUpdate = false)

        % Color of the dividers between horizontal sections
        DividerColor

    end %properties



    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)
        
        % The listbox control
        ListBox (1,1) matlab.ui.control.ListBox

        % Grid
        Grid (1,1) matlab.ui.container.GridLayout

        % The label for each section
        SectionLabel (:,1) matlab.ui.control.Label

        % Buttons used when space is limited
        SectionButton (:,1) matlab.ui.control.Button

        % The dummy section filling up remaining space
        DummySection (:,1) matlab.ui.container.Container

        % Listen to section changes
        SectionChangedListener event.listener

        % Listen to size changes
        SizeChangedListener event.listener

        % Listen to button pushes in sections
        ButtonPushedListener event.listener

    end %properties


    properties (Dependent, NonCopyable, Hidden, SetAccess = private)
        
        % Indicates which sections are open
        SectionIsOpen (:,1) logical

    end %properties


    properties (Constant, Access = private)

        % The down arrow mask for the button icons
        BUTTON_MASK (:,:) logical = sectionButtonIconMask()

    end %properties



    %% Protected methods
    methods (Access = protected)

        function setup(obj)

            % Set default size
            obj.Position(3:4) = [500 90];

            % Construct Grid Layout to Manage Building Blocks
            obj.Grid = uigridlayout(obj);
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.RowSpacing = 2;
            obj.Grid.ColumnSpacing = 2;
            obj.Grid.Padding = 0;

            % Configure style defaults
            obj.FontColor = [1 1 1] * 0.3;
            obj.TitleColor = [1 1 1] * 0.5;

            % Configure grid
            obj.Grid.Padding = [0 0 0 0]; %If changed, modify updateLayout!
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.ColumnWidth = {};
            obj.Grid.ColumnSpacing = 1; %If changed, modify updateLayout!

            obj.Grid.RowHeight = {'1x',15};
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowSpacing = 0;
            obj.Grid.BackgroundColor = [1 1 1]*0.5;

            % Add a dummy section to color the empty space
            obj.DummySection = uicontainer(obj.Grid);
            obj.DummySection.Layout.Row = [1 2];

            % Set Colorable Background Objects
            obj.BackgroundColorableComponents = [obj.DummySection obj.Grid];

            % Listen to size changes
            obj.SizeChangedListener = event.listener(obj,'SizeChanged',...
                @(h,e)obj.updateLayout());

        end %function


        function update(obj)
            % Update the toolbar

            % Update each section
            numSections = numel(obj.Section);
            for sIdx = 1:numSections

                % Do we need to add more labels/buttons?
                if sIdx > numel(obj.SectionLabel)
                    obj.SectionLabel(sIdx) = uilabel(obj.Grid);
                    obj.SectionButton(sIdx) = uibutton(obj.Grid);
                end

                % Update each label
                obj.SectionLabel(sIdx).Parent = obj.Grid;
                obj.SectionLabel(sIdx).Layout.Row = 2;
                obj.SectionLabel(sIdx).Layout.Column = sIdx;
                obj.SectionLabel(sIdx).BackgroundColor = obj.BackgroundColor;
                obj.SectionLabel(sIdx).FontSize = 10;
                obj.SectionLabel(sIdx).HorizontalAlignment = 'center';
                obj.SectionLabel(sIdx).Text = upper(obj.Section(sIdx).Title);

                % Update each button
                obj.SectionButton(sIdx).Parent = obj.Grid;
                obj.SectionButton(sIdx).Layout.Row = [1 2];
                obj.SectionButton(sIdx).Layout.Column = sIdx;
                obj.SectionButton(sIdx).BackgroundColor = obj.BackgroundColor;
                obj.SectionButton(sIdx).FontSize = 10;
                obj.SectionButton(sIdx).IconAlignment = 'bottom';
                obj.SectionButton(sIdx).WordWrap = 'on';
                obj.SectionButton(sIdx).Visible = 'off';
                obj.SectionButton(sIdx).ButtonPushedFcn = @(h,e)obj.onPanelButtonPushed(e);
                obj.SectionButton(sIdx).Text = upper(obj.Section(sIdx).Title);

                % Update each section
                if ~isequal(obj.Section(sIdx).Parent, obj.Grid)
                    obj.Section(sIdx).Parent = obj.Grid;
                end
                obj.Section(sIdx).Layout.Row = 1;
                obj.Section(sIdx).Layout.Column = sIdx;

            end %for

            % Remove any extra labels/buttons
            numLabels = numel(obj.SectionLabel);
            if numLabels > numSections
                idxRemove = (numSections+1):numLabels;
                delete(obj.SectionLabel(idxRemove));
                obj.SectionLabel(idxRemove) = [];
                delete(obj.SectionButton(idxRemove));
                obj.SectionButton(idxRemove) = [];
            end

            % Unparent removed components
            isHorizontalSection = get(obj.Grid.Children,'Type') == ...
                lower("wt.toolbar.HorizontalSection");
            oldSections = obj.Grid.Children(isHorizontalSection);
            removedSections = setdiff(oldSections, obj.Section);
            set(removedSections,'Parent',[]);

            % Update component style lists
            obj.TitleColorableComponents = obj.SectionLabel;
            obj.BackgroundColorableComponents = [
                obj.Section
                obj.SectionButton
                obj.SectionLabel
                obj.DummySection
                ];
            obj.FontStyledComponents = [
                obj.Section
                obj.SectionButton
                ];

            % Update listeners
            obj.ButtonPushedListener = event.listener(obj.Section,...
                'ButtonPushed',@(h,e)obj.onButtonPushed(e));
            obj.SectionChangedListener = event.listener(obj.Section,...
                'PropertyChanged',@(h,e)obj.update());

            % Update the layout
            obj.updateLayout();

        end %function


        function propGroups = getPropertyGroups(obj)
            % Override the ComponentContainer GetPropertyGroups with newly
            % customiziable mixin. This can probably also be specific to each control.

            propGroups = getPropertyGroups@wt.mixin.PropertyViewable(obj);

        end
        

        function updateLayout(obj)
            % Dynamically configure the toolbar based on space

            % Return if no sections to show
            if isempty(obj.Section)
                return
            end

            % How many sections
            numSections = numel(obj.Section);

            % Total available width to place sections
            pos = getpixelposition(obj);
            wAvail = pos(3);

            % Subtract off needed spacing between components
            spacingNeeded = numSections - 1;
            wAvail = wAvail - spacingNeeded;

            % Calculate the width options
            sectionWidth = [obj.Section.TotalWidth];
            minSectionWidth = [obj.Section.MinimizedWidth];
            widthOptions = cumsum(sectionWidth) + ...
                flip( cumsum([0 minSectionWidth(1:end-1)]) ) ;

            % Which panels can be fully shown?
            isFullyShown = widthOptions <= wAvail;

            % Loop on each obj.Section
            for idx = 1:numel(obj.Section)

                % If unparented, fix it
                if ~isequal(obj.Section(idx).Parent, obj.Grid)
                    obj.Section(idx).Parent = obj.Grid;
                    obj.Section(idx).Layout.Row = 1;
                    obj.Section(idx).Layout.Column = idx;
                end

                % Toggle visibilities
                obj.Section(idx).Visible = isFullyShown(idx);
                if idx <= numel(obj.SectionLabel)
                    obj.SectionLabel(idx).Visible = isFullyShown(idx);
                    obj.SectionButton(idx).Visible = ~isFullyShown(idx);
                end

            end %for idx = 1:numel(section)

            % Update the grid column widths
            sectionWidth(~isFullyShown) = minSectionWidth(~isFullyShown);
            obj.Grid.ColumnWidth = [num2cell(sectionWidth),'1x'];

            % Place the dummy section
            obj.DummySection.Layout.Column = numSections + 1;

        end %function


        function updateBackgroundColorableComponents(obj)

            % Update button icons
            obj.updateButtonIcons();

            % Override the default, not setting the Grid background
            hasProp = isprop(obj.BackgroundColorableComponents,'BackgroundColor');
            set(obj.BackgroundColorableComponents(hasProp),...
                "BackgroundColor",obj.BackgroundColor);

        end %function


        function updateFontStyledComponents(obj,varargin)

            % Update button icons
            obj.updateButtonIcons();

            % Call the superclass method
            obj.updateFontStyledComponents@wt.mixin.FontStyled(varargin{:});

        end %function


        function updateButtonIcons(obj)
            % Color the down arrow on buttons the same as font color

            % Has foreground or background color changed?
            if ~isempty(obj.SectionButton) && ( ...
                    ~isequal(obj.SectionButton(end).FontColor, obj.FontColor) || ...
                    ~isequal(obj.DummySection.BackgroundColor, obj.BackgroundColor) )

                % Create the button icon
                icon = cell(1,3);
                for cIdx = 1:3
                    icon{cIdx} = ~obj.BUTTON_MASK * obj.BackgroundColor(cIdx);
                    icon{cIdx}(obj.BUTTON_MASK) = obj.FontColor(cIdx);
                end
                icon = cat(3,icon{:});

                % Update the icon on each button
                set(obj.SectionButton,'Icon',icon);

            end %if

        end %function


        function onButtonPushed(obj,evt)

            % Close any open sections
            obj.updateLayout();

            % Trigger event
            notify(obj,"ButtonPushed",evt);

        end %function

    end %methods



    %% Private methods
    methods (Access = private)


        function onPanelButtonPushed(obj,e)

            % Which panel?
            sectionButton = e.Source;
            isActivatedSection = obj.SectionButton == sectionButton;
            section = obj.Section(isActivatedSection);

            % Get the status of section panels to change
            isOpenSection = obj.SectionIsOpen;
            isDeactivating = isOpenSection & ~isActivatedSection;
            isActivating = ~isOpenSection & isActivatedSection;

            % Is the pushed section already open?
            if any(isOpenSection(isActivatedSection))
                % Instead update the layout to closse it
                obj.updateLayout();
                return
            end

            % Close any deactivated sections
            if any(isDeactivating)
                set(obj.Section(isDeactivating), 'Visible', 'off');
            end

            % Activate the specified section
            if any(isActivating)

                % Where are things located now?
                bPos = getpixelposition(sectionButton, true);
                fig = ancestor(obj,'figure');
                wt.utility.fastSet(fig,"Units","pixels");
                figureWidth = fig.Position(3);

                % Where should the panel go?
                panelX = bPos(1);
                panelWidth = section.TotalWidth;
                panelHeight = bPos(4) - obj.Grid.RowHeight{2} - obj.Grid.RowSpacing;
                panelY = bPos(2) - panelHeight - 1;

                % Adjust panel X position if needed
                panelRightEdge = panelX + panelWidth;
                if panelRightEdge > figureWidth
                    panelX = figureWidth - panelWidth;
                end

                % Now, position and show the panel as a dropdown
                panelPos = [panelX panelY panelWidth panelHeight];

                % Turn it on
                set(obj.Section(isActivating), 'Parent', fig);
                set(obj.Section(isActivating), 'Position', panelPos);
                set(obj.Section(isActivating), 'Visible', 'on');

            end %if any(isActivating)

        end %function

    end %methods


    %% Accessors
    methods

        function value = get.DividerColor(obj)
            value = obj.Grid.BackgroundColor;
        end
        function set.DividerColor(obj,value)
            obj.Grid.BackgroundColor = value;
        end

        function value = get.SectionIsOpen(obj)
            value = false(size(obj.Section));
            for idx = 1:numel(obj.Section)
                value(idx) = ~isequal(obj.Section(idx).Parent, obj.Grid);
            end
        end

    end %methods


end % classdef


%% Helper Functions

function mask = sectionButtonIconMask()

mask = logical([
    0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0
    1 1 1 1 1 1 1 1 1
    0 1 1 1 1 1 1 1 0
    0 0 1 1 1 1 1 0 0
    0 0 0 1 1 1 0 0 0
    0 0 0 0 1 0 0 0 0
    ]);

end %function