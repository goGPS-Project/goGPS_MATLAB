classdef ListSelectorTwoPane < matlab.ui.componentcontainer.ComponentContainer & ...
        wt.mixin.Enableable &...
        wt.mixin.FontStyled & wt.mixin.ButtonColorable &...
        wt.mixin.FieldColorable & wt.mixin.BackgroundColorable & ...
        wt.mixin.PropertyViewable
    
    % Dual lists where selected items are moved from left to right

    % Copyright 2020-2023 The MathWorks Inc.


    %% Events
    events (HasCallbackProperty, NotifyAccess = protected)

        % Triggered when a button is pushed
        ButtonPushed

        % Triggered when the value of the list selection changes
        ValueChanged

        % Triggered when the highlighted value changes
        HighlightedValueChanged

    end %events



    %% Public properties
    properties (AbortSet)

        % List of items to add to the list
        Items (1,:) string = ["Item 1" "Item 2" "Item 3" "Item 4"]

        % Data associated with  items (optional)
        ItemsData (1,:)

        % Indicates whether to allow sort controls
        Sortable (1,1) matlab.lang.OnOffSwitchState = true

    end %properties


    properties (AbortSet, Dependent)

        % Indices of displayed items that are currently added to the list
        SelectedIndex (1,:)

        % The current selection
        Value (1,:)

        % The current highlighted selection
        HighlightedValue (1,:)

    end %properties


    properties (AbortSet, Dependent, UsedInUpdate = false)

        % Width of the buttons
        ButtonWidth

    end %properties



    %% Read-Only properties
    properties (SetAccess = private)

        % Additional user buttons
        UserButtons wt.ButtonGrid

        % The list sorting buttons
        ListButtons wt.ButtonGrid

    end %properties



    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)
        
        % The left listbox control
        LeftList (1,1) matlab.ui.control.ListBox

        % Grid
        Grid (1,1) matlab.ui.container.GridLayout

        % The right listbox control
        RightList (1,1) matlab.ui.control.ListBox

        % Listen to button pushes in sections
        ButtonPushedListener event.listener

    end %properties



    %% Protected methods
    methods (Access = protected)

        function setup(obj)

            % Set default size
            obj.Position(3:4) = [200 160];

            % Construct Grid Layout to Manage Building Blocks
            obj.Grid = uigridlayout(obj);
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.RowSpacing = 2;
            obj.Grid.ColumnSpacing = 2;
            obj.Grid.Padding = 0;

            % Configure grid
            obj.Grid.Padding = 3;
            obj.Grid.ColumnWidth = {'1x',28,'1x'};
            obj.Grid.RowHeight = {'fit','1x'};

            % Create the left Listbox
            obj.LeftList = uilistbox(obj.Grid);
            obj.LeftList.Multiselect = true;
            obj.LeftList.Layout.Column = 1;
            obj.LeftList.Layout.Row = [1 2];
            obj.LeftList.ValueChangedFcn = @(h,e)obj.onLeftSelectionChanged(e);

            % Create the list buttons
            obj.ListButtons = wt.ButtonGrid(obj.Grid);
            obj.ListButtons.Layout.Column = 2;
            obj.ListButtons.Layout.Row = 1;
            obj.ListButtons.Orientation = "vertical";

            % Create an additional button grid for custom buttons
            obj.UserButtons = wt.ButtonGrid(obj.Grid,"Icon",[]);
            obj.UserButtons.Layout.Column = 2;
            obj.UserButtons.Layout.Row = 2;
            obj.UserButtons.Orientation = "vertical";

            % Create the right RightList
            obj.RightList = uilistbox(obj.Grid);
            obj.RightList.Multiselect = true;
            obj.RightList.Layout.Column = 3;
            obj.RightList.Layout.Row = [1 2];
            obj.RightList.ValueChangedFcn = @(h,e)obj.onRightSelectionChanged(e);

            % Update listeners
            obj.ButtonPushedListener = event.listener(...
                [obj.ListButtons obj.UserButtons],...
                'ButtonPushed',@(h,e)obj.onButtonPushed(e) );

            % Update the internal component lists
            obj.BackgroundColorableComponents = [obj.ListButtons, obj.UserButtons obj.Grid];
            obj.FontStyledComponents = [obj.RightList, obj.UserButtons, ...
                obj.ListButtons, obj.LeftList];
            obj.EnableableComponents = [obj.RightList, obj.UserButtons, ...
                obj.ListButtons, obj.LeftList];
            obj.ButtonColorableComponents = [obj.UserButtons obj.ListButtons];
            obj.FieldColorableComponents = [obj.RightList, obj.LeftList];

        end %function


        function update(obj)

            % What is selected?
            selIdx = obj.SelectedIndex;

            % Is the list sortable?
            if obj.Sortable
                obj.ListButtons.Icon = ["right_24.png", "left_24.png", "up_24.png", "down_24.png"];
                obj.ListButtons.ButtonTag = ["Add", "Remove", "Up", "Down"];
                obj.ListButtons.ButtonHeight = {28 28 28 28};
            else
                obj.ListButtons.Icon = ["right_24.png", "left_24.png"];
                obj.ListButtons.ButtonTag = ["Add", "Remove"];
                obj.ListButtons.ButtonHeight = {28 28};
            end

            % Update the list
            obj.RightList.Items = obj.Items(selIdx);
            obj.RightList.ItemsData = selIdx;

            % Update button enable states
            obj.updateEnables();

            % Update the list of choices
            itemIds = 1:numel(obj.Items);
            isNotSelected = ~ismember(itemIds, obj.SelectedIndex);
            obj.LeftList.Items = obj.Items(isNotSelected);
            obj.LeftList.ItemsData = itemIds(isNotSelected);

            % Update button enable states
            obj.updateEnables();

        end %function

        function propGroups = getPropertyGroups(obj)
            % Override the ComponentContainer GetPropertyGroups with newly
            % customiziable mixin. This can probably also be specific to each control.

            propGroups = getPropertyGroups@wt.mixin.PropertyViewable(obj);

        end %function


        function updateEnables(obj)

            % Button enables
            if obj.Enable

                % What is selected?
                selIdx = obj.SelectedIndex;

                % Highlighted selection in list?
                hiliteIdx = obj.getListBoxSelectedIndex();

                % How many items selected into list
                numRows = numel(selIdx);
                numHilite = numel(hiliteIdx);

                % Toggle button enables
                obj.ListButtons.ButtonEnable = [
                    numel(selIdx) < numel(obj.Items) %Add Button
                    ~isempty(hiliteIdx) % Delete Button
                    numHilite && ( hiliteIdx(end) > numHilite ) %Up Button
                    numHilite && ( hiliteIdx(1) <= (numRows - numHilite) ) %Down Button
                    ];

                obj.ListButtons.ButtonEnable(1) = ...
                    ~isempty(obj.LeftList.Value);

            end %if obj.Enable

        end %function


        function onLeftSelectionChanged(obj,~)

            % Update button enable states
            obj.updateEnables();

        end %function


        function onRightSelectionChanged(obj,evt)

            % Get the new and old values
            if isempty(obj.ItemsData)
                itemsData = obj.Items;
            else
                itemsData = obj.ItemsData;
            end

            if isempty(evt.PreviousValue)
                oldValue = itemsData([]);
            else
                oldValue = itemsData(evt.PreviousValue);
            end

            if isempty(evt.Value)
                newValue = itemsData([]);
            else
                newValue = itemsData(evt.Value);
            end

            % Update button enable states
            obj.updateEnables();

            % Trigger event
            evtOut = wt.eventdata.ValueChangedData(newValue, oldValue);
            notify(obj,"HighlightedValueChanged",evtOut);

        end %function


        function onButtonPushed(obj,evt)

            % Which button?
            switch evt.Tag

                case 'Add'

                    % Get the original value
                    oldValue = obj.Value;

                    % Update the selection
                    newSelIdx = [obj.SelectedIndex obj.LeftList.Value];
                    obj.SelectedIndex = newSelIdx;

                    % Trigger event
                    evtOut = wt.eventdata.ValueChangedData(obj.Value, oldValue);
                    notify(obj,"ValueChanged",evtOut);

                case 'Remove'
                    obj.removeListBoxSelection();

                case 'Up'
                    obj.shiftListBoxIndex(-1);

                case 'Down'
                    obj.shiftListBoxIndex(1);

                otherwise
                    % Trigger event for user buttons
                    notify(obj,"ButtonPushed",evt);

            end %switch

            % Request update
            obj.update();

        end %function


        function selIdx = getListBoxSelectedIndex(obj)
            % Get the current selected row indices in the listbox

            warnState = warning('off','MATLAB:structOnObject');
            s = struct(obj.RightList);
            warning(warnState);
            selIdx = s.SelectedIndex;
            if isequal(selIdx, -1)
                selIdx = [];
            end

        end %function


        function removeListBoxSelection(obj)
            % Removes the currently selected items from the listbox

            % What's currently selected?
            idxSel = obj.getListBoxSelectedIndex();

            % Is there something to remove?
            if ~isempty(idxSel)

                % Get the original value
                oldValue = obj.Value;

                % Remove it
                obj.RightList.Items(idxSel) = [];
                obj.RightList.ItemsData(idxSel) = [];

                % Trigger event
                evtOut = wt.eventdata.ValueChangedData(obj.Value, oldValue);
                notify(obj,"ValueChanged",evtOut);

            end %if

        end %function


        function shiftListBoxIndex(obj, shift)
            % Shift selected items up/down within a listbox
            % This assumes ItemsData contains unique values

            % What is the current order?
            selIdx = obj.getListBoxSelectedIndex();

            % Make indices to all items as they are now
            idxNew = 1:numel(obj.RightList.Items);
            idxOld = idxNew;

            % Find the last stable item that doesn't move
            [~,idxStable] = setdiff(idxNew, selIdx, 'stable');
            if ~isempty(idxStable)
                idxFirstStable = idxStable(1);
                idxLastStable = idxStable(end);
            else
                idxFirstStable = inf;
                idxLastStable = 0;
            end

            % Which way do we loop?
            if shift > 0 %Shift to end

                for idxToMove = numel(selIdx):-1:1

                    % Calculate if there's room to move this item
                    idxThisBefore = selIdx(idxToMove);
                    thisShift = max( min(idxLastStable-idxThisBefore, shift), 0 );

                    % Where does this item move from/to
                    idxThisAfter = idxThisBefore + thisShift;

                    % Where do other items move from/to
                    idxOthersBefore = selIdx(idxToMove)+1:1:idxThisAfter;
                    idxOthersAfter = idxOthersBefore - thisShift;

                    % Move the items
                    idxNew([idxThisAfter idxOthersAfter]) = idxNew([idxThisBefore idxOthersBefore]);

                end

            elseif shift < 0 %Shift to start

                for idxToMove = 1:numel(selIdx)

                    % Calculate if there's room to move this item
                    idxThisBefore = selIdx(idxToMove);
                    thisShift = min( max(idxFirstStable-idxThisBefore, shift), 0 );

                    % Where does this item move from/to
                    idxThisAfter = idxThisBefore + thisShift;

                    % Where do other items move from/to
                    idxOthersBefore = idxThisAfter:1:selIdx(idxToMove)-1;
                    idxOthersAfter = idxOthersBefore - thisShift;

                    % Move the items
                    idxNew([idxThisAfter idxOthersAfter]) = idxNew([idxThisBefore idxOthersBefore]);

                end

            end %if shift > 0

            % Was a change made?
            if ~isequal(idxOld, idxNew)

                % Get the original value
                oldValue = obj.Value;

                % Make the shift
                obj.RightList.Items = obj.RightList.Items(idxNew);
                obj.RightList.ItemsData = obj.RightList.ItemsData(idxNew);

                % Trigger event
                evtOut = wt.eventdata.ValueChangedData(obj.Value, oldValue);
                notify(obj,"ValueChanged",evtOut);

            end %if

        end %function

    end %methods



    %% Accessors
    methods

        function value = get.SelectedIndex(obj)
            value = obj.RightList.ItemsData;
        end
        function set.SelectedIndex(obj,value)
            if ~obj.Sortable
                value = sort(value);
            end
            obj.RightList.Items = obj.Items(value);
            obj.RightList.ItemsData = value;
        end

        function value = get.Value(obj)
            if isempty(obj.ItemsData)
                value = obj.Items(:,obj.RightList.ItemsData);
            else
                value = obj.ItemsData(:,obj.RightList.ItemsData);
            end
        end
        function set.Value(obj,value)
            if isempty(value)
                obj.SelectedIndex = [];
            else
                if isempty(obj.ItemsData)
                    [tf, selIdx] = ismember(value, obj.Items);
                else
                    [tf, selIdx] = ismember(value, obj.ItemsData);
                end
                if ~all(tf)
                    warning("widgets:ListSelector:InvalidValue",...
                        "Attempt to set an invalid Value to the list.")
                    selIdx(~tf) = [];
                end
                obj.SelectedIndex = selIdx;
            end
        end

        function value = get.HighlightedValue(obj)
            selIdx = obj.RightList.Value;
            if isempty(selIdx) || ~isnumeric(selIdx)
                selIdx = [];
            end
            if isempty(obj.ItemsData)
                value = obj.Items(:,selIdx);
            else
                value = obj.ItemsData(:,selIdx);
            end
        end        
        function set.HighlightedValue(obj,value)
            if isempty(value)
                return;
            end
            if isempty(obj.ItemsData)
                [~, obj.RightList.Value] = ismember(value, obj.Items);
            else
                [~, obj.RightList.Value] = ismember(value, obj.ItemsData);
            end
        end

        function value = get.ButtonWidth(obj)
            value = obj.Grid.ColumnWidth{2};
        end
        function set.ButtonWidth(obj,value)
            obj.Grid.ColumnWidth{2} = value;
        end

    end %methods

end % classdef