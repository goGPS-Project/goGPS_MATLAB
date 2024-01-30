classdef TaskStatusTable < matlab.ui.componentcontainer.ComponentContainer & ...
        wt.mixin.Enableable & wt.mixin.FontStyled & wt.mixin.Tooltipable & ...
        wt.mixin.ButtonColorable & wt.mixin.BackgroundColorable & ...
        wt.mixin.PropertyViewable

    % A table showing status of multiple tasks
    
    % Copyright 2020-2023 The MathWorks Inc.
    
    
    %% Public properties
    properties (AbortSet)
        
        % Name of each task
        Items (:,1) string = ["Step 1","Step 2","Step 3","Step 4",...
            "Step 5","Step 6","Step 7"]
        
        % Status of task
        Status (:,1) wt.enum.StatusState = [
            wt.enum.StatusState.complete
            wt.enum.StatusState.warning
            wt.enum.StatusState.error
            wt.enum.StatusState.info
            wt.enum.StatusState.running
            wt.enum.StatusState.unknown
            wt.enum.StatusState.none
            ];
        
        % Height of each task row
        RowHeight (1,1) double = 25
        
        % Selected task
        SelectedIndex double {mustBeScalarOrEmpty, mustBeFinite, mustBePositive} = 1
        
        % Selection color
        SelectionColor (1,3) double ...
            {wt.validators.mustBeBetweenZeroAndOne} = [.8 .8 1];
        
        % Enables button and status row to display
        ShowButtonRow (1,1) matlab.lang.OnOffSwitchState = true
        
        % Status message
        StatusMessage (1,1) string = "";
        
        % Enables the task forward button to proceed
        EnableForward (1,1) matlab.lang.OnOffSwitchState = true
        
        % Enables the task back button to proceed
        EnableBack (1,1) matlab.lang.OnOffSwitchState = true
        
    end %properties
    
    
    
    %% Events
    events (HasCallbackProperty, NotifyAccess = protected)
        
        % Triggered on button pushed, has companion callback
        ButtonPushed
        
    end %events
    
    
    
    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)
        
        % Grid for task items
        TaskGrid (1,1) matlab.ui.container.GridLayout

        % Main Grid for Layout
        Grid (1,1) matlab.ui.container.GridLayout
        
        % Task Labels
        Label (1,:) matlab.ui.control.Label
        
        % Edit control
        Icon (1,:) matlab.ui.control.Image
        
        % Back Button
        BackButton (1,1) matlab.ui.control.Button
        
        % Status Label
        StatusLabel (1,1) matlab.ui.control.Label
        
        % Forward Button
        ForwardButton (1,1) matlab.ui.control.Button
        
    end %properties
    
    
    
    %% Protected methods
    methods (Access = protected)
        
        function setup(obj)

            % Set default size
            obj.Position(3:4) = [100 180];
            
            % Construct Grid Layout to Manage Building Blocks
            obj.Grid = uigridlayout(obj);
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.RowSpacing = 2;
            obj.Grid.ColumnSpacing = 2;
            obj.Grid.Padding = 0;   
            
            % Configure Main Grid
            obj.Grid.ColumnWidth = {25,'1x',25};
            obj.Grid.RowHeight = {'1x',25};
            
            % Layout for task grid
            obj.TaskGrid = uigridlayout(obj.Grid, [7 2]);
            obj.TaskGrid.Padding = 2;
            obj.TaskGrid.ColumnSpacing = 5;
            obj.TaskGrid.RowSpacing = 2;
            obj.TaskGrid.Layout.Column = [1 3];
            obj.TaskGrid.ColumnWidth = {16,'1x'};
            obj.TaskGrid.RowHeight = {obj.RowHeight};
            obj.TaskGrid.Scrollable = true;
            
            % Create back button
            obj.BackButton = matlab.ui.control.Button(...
                "Parent",obj.Grid,...
                "Icon","left_24.png",...
                "Text","",...
                "ButtonPushedFcn",@(h,e)obj.onBackButtonPushed(e));
            obj.BackButton.Layout.Row = 2;
            obj.BackButton.Layout.Column = 1;
            
            % Create the status text
            obj.StatusLabel = matlab.ui.control.Label(...
                "Parent",obj.Grid,...
                "Text","",...
                "HorizontalAlignment","center");
            obj.StatusLabel.Layout.Row = 2;
            obj.StatusLabel.Layout.Column = 2;
            
            % Create next button
            obj.ForwardButton = matlab.ui.control.Button(...
                "Parent",obj.Grid,...
                "Icon","right_24.png",...
                "Text","",...
                "ButtonPushedFcn",@(h,e)obj.onForwardButtonPushed(e));
            obj.ForwardButton.Layout.Row = 2;
            obj.ForwardButton.Layout.Column = 3;
            
            % Update the internal component lists
            obj.BackgroundColorableComponents = [obj.TaskGrid obj.Grid];
            obj.ButtonColorableComponents = [obj.BackButton, obj.ForwardButton];
           
        end %function
        
        
        function update(obj)
            
            % How many tasks?
            numOld = numel(obj.Label);
            numNew = numel(obj.Items);
            
            % Update number of rows
            if numNew > numOld
                
                % Add rows
                for idx = (numOld+1):numNew
                    obj.Icon(idx) = uiimage(obj.TaskGrid,...
                        "ImageSource","",...
                        "ScaleMethod","fit");
                    obj.Label(idx) = uilabel(obj.TaskGrid,"Text","");
                end
                
            elseif numOld > numNew
                
                % Remove rows
                delete(obj.Label((numNew+1):end));
                delete(obj.Icon((numNew+1):end));
                obj.Label((numNew+1):end) = [];
                obj.Icon((numNew+1):end) = [];
                
            end %if numNew > numOld
            
            % Update the internal component lists
            if numOld ~= numNew
                obj.FontStyledComponents = [obj.Label, obj.StatusLabel];
                obj.EnableableComponents = [obj.Label, obj.Icon, obj.StatusLabel];
            end
                
            % Update the task names and icons
            status = obj.Status;
            imgFile = string(status) + "_16.png";
            imgFile(status=="running") = "running_16.gif";
            numImg = numel(obj.Icon);
            for idx = 1:numNew
                wt.utility.fastSet(obj.Label(idx),"Text",obj.Items(idx));
                if status(idx) == "none"
                    wt.utility.fastSet(obj.Icon(idx),"Visible","off");
                elseif numImg >= idx
                    if obj.Icon(idx).ImageSource ~= imgFile(idx)
                        if exist(imgFile(idx),'file')
                            obj.Icon(idx).ImageSource = imgFile(idx);
                        else
                            obj.Icon(idx).ImageSource = "";
                        end
                        wt.utility.fastSet(obj.Icon(idx),"Visible","on");
                    end
                else
                    wt.utility.fastSet(obj.Icon(idx),"Visible","off");
                end
            end
            
            % Update row heights
            wt.utility.fastSet(obj.TaskGrid, "RowHeight", repmat({obj.RowHeight},1,numNew));
            
            % Show/hide status and button row
            if obj.ShowButtonRow
                obj.Grid.RowHeight{2} = 25;
            else
                obj.Grid.RowHeight{2} = 0;
            end
            
            % Update status label
            wt.utility.fastSet(obj.StatusLabel,"Text",obj.StatusMessage);
            
            % Update selected task
            selLabel = obj.Label(obj.SelectedIndex);
            nonSelLabel = obj.Label(obj.Label ~= selLabel);
            wt.utility.fastSet(nonSelLabel,"BackgroundColor",'none');
            if ~isempty(selLabel)
                wt.utility.fastSet(selLabel,"BackgroundColor",obj.SelectionColor);
            end
            
            % Update the button enables
            obj.updateEnableableComponents();
            
        end %function
        
        
        function propGroups = getPropertyGroups(obj)
            % Override the ComponentContainer GetPropertyGroups with newly
            % customiziable mixin. This can probably also be specific to each control.

            propGroups = getPropertyGroups@wt.mixin.PropertyViewable(obj);

        end
        

        function updateEnableableComponents(obj)
            
            % Call superclass method first
            obj.updateEnableableComponents@wt.mixin.Enableable();
            
            % Can we go forward or back, based on selected task?
            canGoBack = obj.SelectedIndex > 1;
            canGoForward = obj.SelectedIndex < numel(obj.Items);
            backEnable = obj.Enable && obj.EnableBack && canGoBack;
            forwardEnable = obj.Enable && obj.EnableForward && canGoForward;
            
            % Update button enables
            wt.utility.fastSet(obj.BackButton,"Enable",backEnable);
            wt.utility.fastSet(obj.ForwardButton,"Enable",forwardEnable);
            
        end %function
        
        
        function onBackButtonPushed(obj,evt)
            % Triggered on button pushed
            
            % Go back a step
            if ~isempty(obj.SelectedIndex) && obj.SelectedIndex > 1
                obj.SelectedIndex = obj.SelectedIndex - 1;
            end
            
            % Trigger event
            evtOut = wt.eventdata.ButtonPushedData(evt.Source, "Back");
            notify(obj,"ButtonPushed",evtOut);
            
        end %function
        
        
        function onForwardButtonPushed(obj,evt)
            % Triggered on button pushed
            
            % Go forward a step
            if isempty(obj.SelectedIndex) 
                obj.SelectedIndex = 1;
            elseif obj.SelectedIndex < numel(obj.Items)
                obj.SelectedIndex = obj.SelectedIndex + 1;
            end
            
            % Trigger event
            evtOut = wt.eventdata.ButtonPushedData(evt.Source, "Forward");
            notify(obj,"ButtonPushed",evtOut);
            
        end %function
        
    end %methods
    
    
    %% Accessors
    methods
        
        function value = get.Status(obj)
            numItems = numel(obj.Items);
            value = obj.Status;
            value(end+1:numItems) = "none";
            value(numItems+1:end) = [];
        end
        
        function value = get.SelectedIndex(obj)
            value = obj.SelectedIndex;
            if ~isempty(value) && value > numel(obj.Items)
                value = [];
            end
        end
        
    end % methods
    
    
end % classdef