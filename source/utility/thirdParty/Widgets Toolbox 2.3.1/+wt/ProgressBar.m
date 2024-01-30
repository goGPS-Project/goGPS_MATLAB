classdef ProgressBar < matlab.ui.componentcontainer.ComponentContainer & ...
        wt.mixin.BackgroundColorable & ...
        wt.mixin.FontStyled & wt.mixin.PropertyViewable

    % A progress bar with status and cancel button

    % Copyright 2020-2023 The MathWorks Inc.


    %% Events
    events (HasCallbackProperty, NotifyAccess = protected)
        CancelPressed
    end


    %% Properties
    properties (AbortSet)

        % Shows time remaining during progress
        ShowTimeRemaining  (1,1) matlab.lang.OnOffSwitchState = true

        % Enables a Cancel button during progress
        ShowCancel  (1,1) matlab.lang.OnOffSwitchState = false

        % Indicates that progress is counting
        Indeterminate  (1,1) matlab.lang.OnOffSwitchState = false

    end %properties


    properties (Dependent, UsedInUpdate = false)

        % Bar color
        BarColor (1,3) double {wt.validators.mustBeBetweenZeroAndOne}

    end %properties


    % These properties are internally set by methods, and do not trigger
    % the update method
    properties (SetAccess = protected, UsedInUpdate = false)

        % Current progress value
        Value (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(Value,1)} = 0

        % Current status text displayed
        StatusText (1,1) string

        % Indicates that progress is counting
        Running  (1,1) matlab.lang.OnOffSwitchState = false

        % Indicates that Cancel button was pressed
        CancelRequested  (1,1) matlab.lang.OnOffSwitchState = false

        % Start time
        StartTime (1,1) datetime = NaT

    end %properties


    properties (Dependent, SetAccess = protected, UsedInUpdate = false)

        % Elapsed time
        ElapsedTime (1,1) datetime

        % Remaining time
        RemainingTime (1,1) duration

    end %properties



    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)
        
        % Progress panel
        ProgressPanel (1,1) matlab.ui.container.Panel

        % Grid
        Grid (1,1) matlab.ui.container.GridLayout

        % Indeterminate bar
        IndeterminateBar (1,1) matlab.ui.control.Image

        % Cancel button
        CancelButton (1,1) matlab.ui.control.Image

        % Status text
        StatusTextLabel (1,1) matlab.ui.control.Label

        % Remaining Time text
        RemTimeLabel (1,1) matlab.ui.control.Label

    end %properties



    %% Public Methods
    methods

        function startProgress(obj, text)
            % Begin the progress bar from zero

            arguments
                obj (1,1) wt.ProgressBar
                text (1,1) string = ""
            end

            obj.StartTime = datetime();
            obj.CancelRequested = false;
            obj.Running = true;
            obj.Value = 0;
            obj.StatusText = text;

            obj.update();

            % Process updates immediately
            drawnow

        end %function


        function setProgress(obj, value, text)
            % Set the progress bar and status text

            arguments
                obj (1,1) wt.ProgressBar
                value (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(value,1)} = 0
                text (1,1) string = ""
            end

            obj.Running = true;
            obj.Value = value;
            if ~obj.CancelRequested
                obj.StatusText = text;
            end

            obj.update();

            % Process updates immediately
            drawnow('limitrate')

        end %function


        function finishProgress(obj, text)
            % Finish the progress bar

            arguments
                obj (1,1) wt.ProgressBar
                text (1,1) string = ""
            end

            obj.StartTime = NaT;
            obj.CancelRequested = false;
            obj.Running = false;
            obj.Value = 0;
            obj.StatusText = text;

            obj.update();

            % Process updates immediately
            drawnow

        end %function


        function cancel(obj)
            % Triggered on button pushed

            arguments
                obj (1,1) wt.ProgressBar
            end

            % Set cancel state
            obj.CancelRequested = true;
            obj.StatusText = "Canceling...";

            obj.update();

            % Process updates immediately
            drawnow

            % Trigger event
            evt = wt.eventdata.ButtonPushedData(obj.CancelButton, "Cancel");
            notify(obj,"CancelPressed",evt);

        end %function


        function setStatusText(obj, text)
            % Set the status text

            obj.StatusText = text;

            obj.update();

            % Process updates immediately
            drawnow('limitrate')

        end %function


        function demo(obj)
            % Demonstrates the progress bar moving

            obj.startProgress("Starting");
            pause(0.5)
            n = 20;
            for idx = 1:n
                if obj.CancelRequested
                    break
                end
                pause(0.2);
                obj.setProgress(idx/n, "Running");
            end
            pause(0.5)
            obj.finishProgress("Finished");

        end %function

    end %methods



    %% Protected methods
    methods (Access = protected)

        function setup(obj)

            % Set default size
            obj.Position(3:4) = [200 30];

            % Construct Grid Layout to Manage Building Blocks
            obj.Grid = uigridlayout(obj);
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.RowSpacing = 2;
            obj.Grid.ColumnSpacing = 2;
            obj.Grid.Padding = 2;

            % Configure grid
            obj.Grid.Padding = 3;
            obj.Grid.ColumnWidth = {'0x','1x',25};
            obj.Grid.RowHeight = {'1x'};

            % Progress panel indicator
            obj.ProgressPanel = matlab.ui.container.Panel(...
                "Parent",obj.Grid,...
                "Units","pixels",...
                "BackgroundColor",[.5 .7 1],...
                "BorderType","none");
            obj.ProgressPanel.Layout.Column = 1;
            obj.ProgressPanel.Layout.Row = 1;

            % Indeterminate progress indicator
            obj.IndeterminateBar = matlab.ui.control.Image(...
                "Parent",obj.Grid,...
                "ImageSource","progress_bar_200_10.gif",...
                "ScaleMethod","stretch",...
                "HorizontalAlignment","center",...
                "Visible",false,...
                "ImageClickedFcn",@(h,e)obj.cancel());
            obj.IndeterminateBar.Layout.Column = [1 2];
            obj.IndeterminateBar.Layout.Row = 1;

            % Status text indicator
            obj.StatusTextLabel = matlab.ui.control.Label(...
                "Parent",obj.Grid,...
                "Text","",...
                "HorizontalAlignment","left");
            obj.StatusTextLabel.Layout.Column = [1 2];
            obj.StatusTextLabel.Layout.Row = 1;

            % Remaining time indicator
            obj.RemTimeLabel = matlab.ui.control.Label(...
                "Parent",obj.Grid,...
                "Text","",...
                "HorizontalAlignment","right");
            obj.RemTimeLabel.Layout.Column = [1 2];
            obj.RemTimeLabel.Layout.Row = 1;

            % Cancel Button
            obj.CancelButton = matlab.ui.control.Image(...
                "Parent",obj.Grid,...
                "ImageSource","stop_24.png",...
                "ScaleMethod","scaledown",...
                "HorizontalAlignment","right",...
                "ImageClickedFcn",@(h,e)obj.cancel());
            obj.CancelButton.Layout.Column = 3;
            obj.CancelButton.Layout.Row = 1;

            % Update the internal component lists
            obj.FontStyledComponents = [obj.StatusTextLabel obj.RemTimeLabel];
            obj.BackgroundColorableComponents = obj.Grid;

        end %function


        function update(obj)

            % Toggle visibilities
            remTimeVisible = obj.ShowTimeRemaining;
            wt.utility.fastSet(obj.RemTimeLabel,"Visible",remTimeVisible);

            % Show/hide the indeterminate bar
            indBarVisible = obj.Indeterminate && obj.Running;
            wt.utility.fastSet(obj.IndeterminateBar,"Visible",indBarVisible);

            % Update the progress indication by adjusting the beneath grid
            colWidths = [cellstr([obj.Value 1-obj.Value] + "x") {0}];

            % Update the status text
            if strlength(obj.StatusText)
                statusText = "  " + obj.StatusText;
            else
                statusText = "";
            end
            wt.utility.fastSet(obj.StatusTextLabel,"Text",statusText);

            % Update remaining time
            if obj.ShowTimeRemaining && ~isinf(obj.RemainingTime)
                timeText = string(obj.RemainingTime);
            else
                timeText = "";
            end
            wt.utility.fastSet(obj.RemTimeLabel,"Text",timeText);

            % Update cancel button
            cancelVisible = obj.ShowCancel && obj.Running;
            cancelColor = "none";
            if cancelVisible
                colWidths{3} = 25;
                if obj.CancelRequested
                    cancelColor = "red";
                end
                wt.utility.fastSet(obj.CancelButton,"BackgroundColor",cancelColor);
            end
            wt.utility.fastSet(obj.CancelButton,"Visible",cancelVisible);

            % Update the grid
            obj.Grid.ColumnWidth = colWidths;

        end %function

        
        function propGroups = getPropertyGroups(obj)
            % Override the ComponentContainer GetPropertyGroups with newly
            % customiziable mixin. This can probably also be specific to each control.

            propGroups = getPropertyGroups@wt.mixin.PropertyViewable(obj);

        end % function

    end %methods


    %% Accessors
    methods

        function value = get.ElapsedTime(obj)
            value = datetime - obj.StartTime;
            value = i_formatDuration(value);
        end

        function value = get.RemainingTime(obj)
            value = (1-obj.Value)/obj.Value * obj.ElapsedTime;
            value = i_formatDuration(value);
        end

        function value = get.BarColor(obj)
            value = obj.ProgressPanel.BackgroundColor;
        end
        function set.BarColor(obj,value)
            obj.ProgressPanel.BackgroundColor = value;
        end
    end % methods


end % classdef



%% Helper functions
function value = i_formatDuration(value)

if ismissing(value)
    seconds(0);
end

if value < hours(1)
    value.Format = "mm:ss";
end

end %function