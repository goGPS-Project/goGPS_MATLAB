classdef DatetimeSelector < matlab.ui.componentcontainer.ComponentContainer & ...
        wt.mixin.BackgroundColorable & ...
        wt.mixin.Enableable & wt.mixin.FontStyled & wt.mixin.FieldColorable & ...
        wt.mixin.PropertyViewable

    % Date and time selection control

    % Copyright 2020-2023 The MathWorks Inc.


    %% Events
    events (HasCallbackProperty, NotifyAccess = protected)

        % Triggered on value changed, has companion callback
        ValueChanged

    end %events


    %% Public properties
    properties

        % The current value shown
        Value (1,1) datetime

    end %properties


    properties (AbortSet)

        % The time format
        ShowAMPM (1,1) matlab.lang.OnOffSwitchState = 'on'

        % Show seconds or not
        ShowSeconds (1,1) matlab.lang.OnOffSwitchState = 'off'

        % Show timezone
        ShowTimeZone (1,1) matlab.lang.OnOffSwitchState = 'off'
        %         ShowTimeZone (1,1) wt.enum.TimeZoneDisplayState = ...
        %             wt.enum.TimeZoneDisplayState.none

    end %properties


    properties (AbortSet, Dependent, UsedInUpdate = false)

        % Format for display
        DateFormat

        % Disabled days of week
        DisabledDaysOfWeek

        % DisabledDates
        DisabledDates

        % Limits on the date
        Limits

    end %properties



    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)

        % Button
        DateControl (1,1) matlab.ui.control.DatePicker

        % Grid
        Grid (1,1) matlab.ui.container.GridLayout

        % Hour control
        HourControl (1,1) matlab.ui.control.Spinner

        % Minute control
        MinuteControl (1,1) matlab.ui.control.Spinner

        % Second control
        SecondControl (1,1) matlab.ui.control.Spinner

        % AmPm control
        AmPmControl (1,1) matlab.ui.control.DropDown

        % TimeZone control
        TimeZoneControl (1,1) matlab.ui.control.DropDown

    end %properties


    properties (Constant, Access = protected)

        % Time zone info
        TimeZoneInfo table = getTimeZoneInfo();

    end %properties



    %% Protected methods
    methods (Access = protected)

        function setup(obj)

            % Adjust default size
            obj.Position(3:4) = [270 25];

            % Set default date to today
            obj.Value = datetime("today",...
                "TimeZone","local",...
                "Format","dd-MMM-uuuu hh:mm aa");

            % Construct Default Grid Layout to Manage Building Blocks
            obj.Grid = uigridlayout(obj);
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.RowSpacing = 2;
            obj.Grid.ColumnSpacing = 2;
            obj.Grid.Padding = 0;

            % Configure Grid
            obj.Grid.ColumnWidth = {'9x',5,'4x','4x',0,0,0};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.ColumnSpacing = 0;

            % Create the date control
            obj.DateControl = matlab.ui.control.DatePicker(...
                "Parent",obj.Grid,...
                "ValueChangedFcn",@(h,e)obj.onDateEdited(e));

            % Spacer
            uicontainer(obj.Grid,'Visible','off');

            % Create the time controls
            % obj.HourControl = matlab.ui.control.NumericEditField(...
            obj.HourControl = uispinner(...
                "Parent",obj.Grid,...
                "Limits",[-1 24],...
                "HorizontalAlignment","center",...
                "ValueChangedFcn",@(h,e)obj.onTimeEdited(e));

            %obj.MinuteControl = matlab.ui.control.NumericEditField(...
            obj.MinuteControl = uispinner(...
                "Parent",obj.Grid,...
                "Limits",[-1 60],...
                "ValueDisplayFormat","%02.0f",...
                "HorizontalAlignment","center",...
                "ValueChangedFcn",@(h,e)obj.onTimeEdited(e));

            %obj.SecondControl = matlab.ui.control.NumericEditField(...
            obj.SecondControl = uispinner(...
                "Parent",obj.Grid,...
                "Limits",[-1 60],...
                "ValueDisplayFormat","%02.f",...
                "HorizontalAlignment","center",...
                "ValueChangedFcn",@(h,e)obj.onTimeEdited(e));

            obj.AmPmControl = matlab.ui.control.DropDown(...
                "Parent",obj.Grid,...
                "Items",["AM","PM"],...
                "ValueChangedFcn",@(h,e)obj.onTimeEdited(e));

            obj.TimeZoneControl = uidropdown(obj.Grid);
            obj.TimeZoneControl.ValueChangedFcn = @(h,e)obj.onTimeZoneChanged(e);

            % Update the internal component lists
            allFields = [
                obj.DateControl
                obj.HourControl
                obj.MinuteControl
                obj.SecondControl
                obj.AmPmControl
                obj.TimeZoneControl
                ];
            obj.FontStyledComponents = allFields;
            obj.FieldColorableComponents = allFields;
            obj.EnableableComponents = allFields;
            obj.BackgroundColorableComponents = obj.Grid;

        end %function


        function update(obj)

            % Get the value
            v = obj.Value;

            % Toggle visibilities
            obj.SecondControl.Visible = obj.ShowSeconds;
            obj.AmPmControl.Visible = obj.ShowAMPM;
            obj.TimeZoneControl.Visible = obj.ShowTimeZone;

            if obj.ShowSeconds
                obj.Grid.ColumnWidth{5} = '4x';
            else
                obj.Grid.ColumnWidth{5} = 0;
            end

            if obj.ShowAMPM
                obj.Grid.ColumnWidth{6} = '5x';
            else
                obj.Grid.ColumnWidth{6} = 0;
            end

            if obj.ShowTimeZone
                obj.Grid.ColumnWidth{7} = '8x';
            else
                obj.Grid.ColumnWidth{7} = 0;
            end

            % Update the date control
            valFormat = strtok(v.Format,' ');
            obj.DateControl.Value = v;
            obj.DateControl.DisplayFormat = valFormat;

            % Which hour format?
            if obj.ShowAMPM

                % Update the AM/PM
                if v.Hour >= 12
                    obj.AmPmControl.Value = "PM";
                else
                    obj.AmPmControl.Value = "AM";
                end

                % Update the hour control
                if v.Hour > 12
                    obj.HourControl.Value = v.Hour - 12;
                elseif v.Hour == 0
                    obj.HourControl.Value = 12;
                elseif isnan(v.Hour)
                    obj.HourControl.Value = 12;
                else
                    obj.HourControl.Value = v.Hour;
                end
                %obj.HourControl.Limits = [1 12];
                obj.HourControl.ValueDisplayFormat = "%2.0f";

            else

                % Update the hour control
                %obj.HourControl.Limits = [0 23];
                obj.HourControl.ValueDisplayFormat = "%02.0f";
                if isnan(v.Hour)
                    obj.HourControl.Value = 0;
                else
                    obj.HourControl.Value = v.Hour;
                end

            end %if obj.ShowAMPM

            % Update the minutes
            if isnan(v.Minute)
                obj.MinuteControl.Value = 0;
            else
                obj.MinuteControl.Value = v.Minute;
            end

            % Update the seconds
            if isnan(v.Second)
                obj.SecondControl.Value = 0;
            else
                obj.SecondControl.Value = v.Second;
            end

            % Prepare time zone list and value
            if obj.ShowTimeZone

                % Get the current value
                tzValue = obj.Value.TimeZone;

                % Prepare the list selections
                if startsWith(tzValue,["+","-"])
                    % Current timezone is just a value with no region!

                    allOffsets = unique(obj.TimeZoneInfo.OffsetName,'stable');

                    value = tzoffset(obj.Value);
                    valueStr = string( value );
                    if value >= 0
                        valueStr = "+" + valueStr;
                    end

                    tzItems = allOffsets;
                    tzItemsData = allOffsets;
                    tzValue = valueStr;

                else
                    % Current timezone is a standard region name
                    % Or it is an empty string

                    tzItems = obj.TimeZoneInfo.CombinedName;
                    tzItemsData = obj.TimeZoneInfo.Name;
                    tzValue = obj.Value.TimeZone;

                end

                % Update time zone
                obj.TimeZoneControl.Items = tzItems;
                obj.TimeZoneControl.ItemsData = tzItemsData;
                obj.TimeZoneControl.Value = tzValue;

            end %if obj.ShowTimeZone

        end %function


        function propGroups = getPropertyGroups(obj)
            % Override the ComponentContainer GetPropertyGroups with newly
            % customiziable mixin. This can probably also be specific to each control.

            propGroups = getPropertyGroups@wt.mixin.PropertyViewable(obj);

        end %function


        function onDateEdited(obj,evt)
            % Triggered on edits

            % Get the new date
            newValue = evt.Source.Value;

            % Existing date (if missing, use today)
            dt = obj.Value;
            if ismissing(dt)
                dt = datetime("today",...
                "TimeZone","local",...
                "Format","dd-MMM-uuuu hh:mm aa");
            end

            % Set the date
            dt.Year = newValue.Year;
            dt.Month = newValue.Month;
            dt.Day = newValue.Day;

            % Update value
            obj.Value = dt;

            % Trigger event
            evtOut = wt.eventdata.ValueChangedData(dt);
            notify(obj,"ValueChanged",evtOut);

        end %function


        function onTimeEdited(obj,evt)
            % Triggered on edits

            % Get the new value
            newValue = evt.Source.Value;

            % Existing date (if missing, use today)
            dt = obj.Value;
            if ismissing(dt)
                dt = datetime("today",...
                "TimeZone","local",...
                "Format","dd-MMM-uuuu hh:mm aa");
            end

            % What was changed?
            switch evt.Source

                case obj.HourControl
                    if obj.ShowAMPM && dt.Hour == 0
                        dt.Hour = newValue - 12;
                    elseif obj.ShowAMPM && dt.Hour > 12
                        dt.Hour = newValue + 12;
                    else
                        dt.Hour = newValue;
                    end

                case obj.MinuteControl
                    dt.Minute = newValue;

                case obj.SecondControl
                    dt.Second = newValue;

                case obj.AmPmControl

                    if newValue == "AM" && dt.Hour > 11
                        dt.Hour = dt.Hour - 12;
                    elseif newValue == "PM" && dt.Hour < 12
                        dt.Hour = dt.Hour + 12;
                    end

            end %switch

            % Update value
            obj.Value = dt;

            % Trigger event
            evtOut = wt.eventdata.ValueChangedData(obj.Value);
            notify(obj,"ValueChanged",evtOut);

        end %function


        function onTimeZoneChanged(obj,evt)
            % Triggered on edits

            % Get the new value
            newValue = evt.Value;
            dt = obj.Value;
            dt.TimeZone = newValue;

            % Update value
            obj.Value = dt;

            % Trigger event
            evtOut = wt.eventdata.ValueChangedData(obj.Value);
            notify(obj,"ValueChanged",evtOut);

        end %function

    end % methods



    %% Accessors
    methods

        function value = get.DateFormat(obj)
            value = obj.DateControl.DisplayFormat;
        end
        function set.DateFormat(obj,value)
            obj.DateControl.DisplayFormat = value;
        end

        function value = get.Limits(obj)
            value = obj.DateControl.Limits;
        end
        function set.Limits(obj,value)
            obj.DateControl.Limits = value;
        end

        function value = get.DisabledDaysOfWeek(obj)
            value = obj.DateControl.DisabledDaysOfWeek;
        end
        function set.DisabledDaysOfWeek(obj,value)
            obj.DateControl.DisabledDaysOfWeek = value;
        end

        function value = get.DisabledDates(obj)
            value = obj.DateControl.DisabledDates;
        end
        function set.DisabledDates(obj,value)
            obj.DateControl.DisabledDates = value;
        end

    end %methods

end % classdef


%% Helper Functions
function timeZoneInfo = getTimeZoneInfo()
% Format a table of the necessary time zone info

% Get the timezones
timeZoneInfo = timezones();
timeZoneInfo.Name = string(timeZoneInfo.Name);

% Remove the Etc timezones
% These seem to have UTCOffset that is negative of what the name implies,
% causing issues
timeZoneInfo(timeZoneInfo.Area == "Etc", :) = [];

% Sort by offset from GMT
timeZoneInfo = sortrows(timeZoneInfo, "UTCOffset");

% Prepare the formatted offsets
offsetH = fix(timeZoneInfo.UTCOffset);
offsetM = mod(timeZoneInfo.UTCOffset, 1) * 60;

% Attach items to the table
timeZoneInfo.OffsetName = compose("%+03d:%02d", offsetH, offsetM);
timeZoneInfo.CombinedName = timeZoneInfo.OffsetName + "  " + timeZoneInfo.Name;

% Add a blank row to the top for empty timezone info
timeZoneInfo = vertcat({"","",NaN,NaN,"",""}, timeZoneInfo);

end %function