classdef PasswordField <  matlab.ui.componentcontainer.ComponentContainer & ...
    wt.mixin.BackgroundColorable & wt.mixin.PropertyViewable
    
    % A password entry field
    
    % Copyright 2020-2023 The MathWorks Inc.
    
    
    %% Public properties
    properties (AbortSet)
        
        % The current value shown
        Value (1,1) string
        
                
    end %properties
    
    
    %% Events
    events (HasCallbackProperty, NotifyAccess = protected)
        
        % Triggered on enter key pressed, has companion callback
        ValueChanged
        
        % Triggered on value changing during typing, has companion callback
        ValueChanging
        
    end %events
    
    
    
    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)
        
        % Grid
        Grid (1,1) matlab.ui.container.GridLayout
        
        % Password control
        PasswordControl (1,1) matlab.ui.control.HTML

        % Previous value for callback
        PrevValue (1,1) string
        
    end %properties
    
    
    %% Protected methods
    methods (Access = protected)
        
        function setup(obj)
            
            % Set default size
            obj.Position(3:4) = [100 25];
            
            % Construct Grid Layout to Manage Building Blocks
            obj.Grid = uigridlayout(obj);
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.RowSpacing = 2;
            obj.Grid.ColumnSpacing = 2;
            obj.Grid.Padding = 0;   
            
            % Define the HTML source
            html = ['<input type="password" id="pass" name="password" style="width:100%;height:100%" >',...
                '<script type="text/javascript">',...
                'function setup(htmlComponent) {',...
                '  htmlComponent.addEventListener("DataChanged", function(event) {',... %On uihtml Data prop changed
                '      if (document.getElementById("pass").value !== htmlComponent.Data) {',... %Does JS value not match uihtml Data?
                '        document.getElementById("pass").value = htmlComponent.Data;',... %Update the JS value to uihtml Data
                '      }',...
                '    });',...
                '  document.getElementById("pass").addEventListener("input", function() {',...%On input to html field
                '    htmlComponent.Data = document.getElementById("pass").value;',... %Copy JS value to uihtml data
                '    });',...
                '  document.addEventListener("keyup", function(event) {',... %Listen to key press
                '      if (event.keyCode === 13) {',... %If ENTER key
                '        htmlComponent.Data = document.getElementById("pass").value + ''\n'';',... %Add CR to data to trigger callback
                '      }',...
                '    });',...
                '}',...
                '</script>'];
            
            % Create a html password input
            obj.PasswordControl = uihtml(...
                'Parent',obj.Grid,...
                'HTMLSource',html,...
                'DataChangedFcn',@(h,e)obj.onPasswordChanged(e) );

            % Establish Background Color Listener
            obj.BackgroundColorableComponents = obj.Grid;

        end %function

 
        function update(obj)
            
            % Update the edit control text
            obj.PasswordControl.Data = obj.Value;
            
        end %function
        
        
        function propGroups = getPropertyGroups(obj)
            % Override the ComponentContainer GetPropertyGroups with newly
            % customiziable mixin. This can probably also be specific to each control.

            propGroups = getPropertyGroups@wt.mixin.PropertyViewable(obj);

        end % function

    end %methods
    
    
    
    %% Private methods
    methods (Access = private)
        
        function onPasswordChanged(obj,evt)
            % Triggered on interaction

            % Grab the data in string format
            newValue = string(evt.Data);
            oldValue = obj.Value;

            % Look at the states
            if endsWith(evt.PreviousData, newline)
            % This is needed to ignore double events

                % Clear the newline from the uihtml data
                newValue = erase(newValue, newline);
                obj.PasswordControl.Data = newValue;

            elseif endsWith(newValue, newline)
                % Enter key was pressed in the uihtml component

                % Clear the newline from the new value
                newValue = erase(newValue, newline);

                % Trigger event
                evtOut = wt.eventdata.PropertyChangedData('Value', newValue);
                notify(obj,"ValueChanged",evtOut);

            elseif newValue ~= oldValue

                % Store new result
                obj.Value = newValue;

                % Trigger event
                evtOut = wt.eventdata.PropertyChangedData('Value', ...
                    newValue, oldValue);
                notify(obj,"ValueChanging",evtOut);

            end

        end %function
    
    end %methods
    
end % classdef

