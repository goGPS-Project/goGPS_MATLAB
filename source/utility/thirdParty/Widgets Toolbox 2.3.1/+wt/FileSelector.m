classdef FileSelector < matlab.ui.componentcontainer.ComponentContainer & ...
        wt.mixin.BackgroundColorable & wt.mixin.ErrorHandling & ...
        wt.mixin.Enableable & wt.mixin.FontStyled & wt.mixin.Tooltipable &...
        wt.mixin.FieldColorable & wt.mixin.ButtonColorable & ...
        wt.mixin.PropertyViewable
    
    % File or folder selection control with browse button
    
    % Copyright 2020-2023 The MathWorks Inc.
    
    
    %% Public properties
    properties (AbortSet)
        
        % The current value shown. If a RootDirectory is used, this will be
        % a relative path to the root.
        Value (1,1) string

        % Filter applicable to the uigetfile/uiputfile dialog
        Filter string
        
    end %properties
    
    
    properties (Dependent, SetAccess = protected)
        
        % Absolute path to the file. If RootDirectory is used, this will
        % show the full path combining the root and the Value property.
        FullPath (1,1) string
        
        % Indicates the current value is a valid file path
        ValueIsValidPath (1,1) logical
        
    end %properties
    
    
    properties (AbortSet)
        
        % Selection type: (get)file, folder, putfile
        SelectionType (1,1) wt.enum.FileFolderState = wt.enum.FileFolderState.file
        
        % Optional root directory. If unspecified, Value uses an absolute
        % path (default). If specified, Value will show a relative path to
        % the root directory and Value must be beneath RootDirectory.
        RootDirectory (1,1) string
        
        % Indicates whether to show a dropdown of recent folders
        ShowHistory (1,1) matlab.lang.OnOffSwitchState = false
        
        % List of recently selected folders to display in dropdown
        History (:,1) string
        
    end %properties
    
    
    % These properties do not trigger the update method
    properties (AbortSet, UsedInUpdate = false)
        
        % Optional default directory to start in when Value does not exist
        % and user clicks the browse button.
        DefaultDirectory (1,1) string
        
    end %properties
    
    
    
    %% Events
    events (HasCallbackProperty, NotifyAccess = protected)
        
        % Triggered on value changed, has companion callback
        ValueChanged
        
    end %events
    
    
    
    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)
        
        % Button
        ButtonControl (1,1) matlab.ui.control.Button

        % Grid
        Grid (1,1) matlab.ui.container.GridLayout
        
        % Edit control or dropdown
        EditControl (1,1) matlab.ui.control.EditField
        DropdownControl (1,1) matlab.ui.control.DropDown
        
        % Image for invalid file
        WarnImage (1,1) matlab.ui.control.Image
        
    end %properties
    
    
    
    %% Protected methods
    methods (Access = protected)
        
        function setup(obj)
            
            % Adjust default size
            obj.Position(3:4) = [400 25];
            
            % Construct Grid Layout to Manage Building Blocks
            obj.Grid = uigridlayout(obj);
            obj.Grid.ColumnWidth = {'1x',25,25,25};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.RowSpacing = 2;
            obj.Grid.ColumnSpacing = 2;
            obj.Grid.Padding = 0;   
            
            % Create the standard edit control
            obj.EditControl = uieditfield(obj.Grid);
            obj.EditControl.ValueChangedFcn = @(h,e)obj.onTextChanged(e);
            obj.EditControl.Layout.Column = [1 3];
            obj.EditControl.Layout.Row = 1;
            
            % Create the optional dropdown control (unparented for now)
            obj.DropdownControl = uidropdown("Parent",[]);
            obj.DropdownControl.Editable = true;
            obj.DropdownControl.Value = "";
            obj.DropdownControl.ValueChangedFcn = @(h,e)obj.onTextChanged(e);
            
            % Create Button
            obj.ButtonControl = uibutton(obj.Grid);
            obj.ButtonControl.Text = "";
            obj.ButtonControl.ButtonPushedFcn = @(h,e)obj.onButtonPushed(e);
            obj.ButtonControl.Layout.Column = 4;
            obj.ButtonControl.Layout.Row = 1;
            obj.updateButtonIcon();

            % Create overlay
            obj.WarnImage = uiimage(obj.Grid);
            obj.WarnImage.ScaleMethod = "none";
            obj.WarnImage.Visible = "off";
            obj.WarnImage.ImageSource = "warning_16.png";
            obj.WarnImage.Layout.Column = 3;
            obj.WarnImage.Layout.Row = 1;
            
            % Update the internal component lists
            obj.FontStyledComponents = [obj.EditControl, obj.DropdownControl];
            obj.FieldColorableComponents = [obj.EditControl, obj.DropdownControl];
            obj.EnableableComponents = [obj.EditControl, obj.DropdownControl, obj.ButtonControl];
            obj.ButtonColorableComponents = obj.ButtonControl;
            obj.TooltipableComponents = [obj.EditControl, obj.DropdownControl, obj.ButtonControl];
            obj.BackgroundColorableComponents = obj.Grid;
            
        end %function
        
        
        function update(obj)
            
            % Is history being shown?
            if obj.ShowHistory
                % YES - Using dropdown control
                
                % Get the dropdown (history) items
                histItems = obj.History;
                
                % Current value must be in the list
                if ~isempty(obj.Value) && ~any(histItems == obj.Value)
                    histItems = vertcat(obj.Value,histItems);
                end
                
                % Update the items and value
                obj.DropdownControl.Items = histItems;
                obj.DropdownControl.Value = obj.Value;
                
            else
                % NO - Using edit control
                
                % Update the edit control text
                obj.EditControl.Value = obj.Value;
                
            end %if obj.ShowHistory

            % Show the warning icon?
            showWarn = strlength(obj.Value) && ~obj.ValueIsValidPath;
            obj.WarnImage.Visible = showWarn;
            
        end %function
        
        
        function propGroups = getPropertyGroups(obj)
            % Override the ComponentContainer GetPropertyGroups with newly
            % customiziable mixin. This can probably also be specific to each control.

            propGroups = getPropertyGroups@wt.mixin.PropertyViewable(obj);

        end %function

        
        function updateButtonIcon(obj)
            
            % Update the button icon
            if obj.SelectionType == "file"
                obj.ButtonControl.Icon = "folder_file_24.png";
            elseif obj.SelectionType == "putfile"
                obj.ButtonControl.Icon = "folder_file_24.png";
            else
                obj.ButtonControl.Icon = "folder_24.png";
            end
            
        end %function
        
        
        function updateControlType(obj)
            
            % Is history being shown? If so, update history and items
            if obj.ShowHistory
                % Using dropdown control
                
                if isempty(obj.DropdownControl.Parent)
                    obj.EditControl.Parent = [];
                    obj.DropdownControl.Parent = obj.Grid;
                    obj.DropdownControl.Layout.Column = [1 3];
                    obj.DropdownControl.Layout.Row = 1;
                    obj.WarnImage.Parent = [];
                    obj.WarnImage.Parent = obj.Grid;
                    obj.WarnImage.Layout.Column = 2;
                end
            else
                % Using edit control
                
                if isempty(obj.EditControl.Parent)
                    obj.DropdownControl.Parent = [];
                    obj.EditControl.Parent = obj.Grid;
                    obj.EditControl.Layout.Column = [1 3];
                    obj.EditControl.Layout.Row = 1;
                    obj.WarnImage.Parent = [];
                    obj.WarnImage.Parent = obj.Grid;
                    obj.WarnImage.Layout.Column = 3;
                end
            end %if obj.ShowHistory
            
        end %function
        
        
        function onButtonPushed(obj,~)
            % Triggered on button pushed
            
            % What folder should the prompt start at?
            if obj.ValueIsValidPath
                initialPath = obj.FullPath;
            elseif exist(obj.DefaultDirectory,"dir")
                initialPath = obj.DefaultDirectory;
            else
                initialPath = pwd;
            end

            % Get the filter
            filter = obj.Filter;
            if isempty(filter)
                filter = "";
            end

            % Get the figure (to restore focus)
            fig = ancestor(obj,"figure");
            
            % Prompt user for the path
            if obj.SelectionType == "file"
                [fileName,pathName] = uigetfile(filter,"Select a file",initialPath);
            elseif obj.SelectionType == "putfile"
                [fileName,pathName] = uiputfile(filter,"Specify an output file",initialPath);
            else
                pathName = uigetdir(initialPath, "Select a folder");
                fileName = "";
            end

            % Restore figure focus
            if isscalar(fig) && isvalid(fig)
                figure(fig)
            end
            
            % Proceed if user didn't cancel
            if ~isequal(pathName,0)
                obj.setValueFromFullPath( fullfile(pathName,fileName) );
            end %if ~isequal(pathName,0)
            
        end %function
        
        
        function onTextChanged(obj,evt)
            % Triggered on text interaction
            
            % Prepare event data
            evtOut = wt.eventdata.PropertyChangedData('Value', evt.Value, obj.Value);
            
            % Store new result
            obj.Value = evt.Value;
            
            % Trigger event
            notify(obj,"ValueChanged",evtOut);
            
        end %function
        
        
        function addToHistory(obj,value)
            % Add new item to history and clean up
            
            % Add this to the top of the history
            histFiles = obj.History;
            histFiles = vertcat(value,histFiles);
            
            % Clean up the history
            histFiles = wt.utility.cleanPath(histFiles);
            histFiles = unique(histFiles,"stable");
            
            % Filter to valid paths
            if obj.SelectionType == "file"
                fcn = @(x)gt(exist(fullfile(obj.RootDirectory, x),"file"), 0);
            elseif obj.SelectionType == "putfile"
                fcn = @(x)gt(exist(fullfile(obj.RootDirectory, x),"file"), 0);
            else
                fcn = @(x)eq(exist(fullfile(obj.RootDirectory, x),"dir"), 7);
            end
            isValidPath = arrayfun(fcn,histFiles);
            histFiles(~isValidPath) = [];
            
            % Limit history length
            histFiles(11:end) = [];
            
            % Store the result
            obj.History = histFiles;
            
        end %function
        
        
        function setValueFromFullPath(obj,fullPath)
            
            arguments
                obj (1,1)
                fullPath (1,1) string
            end %arguments
            
            fullPath = wt.utility.cleanPath(fullPath);
            
            if strlength(obj.RootDirectory) && strlength(fullPath)
                % Calculate the relative path for the value
                
                % Split the paths apart
                rootParts = strsplit(obj.RootDirectory,filesep);
                fullParts = strsplit(fullPath,filesep);
                
                % Find where the paths diverge
                idx = 1;
                smallestPath = min(numel(rootParts), numel(fullParts));
                while idx<=smallestPath && strcmpi(rootParts(idx),fullParts(idx))
                    idx = idx+1;
                end
                
                % Is the specified path outside of the root directory?
                numAbove = max(numel(rootParts) - idx + 1, 0);
                if numAbove>0
                    
                    msg = "Select a path beneath the root:" + ...
                        newline + newline + obj.RootDirectory;
                    obj.throwError(msg);
                    
                else
                    
                    % In case full path is above the RootPath, add ".." paths
                    parentPaths = string(repmat(['..' filesep],1,numAbove));
                    
                    % Form the relative path
                    relPath = filesep + fullfile(parentPaths, fullParts{idx:end});
                    
                    % What if paths are still the same?
                    if isempty(relPath)
                        relPath = "." + filesep;
                    end
                    
                    % Prepare event data
                    evtOut = wt.eventdata.PropertyChangedData('Value',relPath, obj.Value);
                    
                    % Store new result
                    obj.Value = relPath;
                    
                    % Trigger event
                    notify(obj,"ValueChanged",evtOut);
                    
                end
                
            else
                    
                    % Prepare event data
                    evtOut = wt.eventdata.PropertyChangedData('Value',fullPath, obj.Value);
                    
                    % Store new result
                    obj.Value = fullPath;
                    
                    % Trigger event
                    notify(obj,"ValueChanged",evtOut);
                
            end %if
            
        end %function
        
    end % methods
    
    
    %% Accessors
    methods
        
        function set.Value(obj,value)
            value = wt.utility.cleanPath(value);
            obj.Value = value;
            obj.addToHistory(value)
        end
        
        function set.RootDirectory(obj,value)
            value = wt.utility.cleanPath(value);
            obj.RootDirectory = value;
            obj.addToHistory(value)
        end
        
        function set.SelectionType(obj,value)
            obj.SelectionType = value;
            obj.updateButtonIcon()
        end
        
        function set.ShowHistory(obj,value)
            obj.ShowHistory = value;
            obj.updateControlType()
        end
        
        function value = get.FullPath(obj)
            value = fullfile(obj.RootDirectory, obj.Value);
        end
        
        function value = get.ValueIsValidPath(obj)
            filePath = fullfile(obj.RootDirectory, obj.Value);
            value = ( obj.SelectionType == "file" && isfile(filePath) ) || ...
                ( obj.SelectionType == "putfile" && isfolder(fileparts(filePath)) ) || ...
                ( obj.SelectionType == "folder" && isfolder(filePath) );
        end
        
    end % methods
    
    
end % classdef