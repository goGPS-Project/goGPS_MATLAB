classdef BaseApp < matlab.apps.AppBase & matlab.mixin.SetGetExactNames & ...
        wt.mixin.ErrorHandling
    % Base class for Widgets Toolbox apps
    
    % Copyright 2020-2023 The MathWorks, Inc.
    
    
    %% Properties
    properties (AbortSet)
        
        % Name of the app
        Name (1,1) string = "My App"
        
    end %properties
    
    
    properties (SetAccess = protected)
        
        % Model class for App preferences
        %(may subclass wt.model.Preferences to add more prefs)
        Preferences (1,1) wt.model.Preferences
        
        % Name of group to store preferences (defaults to class name)
        PreferenceGroup (1,1) string
        
    end %properties
    
    
    properties (AbortSet, Dependent)
        
        % Position of the app window
        Position
        
        % Visibility of the app window
        Visible
        
        % State of the app window
        WindowState
        
    end %properties
    
    
    
    %% Internal properties
    properties (Hidden, Transient, NonCopyable, SetAccess = immutable)
        
        % Figure window of the app
        Figure matlab.ui.Figure
        
        % Primary grid to place contents
        Grid matlab.ui.container.GridLayout
        
    end %properties
    
    
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)
        
        % Last used folder (for file operations)
        LastFolder (1,1) string = pwd
        
        % Is setup complete?
        SetupComplete (1,1) logical = false;
        
    end %properties
    
    
    
    %% Abstract methods (subclass must implement these)
    methods (Abstract, Access = protected)
        
        setup(app)
        
        update(app)
        
    end %methods
    
    
    
    %% Debugging Methods
    methods
        
        function forceUpdate(app)
            % Forces update to run (For debugging only!)
            
            app.update();
            
        end %function
        
    end %methods
    
    
    
    %% Constructor / destructor
    methods (Access = public)
        
        function app = BaseApp(varargin)
            % Constructor
            
            % Create the figure and hide until components are created
            app.Figure = uifigure( ...
                'AutoResizeChildren','off',...
                'Units','pixels', ...
                'DeleteFcn',@(h,e)delete(app), ...
                'CloseRequestFcn',@(h,e)close(app), ...
                'Visible','off');
            
            % Create MainGridLayout
            app.Grid = uigridlayout(app.Figure,[1 1]);
            app.Grid.Padding = [0 0 0 0];
            
            % Check for preference input and assign it first, in case
            % Preferences was subclassed
            % [splitArgs,varargin] = app.splitArgs('Preferences', varargin{:});
            % if ~isempty(splitArgs)
            %     app.Preferences = splitArgs{2};
            % end
            
            % Retrieve preferences
            app.loadPreferences();
            
            % Load last figure position
            % Note app.Figure.Position is inner position, where
            % app.Position is outer position. Inner position is the same
            % regardless of any menubar or toolbar, while outer position
            % changes if either is added. Use inner position for this
            % purpose so it does not depend on if/when any menubar or
            % toolbar is added to the figure.
            app.Figure.Position = app.getPreference('Position',[100 100 1000 700]);
            
            % Set up components
            app.setup_internal();
            app.setup();
            
            % Set any P-V pairs
            if ~isempty(varargin)
                set(app, varargin{:});
            end
            
            % Register the app with App Designer
            registerApp(app, app.Figure)
            
            % Ensure it's on screen
            app.moveOnScreen();
            
            % Mark the setup complete
            app.SetupComplete = true;
            
            % Update the app
            app.update();
            
            % Update the title
            app.updateTitle();
            
            % Force drawing to finish
            drawnow
            
            % Now, make it visible
            app.Figure.Visible = 'on';
            
        end %function
        
        
        function delete(app)
            % Destructor
            
            % Store last position in preferences
            if isscalar(app.Figure) && isvalid(app.Figure)
                app.setPreference('Position',app.Figure.Position)
            end
            
            % Save preferences
            app.savePreferences();
            
            % Now, delete the figure
            delete(app.Figure)
            
        end %function
        
    end %methods
    
    
    
    %% Public Methods
    methods
        
        function close(app)
            % Triggered on figure closed
            
            app.delete();
            
        end %function
        
        
        function selection = promptYesNoCancel(app, message, title, default, icon)
            % Prompt the user with a yes/no/cancel selection
            
            % Define arguments
            arguments
                app (1,1) wt.apps.BaseApp
                message (1,1) string = "Are you sure?"
                title (1,1) string = ""
                default (1,1) string = "Cancel"
                icon (1,1) string = "question"
            end
            
            % Launch the prompt
            selection = uiconfirm(app.Figure, message, title,...
                "Options",["Yes","No","Cancel"],...
                "DefaultOption",default,...
                "CancelOption","Cancel",...
                "Icon",icon);
            
        end %function
        
        
        function filePath = promptToSaveAs(app, filePath, filter, title)
            % Prompt the user to save a file
            
            % Define arguments
            arguments
                app (1,1) wt.apps.BaseApp
                filePath (1,1) string = pwd
                filter = ["*.mat","MATLAB MAT File"];
                title (1,1) string = "Save as"
            end
            
            % Prompt for the file
            [fileName,pathName] = uiputfile(filter, title, filePath);
            
            % Did the user cancel?
            if isequal(fileName,0)
                filePath = string.empty(0);
            else
                filePath = fullfile(pathName,fileName);
                app.LastFolder = pathName;
            end %if isequal(fileName,0)
            
        end %function
        
        
        function filePath = promptToLoad(app, filter, title)
            % Prompt the user to load a file
            
            % Define arguments
            arguments
                app (1,1) wt.apps.BaseApp
                filter = ["*.mat","MATLAB MAT File"];
                title (1,1) string = "Open"
            end
            
            % Prompt for the file
            [fileName,pathName] = uigetfile(filter, title, app.LastFolder);
            
            % Did the user cancel?
            if isequal(fileName,0)
                filePath = string.empty(0);
            else
                filePath = fullfile(pathName,fileName);
                app.LastFolder = pathName;
            end %if isequal(fileName,0)
            
        end %function
        
    end %methods
    
    
    
    %% Sealed Public methods
    methods (Sealed)
        
        function value = getPreference(app,propName,defaultValue)
            % Get an app preference from the Preferences object
            
            if isprop(app.Preferences, propName)
                value = app.Preferences.(propName);
            elseif nargin>2
                value = defaultValue;
            else
                value = [];
            end
            
        end %function
        
        
        function setPreference(app,propName,value)
            % Set an app preference in the Preferences object
            
            if ~isprop(app.Preferences,propName)
                addprop(app.Preferences,propName);
            end
            app.Preferences.(propName) = value;
            
        end %function
        
        
        function moveOnScreen(app)
            % Ensure the figure is placed on screen
            
            if strcmp(app.Figure.Units,'pixels')
                
                % Get the corners of each screen
                g = groot;
                screenPos = g.MonitorPositions;
                screenCornerA = screenPos(:,1:2);
                screenCornerB = screenPos(:,1:2) + screenPos(:,3:4) - 1;
                
                % Buffer for title bar
                titleBarHeight = 30;
                
                % Get the corners of the figure (bottom left and top right)
                figPos = app.Figure.OuterPosition;
                figCornerA = figPos(1:2);
                figCornerB = figPos(1:2) + figPos(:,3:4) - 1;
                
                % Are the corners on any screen?
                aIsOnScreen = all( figCornerA >= screenCornerA & ...
                    figCornerA <= screenCornerB, 2 );
                bIsOnScreen = all( figCornerB >= screenCornerA & ...
                    figCornerB <= screenCornerB, 2);
                
                % Are corners on a screen?
                
                % Are both corners fully on any screen?
                if any(aIsOnScreen) && any(bIsOnScreen)
                    % Yes - do nothing
                    
                elseif any(bIsOnScreen)
                    % No - only upper right corner is on a screen
                    
                    % Calculate the adjustment needed, and make it
                    figAdjust = max(figCornerA, screenCornerA(bIsOnScreen,:)) ...
                        - figCornerA;
                    figPos(1:2) = figPos(1:2) + figAdjust;
                    
                    % Ensure the upper right corner still fits
                    figPos(3:4) = min(figPos(3:4), ...
                        screenCornerB(bIsOnScreen,:) - figPos(1:2) - [0 titleBarHeight] + 1);
                    
                    % Move the figure
                    app.Figure.Position = figPos;
                    
                elseif any(aIsOnScreen)
                    % No - only lower left corner is on a screen
                    
                    % Calculate the adjustment needed, and make it
                    figAdjust = min(figCornerB, screenCornerB(aIsOnScreen,:)) ...
                        - figCornerB;
                    figPos(1:2) = max( screenCornerA(aIsOnScreen,:),...
                        figPos(1:2) + figAdjust );
                    
                    % Ensure the upper right corner still fits
                    figPos(3:4) = min(figPos(3:4), ...
                        screenCornerB(aIsOnScreen,:) - figPos(1:2) - [0 titleBarHeight] + 1);
                    
                    % Move the figure
                    app.Figure.Position = figPos;
                    
                else
                    % No - Not on any screen
                    
                    % This is slower, but uncommon anyway
                    movegui(app.Figure,'onscreen');
                    
                end %if any( all(aIsOnScreen,2) & all(bIsOnScreen,2) )
                
            else
                
                % This is slower, but uncommon anyway
                movegui(app.Figure,'onscreen');
                
            end %if strcmp(app.Figure.Units,'pixels')
            
        end %function
        
    end %methods
    
    
    
    %% Protected Methods
    methods (Access = protected)
        
        function setup_internal(~)
            % Preform internal pre-setup necessary
            
            % This is used for session managed apps
            
        end %function
        
        
        function loadPreferences(app)
            % Load stored preferences
            
            app.Preferences.load(app.PreferenceGroup);
            
        end %function
        
        
        function savePreferences(app)
            % Save preferences
            
            app.Preferences.save(app.PreferenceGroup);
            
        end %function
        
        
        function updateTitle(app)
            % Update the figure title
            
            app.Figure.Name = app.Name;
            
        end %function
        
    end %methods
    
    
    
    %% Accessors
    methods
        
        function set.Name(app,value)
            app.Name = value;
            app.updateTitle();
        end
        
        
        function value = get.PreferenceGroup(app)
            value = app.PreferenceGroup;
            if isempty(value)
                value = class(app);
            end
            value = matlab.lang.makeValidName(value);
        end
        
        
        function value = get.Position(app)
            value = app.Figure.Position;
        end
        
        function set.Position(app,value)
            app.Figure.Position = value;
        end
        
        
        function value = get.Visible(app)
            value = app.Figure.Visible;
        end
        
        function set.Visible(app,value)
            app.Figure.Visible = value;
        end
        
        
        function value = get.WindowState(app)
            value = app.Figure.WindowState;
        end
        
        function set.WindowState(app,value)
            app.Figure.WindowState = value;
        end
        
    end %methods
    
    
end % classdef

