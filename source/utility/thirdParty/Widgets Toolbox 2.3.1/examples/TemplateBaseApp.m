classdef TemplateBaseApp < wt.apps.BaseApp
    % Implements a template for a BaseApp
    
    % Copyright 2020-2023 The MathWorks Inc.


    %% Internal Components
    %   Create properties here for each control, layout, or view component
    %   that will be placed directly into the main app window.
    properties ( Transient, NonCopyable, Hidden, SetAccess = protected )

        % Grid Layouts
        % (BaseApp already brings "Grid", the main Grid layout in the window)
        Tab1Grid matlab.ui.container.GridLayout
        Tab2Grid matlab.ui.container.GridLayout
        Panel1Grid matlab.ui.container.GridLayout
        Panel2Grid matlab.ui.container.GridLayout

        % Panels
        Panel1 matlab.ui.container.Panel
        Panel2 matlab.ui.container.Panel

        % Tabs
        TabGroup matlab.ui.container.TabGroup
        Tab1 matlab.ui.container.Tab
        Tab2 matlab.ui.container.Tab

        % View Components
        %View1 packageName.ClassName
        %View2 packageName.ClassName

        % Temporary label components
        Panel1Text matlab.ui.control.Label
        Panel2Text matlab.ui.control.Label
        Tab1Text matlab.ui.control.Label
        Tab2Text matlab.ui.control.Label

    end %properties



    %% Setup and Configuration of the App
    methods  ( Access = protected )

        function setup(app)
            % Runs once on instantiation of the app

            % Set the name
            app.Name = "My App";

            % Configure the main grid
            app.Grid.ColumnWidth = {300,'1x'};
            app.Grid.RowHeight = {'1x',150};
            app.Grid.Padding = 5;

            % Create a panel
            app.Panel1 = uipanel(app.Grid);
            app.Panel1.Title = "Panel 1";
            app.Panel1.Layout.Row = [1 2];
            app.Panel1.Layout.Column = 1;
            app.Panel1.BackgroundColor = "magenta";

            % Create a panel
            app.Panel2 = uipanel(app.Grid);
            app.Panel2.Title = "Panel 2";
            app.Panel2.Layout.Row = 2;
            app.Panel2.Layout.Column = 2;
            app.Panel2.BackgroundColor = "white";

            % Create a tab group
            app.TabGroup = uitabgroup(app.Grid);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 2;
            app.TabGroup.SelectionChangedFcn = @(h,e)onTabChanged(app,e);

            % Create Tabs, each with a grid
            app.Tab1 = uitab(app.TabGroup);
            app.Tab1.Title = 'Tab 1';

            app.Tab2 = uitab(app.TabGroup);
            app.Tab2.Title = 'Tab 2';

            % Create grid layouts to position content inside each container
            app.Panel1Grid = uigridlayout(app.Panel1, [1,1],"Padding", 0);
            app.Panel2Grid = uigridlayout(app.Panel2, [1,1],"Padding", 0);
            app.Tab1Grid = uigridlayout(app.Tab1, [1,1],"Padding", 0);
            app.Tab2Grid = uigridlayout(app.Tab2, [1,1],"Padding", 0);

            % Place some temporary content in each container
            app.Panel1Text = uilabel(app.Panel1Grid);
            app.Panel1Text.Text = "Panel 1 Contents";
            app.Panel1Text.Layout.Row = 1;
            app.Panel1Text.Layout.Column = 1;
            app.Panel1Text.BackgroundColor = "red";

            app.Panel2Text = uilabel(app.Panel2Grid);
            app.Panel2Text.Text = "Panel 2 Contents";
            app.Panel2Text.Layout.Row = 1;
            app.Panel2Text.Layout.Column = 1;
            app.Panel2Text.BackgroundColor = "green";

            app.Tab1Text = uilabel(app.Tab1Grid);
            app.Tab1Text.Text = "Tab 1 Contents";
            app.Tab1Text.Layout.Row = 1;
            app.Tab1Text.Layout.Column = 1;
            app.Tab1Text.BackgroundColor = "cyan";

            app.Tab2Text = uilabel(app.Tab2Grid);
            app.Tab2Text.Text = "Tab 2 Contents";
            app.Tab2Text.Layout.Row = 1;
            app.Tab2Text.Layout.Column = 1;
            app.Tab2Text.BackgroundColor = "blue";

            % Additional examples:
            % (add other views, layouts, and components here as needed)
            %app.View1 = packageName.ClassName( app.Tab1Grid );
            %app.View2 = packageName.ClassName( app.Tab2Grid );

        end %function

    end %methods

    %% Update
    methods  ( Access = protected )

        function update(app)
            % Update the display of the app
            % For the main app, app.update() must be called explicitly
            % during callbacks or other changes that require contents to
            % refresh.

        end %function

    end %methods


    %% Callbacks
    methods

        function onTabChanged(app,evt)
            % Triggered on changing tab

            % (this method is optional if using tabs)

            newTab = evt.NewValue;
            disp("Selected Tab: " + newTab.Title);

        end %function

    end %methods



    %% Private Methods
    methods ( Access = private )

        % (optionally add any internal methods here)

    end %methods

end %classdef