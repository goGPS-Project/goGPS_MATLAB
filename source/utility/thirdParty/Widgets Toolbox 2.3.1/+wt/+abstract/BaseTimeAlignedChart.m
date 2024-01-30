classdef BaseTimeAlignedChart < matlab.graphics.chartcontainer.ChartContainer & ...
        wt.mixin.FontStyled & ...
        wt.mixin.ErrorHandling
    % Base class for a chart with time-aligned axes

    % Copyright 2022-2023 The MathWorks Inc.


    %% Public Properties
    properties (AbortSet, Access = public)

        % How many axes to display
        NumAxes (1,1) double {mustBePositive,mustBeInteger} = 1

        % Which axes are selected?
        SelectedAxes (1,1) double {mustBePositive,mustBeInteger} = 1

        % Show grid on each axes?
        ShowGrid (1,1) logical = true;

        % Show legend on each axes?
        ShowLegend (1,1) logical = false;

        % Background color of axes components
        AxesColor (1,3) double ...
            {mustBeInRange(AxesColor,0,1)} = [1 1 1]

        % Selected color of axes components
        AxesSelectedColor (1,3) double ...
            {mustBeInRange(AxesSelectedColor,0,1)} = [.9 .9 1]

        % Grid color of axes components
        AxesGridColor (1,3) double ...
            {mustBeInRange(AxesGridColor,0,1)} = [.15 .15 .15]

    end %properties


    properties (AbortSet, Dependent, UsedInUpdate = false)

        % Axes Y-Limits
        YLim (:,1) cell

        % Axes Y-Limit Mode
        YLimMode (:,1) string

        % X-Axis Labels
        XLabel (:,1) string

        % Y-Axis Labels
        YLabel (:,1) string

        % Axes Titles
        Title (:,1) string

    end %properties


    % Accessors
    methods

        function set.NumAxes(obj,value)
            obj.NumAxes = value;
            obj.AxesLayoutDirty = true; %#ok<MCSUP>
        end

        function set.ShowLegend(obj,value)
            obj.ShowLegend = value;
            obj.AxesLayoutDirty = true; %#ok<MCSUP>
        end

        function set.AxesColor(obj,value)
            obj.AxesColor = value;
            obj.updateAxesColors();
        end

        function set.AxesGridColor(obj,value)
            obj.AxesGridColor = value;
            obj.updateAxesColors();
        end

        function value = get.YLim(obj)
            if isempty(obj.Axes)
                value = string.empty(0,1);
            else
                value = {obj.Axes.YLim}';
            end
        end

        function set.YLim(obj, value)
            numSet = min(numel(obj.Axes), numel(value));
            for idx = 1:numSet
                if ~isequal(obj.Axes(idx).YLim, value{idx,:})
                    obj.Axes(idx).YLim = value{idx,:};
                end
            end
        end

        function value = get.YLimMode(obj)
            if isempty(obj.Axes)
                value = string.empty(0,1);
            else
                value = string({obj.Axes.YLimMode}');
            end
        end

        function set.YLimMode(obj, value)
            numSet = min(numel(obj.Axes), numel(value));
            for idx = 1:numSet
                obj.Axes(idx).YLimMode = value(idx);
            end
        end

        function value = get.XLabel(obj)
            if isempty(obj.Axes)
                value = string.empty(0,1);
            else
                allXLabel = [obj.Axes.XLabel];
                value = string({allXLabel.String}');
            end
        end

        function set.XLabel(obj, value)
            numSet = min(numel(obj.Axes), numel(value));
            for idx = 1:numSet
                obj.Axes(idx).XLabel.String = value(idx);
            end
        end

        function value = get.YLabel(obj)
            if isempty(obj.Axes)
                value = string.empty(0,1);
            else
                allYLabel = [obj.Axes.YLabel];
                value = string({allYLabel.String}');
            end
        end

        function set.YLabel(obj, value)
            numSet = min(numel(obj.Axes), numel(value));
            for idx = 1:numSet
                obj.Axes(idx).YLabel.String = value(idx);
            end
        end

        function value = get.Title(obj)
            if isempty(obj.Axes)
                value = string.empty(0,1);
            else
                allTitle = [obj.Axes.Title];
                value = string({allTitle.String}');
            end
        end

        function set.Title(obj, value)
            numSet = min(numel(obj.Axes), numel(value));
            for idx = 1:numSet
                obj.Axes(idx).Title.String = value(idx);
            end
        end

    end %methods



    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls that BackgroundColor should apply to
        BackgroundColorableComponents (:,1) matlab.graphics.Graphics

        % TiledLayout for axes
        TiledLayout matlab.graphics.layout.TiledChartLayout

        % Axes to display the signal
        Axes (1,:) matlab.graphics.axis.Axes

        % Legend of each axes
        Legend (1,:) matlab.graphics.illustration.Legend

    end %properties


    properties (Transient, NonCopyable, UsedInUpdate, Access = private)

        % Used internally to request update call
        RequestUpdate_I (1,1) logical = false

    end %properties


    properties (Transient, NonCopyable, Access = private)

        % State of axes layout
        AxesLayoutDirty (1,1) logical = true

    end %properties


    % Hidden
    properties (Hidden)

        % Enable clicking axes to select?
        EnableSelection (1,1) logical = false

    end %properties


    %% Debugging Methods
    methods

        function forceUpdateChart(obj)
            % Forces update to run (For debugging only!)

            obj.update();

        end %function

    end %methods


    %% Protected methods
    methods (Access = protected)

        function setup(obj)
            % Create the underlying components

            % Configure Layout
            obj.TiledLayout = getLayout(obj);
            obj.TiledLayout.Padding = "compact";
            obj.TiledLayout.TileSpacing = "compact";

        end %function


        function update(obj)
            % Update the underlying components

            % Do we need to recreate the axes?
            if obj.AxesLayoutDirty
                obj.recreateContent();
            end

            % If an axes is selected, update colors to highlight it
            if obj.EnableSelection
                obj.updateAxesColors()
            end

            % Toggle legends and grids
            set(obj.Axes,"XGrid",obj.ShowGrid)
            set(obj.Axes,"YGrid",obj.ShowGrid)

        end %function


        function requestUpdate(obj)
            % Request update method to run

            % Trigger set of a UsedInUpdate property to request update
            % during next drawnow. (for optimal efficiency)
            obj.RequestUpdate_I = ~obj.RequestUpdate_I;

        end %function


        function recreateContent(obj)
            % Create / recreate all the content

            % Delete existing axes
            % This is required to change the tiledlayout size
            delete(obj.Axes);
            delete(obj.TiledLayout.Children)
            obj.Axes(:) = [];

            % Update the TiledLayout
            obj.TiledLayout.GridSize = [obj.NumAxes 1];

            % Create the axes
            ax = gobjects(1,obj.NumAxes);
            lgnd = matlab.graphics.illustration.Legend.empty(1,0);
            for idx = 1:obj.NumAxes

                % Create the axes
                ax(idx) = nexttile(obj.TiledLayout);

                % Configure axes
                ax(idx).NextPlot = "add";
                ax(idx).XAxis = matlab.graphics.axis.decorator.DurationRuler();

                % Keep X ticks only on the last axes
                if idx < obj.NumAxes
                    ax(idx).XTickLabel = {};
                end

                % Configure Interpreters
                ax(idx).Title.Interpreter = "none";
                ax(idx).XLabel.Interpreter = "none";
                ax(idx).YLabel.Interpreter = "none";
                ax(idx).YAxis.TickLabelInterpreter = "none";

                % Start legend
                if obj.ShowLegend
                    lgnd(idx) = legend(ax(idx),'Location',"northwest");
                    lgnd(idx).Interpreter = "none";
                end

            end %for

            % Set modes
            set(ax, "ZLimMode", "manual")

            % Configure callback for clicking on axes
            set(ax,"ButtonDownFcn",@(src,evt)onAxesButtonDown(obj,evt));

            % Link all the axes in X
            if numel(ax) > 1  && all(isvalid(ax))
                try
                    linkaxes(ax, 'x');
                catch err
                    warning("wtslrt:BaseTimeAlignedChart:LinkAxesError",...
                        "Unable to link axes: %s", err.message);
                    %RJ - Occasionally gets here, unsure why.
                    keyboard
                end
            end

            % Store new objects
            obj.Axes = ax;
            obj.Legend = lgnd;

            % Update styles
            obj.updateFontStyledComponents();
            obj.updateAxesColors();

            % Mark axes layout clean
            obj.AxesLayoutDirty = false;

        end %function

        function updateBackgroundColorableComponents(obj)
            % Update components that are affected by BackgroundColor

            obj.Grid.BackgroundColor = obj.BackgroundColor;
            hasProp = isprop(obj.BackgroundColorableComponents,'BackgroundColor');
            set(obj.BackgroundColorableComponents(hasProp),...
                "BackgroundColor",obj.BackgroundColor);

        end %function


        function updateFontStyledComponents(obj,prop,value)

            % All components together
            axesAndLegend = horzcat(obj.Axes, obj.Legend);
            allComponents = horzcat(axesAndLegend, obj.Axes.Title);

            if nargin < 2
                % Update all

                set(allComponents, "FontName", obj.FontName)
                set(allComponents, "FontSize", obj.FontSize)
                set(axesAndLegend, "FontWeight", obj.FontWeight)
                set(allComponents, "FontAngle", obj.FontAngle)

                set(obj.Axes, "XColor", obj.FontColor)
                set(obj.Axes, "YColor", obj.FontColor)
                set(obj.Axes, "ZColor", obj.FontColor)
                set(obj.Legend,"TextColor", obj.FontColor)
                set([obj.Axes.Title],"Color", obj.FontColor)

            elseif prop == "FontColor"

                set(obj.Axes, "XColor", value)
                set(obj.Axes, "YColor", value)
                set(obj.Axes, "ZColor", obj.FontColor)
                set(obj.Legend,"TextColor", value)
                set([obj.Axes.Title],"Color", obj.FontColor)

            elseif prop == "FontWeight"

                set(axesAndLegend, "FontWeight", obj.FontWeight)

            else
                % Update other specific property

                wt.utility.fastSet(allComponents,prop,value);

            end %if

        end %function


        function updateAxesColors(obj)
            % Update axes colors

            % All components together
            allComponents = horzcat(obj.Axes, obj.Legend);

            % Set axes and legend background
            set(allComponents,"Color", obj.AxesColor)

            % Update selected axes color
            if obj.EnableSelection && obj.SelectedAxes <= numel(obj.Axes)
                set(obj.Axes(obj.SelectedAxes),"Color", obj.AxesSelectedColor)
            end

            % Set grid color
            set(obj.Axes, "GridColor", obj.AxesGridColor)

        end %function


        function onAxesButtonDown(obj, evt)
            % Triggered on button pushed on axes

            if obj.EnableSelection

                % What was clicked?
                clickedAxes = evt.Source;

                % Find the selected index
                axIdx = find(clickedAxes == obj.Axes, 1);

                % Set the selection
                % This will trigger update
                obj.SelectedAxes = axIdx;

            end %if

        end %function

    end %methods

end %classdef