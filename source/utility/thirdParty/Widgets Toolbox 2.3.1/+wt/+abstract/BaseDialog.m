classdef BaseDialog < wt.abstract.BaseWidget & ...
        wt.mixin.FontStyled & ...
        wt.mixin.ButtonColorable & ...
        wt.mixin.FieldColorable
    % Base class for a dialog panel

    % Copyright 2022-2023 The MathWorks Inc.
    

    %% Public Properties
    properties (AbortSet, Dependent, Access = public)

        % Dialog Title
        Title

        % Dialog Size
        Size (1,2) double {mustBePositive}

    end %properties


    % Accessors
    methods

        function value = get.Title(obj)
            value = obj.TitleBar.Text;
        end

        function set.Title(obj, value)
            obj.TitleBar.Text = value;
        end

        function value = get.Size(obj)
            value = obj.Position(3:4);
        end

        function set.Size(obj, value)
            obj.Position(3:4) = value;
        end

    end %methods
   

    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = private)

        % Dialog's grid managing the outer container
        DialogGrid matlab.ui.container.GridLayout
       
        % Title text
        TitleBar matlab.ui.control.Label

        % Close button
        CloseButton matlab.ui.control.Button

        % Listeners to reference/parent objects to trigger dialog delete
        LifecycleListeners (1,:) event.listener

    end %properties



    %% Protected methods
    methods (Access = protected)
        
        function setup(obj)
            % Configure the widget

            % Configure component
            %obj.Position = [10 10 300 300];
            obj.BackgroundColor = [0 0 0];

            % Call superclass method
            obj.setup@wt.abstract.BaseWidget();
            
            % Outer Grid to manage building blocks
            obj.DialogGrid = uigridlayout(obj);
            obj.DialogGrid.RowHeight = {'fit','1x'};
            obj.DialogGrid.ColumnWidth = {'1x',25};
            obj.DialogGrid.RowSpacing = 1;
            obj.DialogGrid.ColumnSpacing = 1;
            obj.DialogGrid.Padding = 1;
            obj.DialogGrid.BackgroundColor = [.5 .5 .5];

            % Title Bar
            obj.TitleBar = uilabel(obj.DialogGrid);
            obj.TitleBar.Text = "  Configuration";
            obj.TitleBar.BackgroundColor = [1 1 1];
            obj.TitleBar.Layout.Row = 1;
            obj.TitleBar.Layout.Column = 1;

            % Close Button
            obj.CloseButton = uibutton(obj.DialogGrid);
            obj.CloseButton.Text = "X";
            obj.CloseButton.ButtonPushedFcn = @(src,evt)obj.onClosePressed();
            obj.CloseButton.Layout.Row = 1;
            obj.CloseButton.Layout.Column = 2;
            
            % Inner Grid to manage building blocks
            obj.Grid.Parent = obj.DialogGrid;
            obj.Grid.Padding = 10;
            obj.Grid.Scrollable = true;
            obj.Grid.Layout.Row = 2;
            obj.Grid.Layout.Column = [1 2];
            
        end %function

        
        function updateBackgroundColorableComponents(obj)
            % Update components that are affected by BackgroundColor
            
            % Update dialog components
            obj.Grid.BackgroundColor = obj.BackgroundColor;
            obj.TitleBar.BackgroundColor = obj.BackgroundColor;

            % Update other components
            hasProp = isprop(obj.BackgroundColorableComponents,'BackgroundColor');
            set(obj.BackgroundColorableComponents(hasProp),...
                "BackgroundColor",obj.BackgroundColor);
            
        end %function

        
        function updateButtonColorableComponents(obj)
            % Update components that are affected by ButtonColor
            
            % Update dialog components
            obj.CloseButton.BackgroundColor = obj.ButtonColor;

            % Update other components
            obj.updateButtonColorableComponents@wt.mixin.ButtonColorable()

        end %function


        function updateFontStyledComponents(obj,varargin)

            % Handle the dialog components
            if nargin == 3

                % Update title bar and close (X) button with specific font
                % properties only
                props = ["FontColor","FontName","FontSize"];
                if matches(varargin{1}, props)

                    set(obj.TitleBar,varargin{1},varargin{2})
                    set(obj.CloseButton,varargin{1},varargin{2})

                end %if

            end %if

            % Call superclass to update the rest
            obj.updateFontStyledComponents@wt.mixin.FontStyled(varargin{:});

        end %function


        function positionWithin(obj, fig, referenceWidget)
            % Positions within the figure relative to the referenceWidget

            % Must be within a figure
            if isempty(obj.Parent) || ...
                    ~isvalid(obj.Parent) || ...
                    ~isa(obj.Parent,"matlab.ui.Figure") || ...
                    obj.Units ~= "pixels"
                return
            end

            % What size is the dialog?
            dlgSize = obj.Position(3:4);

            % Reference widget position
            refPos = getpixelposition(referenceWidget, true);
            refSize = refPos(3:4);
            refCornerA = refPos(1:2);
            %refCornerB = refPos(1:2) + refPos(:,3:4) - 1;

            % Does it fit entirely within the reference component?
            if all(refSize >= dlgSize)
                % Yes - center it over the component

                % Calculate lower-left corner
                dlgPos = floor((refSize - dlgSize) / 2) + refCornerA;

            else
                % NO - position within the figure

                % Get the corners of the figure (bottom left and top right)
                figPos = getpixelposition(fig);
                figSize = figPos(3:4);

                % Start with dialog position in lower-left of widget
                dlgPos = refCornerA;
                dlgCornerB = dlgPos + dlgSize;

                % Move left and down as needed to fit in figure
                adj = figSize - dlgCornerB;
                adj(adj>0) = 0;
                dlgPos = max(dlgPos + adj, [1 1]);
                dlgCornerB = dlgPos + dlgSize;

                % If it doesn't fit in the figure, shrink it
                adj = figSize - dlgCornerB;
                adj(adj>0) = 0;
                dlgSize = dlgSize + adj;

            end %if

            % Set final position
            obj.Position = [dlgPos dlgSize];
            
        end %function


        function labels = addRowLabels(obj, names, parent, column, startRow)
            % Add a group of standard row labels to the grid (or specified
            % grid)

            arguments
                obj %#ok<INUSA> 
                names (:,1) string
                parent = obj.Grid
                column = 1
                startRow = 1
            end

            numRows = numel(names);
            labels = gobjects(1,numRows);
            hasText = false(1,numRows);
            for idx = 1:numel(names)
                thisName = names(idx);
                hasText(idx) = strlength(thisName) > 0;
                if hasText(idx)
                    h = uilabel(parent);
                    h.HorizontalAlignment = "right";
                    h.Text = thisName;
                    h.Layout.Column = column;
                    h.Layout.Row = idx + startRow - 1;
                    labels(idx) = h;
                end
            end

            % Remove the empty spaces
            labels(~hasText) = [];
            
        end %function


        function onClosePressed(obj)
            % Delete the dialog on close button pressed

            delete(obj)

        end %function


        function attachLifecycleListeners(obj, listenObj)
            % Delete the dialog upon destruction of the specified objects

            arguments
                obj (1,1) wt.abstract.BaseDialog
                listenObj handle   
            end

            % Create listeners
            % The dialog will be deleted if the listenObj is deleted
            newListeners = listener(listenObj, "ObjectBeingDestroyed",...
                @(src,evt)delete(obj));

            % Add to any existing listeners
            obj.LifecycleListeners = horzcat(obj.LifecycleListeners, newListeners);

        end %function


        function copyStylesFromReference(obj, reference)
            % Copy styles from reference to this dialog where possible

            % Style fields to check
            styleFields = [
                "BackgroundColor"
                "FontName"
                "FontSize"
                "FontWeight"
                "FontAngle"
                "FontColor"
                "ButtonColor"
                "FieldColor"
                ];

            % Loop on each style field
            for thisField = styleFields'

                % Does the reference component have this style fielid?
                if isprop(reference, thisField)

                    % Attempt to apply it to this dialog
                    try
                        thisValue = reference.(thisField);
                        obj.(thisField) = thisValue;
                    catch ME
                        warning("wt:BaseDialog:CopyStylesError",...
                            "Unable to apply style '%s' to dialog '%s'.\n%s",...
                            thisField, class(reference), ME.message);
                    end %try

                end %if

            end %for


        end %function
        
    end %methods

end %classdef