classdef (Abstract) BaseWidget < ...
        matlab.ui.componentcontainer.ComponentContainer & ...
        wt.mixin.BackgroundColorable & ...
        wt.mixin.PropertyViewable & ...
        wt.mixin.ErrorHandling
    % Base class for a graphical widget

    % Copyright 2020-2023 The MathWorks Inc.
    

    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)
        
        % The internal grid to manage contents
        Grid matlab.ui.container.GridLayout
        
    end %properties


    properties (Transient, NonCopyable, Hidden, SetAccess = private)

        % Internal flag to confirm setup has finished
        SetupFinished (1,1) logical = false
        
    end %properties


    properties (Transient, NonCopyable, UsedInUpdate = true, ...
            GetAccess = private, SetAccess = protected)
        
        % Internal flag to trigger an update call
        Dirty (1,1) logical = false
        
    end %properties


    %% Special support for pre-R2021a releases (g2282435)
    properties (Transient, Access = private)
        
        % Flag to indicate earlier release (to be updated during setup)
        IsPre21a (1,1) logical = verLessThan('matlab','9.10')
        
    end %properties

    
    %% Debugging Methods
    methods
        
        function forceUpdate(obj)
            % Forces update to run (For debugging only!)

            disp("DEBUG: Forcing update for " + class(obj));
            obj.update();

        end %function
        
    end %methods
    

    %% Constructor
    methods

        function obj = BaseWidget(varargin)

            % Attach internal postSetup callback
            args = horzcat(varargin, {"CreateFcn",  @(src,evt)postSetup_I(src)});

            % Call superclass constructor
            obj = obj@matlab.ui.componentcontainer.ComponentContainer(args{:});

        end %function

    end %methods
    

    %% Protected Methods
    methods (Access = protected)
        
        function setup(obj)
            % Configure the widget

            % Grid Layout to manage building blocks
            obj.Grid = uigridlayout(obj);
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowHeight = {'1x'};
            obj.Grid.RowSpacing = 2;
            obj.Grid.ColumnSpacing = 2;
            obj.Grid.Padding = [0 0 0 0];
            
        end %function


        function postSetup(~)
            % Optional post-setup method 
            % (after setup and input arguments set, before update)
            
        end %function
        
        
        function requestUpdate(obj)
            % Request update method to run 
            % (This is for backward compatibility. 
            % Alternatively just set "obj.Dirty = true")

            % Support for initial release before R2021a
            if obj.IsPre21a && obj.SetupFinished

                % g2282435 Force a call to update
                obj.update();

            else

                % Trigger property to request update during next drawnow
                % (for optimal efficiency)
                obj.Dirty = true;

            end %if
            
        end %function
        

        function updateBackgroundColorableComponents(obj)
            % Update components that are affected by BackgroundColor
            % (overrides the superclass method)
            
            % Update grid color
            obj.Grid.BackgroundColor = obj.BackgroundColor;

            % Call superclass method
            obj.updateBackgroundColorableComponents@wt.mixin.BackgroundColorable();
            
        end %function


        function groups = getPropertyGroups(obj)
            % Customize the property display
            % (override to use the mixin implementation, since multiple
            % superclasses have competing implementations)

            groups = getPropertyGroups@wt.mixin.PropertyViewable(obj);
            
        end %function
        
    end %methods


    %% Private Methods
    methods (Access = private)
            
        function postSetup_I(obj)
            % Indicate setup is complete
            
            obj.SetupFinished = true;
            obj.CreateFcn = '';

            % Call any custom postSetup method
            obj.postSetup();
            
        end %function
        
    end %methods

    
end %classdef