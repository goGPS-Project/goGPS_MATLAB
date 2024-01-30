classdef BaseComponentDialog < wt.abstract.BaseDialog
    % Base dialog attached to a ComponentContainer component

    % Copyright 2022-2023 The MathWorks Inc.


    %% Read-Only Properties
    properties (SetAccess = {?matlab.ui.componentcontainer.ComponentContainer})

        % The parent component (must be provided to constructor)
        ParentComponent matlab.ui.componentcontainer.ComponentContainer

    end %properties


    %% Constructor / Destructor
    methods
        function obj = BaseComponentDialog(parentComponent)
            % Construct a dialog

            arguments
                parentComponent (1,1) matlab.ui.componentcontainer.ComponentContainer
            end

            % Get ancestor figure, if possible
            % If no figure available, place dialog in a new figure
            parent = ancestor(parentComponent,'figure');
            if isempty(parent)
                fig = uifigure;
                parent = uigridlayout(fig,[1 1]);
            end

            % Call superclass
            obj@wt.abstract.BaseDialog(parent, ...
                "ParentComponent", parentComponent);

            % Position within figure
            obj.positionWithin(parent, parentComponent)

            % Set dialog lifecycle to that of the parent component
            % The dialog will be deleted if the parent component is deleted
            obj.attachLifecycleListeners(parentComponent);

            % Pass along any fonts from the reference component
            obj.copyStylesFromReference(parentComponent);


        end %function
        
    end %methods


    %% Protected methods
    methods (Access = protected)
        
        function setup(obj)
            % Configure the widget

            % Call superclass method
            obj.setup@wt.abstract.BaseDialog();
            
        end %function
        
    end %methods


    %% Public Methods
    methods

        function positionWithinComponent(obj)

            % Position within figure parent
            obj.positionWithin(obj.Parent, obj.ParentComponent)

        end %function
        
    end %methods

end %classdef