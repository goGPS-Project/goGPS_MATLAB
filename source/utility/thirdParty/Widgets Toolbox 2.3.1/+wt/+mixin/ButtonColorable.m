classdef ButtonColorable < handle
    % Mixin to add styles to a component

    % Copyright 2020-2023 The MathWorks Inc.


    %% Properties
    properties (AbortSet)

        % Button Color
        ButtonColor (1,3) double {wt.validators.mustBeBetweenZeroAndOne} = [1 1 1] * 0.96

    end %properties



    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls to apply to
        ButtonColorableComponents (:,1) matlab.graphics.Graphics

    end %properties



    %% Accessors
    methods

        function set.ButtonColor(obj,value)
            obj.ButtonColor = value;
            obj.updateButtonColorableComponents()
        end

        function set.ButtonColorableComponents(obj,value)
            obj.ButtonColorableComponents = value;
            obj.updateButtonColorableComponents()
        end

    end %methods



    %% Methods
    methods (Access = protected)

        function updateButtonColorableComponents(obj)

            % What needs to be updated?
            comps = obj.ButtonColorableComponents;
            newValue = obj.ButtonColor;
            propNames = ["ButtonColor","BackgroundColor"];

            % Set the subcomponent properties in prioritized order
            wt.utility.setStylePropsInPriority(comps, propNames, newValue);

        end %function

    end %methods

end %classdef