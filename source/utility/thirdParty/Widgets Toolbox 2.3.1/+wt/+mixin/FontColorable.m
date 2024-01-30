classdef FontColorable < handle
    % Mixin to add styles to a component

    % Copyright 2020-2023 The MathWorks Inc.


    %% Properties
    properties (AbortSet)

        % Font color
        FontColor (1,3) double {wt.validators.mustBeBetweenZeroAndOne} = [0 0 0]

    end %properties



    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls to apply to
        FontColorableComponents (:,1) matlab.graphics.Graphics

    end %properties



    %% Accessors
    methods

        function set.FontColor(obj,value)
            obj.FontColor = value;
            obj.updateFontColorableComponents()
        end

        function set.FontColorableComponents(obj,value)
            obj.FontColorableComponents = value;
            obj.updateFontColorableComponents()
        end

    end %methods



    %% Methods
    methods (Access = protected)

        function updateFontColorableComponents(obj)

            % What needs to be updated?
            comps = obj.FontColorableComponents;
            newValue = obj.FontColor;
            propNames = ["FontColor","ForegroundColor"];

            % Set the subcomponent properties in prioritized order
            wt.utility.setStylePropsInPriority(comps, propNames, newValue);

        end %function

    end %methods

end %classdef