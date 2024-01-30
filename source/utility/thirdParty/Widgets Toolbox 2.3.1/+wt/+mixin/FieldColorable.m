classdef FieldColorable < handle
    % Mixin to add styles to a component

    % Copyright 2020-2023 The MathWorks Inc.


    %% Properties
    properties (AbortSet)

        % Field Color
        FieldColor (1,3) double {wt.validators.mustBeBetweenZeroAndOne} = [1 1 1]

    end %properties



    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls to apply to
        FieldColorableComponents (:,1) matlab.graphics.Graphics

    end %properties



    %% Accessors
    methods

        function set.FieldColor(obj,value)
            obj.FieldColor = value;
            obj.updateFieldColorableComponents()
        end

        function set.FieldColorableComponents(obj,value)
            obj.FieldColorableComponents = value;
            obj.updateFieldColorableComponents()
        end

    end %methods



    %% Methods
    methods (Access = protected)

        function updateFieldColorableComponents(obj)

            % What needs to be updated?
            comps = obj.FieldColorableComponents;
            newValue = obj.FieldColor;
            propNames = ["FieldColor","BackgroundColor","Color"];

            % Set the subcomponent properties in prioritized order
            wt.utility.setStylePropsInPriority(comps, propNames, newValue);

        end %function

    end %methods

end %classdef