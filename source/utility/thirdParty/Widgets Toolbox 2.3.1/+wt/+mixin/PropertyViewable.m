classdef PropertyViewable < handle
    % Mixin for component to organize the display of properties in the
    % command window.
    
    % Copyright 2020-2023 The MathWorks Inc.

    
    %% Protected Methods
    methods (Access = protected)

        function groups = getPropertyGroups(obj)
            % Customize how the properties are displayed

            % Ignore most superclass properties for default display
            persistent superProps
            if isempty(superProps)
                superProps = properties('matlab.ui.componentcontainer.ComponentContainer');
            end

            % Get the relevant properties (ignore Superclass properties)
            allProps = properties(obj);
            propNames = setdiff(allProps, superProps, 'stable');

            % Remove properties we don't need to see
            propNames(startsWith(propNames, "Font")) = [];
            propNames(matches(propNames, "Enable")) = [];

            % Split callbacks, fonts, and colorizations
            isCallback = endsWith(propNames, "Fcn");
            isColor = endsWith(propNames, "Color");
            normalProps = propNames(~isCallback);
            callbackProps = propNames(isCallback & ~isColor);

            % Define the property gorups
            groups = [
                matlab.mixin.util.PropertyGroup(callbackProps)
                matlab.mixin.util.PropertyGroup(normalProps)
                matlab.mixin.util.PropertyGroup(["Position", "Units"])
                ];

            % Ignore Empty Groups
            groups(~[groups.NumProperties]) = [];

        end %function
        
    end %methods
    
end %classdef