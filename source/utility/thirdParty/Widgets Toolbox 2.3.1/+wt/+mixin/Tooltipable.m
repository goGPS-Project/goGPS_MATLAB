classdef Tooltipable < handle
    % Mixin for component with Tooltip property

    % Copyright 2020-2023 The MathWorks Inc.
    
    
    %% Properties
    properties (AbortSet)
        
        % Tooltip of the component
        Tooltip (1,1) string
        
    end %properties
    
    
    
    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls to apply to
        TooltipableComponents (:,1) matlab.graphics.Graphics
        
    end %properties
    
    
    
    %% Accessors
    methods
        
        function set.Tooltip(obj,value)
            obj.Tooltip = value;
            obj.updateTooltipableComponents()
        end
        
        function set.TooltipableComponents(obj,value)
            obj.TooltipableComponents = value;
            obj.updateTooltipableComponents()
        end
        
    end %methods
    
    
    
    %% Methods
    methods (Access = protected)
        
        function updateTooltipableComponents(obj)

            % What needs to be updated?
            comps = obj.TooltipableComponents;
            newValue = obj.Tooltip;
            propNames = "Tooltip";

            % Set the subcomponent properties in prioritized order
            wt.utility.setStylePropsInPriority(comps, propNames, newValue);
            
        end %function
        
    end %methods
    
end %classdef