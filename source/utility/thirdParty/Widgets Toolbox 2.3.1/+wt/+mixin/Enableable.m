classdef Enableable < handle
    % Mixin to add styles to a component

    % Copyright 2020-2023 The MathWorks Inc.
    
    
    %% Properties
    properties (AbortSet)
        
        % Enable of the component
        Enable (1,1) matlab.lang.OnOffSwitchState = 'on'
        
    end %properties
    
    
    
    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls to apply to
        EnableableComponents (:,1) matlab.graphics.Graphics
        
    end %properties
    
    
    
    %% Accessors
    methods
        
        function set.Enable(obj,value)
            obj.Enable = value;
            obj.updateEnableableComponents()
        end
        
        function set.EnableableComponents(obj,value)
            obj.EnableableComponents = value;
            obj.updateEnableableComponents()
        end
        
    end %methods
    
    
    
    %% Methods
    methods (Access = protected)
        
        function updateEnableableComponents(obj)

            % What needs to be updated?
            comps = obj.EnableableComponents;
            newValue = obj.Enable;
            propNames = "Enable";

            % Set the subcomponent properties in prioritized order
            wt.utility.setStylePropsInPriority(comps, propNames, newValue);
            
        end %function
        
    end %methods
    
end %classdef