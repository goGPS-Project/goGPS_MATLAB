classdef TitleColorable < handle
    % Mixin for component with Title font color

    % Copyright 2020-2023 The MathWorks Inc.
    
    
    %% Properties
    properties (AbortSet)
        
        % Title color
        TitleColor (1,3) double {wt.validators.mustBeBetweenZeroAndOne} = [0 0 0]
        
    end %properties
    
    
    
    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls to apply to
        TitleColorableComponents (:,1) matlab.graphics.Graphics
        
    end %properties
    
    
    
    %% Accessors
    methods
        
        function set.TitleColor(obj,value)
            obj.TitleColor = value;
            obj.updateTitleColorableComponents()
        end
        
        function set.TitleColorableComponents(obj,value)
            obj.TitleColorableComponents = value;
            obj.updateTitleColorableComponents()
        end
        
    end %methods
    
    
    
    %% Methods
    methods (Access = protected)
        
        function updateTitleColorableComponents(obj)

            % What needs to be updated?
            comps = obj.TitleColorableComponents;
            newValue = obj.TitleColor;
            propNames = ["FontColor","ForegroundColor"];

            % Set the subcomponent properties in prioritized order
            wt.utility.setStylePropsInPriority(comps, propNames, newValue);
            
        end %function
        
    end %methods
    
end %classdef