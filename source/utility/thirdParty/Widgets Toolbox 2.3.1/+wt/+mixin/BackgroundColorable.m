classdef BackgroundColorable < handle
    % Mixin to add styles to a component

    % Copyright 2020-2023 The MathWorks Inc.


    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls to apply to
        BackgroundColorableComponents (:,1) matlab.graphics.Graphics

        % Listener to background color changes
        BackgroundColorListener event.proplistener

    end %properties


    %% Accessors
    methods

        function set.BackgroundColorableComponents(obj,value)
            obj.BackgroundColorableComponents = value;
            obj.updateBackgroundColorableComponents()
            obj.listenForBackgroundChange();
        end

    end %methods



    %% Methods
    methods (Access = protected)

        function updateBackgroundColorableComponents(obj)

            % What needs to be updated?
            comps = obj.BackgroundColorableComponents;
            newValue = obj.BackgroundColor; %#ok<MCNPN> 
            propNames = ["BackgroundColor","Color"];

            % Set the subcomponent properties in prioritized order
            wt.utility.setStylePropsInPriority(comps, propNames, newValue);

        end %function


        function listenForBackgroundChange(obj)

            % Establish Listener for Background Color Change
            if isempty(obj.BackgroundColorListener)
                obj.BackgroundColorListener = ...
                    addlistener(obj,'BackgroundColor','PostSet',...
                    @(h,e)obj.updateBackgroundColorableComponents());
            end

        end %function

    end %methods

end %classdef