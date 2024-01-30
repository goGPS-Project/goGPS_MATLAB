classdef FontStyled < handle
    % Mixin for component with Font properties

    % Copyright 2020-2023 The MathWorks Inc.


    %% Properties
    properties (AbortSet)

        % Font name
        FontName char {mustBeNonempty} = 'Helvetica'

        % Font size in points
        FontSize (1,1) double {mustBePositive,mustBeFinite} = 12

        % Font weight (normal/bold)
        FontWeight {mustBeMember(FontWeight,{'normal','bold'})} = 'normal'

        % Font angle (normal/italic)
        FontAngle {mustBeMember(FontAngle,{'normal','italic'})} = 'normal'

        % Font color
        FontColor (1,3) double {wt.validators.mustBeBetweenZeroAndOne} = [0 0 0]

    end %properties



    %% Internal properties
    properties (AbortSet, Transient, NonCopyable, Hidden, SetAccess = protected)

        % List of graphics controls to apply to
        FontStyledComponents (:,1) matlab.graphics.Graphics

    end %properties



    %% Accessors
    methods

        function set.FontName(obj,value)
            obj.FontName = value;
            obj.updateFontStyledComponents("FontName",value)
        end

        function set.FontSize(obj,value)
            obj.FontSize = value;
            obj.updateFontStyledComponents("FontSize",value)
        end

        function set.FontWeight(obj,value)
            obj.FontWeight = value;
            obj.updateFontStyledComponents("FontWeight",value)
        end

        function set.FontAngle(obj,value)
            obj.FontAngle = value;
            obj.updateFontStyledComponents("FontAngle",value)
        end

        function set.FontColor(obj,value)
            obj.FontColor = value;
            obj.updateFontStyledComponents("FontColor",value)
        end

        function set.FontStyledComponents(obj,value)
            obj.FontStyledComponents = value;
            obj.updateFontStyledComponents();
        end

    end %methods



    %% Methods
    methods (Access = protected)

        function updateFontStyledComponents(obj,prop,value)

            % Get the components
            comps = obj.FontStyledComponents;

            % Updating all or a specific property?
            if nargin < 3

                % Font color properties in prioritized order
                colorProps = ["FontColor","ForegroundColor"];

                % Set all subcomponent properties
                wt.utility.setStylePropsInPriority(comps,"FontName",obj.FontName)
                wt.utility.setStylePropsInPriority(comps,"FontSize",obj.FontSize)
                wt.utility.setStylePropsInPriority(comps,"FontWeight",obj.FontWeight)
                wt.utility.setStylePropsInPriority(comps,"FontAngle",obj.FontAngle)
                wt.utility.setStylePropsInPriority(comps,colorProps, obj.FontColor);

            elseif prop == "FontColor"
                % Update just the FontColor property

                % Set the subcomponent property
                wt.utility.setStylePropsInPriority(comps,...
                    ["FontColor","ForegroundColor"], value);

            else

                % Set the subcomponent property
                wt.utility.setStylePropsInPriority(comps,prop,value);

            end %if

        end %function

    end %methods

end %classdef