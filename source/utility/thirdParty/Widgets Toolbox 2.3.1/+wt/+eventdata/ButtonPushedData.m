classdef ButtonPushedData < event.EventData & dynamicprops
    % Event data for widget button push
    %
    % Syntax:
    %           obj = wt.eventdata.ButtonPushedData(button)
    %           obj = wt.eventdata.ButtonPushedData(button, value)
    %           obj = wt.eventdata.ButtonPushedData(eventData)
    %
    
    % Copyright 2020-2021 The MathWorks, Inc.

    %% Properties
    properties (SetAccess = protected)
        Button matlab.graphics.Graphics
        Text char
        Tag char
        Value
    end %properties


    %% Constructor / destructor
    methods
        function obj = ButtonPushedData(buttonOrEvt, value)
            
            % Was eventdata provided?
            if isa(buttonOrEvt,'matlab.ui.eventdata.ButtonPushedData')
                % Push button eventdata
                
                obj.Button = buttonOrEvt.Source;
                
            elseif isa(buttonOrEvt,'matlab.ui.eventdata.ValueChangedData')
                % State button eventdata
                
                obj.Button = buttonOrEvt.Source;
                obj.Value = value;
                
            else
                
                obj.Button = buttonOrEvt;
                
            end
            
            % Populate other data
            if isprop(obj.Button,'Text')
                obj.Text = obj.Button.Text;
            end
            obj.Tag = obj.Button.Tag;
            if nargin >= 2
                obj.Value = value;
            end

        end %constructor
    end %methods

end % classdef