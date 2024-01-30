classdef ValueChangedData < event.EventData & dynamicprops
    % Event data for widget value changes
    %
    % Syntax:
    %           obj = wt.eventdata.ValueChangedData(newValue,oldValue)
    %           obj = wt.eventdata.ValueChangedData(valueChangedEvent)
    %           obj = wt.eventdata.ValueChangedData(...,'p1',v1,...)
    %
    
    % Copyright 2020-2021 The MathWorks, Inc.

    %% Properties
    properties (SetAccess = protected)
        Value
        PreviousValue
    end %properties


    %% Constructor / destructor
    methods
        function obj = ValueChangedData(newValue, previousValue)
            
            arguments
                newValue
                previousValue = []
            end
            
            % Is input a MATLAB eventdata?
            if isa(newValue,'matlab.ui.eventdata.ValueChangedData')
                obj.Value = newValue.Value;
                obj.PreviousValue = newValue.PreviousValue;
            elseif isa(newValue,'matlab.ui.eventdata.ValueChangingData')
                obj.Value = newValue.Value;
            else
                % No - use the value directly
                obj.Value = newValue;
                obj.PreviousValue = previousValue;
            end
            
        end %constructor
    end %methods

end % classdef