classdef SliderCheckboxChangedData < event.EventData & dynamicprops
    % Event data for widget property value changes
    %
    % Syntax:
    %           obj = wt.eventdata.SliderCheckboxChangedData(name,index,prop,state,value)
    %
    
    % Copyright 2020-2021 The MathWorks, Inc.

    %% Properties
    properties (SetAccess = protected)
        Name (1,1) string
        Index (1,1) double
        Property (1,1) string
        State (1,1) logical
        Value (1,1) double
    end %properties


    %% Constructor / destructor
    methods
        function obj = SliderCheckboxChangedData(name,index,prop,state,value)
            
            % Set the changed properties
            obj.Name = name;
            obj.Index = index;
            obj.Property = prop;
            obj.State = state;
            obj.Value = value;

        end %constructor
    end %methods

end % classdef