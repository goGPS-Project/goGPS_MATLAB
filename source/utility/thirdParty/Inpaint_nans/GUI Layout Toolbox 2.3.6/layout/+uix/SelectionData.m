classdef( Hidden, Sealed ) SelectionData < event.EventData
    %uix.SelectionData  Event data for selection event
    %
    %  e = uix.SelectionData(o,n) creates event data including the old
    %  value o and the new value n.
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    properties( SetAccess = private )
        OldValue % old value
        NewValue % newValue
    end
    
    methods
        
        function obj = SelectionData( oldValue, newValue )
            %uix.SelectionData  Event data for selection event
            %
            %  e = uix.SelectionData(o,n) creates event data including the
            %  old value o and the new value n.
            
            % Set properties
            obj.OldValue = oldValue;
            obj.NewValue = newValue;
            
        end % constructor
        
    end % structors
    
end % classdef