classdef ( Hidden, Sealed ) FigureData < event.EventData
    %uix.FigureData  Event data for FigureChanged on uix.FigureObserver
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    properties( SetAccess = private )
        OldFigure % old figure
        NewFigure % new figure
    end
    
    methods( Access = ?uix.FigureObserver )
        
        function obj = FigureData( oldFigure, newFigure )
            %uix.FigureData  Create event data
            %
            %  d = uix.FigureData(oldFigure,newFigure)
            
            obj.OldFigure = oldFigure;
            obj.NewFigure = newFigure;
            
        end % constructor
        
    end % methods
    
end % classdef