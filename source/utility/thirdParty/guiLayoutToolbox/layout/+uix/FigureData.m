classdef ( Hidden, Sealed ) FigureData < event.EventData
    %uix.FigureData  Event data for FigureChanged on uix.FigureObserver
    
    %  Copyright 2014-2015 The MathWorks, Inc.
    %  $Revision: 1435 $ $Date: 2016-11-17 17:50:34 +0000 (Thu, 17 Nov 2016) $
    
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