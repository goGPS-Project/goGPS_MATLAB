classdef ( Hidden, Sealed ) FigureObserver < handle
    %uix.FigureObserver  Figure observer
    %
    %  A figure observer raises an event FigureChanged when the figure
    %  ancestor of a subject changes.
    
    %  Copyright 2014-2015 The MathWorks, Inc.
    %  $Revision: 1435 $ $Date: 2016-11-17 17:50:34 +0000 (Thu, 17 Nov 2016) $
    
    properties( SetAccess = private )
        Subject % subject
        Figure % figure ancestor
    end
    
    properties( Access = private )
        PreSetListeners % listeners to Parent PreGet
        PostSetListeners % listeners to Parent PreGet
        OldFigure = gobjects( 0 ) % previous figure ancestor
    end
    
    events( NotifyAccess = private )
        FigureChanged
    end
    
    methods
        
        function obj = FigureObserver( subject )
            %uix.FigureObserver  Create figure observer
            %
            %  o = uix.FigureObserver(s) creates a figure observer for the
            %  subject s.
            
            % Check
            validateattributes( subject, {'matlab.graphics.Graphics'}, ...
                {'scalar'}, '', 'subject' )
            
            % Store subject
            obj.Subject = subject;
            
            % Set up object
            obj.update()
            
        end % constructor
        
    end % structors
    
    methods( Access = private )
        
        function update( obj )
            %update  Update listeners and Figure property
            
            % Create fresh listeners
            obj.PreSetListeners = event.proplistener.empty( [1 0] ); % clear
            obj.PostSetListeners = event.proplistener.empty( [1 0] ); % clear
            o = obj.Subject;
            while ~isempty( o ) && ~isa( o, 'matlab.ui.Figure' )
                obj.PreSetListeners(end+1) = event.proplistener( o, ...
                    findprop( o, 'Parent' ), 'PreSet', @obj.onParentPreSet );
                obj.PostSetListeners(end+1) = event.proplistener( o, ...
                    findprop( o, 'Parent' ), 'PostSet', @obj.onParentPostSet );
                o = o.Parent;
            end
            
            % Store figure
            obj.Figure = o;
            
        end % update
        
        function onParentPreSet( obj, ~, ~ )
            %onParentPreSet  Event handler
            
            % Store old figure
            obj.OldFigure = obj.Figure;
            
        end % onParentPreSet
        
        function onParentPostSet( obj, ~, ~ )
            %onParentPostSet  Event handler
            
            % Update object
            obj.update()
            
            % Raise event
            oldFigure = obj.OldFigure;
            newFigure = obj.Figure;
            if ~isequal( oldFigure, newFigure )
                notify( obj, 'FigureChanged', ...
                    uix.FigureData( oldFigure, newFigure ) )
            end
            
            % Clear old figure
            obj.OldFigure = gobjects( 0 );
            
        end % onParentPostSet
        
    end % private methods
    
end % classdef