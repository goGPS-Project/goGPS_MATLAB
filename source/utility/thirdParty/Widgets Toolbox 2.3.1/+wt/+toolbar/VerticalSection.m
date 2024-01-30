classdef (Sealed, Hidden) VerticalSection < wt.toolbar.BaseSection
    % A column of vertically stacked controls for a toolbar
    
    % Copyright 2020-2021 The MathWorks Inc.
    
    
    %% Public properties
    properties (AbortSet)
        
        % Components that are part of this section
        Component
        
    end %properties
    
    
    
    %% Public methods
    methods
        
        
        function button = addButton(obj,icon,text)
            
            % Call superclass method first
            button = obj.addButton@wt.toolbar.BaseSection(icon,text);
            
            % Customize for vertical stack
            button.IconAlignment = 'left';
            button.WordWrap = 'off';
            
        end %function
        
        
        function button = addStateButton(obj,icon,text)
            
            % Call superclass method first
            button = obj.addStateButton@wt.toolbar.BaseSection(icon,text);
            
            % Customize for vertical stack
            button.IconAlignment = 'left';
            button.WordWrap = 'off';
            
        end %function
        
    end %methods
    
    
    %% Protected methods
    methods (Access = protected)
        
        function setup(obj)
            
            % Call superclass setup first to establish the grid
            obj.setup@wt.abstract.BaseWidget();
            
            % Configure grid
            obj.Grid.RowHeight = {};
            obj.Grid.ColumnWidth = {'1x'};
            
        end %function
        
        
        function update(obj)
            
            % Call superclass method
            obj.update@wt.toolbar.BaseSection();
            
            % Configure the layout
            obj.Grid.ColumnWidth = {'1x'};
            if numel(obj.Grid.RowHeight) ~= numel(obj.Component)
                obj.Grid.RowHeight = repmat({'fit'},1,numel(obj.Component));
            end
            for idx = 1:numel(obj.Component)
                obj.Component(idx).Layout.Column = 1;
                obj.Component(idx).Layout.Row = idx;
            end
            
        end %function
        
        
        function onButtonPushed(obj,evt)
            
            % Use custom event data
            if isa(evt,'matlab.ui.eventdata.ButtonPushedData')
                % Push button
                evt = wt.eventdata.ButtonPushedData(evt.Source);
            elseif isa(evt,'matlab.ui.eventdata.ValueChangedData')
                % State button
                evt = wt.eventdata.ButtonPushedData(evt.Source, evt.Source.Value);
            end
            
            % Trigger event
            notify(obj,"ButtonPushed",evt);
            
            % Trigger callback
            if ~isempty(obj.ButtonPushedFcn)
                feval(obj.ButtonPushedFcn, obj, evt);
            end
            
        end %function
        
    end %methods
    
    
end % classdef