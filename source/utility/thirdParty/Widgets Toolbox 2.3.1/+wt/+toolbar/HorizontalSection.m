classdef (Sealed, Hidden) HorizontalSection < wt.toolbar.BaseSection
    % A row of horizontally stacked controls for a toolbar
    
    % Copyright 2020-2021 The MathWorks Inc.
    
    
    %% Public properties
    properties (AbortSet)
        
        % The title to display beneath the section
        Title (1,1) string = "NEW SECTION"
        
        % Components that are part of this section
        Component
        
        % Width of each component
        ComponentWidth (1,:) double
        
        % Width of each section button
        MinimizedWidth (1,1) double = 60
        
    end %properties
    
    
    properties (Dependent, SetAccess = protected)
        
        % Width of the whole section
        TotalWidth (1,1) double = 0
        
    end %properties
    
    
    
    %% Read-Only Properties
    properties (SetAccess = protected)
        
        % Listen to button pushes in sections
        ButtonPushedListener event.listener
        
    end %properties
    
    
    properties (Dependent, SetAccess = protected)
        
        % Indicates which components are vertical sections
        IsVerticalSection (:,1) logical
        
    end %properties
    
    
    
    %% Public methods
    methods
        
        function section = addVerticalSection(obj)
            % Adds a vertical section to the toolbar
            
            % Add a vertical subsection
            section = wt.toolbar.VerticalSection();
            
            % Add the component
            obj.Component(end+1) = section;
            
        end %function
        
    end %methods
    
    
    
    %% Protected methods
    methods (Access = protected)
        
        function setup(obj)
            
            % Call superclass setup first to establish the grid
            obj.setup@wt.toolbar.BaseSection();
            
            % Configure grid
            obj.Grid.Padding = [8 4 8 4];
            obj.Grid.ColumnSpacing = 4;
            
        end %function
        
        
        function update(obj)
            
            % Call superclass method
            obj.update@wt.toolbar.BaseSection();
            
            % Configure the layout
            obj.Grid.ColumnWidth = num2cell(obj.ComponentWidth);
            obj.Grid.RowHeight = {'1x'};
            for idx = 1:numel(obj.Component)
                obj.Component(idx).Layout.Row = 1;
                obj.Component(idx).Layout.Column = idx;
            end
            
            % Listen to any vertical sections
            isVerticalSection = arrayfun(@(x)isa(x,'wt.toolbar.VerticalSection'),obj.Component);
            obj.ButtonPushedListener = event.listener( ...
                obj.Component(isVerticalSection), 'ButtonPushed', ...
                @(h,e)onButtonPushed(obj,e) );
            
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
            
            % Add section name to event data
            addprop(evt,'Section');
            evt.Section = obj.Title;
            
            % Trigger event
            notify(obj,"ButtonPushed",evt);
            
            % Trigger callback
            if ~isempty(obj.ButtonPushedFcn)
                feval(obj.ButtonPushedFcn, obj, evt);
            end
            
        end %function
        
    end %methods
    
    
    
    %% Accessors
    methods
        
        function value = get.ComponentWidth(obj)
            % Expand/contract length to the number of components
            numComp = numel(obj.Component);
            value = repmat(50,1,numComp);
            value(obj.IsVerticalSection) = 90;
            numSet = min(numComp, numel(obj.ComponentWidth));
            value(1:numSet) = obj.ComponentWidth(1:numSet);
        end
        
        function set.ComponentWidth(obj,value)
            evt = wt.eventdata.PropertyChangedData('ComponentWidth',value,obj.ComponentWidth);
            obj.ComponentWidth = value;
            notify(obj,'PropertyChanged',evt)
        end
        
        
        function value = get.TotalWidth(obj)
            numComp = numel(obj.ComponentWidth);
            value = sum(obj.Grid.Padding([1 3])) + ...
                sum(obj.ComponentWidth) + ...
                obj.Grid.ColumnSpacing * (numComp - 1);
        end
        
        
        function value = get.IsVerticalSection(obj)
            if isempty(obj.Component)
                value = false(0,1);
            else
                types = get(obj.Component,'Type');
                value = types == lower("wt.toolbar.VerticalSection");
            end
        end
        
        
        function set.Title(obj,value)
            evt = wt.eventdata.PropertyChangedData('Title',value,obj.Title);
            obj.Title = value;
            notify(obj,'PropertyChanged',evt)
        end
        
        
        function set.MinimizedWidth(obj,value)
            evt = wt.eventdata.PropertyChangedData('MinimizedWidth',value,obj.MinimizedWidth);
            obj.MinimizedWidth = value;
            notify(obj,'PropertyChanged',evt)
        end
        
    end %methods
    
end % classdef