classdef (Hidden, AllowedSubclasses = ...
        {?wt.toolbar.VerticalSection, ?wt.toolbar.HorizontalSection} ) ...
        BaseSection < wt.abstract.BaseWidget & wt.mixin.FontStyled
    % Base class for a toolbar section
    
    % Copyright 2020-2021 The MathWorks Inc.
    
    
    %% Events
    %events (HasCallbackProperty, NotifyAccess = protected)
    %RAJ - manually added callback due to g2405728
    events (NotifyAccess = protected)
        
        % Event triggered when a button is pushed
        ButtonPushed
        
    end %events
        
    
    events (ListenAccess = ?wt.Toolbar, NotifyAccess = protected)
        
        % Triggered when properties that affect toolbar display are changed
        PropertyChanged 
        
    end %events
    
    
    
    %% Public properties
    properties (Abstract, AbortSet)
        
        % Components that are part of this section
        Component (:,1) matlab.graphics.Graphics
        
    end %properties
    
    
    properties (UsedInUpdate = false)
        
        % Callback
        ButtonPushedFcn function_handle {mustBeScalarOrEmpty}
        
    end %properties
    
    
    %% Constructor
    methods
        
        function obj = BaseSection(varargin)
            % Construct the widget
            
            % If no inputs provided, leave widget unparented
            if ~nargin
                varargin = {'Parent',[]};
            end
            
            % Call superclass constructor
            obj@wt.abstract.BaseWidget(varargin{:});
            
        end %function
        
    end %constructor/destructor
    
    
    
    %% Abstract methods
    methods (Abstract, Access = protected)
        
        onButtonPushed(obj,evt)
        
    end %methods
    
    

    %% Public methods
    methods
        
        function button = addButton(obj,icon,text)
            % Adds a push button to the toolbar
        
            arguments
                obj (1,1) wt.toolbar.BaseSection
                icon char = ''
                text char = ''
            end
                
            % Create the button
            button = uibutton('push','Parent',[]);
            button.ButtonPushedFcn = @(h,e)onButtonPushed(obj,e);
            button.Text = text;
            button.Icon = icon;
            button.IconAlignment = 'top';
            button.WordWrap = 'on';
            
            obj.Component(end+1) = button;
            
        end %function
        
        
        function button = addStateButton(obj,icon,text)
            % Adds a state button to the toolbar
        
            arguments
                obj (1,1) wt.toolbar.BaseSection
                icon char = ''
                text char = ''
            end
                
            % Create the button
            button = uibutton('state','Parent',[]);
            button.ValueChangedFcn = @(h,e)onButtonPushed(obj,e);
            button.Text = text;
            button.Icon = icon;
            button.IconAlignment = 'top';
            button.WordWrap = 'on';
            
            % Add the component
            obj.Component(end+1) = button;
            
        end %function
        
    end %methods
    
    
    
    %% Protected methods
    methods (Access = protected)
        
        function update(obj)
            
            % Unparent removed components
            removedComponents = setdiff(obj.Grid.Children, obj.Component);
            removedComponents = removedComponents(isvalid(removedComponents));
            set(removedComponents,'Parent',[]);
            
            % Link styles
            obj.FontStyledComponents = obj.Component;
            obj.BackgroundColorableComponents = obj.Component;
            
            % Set column width
            obj.Grid.ColumnWidth = {'1x'};
            obj.Grid.RowHeight = repmat({'fit'},1,numel(obj.Component));

            % Parent new components
            set(obj.Component,'Parent',obj.Grid);
            
        end %function
        
    end %methods
    
    
end % classdef