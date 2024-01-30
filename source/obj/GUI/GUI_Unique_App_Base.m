classdef GUI_Unique_App_Base < matlab.apps.AppBase
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Access = protected)

        % Constructor
        function app = getUIFigure(app)
            
            app_instance = app.getAppInstance();
            % Check if there's already an instance running
            if ~isempty(app_instance) && isvalid(app_instance)
                % Bring the running instance to the front
                app = app_instance;
            else
                % No instance running, so create a new one
                createComponents(app);
                % Register the app with App Designer
                registerApp(app, app.UIFigure)
            end
            % Set the instance
            app.getAppInstance(app);
        end
    end

    %% METHODS UTILITY
    % ==================================================================================================================================================
    methods (Abstract)
        % Create UIFigure and components
        createComponents(app)
    end    
    methods (Static)
        function app = getAppInstance(app)
            persistent app_instance;
            if nargin == 1
                app_instance = app;
            end
            app = app_instance;
        end
    end

    methods
        function bringToFront(app)
            % Bring th app to front
            try
                % Not optimal solution suppose the name of the figure to be UIFigure

                % Get all properties of the app
                props = properties(app);

                % Loop through properties to find visible UIFigures
                for i = 1:length(props)
                    fig = app.(props{i});
                    if ~isempty(fig) && isa(fig, 'matlab.ui.Figure') && strcmp(fig.Visible, 'on')
                        figure(fig);
                    end                    
                end
            catch ex
                Core_Utils.printEx(ex);
            end
        end
    end

end
