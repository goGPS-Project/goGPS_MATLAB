classdef Preferences < wt.model.BaseModel
    % Base class for managing app preferences
    %
    % Notes:
    %   The preferences may be set to automatically load upon creating an
    %   instance, and save upon destruction of an instance of this class.
    
    % Copyright 2020-2021 The MathWorks Inc.
    
    
    
    %% Preference Properties
    properties (AbortSet)
        
        % Starting window position
        Position (1,4) double = [100 100 1000 700] 
        
        % Last folder that was used when opening a file
        LastFolder (1,1) string {mustBeFolder} = pwd 
        
        % List of recent session files
        RecentSessionPaths (:,1) string {mustBeFile}
        
    end %properties
    
    
    
    %% Saving and Loading Preferences
    methods (Sealed)
        
        function load(obj,prefGroup)
            % Load stored preferences
            
            % Get a list of preference properties
            propList = getPreferenceProperties(obj);
            
            % Get the preference values
            if ispref(prefGroup)
                prefValues = getpref(prefGroup);
            else
                prefValues = struct();
            end
            
            % Loop on each preference field
            for idx = 1:numel(propList)
                
                % Get the current preference
                thisProp = propList(idx).Name;
                thisDefault = propList(idx).DefaultValue;
                if isfield(prefValues,thisProp)
                    thisValue = prefValues.(thisProp);
                else
                    thisValue = thisDefault;
                end
                
                % Get the preference
                try
                    obj.(thisProp) = thisValue;
                catch err
                    warning('wt:Preferences:LoadFailure',...
                        ['Unable to load stored preference for %s:%s. '...
                        'Reverting to default value.\n%s'],...
                        prefGroup,thisProp,err.message);
                    obj.(thisProp) = thisDefault;
                end
                
            end %for idx = 1:numel(propList)
            
        end %function
        
        
        function save(obj,prefGroup)
            % Save preferences
            
            % Get a list of preference properties
            propList = getPreferenceProperties(obj);
            
            % Remove existing preferences first (in case the props changed)
            if ispref(prefGroup)
                rmpref(prefGroup);
            end
            
            % Loop on each preference field
            if ~isempty(prefGroup)
                for idx = 1:numel(propList)
                    
                    % Get the current preference
                    thisProp = propList(idx).Name;
                    thisValue = obj.(thisProp);
                    
                    % Set the preference
                    setpref(prefGroup,thisProp,thisValue);
                    
                end %for idx = 1:numel(propList)
            end %if ~isempty(prefGroup)
            
        end %function
        
    end %methods
    
    
    
    %% Private methods
    methods (Access=private)
        
        function propList = getPreferenceProperties(obj)
            % Get the info for preference properties
            
            %mc = metaclass(obj);
            %propList = mc.PropertyList;
            propNames = properties(obj);
            for idx = numel(propNames):-1:1
                propList(idx) = findprop(obj,propNames{idx});
            end
            propsToKeep = ~[propList.Constant] & ...
                ~[propList.Transient] & ...
                ~[propList.Hidden] & ...
                [propList.HasDefault] & ...
                strcmp({propList.GetAccess},'public') & ...
                strcmp({propList.SetAccess},'public');
            propList(~propsToKeep) = [];
            
        end %function
        
    end %methods
    
    
end % classdef