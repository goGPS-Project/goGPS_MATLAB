classdef (Abstract) BaseModel < handle & ...
        matlab.mixin.SetGetExactNames & ...
        wt.mixin.DisplayNonScalarObjectAsTable
    % Base model class for apps that provides:
    %   PV pairs assignment on construction
    %   Display nonscalar arrays in table format
    %   Public PropertyChanged event for properties tagged as SetObservable
    %     to enable apps/widgets to listen to model changes
    %    

    % Copyright 2020-2023 The MathWorks, Inc.
    

    
    %% Events
    events
        
        % Triggered when SetObservable properties are changed
        PropertyChanged 
        
    end %events
    
    
    
    %% Internal Properties
    properties (Transient, NonCopyable, Hidden, SetAccess = private)
        
        % Listeners to public properties
        PropListeners
        
    end %properties
    
    
    
    %% Constructor
    methods
        function obj = BaseModel(varargin)
            % Constructor
            
            % Populate public properties from P-V input pairs
            if nargin
                obj.assignPVPairs(varargin{:});
            end
            
            % Create listeners to public properties
            obj.createPropListeners();
            
        end %function obj = ObjectModel(varargin)
    end %methods
    
    
    
    %% Static methods
    methods (Static)
        
        function obj = loadobj(obj)
            % Customize loading from file
            
            if isstruct(obj)
                error('Unable to load object.');
            end
            
            % Need to recreate listeners
            obj.createPropListeners();
            
        end %function
        
    end %methods
    
    
    
    %% Protected Methods
    
    methods (Access=protected)
        
        % This method is similar to
        % matlab.io.internal.mixin.HasPropertiesAsNVPairs, but this one
        % generally performs faster.
        
        function varargout = assignPVPairs(obj,varargin)
            % Assign the specified property-value pairs
            
            if nargin > 1
                
                % Get a singleton parser for this class
                keepUnmatched = nargout > 0;
                p = getParser(obj, keepUnmatched);
                
                % Parse the P-V pairs
                p.parse(varargin{:});
                
                % Set just the parameters the user passed in
                ParamNamesToSet = varargin(1:2:end);
                ParamsInResults = fieldnames(p.Results);
                
                % Assign properties
                for ThisName = ParamNamesToSet
                    isSettable = any(strcmpi(ThisName,ParamsInResults));
                    if isSettable && ~any(strcmpi(ThisName,p.UsingDefaults))
                        obj.(ThisName{1}) = p.Results.(ThisName{1});
                    end
                end
                
                % Return unmatched pairs
                if nargout
                    varargout{1} = p.Unmatched;
                end
                
            elseif nargout
                
                varargout{1} = struct;
                
            end %if nargin > 1
            
        end %function
        
        
        function createPropListeners(obj)
            
            for idx = 1:numel(obj)
                mc = metaclass(obj(idx));
                isObservable = [mc.PropertyList.SetObservable];
                props = mc.PropertyList(isObservable);
                obj(idx).PropListeners = event.proplistener(obj(idx),props,...
                    'PostSet',@(h,e)onPropChanged(obj(idx),e) );
            end %for
            
        end %function
        
        
        function onPropChanged(obj,e)
            
            evt = wt.eventdata.PropertyChangedData(e.Source.Name, obj.(e.Source.Name));
            obj.notify('PropertyChanged',evt)
            
        end %function
        
    end %methods
    
    
    
    %% Private methods
    methods (Access=private)
        
        function thisParser = getParser(obj,keepUnmatched)
            
            % What class is this?
            className = class(obj);
            
            % Keep a list of reusable parsers for each class
            persistent allParsers
            if isempty(allParsers)
                allParsers = containers.Map('KeyType','char','ValueType','any');
            end
            
            % Get or make a custom parser for this class
            try
                
                thisParser = allParsers(className);
                
            catch
                
                % Get a list of public properties
                metaObj = metaclass(obj);
                isSettableProp = strcmp({metaObj.PropertyList.SetAccess}','public');
                settableProps = metaObj.PropertyList(isSettableProp);
                publicPropNames = {settableProps.Name}';
                hasDefault = [settableProps.HasDefault]';
                defaultValues = repmat({[]},size(hasDefault));
                defaultValues(hasDefault) = {settableProps(hasDefault).DefaultValue};
                
                % Create custom parser for this class
                thisParser = inputParser;
                thisParser.KeepUnmatched = keepUnmatched;
                thisParser.FunctionName = className;
                
                % Add each public property to the parser
                for pIdx = 1:numel(publicPropNames)
                    thisParser.addParameter(publicPropNames{pIdx}, defaultValues{pIdx});
                end
                
                % Add this parser to the map
                allParsers(className) = thisParser;
                
            end %if allParsers.isKey(className)
            
        end %function
        
    end %methods
        
    
end %classdef
