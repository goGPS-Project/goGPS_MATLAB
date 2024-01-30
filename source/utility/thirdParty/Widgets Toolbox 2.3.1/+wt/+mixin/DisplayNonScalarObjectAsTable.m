classdef (HandleCompatible) DisplayNonScalarObjectAsTable < matlab.mixin.CustomDisplay
    % DisplayNonScalarObjectAsTable -
    % ---------------------------------------------------------------------
    % Abstract: This mixin class overrides displayNonScalarObject() to
    % display an array of objects as a table.
    %
    
    %   Copyright 2018-2019 The MathWorks, Inc.
    %
    % Auth/Revision:
    %   MathWorks Consulting
    %   $Author: rjackey $
    %   $Revision: 348 $  $Date: 2018-03-02 15:51:54 -0500 (Fri, 02 Mar 2018) $
    % ---------------------------------------------------------------------
    
    
    %% Public Methods
    methods (Sealed)
        
        function t = toDisplayTable(obj)
            % Convert the object to a table
            
            % Prepare indices as row names
            numObj = numel(obj);
            if isrow(obj)
                indices = string(1:numObj)';
            else
                varargout = cell(1,ndims(obj));
                [varargout{:}] = ind2sub(size(obj),1:numObj);
                indices = string(vertcat(varargout{:})');
                indices = join(indices,',',2);
            end
            rowNames = "(" + indices + ")";
            
            % Check for deleted handle objects
            if isa(obj,'handle')
                isDeleted = ~isvalid(obj);
                obj(isDeleted) = [];
                rowNames(isDeleted) = [];
            end
            
            % Gather property info
            props = string( properties(obj) )';
            
            % Preallocate a cell to store values
            numProps = numel(props);
            numObj = numel(obj);
            values = cell(numObj,numProps);
            
            % Populate cell with property values
            for pIdx = 1:numProps
                thisProp = props(pIdx);
                values(:,pIdx) = {obj.(thisProp)}';
            end %for
            
            % Convert to table
            t = cell2table(values,'VariableNames',props,'RowNames',rowNames);
            
        end %function
        
    end %methods
    
    
    %% Protected Methods
    methods (Sealed, Access=protected)
        
        function displayNonScalarObject(obj)
            
            % Format text to display
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            dimStr = matlab.mixin.CustomDisplay.convertDimensionsToString(obj);
            
            % Display the header
            if isa(obj,'matlab.mixin.Heterogeneous')
                fprintf('  %s Heterogeneous %s with common properties:\n\n',dimStr,className);
            else
                fprintf('  %s %s with properties:\n\n',dimStr,className);
            end
            
            % Show the group list in a table
            disp( obj.toDisplayTable() );
            
            if isa(obj,'handle') && any(~isvalid(obj))
               fprintf('  Object array contains deleted handles.\n'); 
            end
            
        end %function
        
    end %methods
    
    
end % classdef
