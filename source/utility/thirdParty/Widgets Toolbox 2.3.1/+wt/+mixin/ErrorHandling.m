classdef ErrorHandling < handle
    %ErrorHandling Error handling methods
    
    % Copyright 2020-2023 The MathWorks Inc.
    
    methods ( Access = protected )
        
        function throwError(obj,err,title)
            % Throws an error to a dialog
            
            % Validate arguments
            arguments
                obj (1,1) wt.mixin.ErrorHandling
                err % string or MException
                title (1,1) string = "Error"
            end
            
            % Prepare the message
            % Was an exception provided?
            if isa(err,'MException')
                message = err.message;
            else
                message = string(err);
            end
            
            % Locate ancestor figure
            if isprop(obj,"Figure")
                fig = obj.Figure;
            else
                fig = ancestor(obj,'figure');
            end
            
            % Place in a dialog if possible
            if ~isempty(fig)
                uialert(fig,message,title);
            elseif isa(err,'MException')
                err.throwAsCaller();
            else
                error(message);
            end
            
        end %function
        
        
        function dlg = showProgress(obj,title,message,cancelOn)
            % Places a progress dialog in the widget's figure
            
            % Validate arguments
            arguments
                obj (1,1) wt.mixin.ErrorHandling
                title (1,1) string = "Please Wait"
                message (1,1) string = ""
                cancelOn (1,1) logical = false
            end
            
            % Locate ancestor figure
            if isprop(obj,"Figure")
                fig = obj.Figure;
            else
                fig = ancestor(obj,'figure');
            end
            
            % Place in a dialog if possible
            if isempty(fig)
                dlg = matlab.ui.dialog.ProgressDialog.empty(0,0);
            else
                dlg = uiprogressdlg(fig,...
                    "Title",title,...
                    "Message",message,...
                    "Cancelable",cancelOn);
            end
            
        end %function
        
        
        function dlg = showIndeterminateProgress(obj,title,message,cancelOn)
            % Places an indeterminate progress dialog in the widget's figure
            
            dlg = showProgress(obj,title,message,cancelOn);
            dlg.Indeterminate = true;
            
        end %function
        
    end%methods
    
end %classdef 

