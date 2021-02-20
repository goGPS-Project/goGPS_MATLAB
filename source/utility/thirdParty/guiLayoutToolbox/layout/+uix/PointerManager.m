classdef ( Hidden, Sealed ) PointerManager < handle
    %uix.PointerManager  Pointer manager
    
    %  Copyright 2016 The MathWorks, Inc.
    %  $Revision: 1435 $ $Date: 2016-11-17 17:50:34 +0000 (Thu, 17 Nov 2016) $
    
    properties( SetAccess = private )
        Figure % figure
    end
    
    properties( Access = private )
        Tokens % tokens
        Pointers % pointers
        NextToken % next token
        PointerListener % listener
    end
    
    methods( Access = private )
        
        function obj = PointerManager( figure )
            %uix.PointerManager  Create pointer manager
            %
            %  m = uix.PointerManager(f) creates a pointer manager for the
            %  figure f.
            
            obj.Figure = figure;
            obj.Tokens = 0;
            obj.Pointers = {figure.Pointer};
            obj.NextToken = 1;
            obj.PointerListener = event.proplistener( figure, ...
                findprop( figure, 'Pointer' ), 'PostSet', ...
                @obj.onPointerChanged );
            
        end % constructor
        
    end % structors
    
    methods( Access = private )
        
        function doSetPointer( obj, token, pointer )
            %doSetPointer  Set pointer
            %
            %  m.doSetPointer(t,p) sets the pointer to p with the token t.
            
            % Remove old entry
            tf = obj.Tokens == token;
            obj.Tokens(tf) = [];
            obj.Pointers(tf) = [];
            
            % Add new entry
            obj.Tokens(end+1) = token;
            obj.Pointers{end+1} = pointer;
            
            % Set pointer
            obj.PointerListener.Enabled = false;
            obj.Figure.Pointer = pointer;
            obj.PointerListener.Enabled = true;
            
        end % doSetPointer
        
        function doUnsetPointer( obj, token )
            %doUnsetPointer  Unset pointer
            %
            %  m.doUnsetPointer(s) unsets the pointer with the token t.
            
            % Remove old entry
            tf = obj.Tokens == token;
            obj.Tokens(tf) = [];
            obj.Pointers(tf) = [];
            
            % Update pointer
            obj.PointerListener.Enabled = false;
            obj.Figure.Pointer = obj.Pointers{end};
            obj.PointerListener.Enabled = true;
            
        end % doUnsetPointer
        
    end % private methods
    
    methods
        
        function onPointerChanged( obj, ~, ~ )
            %onPointerChanged  Event handler
            
            % Log as unknown setter
            obj.doSetPointer( 0, obj.Figure.Pointer )
            
        end % onPointerChanged
        
    end % event handlers
    
    methods( Static )
        
        function token = setPointer( figure, pointer )
            %setPointer  Set pointer
            %
            %  t = uix.PointerManager.setPointer(f,p) sets the pointer of
            %  the figure f to p.  The returned token t can be used
            %  subsequently to unset the pointer.
            
            % Get pointer manager
            obj = uix.PointerManager.getInstance( figure );
            
            % Retrieve token
            token = obj.NextToken;
            
            % Set
            obj.doSetPointer( token, pointer )
            
            % Increment token
            obj.NextToken = token + 1;
            
        end % setPointer
        
        function unsetPointer( figure, token )
            %unsetPointer  Unset pointer
            %
            %  uix.PointerManager.unsetPointer(f,t) unsets the pointer of
            %  the figure f using the token t.
            
            % Check ID
            validateattributes( token, {'numeric'}, {'scalar','integer','>',0} )
            
            % Get pointer manager
            obj = uix.PointerManager.getInstance( figure );
            
            % Unset
            obj.doUnsetPointer( token )
            
        end % unsetPointer
        
        function obj = getInstance( figure )
            %getInstance  Get pointer manager
            %
            %  m = uix.PointerManager.getInstance(f) gets the pointer
            %  manager for the figure f.
            
            % Get pointer manager
            name = 'UIxPointerManager';
            if isprop( figure, name ) % existing, retrieve
                obj = figure.( name );
            else % new, create and store
                obj = uix.PointerManager( figure );
                p = addprop( figure, name );
                p.Hidden = true;
                figure.( name ) = obj;
            end
            
        end % getInstance
        
    end % static methods
    
end % classdef