classdef Flex < handle
    %uix.mixin.Flex  Flex mixin
    %
    %  uix.mixin.Flex is a mixin class used by flex containers to provide
    %  various properties and helper methods.
    
    %  Copyright 2016 The MathWorks, Inc.
    %  $Revision: 1435 $ $Date: 2016-11-17 17:50:34 +0000 (Thu, 17 Nov 2016) $
    
    properties( GetAccess = protected, SetAccess = private )
        Pointer = 'unset' % mouse pointer
    end
    
    properties( Access = private )
        Figure = gobjects( 0 ); % mouse pointer figure
        Token = -1 % mouse pointer token
    end
    
    methods
        
        function delete( obj )
            %delete  Destructor
            
            % Clean up
            if ~strcmp( obj.Pointer, 'unset' )
                obj.unsetPointer()
            end
            
        end % destructor
        
    end % structors
    
    methods( Access = protected )
        
        function setPointer( obj, figure, pointer )
            %setPointer  Set pointer
            %
            %  c.setPointer(f,p) sets the pointer for the figure f to p.
            
            % If set, unset
            if obj.Token ~= -1
                obj.unsetPointer()
            end
            
            % Set
            obj.Token = uix.PointerManager.setPointer( figure, pointer );
            obj.Figure = figure;
            obj.Pointer = pointer;
            
        end % setPointer
        
        function unsetPointer( obj )
            %unsetPointer  Unset pointer
            %
            %  c.unsetPointer() undoes the previous pointer set.
            
            % Check
            assert( obj.Token ~= -1, 'uix:InvalidOperation', ...
                'Pointer is already unset.' )
            
            % Unset
            uix.PointerManager.unsetPointer( obj.Figure, obj.Token );
            obj.Figure = gobjects( 0 );
            obj.Pointer = 'unset';
            obj.Token = -1;
            
        end % unsetPointer
        
    end % helper methods
    
end % classdef