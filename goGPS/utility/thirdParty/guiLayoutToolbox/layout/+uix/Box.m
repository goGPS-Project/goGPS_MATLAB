classdef Box < uix.Container & uix.mixin.Container
    %uix.Box  Box and grid base class
    %
    %  uix.Box is a base class for containers with spacing between
    %  contents.
    
    %  Copyright 2009-2015 The MathWorks, Inc.
    %  $Revision: 1594 $ $Date: 2018-03-28 02:27:52 +1100 (Wed, 28 Mar 2018) $
    
    properties( Access = public, Dependent, AbortSet )
        Spacing = 0 % space between contents, in pixels
    end
    
    properties( Access = protected )
        Spacing_ = 0 % backing for Spacing
    end
    
    methods
        
        function value = get.Spacing( obj )
            
            value = obj.Spacing_;
            
        end % get.Spacing
        
        function set.Spacing( obj, value )
            
            % Check
            assert( isa( value, 'double' ) && isscalar( value ) && ...
                isreal( value ) && ~isinf( value ) && ...
                ~isnan( value ) && value >= 0, ...
                'uix:InvalidPropertyValue', ...
                'Property ''Spacing'' must be a non-negative scalar.' )
            
            % Set
            obj.Spacing_ = value;
            
            % Mark as dirty
            obj.Dirty = true;
            
        end % set.Spacing
        
    end % accessors
    
end % classdef