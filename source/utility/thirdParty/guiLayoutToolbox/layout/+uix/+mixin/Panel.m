classdef Panel < uix.mixin.Container
    %uix.mixin.Panel  Panel mixin
    %
    %  uix.mixin.Panel is a mixin class used by panels to provide various
    %  properties and template methods.
    
    %  Copyright 2009-2015 The MathWorks, Inc.
    %  $Revision: 1435 $ $Date: 2016-11-17 17:50:34 +0000 (Thu, 17 Nov 2016) $
    
    properties( Access = public, Dependent, AbortSet )
        Selection % selected contents
    end
    
    properties( Access = protected )
        Selection_ = 0 % backing for Selection
    end
    
    properties( Access = protected )
        G1218142 = false % bug flag
    end
    
    events( NotifyAccess = protected )
        SelectionChanged % selection changed
    end
    
    methods
        
        function value = get.Selection( obj )
            
            value = obj.Selection_;
            
        end % get.Selection
        
        function set.Selection( obj, value )
            
            % Check
            assert( isa( value, 'double' ), 'uix:InvalidPropertyValue', ...
                'Property ''Selection'' must be of type double.' )
            assert( isequal( size( value ), [1 1] ), ...
                'uix:InvalidPropertyValue', ...
                'Property ''Selection'' must be scalar.' )
            assert( isreal( value ) && rem( value, 1 ) == 0, ...
                'uix:InvalidPropertyValue', ...
                'Property ''Selection'' must be an integer.' )
            n = numel( obj.Contents_ );
            if n == 0
                assert( value == 0, 'uix:InvalidPropertyValue', ...
                    'Property ''Selection'' must be 0 for a container with no children.' )
            else
                assert( value >= 1 && value <= n, 'uix:InvalidPropertyValue', ...
                    'Property ''Selection'' must be between 1 and the number of children.' )
            end
            
            % Set
            oldSelection = obj.Selection_;
            newSelection = value;
            obj.Selection_ = newSelection;
            
            % Show selected child
            obj.showSelection()
            
            % Mark as dirty
            obj.Dirty = true;
            
            % Raise event
            notify( obj, 'SelectionChanged', ...
                uix.SelectionData( oldSelection, newSelection ) )
            
        end % set.Selection
        
    end % accessors
    
    methods( Access = protected )
        
        function addChild( obj, child )
            
            % Check for bug
            if verLessThan( 'MATLAB', '8.5' ) && strcmp( child.Visible, 'off' )
                obj.G1218142 = true;
            end
            
            % Select new content
            oldSelection = obj.Selection_;
            newSelection = numel( obj.Contents_ ) + 1;
            obj.Selection_ = newSelection;
            
            % Call superclass method
            addChild@uix.mixin.Container( obj, child )
            
            % Show selected child
            obj.showSelection()
            
            % Notify selection change
            obj.notify( 'SelectionChanged', ...
                uix.SelectionData( oldSelection, newSelection ) )
            
        end % addChild
        
        function removeChild( obj, child )
            
            % Adjust selection if required
            contents = obj.Contents_;
            index = find( contents == child );
            oldSelection = obj.Selection_;
            if index < oldSelection
                newSelection = oldSelection - 1;
            elseif index == oldSelection
                newSelection = min( oldSelection, numel( contents ) - 1 );
            else % index > oldSelection
                newSelection = oldSelection;
            end
            obj.Selection_ = newSelection;
            
            % Call superclass method
            removeChild@uix.mixin.Container( obj, child )
            
            % Show selected child
            obj.showSelection()
            
            % Notify selection change
            if oldSelection ~= newSelection
                obj.notify( 'SelectionChanged', ...
                    uix.SelectionData( oldSelection, newSelection ) )
            end
            
        end % removeChild
        
        function reorder( obj, indices )
            %reorder  Reorder contents
            %
            %  c.reorder(i) reorders the container contents using indices
            %  i, c.Contents = c.Contents(i).
            
            % Reorder
            selection = obj.Selection_;
            if selection ~= 0
                obj.Selection_ = find( indices == selection );
            end
            
            % Call superclass method
            reorder@uix.mixin.Container( obj, indices )
            
        end % reorder
        
        function showSelection( obj )
            %showSelection  Show selected child, hide the others
            %
            %  c.showSelection() shows the selected child of the container
            %  c, and hides the others.
            
            % Set positions and visibility
            selection = obj.Selection_;
            children = obj.Contents_;
            for ii = 1:numel( children )
                child = children(ii);
                if ii == selection
                    if obj.G1218142
                        warning( 'uix:G1218142', ...
                            'Selected child of %s is not visible due to bug G1218142.  The child will become visible at the next redraw.', ...
                            class( obj ) )
                        obj.G1218142 = false;
                    else
                        child.Visible = 'on';
                    end
                    if isa( child, 'matlab.graphics.axis.Axes' )
                        child.ContentsVisible = 'on';
                    end
                else
                    child.Visible = 'off';
                    if isa( child, 'matlab.graphics.axis.Axes' )
                        child.ContentsVisible = 'off';
                    end
                    % As a remedy for g1100294, move off-screen too
                    margin = 1000;
                    if isa( child, 'matlab.graphics.axis.Axes' ) ...
                            && strcmp(child.ActivePositionProperty, 'outerposition' )
                        child.OuterPosition(1) = -child.OuterPosition(3)-margin;
                    else
                        child.Position(1) = -child.Position(3)-margin;
                    end
                end
            end
            
        end % showSelection
        
    end % template methods
    
end % classdef