classdef ( Hidden, Sealed ) ChildObserver < handle
    %uix.ChildObserver  Child observer
    %
    %  co = uix.ChildObserver(o) creates a child observer for the graphics
    %  object o.  A child observer raises events when objects are added to
    %  and removed from the property Children of o.
    %
    %  See also: uix.Node
    
    %  Copyright 2009-2016 The MathWorks, Inc.
    %  $Revision: 1763 $ $Date: 2018-10-25 15:56:50 +0100 (Thu, 25 Oct 2018) $
    
    properties( Access = private )
        Root % root node
    end
    
    events( NotifyAccess = private )
        ChildAdded % child added
        ChildRemoved % child removed
    end
    
    methods
        
        function obj = ChildObserver( oRoot )
            %uix.ChildObserver  Child observer
            %
            %  co = uix.ChildObserver(o) creates a child observer for the
            %  graphics object o.  A child observer raises events when
            %  objects are added to and removed from the property Children
            %  of o.
            
            % Check
            assert( ispositionable( oRoot ) && ...
                isequal( size( oRoot ), [1 1] ), 'uix.InvalidArgument', ...
                'Object must be a graphics object.' )
            
            % Create root node
            nRoot = uix.Node( oRoot );
            childAddedListener = event.listener( oRoot, ...
                'ObjectChildAdded', ...
                @(~,e)obj.addChild(nRoot,e.Child) );
            childAddedListener.Recursive = true;
            nRoot.addprop( 'ChildAddedListener' );
            nRoot.ChildAddedListener = childAddedListener;
            childRemovedListener = event.listener( oRoot, ...
                'ObjectChildRemoved', ...
                @(~,e)obj.removeChild(nRoot,e.Child) );
            childRemovedListener.Recursive = true;
            nRoot.addprop( 'ChildRemovedListener' );
            nRoot.ChildRemovedListener = childRemovedListener;
            
            % Add children
            oChildren = hgGetTrueChildren( oRoot );
            for ii = 1:numel( oChildren )
                obj.addChild( nRoot, oChildren(ii) )
            end
            
            % Store properties
            obj.Root = nRoot;
            
        end % constructor
        
    end % structors
    
    methods( Access = private )
        
        function addChild( obj, nParent, oChild )
            %addChild  Add child object to parent node
            %
            %  co.addChild(np,oc) adds the child object oc to the parent
            %  node np, either as part of construction of the child
            %  observer co, or in response to an ObjectChildAdded event on
            %  an object of interest to co.  This may lead to ChildAdded
            %  events being raised on co.
            
            % Create child node
            nChild = uix.Node( oChild );
            nParent.addChild( nChild )
            positionable = ispositionable( oChild );
            if positionable == true
                % Add Internal PreSet property listener
                internalPreSetListener = event.proplistener( oChild, ...
                    findprop( oChild, 'Internal' ), 'PreSet', ...
                    @(~,~)obj.preSetInternal(nChild) );
                nChild.addprop( 'InternalPreSetListener' );
                nChild.InternalPreSetListener = internalPreSetListener;
                % Add Internal PostSet property listener
                internalPostSetListener = event.proplistener( oChild, ...
                    findprop( oChild, 'Internal' ), 'PostSet', ...
                    @(~,~)obj.postSetInternal(nChild) );
                nChild.addprop( 'InternalPostSetListener' );
                nChild.InternalPostSetListener = internalPostSetListener;
            else
                % Add ObjectChildAdded listener
                childAddedListener = event.listener( oChild, ...
                    'ObjectChildAdded', ...
                    @(~,e)obj.addChild(nChild,e.Child) );
                nChild.addprop( 'ChildAddedListener' );
                nChild.ChildAddedListener = childAddedListener;
                % Add ObjectChildRemoved listener
                childRemovedListener = event.listener( oChild, ...
                    'ObjectChildRemoved', ...
                    @(~,e)obj.removeChild(nChild,e.Child) );
                nChild.addprop( 'ChildRemovedListener' );
                nChild.ChildRemovedListener = childRemovedListener;
            end
            
            % Raise ChildAdded event
            if positionable == true && oChild.Internal == false
                notify( obj, 'ChildAdded', uix.ChildEvent( oChild ) )
            end
            
            % Add grandchildren
            if positionable == false && isblacklisted( oChild ) == false
                oGrandchildren = hgGetTrueChildren( oChild );
                for ii = 1:numel( oGrandchildren )
                    obj.addChild( nChild, oGrandchildren(ii) )
                end
            end
            
        end % addChild
        
        function removeChild( obj, nParent, oChild )
            %removeChild  Remove child object from parent node
            %
            %  co.removeChild(np,oc) removes the child object oc from the
            %  parent node np, in response to an ObjectChildRemoved event
            %  on an object of interest to co.  This may lead to
            %  ChildRemoved events being raised on co.
            
            % Get child node
            nChildren = nParent.Children;
            tf = oChild == [nChildren.Object];
            nChild = nChildren(tf);
            
            % Raise ChildRemoved event(s)
            notifyChildRemoved( nChild )
            
            % Delete child node
            delete( nChild )
            
            function notifyChildRemoved( nc )
                
                % Process child nodes
                ngc = nc.Children;
                for ii = 1:numel( ngc )
                    notifyChildRemoved( ngc(ii) )
                end
                
                % Process this node
                oc = nc.Object;
                if ispositionable( oc ) == true && oc.Internal == false
                    notify( obj, 'ChildRemoved', uix.ChildEvent( oc ) )
                end
                
            end % notifyChildRemoved
            
        end % removeChild
        
        function preSetInternal( ~, nChild )
            %preSetInternal  Perform property PreSet tasks
            %
            %  co.preSetInternal(n) caches the previous value of the
            %  property Internal of the object referenced by the node n, to
            %  enable PostSet tasks to identify whether the value changed.
            %  This is necessary since Internal AbortSet is false.
            
            oldInternal = nChild.Object.Internal;
            nChild.addprop( 'OldInternal' );
            nChild.OldInternal = oldInternal;
            
        end % preSetInternal
        
        function postSetInternal( obj, nChild )
            %postSetInternal  Perform property PostSet tasks
            %
            %  co.postSetInternal(n) raises a ChildAdded or ChildRemoved
            %  event on the child observer co in response to a change of
            %  the value of the property Internal of the object referenced
            %  by the node n.
            
            % Retrieve old and new values
            oChild = nChild.Object;
            newInternal = oChild.Internal;
            oldInternal = nChild.OldInternal;
            
            % Clean up node
            delete( findprop( nChild, 'OldInternal' ) )
            
            % Raise event
            switch newInternal
                case oldInternal % no change
                    % no event
                case true % false to true
                    notify( obj, 'ChildRemoved', uix.ChildEvent( oChild ) )
                case false % true to false
                    notify( obj, 'ChildAdded', uix.ChildEvent( oChild ) )
            end
            
        end % postSetInternal
        
    end % event handlers
    
end % classdef

function tf = ispositionable( o )
%ispositionable  True for positionable graphics

p = findprop( o, 'Position' );
tf = isgraphics( o ) && ~isempty( p ) && ...
    isequal( p.GetAccess, 'public' ) && ...
    isequal( p.SetAccess, 'public' ) && ...
    isequal( size( o.Position ), [1 4] );

end % ispositionable

function tf = isblacklisted( o )
%isblacklisted  True for objects that never have positionable graphics

tf = isa( o, 'matlab.ui.container.Menu' ) || ...
    isa( o, 'matlab.ui.container.Toolbar' ) || ...
    isa( o, 'matlab.graphics.shape.internal.AnnotationPane' );

end % isblacklisted