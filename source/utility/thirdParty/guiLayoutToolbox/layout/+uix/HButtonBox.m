classdef HButtonBox < uix.ButtonBox
    %uix.HButtonBox  Horizontal button box
    %
    %  b = uix.HButtonBox(p1,v1,p2,v2,...) constructs a horizontal button
    %  box and sets parameter p1 to value v1, etc.
    %
    %  A horizontal button box lays out equally sized buttons from left to
    %  right.
    %
    %  See also: uix.VButtonBox
    
    %  Copyright 2009-2016 The MathWorks, Inc.
    %  $Revision: 1480 $ $Date: 2017-02-15 16:56:13 +0100 (Wed, 15 Feb 2017) $
    
    methods
        
        function obj = HButtonBox( varargin )
            %uix.HButtonBox  Horizontal button box constructor
            %
            %  b = uix.HButtonBox() constructs a horizontal button box.
            %
            %  b = uix.HButtonBox(p1,v1,p2,v2,...) sets parameter p1 to
            %  value v1, etc.
            
            % Set properties
            try
                uix.set( obj, varargin{:} )
            catch e
                delete( obj )
                e.throwAsCaller()
            end
            
        end % constructor
        
    end % structors
    
    methods( Access = protected )
        
        function redraw( obj )
            
            % Compute positions
            bounds = hgconvertunits( ancestor( obj, 'figure' ), ...
                [0 0 1 1], 'normalized', 'pixels', obj );
            buttonSize = obj.ButtonSize_;
            padding = obj.Padding_;
            spacing = obj.Spacing_;
            c = numel( obj.Contents_ );
            if 2 * padding + (c-1) * spacing + c * buttonSize(1) > bounds(3)
                xSizes = uix.calcPixelSizes( bounds(3), -ones( [c 1] ), ...
                    ones( [c 1] ), padding, spacing ); % shrink to fit
            else
                xSizes = repmat( buttonSize(1), [c 1] );
            end
            switch obj.HorizontalAlignment
                case 'left'
                    xPositions = [cumsum( [0; xSizes(1:c-1,:)] ) + ...
                        padding + spacing * transpose( 0:c-1 ) + 1, xSizes];
                case 'center'
                    xPositions = [cumsum( [0; xSizes(1:c-1,:)] ) + ...
                        spacing * transpose( 0:c-1 ) + bounds(3) / 2 - ...
                        sum( xSizes ) / 2 - spacing * (c-1) / 2 + 1, ...
                        xSizes];
                case 'right'
                    xPositions = [cumsum( [0; xSizes(1:c-1,:)] ) + ...
                        spacing * transpose( 0:c-1 ) + bounds(3) - ...
                        sum( xSizes ) - spacing * (c-1) - padding + 1, ...
                        xSizes];
            end
            if 2 * padding + buttonSize(2) > bounds(4)
                ySizes = repmat( uix.calcPixelSizes( bounds(4), -1, 1, ...
                    padding, spacing ), [c 1] ); % shrink to fit
            else
                ySizes = repmat( buttonSize(2), [c 1] );
            end
            switch obj.VerticalAlignment
                case 'top'
                    yPositions = [bounds(4) - ySizes - padding + 1, ySizes];
                case 'middle'
                    yPositions = [(bounds(4) - ySizes) / 2 + 1, ySizes];
                case 'bottom'
                    yPositions = [repmat( padding, [c 1] ) + 1, ySizes];
            end
            positions = [xPositions(:,1), yPositions(:,1), ...
                xPositions(:,2), yPositions(:,2)];
            
            % Set positions
            children = obj.Contents_;
            for ii = 1:numel( children )
                uix.setPosition( children(ii), positions(ii,:), 'pixels' )
            end
            
        end % redraw
        
    end % template methods
    
end % classdef