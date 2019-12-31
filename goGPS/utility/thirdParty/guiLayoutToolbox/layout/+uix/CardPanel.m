classdef CardPanel < uix.Container & uix.mixin.Panel
    %uix.CardPanel  Card panel
    %
    %  b = uix.CardPanel(p1,v1,p2,v2,...) constructs a card panel and sets
    %  parameter p1 to value v1, etc.
    %
    %  A card panel is a standard container (uicontainer) that shows one
    %  its contents and hides the others.
    %
    %  See also: uix.Panel, uix.BoxPanel, uix.TabPanel, uicontainer
    
    %  Copyright 2009-2016 The MathWorks, Inc.
    %  $Revision: 1480 $ $Date: 2017-02-15 16:56:13 +0100 (Wed, 15 Feb 2017) $
    
    methods
        
        function obj = CardPanel( varargin )
            %uix.CardPanel  Card panel constructor
            %
            %  p = uix.CardPanel() constructs a card panel.
            %
            %  p = uix.CardPanel(p1,v1,p2,v2,...) sets parameter p1 to
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
            %redraw  Redraw
            
            % Compute positions
            bounds = hgconvertunits( ancestor( obj, 'figure' ), ...
                [0 0 1 1], 'normalized', 'pixels', obj );
            padding = obj.Padding_;
            xSizes = uix.calcPixelSizes( bounds(3), -1, 1, padding, 0 );
            ySizes = uix.calcPixelSizes( bounds(4), -1, 1, padding, 0 );
            position = [padding+1 padding+1 xSizes ySizes];
            
            % Redraw contents
            selection = obj.Selection_;
            if selection ~= 0
                uix.setPosition( obj.Contents_(selection), position, 'pixels' )
            end
            
        end % redraw
        
    end % template methods
    
end % classdef