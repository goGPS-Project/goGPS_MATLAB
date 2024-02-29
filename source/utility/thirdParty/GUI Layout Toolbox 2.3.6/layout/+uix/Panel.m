classdef Panel < matlab.ui.container.Panel & uix.mixin.Panel
    %uix.Panel  Standard panel
    %
    %  b = uix.Panel(p1,v1,p2,v2,...) constructs a standard panel and sets
    %  parameter p1 to value v1, etc.
    %
    %  A card panel is a standard panel (uipanel) that shows one its
    %  contents and hides the others.
    %
    %  See also: uix.CardPanel, uix.BoxPanel, uipanel
    
    %  Copyright 2009-2020 The MathWorks, Inc.
    
    methods
        
        function obj = Panel( varargin )
            %uix.Panel  Standard panel constructor
            %
            %  p = uix.Panel() constructs a standard panel.
            %
            %  p = uix.Panel(p1,v1,p2,v2,...) sets parameter p1 to value
            %  v1, etc.           
            
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