classdef HButtonBox < uix.HButtonBox
    %uiextras.HButtonBox  Arrange buttons horizontally in a single row
    %
    %   obj = uiextras.HButtonBox() is a type of HBox specialised for
    %   arranging a row of buttons, check-boxes or similar graphical
    %   elements. All buttons are given equal size and by default are
    %   centered in the drawing area. The justification can be changed as
    %   required.
    %
    %   obj = uiextras.HButtonBox(param,value,...) also sets one or more
    %   parameter values.
    %
    %   See the <a href="matlab:doc uiextras.HButtonBox">documentation</a> for more detail and the list of properties.
    %
    %   Examples:
    %   >> f = figure();
    %   >> b = uiextras.HButtonBox( 'Parent', f );
    %   >> uicontrol( 'Parent', b, 'String', 'One' );
    %   >> uicontrol( 'Parent', b, 'String', 'Two' );
    %   >> uicontrol( 'Parent', b, 'String', 'Three' );
    %   >> set( b, 'ButtonSize', [130 35], 'Spacing', 5 );
    %
    %   See also: uiextras.VButtonBox
    %             uiextras.HBox
    
    %  Copyright 2009-2014 The MathWorks, Inc.
    %  $Revision: 1077 $ $Date: 2015-03-19 16:44:14 +0000 (Thu, 19 Mar 2015) $
    
    methods
        
        function obj = HButtonBox( varargin )
            %uiextras.HButtonBox  Create a new horizontal button box
            
            % Call uix constructor
            obj@uix.HButtonBox( varargin{:} )
            
            % Auto-parent
            if ~ismember( 'Parent', varargin(1:2:end) )
                obj.Parent = gcf();
            end
            
        end % constructor
        
    end % structors
    
end % classdef