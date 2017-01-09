%   CLASS PP_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the Pre Processing parameters
%
% EXAMPLE
%   settings = PP_Settings();
%
% FOR A LIST OF CONSTANTs and METHODS use doc PP_Settings

%----------------------------------------------------------------------------------------------
%                           goGPS v0.9.1
% Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
% Written by:       Gatti Andrea
% Contributors:     Gatti Andrea, ...
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------
classdef PP_Settings < Settings_Interface
    
    properties (Constant)
    end
    
    properties
    end
    
    methods
        function obj = PP_Settings()            
        end
    end
    
    methods 
        function import(obj, settings)
            % This function import PP (only) settings from another setting object
        end
        
        function str = toString(obj, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end
            
            str = [str '---- POST PROCESSING -------------------------------------------' 10 10];
            str = [str ' nothing to report ' 10 10];
        end
        
        function str_cell = export(obj, str_cell)            
            % Conversion to string ini format of the minimal information needed to reconstruct the obj            
            if (nargin == 1)
                str_cell = {};
            end
        end
    end        
end
