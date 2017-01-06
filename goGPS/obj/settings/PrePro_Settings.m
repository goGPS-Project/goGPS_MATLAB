%   CLASS PrePro_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the Pre Processing parameters
%
% EXAMPLE
%   settings = PrePro_Settings();
%
% FOR A LIST OF CONSTANTs and METHODS use doc PrePro_Settings

%----------------------------------------------------------------------------------------------
%                           goGPS v0.5.9
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
classdef PrePro_Settings < Settings_Interface
    
    properties (Constant)
    end
    
    properties

        % PRE PROCESSING --------------------------------------------------
        
        % Cycle slip threshold (pre-processing) [cycles]
        cs_thr_pre_pro = 1;
    end
    
    methods
        function obj = PrePro_Settings()            
        end
    end
    
    methods 
        function copyFrom(obj, settings)
            % This function import Post Processing (only) settings from another setting object
            obj.cs_thr_pre_pro = settings.cs_thr_pre_pro;
        end
        
        function str = toString(obj, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end
            
            str = [str '---- PRE PROCESSING ------------------------------------------------------' 10 10];
            str = [str sprintf(' Cycle slip threshold [cycles]            %g\n\n', obj.cs_thr_pre_pro)];
        end
        
        function str_cell = toIniString(obj, str_cell)            
            % Conversion to string of the minimal information needed to reconstruct the obj            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniStringSection('PRE_PROCESSING', str_cell);
            str_cell = Ini_Manager.toIniString('cs_thr_pre_pro', obj.cs_thr_pre_pro, str_cell);
        end
    end        
end
