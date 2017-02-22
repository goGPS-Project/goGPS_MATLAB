%   CLASS Settings_Interface
% =========================================================================
%
% DESCRIPTION
%   Abstract class with basic properties and methods that a setting class
%   must have
%
% COMMENTS
% Settings have been build with multiple inheritance
% A standard abstract interface have been created: Setting_Interface
% it add to each subclass the object logger
% force the subclasses to implement three basic methods:
%  - import:       when a setting object of identical class, or that inherits from the same class is passed to this function, the relative parameters are copied in the calling object
%                  when an ini file is the input of the function, the object is updated with the settings contained into the file
%  - toString:     display the content of the object, in a human readable way, a goGPS user can "ask" for the value of the settings on screen
%  - export:       create a cell array of strings containing the settings in plain text ini format. The variable it's the raw data format of Ini_Manager
% 

%--------------------------------------------------------------------------
%               ___ ___ ___ 
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.9.1
% 
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011 
%--------------------------------------------------------------------------

classdef Settings_Interface < handle
    properties (SetAccess = protected, GetAccess = protected)
        logger = Logger.getInstance(); % Handler to the logger object
    end
    
    properties (Abstract)
    end
        
    methods  (Abstract)        
        import(obj, settings)
        % This function import the settings of the current class from another setting object having the same properties, or an ini file
        
        str = toString(obj, str)
        % Display the content of the class in a human readable format

        str_cell = export(obj, str_cell)
        % Conversion to string ini format of the minimal information needed to reconstruct the obj
    end
    
    methods (Access = 'protected')
        function testInterfaceRoutines(this)
            % test the class (Interface Routines)
            
            try
                vl = this.logger.getVerbosityLev();
                this.logger.setVerbosityLev(1e3);
                test = this;
                raw_data = test.export();
                
                ini = Ini_Manager('test.ini', raw_data);
                ini.showData();
                
                test.import(ini);
                test_copy = repmat(test,2,1); test_copy = test_copy(2);
                test.import(test_copy);
                clear test_copy;
                fprintf('\n');
                disp(test.toString());
            catch ex
                this.logger.addError(['Test failed: ' ex.message]);
            end
            this.logger.setVerbosityLev(vl);
        end
    end
end
