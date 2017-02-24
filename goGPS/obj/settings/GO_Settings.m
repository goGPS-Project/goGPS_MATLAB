%   CLASS GO_Settings
% =========================================================================
%
% DESCRIPTION
%   Collector of settings to manage a useful parameters of goGPS
%   This singleton class collects multiple objects containing various
%   parameters
%
% EXAMPLE
%   settings = GO_Settings.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc GO_Settings

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

classdef GO_Settings < Settings_Interface
    
    properties (Constant)
        V_LIGHT = 299792458;                % Velocity of light in the void [m/s]
        
        % Values as defined by standards
        
        PI_ORBIT = 3.1415926535898;         % pi as from standards
        CIRCLE_RAD = 6.2831853071796;       % Circle as from standards
        
        % Standard atmosphere parameters (struct PRES, STD_TEMP. STD_HUMI) - Berg, 1948
        ATM = struct('PRES', 1013.25, ...      % pressure [mbar]
                     'STD_TEMP', 291.15, ...   % temperature [K] (18? C)
                     'STD_HUMI', 50.0);        % humidity [%]               
    end
    
    properties % Public Access
        cur_settings = Main_Settings();        % Processing settings
    end
    
    % =========================================================================
    %  INIT
    % =========================================================================
    methods (Access = private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function obj = GO_Settings()
        end
    end
    
    % =========================================================================
    %  SINGLETON GETTERS
    % =========================================================================
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function obj = getInstance()
            persistent unique_instance_settings__
            if isempty(unique_instance_settings__)
                obj = GO_Settings();
                unique_instance_settings__ = obj;
            else
                obj = unique_instance_settings__;
            end
        end
        
        function cur_settings = getCurrentSettings()
            this = GO_Settings.getInstance();
            % Return the handler to the object containing the current settings
            cur_settings = handle(this.cur_settings);
        end
        
    end
    
    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods % Public Access
        function import(obj, settings)
            % This function try to import settings from another setting object
            if isprop(settings, 'ps')
                obj.cur_settings.import(settings.cur_settings);
            else
                try
                    obj.cur_settings.import(settings);
                catch ex
                    obj.logger.addWarning(['GO_Settings.import failed to import settings (invalid input settings) ', ex.message()]);
                end
            end            
        end
        
        function str = toString(obj, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end            
            str = [str '---- CONSTANTS -----------------------------------------------------------' 10 10];
            str = [str sprintf(' VLIGHT:                                           %g\n', obj.V_LIGHT)];
            str = [str sprintf(' PI_ORBIT:                                         %g\n', obj.PI_ORBIT)];
            str = [str sprintf(' CIRCLE_RAD:                                       %g\n', obj.CIRCLE_RAD)];
            str = [str sprintf(' STANDARD ATMOSPHERE (Berg, 1948):\n  - PRESSURE [mBar]                                %g\n  - TEMPERATURE [K]                                %g\n  - HUMIDITY [%%]                                   %g\n\n', obj.ATM.PRES, obj.ATM.STD_TEMP, obj.ATM.STD_HUMI)];
            str = obj.cur_settings.toString(str);
        end
        
        function str_cell = export(obj, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the obj            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = obj.cur_settings.export(str_cell);
        end        
    end
    
    % =========================================================================
    %  TEST
    % =========================================================================
    methods (Static, Access = 'public')
        function test()      
            % test the class
            s = GO_Settings.getInstance();
            s.testInterfaceRoutines();
        end
    end    
end
