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
classdef GO_Settings < Settings_Interface
    
    properties (Constant)
        V_LIGHT = 299792458;                % Velocity of light in the void [m/s]
        
        % Values as defined by standards
        
        PI_ORBIT = 3.1415926535898;         % pi as from standards
        CIRCLE_RAD = 6.2831853071796;       % Circle as from standards
        
        % Standard atmosphere parameters (struct PRES, STD_TEMP. STD_HUMI) - Berg, 1948
        ATM = struct('PRES', 1013.25, ...      % pressure [mbar]
                     'STD_TEMP', 291.15, ...   % temperature [K]
                     'STD_HUMI', 50.0);        % humidity [%]               
    end
    
    properties % Public Access
        
        processing = Processing_Settings();        
    end
    
    methods (Access = private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function obj = GO_settings()
        end
    end
    
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function obj = getInstance()
            persistent uniqueInstance
            if isempty(uniqueInstance)
                obj = GO_Settings();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
    
    %*** Define your own methods for SingletonImpl.
    methods % Public Access

            
    end
    
    methods % Public Access (Legacy support)
    end
    
end
