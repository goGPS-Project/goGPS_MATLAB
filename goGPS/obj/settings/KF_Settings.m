%   CLASS KF_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the Kalman Filter Settings
%
% EXAMPLE
%   settings = KF_Settings();
%
% FOR A LIST OF CONSTANTs and METHODS use doc KF_Settings

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
classdef KF_Settings < Settings_Interface
    
    properties (Constant)
    end
    
    properties (SetAccess = protected, GetAccess = public)
        
        %------------------------------------------------------------------
        % RECEIVER POSITION / MOTION 
        %------------------------------------------------------------------
        
        % Std of initial state [m]
        sigma0_pos = 1; 
        
        % Std of velocity ENU coordinates [m/s]
        sigma_vel_ENU = struct('E', 0.5, 'N', 0.5, 'U', 0.1);
        % Std of 3D velocity modulus [m/s]
        sigma_vel_mod = 0.1;
                        
        %------------------------------------------------------------------
        % ATHMOSPHERE 
        %------------------------------------------------------------------

        % Std of apriori tropospheric delay
        sigma0_tropo = 0.1;
        % Std of tropospheric delay
        sigma_tropo = 4.6e-4;
        
        %------------------------------------------------------------------
        % KF
        %------------------------------------------------------------------
        
        % Minimum number of satellites to be used in the Kalman filter
        kf_min_n_sat = 2;
        
        % Order of the dynamic model polynomial
        kf_order = 1;

    end
    
    methods
        function obj = KF_Settings()
        end
    end
    
    methods 
        function copyFrom(obj, kf_settings)
            % This function import KF (only) settings from another setting object
            obj.sigma0_pos    = kf_settings.sigma0_pos;
            obj.sigma_vel_ENU = kf_settings.sigma_vel_ENU;
            obj.sigma_vel_mod = kf_settings.sigma_vel_mod;
            obj.sigma0_tropo  = kf_settings.sigma0_tropo;
            obj.sigma_tropo   = kf_settings.sigma_tropo;
            obj.kf_min_n_sat  = kf_settings.kf_min_n_sat;
            obj.kf_order      = kf_settings.kf_order;
        end
        
        function str = toString(obj, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end            
            str = [str '---- KALMAN FILTER PARAMETERS --------------------------------------------' 10 10];
            str = [str sprintf(' STD of initial state [m]:                         %g\n', obj.sigma0_pos)];
            str = [str sprintf(' STD of ENU velocity [m]:                          %g %g %g\n', struct2array(obj.sigma_vel_ENU))];
            str = [str sprintf(' STD of 3D velocity modulus [m]:                   %g\n\n', obj.sigma_vel_mod)];
            str = [str sprintf(' STD of apriori tropospheric delay:                %g\n', obj.sigma0_tropo)];
            str = [str sprintf(' STD of tropospheric delay:                        %g\n\n', obj.sigma_tropo)];
            str = [str sprintf(' Minimum number of satellite per epoch:            %d\n', obj.kf_min_n_sat)];
            str = [str sprintf(' Oreder of the KF:                                 %d\n\n', obj.kf_order)];
        end
        
        function str_cell = toIniString(obj, str_cell)            
            % Conversion to string of the minimal information needed to reconstruct the obj            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniStringSection('KALMAN_FILTER', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_pos', obj.sigma0_pos, str_cell);
            str_cell = Ini_Manager.toIniString('sigma_vel_ENU', struct2array(obj.sigma_vel_ENU), str_cell);
            str_cell = Ini_Manager.toIniString('sigma_vel_mod', obj.sigma_vel_mod, str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_tropo', obj.sigma0_tropo, str_cell);
            str_cell = Ini_Manager.toIniString('sigma_tropo', obj.sigma_tropo, str_cell);
            str_cell = Ini_Manager.toIniString('kf_min_n_sat', obj.kf_min_n_sat, str_cell);
            str_cell = Ini_Manager.toIniString('kf_order', obj.kf_order, str_cell);
        end
        
    end        
end
