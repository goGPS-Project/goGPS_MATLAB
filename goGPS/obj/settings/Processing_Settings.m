%   CLASS Processing_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the Processing parameters
%
% EXAMPLE
%   settings = Processing_Settings();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Processing_Settings

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
classdef Processing_Settings < Settings_Interface & IO_Settings & PrePro_Settings & KF_Settings
    
    properties (Constant)
    end
    
    %  Processing parameters
    % ------------------------------
    % note that some of the processing parameters are actually properties of other objects 
    % e.g. sigma_vel_ENU depend on the receiver motion
    % However for various reasons I can decide to ignore the "correct"
    % values and redefine them. So the values are also
    % stored here as parameters used for a specific processing
    properties

        %------------------------------------------------------------------
        % RECEIVER 
        %------------------------------------------------------------------
                
        % Std of phase observations [m]
        sigma_ph = 0.003; % (maximize to obtain a code-based solution)
        % Std of iono-free phase observations [m]
        sigmaq_ph_if = 0.009;
        
        % Std of a priori ambiguity combinations [cycles]
        sigma0_N = 1000;
        
        % Std of apriori receiver clock
        sigma0_clock = 4.5e-09        
        % Std of receiver clock
        sigma0_r_clock = 1e3;
        
        % Signal-to-noise ratio threshold [dB]
        snr_thr = 0;
        
        %------------------------------------------------------------------
        % PROCESSING PARAMETERS 
        %------------------------------------------------------------------
        
        % Cycle slip threshold (processing) [cycles]
        cs_thr = 1;

        % Parameter used to select the weightening mode for GPS observations
        w_mode = 1;
        %  - weights = 0: same weight for all the observations
        %  - weights = 1: weight based on satellite elevation (sin)
        %  - weights = 2: weight based on signal-to-noise ratio
        %  - weights = 3: weight based on combined elevation and signal-to-noise ratio
        %  - weights = 4: weight based on satellite elevation (exp)
        
        % Weight function parameters (when based on SNR)
        w_snr = struct('a', 30, 'zero', 10, 'one', 50, 'A', 30);
        
        % Cut-off [degrees]
        cutoff = 10;

        %------------------------------------------------------------------
        % INTEGER AMBIGUITY RESOLUTION
        %------------------------------------------------------------------

        % Ambiguity restart mode
        iar_restart_mode = 2;

        % Parameter used to select Integer Least Squares estimator
        iar_mode = 1; % ILS method with shrinking ellipsoid during search (LAMBDA3)
        % - iar_mode = 0; % ILS method with numeration in search (LAMBDA2)
        % - iar_mode = 1; % ILS method with shrinking ellipsoid during search (LAMBDA3)
        % - iar_mode = 2; % ILS method with numeration in search (LAMBDA3)
        % - iar_mode = 3; % integer rounding method (LAMBDA3)
        % - iar_mode = 4; % integer bootstrapping method (LAMBDA3)
        % - iar_mode = 5; % Partial Ambiguity Resolution (PAR) (LAMBDA3)
        
        % User defined fixed failure rate (for methods 1,2) or minimum required success rate (for method 5)
        iar_P0 = 0.001;
        
        % User defined threshold for ratio test
        iar_mu = 0.5;
        
        % Flag for enabling the automatic determination of mu
        iar_flag_auto_mu = true;
        
        % Flag for enabling the default value for P0
        iar_flag_default_P0 = true;

        %------------------------------------------------------------------
        % DEM 
        %------------------------------------------------------------------
        
        % Std of DEM height [m]
        sigma_dtm = 0.03  % (maximize to disable DEM usage: e.g. 1e30)
        
        %------------------------------------------------------------------
        % ATHMOSPHERE
        %------------------------------------------------------------------
        
        % Ionospheric model to be used (0: none, 1: Geckle and Feen, 2: Klobuchar, 3: SBAS)
        iono_model = 2; % Klobuchar model
        % - iono_model = 0: no model
        % - iono_model = 1: Geckle and Feen model
        % - iono_model = 2: Klobuchar model
        % - iono_model = 3: SBAS grid
        
        % Tropospheric model to be used (0: none, 1: Saastamoinen std parameters, 2: Saastamoinen global pararameters)
        tropo_model = 1; % Saastamoinen model (with standard atmosphere parameters)
        % - tropo_model = 0: no model
        % - tropo_model = 1: Saastamoinen model (with standard atmosphere parameters)
        % - tropo_model = 2: Saastamoinen model (with Global Pressure Temperature model)
        
        %-------------------------------------------------------------------------------
        % THRESHOLDS
        %-------------------------------------------------------------------------------
        
        % Threshold on the code point-positioning least squares estimation error [m]
        pp_spp_thr = 4;                 
        % Threshold on the maximum residual of code observations [m]
        pp_max_code_err_thr = 30;            
        % Threshold on the maximum residual of phase observations [m]
        pp_max_phase_err_thr = 0.2;  
        
    end
    
    methods
        function obj = Processing_Settings()
        end
    end
    
    methods 
        function copyFrom(obj, settings)
            % This function import Processing settings from another setting object
            obj.copyFrom@IO_Settings(settings);
            obj.copyFrom@PrePro_Settings(settings);
            obj.copyFrom@KF_Settings(settings);
            obj.sigma0_pos = kf_settings.sigma0_pos;
        end
                
        function str = toString(obj, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end
            
            str = obj.toString@IO_Settings(str);
            str = obj.toString@PrePro_Settings(str);
            str = obj.toString@KF_Settings(str);
        end
        
        function str_cell = toIniString(obj, str_cell)            
            % Conversion to string of the minimal information needed to reconstruct the obj            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = obj.toIniString@IO_Settings(str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = obj.toIniString@PrePro_Settings(str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = obj.toIniString@KF_Settings(str_cell);
        end

    end        
end
