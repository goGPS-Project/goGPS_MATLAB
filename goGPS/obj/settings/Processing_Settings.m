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
    
    properties (Constant, Access = 'protected')
        IAR_MODE = {'0 - ILS method with numeration in search (LAMBDA2)', ...
                    '1 - ILS method with shrinking ellipsoid during search (LAMBDA3)' ...
                    '2 - ILS method with numeration in search (LAMBDA3)', ...
                    '3 - integer rounding method (LAMBDA3)', ...
                    '4 - iar_mode = 4; % integer bootstrapping method (LAMBDA3)', ...
                    '5 - Partial Ambiguity Resolution (PAR) (LAMBDA3)'}
       
        W_MODE = {'same weight for all the observations', ...
                  'weight based on satellite elevation (sin)' ...
                  'weight based on signal-to-noise ratio', ...
                  'weight based on combined elevation and signal-to-noise ratio', ...
                  'weight based on satellite elevation (exp)'}
              
        IONO_MODE = {'0 - no model', ...
                     '1 - Geckle and Feen model' ...
                     '2 - Klobuchar model', ...
                     '3 - SBAS grid'}
                 
        TROPO_MODE = {'0 - no model', ...
                      '1 - Saastamoinen model (with standard atmosphere parameters)' ...
                      '2 - Saastamoinen model (with Global Pressure Temperature model)'}
         

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
        
        % Std of apriori receiver clock
        sigma0_clock = 4.5e-09        
        % Std of receiver clock
        sigma0_r_clock = 1e3;
        
        % Signal-to-noise ratio threshold [dB]
        snr_thr = 0;

        % Cut-off [degrees]
        cutoff = 10;
                
        %------------------------------------------------------------------
        % PROCESSING PARAMETERS 
        %------------------------------------------------------------------
        
        % Parameter used to select the weightening mode for GPS observations
        w_mode = 1;
        %  - weights = 0: same weight for all the observations
        %  - weights = 1: weight based on satellite elevation (sin)
        %  - weights = 2: weight based on signal-to-noise ratio
        %  - weights = 3: weight based on combined elevation and signal-to-noise ratio
        %  - weights = 4: weight based on satellite elevation (exp)
        
        % Cycle slip threshold (processing) [cycles]
        cs_thr = 1;

        % Weight function parameters (when based on SNR)
        w_snr = struct('a', 30, 'zero', 10, 'one', 50, 'A', 30);

        %-------------------------------------------------------------------------------
        % THRESHOLDS
        %-------------------------------------------------------------------------------
        
        % Threshold on the code point-positioning least squares estimation error [m]
        pp_spp_thr = 4;                 
        % Threshold on the maximum residual of code observations [m]
        pp_max_code_err_thr = 30;            
        % Threshold on the maximum residual of phase observations [m]
        pp_max_phase_err_thr = 0.2;  
        
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

        % Std of apriori ambiguity combinations [cycles]
        sigma0_N = 1000;

        % User defined threshold for ratio test
        iar_mu = 0.5;
        
        % Flag for enabling the automatic determination of mu
        iar_flag_auto_mu = true;
        
        % Flag for enabling the default value for P0
        iar_flag_default_P0 = true;
        
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
        
        %------------------------------------------------------------------
        % DEM 
        %------------------------------------------------------------------
        
        % Std of DEM height [m]
        sigma_dtm = 0.03  % (maximize to disable DEM usage: e.g. 1e30)
        
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
            str = [str '---- RECEIVERS -----------------------------------------------------------' 10 10];
            str = [str sprintf(' STD of phase observations [m]:                    %g\n', obj.sigma_ph)];
            str = [str sprintf(' STD of iono-free phase observations [m]:          %g\n\n', obj.sigmaq_ph_if)];
            str = [str sprintf(' STD of apriori receiver clock:                    %g\n', obj.sigma0_clock)];
            str = [str sprintf(' STD of receiver clock:                            %g\n\n', obj.sigma0_r_clock)];
            str = [str sprintf(' Signal-to-noise ratio threshold [dB]:             %d\n\n', obj.snr_thr)];
            str = [str sprintf(' Cut-off [degrees]:                                %d\n\n', obj.cutoff)];
            str = obj.toString@PrePro_Settings(str);
            str = [str '---- PROCESSING PARAMETERS -----------------------------------------------' 10 10];
            str = [str sprintf(' Processing using %s\n\n', obj.W_MODE{obj.w_mode+1})];
            str = [str sprintf(' Cycle slip threshold (processing) [cycles]:       %d\n\n', obj.cs_thr)];
            str = [str sprintf(' Weight function parameters (when based on SNR): \n')];
            str = [str sprintf('   - w.a:     %d\n', obj.w_snr.a)];
            str = [str sprintf('   - w.zero:  %d\n', obj.w_snr.zero)];
            str = [str sprintf('   - w.one:   %d\n', obj.w_snr.one)];
            str = [str sprintf('   - w.A:     %d\n\n', obj.w_snr.A)];            
            str = [str '---- THRESHOLDS ----------------------------------------------------------' 10 10];
            str = [str sprintf(' Threshold on code LS estimation error [m]:        %g\n', obj.pp_spp_thr)];
            str = [str sprintf(' Threshold on maximum residual of code obs [m]:    %g\n', obj.pp_max_code_err_thr)];
            str = [str sprintf(' Threshold on maximum residual of phase obs [m]:   %g\n\n', obj.pp_max_phase_err_thr)];
            str = obj.toString@KF_Settings(str);
            str = [str '---- ABIGUITY ------------------------------------------------------------' 10 10];
            str = [str sprintf(' Ambiguity restart mode:                           %d\n\n', obj.iar_restart_mode)];
            str = [str sprintf(' Using method: %s\n\n', obj.IAR_MODE{obj.iar_mode+1})];
            str = [str sprintf(' User defined fixed failure rate (methods 1,2):    %g\n', obj.iar_P0)];
            str = [str sprintf(' User defined minimum success rate (for method 5): %g\n', obj.iar_P0)];
            str = [str sprintf(' STD of apriori ambiguity combinations [cycles]:   %d\n\n', obj.sigma0_N)];
            str = [str sprintf(' User defined threshold for ratio test:            %g\n', obj.iar_mu)];
            str = [str sprintf(' Automatic determination of mu:                    %d\n', obj.iar_flag_auto_mu)];
            str = [str sprintf(' Use default value for P0:                         %d\n\n', obj.iar_flag_default_P0)];
            str = [str '---- ATHMOSPHERE ---------------------------------------------------------' 10 10];
            str = [str sprintf(' Ionospheric model:  %s\n', obj.IONO_MODE{obj.iono_model+1})];
            str = [str sprintf(' Tropospheric model: %s\n\n', obj.TROPO_MODE{obj.tropo_model+1})];
            str = [str '---- DEM -----------------------------------------------------------------' 10 10];
            str = [str sprintf(' STD of DEM model [m]:                             %g\n', obj.iar_mu)];
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
