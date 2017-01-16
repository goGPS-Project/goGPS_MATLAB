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
classdef Processing_Settings < Settings_Interface & IO_Settings
    
    properties (Constant, Access = 'protected')
        IAR_MODE = {'0: ILS method with numeration in search (LAMBDA2)', ...
                    '1: ILS method with shrinking ellipsoid during search (LAMBDA3)' ...
                    '2: ILS method with numeration in search (LAMBDA3)', ...
                    '3: integer rounding method (LAMBDA3)', ...
                    '4: iar_mode = 4; % integer bootstrapping method (LAMBDA3)', ...
                    '5: Partial Ambiguity Resolution (PAR) (LAMBDA3)'}
       
        W_MODE = {'same weight for all the observations', ...
                  'weight based on satellite elevation (sin)' ...
                  'weight based on signal-to-noise ratio', ...
                  'weight based on combined elevation and signal-to-noise ratio', ...
                  'weight based on satellite elevation (exp)'}
              
        IONO_MODE = {'0: no model', ...
                     '1: Geckle and Feen model' ...
                     '2: Klobuchar model', ...
                     '3: SBAS grid'}
                 
        TROPO_MODE = {'0: no model', ...
                      '1: Saastamoinen model (with standard atmosphere parameters)' ...
                      '2: Saastamoinen model (with Global Pressure Temperature model)'}         

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
        sigma_ph_if = 0.009;
        
        % Std of apriori receiver clock
        sigma0_clock = 4.5e-09        
        % Std of receiver clock
        sigma0_r_clock = 1e3;
        
        % Signal-to-noise ratio threshold [dB]
        snr_thr = 0;

        % Cut-off [degrees]
        cutoff = 10;
        
        %------------------------------------------------------------------
        % PRE PROCESSING 
        %------------------------------------------------------------------
        
        % Cycle slip threshold (pre-processing) [cycles]
        cs_thr_pre_pro = 1;
                
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
        flag_iar_auto_mu = true;
        
        % Flag for enabling the default value for P0
        flag_iar_default_P0 = true;
        
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
        
        
        % DTM flag (use / not use) DTM
        flag_dtm = false;

        % Folder containing DTM files
        dtm_dir = '../data/dtm';

        % Std of DEM height [m]
        sigma_dtm = 0.03  % (maximize to disable DEM usage: e.g. 1e30)
        
        % Parameters common to all DTM tiles
        dtm_tile_header = struct('nrows', 0, 'ncols', 0, 'cellsize', 0, 'nodata', 0);
        
        % Parameters used to georeference every DTM tile
        dtm_tile_georef = zeros(1,1,4);        
    end
    
    methods
        function obj = Processing_Settings()
            % Creator
            obj.postImportInit();
        end
    end
    
    methods (Access = 'protected')
        function postImportInit(obj)
            % Operations to run after the import of new parameters
            obj.init_dtm();
        end
    end
    
    methods 
        function import(obj, settings)
            % This function import processing settings from another setting object or ini file
            if isa(settings, 'Ini_Manager')                
                obj.sigma_ph = settings.getData('sigma_ph');
                obj.sigma_ph_if = settings.getData('sigma_ph_if');
                obj.sigma0_clock = settings.getData('sigma0_clock');
                obj.sigma0_r_clock = settings.getData('sigma0_r_clock');
                obj.snr_thr = settings.getData('snr_thr');

                obj.cs_thr_pre_pro = settings.getData('cs_thr_pre_pro');

                obj.cutoff = settings.getData('cutoff');
                obj.w_mode = settings.getData('w_mode');
                obj.cs_thr = settings.getData('cs_thr');
                tmp = settings.getData('w_snr');
                obj.w_snr = struct('a', tmp(1), 'zero', tmp(2), 'one', tmp(2), 'A', tmp(4));
                obj.pp_spp_thr = settings.getData('pp_spp_thr');
                obj.pp_max_code_err_thr = settings.getData('pp_max_code_err_thr');
                obj.pp_max_phase_err_thr = settings.getData('pp_max_phase_err_thr');
                
                obj.sigma0_pos    = settings.getData('sigma0_pos');
                obj.sigma_vel_ENU = settings.getData('sigma_vel_ENU');
                obj.sigma_vel_mod = settings.getData('sigma_vel_mod');
                obj.sigma0_tropo  = settings.getData('sigma0_tropo');
                obj.sigma_tropo   = settings.getData('sigma_tropo');
                obj.kf_min_n_sat  = settings.getData('kf_min_n_sat');
                obj.kf_order      = settings.getData('kf_order');          

                obj.iar_restart_mode = settings.getData('iar_restart_mode');
                obj.iar_mode = settings.getData('iar_mode');
                obj.iar_P0 = settings.getData('iar_P0');
                obj.sigma0_N = settings.getData('sigma0_N');
                obj.iar_mu = settings.getData('iar_mu');
                obj.flag_iar_auto_mu = settings.getData('flag_iar_auto_mu');
                obj.flag_iar_default_P0 = settings.getData('flag_iar_default_P0');
                obj.iono_model = settings.getData('iono_model');
                obj.tropo_model = settings.getData('tropo_model');
                obj.dtm_dir    = settings.getData('dtm_dir');
                obj.flag_dtm = settings.getData('flag_dtm');
                obj.sigma_dtm = settings.getData('sigma_dtm');
            else                
                obj.sigma_ph = settings.sigma_ph;
                obj.sigma_ph_if = settings.sigma_ph_if;
                obj.sigma0_clock = settings.sigma0_clock;
                obj.sigma0_r_clock = settings.sigma0_r_clock;
                obj.snr_thr = settings.snr_thr;
                obj.cutoff = settings.cutoff;

                obj.cs_thr_pre_pro = settings.cs_thr_pre_pro;

                obj.w_mode = settings.w_mode;
                obj.cs_thr = settings.cs_thr;
                obj.w_snr = settings.w_snr;
                obj.pp_spp_thr = settings.pp_spp_thr;
                obj.pp_max_code_err_thr = settings.pp_max_code_err_thr;
                obj.pp_max_phase_err_thr = settings.pp_max_phase_err_thr;
                
                obj.sigma0_pos    = settings.sigma0_pos;
                obj.sigma_vel_ENU = settings.sigma_vel_ENU;
                obj.sigma_vel_mod = settings.sigma_vel_mod;
                obj.sigma0_tropo  = settings.sigma0_tropo;
                obj.sigma_tropo   = settings.sigma_tropo;
                obj.kf_min_n_sat  = settings.kf_min_n_sat;
                obj.kf_order      = settings.kf_order;
                
                obj.iar_restart_mode = settings.iar_restart_mode;
                obj.iar_mode = settings.iar_mode;
                obj.iar_P0 = settings.iar_P0;
                obj.sigma0_N = settings.sigma0_N;
                obj.iar_mu = settings.iar_mu;
                obj.flag_iar_auto_mu = settings.flag_iar_auto_mu;
                obj.flag_iar_default_P0 = settings.flag_iar_default_P0;
                obj.iono_model = settings.iono_model;
                obj.tropo_model = settings.tropo_model;
                obj.dtm_dir    = settings.dtm_dir;
                obj.flag_dtm = settings.flag_dtm;
                obj.sigma_dtm = settings.sigma_dtm;
            end
            % Call to Super Methods
            obj.import@IO_Settings(settings);
            obj.import@PrePro_Settings(settings);
            obj.import@KF_Settings(settings);
            
            obj.postImportInit();
        end
        
        function str = toString(obj, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end
            
            str = obj.toString@IO_Settings(str);
            str = [str '---- RECEIVERS -----------------------------------------------------------' 10 10];
            str = [str sprintf(' STD of phase observations [m]:                    %g\n', obj.sigma_ph)];
            str = [str sprintf(' STD of iono-free phase observations [m]:          %g\n\n', obj.sigma_ph_if)];
            str = [str sprintf(' STD of apriori receiver clock:                    %g\n', obj.sigma0_clock)];
            str = [str sprintf(' STD of receiver clock:                            %g\n\n', obj.sigma0_r_clock)];
            str = [str sprintf(' Signal-to-noise ratio threshold [dB]:             %d\n\n', obj.snr_thr)];
            str = [str sprintf(' Cut-off [degrees]:                                %d\n\n', obj.cutoff)];
            
            str = [str '---- PRE PROCESSING ------------------------------------------------------' 10 10];
            str = [str sprintf(' Cycle slip threshold [cycles]                     %g\n\n', obj.cs_thr_pre_pro)];

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

            str = [str '---- KALMAN FILTER PARAMETERS --------------------------------------------' 10 10];
            str = [str sprintf(' STD of initial state [m]:                         %g\n', obj.sigma0_pos)];
            str = [str sprintf(' STD of ENU velocity [m]:                          %g %g %g\n', struct2array(obj.sigma_vel_ENU))];
            str = [str sprintf(' STD of 3D velocity modulus [m]:                   %g\n\n', obj.sigma_vel_mod)];
            str = [str sprintf(' STD of apriori tropospheric delay:                %g\n', obj.sigma0_tropo)];
            str = [str sprintf(' STD of tropospheric delay:                        %g\n\n', obj.sigma_tropo)];
            str = [str sprintf(' Minimum number of satellite per epoch:            %d\n', obj.kf_min_n_sat)];
            str = [str sprintf(' Oreder of the KF:                                 %d\n\n', obj.kf_order)];
            
            str = [str '---- ABIGUITY (IAR) -------------------------------------------------------' 10 10];
            str = [str sprintf(' Ambiguity restart mode:                           %d\n\n', obj.iar_restart_mode)];
            str = [str sprintf(' Using method: %s\n\n', obj.IAR_MODE{obj.iar_mode+1})];
            str = [str sprintf(' User defined fixed failure rate (methods 1,2):    %g\n', obj.iar_P0)];
            str = [str sprintf(' User defined minimum success rate (for method 5): %g\n', obj.iar_P0)];
            str = [str sprintf(' STD of apriori ambiguity combinations [cycles]:   %d\n\n', obj.sigma0_N)];
            str = [str sprintf(' User defined threshold for ratio test:            %g\n', obj.iar_mu)];
            str = [str sprintf(' Automatic determination of mu:                    %d\n', obj.flag_iar_auto_mu)];
            str = [str sprintf(' Use default value for P0:                         %d\n\n', obj.flag_iar_default_P0)];
            str = [str '---- ATHMOSPHERE ---------------------------------------------------------' 10 10];
            str = [str sprintf(' Ionospheric model:  %s\n', obj.IONO_MODE{obj.iono_model+1})];
            str = [str sprintf(' Tropospheric model: %s\n\n', obj.TROPO_MODE{obj.tropo_model+1})];
            str = [str '---- DEM -----------------------------------------------------------------' 10 10];
            str = [str sprintf(' Use DTM:                                          %d\n', obj.flag_dtm)];
            str = [str sprintf(' Folder containing DTM data:                       %s\n', obj.dtm_dir)];
            str = [str sprintf(' STD of DEM model [m]:                             %g\n', obj.sigma_dtm)];
        end
        
        function str_cell = export(obj, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the obj
            if (nargin == 1)
                str_cell = {};
            end
                        
            str_cell = obj.export@IO_Settings(str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            str_cell = Ini_Manager.toIniStringSection('RECEIVERS', str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of phase observations [m]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma_ph', obj.sigma_ph, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of iono-free phase observations [m', str_cell);
            str_cell = Ini_Manager.toIniString('sigma_ph_if', obj.sigma_ph_if, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of apriori receiver clock', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_clock', obj.sigma0_clock, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of receiver clock', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_r_clock', obj.sigma0_r_clock, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Signal-to-noise ratio threshold [dB]', str_cell);
            str_cell = Ini_Manager.toIniString('snr_thr', obj.snr_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Cut-off [degrees]', str_cell);
            str_cell = Ini_Manager.toIniString('cutoff', obj.cutoff, str_cell);            
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            str_cell = Ini_Manager.toIniStringSection('PRE_PROCESSING', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Cycle slip threshold [cycles]', str_cell);
            str_cell = Ini_Manager.toIniString('cs_thr_pre_pro', obj.cs_thr_pre_pro, str_cell);

            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            str_cell = Ini_Manager.toIniStringSection('PROCESSING', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Processing using weighting mode:', str_cell);
            str_cell = Ini_Manager.toIniString('w_mode', obj.w_mode, str_cell);
            for (i = 1 : numel(obj.W_MODE))
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %d: %s', i - 1, obj.W_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Cycle slip threshold (processing) [cycles]', str_cell);
            str_cell = Ini_Manager.toIniString('cs_thr', obj.cs_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Weight function parameters (when based on SNR): a / 0 / 1 / A', str_cell);
            str_cell = Ini_Manager.toIniString('w_snr', struct2array(obj.w_snr), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            str_cell = Ini_Manager.toIniStringSection('THRESHOLDS', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Threshold on code LS estimation error [m]', str_cell);
            str_cell = Ini_Manager.toIniString('pp_spp_thr', obj.pp_spp_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Threshold on maximum residual of code obs [m]', str_cell);
            str_cell = Ini_Manager.toIniString('pp_max_code_err_thr', obj.pp_max_code_err_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Threshold on maximum residual of phase obs [m]', str_cell);
            str_cell = Ini_Manager.toIniString('pp_max_phase_err_thr', obj.pp_max_phase_err_thr, str_cell);
            
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringSection('KALMAN_FILTER', str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of initial state [m]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_pos', obj.sigma0_pos, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of ENU velocity [m]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma_vel_ENU', struct2array(obj.sigma_vel_ENU), str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of 3D velocity modulus [m]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma_vel_mod', obj.sigma_vel_mod, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of apriori tropospheric delay', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_tropo', obj.sigma0_tropo, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of tropospheric delay', str_cell);
            str_cell = Ini_Manager.toIniString('sigma_tropo', obj.sigma_tropo, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Minimum number of satellite per epoch', str_cell);
            str_cell = Ini_Manager.toIniString('kf_min_n_sat', obj.kf_min_n_sat, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Oreder of the KF', str_cell);
            str_cell = Ini_Manager.toIniString('kf_order', obj.kf_order, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            str_cell = Ini_Manager.toIniStringSection('ABIGUITY', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Ambiguity restart mode', str_cell);
            str_cell = Ini_Manager.toIniString('iar_restart_mode', obj.iar_restart_mode, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Ambiguity detection mode', str_cell);
            str_cell = Ini_Manager.toIniString('iar_mode', obj.iar_mode, str_cell);
            for (i = 1 : numel(obj.IAR_MODE))
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', obj.IAR_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('User defined fixed failure rate (methods 1,2) / user defined minimum success rate (for method 5)', str_cell);
            str_cell = Ini_Manager.toIniString('iar_P0', obj.iar_P0, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of apriori ambiguity combinations [cycles]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_N', obj.sigma0_N, str_cell);
            str_cell = Ini_Manager.toIniStringComment('User defined threshold for ratio test', str_cell);
            str_cell = Ini_Manager.toIniString('iar_mu', obj.iar_mu, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Automatic determination of mu', str_cell);
            str_cell = Ini_Manager.toIniString('flag_iar_auto_mu', obj.flag_iar_auto_mu, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use default value for P0', str_cell);
            str_cell = Ini_Manager.toIniString('flag_iar_default_P0', obj.flag_iar_default_P0, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            str_cell = Ini_Manager.toIniStringSection('ATHMOSPHERE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Ionospheric model', str_cell);
            str_cell = Ini_Manager.toIniString('iono_model', obj.iono_model, str_cell);
            for (i = 1 : numel(obj.IONO_MODE))
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', obj.IONO_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Tropospheric model', str_cell);
            str_cell = Ini_Manager.toIniString('tropo_model', obj.tropo_model, str_cell);
            for (i = 1 : numel(obj.IONO_MODE))
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', obj.IONO_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            str_cell = Ini_Manager.toIniStringSection('DEM', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use DTM (true/false)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_dtm', obj.flag_dtm, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Folder containing DTM data', str_cell);
            str_cell = Ini_Manager.toIniString('dtm_dir', obj.dtm_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of DEM model [m]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma_dtm', obj.sigma_dtm, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end
        
        function init_dtm(obj, dtm_dir)
            % Try to load default DTM values (if flag_dtm is on)
            if obj.flag_dtm
                if nargin == 1
                    dtm_dir = obj.dtm_dir; % from superclass IO_Settings
                end
                try
                    load([dtm_dir filesep 'tiles' filesep 'tile_header'], 'tile_header');
                    obj.tile_header = tile_header;
                    obj.logger.addMessage(sprintf(' - DTM tile header in %s have been read', [dtm_dir filesep 'tiles' filesep 'tile_header']));
                    load([dtm_dir filesep 'tiles' filesep 'tile_georef'], 'tile_georef');
                    obj.tile_georef = tile_georef;
                    obj.logger.addMessage(sprintf(' - DTM tile georef in %s have been read', [dtm_dir filesep 'tiles' filesep 'tile_georef']));
                catch
                    obj.logger.addWarning(sprintf('Failed to read DTM stored in %s', [dtm_dir '/tiles/']));
                    % use default zeroes values
                end
            end
        end
        
    end
end
