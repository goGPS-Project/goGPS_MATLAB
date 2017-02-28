%   CLASS Main_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the Processing parameters
%
% EXAMPLE
%   state = Main_Settings();
%
% SEE ALSO
%   - GO_Settings     
%         it is a singleton class that store the status of the
%         Main_Settings to be used in one goGPS session
%
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___ 
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.0
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

classdef Main_Settings < Settings_Interface & IO_Settings & Mode_Settings
    
    % Default values for each field - useful to restore corrupted fields
    properties (Constant, Access = 'private')
        % RECEIVER DEFAULT PARAMETERS
        STD_CODE = 3;                                   % Std of code observations [m]
        STD_PHASE = 0.03;                               % Std of phase observations [m]
        STD_PHASE_IF = 0.009;                           % Std of iono-free phase observations [m]
        SIGMA0_CLOCK = 4.47e-9;                         % Std of a priori receiver clock
        SIGMA0_R_CLOCK = 31                             % Std of receiver clock
        FLAG_RINEX_MPOS = true;                         % Flag to read the position of the master form the RINEX file (deprecate)
        MPOS = struct('X', 0, 'Y', 0, 'Z', 0);          % Default master position (these are overrided when the coordinates are specified elsewhere) (deprecate)
        
        % DATA SELECTION        
        CC = Constellation_Collector('G');              % object containing info on the activated constellations
        P_RATE = 1;                                     % Minimum processing rate [s]
        MIN_N_SAT = 2;                                  % Minimum number of satellites to be used in the Kalman filter
        CUT_OFF = 10;                                   % Cut-off [degrees]
        SNR_THR = 0;                                    % Signal-to-noise ratio threshold [dB]
        FLAG_OCEAN = false;                             % Flag for enabling the usage of ocean tides modeling
        MIN_ARC = 10;                                   % Minimum length an arc (a satellite to be used must be seen for a number of consecutive epochs greater than this value)
        
        % PRE PROCESSING
        FLAG_PRE_PRO = true;                            % Flag for enabling pre-processing
        CS_THR_PRE_PRO = 3;                             % Cycle slip threshold (pre-processing) [cycles]
        
        % OUTLIER DETECTION
        FLAG_OUTLIER = true;                            % Flag for enabling outlier detection
        PP_SPP_THR = 4;                                 % Threshold on the code point-positioning least squares estimation error [m]
        PP_MAX_CODE_ERR_THR = 30;                       % Threshold on the maximum residual of code observations [m]
        PP_MAX_PHASE_ERR_THR = 0.2;                     % Threshold on the maximum residual of phase observations [m]
        
        % PROCESSING PARAMETERS
        FLAG_TROPO = false;                             % Flag for enabling the computation of thropospheric derived informations
        W_MODE = 1;                                     % Parameter used to select the weightening mode for GPS observations
                                                        %  - weights = 0: same weight for all the observations
                                                        %  - weights = 1: weight based on satellite elevation (sin)
                                                        %  - weights = 2: weight based on signal-to-noise ratio
                                                        %  - weights = 3: weight based on combined elevation and signal-to-noise ratio
                                                        %  - weights = 4: weight based on satellite elevation (exp)
        W_SNR = struct('a', 30, 'zero', 10, 'one', 50, 'A', 30); % Weight function parameters (when based on SNR)
        CS_THR = 1;                                     % Cycle slip threshold (processing) [cycles]
        FLAG_IONOFREE = false;                          % Flag for enabling the usage of iono-free combination
        CONSTRAIN = false;                              % Constrain the solution using a reference path
        STOP_GO_STOP = false;                           % This flag add the possibility to process in stop go stop mode        

        % INTEGER AMBIGUITY RESOLUTION
        FLAG_IAR = 0;                                   % Flag for enabling the automatic detection of cycle sleep
        IAR_RESTART_MODE = 2;                           % Ambiguity restart mode
                                                        % - iar_restart_mode = 0; % Observed code - phase difference
                                                        % - iar_restart_mode = 1; % Kalman-predicted code - phase difference
                                                        % - iar_restart_mode = 2; % Least squares adjustment
        IAR_MODE = 1;                                   % Parameter used to select Integer Least Squares estimator
                                                        % - iar_mode = 0: % ILS method with numeration in search (LAMBDA2)
                                                        % - iar_mode = 1: % ILS method with shrinking ellipsoid during search (LAMBDA3)
                                                        % - iar_mode = 2: % ILS method with numeration in search (LAMBDA3)
                                                        % - iar_mode = 3: % integer rounding method (LAMBDA3)
                                                        % - iar_mode = 4: % integer bootstrapping method (LAMBDA3)
                                                        % - iar_mode = 5: % Partial Ambiguity Resolution (PAR) (LAMBDA3)
        IAR_P0 = 0.001;                                 % User defined fixed failure rate (for methods 1,2) or minimum required success rate (for method 5)
        SIGMA0_N = 31;                                  % Std of a priori ambiguity combinations [cycles]
        IAR_MU = 0.5;                                   % User defined threshold for ratio test
        FLAG_IAR_AUTO_MU = true;                        % Flag for enabling the automatic determination of mu
        FLAG_IAR_DEFAULT_P0 = true;                     % Flag for enabling the default value for P0
        FLAG_DOPPLER = false;                           % Flag for using doppler-predicted phase range for detecting cycle slips
        
        % KF
        KF_ORDER = 1;                                   % Order of the dynamic model polynomial
                                                        % - kf_order = 0; static
                                                        % - kf_order = 1; constant velocity
                                                        % - kf_order = 2; constant acceleration
                                                        % - kf_order = 3; variable (stop-go-stop)

        % RECEIVER POSITION / MOTION 
        SIGMA0_K_POS = 1;                               % Std of initial state [m]
        STD_K_ENU = struct('E', 0.5, 'N', 0.5, 'U', 0.1); % Std of ENU coordinates variation [m] / [m/s] / [m/s^2]
        STD_K_VEL_MOD = 0.1                             % Std of 3D modulus variation [m] / [m/s] / [m/s^2]
                                                        
        % ATMOSPHERE
        SIGMA0_TROPO = 0.1;                             % Std of a priori tropospheric delay
        STD_TROPO = 4.5644e-4;                          % Std of tropospheric delay
        IONO_MODEL = 2;                                 % Ionospheric model to be used (0: none, 1: Geckle and Feen, 2: Klobuchar, 3: SBAS)
                                                        % - iono_model = 0: no model
                                                        % - iono_model = 1: Geckle and Feen model
                                                        % - iono_model = 2: Klobuchar model
                                                        % - iono_model = 3: SBAS grid
        TROPO_MODEL = 1;                                % Tropospheric model to be used (0: none, 1: Saastamoinen std parameters, 2: Saastamoinen global pararameters)
                                                        % - tropo_model = 0: no model
                                                        % - tropo_model = 1: Saastamoinen model (with standard atmosphere parameters)
                                                        % - tropo_model = 2: Saastamoinen model (with Global Pressure Temperature model)

        % DTM
        FLAG_DTM = false;                               % DTM flag (use / do not use)
        STD_DTM = 0.03;                                 % Std of DEM height [m]
        ANTENNA_H = 0;                                  % Elevation of the antenna above ground
        DTM_TILE_HEADER = struct('nrows', 0, 'ncols', 0, 'cellsize', 0, 'nodata', 0); % Parameters common to all DTM tiles
        DTM_TILE_GEOREF = zeros(1,1,4);                 % Parameters used to georeference every DTM tile

        % GUI
        PLOT_PROC = true;                               % plot during processing
        PLOT_REF_PATH = false;                          % plot ref during processing
        PLOT_SKYPLOT_SNR = false;                       % plot sky plot during processing
        PLOT_ERR_ELLIPSE = false;                       % plot error_ellipse
        PLOT_AMBIGUITIES = false;                       % plot ambiguities
        PLOT_MASTER = false;                            % plot master station
        PLOT_GOOGLE_EARTH = false;                      % plot on google earth
        
        % CAPTURE
        C_N_RECEIVERS = 1;                              % Number of receiver to use for capturing data
        C_RATE = 1;                                     % Capture rate [s]
        C_PRTC = 1;                                     % Array with the size of c_n_receivers
                                                        % - c_prtc = 1: UBX (u-blox)
                                                        % - c_prtc = 2: iTalk (Fastrax)
                                                        % - c_prtc = 3: SkyTraq
                                                        % - c_prtc = 4: BINR (NVS)
        C_COM_ADDR = {'/dev/tty.lpss-serial1'};         % Cell array with the com address of each receiver to be used
        
        % NTRIP
        FLAG_NTRIP = false;                             % NTRIP flag (use / do not use) 
        NTRIP = struct('ip_addr', '127.0.0.1', ...      % Struct containing NTRIP parameters:
                       'port', '2101', ...
                       'mountpoint', '/', ...
                       'username', 'user', ...
                       'password', 'pswd', ...
                       'approx_position', struct('lat', 0, 'lon', 0, 'h', 0));
    end
    
    properties (Constant, Access = 'protected')
        % id to string of Kalman Filter dynamic modes
        DYN_SMODE_RT = {'0: constant', ...
                        '1: variable'}
        DYN_SMODE_PP = {'0: static', ...
                        '1: constant velocity' ...
                        '2: constant acceleration', ...
                        '3: variable (stop-go-stop)'}

        % id to string of IAR modes
        IAR_SMODE = {'0: ILS method with numeration in search (LAMBDA2)', ...
                     '1: ILS method with shrinking ellipsoid during search (LAMBDA3)' ...
                     '2: ILS method with numeration in search (LAMBDA3)', ...
                     '3: integer rounding method (LAMBDA3)', ...
                     '4: integer bootstrapping method (LAMBDA3)', ...
                     '5: Partial Ambiguity Resolution (PAR) (LAMBDA3)'}
                
        % id to string of IAR restart function
        IAR_SRESTART = {'0: Observed code - phase difference', ...
                        '1: Kalman-predicted code - phase difference', ...
                        '2: Least squares adjustment'}
       
        % id to string of weight functions
        W_SMODE = {'same weight for all the observations', ...
                   'weight based on satellite elevation (sin)' ...
                   'weight based on satellite elevation (exp)', ...
                   'weight based on signal-to-noise ratio', ...
                   'weight based on combined elevation and signal-to-noise ratio'}
              
        % id to string of ionospheric models
        IONO_SMODE = {'0: no model', ...
                      '1: Geckle and Feen model' ...
                      '2: Klobuchar model', ...
                      '3: SBAS grid'}
                 
        % id to string of tropospheric models 
        TROPO_SMODE = {'0: no model', ...
                       '1: Saastamoinen model (with standard atmosphere parameters)' ...
                       '2: Saastamoinen model (with Global Pressure Temperature model)'} 
                          
        % id to string of capturing protocols
        C_SPROTOCOL = {'1: UBX (u-blox)', ...
                       '2: iTalk (Fastrax)', ...
                       '3: SkyTraq', ...
                       '4: BINR (NVS)'}

        % processing rates used in UI -> to convert UI to settings format
        UI_P_SRATE = goGUIclass.UI_P_SRATE;
        
        % capture rates used in UI -> to convert UI to settings format
        UI_C_SRATE = goGUIclass.UI_C_SRATE;
    end
    
    %  Processing parameters
    % ------------------------------
    % note that some of the processing parameters are actually properties of other objects 
    % e.g. std_k_ENU depend on the receiver motion
    % However for various reasons I can decide to ignore the "correct"
    % values and redefine them. So the values are also
    % stored here as parameters used for a specific processing
    properties (SetAccess = private, GetAccess = public)
        
        %------------------------------------------------------------------
        % RECEIVER DEFAULT PARAMETERS
        %------------------------------------------------------------------
        % These values are kept if not specified elsewhere
        
        % Std of code observations [m]
        std_code = Main_Settings.STD_CODE; 

        % Std of phase observations [m]
        std_phase = Main_Settings.STD_PHASE; % (maximize to obtain a code-based solution)
        % Std of iono-free phase observations [m]
        std_phase_if = Main_Settings.STD_PHASE_IF;
        
        % Std of a priori receiver clock
        sigma0_clock = Main_Settings.SIGMA0_CLOCK;       
        % Std of receiver clock
        sigma0_r_clock = Main_Settings.SIGMA0_R_CLOCK;
        
        % Flag to read the position of the master form the RINEX file (deprecate)
        flag_rinex_mpos = Main_Settings.FLAG_RINEX_MPOS;
        % Default master position (these are overrided when the coordinates are specified elsewhere) (deprecate)
        mpos = Main_Settings.MPOS;
        
        %------------------------------------------------------------------
        % DATA SELECTION
        %------------------------------------------------------------------
        
        % object containing info on the activated constellations
        cc =  Main_Settings.CC;
        % Minimum processing rate [s]
        p_rate = Main_Settings.P_RATE;
        % Minimum number of satellites to be used in the Kalman filter
        min_n_sat = Main_Settings.MIN_N_SAT;        
        % Cut-off [degrees]
        cut_off = Main_Settings.CUT_OFF;
        % Signal-to-noise ratio threshold [dB]
        snr_thr = Main_Settings.SNR_THR;                
        % Flag for enabling the usage of ocean tides modeling
        flag_ocean = Main_Settings.FLAG_OCEAN;        
        % Minimum length an arc (a satellite to be used must be seen for a number of consecutive epochs greater than this value)
        min_arc = Main_Settings.MIN_ARC;

        %------------------------------------------------------------------
        % PRE PROCESSING 
        %------------------------------------------------------------------
        
        % Flag for enabling pre-processing
        flag_pre_pro = Main_Settings.FLAG_PRE_PRO;        
        % Cycle slip threshold (pre-processing) [cycles]
        cs_thr_pre_pro = Main_Settings.CS_THR_PRE_PRO;

        %------------------------------------------------------------------
        % OUTLIER DETECTION
        %------------------------------------------------------------------

        % Flag for enabling outlier detection
        flag_outlier = Main_Settings.FLAG_OUTLIER;
                
        % Threshold on the code point-positioning least squares estimation error [m]
        pp_spp_thr = Main_Settings.PP_SPP_THR;
        % Threshold on the maximum residual of code observations [m]
        pp_max_code_err_thr = Main_Settings.PP_MAX_CODE_ERR_THR;            
        % Threshold on the maximum residual of phase observations [m]
        pp_max_phase_err_thr = Main_Settings.PP_MAX_PHASE_ERR_THR;  
        
        %------------------------------------------------------------------
        % PROCESSING PARAMETERS 
        %------------------------------------------------------------------
        
        % Flag for enabling the computation of thropospheric derived informations
        flag_tropo = Main_Settings.FLAG_TROPO;
        
        % Parameter used to select the weightening mode for GPS observations
        w_mode = Main_Settings.W_MODE;
        %  - weights = 0: same weight for all the observations
        %  - weights = 1: weight based on satellite elevation (sin)
        %  - weights = 2: weight based on signal-to-noise ratio
        %  - weights = 3: weight based on combined elevation and signal-to-noise ratio
        %  - weights = 4: weight based on satellite elevation (exp)
                
        % Weight function parameters (when based on SNR)
        w_snr = Main_Settings.W_SNR;
        
        % Cycle slip threshold (processing) [cycles]
        cs_thr = Main_Settings.CS_THR;

        % Flag for enabling the usage of iono-free combination
        flag_ionofree = Main_Settings.FLAG_IONOFREE;

        % Constrain the solution using a reference path
        constrain = Main_Settings.CONSTRAIN; 
        
        % This flag add the possibility to process in stop go stop mode
        % static convergence / movement / static convergence
        stop_go_stop = Main_Settings.STOP_GO_STOP;
        
        %------------------------------------------------------------------
        % INTEGER AMBIGUITY RESOLUTION
        %------------------------------------------------------------------
        
        % Flag for enabling the automatic detection of cycle sleep
        flag_iar = Main_Settings.FLAG_IAR;

        % Ambiguity restart mode
        iar_restart_mode = Main_Settings.IAR_RESTART_MODE;        
        % - iar_restart_mode = 0; % Observed code - phase difference
        % - iar_restart_mode = 1; % Kalman-predicted code - phase difference
        % - iar_restart_mode = 2; % Least squares adjustment
        
        % Parameter used to select Integer Least Squares estimator
        iar_mode = Main_Settings.IAR_MODE;
        % - iar_mode = 0: % ILS method with numeration in search (LAMBDA2)
        % - iar_mode = 1: % ILS method with shrinking ellipsoid during search (LAMBDA3)
        % - iar_mode = 2: % ILS method with numeration in search (LAMBDA3)
        % - iar_mode = 3: % integer rounding method (LAMBDA3)
        % - iar_mode = 4: % integer bootstrapping method (LAMBDA3)
        % - iar_mode = 5: % Partial Ambiguity Resolution (PAR) (LAMBDA3)

        % User defined fixed failure rate (for methods 1,2) or minimum required success rate (for method 5)
        iar_p0 = Main_Settings.IAR_P0;

        % Std of a priori ambiguity combinations [cycles]
        sigma0_N = Main_Settings.SIGMA0_N;

        % User defined threshold for ratio test
        iar_mu = Main_Settings.IAR_MU;
        
        % Flag for enabling the automatic determination of mu
        flag_iar_auto_mu = Main_Settings.FLAG_IAR_AUTO_MU;
        
        % Flag for enabling the default value for P0
        flag_iar_default_p0 = Main_Settings.FLAG_IAR_DEFAULT_P0;
        
        % Flag for using doppler-predicted phase range for detecting cycle slips
        flag_doppler = Main_Settings.FLAG_DOPPLER;
        
        %------------------------------------------------------------------
        % KF
        %------------------------------------------------------------------
        
        % Order of the dynamic model polynomial
        kf_order = Main_Settings.KF_ORDER;
        % - kf_order = 0; static
        % - kf_order = 1; constant velocity
        % - kf_order = 2; constant acceleration
        % - kf_order = 3; variable (stop-go-stop)
        
        %------------------------------------------------------------------
        % RECEIVER POSITION / MOTION 
        %------------------------------------------------------------------
        
        % Std of initial state [m]
        sigma0_k_pos = Main_Settings.SIGMA0_K_POS;         
        % Std of ENU coordinates variation [m] / [m/s] / [m/s^2]
        std_k_ENU = Main_Settings.STD_K_ENU;
        % Std of 3D modulus variation [m] / [m/s] / [m/s^2]
        std_k_vel_mod = Main_Settings.STD_K_VEL_MOD;
                                
        %------------------------------------------------------------------
        % ATMOSPHERE
        %------------------------------------------------------------------

        % Std of a priori tropospheric delay
        sigma0_tropo = Main_Settings.SIGMA0_TROPO;
        % Std of tropospheric delay
        std_tropo = Main_Settings.STD_TROPO;
                
        % Ionospheric model to be used (0: none, 1: Geckle and Feen, 2: Klobuchar, 3: SBAS)
        iono_model = Main_Settings.IONO_MODEL;
        % - iono_model = 0: no model
        % - iono_model = 1: Geckle and Feen model
        % - iono_model = 2: Klobuchar model
        % - iono_model = 3: SBAS grid
        
        % Tropospheric model to be used (0: none, 1: Saastamoinen std parameters, 2: Saastamoinen global pararameters)
        tropo_model = Main_Settings.TROPO_MODEL;
        % - tropo_model = 0: no model
        % - tropo_model = 1: Saastamoinen model (with standard atmosphere parameters)
        % - tropo_model = 2: Saastamoinen model (with Global Pressure Temperature model)
        
        %------------------------------------------------------------------
        % DTM 
        %------------------------------------------------------------------
                
        % DTM flag (use / do not use)
        flag_dtm = Main_Settings.FLAG_DTM;

        % Folder containing DTM files
        % dtm_dir = '../data/dtm'; -> defined in superclass IO_Settings

        % Std of DEM height [m]
        std_dtm = Main_Settings.STD_DTM  % (maximize to disable DEM usage: e.g. 1e30)
        
        % Elevation of the antenna above ground
        antenna_h = Main_Settings.ANTENNA_H; % when fixing 
        
        % Parameters common to all DTM tiles
        dtm_tile_header = Main_Settings.DTM_TILE_HEADER;
        
        % Parameters used to georeference every DTM tile
        dtm_tile_georef = Main_Settings.DTM_TILE_GEOREF;        
        
        %------------------------------------------------------------------
        % GUI
        %------------------------------------------------------------------
        
        % plot during processing
        plot_proc = Main_Settings.PLOT_PROC;        
        % plot ref during processing
        plot_ref_path = Main_Settings.PLOT_REF_PATH;        
        % plot sky plot during processing
        plot_skyplot_snr = Main_Settings.PLOT_SKYPLOT_SNR;        
        % plot error_ellipse
        plot_err_ellipse = Main_Settings.PLOT_ERR_ELLIPSE;        
        % plot ambiguities
        plot_ambiguities = Main_Settings.PLOT_AMBIGUITIES;        
        % plot master station
        plot_master = Main_Settings.PLOT_MASTER;        
        % plot on google earth
        plot_google_earth = Main_Settings.PLOT_GOOGLE_EARTH;     
                
        %------------------------------------------------------------------
        % CAPTURE
        %------------------------------------------------------------------
        
        % Number of receiver to use for capturing data
        c_n_receivers = Main_Settings.C_N_RECEIVERS;
        
        % Capture rate [s]
        c_rate = Main_Settings.C_RATE;
        
        % Array with the size of c_n_receivers
        c_prtc = Main_Settings.C_PRTC
        % - c_prtc = 1: UBX (u-blox)
        % - c_prtc = 2: iTalk (Fastrax)
        % - c_prtc = 3: SkyTraq
        % - c_prtc = 4: BINR (NVS)
        
        % Cell array with the com address of each receiver to be used
        c_com_addr = Main_Settings.C_COM_ADDR;
        
        %------------------------------------------------------------------
        % NTRIP
        %------------------------------------------------------------------
        
        % NTRIP flag (use / do not use) 
        flag_ntrip = Main_Settings.FLAG_NTRIP;
        
        % Struct containing NTRIP parameters:
        ntrip = Main_Settings.NTRIP;
    end
    
    % =========================================================================
    %  INIT
    % =========================================================================    
    methods
        function this = Main_Settings()
            % Creator
            % SYNTAX: s_obj = Main_Settings();
            this.postImportInit();
            if exist(this.LAST_SETTINGS, 'file')
                this.logger.addMessage('Importing settings from last settings file')
                this.importIniFile(this.LAST_SETTINGS);
            end
        end
    end
    
    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================    
    methods (Access = 'public')
        function import(this, state)
            % This function import processing settings from another setting object or ini file
            % SYNTAX: s_obj.import(state)
            
            if isa(state, 'Ini_Manager')         
                % RECEIVER DEFAULT PARAMETERS
                this.std_code = state.getData('std_code');
                this.std_phase = state.getData('std_phase');
                this.std_phase_if = state.getData('std_phase_if');
                this.sigma0_clock = state.getData('sigma0_clock');
                this.sigma0_r_clock = state.getData('sigma0_r_clock');

                this.flag_rinex_mpos = state.getData('flag_rinex_mpos');                
                tmp = state.getData('mpos_XYZ');
                this.mpos = struct('X', tmp(1), 'Y', tmp(2), 'Z', tmp(3));
                
                % DATA SELECTION
                this.cc.import(state);
                this.p_rate = state.getData('p_rate');
                this.min_n_sat  = state.getData('min_n_sat');
                this.cut_off = state.getData('cut_off');
                this.snr_thr = state.getData('snr_thr');                
                this.flag_ocean = state.getData('flag_ocean');
                this.min_arc = state.getData('min_arc');
                
                % PRE PROCESSING                
                this.flag_pre_pro = state.getData('flag_pre_pro');
                this.cs_thr_pre_pro = state.getData('cs_thr_pre_pro');
                
                % OUTLIER DETECTION
                this.flag_outlier = state.getData('flag_outlier');
                this.pp_spp_thr = state.getData('pp_spp_thr');
                this.pp_max_code_err_thr = state.getData('pp_max_code_err_thr');
                this.pp_max_phase_err_thr = state.getData('pp_max_phase_err_thr');

                % PROCESSING PARAMETERS
                this.flag_tropo = state.getData('flag_tropo');
                this.w_mode = state.getData('w_mode');
                tmp = state.getData('w_snr');
                this.w_snr = struct('a', tmp(1), 'zero', tmp(2), 'one', tmp(2), 'A', tmp(4));
                this.cs_thr = state.getData('cs_thr');
                this.flag_ionofree = state.getData('flag_ionofree');
                this.constrain = state.getData('constrain');
                this.stop_go_stop = state.getData('stop_go_stop');
                
                % INTEGER AMBIGUITY RESOLUTION
                this.flag_iar = state.getData('flag_iar');
                this.iar_restart_mode = state.getData('iar_restart_mode');
                this.iar_mode = state.getData('iar_mode');
                this.iar_p0 = state.getData('iar_p0');
                this.sigma0_N = state.getData('sigma0_N');
                this.iar_mu = state.getData('iar_mu');
                this.flag_iar_auto_mu = state.getData('flag_iar_auto_mu');
                this.flag_iar_default_p0 = state.getData('flag_iar_default_p0');
                this.flag_doppler = state.getData('flag_doppler');
                
                % KF
                this.kf_order = state.getData('kf_order');

                % RECEIVER POSITION / MOTION
                this.sigma0_k_pos    = state.getData('sigma0_k_pos');
                tmp = state.getData('std_k_ENU');
                this.std_k_ENU = struct('E', tmp(1), 'N', tmp(2), 'U', tmp(3));
                this.std_k_vel_mod = state.getData('std_k_vel_mod');
                                
                % ATMOSPHERE
                this.sigma0_tropo  = state.getData('sigma0_tropo');
                this.std_tropo   = state.getData('std_tropo');

                this.iono_model = state.getData('iono_model');
                this.tropo_model = state.getData('tropo_model');
                
                % DTM
                %this.dtm_dir    = state.getData('dtm_dir');
                this.flag_dtm = state.getData('flag_dtm');
                this.std_dtm = state.getData('std_dtm');
                this.antenna_h = state.getData('antenna_h');
                
                % GUI
                this.plot_proc = state.getData('plot_proc');
                this.plot_ref_path = state.getData('plot_ref_path');
                this.plot_skyplot_snr = state.getData('plot_skyplot_snr');
                this.plot_err_ellipse = state.getData('plot_err_ellipse');
                this.plot_ambiguities = state.getData('plot_ambiguities');
                this.plot_master = state.getData('plot_master');
                this.plot_google_earth = state.getData('plot_google_earth');
                
                % CAPTURE
                this.c_n_receivers = state.getData('c_n_receivers');
                this.c_rate = state.getData('c_rate');
                for r = 1 : this.c_n_receivers
                    this.c_prtc(r) = state.getData(sprintf('c_prtc_%02d', r));
                    this.c_com_addr{r} = state.getData(sprintf('c_com_addr_%02d', r));
                end
                
                % NTRIP
                this.flag_ntrip = state.getData('flag_ntrip');
                this.ntrip = struct('ip_addr', state.getData('NTRIP','ip_addr'), ...
                       'port', state.getData('NTRIP','port'), ...
                       'mountpoint', state.getData('NTRIP','mountpoint'), ...
                       'username', state.getData('NTRIP','username'), ...
                       'password', state.getData('NTRIP','password'), ...
                       'approx_position', struct('lat', state.getData('NTRIP','ntrip_lat'), 'lon', state.getData('NTRIP','ntrip_lon'), 'h', state.getData('NTRIP','ntrip_h')));                
            else
                % RECEIVER DEFAULT PARAMETERS
                this.std_code = state.std_code;
                this.std_phase = state.std_phase;
                this.std_phase_if = state.std_phase_if;
                this.sigma0_clock = state.sigma0_clock;
                this.sigma0_r_clock = state.sigma0_r_clock;
                this.flag_rinex_mpos = state.flag_rinex_mpos;
                this.mpos = state.mpos;
                
                % DATA SELECTION
                this.cc.import(state.cc);
                this.p_rate = state.p_rate;
                this.min_n_sat  = state.min_n_sat;
                this.cut_off = state.cut_off;
                this.snr_thr = state.snr_thr;
                this.flag_ocean = state.flag_ocean;
                this.min_arc = state.min_arc;
                
                % PRE PROCESSING
                this.flag_pre_pro = state.flag_pre_pro;
                this.cs_thr_pre_pro = state.cs_thr_pre_pro;
                
                % OUTLIER DETECTION                
                this.flag_outlier = state.flag_outlier;
                this.pp_spp_thr = state.pp_spp_thr;
                this.pp_max_code_err_thr = state.pp_max_code_err_thr;
                this.pp_max_phase_err_thr = state.pp_max_phase_err_thr;

                % PROCESSING PARAMETERS                
                this.flag_tropo = state.flag_tropo;
                this.w_mode = state.w_mode;
                this.w_snr = state.w_snr;
                this.cs_thr = state.cs_thr;
                this.flag_ionofree = state.flag_ionofree;
                this.constrain = state.constrain;
                this.stop_go_stop = state.stop_go_stop;
                
                % INTEGER AMBIGUITY RESOLUTION
                this.flag_iar = state.flag_iar;
                this.iar_restart_mode = state.iar_restart_mode;
                this.iar_mode = state.iar_mode;
                this.iar_p0 = state.iar_p0;
                this.sigma0_N = state.sigma0_N;
                this.iar_mu = state.iar_mu;
                this.flag_iar_auto_mu = state.flag_iar_auto_mu;
                this.flag_iar_default_p0 = state.flag_iar_default_p0;
                this.flag_doppler = state.flag_doppler;                                
                
                % KF
                this.kf_order = state.kf_order;

                % RECEIVER POSITION / MOTION 
                this.sigma0_k_pos = state.sigma0_k_pos;
                this.std_k_ENU = state.std_k_ENU;
                this.std_k_vel_mod = state.std_k_vel_mod;
                       
                % ATMOSPHERE
                this.sigma0_tropo = state.sigma0_tropo;
                this.std_tropo = state.std_tropo;                                
                this.iono_model = state.iono_model;
                this.tropo_model = state.tropo_model;
                        
                % DTM 
                this.flag_dtm = state.flag_dtm;
                %this.dtm_dir    = state.dtm_dir;
                this.std_dtm = state.std_dtm;
                this.antenna_h = state.antenna_h;
                
                % GUI
                this.plot_proc = state.plot_proc;
                this.plot_ref_path = state.plot_ref_path;
                this.plot_skyplot_snr = state.plot_skyplot_snr;
                this.plot_err_ellipse = state.plot_err_ellipse;
                this.plot_ambiguities = state.plot_ambiguities;
                this.plot_master = state.plot_master;
                this.plot_google_earth = state.plot_google_earth;
                
                % CAPTURE
                this.c_n_receivers = state.c_n_receivers;
                this.c_rate = state.c_rate;
                this.c_prtc = state.c_prtc;
                this.c_com_addr = state.c_com_addr;
                
                % NTRIP
                this.flag_ntrip = state.flag_ntrip;
                this.ntrip = state.ntrip;
            end

            % Call to Super Methods
            this.import@Mode_Settings(state);
            this.import@IO_Settings(state);
            
            this.check(); % check after import
            this.postImportInit();
        end
        
        function str = toString(this, str)
            % Display the satellite system in use
            % SYNTAX: s_obj.toString(str)

            if (nargin == 1)
                str = '';
            end
            
            str = this.toString@IO_Settings(str);
            str = [str '---- RECEIVERS -----------------------------------------------------------' 10 10];
            str = [str sprintf(' Default STD of code observations [m]:             %g\n', this.std_code)];
            str = [str sprintf(' Default STD of phase observations [m]:            %g\n', this.std_phase)];
            str = [str sprintf(' Default STD of iono-free phase observations [m]:  %g\n', this.std_phase_if)];
            str = [str sprintf(' Default STD of a priori receiver clock:           %g\n', this.sigma0_clock)];
            str = [str sprintf(' Default STD of receiver clock:                    %g\n', this.sigma0_r_clock)];
            str = [str sprintf(' Read master position from RINEX:                  %d\n', this.flag_rinex_mpos)];
            str = [str sprintf(' Default Master position (when not ovewrrided):\n')];
            str = [str sprintf('   - X:  %12.4f [m]\n', this.mpos.X)];            
            str = [str sprintf('   - Y:  %12.4f [m]\n', this.mpos.Y)];            
            str = [str sprintf('   - Z:  %12.4f [m]\n\n', this.mpos.Z)];            
            
            str = [str '---- DATA SELECTION ------------------------------------------------------' 10 10];
            str = this.cc.toString(str);
            str = [str sprintf(' Using a minimum rate of [s]                       %d\n', this.p_rate)];
            str = [str sprintf(' Minimum number of satellite per epoch:            %d\n', this.min_n_sat)];
            str = [str sprintf(' Cut-off [degrees]:                                %d\n', this.cut_off)];
            str = [str sprintf(' Signal-to-noise ratio threshold [dB]:             %d\n', this.snr_thr)];
            str = [str sprintf(' Use ocean tides:                                  %d\n', this.flag_ocean)];
            str = [str sprintf(' Minimum number of epoch in an arc of observations %d\n\n', this.min_arc)];

            str = [str '---- PRE PROCESSING ------------------------------------------------------' 10 10];
            str = [str sprintf(' Enable Pre Processing                             %d\n', this.flag_pre_pro)];
            str = [str sprintf(' Cycle slip threshold [cycles]                     %g\n\n', this.cs_thr_pre_pro)];

            str = [str '---- OUTLIER DETECTION ---------------------------------------------------' 10 10];
            str = [str sprintf(' Enable Outlier detection                          %d\n', this.flag_outlier)];
            str = [str sprintf(' Threshold on code LS estimation error [m]:        %g\n', this.pp_spp_thr)];
            str = [str sprintf(' Threshold on maximum residual of code obs [m]:    %g\n', this.pp_max_code_err_thr)];
            str = [str sprintf(' Threshold on maximum residual of phase obs [m]:   %g\n\n', this.pp_max_phase_err_thr)];

            str = [str '---- PROCESSING PARAMETERS -----------------------------------------------' 10 10];
            str = this.toString@Mode_Settings(str);
            str = [str sprintf(' Compute tropospheric indicators                   %d\n\n', this.flag_tropo)];
            str = [str sprintf(' Using %s\n\n', this.W_SMODE{this.w_mode+1})];
            str = [str sprintf(' Weight function parameters (when based on SNR): \n')];
            str = [str sprintf('   - w.a:     %d\n', this.w_snr.a)];
            str = [str sprintf('   - w.zero:  %d\n', this.w_snr.zero)];
            str = [str sprintf('   - w.one:   %d\n', this.w_snr.one)];
            str = [str sprintf('   - w.A:     %d\n\n', this.w_snr.A)];            
            str = [str sprintf(' Cycle slip threshold (processing) [cycles]:       %d\n', this.cs_thr)];
            str = [str sprintf(' Enable iono free combination:                     %d\n', this.flag_ionofree)];
            str = [str sprintf(' Constrain the solution using a reference path:    %d\n', this.constrain)];
            str = [str sprintf(' Stop go stop option:                              %d\n\n', this.stop_go_stop)];

            str = [str '---- AMBIGUITY (IAR) ------------------------------------------------------' 10 10];
            str = [str sprintf(' Use ambiguity fix resolution:                     %d\n\n', this.flag_iar)];
            str = [str sprintf(' Ambiguity restart mode: %s\n\n', this.IAR_SRESTART{this.iar_restart_mode+1})];
            str = [str sprintf(' Using method %s\n\n', this.IAR_SMODE{this.iar_mode+1})];
            str = [str sprintf(' User defined fixed failure rate (methods 1,2):    %g\n', this.iar_p0)];
            str = [str sprintf(' User defined minimum success rate (for method 5): %g\n', this.iar_p0)];
            str = [str sprintf(' STD of a priori ambiguity combinations [cycles]:   %d\n\n', this.sigma0_N)];
            str = [str sprintf(' User defined threshold for ratio test:            %g\n', this.iar_mu)];
            str = [str sprintf(' Automatic determination of mu:                    %d\n', this.flag_iar_auto_mu)];
            str = [str sprintf(' Use default value for P0:                         %d\n', this.flag_iar_default_p0)];
            str = [str sprintf(' Use doppler predicted phase range:                %d\n\n', this.flag_doppler)];            
            
            str = [str '---- KALMAN FILTER PARAMETERS --------------------------------------------' 10 10];
            if this.isPP(this.p_mode)
                str = [str sprintf(' Order of the KF %s\n\n', this.DYN_SMODE_PP{this.kf_order+1})];
            else
                str = [str sprintf(' Order of the KF %s\n\n', this.DYN_SMODE_RT{this.kf_order+1})];
            end
            str = [str sprintf(' STD of initial state:                             %g\n', this.sigma0_k_pos)];
            str = [str sprintf(' STD of ENU variation:                             %g %g %g\n', struct2array(this.std_k_ENU))];
            str = [str sprintf(' STD of 3D modulus variation:                      %g\n\n', this.std_k_vel_mod)];
            str = [str sprintf(' STD of a priori tropospheric delay:               %g\n', this.sigma0_tropo)];
            str = [str sprintf(' STD of tropospheric delay:                        %g\n\n', this.std_tropo)];
            
            str = [str '---- ATMOSPHERE ----------------------------------------------------------' 10 10];
            str = [str sprintf(' Ionospheric model  %s\n', this.IONO_SMODE{this.iono_model+1})];
            str = [str sprintf(' Tropospheric model %s\n\n', this.TROPO_SMODE{this.tropo_model+1})];
            
            str = [str '---- DTM -----------------------------------------------------------------' 10 10];
            str = [str sprintf(' Use DTM:                                          %d\n', this.flag_dtm)];
            str = [str sprintf(' Folder containing DTM data:                       %s\n', this.dtm_dir)];
            str = [str sprintf(' STD of DEM model [m]:                             %g\n', this.std_dtm)];
            str = [str sprintf(' Height of the antenna above ground [m]:           %g\n\n', this.antenna_h)];
            
            str = [str '---- UI ------------------------------------------------------------------' 10 10];
            str = [str sprintf(' Plot during processing:                           %d\n', this.plot_proc)];
            str = [str sprintf(' Plot ref during processing:                       %d\n', this.plot_ref_path)];
            str = [str sprintf(' Plot sky plot during processing:                  %d\n', this.plot_skyplot_snr)];
            str = [str sprintf(' Plot error_ellipse:                               %d\n', this.plot_err_ellipse)];
            str = [str sprintf(' Plot ambiguities:                                 %d\n', this.plot_ambiguities)];
            str = [str sprintf(' Plot master station:                              %d\n', this.plot_master)];
            str = [str sprintf(' Plot on google earth:                             %d\n\n', this.plot_google_earth)];            
            
            str = [str '---- CAPTURE -------------------------------------------------------------' 10 10];
            str = [str sprintf(' Number of receivers for capturing data:           %d\n', this.c_n_receivers)];
            str = [str sprintf(' Capture rate [s]:                                 %d\n\n', this.c_rate)];
            for r = 1 : this.c_n_receivers
                str = [str sprintf(' Receiver number %d\n', r)];  %#ok<AGROW>
                str = [str sprintf('  - Protocol in use %s\n', this.C_SPROTOCOL{this.c_prtc(r)})]; %#ok<AGROW>
                str = [str sprintf('  - COM address: %s\n\n', this.c_com_addr{r})]; %#ok<AGROW>
            end
            
            str = [str '---- NTRIP ---------------------------------------------------------------' 10 10];
            str = [str sprintf(' Use NTRIP protocol:                               %d\n', this.flag_ntrip)];
            str = [str sprintf(' ip address of the server:                         %s\n', this.ntrip.ip_addr)];
            str = [str sprintf(' Port number of the server:                        %s\n', this.ntrip.port)];
            str = [str sprintf(' Mountpoint:                                       %s\n', this.ntrip.mountpoint)];
            str = [str sprintf(' Username:                                         %s\n', this.ntrip.username)];
            str = [str sprintf(' Password:                                         %s\n', this.ntrip.password)];
            str = [str sprintf(' Approximate position - lat [deg]:                 %.9g\n', this.ntrip.approx_position.lat)];
            str = [str sprintf(' Approximate position - lon [deg]:                 %.9g\n', this.ntrip.approx_position.lon)];
            str = [str sprintf(' Approximate position - h   [m]:                   %.9g\n\n', this.ntrip.approx_position.h)];
        end
        
        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the this
            % SYNTAX: s_obj.export(str_cell)
            
            if (nargin == 1)
                str_cell = {};
            end
            
            this.check(); % check before export

            str_cell = this.export@IO_Settings(str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            % RECEIVER DEFAULT PARAMETERS
            str_cell = Ini_Manager.toIniStringSection('RECEIVERS', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default STD of code observations [m]', str_cell);
            str_cell = Ini_Manager.toIniString('std_code', this.std_code, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default STD of phase observations [m]', str_cell);
            str_cell = Ini_Manager.toIniString('std_phase', this.std_phase, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default STD of iono-free phase observations [m', str_cell);
            str_cell = Ini_Manager.toIniString('std_phase_if', this.std_phase_if, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default STD of a priori receiver clock', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_clock', this.sigma0_clock, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default STD of receiver clock', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_r_clock', this.sigma0_r_clock, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Read master position from RINEX (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_rinex_mpos', this.flag_rinex_mpos, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default Master position (this values are read when not specified elsewhere)', str_cell);
            str_cell = Ini_Manager.toIniString('mpos_XYZ', struct2array(this.mpos), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            % DATA SELECTION
            str_cell = Ini_Manager.toIniStringSection('DATA_SELECTION', str_cell);
            str_cell = this.cc.export(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Processing using a minimum rate of [s]:', str_cell);
            str_cell = Ini_Manager.toIniString('p_rate', this.p_rate, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Minimum number of satellite per epoch', str_cell);
            str_cell = Ini_Manager.toIniString('min_n_sat', this.min_n_sat, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Cut-off [degrees]', str_cell);
            str_cell = Ini_Manager.toIniString('cut_off', this.cut_off, str_cell);            
            str_cell = Ini_Manager.toIniStringComment('Signal-to-noise ratio threshold [dB]', str_cell);
            str_cell = Ini_Manager.toIniString('snr_thr', this.snr_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable ocean tides modeling (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_ocean', this.flag_ocean, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Minimum length an arc (a satellite to be used must be seen for a number of consecutive epochs equal or greater than this value)', str_cell);
            str_cell = Ini_Manager.toIniString('min_arc', this.min_arc, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            % PRE PROCESSING
            str_cell = Ini_Manager.toIniStringSection('PRE_PROCESSING', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable pre-processing (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_pre_pro', this.flag_pre_pro, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Cycle slip threshold [cycles]', str_cell);
            str_cell = Ini_Manager.toIniString('cs_thr_pre_pro', this.cs_thr_pre_pro, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            % OUTLIER DETECTION
            str_cell = Ini_Manager.toIniStringSection('OUTLIER_DETECTION', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable outlier detection (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_outlier', this.flag_outlier, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Threshold on code LS estimation error [m]', str_cell);
            str_cell = Ini_Manager.toIniString('pp_spp_thr', this.pp_spp_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Threshold on maximum residual of code obs [m]', str_cell);
            str_cell = Ini_Manager.toIniString('pp_max_code_err_thr', this.pp_max_code_err_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Threshold on maximum residual of phase obs [m]', str_cell);
            str_cell = Ini_Manager.toIniString('pp_max_phase_err_thr', this.pp_max_phase_err_thr, str_cell);            
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            % PROCESSING PARAMETERS
            str_cell = Ini_Manager.toIniStringSection('PROCESSING', str_cell);
            str_cell = this.export@Mode_Settings(str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Compute tropospheric indicators (e.g. ZTD):', str_cell);
            str_cell = Ini_Manager.toIniString('flag_tropo', this.flag_tropo, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Processing using weighting mode:', str_cell);
            str_cell = Ini_Manager.toIniString('w_mode', this.w_mode, str_cell);
            for i = 1 : numel(this.W_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %d: %s', i - 1, this.W_SMODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Weight function parameters (when based on SNR): a / 0 / 1 / A', str_cell);
            str_cell = Ini_Manager.toIniString('w_snr', struct2array(this.w_snr), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Cycle slip threshold (processing) [cycles]', str_cell);
            str_cell = Ini_Manager.toIniString('cs_thr', this.cs_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable usage of iono-free combination in PPP (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_ionofree', this.flag_ionofree, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Constrain the solution using a reference path', str_cell);
            str_cell = Ini_Manager.toIniString('constrain', this.constrain, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable / Disable stop go stop mode option (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('stop_go_stop', this.stop_go_stop, str_cell);            
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            % INTEGER AMBIGUITY RESOLUTION
            str_cell = Ini_Manager.toIniStringSection('AMBIGUITY', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use integer ambiguity resolution (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_iar', this.flag_iar, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Ambiguity restart mode', str_cell);
            str_cell = Ini_Manager.toIniString('iar_restart_mode', this.iar_restart_mode, str_cell);
            for i = 1 : numel(this.IAR_SRESTART)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.IAR_SRESTART{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringComment('Ambiguity detection mode', str_cell);
            str_cell = Ini_Manager.toIniString('iar_mode', this.iar_mode, str_cell);
            for i = 1 : numel(this.IAR_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.IAR_SMODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('User defined fixed failure rate (methods 1,2) / user defined minimum success rate (for method 5)', str_cell);
            str_cell = Ini_Manager.toIniString('iar_p0', this.iar_p0, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of a priori ambiguity combinations [cycles]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_N', this.sigma0_N, str_cell);
            str_cell = Ini_Manager.toIniStringComment('User defined threshold for ratio test', str_cell);
            str_cell = Ini_Manager.toIniString('iar_mu', this.iar_mu, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Automatic determination of mu (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_iar_auto_mu', this.flag_iar_auto_mu, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use default value for P0 (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_iar_default_p0', this.flag_iar_default_p0, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use Doppler-predicted phase range for detecting cycle slips (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_doppler', this.flag_doppler, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);            
            
            % KF
            str_cell = Ini_Manager.toIniStringSection('KALMAN_FILTER', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Order of the KF', str_cell);
            str_cell = Ini_Manager.toIniString('kf_order', this.kf_order, str_cell);
            str_cell = Ini_Manager.toIniStringComment('When capture/monitor modes are in use', str_cell);
            for i = 1 : numel(this.DYN_SMODE_RT)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.DYN_SMODE_RT{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringComment('When post processing is in use:', str_cell);
            for i = 1 : numel(this.DYN_SMODE_PP)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.DYN_SMODE_PP{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of initial state [m]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_k_pos', this.sigma0_k_pos, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of ENU variation [m] / [m/s] / [m/s^2]', str_cell);
            str_cell = Ini_Manager.toIniString('std_k_ENU', struct2array(this.std_k_ENU), str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of 3D modulus variation [m] / [m/s] / [m/s^2]', str_cell);
            str_cell = Ini_Manager.toIniString('std_k_vel_mod', this.std_k_vel_mod, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of a priori tropospheric delay', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_tropo', this.sigma0_tropo, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of tropospheric delay', str_cell);
            str_cell = Ini_Manager.toIniString('std_tropo', this.std_tropo, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            % ATMOSPHERE
            str_cell = Ini_Manager.toIniStringSection('ATMOSPHERE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Ionospheric model', str_cell);
            str_cell = Ini_Manager.toIniString('iono_model', this.iono_model, str_cell);
            for i = 1 : numel(this.IONO_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.IONO_SMODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Tropospheric model', str_cell);
            str_cell = Ini_Manager.toIniString('tropo_model', this.tropo_model, str_cell);
            for i = 1 : numel(this.TROPO_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.TROPO_SMODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
                        
            % DTM
            str_cell = Ini_Manager.toIniStringSection('DTM', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use DTM (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_dtm', this.flag_dtm, str_cell);
            %str_cell = Ini_Manager.toIniStringComment('Folder containing DTM data', str_cell);
            %str_cell = Ini_Manager.toIniString('dtm_dir', this.dtm_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of DEM model [m]', str_cell);
            str_cell = Ini_Manager.toIniString('std_dtm', this.std_dtm, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Elevation of the antenna above ground [m]', str_cell);
            str_cell = Ini_Manager.toIniString('antenna_h', this.antenna_h, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            % GUI
            str_cell = Ini_Manager.toIniStringSection('UI', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Plot during processing (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('plot_proc', this.plot_proc, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Plot reference during processing (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('plot_ref_path', this.plot_ref_path, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Plot sky plot during processing (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('plot_skyplot_snr', this.plot_skyplot_snr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Plot error_ellipse (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('plot_err_ellipse', this.plot_err_ellipse, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Plot ambiguities (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('plot_ambiguities', this.plot_ambiguities, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Plot master station (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('plot_master', this.plot_master, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Plot on google earth (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('plot_google_earth', this.plot_google_earth, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            % CAPTURE
            str_cell = Ini_Manager.toIniStringSection('CAPTURE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Number of receivers for capturing data', str_cell);
            str_cell = Ini_Manager.toIniString('c_n_receivers', this.c_n_receivers, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Capture rate', str_cell);
            str_cell = Ini_Manager.toIniString('c_rate', this.c_rate, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            for r = 1 : this.c_n_receivers
                str_cell = Ini_Manager.toIniStringComment(sprintf('Protocol for receiver %d', r), str_cell);
                str_cell = Ini_Manager.toIniString(sprintf('c_prtc_%02d', r), this.c_prtc(r), str_cell);
                for i = 1 : numel(this.C_SPROTOCOL)
                    str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.C_SPROTOCOL{i}), str_cell);
                end
                str_cell = Ini_Manager.toIniStringNewLine(str_cell);
                str_cell = Ini_Manager.toIniStringComment(sprintf('COM address for receiver %d', r), str_cell);
                str_cell = Ini_Manager.toIniString(sprintf('c_com_addr_%02d', r), this.c_com_addr{r}, str_cell);                str_cell = Ini_Manager.toIniStringNewLine(str_cell);
                str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            end
            
            % NTRIP
            str_cell = Ini_Manager.toIniStringSection('NTRIP', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use NTRIP protocol (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_ntrip', this.flag_ntrip, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Here the NTRIP server parameters will follow (ip_addr, port, mountpoint, user, password, approximate_position):', str_cell);
            str_cell = Ini_Manager.toIniString('ip_addr', this.ntrip.ip_addr, str_cell);
            str_cell = Ini_Manager.toIniString('port', this.ntrip.port, str_cell);
            str_cell = Ini_Manager.toIniString('mountpoint', this.ntrip.mountpoint, str_cell);
            str_cell = Ini_Manager.toIniString('username', this.ntrip.username, str_cell);
            str_cell = Ini_Manager.toIniString('password', this.ntrip.password, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Approximate position [degrees / degrees / m]:', str_cell);
            str_cell = Ini_Manager.toIniString('ntrip_lat', this.ntrip.approx_position.lat, str_cell);
            str_cell = Ini_Manager.toIniString('ntrip_lon', this.ntrip.approx_position.lon, str_cell);
            str_cell = Ini_Manager.toIniString('ntrip_h', this.ntrip.approx_position.h, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end
    end
    
    % =========================================================================
    %  LEGACY IMPORT
    % =========================================================================
    methods (Access = 'public')
        function legacyImport(this, state)
            % import from the state variable (saved into the old interface mat file of goGPS)
            % If a group of imports fails display a warning but continue the
            % import of other groups
            % SYNTAX: this.legacyImport(state)
            
            this.legacyImport@IO_Settings(state);
            
            % RECEIVER DEFAULT PARAMETERS ---------------------------------            
            try 
                this.std_code = str2double(state.std_code);
                if (state.toggle_std_phase)
                    this.std_phase = str2double(state.std_phase);
                else
                    this.std_phase = 1e30;
                end
            catch ex
                this.logger.addWarning(['Legacy import "Receiver defaults" failed - ', ex.message])
            end
            
            % DATA SELECTION  ---------------------------------------------
            try
                this.cc.legacyImport(state);
                if (isfield(state,'srate'))
                    rates = this.UI_P_SRATE;
                    this.p_rate = rates(state.srate);
                end
                this.min_n_sat = str2double(state.min_sat);
                this.cut_off = str2double(state.cut_off);
                this.snr_thr = str2double(state.snr_thres);                
                if (isfield(state,'flag_ocean'))
                    this.flag_ocean = state.ocean;
                end
                if (isfield(state,'min_arc'))
                    this.min_arc = str2double(state.min_arc);
                end
                
            catch ex
                this.logger.addWarning(['Legacy import "Data Selection" failed - ', ex.message])
            end
            
            % PRE PROCESSING ----------------------------------------------
            try
                if (isfield(state,'pre_pro'))
                    this.flag_pre_pro = state.pre_pro;
                end
                if (isfield(state,'cs_thresh'))
                    this.cs_thr_pre_pro = str2double(state.cs_thresh);
                end
            catch ex
                this.logger.addWarning(['Legacy import "Pre Processing" failed - ', ex.message])
            end                              
            
            % OUTLIER DETECTION -------------------------------------------
            try
                if (isfield(state,'outlier'))
                    this.flag_outlier = state.outlier;
                end
                if (isfield(state,'spp_thr'))
                    this.pp_spp_thr = str2double(state.spp_thr);
                end
                if (isfield(state,'code_thr'))
                    this.pp_max_code_err_thr = str2double(state.code_thr);
                end
                if (isfield(state,'phase_thr'))
                    this.pp_max_phase_err_thr = str2double(state.phase_thr);
                end
            catch ex
                this.logger.addWarning(['Legacy import "Outlier Detection" failed - ', ex.message])
            end      
            
            % PROCESSING --------------------------------------------------
            try
                if (isfield(state,'tropo'))
                    this.flag_tropo = state.tropo;
                end
                this.legacyImport@Mode_Settings(state);
                this.constrain = state.constraint;
                this.cs_thr = str2double(state.cs_thresh);
                if (isfield(state,'wModel'))
                    this.w_mode = state.wModel - 1;
                elseif (isfield(state,'weight_0'))
                    this.w_mode =  max(4, 1 * state.weight_1 + 2 * state.weight_2 + 3 * state.weight_3 + 4 * state.weight_4);
                end  
                if (isfield(state,'obs_comb'))
                    this.flag_ionofree = (state.obs_comb == 2);
                end
                this.stop_go_stop = state.stopGOstop;
            catch ex
                this.logger.addWarning(['Legacy import "Processing" failed - ', ex.message])
            end
                        
            % MASTER STATION DEFAULT PARAMETERS ---------------------------
            try 
                this.flag_rinex_mpos = state.master_pos;   
                if (state.crs == 1)
                    X = str2double(state.master_X); if isempty(X); X = 0; end
                    Y = str2double(state.master_Y); if isempty(Y); Y = 0; end
                    Z = str2double(state.master_Z); if isempty(Z); Z = 0; end
                    this.mpos = struct('X', X, 'Y', Y, 'Z', Z);
                else
                    [X, Y, Z] = geod2cart(str2double(state.master_lat), str2double(state.master_lon), str2double(state.master_h));
                    this.mpos = struct('X', X, 'Y', Y, 'Z', Z);
                end
                
                % master_pos, crs, master_X, master_Y, master_Z,
                % master_lat, master_lon, master_h are no longer supported
                % since they are considered data and not a settings they
                % will not be included into the settings parameters
                if (state.master_pos == 0) % If the coordinates of the master are not read from RINEX
                    this.logger.addWarning(['Master position should be set into the RINEX file!', ex.message])
                    this.logger.addWarning(['the coordinates are still read from UI but soon they won''t be imported anymore', ex.message])
                end
            catch ex
                this.logger.addWarning(['Legacy import "Master position" failed - ', ex.message])
            end                
            
            % INTEGER AMBIGUITY RESOLUTION --------------------------------            
            try
                this.iar_restart_mode = state.amb_select - 1;
                this.flag_iar = state.use_lambda;
                this.iar_mode = state.lambda_method-1;
                this.iar_p0 = str2double(state.lambda_P0);
                this.iar_mu = str2double(state.lambda_mu);
                this.flag_iar_auto_mu = state.lambda_auto_mu;                
                this.flag_iar_default_p0 = state.lambda_default_P0;                
                this.flag_doppler = state.flag_doppler;
            catch ex
                this.logger.addWarning(['Legacy import "iar" failed - ', ex.message])
            end
        
            % KALMAN FILTER PARAMETERS ------------------------------------
            try                
                if this.isPP(this.p_mode)
                    interface2settings = [1 2 0 3];
                    this.kf_order = interface2settings(state.dyn_mod);
                else
                    this.kf_order = state.dyn_mod - 1;
                end
                this.sigma0_k_pos = str2double(state.std_init);
                this.std_k_ENU = struct('E', str2double(state.std_X), 'N', str2double(state.std_Y), 'U', str2double(state.std_Z));
                this.std_k_vel_mod = str2double(state.std_vel);
                
                if (state.toggle_std_dtm)
                    this.std_dtm = str2double(state.std_dtm);
                else
                    this.std_dtm = 1e30;
                end
                
                this.antenna_h = str2double(state.antenna_h);
            catch ex                
                this.logger.addWarning(['Legacy import "Kalman filter parameters" failed - ', ex.message])
            end
            
            % ATMOSPHERE ------------------------------------
            try
                if (isfield(state,'ionoModel'))
                    this.iono_model = state.ionoModel - 1;
                end
                if (isfield(state,'tropoModel'))
                    this.tropo_model = state.tropoModel - 1;
                end
            catch ex                
                this.logger.addWarning(['Legacy import "Atmosphere models" failed - ', ex.message])
            end
            
            % UI ----------------------------------------------------------            
            try
                this.plot_proc = state.plotproc;                
                this.plot_ref_path = state.ref_path;                
                this.plot_skyplot_snr = ~(state.no_skyplot_snr);
                this.plot_err_ellipse = state.err_ellipse;
                this.plot_ambiguities = state.plot_amb;
                this.plot_master = state.plot_master;
                this.plot_google_earth = state.google_earth;
            catch ex
                this.logger.addWarning(['Legacy import "UI Plotting flags" failed - ', ex.message])
            end
            
            % NTRIP -------------------------------------------------------
            try
                this.flag_ntrip = state.use_ntrip;
                this.ntrip = struct('ip_addr', state.IP_address, ...
                                    'port', state.port, ...
                                    'mountpoint', state.mountpoint, ...
                                    'username', state.username, ...
                                    'password', state.password, ...
                                    'approx_position', struct('lat', str2double(state.approx_lat), 'lon', str2double(state.approx_lon), 'h', str2double(state.approx_h)));
            catch ex
                this.logger.addWarning(['Legacy import "NTRIP parameters" failed - ', ex.message])
            end   
            
            % CAPTURE -----------------------------------------------------
            try
                this.c_n_receivers = state.num_receivers;
                rates = this.UI_C_SRATE;
                this.c_rate = rates(state.captureRate);
                for r = 1 : this.c_n_receivers
                    this.c_prtc(r) = state.(sprintf('protocol_select_%d', r));
                    this.c_com_addr{r} = state.(sprintf('com_select_%d', r));
                end
            catch ex
                this.logger.addWarning(['Legacy import "Capture parameters" failed - ', ex.message])
            end
            
            this.check(); % check after import
        end
    end    
    
    % =========================================================================
    %  ADDITIONAL PROTECTED
    % =========================================================================    
    methods (Access = 'protected')
        function postImportInit(this)
            % Operations to run after the import of new parameters
            % SYNTAX: this.postImportInit
            this.init_dtm();
        end
    end
        
    % =========================================================================
    %  ADDITIONAL PUBLIC METHODS
    % =========================================================================    
    methods (Access = 'public')
        function ini = save(this, file_path)
            % Save to a file (in INI fomat) the content of the settings object
            % SYNTAX: <ini> = this.save(<file_path>);
            %
            % when file_path is not specified settings are saved on the
            % current settings file stored whose location is stored into the
            % property "cur_ini" defined in the superclass IO_Settings
            % return optionally the ini manager object used by the save function
            
            if (nargin == 1)
                file_path = this.cur_ini;
            end
            [dir_path, ~, ~] = fileparts(file_path);
            if not(exist(dir_path, 'dir'))
                mkdir(dir_path);
            end
            this.check(); % check before saving
            ini = this.save@Settings_Interface(file_path);
        end
        
        function importIniFile(this, file_path)
            % Import from an INI file the content of the settings object
            % SYNTAX: this.importIniFile(<file_path>);
            % when file_path is not specified settings are saved on the
            % current settings file stored whose location is stored into the
            % property "cur_ini" defined in the superclass IO_Settings            if (nargin == 1)
            if (nargin == 1)
                file_path = this.cur_ini;
            end
            this.importIniFile@Settings_Interface(file_path);
        end
        
        function importLegacyFile(this, file_path)
            % Import from an INI file the content of the Settings object
            % SYNTAX: this.importIniFile(file_path);
            try
                load(file_path, 'state');
                this.legacyImport(state);
            catch ex
                this.logger.addError(sprintf('Failed to load state variable from legacy ".mat" file - %s', ex.message))
            end
        end
        
        function init_dtm(this, dtm_dir)
            % Try to load default DTM values (if flag_dtm is on)?
            % SYNTAX: s_obj.init_dtm(<dtm_dir>)

            if this.flag_dtm
                if nargin == 1
                    dtm_dir = this.dtm_dir; % from superclass IO_Settings
                end
                try
                    load([dtm_dir filesep 'tiles' filesep 'tile_header'], 'tile_header');
                    this.tile_header = tile_header;
                    this.logger.addMessage(sprintf(' - DTM tile header in %s have been read', [dtm_dir filesep 'tiles' filesep 'tile_header']));
                    load([dtm_dir filesep 'tiles' filesep 'tile_georef'], 'tile_georef');
                    this.tile_georef = tile_georef;
                    this.logger.addMessage(sprintf(' - DTM tile georef in %s have been read', [dtm_dir filesep 'tiles' filesep 'tile_georef']));
                catch
                    this.logger.addWarning(sprintf('Failed to read DTM stored in %s', [dtm_dir '/tiles/']));
                    % use default zeroes values
                end
            end
        end
    end

    % =========================================================================
    %  TEST PARAMETERS VALIDITY
    % =========================================================================    
    methods (Access = 'protected')
        
        function checkLogicalField(this, field_name)
            % Check if a logical field of the object is a valid logical number
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.checkLogicalField(string_field_name);
            this.(field_name) = this.checkLogical(field_name, this.(field_name), this.(upper(field_name)));
        end

        function checkStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.checkStringField(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            
            switch nargin
                case 2, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid);
                case 4, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid, check_existence);
                otherwise, error('Settings checkStringField called with the wrong number of parameters');
            end
        end
        
        function checkNumericField(this, field_name, limits, valid_val)
            % Check if a numeric field of the object is valid
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.checkNumericField(string_field_name, <limits>, <valid_values>);
            switch nargin
                case 2, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits);
                case 4, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits, valid_val);
                otherwise, error('Settings checkNumericField called with the wrong number of parameters');
            end
        end
    end
    
    % =========================================================================
    %  TEST PARAMETERS VALIDITY
    % =========================================================================       
    
    methods (Access = 'public')
        function check(this)    
            % Check the validity of the fields
            % SYNTAX: this.check();
            
            % RECEIVER DEFAULT PARAMETERS
            this.checkNumericField('std_code',[0 1e50]);
            this.checkNumericField('std_phase',[0 1e50]);
            this.checkNumericField('std_phase_if',[0 1e50]);
            this.checkNumericField('sigma0_clock',[0 1e50]);
            this.checkNumericField('sigma0_r_clock',[0 1e50]);            
            this.checkLogicalField('flag_rinex_mpos');
           
            this.mpos.X = this.checkNumber('mpos.X', this.mpos.X, this.MPOS.X, [-1e9 1e9]);
            this.mpos.Y = this.checkNumber('mpos.Y', this.mpos.Y, this.MPOS.Y, [-1e9 1e9]);
            this.mpos.Z = this.checkNumber('mpos.Z', this.mpos.Z, this.MPOS.Z, [-1e9 1e9]);
    
            % DATA SELECTION            
            this.checkNumericField('p_rate',[0.0001 1800]);
            this.checkNumericField('min_n_sat',[1 300]);
            this.checkNumericField('cut_off',[0 90]);
            this.checkNumericField('snr_thr',[0 70]);
            this.checkLogicalField('flag_ocean');
            this.checkNumericField('min_arc',[1 1800]);

            % PRE PROCESSING 
            this.checkLogicalField('flag_pre_pro');
            this.checkNumericField('cs_thr_pre_pro',[0 1e50]);
    
            % OUTLIER DETECTION
            this.checkLogicalField('flag_outlier');
            this.checkNumericField('pp_spp_thr',[0.001 1e50]);
            this.checkNumericField('pp_max_code_err_thr',[0.001 1e50]);
            this.checkNumericField('pp_max_phase_err_thr',[0.001 1e50]);
    
            % PROCESSING PARAMETERS
            this.checkLogicalField('flag_tropo');
            this.checkNumericField('w_mode',[0 numel(this.W_SMODE)-1]);

            this.w_snr.a = this.checkNumber('w_snr.a', this.w_snr.a, this.W_SNR.a, [0.001 100]);
            this.w_snr.zero = this.checkNumber('w_snr.zero', this.w_snr.zero, this.W_SNR.zero, [0.001 100]);
            this.w_snr.one = this.checkNumber('w_snr.one', this.w_snr.one, this.W_SNR.one, [0.001 100]);
            this.w_snr.A = this.checkNumber('w_snr.A', this.w_snr.A, this.W_SNR.A, [0.001 100]);            
            
            this.checkNumericField('cs_thr',[0 1e50]);
            this.checkLogicalField('flag_ionofree');
            this.checkLogicalField('constrain');
            this.checkLogicalField('stop_go_stop');
    
            % INTEGER AMBIGUITY RESOLUTION
            this.checkLogicalField('flag_iar');
            this.checkNumericField('iar_restart_mode',[0 numel(this.IAR_SRESTART)-1]);
            this.checkNumericField('iar_mode',[0 numel(this.IAR_SMODE)-1]);
            this.checkNumericField('iar_p0',[0.00001 1]);
            this.checkNumericField('sigma0_N',[1 100]);
            this.checkNumericField('iar_mu',[0.001 1]);
            this.checkLogicalField('flag_iar_auto_mu');
            this.checkLogicalField('flag_iar_default_p0');
            this.checkLogicalField('flag_doppler');
    
            % KF
            if this.isPP(this.p_mode)
                this.checkNumericField('kf_order',[0 numel(this.DYN_SMODE_PP)-1]);
            else
                this.checkNumericField('kf_order',[0 numel(this.DYN_SMODE_RT)-1]);
            end
        
            % RECEIVER POSITION / MOTION 
            this.checkNumericField('sigma0_k_pos',[0 1e3]);
            this.std_k_ENU.E = this.checkNumber('std_k_ENU.E', this.std_k_ENU.E, this.STD_K_ENU.E, [0 1e9]);
            this.std_k_ENU.N = this.checkNumber('std_k_ENU.N', this.std_k_ENU.N, this.STD_K_ENU.N, [0 1e9]);
            this.std_k_ENU.U = this.checkNumber('std_k_ENU.U', this.std_k_ENU.U, this.STD_K_ENU.U, [0 1e9]);
            this.checkNumericField('std_k_vel_mod',[0 1e9]);
                                
            % ATMOSPHERE
            this.checkNumericField('sigma0_tropo',[1e-6 10]);
            this.checkNumericField('std_tropo',[1e-6 1]);
            this.checkNumericField('iono_model',[0 numel(this.IONO_SMODE)-1]);
            this.checkNumericField('tropo_model',[0 numel(this.TROPO_SMODE)-1]);
       
            % DTM 
            this.checkLogicalField('flag_dtm');
            this.checkNumericField('std_dtm',[0 1e50]);
            this.checkNumericField('antenna_h',[0 1e6]);
        
            % GUI        
            this.checkLogicalField('plot_proc');
            this.checkLogicalField('plot_ref_path');
            this.checkLogicalField('plot_skyplot_snr');
            this.checkLogicalField('plot_err_ellipse');
            this.checkLogicalField('plot_ambiguities');
            this.checkLogicalField('plot_master');
            this.checkLogicalField('plot_google_earth');   
                
            % CAPTURE
            this.checkNumericField('c_n_receivers',[0 100]);
            this.checkNumericField('c_rate',[0.0001 1800]);
            this.checkNumericField('c_prtc',[1 numel(this.C_SPROTOCOL)]);

            % NTRIP        
            this.checkLogicalField('flag_ntrip');
            this.ntrip.ip_addr = this.checkString('ntrip.ip_addr', this.ntrip.ip_addr, this.NTRIP.ip_addr, true);
            this.ntrip.port = this.checkString('ntrip.port', this.ntrip.port, this.NTRIP.port, true);
            this.ntrip.mountpoint = this.checkString('ntrip.mountpoint', this.ntrip.mountpoint, this.NTRIP.mountpoint, true);
            this.ntrip.username = this.checkString('ntrip.username', this.ntrip.username, this.NTRIP.username, true);
            this.ntrip.password = this.checkString('ntrip.password', this.ntrip.password, this.NTRIP.password, true);
            this.ntrip.approx_position.lat = this.checkNumber('ntrip.approx_position.lat', this.ntrip.approx_position.lat, this.NTRIP.approx_position.lat, [-90 90]);
            this.ntrip.approx_position.lon = this.checkNumber('ntrip.approx_position.lon', this.ntrip.approx_position.lon, this.NTRIP.approx_position.lon, [-90 90]);
            this.ntrip.approx_position.h = this.checkNumber('ntrip.approx_position.h', this.ntrip.approx_position.h, this.NTRIP.approx_position.h, [-1e4 1e6]);            
            this.check@IO_Settings();
        end
    end
    
    % =========================================================================
    %  TEST
    % =========================================================================    
    methods (Static, Access = 'public')        
        function test()      
            % Test the class
            % SYNTAX: state.test()            
            s = Main_Settings();
            s.testInterfaceRoutines();
        end
    end    
end
