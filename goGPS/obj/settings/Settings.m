%   CLASS Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the Processing parameters
%
% EXAMPLE
%   settings = Settings();
%
% SEE ALSO
%   - GO_Settings     
%         it is a singleton class that store the status of the
%         settings to be used in one goGPS session
%
% FOR A LIST OF CONSTANTs and METHODS use doc Settings

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

classdef Settings < Settings_Interface & IO_Settings & Mode_Settings
    
    properties (Constant, Access = 'protected')
        % id to string of IAR modes
        IAR_MODE = {'0: ILS method with numeration in search (LAMBDA2)', ...
                    '1: ILS method with shrinking ellipsoid during search (LAMBDA3)' ...
                    '2: ILS method with numeration in search (LAMBDA3)', ...
                    '3: integer rounding method (LAMBDA3)', ...
                    '4: iar_mode = 4; % integer bootstrapping method (LAMBDA3)', ...
                    '5: Partial Ambiguity Resolution (PAR) (LAMBDA3)'}
                
        % id to string of IAR restart function
        IAR_RESTART = {'0: Observed code - phase difference', ...
                       '1: Kalman-predicted code - phase difference', ...
                       '2: Least squares adjustment'}
       
        % id to string of weight functions
        W_MODE = {'same weight for all the observations', ...
                  'weight based on satellite elevation (sin)' ...
                  'weight based on signal-to-noise ratio', ...
                  'weight based on combined elevation and signal-to-noise ratio', ...
                  'weight based on satellite elevation (exp)'}
              
        % id to string of ionospheric models
        IONO_MODE = {'0: no model', ...
                     '1: Geckle and Feen model' ...
                     '2: Klobuchar model', ...
                     '3: SBAS grid'}
                 
        % id to string of tropospheric models 
        TROPO_MODE = {'0: no model', ...
                      '1: Saastamoinen model (with standard atmosphere parameters)' ...
                      '2: Saastamoinen model (with Global Pressure Temperature model)'} 
                          
        % id to string of capturing protocols
        C_PROTOCOL = {'1: UBX (u-blox)', ...
                      '2: iTalk (Fastrax)', ...
                      '3: SkyTraq', ...
                      '4: BINR (NVS)'}
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
        std_code = 0.3; 

        % Std of phase observations [m]
        std_ph = 0.003; % (maximize to obtain a code-based solution)
        % Std of iono-free phase observations [m]
        std_ph_if = 0.003;
        
        % Std of a priori receiver clock
        sigma0_clock = 4.5e-09        
        % Std of receiver clock
        sigma0_r_clock = 1e3;
        
        % Flag to read the position of the master form the RINEX file (deprecate)
        flag_rinex_mpos = true;
        % Default master position (these are overrided when the coordinates are specified elsewhere) (deprecate)
        mpos = struct('X', 0, 'Y', 0, 'Z', 0);
        
        %------------------------------------------------------------------
        % DATA SELECTION
        %------------------------------------------------------------------
        
        % object containing info on the activated constellations
        cc =  Constellation_Collector('GRECJ');
        % Minimum processing rate [s]
        p_rate = 1;
        % Minimum number of satellites to be used in the Kalman filter
        min_n_sat = 2;        
        % Cut-off [degrees]
        cut_off = 10;
        % Signal-to-noise ratio threshold [dB]
        snr_thr = 0;                
        % Flag for enabling the usage of ocean tides modeling
        flag_ocean = true;        
        % Minimum length an arc (a satellite to be used must be seen for a number of consecutive epochs greater than this value)
        min_arc = 10;

        %------------------------------------------------------------------
        % PRE PROCESSING 
        %------------------------------------------------------------------
        
        % Flag for enabling pre-processing
        flag_pre_pro = true;        
        % Cycle slip threshold (pre-processing) [cycles]
        cs_thr_pre_pro = 1;

        %------------------------------------------------------------------
        % OUTLIER DETECTION
        %------------------------------------------------------------------

        % Flag for enabling outlier detection
        flag_outlier = true;
                
        % Threshold on the code point-positioning least squares estimation error [m]
        pp_spp_thr = 4;                 
        % Threshold on the maximum residual of code observations [m]
        pp_max_code_err_thr = 30;            
        % Threshold on the maximum residual of phase observations [m]
        pp_max_phase_err_thr = 0.2;  
        
        %------------------------------------------------------------------
        % PROCESSING PARAMETERS 
        %------------------------------------------------------------------
        
        % Flag for enabling the computation of thropospheric derived informations
        flag_tropo = true;
        
        % Parameter used to select the weightening mode for GPS observations
        w_mode = 1;
        %  - weights = 0: same weight for all the observations
        %  - weights = 1: weight based on satellite elevation (sin)
        %  - weights = 2: weight based on signal-to-noise ratio
        %  - weights = 3: weight based on combined elevation and signal-to-noise ratio
        %  - weights = 4: weight based on satellite elevation (exp)
                
        % Weight function parameters (when based on SNR)
        w_snr = struct('a', 30, 'zero', 10, 'one', 50, 'A', 30);
        
        % Cycle slip threshold (processing) [cycles]
        cs_thr = 1;

        % Flag for enabling the usage of iono-free combination
        flag_ionofree = true;

        % Constrain the solution using a reference path
        constrain = false; 
        
        % This flag add the possibility to process in stop go stop mode
        % static convergence / movement / static convergence
        stop_go_stop = false;
        
        %------------------------------------------------------------------
        % INTEGER AMBIGUITY RESOLUTION
        %------------------------------------------------------------------
        
        % Flag for enabling the automatic detection of cycle sleep
        flag_iar = 0;

        % Ambiguity restart mode
        iar_restart_mode = 2;        
        % - iar_restart_mode = 0; % Observed code - phase difference
        % - iar_restart_mode = 1; % Kalman-predicted code - phase difference
        % - iar_restart_mode = 2; % Least squares adjustment
        
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

        % Std of a priori ambiguity combinations [cycles]
        sigma0_N = 1000;

        % User defined threshold for ratio test
        iar_mu = 0.5;
        
        % Flag for enabling the automatic determination of mu
        flag_iar_auto_mu = true;
        
        % Flag for enabling the default value for P0
        flag_iar_default_P0 = true;
        
        % Flag for using doppler-predicted phase range for detecting cycle slips
        flag_doppler = false;
        
        %------------------------------------------------------------------
        % KF
        %------------------------------------------------------------------
        
        % Order of the dynamic model polynomial
        kf_order = 1;

        %------------------------------------------------------------------
        % RECEIVER POSITION / MOTION 
        %------------------------------------------------------------------
        
        % Std of initial state [m]
        sigma0_k_pos = 1;         
        % Std of ENU coordinates variation [m] / [m/s] / [m/s^2]
        std_k_ENU = struct('E', 0.5, 'N', 0.5, 'U', 0.1);
        % Std of 3D modulus variation [m] / [m/s] / [m/s^2]
        std_k_vel_mod = 0.1;
                                
        %------------------------------------------------------------------
        % ATMOSPHERE
        %------------------------------------------------------------------

        % Std of a priori tropospheric delay
        sigma0_tropo = 0.1;
        % Std of tropospheric delay
        std_tropo = 4.6e-4;
                
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
        % DTM 
        %------------------------------------------------------------------
                
        % DTM flag (use / do not use)
        flag_dtm = false;

        % Folder containing DTM files
        dtm_dir = '../data/dtm';

        % Std of DEM height [m]
        std_dtm = 0.03  % (maximize to disable DEM usage: e.g. 1e30)
        
        % Elevation of the antenna above ground
        antenna_h = 0; % when fixing 
        
        % Parameters common to all DTM tiles
        dtm_tile_header = struct('nrows', 0, 'ncols', 0, 'cellsize', 0, 'nodata', 0);
        
        % Parameters used to georeference every DTM tile
        dtm_tile_georef = zeros(1,1,4);        
        
        %------------------------------------------------------------------
        % GUI
        %------------------------------------------------------------------
        
        % plot during processing
        plot_proc = false;        
        % plot ref during processing
        plot_ref_path = false;        
        % plot sky plot during processing
        plot_skyplot_snr = false;        
        % plot error_ellipse
        plot_err_ellipse = false;        
        % plot ambiguities
        plot_ambiguities = false;        
        % plot master station
        plot_master = false;        
        % plot on google earth
        plot_google_earth = false;     
                
        %------------------------------------------------------------------
        % CAPTURE
        %------------------------------------------------------------------
        
        % Number of receiver to use for capturing data
        c_n_receivers = 1;
        
        % Capture rate [s]
        c_rate = 1;
        
        % Array with the size of c_n_receivers
        c_prtc = 1
        % - c_prtc = 1: UBX (u-blox)
        % - c_prtc = 2: iTalk (Fastrax)
        % - c_prtc = 3: SkyTraq
        % - c_prtc = 4: BINR (NVS)
        
        % Cell array with the com address of each receiver to be used
        c_com_addr = {'/dev/tty.lpss-serial1'};
        
        %------------------------------------------------------------------
        % NTRIP
        %------------------------------------------------------------------
        
        % NTRIP flag (use / do not use) 
        flag_ntrip = false;
        
        % struct containing NTRIP parameters:
        ntrip = struct('ip_addr', '127.0.0.1', ...
                       'port', '2101', ...
                       'mountpoint', '/', ...
                       'username', 'user', ...
                       'password', 'pswd', ...
                       'approx_position', struct('lat', 0, 'lon', 0, 'h', 0));        
    end
    
    % =========================================================================
    %  INIT
    % =========================================================================    
    methods
        function this = Settings()
            % Creator
            % SYNTAX: s_obj = Settings();
            this.postImportInit();
        end
    end
    
    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================    
    methods (Access = 'public')
        function import(this, settings)
            % This function import processing settings from another setting object or ini file
            % SYNTAX: s_obj.import(settings)
            
            if isa(settings, 'Ini_Manager')         
                % RECEIVER DEFAULT PARAMETERS
                this.std_code = settings.getData('std_code');
                this.std_ph = settings.getData('std_ph');
                this.std_ph_if = settings.getData('std_ph_if');
                this.sigma0_clock = settings.getData('sigma0_clock');
                this.sigma0_r_clock = settings.getData('sigma0_r_clock');

                this.flag_rinex_mpos = settings.getData('flag_rinex_mpos');                
                tmp = settings.getData('mpos_XYZ');
                this.mpos = struct('X', tmp(1), 'Y', tmp(2), 'Z', tmp(3));
                
                % DATA SELECTION
                this.cc.import(settings);
                this.p_rate = settings.getData('p_rate');
                this.min_n_sat  = settings.getData('min_n_sat');
                this.cut_off = settings.getData('cut_off');
                this.snr_thr = settings.getData('snr_thr');                
                this.flag_ocean = settings.getData('flag_ocean');
                this.min_arc = settings.getData('min_arc');
                
                % PRE PROCESSING                
                this.flag_pre_pro = settings.getData('flag_pre_pro');
                this.cs_thr_pre_pro = settings.getData('cs_thr_pre_pro');
                
                % OUTLIER DETECTION
                this.flag_outlier = settings.getData('flag_outlier');
                this.pp_spp_thr = settings.getData('pp_spp_thr');
                this.pp_max_code_err_thr = settings.getData('pp_max_code_err_thr');
                this.pp_max_phase_err_thr = settings.getData('pp_max_phase_err_thr');

                % PROCESSING PARAMETERS
                this.flag_tropo = settings.getData('flag_tropo');
                this.w_mode = settings.getData('w_mode');
                tmp = settings.getData('w_snr');
                this.w_snr = struct('a', tmp(1), 'zero', tmp(2), 'one', tmp(2), 'A', tmp(4));
                this.cs_thr = settings.getData('cs_thr');
                this.flag_ionofree = settings.getData('flag_ionofree');
                this.constrain = settings.getData('constrain');
                this.stop_go_stop = settings.getData('stop_go_stop');
                
                % INTEGER AMBIGUITY RESOLUTION
                this.flag_iar = settings.getData('flag_iar');
                this.iar_restart_mode = settings.getData('iar_restart_mode');
                this.iar_mode = settings.getData('iar_mode');
                this.iar_P0 = settings.getData('iar_P0');
                this.sigma0_N = settings.getData('sigma0_N');
                this.iar_mu = settings.getData('iar_mu');
                this.flag_iar_auto_mu = settings.getData('flag_iar_auto_mu');
                this.flag_iar_default_P0 = settings.getData('flag_iar_default_P0');
                this.flag_doppler = settings.getData('flag_doppler');
                
                % KF
                this.kf_order      = settings.getData('kf_order');

                % RECEIVER POSITION / MOTION
                this.sigma0_k_pos    = settings.getData('sigma0_k_pos');
                tmp = settings.getData('std_k_ENU');
                this.std_k_ENU = struct('E', tmp(1), 'N', tmp(2), 'U', tmp(3));
                this.std_k_vel_mod = settings.getData('std_k_vel_mod');
                                
                % ATMOSPHERE
                this.sigma0_tropo  = settings.getData('sigma0_tropo');
                this.std_tropo   = settings.getData('std_tropo');

                this.iono_model = settings.getData('iono_model');
                this.tropo_model = settings.getData('tropo_model');
                
                % DTM
                this.dtm_dir    = settings.getData('dtm_dir');
                this.flag_dtm = settings.getData('flag_dtm');
                this.std_dtm = settings.getData('std_dtm');
                this.antenna_h = settings.getData('antenna_h');
                
                % GUI
                this.plot_proc = settings.getData('plot_proc');
                this.plot_ref_path = settings.getData('plot_ref_path');
                this.plot_skyplot_snr = settings.getData('plot_skyplot_snr');
                this.plot_err_ellipse = settings.getData('plot_err_ellipse');
                this.plot_ambiguities = settings.getData('plot_ambiguities');
                this.plot_master = settings.getData('plot_master');
                this.plot_google_earth = settings.getData('plot_google_earth');
                
                % CAPTURE
                this.c_n_receivers = settings.getData('c_n_receivers');
                this.c_rate = settings.getData('c_rate');
                for r = 1 : this.c_n_receivers
                    this.c_prtc(r) = settings.getData(sprintf('c_prtc_%02d', r));
                    this.c_com_addr{r} = settings.getData(sprintf('c_com_addr_%02d', r));
                end
                
                % NTRIP
                this.flag_ntrip = settings.getData('flag_ntrip');
                this.ntrip = struct('ip_addr', settings.getData('NTRIP','ip_addr'), ...
                       'port', settings.getData('NTRIP','port'), ...
                       'mountpoint', settings.getData('NTRIP','mountpoint'), ...
                       'username', settings.getData('NTRIP','username'), ...
                       'password', settings.getData('NTRIP','password'), ...
                       'approx_position', struct('lat', settings.getData('NTRIP','ntrip_lat'), 'lon', settings.getData('NTRIP','ntrip_lon'), 'h', settings.getData('NTRIP','ntrip_h')));                
            else
                % RECEIVER DEFAULT PARAMETERS
                this.std_code = settings.std_code;
                this.std_ph = settings.std_ph;
                this.std_ph_if = settings.std_ph_if;
                this.sigma0_clock = settings.sigma0_clock;
                this.sigma0_r_clock = settings.sigma0_r_clock;
                this.flag_rinex_mpos = settings.flag_rinex_mpos;
                this.mpos = settings.mpos;
                
                % DATA SELECTION
                this.cc.import(settings.cc);
                this.p_rate = settings.p_rate;
                this.min_n_sat  = settings.min_n_sat;
                this.cut_off = settings.cut_off;
                this.snr_thr = settings.snr_thr;
                this.flag_ocean = settings.flag_ocean;
                this.min_arc = settings.min_arc;
                
                % PRE PROCESSING
                this.flag_pre_pro = settings.flag_pre_pro;
                this.cs_thr_pre_pro = settings.cs_thr_pre_pro;
                
                % OUTLIER DETECTION                
                this.flag_outlier = settings.flag_outlier;
                this.pp_spp_thr = settings.pp_spp_thr;
                this.pp_max_code_err_thr = settings.pp_max_code_err_thr;
                this.pp_max_phase_err_thr = settings.pp_max_phase_err_thr;

                % PROCESSING PARAMETERS                
                this.flag_tropo = settings.flag_tropo;
                this.w_mode = settings.w_mode;
                this.w_snr = settings.w_snr;
                this.cs_thr = settings.cs_thr;
                this.flag_ionofree = settings.flag_ionofree;
                this.constrain = settings.constrain;
                this.stop_go_stop = settings.stop_go_stop;
                
                % INTEGER AMBIGUITY RESOLUTION
                this.flag_iar = settings.flag_iar;
                this.iar_restart_mode = settings.iar_restart_mode;
                this.iar_mode = settings.iar_mode;
                this.iar_P0 = settings.iar_P0;
                this.sigma0_N = settings.sigma0_N;
                this.iar_mu = settings.iar_mu;
                this.flag_iar_auto_mu = settings.flag_iar_auto_mu;
                this.flag_iar_default_P0 = settings.flag_iar_default_P0;
                this.flag_doppler = settings.flag_doppler;                                
                
                % KF
                this.kf_order = settings.kf_order;

                % RECEIVER POSITION / MOTION 
                this.sigma0_k_pos = settings.sigma0_k_pos;
                this.std_k_ENU = settings.std_k_ENU;
                this.std_k_vel_mod = settings.std_k_vel_mod;
                       
                % ATMOSPHERE
                this.sigma0_tropo = settings.sigma0_tropo;
                this.std_tropo = settings.std_tropo;                                
                this.iono_model = settings.iono_model;
                this.tropo_model = settings.tropo_model;
                        
                % DTM 
                this.flag_dtm = settings.flag_dtm;
                this.dtm_dir    = settings.dtm_dir;
                this.std_dtm = settings.std_dtm;
                this.antenna_h = settings.antenna_h;
                
                % GUI
                this.plot_proc = settings.plot_proc;
                this.plot_ref_path = settings.plot_ref_path;
                this.plot_skyplot_snr = settings.plot_skyplot_snr;
                this.plot_err_ellipse = settings.plot_err_ellipse;
                this.plot_ambiguities = settings.plot_ambiguities;
                this.plot_master = settings.plot_master;
                this.plot_google_earth = settings.plot_google_earth;
                
                % CAPTURE
                this.c_n_receivers = settings.c_n_receivers;
                this.c_rate = settings.c_rate;
                this.c_prtc = settings.c_prtc;
                this.c_com_addr = settings.c_com_addr;
                
                % NTRIP
                this.flag_ntrip = settings.flag_ntrip;
                this.ntrip = settings.ntrip;
            end

            % Call to Super Methods
            this.import@Mode_Settings(settings);
            this.import@IO_Settings(settings);
            
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
            str = [str sprintf(' Default STD of phase observations [m]:            %g\n', this.std_ph)];
            str = [str sprintf(' Default STD of iono-free phase observations [m]:  %g\n', this.std_ph_if)];
            str = [str sprintf(' Default STD of a priori receiver clock:            %g\n', this.sigma0_clock)];
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
            str = [str sprintf(' Using %s\n\n', this.W_MODE{this.w_mode+1})];
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
            str = [str sprintf(' Ambiguity restart mode: %s\n\n', this.IAR_RESTART{this.iar_restart_mode+1})];
            str = [str sprintf(' Using method: %s\n\n', this.IAR_MODE{this.iar_mode+1})];
            str = [str sprintf(' User defined fixed failure rate (methods 1,2):    %g\n', this.iar_P0)];
            str = [str sprintf(' User defined minimum success rate (for method 5): %g\n', this.iar_P0)];
            str = [str sprintf(' STD of a priori ambiguity combinations [cycles]:   %d\n\n', this.sigma0_N)];
            str = [str sprintf(' User defined threshold for ratio test:            %g\n', this.iar_mu)];
            str = [str sprintf(' Automatic determination of mu:                    %d\n', this.flag_iar_auto_mu)];
            str = [str sprintf(' Use default value for P0:                         %d\n', this.flag_iar_default_P0)];
            str = [str sprintf(' Use doppler predicted phase range:                %d\n\n', this.flag_doppler)];            
            
            str = [str '---- KALMAN FILTER PARAMETERS --------------------------------------------' 10 10];
            str = [str sprintf(' Order of the KF:                                  %d\n\n', this.kf_order)];
            str = [str sprintf(' STD of initial state:                             %g\n', this.sigma0_k_pos)];
            str = [str sprintf(' STD of ENU variation:                             %g %g %g\n', struct2array(this.std_k_ENU))];
            str = [str sprintf(' STD of 3D modulus variation:                      %g\n\n', this.std_k_vel_mod)];
            str = [str sprintf(' STD of a priori tropospheric delay:                %g\n', this.sigma0_tropo)];
            str = [str sprintf(' STD of tropospheric delay:                        %g\n\n', this.std_tropo)];
            
            str = [str '---- ATMOSPHERE ----------------------------------------------------------' 10 10];
            str = [str sprintf(' Ionospheric model:  %s\n', this.IONO_MODE{this.iono_model+1})];
            str = [str sprintf(' Tropospheric model: %s\n\n', this.TROPO_MODE{this.tropo_model+1})];
            
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
                str = [str sprintf('  - Protocol in use %s\n', this.C_PROTOCOL{this.c_prtc(r)})]; %#ok<AGROW>
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
                        
            str_cell = this.export@IO_Settings(str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            % RECEIVER DEFAULT PARAMETERS
            str_cell = Ini_Manager.toIniStringSection('RECEIVERS', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default STD of code observations [m]', str_cell);
            str_cell = Ini_Manager.toIniString('std_code', this.std_code, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default STD of phase observations [m]', str_cell);
            str_cell = Ini_Manager.toIniString('std_ph', this.std_ph, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Default STD of iono-free phase observations [m', str_cell);
            str_cell = Ini_Manager.toIniString('std_ph_if', this.std_ph_if, str_cell);
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
            for i = 1 : numel(this.W_MODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %d: %s', i - 1, this.W_MODE{i}), str_cell);
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
            for i = 1 : numel(this.IAR_RESTART)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.IAR_RESTART{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringComment('Ambiguity detection mode', str_cell);
            str_cell = Ini_Manager.toIniString('iar_mode', this.iar_mode, str_cell);
            for i = 1 : numel(this.IAR_MODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.IAR_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('User defined fixed failure rate (methods 1,2) / user defined minimum success rate (for method 5)', str_cell);
            str_cell = Ini_Manager.toIniString('iar_P0', this.iar_P0, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of a priori ambiguity combinations [cycles]', str_cell);
            str_cell = Ini_Manager.toIniString('sigma0_N', this.sigma0_N, str_cell);
            str_cell = Ini_Manager.toIniStringComment('User defined threshold for ratio test', str_cell);
            str_cell = Ini_Manager.toIniString('iar_mu', this.iar_mu, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Automatic determination of mu (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_iar_auto_mu', this.flag_iar_auto_mu, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use default value for P0 (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_iar_default_P0', this.flag_iar_default_P0, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use Doppler-predicted phase range for detecting cycle slips (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_doppler', this.flag_doppler, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            
            
            % KF
            str_cell = Ini_Manager.toIniStringSection('KALMAN_FILTER', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Order of the KF', str_cell);
            str_cell = Ini_Manager.toIniString('kf_order', this.kf_order, str_cell);
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
            for i = 1 : numel(this.IONO_MODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.IONO_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Tropospheric model', str_cell);
            str_cell = Ini_Manager.toIniString('tropo_model', this.tropo_model, str_cell);
            for i = 1 : numel(this.TROPO_MODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.TROPO_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
                        
            % DTM
            str_cell = Ini_Manager.toIniStringSection('DTM', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Use DTM (0/1)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_dtm', this.flag_dtm, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Folder containing DTM data', str_cell);
            str_cell = Ini_Manager.toIniString('dtm_dir', this.dtm_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('STD of DEM model [m]', str_cell);
            str_cell = Ini_Manager.toIniString('std_dtm', this.std_dtm, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Elevation of the antenna above ground [m]', str_cell);
            str_cell = Ini_Manager.toIniString('antenna_h', this.std_dtm, str_cell);
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
                for i = 1 : numel(this.C_PROTOCOL)
                    str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.C_PROTOCOL{i}), str_cell);
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
            % SYNTAX: s_obj.legacyImport(state)

            
            % RECEIVER DEFAULT PARAMETERS ---------------------------------            
            try 
                this.std_code = str2double(state.std_code);
                if (state.toggle_std_phase)
                    this.std_ph = str2double(state.std_phase);
                else
                    this.std_ph = 1e30;
                end
            catch ex
                this.logger.addWarning(['Legacy import "Receiver defaults" failed - ', ex.message])
            end
            
            % DATA SELECTION  ---------------------------------------------
            try
                this.cc.legacyImport(state);
                if (isfield(state,'srate'))
                    rates = [1/10 1/5 1/2 1 5 15 30];
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
                    this.mpos = struct('X', str2double(state.master_X), 'Y', str2double(state.master_Y), 'Z', str2double(state.master_Z));
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
                this.iar_restart_mode = state.amb_select;
                this.flag_iar = state.use_lambda;
                this.iar_mode = state.lambda_method;
                this.iar_P0 = str2double(state.lambda_P0);
                this.iar_mu = str2double(state.lambda_mu);
                this.flag_iar_auto_mu = state.lambda_auto_mu;                
                this.flag_iar_default_P0 = state.lambda_default_P0;                
                this.flag_doppler = state.flag_doppler;
            catch ex
                this.logger.addWarning(['Legacy import "iar" failed - ', ex.message])
            end
        
            % KALMAN FILTER PARAMETERS ------------------------------------
            try
                this.kf_order = state.dyn_mod;
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
                this.c_rate = state.captureRate;
                for r = 1 : this.c_n_receivers      
                    this.c_prtc(r) = state.(sprintf('protocol_select_%d', r));
                    this.c_com_addr{r} = state.(sprintf('com_select_%d', r));
                end
            catch ex
                this.logger.addWarning(['Legacy import "Capture parameters" failed - ', ex.message])
            end   

        end
    end
    
    % =========================================================================
    %  Additional Protected
    % =========================================================================    
    methods (Access = 'protected')
        function postImportInit(this)
            % Operations to run after the import of new parameters
            % SYNTAX: s_obj.postImportInit
            this.init_dtm();
        end        
    end
    
    % =========================================================================
    %  Additional Public methods
    % =========================================================================    
    methods (Access = 'public')
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
    %  TEST
    % =========================================================================    
    methods (Static, Access = 'public')        
        function test()      
            % Test the class
            % SYNTAX: Settings.test()            
            s = Settings();
            s.testInterfaceRoutines();
        end
    end    
end
