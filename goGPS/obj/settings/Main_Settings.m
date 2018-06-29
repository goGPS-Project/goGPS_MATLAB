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
%   - Global_Configuration
%         it is a singleton class that store the status of the
%         Main_Settings to be used in one goGPS session
%
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 3 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, Giulio Taliaferro, ...
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

classdef Main_Settings < Settings_Interface & Command_Settings

    properties (Constant, Access = 'protected')
        % id to string of out modes
        DEFAULT_DIR_IN = ['..' filesep '..' filesep ];
        DEFAULT_DIR_OUT = ['..' filesep '..' filesep];
    end
    
    % Real constant
    properties(Constant, Access = 'private')
        % PRE PROCESSING
        CS_THR_PRE_PRO = 1;                             % Cycle slip threshold (pre-processing) [cycles]
    end
    
    properties (Constant, Access = 'public')
        % Location of the latest project (the ini contains just a reference to the default project ini file - that is actually a settings file
        LAST_SETTINGS = 'last_settings.ini';
    end

    % Default values for each field - useful to restore corrupted fields
    properties (Constant, Access = 'private')
        
        % IO ---------------------------------------------------------------------------------------------------------------------------------------------------
        % PROJECT
        PRJ_NAME = 'Default PPP project';  % Name of the project
        PRJ_HOME = [fileparts(which('goGPS.m')) filesep '..' filesep 'data' filesep 'project' filesep 'default_PPP' filesep]; % Location of the project <relative path from goGPS folder>
        CUR_INI = [Main_Settings.PRJ_HOME 'Config' filesep 'settings.ini']; % Location of the current ini file

        % SESSION
        SSS_DATE_START = GPS_Time('2015-08-23'); % Start of the processing session
        SSS_DATE_STOP = GPS_Time('2015-08-23');  % End of the processing session
        SSS_ID_LIST = '0';   % id character sequence to be use for the session $(S) special keyword
        SSS_ID_START = '0';  % first session id (char of sss_id_list)
        SSS_ID_STOP = '0';   % last session id (char of sss_id_list)
        
        SSS_FILE_BASED = true;         % is the session management file based
        SSS_DURATION = 86400;          % session duration in seconds
        SSS_BUFFER = [3600*3 3600*3]; % session overlap in seconds [left right]
        
        
        FLAG_KEEP_REC_LIST = true; % Flag to store the receivers for all the sessions
        
        % STATIONS
        OBS_DIR = 'RINEX';
        OBS_NAME = {'ZIMM${DOY}${S}.${YY}O'};
        REC_DYN_MODE = 0;    % Array for each receiver specify the kind of station
                             % - rec_mode = 0; static
                             % - rec_mode = 1; constant velocity
                             % - rec_mode = 2; constant acceleration
                             % - rec_mode = 3; variable (stop-go-stop)

        % id to string of Kalman Filter dynamic modes
        DYN_MODE = {'0: static', ...
        '1: constant velocity' ...
        '2: constant acceleration', ...
        '3: variable (stop-go-stop)'}

        REC_TARGET = 0;
        REC_MASTER = 1;
        REC_REFERENCE = 2;
        OBS_TYPE_LIST = {'Target', ...
                         'Master', ...
                         'SEID reference'};

        CRD_DIR = [Main_Settings.DEFAULT_DIR_IN 'station' filesep 'CRD' filesep]; % Path to Ephemeris files folder
        CRD_NAME = '';    % Location of the stations coordinate file
        MET_DIR = [Main_Settings.DEFAULT_DIR_IN 'station' filesep 'MET' filesep]; % Path to Clock Offset files folder
        MET_NAME = '';    % Location of the meteorological file
        OCEAN_DIR = [Main_Settings.DEFAULT_DIR_IN 'station' filesep 'ocean' filesep]; % Path to CRX folder containing files of Satellites problems
        OCEAN_NAME = '';  % Location of the ocean loading file

        % REFERENCE
        %REMOTE_RES_CONF_DIR = [Main_Settings.DEFAULT_DIR_IN filesep 'goGPSconfig' filesep];
        REMOTE_RES_CONF_DIR = '';
        ERP_DIR = [Main_Settings.DEFAULT_DIR_IN 'reference' filesep 'ERP' filesep]; % Earth Rotation Parameters
        ERP_NAME = ''; % Name of ERP files
        IGRF_DIR = [Main_Settings.DEFAULT_DIR_IN 'reference' filesep 'IGRF' filesep]; % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        IGRF_NAME = 'igrf12coeff.txt';

        GEOID_DIR = [Main_Settings.DEFAULT_DIR_IN 'reference' filesep 'geoid' filesep]; % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        GEOID_NAME = 'geoid_EGM2008_05.mat'; % File name of the Geoid containing the geoid to be used for the computation of hortometric heighs
        IONO_DIR = [Main_Settings.DEFAULT_DIR_IN 'reference' filesep 'IONO' filesep];
        IONO_NAME = '';
        ATM_LOAD_DIR = [Main_Settings.DEFAULT_DIR_IN 'reference' filesep 'ATM_LOAD' filesep];
        ATM_LOAD_NAME_NT = '';
        ATM_LOAD_NAME_T = 's1_s2_s3_cm_noib_grid.dat';
        VMF_DIR = [Main_Settings.DEFAULT_DIR_IN 'reference' filesep 'VMF' filesep '${YYYY}' filesep];
        VMF_NAME = '';

        % COMPUTATION CENTERS
        % With official products for orbits and clocks
        PREFERRED_EPH = {'final', 'rapid', 'ultra', 'broadcast'}
        PREFERRED_CENTER = {'default'}
        PREFERRED_IONO = {'final', 'predicted1', 'predicted2', 'broadcast'}
        FLAG_CHECK_REMOTE = true;

        % SATELLITES
        EPH_DIR = [Main_Settings.DEFAULT_DIR_IN 'satellite' filesep 'EPH' filesep '${WWWW}' filesep]; % Path to Ephemeris files folder
        EPH_NAME = ''; % Name for Ephemeris files
        CLK_DIR = [Main_Settings.DEFAULT_DIR_IN 'satellite' filesep 'CLK' filesep '${WWWW}' filesep]; % Path to Clock Offset files folder
        CLK_NAME = ''; % Name of Clock Offset files
        CRX_DIR = [Main_Settings.DEFAULT_DIR_IN 'satellite' filesep 'CRX' filesep]; % Path to CRX folder containing files of Satellites problems
        CRX_NAME = 'SAT_${YYYY}.CRX';
        DCB_DIR = [Main_Settings.DEFAULT_DIR_IN 'satellite' filesep 'DCB' filesep]; % Path to DCB folder containing files of Differential Code Biases
        EMS_DIR = [Main_Settings.DEFAULT_DIR_IN 'satellite' filesep 'SBAS' filesep 'EMS' filesep]; % Path to EMS folder containing files of EGNOS Message Server.

        % ANTENNA
        ATX_DIR = [Main_Settings.DEFAULT_DIR_IN 'antenna' filesep 'ATX' filesep]; % Location of the antex files
        ATX_NAME = 'igs14.atx';    % Name antex file

        % OUT PATH
        OUT_DIR = [Main_Settings.DEFAULT_DIR_OUT  'project' filesep 'default_DD' filesep 'out' filesep]; % Directory containing the output of the project
        OUT_PREFIX = 'out';  % Every time a solution is computed a folder with prefix followed by the run number is created
        RUN_COUNTER = [];     % This parameter store the current run number

        % STD PAR ----------------------------------------------------------------------------------------------------------------------------------------------
        
        % ADV RECEIVER DEFAULT PARAMETERS
        STD_CODE = 3;                                   % Std of code observations [m]
        STD_PHASE = 0.03;                               % Std of phase observations [m]
        STD_PHASE_IF = 0.009;                           % Std of iono-free phase observations [m]
        SIGMA0_CLOCK = 4.47e-9;                         % Std of a priori receiver clock
        SIGMA0_R_CLOCK = 31                             % Std of receiver clock

        % DATA SELECTION
        CC = Constellation_Collector('G');              % object containing info on the activated constellations
        MIN_N_SAT = 2;                                  % Minimum number of satellites needed to process a valid epoch
        CUT_OFF = 10;                                   % Cut-off [degrees]
        SNR_THR = 0;                                    % Signal-to-noise ratio threshold [dB]
        MIN_ARC = 10;                                   % Minimum length of an arc (a satellite to be used must be seen for a number of consecutive epochs greater than this value)

        % ADV DATA SELECTION
        FLAG_OUTLIER = true;                            % Flag for enabling outlier detection
        PP_SPP_THR = 100;                               % Threshold on the code point-positioning least squares estimation error [m]
        PP_MAX_CODE_ERR_THR = 30;                       % Threshold on the maximum residual of code observations [m]
        PP_MAX_PHASE_ERR_THR = 0.2;                     % Threshold on the maximum residual of phase observations [m]

        % PROCESSING PARAMETERS
        W_MODE = 1;                                     % Parameter used to select the weightening mode for GPS observations
                                                        %  - weights = 0: same weight for all the observations
                                                        %  - weights = 1: weight based on satellite elevation (sin)
                                                        %  - weights = 2: weight based on signal-to-noise ratio
                                                        
        FLAG_AMB_FIX = false;                           % try to fix ambiguity
        FLAG_SMOOTH_TROPO_OUT = true;                   % smooth the output paramtere at bounadries
        
        FLAG_SOLID_EARTH = true;                        % Flag to enable solid eearth tide corrections
        FLAG_POLE_TIDE = true;                          % Falg to enable pole earth tide corrections
        FLAG_PHASE_WIND = true;                         % Flag to enable pahse wrap up correction
        FLAG_SHAPIRO = true;                            % Flag to enable shapiro delay corrections
        FLAG_OCEAN_LOAD = true;                         % Flag to enable Ocean Loading Corrections
        FLAG_ATM_LOAD = false;                          % FAlg to enable Atmospheric Loading Corrections
        FLAG_HOI = false;                               % Flag to enable High Order Ionospherich effects and bendigs
        FLAG_REC_PCV = true;                            % Flag to enable receiver pcv corrections

        % ATMOSPHERE
        FLAG_TROPO = false;                             % Flag for enabling the estimation of tropospheric delay
        FLAG_TROPO_GRADIENT = false;                    % Flag for enabling the estimation of tropospheric delay gradient
        IONO_MANAGEMENT = 1;                            % Flag for enabling the usage of iono-free combination
        IONO_MODEL = 2;                                 % Ionospheric model to be used (1: none, 2: Klobuchar, 3: SBAS)
                                                        % - iono_model = 1: no model
                                                        % - iono_model = 2: Klobuchar model
                                                        % - iono_model = 3: SBAS grid
                                                        
        ZD_MODEL  = 2;                                  % A-priori Tropospheric Zenith delay model to be used (0: none, 1: Saastamoinen , 2: Vienna mapping function gridded delays)
                                                        % - zd_model = 1: no model
                                                        % - zd_model = 2: Saastamoinen
                                                        % - zd_model = 3:  Vienna mapping function gridded delay
        MAPPING_FUNCTION = 1                            % Mapping function to be used
                                                        % 1 : GMF
                                                        % 2 : VMF gridded
                                                        % ADV ATMOSPHERE
        METEO_DATA = 2;                                 % Meteo data to be used
                                                        % 1: standard atmopshere
                                                        % 2: GPT
                                                        % 3: MET File
        %SIGMA0_TROPO = 0.2;                             % Std of a priori tropospheric delay
        %SIGMA0_TROPO_GRADIENT = 1;                      % Std of a priori tropospheric gradient

        STD_TROPO = 0.0015;                             % Std of tropospheric delay [m/h]
        STD_TROPO_GRADIENT = 0.0004;                    % Std of tropospheric gradient [m/h]
    end

    properties (Constant, Access = 'public')

        % id to string of weight functions
        W_SMODE = {'same weight for all the observations', ...
                   'weight based on satellite elevation (sin)' ...
                   'weight based on satellite elevation (exp)'}

        % id to string of ionospheric models
        IONO_SMODE = {'1: no model', ...
                      '2: Klobuchar model', ...
                      '3: IONEX'}
        IONO_LABEL = {'No model', ...
                      'Klobuchar model', ...
                      'IONEX'}
        % id to string of tropospheric models
        ZD_SMODE = {'1: Saastamoinen model' ...
            '2: Vienna Mapping Function gridded'}
        ZD_LABEL =  {'Saastamoinen','VMF_GRD'}
        % id to string of mappig functions
        MF_SMODE = {'1: Global Mapping Function', ...
            '2: Vienna Mapping Function gridded'}
        MF_LABEL = {'GMF','VMF_GRD'}
        % id to string of meteo dtata
        MD_SMODE = {'1: standard atmosphere', ...
                       '2: Global Pressure Temperature Model' ...
                       '3: MET file'}
        MD_LABEL = {'Standard Atmosphere','GPT','MET file'}
        IE_SMODE = {'1: Iono free', ...
                    '2: smoothed geometry free re-applyed to observables',...
                    '3: external model'}
        IE_LABEL = {'Iono free','Smooth GF','External model'}
    end

    properties (SetAccess = protected, GetAccess = protected)
        % Location of the current ini file
        cur_ini = Main_Settings.CUR_INI;
    end

    %%  Processing parameters
    % ------------------------------
    % note that some of the processing parameters are actually properties of other objects
    % e.g. std_k_ENU depend on the receiver motion
    % However for various reasons I can decide to ignore the "correct"
    % values and redefine them. So the values are also
    % stored here as parameters used for a specific processing
    properties (SetAccess = public, GetAccess = public)
        %------------------------------------------------------------------
        % PROJECT
        %------------------------------------------------------------------

        % Name of the project
        prj_name = Main_Settings.PRJ_NAME;

        % Location of the project <relative path from goGPS folder>
        prj_home = Main_Settings.PRJ_HOME;

        %------------------------------------------------------------------
        % SESSION
        %------------------------------------------------------------------

        sss_date_start = Main_Settings.SSS_DATE_START;    % start of the processing session
        sss_date_stop =  Main_Settings.SSS_DATE_STOP;     % end of the processing session
        sss_id_list =    Main_Settings.SSS_ID_LIST;       % id character sequence to be use for the session $(S) special keyworc
        sss_id_start =   Main_Settings.SSS_ID_START;      % first session id (char of sss_id_list)
        sss_id_stop =    Main_Settings.SSS_ID_STOP;       % last session id (char of sss_id_list)
        
        sss_file_based = Main_Settings.SSS_FILE_BASED;
        sss_duration   = Main_Settings.SSS_DURATION;
        sss_buffer    = Main_Settings.SSS_BUFFER;
        
        flag_keep_rec_list = Main_Settings.FLAG_KEEP_REC_LIST; % Flag to store the receivers for all the sessions
        flag_smooth_tropo_out = Main_Settings.FLAG_SMOOTH_TROPO_OUT;

        %------------------------------------------------------------------
        % STATIONS
        %------------------------------------------------------------------

        % Observation files of the Receivers
        % reference receivers (e.g. master, SEID reference)

        obs_dir = Main_Settings.OBS_DIR;    % Directory containing the data (static)
        obs_name = Main_Settings.OBS_NAME;  % File name of the receivers (can contain special keywords)
        obs_full_name;                    % Full name of the observations generated during runtime from the provided parameters

        rec_dyn_mode = Main_Settings.REC_DYN_MODE; % Array for each receiver specify the kind of station
                                          % - rec_dyn_mode = 0; static
                                          % - rec_dyn_mode = 1; constant velocity
                                          % - rec_dyn_mode = 2; constant acceleration
                                          % - rec_dyn_mode = 3; variable (stop-go-stop)

        % Path to stations coordinates files
        crd_dir = Main_Settings.CRD_DIR;
        % Location of the stations coordinate file
        crd_name = Main_Settings.CRD_NAME;

        % Path to stations meteorological files
        met_dir = Main_Settings.MET_DIR;
        % Location of the meteorological file
        met_name =  Main_Settings.MET_NAME;
        met_full_name; % Full name of the met file generated during runtime from the provided parameters

        % Path to stations ocean loading files
        ocean_dir = Main_Settings.OCEAN_DIR;

        %------------------------------------------------------------------
        % COMPUTATION CENTERS
        %------------------------------------------------------------------
        % Centers for computation of orbits and other related parameters
        preferred_eph = Main_Settings.PREFERRED_EPH;          % kind of orbits to prefer
        preferred_iono = Main_Settings.PREFERRED_IONO;
        preferred_center = Main_Settings.PREFERRED_CENTER;
        flag_check_remote = Main_Settings.FLAG_CHECK_REMOTE;

        %------------------------------------------------------------------
        % SATELLITES
        %------------------------------------------------------------------

        eph_dir = Main_Settings.EPH_DIR;    % Path to Ephemeris files folder
        eph_name = Main_Settings.EPH_NAME;  % File name of ephemeris
        eph_full_name;                    % Full name of the ephemeris generated during runtime from the provided parameters

        clk_dir = Main_Settings.CLK_DIR;    % Path to Clock Offset files folder
        clk_name = Main_Settings.CLK_NAME;  % File name of clock offsets
        clk_full_name;                    % Full name of the clock offsets generated during runtime from the provided parameters

        % Path to CRX folder containing files of Satellites problems
        crx_dir = Main_Settings.CRX_DIR;
        crx_name = Main_Settings.CRX_NAME;
        % Path to DCB folder containing files of Differential Code Biases
        dcb_dir = Main_Settings.DCB_DIR;
        dcb_name = []; % setted in File_Wizard.ConjureDCB
        % Path to EMS folder containing files of EGNOS Message Server.
        ems_dir = Main_Settings.EMS_DIR;

        %------------------------------------------------------------------
        % ANTENNAS
        %------------------------------------------------------------------

        atx_dir = Main_Settings.ATX_DIR;    % Location of the antex file
        atx_name = Main_Settings.ATX_NAME;  % Location of the antex file

        %------------------------------------------------------------------
        % REFERENCE
        %------------------------------------------------------------------

        % Path to file containing the reference path
        erp_dir = Main_Settings.ERP_DIR;     % Path to ERP files folder
        erp_name = Main_Settings.ERP_NAME;   % File name of ERP
        erp_full_name;                     % Full name of ERPs generated during runtime from the provided parameters

        iono_dir = Main_Settings.IONO_DIR;   % Path to IONO files folder
        iono_name = Main_Settings.IONO_NAME; % Path to IONO files folder
        atm_load_dir = Main_Settings.ATM_LOAD_DIR;   % Path to IONO files folder
        atm_load_name_nt = Main_Settings.ATM_LOAD_NAME_NT; % Path to IONO files folder
        atm_load_name_t = Main_Settings.ATM_LOAD_NAME_T;
        vmf_dir = Main_Settings.VMF_DIR;   % Path to IONO files folder
        vmf_name = Main_Settings.VMF_NAME; % Path to IONO files folder
        iono_full_name;                    % Full name of ERPs generated during runtime from the provided parameters
        remote_res_conf_dir = Main_Settings.REMOTE_RES_CONF_DIR;

        igrf_dir = Main_Settings.IGRF_DIR;   % Path to IGRF files folder
        igrf_name = Main_Settings.IGRF_NAME;

        geoid_dir = Main_Settings.GEOID_DIR;   % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        geoid_name = Main_Settings.GEOID_NAME; % Name of the Geoid file containing the geoid to be used for the computation of hortometric heighs

        ocean_name =  Main_Settings.OCEAN_NAME; % Location of the ocean loading file

        %------------------------------------------------------------------
        % OUTPUT
        %------------------------------------------------------------------

        out_dir = Main_Settings.OUT_DIR;        % Directory containing the output of the project
        out_prefix = Main_Settings.OUT_PREFIX;  % Every time a solution is computed a folder with prefix followed by the run number is created
        out_full_path;                        % Full prefix of the putput files generated during runtime from the provided parameters

        % This parameter store the current run number
        run_counter = Main_Settings.RUN_COUNTER;
        run_counter_is_set = false; % When importing the run counter, check if is set -> when set overwrite output
                
        %------------------------------------------------------------------
        % ADVANCED RECEIVER DEFAULT PARAMETERS
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

        %------------------------------------------------------------------
        % DATA SELECTION
        %------------------------------------------------------------------

        % object containing info on the activated constellations
        cc =  Main_Settings.CC;
        % Minimum number of satellites needed to process a valid epoch
        min_n_sat = Main_Settings.MIN_N_SAT;
        % Cut-off [degrees]
        cut_off = Main_Settings.CUT_OFF;
        % Signal-to-noise ratio threshold [dB]
        snr_thr = Main_Settings.SNR_THR;
        % Minimum length of an arc (a satellite to be used must be seen for a number of consecutive epochs greater than this value)
        min_arc = Main_Settings.MIN_ARC;

        %------------------------------------------------------------------
        % ADV DATA SELECTION
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

        % Parameter used to select the weightening mode for GPS observations
        w_mode = Main_Settings.W_MODE;
        %  - weights = 0: same weight for all the observations
        %  - weights = 1: weight based on satellite elevation (sin)
        %  - weights = 2: weight based on signal-to-noise ratio
        %  - weights = 3: weight based on combined elevation and signal-to-noise ratio
        %  - weights = 4: weight based on satellite elevation (exp)

        flag_amb_fix = Main_Settings.FLAG_AMB_FIX;
        
        
        
        % Flag for enabling the usage of iono-free combination
        iono_management  = Main_Settings.IONO_MANAGEMENT;
        flag_solid_earth = Main_Settings.FLAG_SOLID_EARTH;
        flag_pole_tide   = Main_Settings.FLAG_POLE_TIDE;
        flag_phase_wind  = Main_Settings.FLAG_PHASE_WIND;
        flag_shapiro     = Main_Settings.FLAG_SHAPIRO;
        flag_ocean_load  = Main_Settings.FLAG_OCEAN_LOAD;
        flag_atm_load    = Main_Settings.FLAG_ATM_LOAD;
        flag_hoi         = Main_Settings.FLAG_HOI;
        flag_rec_pcv     = Main_Settings.FLAG_REC_PCV;

        %------------------------------------------------------------------
        % ATMOSPHERE
        %------------------------------------------------------------------

        % Flag for enabling the estimation of tropospheric delay
        flag_tropo = Main_Settings.FLAG_TROPO;
        % Flag for enabling the estimation of tropospheric delay
        flag_tropo_gradient = Main_Settings.FLAG_TROPO_GRADIENT;

        % Ionospheric model to be used (0: none, 1: Geckle and Feen, 2: Klobuchar, 3: SBAS)
        iono_model = Main_Settings.IONO_MODEL;
        % - iono_model = 0: no model
        % - iono_model = 1: Geckle and Feen model
        % - iono_model = 2: Klobuchar model
        % - iono_model = 3: SBAS grid

        % A-priori Tropospheric model to be used (0: none, 1: Saastamoinen std parameters, 2: Saastamoinen global pararameters)
        zd_model = Main_Settings.ZD_MODEL;
        mapping_function = Main_Settings.MAPPING_FUNCTION;
        meteo_data = Main_Settings.METEO_DATA;

        %------------------------------------------------------------------
        % ADV ATMOSPHERE
        %------------------------------------------------------------------
        % Std of a priori tropospheric delay
        %sigma0_tropo = Main_Settings.SIGMA0_TROPO;
        %sigma0_tropo_gradient = Main_Settings.SIGMA0_TROPO_GRADIENT;
        % Std of tropospheric delay [m / h]
        std_tropo = Main_Settings.STD_TROPO;
        std_tropo_gradient = Main_Settings.STD_TROPO;
    end

    % =========================================================================
    %%  INIT
    % =========================================================================
    methods
        function this = Main_Settings(ini_settings_file)
            % Creator
            %
            % SYNTAX
            %   s_obj = Main_Settings(<ini_settings_file>);

            this.initLogger();
            this.log.addMarkedMessage('Building settings object...');
            this.log.newLine();
            if (nargin == 1)
                if ~exist(ini_settings_file, 'file')
                    this.log.addWarning(sprintf('File "%s" not found!', ini_settings_file));
                    ini_settings_file = this.LAST_SETTINGS;
                end
            else
                ini_settings_file = this.LAST_SETTINGS;
            end

            if (exist(ini_settings_file, 'file') == 2)
                this.importIniFile(ini_settings_file);
            else
                this.log.addMarkedMessage('Using default settings');
                this.log.newLine();
                this.postImportInit();
            end
        end
    end

    % =========================================================================
    %%  INTERFACE REQUIREMENTS
    % =========================================================================
    methods (Access = 'public')
        function import(this, state)
            % This function import processing settings from another setting object or ini file
            %
            % SYNTAX
            %   s_obj.import(state)

            fnp = File_Name_Processor;
            
            if isa(state, 'Ini_Manager')
                % PROJECT
                this.prj_name   = fnp.checkPath(state.getData('prj_name'));
                if isempty(fnp.getPath(state.file_name))
                    dir_fallback = fnp.getFullDirPath([fileparts(which('goGPS.m')) filesep '..' filesep 'data' filesep 'project' filesep 'default_DD' filesep]);
                else
                    dir_fallback = fnp.getRelDirPath([fnp.getPath(state.file_name) filesep '..'], pwd);
                end
                if isempty(state.getData('prj_home'))
                    this.prj_home  = dir_fallback;
                else
                   this.prj_home = fnp.getFullDirPath(state.getData('prj_home'),  pwd, dir_fallback);
                end
                if ~exist(this.prj_home, 'dir')
                    this.log.addWarning(sprintf('Project home "%s" does not exist\nusing prj_home = "%s"', this.prj_home, dir_fallback));
                    this.prj_home = dir_fallback;
                end

                % SESSION
                if ~(isempty(state.getData('sss_date_start')))
                    this.sss_date_start = GPS_Time(datenum(state.getData('sss_date_start')));
                end
                if ~(isempty(state.getData('sss_date_stop')))
                    this.sss_date_stop = GPS_Time(datenum(state.getData('sss_date_stop')));
                end
                this.sss_id_list = state.getData('sss_id_list');
                this.sss_id_start = state.getData('sss_id_start');
                this.sss_id_stop = state.getData('sss_id_stop');
                
                this.sss_file_based = state.getData('sss_file_based');
                this.sss_duration   = state.getData('sss_duration');
                this.sss_buffer    = state.getData('sss_buffer');

                this.flag_keep_rec_list = state.getData('flag_keep_rec_list');
                this.flag_smooth_tropo_out  = state.getData('flag_smooth_tropo_out');
                
                % STATIONS
                this.obs_dir  = fnp.getFullDirPath(state.getData('obs_dir'), this.prj_home, pwd);
                this.obs_name = fnp.checkPath(state.getData('obs_name'));
                this.obs_full_name = {};
                
                this.rec_dyn_mode = state.getData('rec_dyn_mode');

                this.crd_dir    = fnp.getFullDirPath(state.getData('crd_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('crd_dir')), this.prj_home));
                this.crd_name   = fnp.checkPath(state.getData('crd_name'));
                this.met_dir    = fnp.getFullDirPath(state.getData('met_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('met_dir')), this.prj_home));
                this.met_name   = fnp.checkPath(state.getData('met_name'));
                this.met_full_name = {};
                this.ocean_dir  = fnp.getFullDirPath(state.getData('ocean_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('ocean_dir')), this.prj_home));
                this.ocean_name = fnp.checkPath(state.getData('ocean_name'));

                % REFERENCE
                %this.remote_res_conf_dir = fnp.getFullDirPath(settings.getData('remote_res_conf_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('remote_res_conf_dir')), this.prj_home));
                this.igrf_dir   = fnp.getFullDirPath(state.getData('igrf_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('igrf_dir')), this.prj_home));
                this.igrf_name  = fnp.checkPath(state.getData('igrf_name'));
                this.erp_dir    = fnp.getFullDirPath(state.getData('erp_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('erp_dir')), this.prj_home));
                this.geoid_dir  = fnp.getFullDirPath(state.getData('geoid_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('geoid_dir')), this.prj_home));
                this.geoid_name = fnp.checkPath(state.getData('geoid_name'));
                this.iono_dir   = fnp.getFullDirPath(state.getData('iono_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('iono_dir')), this.prj_home));
                this.atm_load_dir   = fnp.getFullDirPath(state.getData('atm_load_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('atm_load_dir')), this.prj_home));
                this.vmf_dir   = fnp.getFullDirPath(state.getData('vmf_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('vmf_dir')), this.prj_home));

                % COMPUTATION CENTERS
                this.preferred_eph = fnp.checkPath(state.getData('preferred_eph'));
                this.preferred_iono = fnp.checkPath(state.getData('preferred_iono'));
                this.preferred_center = fnp.checkPath(state.getData('preferred_center'));

                % SATELLITES
                this.eph_dir    = fnp.getFullDirPath(state.getData('eph_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('eph_dir')), this.prj_home));
                this.clk_dir    = fnp.getFullDirPath(state.getData('clk_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('clk_dir')), this.prj_home));
                this.crx_dir    = fnp.getFullDirPath(state.getData('crx_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('crx_dir')), this.prj_home));
                this.dcb_dir    = fnp.getFullDirPath(state.getData('dcb_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('dcb_dir')), this.prj_home));
                this.ems_dir    = fnp.getFullDirPath(state.getData('ems_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('ems_dir')), this.prj_home));

                % ANTENNAS
                this.atx_dir    = fnp.getFullDirPath(state.getData('atx_dir'), this.prj_home, pwd);
                this.atx_name   = fnp.checkPath(state.getData('atx_name'));

                % OUTPUT
                this.out_dir = fnp.getFullDirPath(state.getData('out_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('out_dir')), this.prj_home));
                if ~exist(this.out_dir, 'dir')
                    % fallback of fallback
                    this.out_dir = fnp.getFullDirPath(state.getData('out_dir'), this.prj_home);
                end
                this.out_prefix = fnp.checkPath(state.getData('out_prefix'));
                this.run_counter = state.getData('run_counter');
                this.run_counter_is_set = ~isempty(this.run_counter);
                
                
                
                % RECEIVER DEFAULT PARAMETERS
                this.std_code = state.getData('std_code');
                this.std_phase = state.getData('std_phase');
                this.std_phase_if = state.getData('std_phase_if');
                this.sigma0_clock = state.getData('sigma0_clock');
                this.sigma0_r_clock = state.getData('sigma0_r_clock');

                % DATA SELECTION
                this.cc.import(state);
                this.min_n_sat  = state.getData('min_n_sat');
                this.cut_off = state.getData('cut_off');
                this.snr_thr = state.getData('snr_thr');
                this.min_arc = state.getData('min_arc');

                % ADV DATA SELECTION
                this.flag_outlier = state.getData('flag_outlier');
                this.pp_spp_thr = state.getData('pp_spp_thr');
                this.pp_max_code_err_thr = state.getData('pp_max_code_err_thr');
                this.pp_max_phase_err_thr = state.getData('pp_max_phase_err_thr');

                % PROCESSING PARAMETERS
                this.w_mode = state.getData('w_mode');
                this.flag_solid_earth = state.getData('flag_solid_earth');
                this.flag_pole_tide = state.getData('flag_pole_tide');
                this.flag_phase_wind = state.getData('flag_phase_wind');
                this.flag_shapiro = state.getData('flag_shapiro');
                this.flag_ocean_load = state.getData('flag_ocean_load');
                this.flag_atm_load = state.getData('flag_atm_load');
                this.flag_hoi = state.getData('flag_hoi');
                this.flag_rec_pcv = state.getData('flag_rec_pcv');

                % ATMOSPHERE
                this.flag_tropo = state.getData('flag_tropo');
                this.flag_tropo_gradient = state.getData('flag_tropo_gradient');
                this.iono_management = state.getData('iono_management');
                this.iono_model = state.getData('iono_model');
                this.zd_model = state.getData('zd_model');
                this.mapping_function = state.getData('mapping_function');
                this.meteo_data = state.getData('meteo_data');

                % ADV ATMOSPHERE
                %this.sigma0_tropo  = state.getData('sigma0_tropo');
                %this.sigma0_tropo_gradient  = state.getData('sigma0_tropo_gradient');
                this.std_tropo   = state.getData('std_tropo');
                this.std_tropo_gradient   = state.getData('std_tropo_gradient');
            else
                % PROJECT
                this.prj_name   = state.prj_name;
                this.prj_home   = state.prj_home;
                %this.cur_ini   = settings.cur_ini;

                % SESSION
                this.sss_date_start = state.sss_date_start;
                this.sss_date_stop = state.sss_date_stop;
                this.sss_id_list = state.sss_id_list;
                this.sss_id_start = state.sss_id_start;
                this.sss_id_stop = state.sss_id_stop;
                
                this.sss_file_based = state.sss_file_based;
                this.sss_duration   = state.sss_duration;
                this.sss_buffer    = state.sss_buffer;

                this.flag_keep_rec_list = state.flag_keep_rec_list;
                this.flag_smooth_tropo_out = state.flag_smooth_tropo_out;

                % STATIONS
                this.obs_dir = state.obs_dir;
                this.obs_name = state.obs_name;
                this.obs_full_name = {};

                this.rec_dyn_mode = state.rec_dyn_mode;
                
                this.crd_dir     = state.crd_dir;
                this.crd_name    = state.crd_name;
                this.met_dir     = state.met_dir;
                this.met_name    = state.met_name;
                this.met_full_name = {};
                this.ocean_dir   = state.ocean_dir;
                this.ocean_name  = state.ocean_name;

                % REFERENCE
                %this.remote_res_conf_dir = settings.remote_res_conf_dir;
                this.igrf_dir  = state.igrf_dir;
                this.igrf_name = state.igrf_name;
                this.erp_dir    = state.erp_dir;
                this.geoid_dir  = state.geoid_dir;
                this.geoid_name = state.geoid_name;
                this.iono_dir   = state.iono_dir;
                this.atm_load_dir   = state.atm_load_dir;
                this.vmf_dir   = state.vmf_dir;

                % COMPUTATION CENTERS
                this.preferred_eph = state.preferred_eph;
                this.preferred_iono = state.preferred_iono;
                this.preferred_center = state.preferred_center;

                % SATELLITES
                this.eph_dir     = state.eph_dir;
                this.eph_name    = state.eph_name;
                this.clk_dir     = state.clk_dir;
                this.clk_name    = state.clk_name;
                this.crx_dir     = state.crx_dir;
                this.dcb_dir     = state.dcb_dir;
                this.ems_dir     = state.ems_dir;

                % ANTENNA
                this.atx_dir     = state.atx_dir;
                this.atx_name    = state.atx_name;

                % OUTPUT
                this.out_dir = state.out_dir;
                this.out_prefix = state.out_prefix;
                this.run_counter = state.run_counter;
                this.run_counter_is_set = ~isempty(this.run_counter);
                
                
                % ADV RECEIVER DEFAULT PARAMETERS
                this.std_code = state.std_code;
                this.std_phase = state.std_phase;
                this.std_phase_if = state.std_phase_if;
                this.sigma0_clock = state.sigma0_clock;
                this.sigma0_r_clock = state.sigma0_r_clock;

                % DATA SELECTION
                this.cc.import(state.cc);
                this.min_n_sat  = state.min_n_sat;
                this.cut_off = state.cut_off;
                this.snr_thr = state.snr_thr;
                this.min_arc = state.min_arc;

                this.flag_outlier = state.flag_outlier;
                this.pp_spp_thr = state.pp_spp_thr;
                this.pp_max_code_err_thr = state.pp_max_code_err_thr;
                this.pp_max_phase_err_thr = state.pp_max_phase_err_thr;

                % PROCESSING PARAMETERS
                this.w_mode = state.w_mode;
                this.flag_solid_earth = state.flag_solid_earth;
                this.flag_pole_tide = state.flag_pole_tide;
                this.flag_phase_wind = state.flag_phase_wind;
                this.flag_shapiro = state.flag_shapiro;
                this.flag_ocean_load = state.flag_ocean_load;
                this.flag_atm_load = state.flag_atm_load;
                this.flag_hoi = state.flag_hoi;
                this.flag_rec_pcv = state.flag_rec_pcv;
               
                % ATMOSPHERE
                this.flag_tropo = state.flag_tropo;
                this.flag_tropo_gradient = state.flag_tropo_gradient;
                this.iono_management = state.iono_management;
                this.iono_model = state.iono_model;
                this.zd_model = state.zd_model;
                this.mapping_function = state.mapping_function;
                this.meteo_data = state.meteo_data;

                % ADV ATMOSPHERE
                %this.sigma0_tropo = state.sigma0_tropo;
                %this.sigma0_tropo_gradient = state.sigma0_tropo_gradient;
                this.std_tropo = state.std_tropo;
                this.std_tropo_gradient = state.std_tropo_gradient;
            end

            this.check(); % check after import
            this.eph_full_name = '';
            this.clk_full_name = '';
            this.erp_full_name = '';
            this.updateObsFileName();

            % Call to Super Methods
            this.import@Command_Settings(state);

            this.postImportInit();
        end

        function str = toString(this, str)
            % Display the satellite system in use
            %
            % SYNTAX
            %   s_obj.toString(str)

            if (nargin == 1)
                str = '';
            end

            fnp = File_Name_Processor();
            
            str = [str '---- PROJECT --------------------------------------------------------------' 10 10];
            str = [str sprintf(' Project name:                                     %s\n', this.prj_name)];
            str = [str sprintf(' Project home:                                     %s\n', fnp.getRelDirPath(this.prj_home, pwd))];
            str = [str sprintf(' Path to the current project ini file:             %s\n\n', this.cur_ini)];
            str = [str '---- SESSION --------------------------------------------------------------' 10 10];
            str = [str sprintf(' Definition of the file names to be parsed\n')];
            if ~(this.sss_date_start.isempty)
                str = [str sprintf(' Session start at                              %s \n', this.sss_date_start.toString())];
                str = [str sprintf(' Session end at                                %s \n', this.sss_date_stop.toString())];
            end
            str = [str sprintf(' Character sequence to be used for the sessions    %s \n', this.sss_id_list)];
            str = [str sprintf(' First session char                                %c \n', this.sss_id_start)];
            str = [str sprintf(' Last session char                                 %c \n', this.sss_id_stop)];
            
            str = [str sprintf(' Session management file based                     %d \n', this.sss_file_based)];
            str = [str sprintf(' Duration of the session                           %d \n', this.sss_duration  )];
            str = [str sprintf(' Overlap of the sessions [left right]              %d %d \n', this.sss_buffer(1), this.sss_buffer(end))];
            
            str = [str '---- INPUT: STATIONS  -----------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of the observation files                %s \n', fnp.getRelDirPath(this.obs_dir, this.prj_home))];
            str = [str sprintf(' Name of the observation files                     %s \n', strCell2Str(this.obs_name, ', '))];
            str = [str sprintf(' Directory of coordinates file:                    %s\n\n', fnp.getRelDirPath(this.crd_dir, this.prj_home))];
            str = [str sprintf(' Name of coordinate (CRD) file:                    %s\n', this.crd_name)];
            str = [str sprintf(' Name of the observation files                     %s\n', strCell2Str(this.obs_name, ', '))];
            str = [str sprintf(' Directory of meteorological data:                 %s\n', fnp.getRelDirPath(this.met_dir, this.prj_home))];
            str = [str sprintf(' Name of meteorological (met) files:               %s\n', strCell2Str(this.met_name))];
            str = [str sprintf(' Directory of ocean loading files:                 %s\n', fnp.getRelDirPath(this.ocean_dir, this.prj_home))];
            str = [str sprintf(' Name of ocean loading file:                       %s\n\n', this.ocean_name)];
            str = [str sprintf(' Keep all the receiver objects:                    %d\n\n', this.flag_keep_rec_list)];
            str = [str sprintf(' Smooth tropospheric outputas:                     %d\n\n', this.flag_smooth_tropo_out)];
            
            str = [str '---- INPUT: REFERENCE ------------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of ERP files:                           %s\n', fnp.getRelDirPath(this.erp_dir, this.prj_home))];
            str = [str sprintf(' Name of ERP files:                                %s\n', this.erp_name)];
            str = [str sprintf(' Directory of IGRF files:                          %s\n', fnp.getRelDirPath(this.igrf_dir, this.prj_home))];
            str = [str sprintf(' Name of IGRF files:                               %s\n', this.igrf_name)];
            str = [str sprintf(' Directory of Geoid models:                        %s\n', fnp.getRelDirPath(this.geoid_dir, this.prj_home))];
            str = [str sprintf(' Name of the Geoid map file:                       %s\n', this.geoid_name)];
            str = [str sprintf(' Directory of Iono models:                         %s\n', fnp.getRelDirPath(this.iono_dir, this.prj_home))];
            str = [str sprintf(' Name of the iono mnodels/maps files:              %s\n', this.iono_name)];
            str = [str sprintf(' Directory of Atmospehric Loading models:          %s\n', fnp.getRelDirPath(this.atm_load_dir, this.prj_home))];
            str = [str sprintf(' Name of the non tidal Atmospehric Loading files:  %s\n', this.atm_load_name_nt)];
            str = [str sprintf(' Name of the tidal Atmospehric Loading files:      %s\n', this.atm_load_name_t)];
            str = [str sprintf(' Directory of Vienna Mapping Function coefficients:%s\n', fnp.getRelDirPath(this.vmf_dir, this.prj_home))];
            str = [str sprintf(' Name of Vienna Mapping Function coefficients:     %s\n\n', this.vmf_name)];
            
            str = [str '---- COMPUTATION CENTER ---------------------------------------------------' 10 10];
            str = [str sprintf(' List of server to be used for downloading ephemeris\n')];
            str = [str sprintf(' Preferred order for orbits products:              %s\n', strCell2Str(this.preferred_eph))];
            str = [str sprintf(' Preferred order for iono products:                %s\n', strCell2Str(this.preferred_iono))];
            str = [str sprintf(' Preferred center:                                 %s\n\n', strCell2Str(this.preferred_center))];
            
            str = [str '---- INPUT: SATELLITE ------------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of Ephemeris files:                     %s\n', fnp.getRelDirPath(this.eph_dir, this.prj_home))];
            str = [str sprintf(' Name of Ephemeris files:                          %s\n', this.eph_name)];
            str = [str sprintf(' Directory of Satellite clock offsets:             %s\n', fnp.getRelDirPath(this.clk_dir, this.prj_home))];
            str = [str sprintf(' Name of Satellite clock offsets:                  %s\n', this.clk_name)];
            str = [str sprintf(' Directory of CRX (satellite problems):            %s\n', fnp.getRelDirPath(this.crx_dir, this.prj_home))];
            str = [str sprintf(' Directory of DCB (Differential Code Biases):      %s\n', fnp.getRelDirPath(this.dcb_dir, this.prj_home))];
            str = [str sprintf(' Directory of EMS (EGNOS Message Server):          %s\n\n', fnp.getRelDirPath(this.ems_dir, this.prj_home))];
            
            str = [str '---- INPUT: ANTENNAS -------------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of antennas (atx) files                 %s \n', this.atx_dir)];
            str = [str sprintf(' Antenna antex (ATX) file                          %s \n\n', this.atx_name)];
            
            str = [str '---- OUTPUT SETTINGS ------------------------------------------------------' 10 10];
            str = [str sprintf(' Directory containing the output of the project:   %s\n', fnp.getRelDirPath(this.out_dir, this.prj_home))];
            str = [str sprintf(' Prefix of each run:                               %s\n', this.out_prefix)];
            str = [str sprintf(' Run counter:                                      %d\n', this.run_counter)];
            if (this.run_counter_is_set)
                str = [str sprintf(' Run counter has been set manually => overwriting output\n\n')];
            else
                str = [str sprintf(' Run counter has not been previously set \n => it will be set automatically to avoid overwriting of the oputputs\n\n')];
            end
            
            str = [str '---- ADV RECEIVERS -------------------------------------------------------' 10 10];
            str = [str sprintf(' Default STD of code observations [m]:             %g\n', this.std_code)];
            str = [str sprintf(' Default STD of phase observations [m]:            %g\n', this.std_phase)];
            str = [str sprintf(' Default STD of iono-free phase observations [m]:  %g\n', this.std_phase_if)];
            str = [str sprintf(' Default STD of a priori receiver clock:           %g\n', this.sigma0_clock)];
            str = [str sprintf(' Default STD of receiver clock:                    %g\n', this.sigma0_r_clock)];

            str = [str '---- DATA SELECTION ------------------------------------------------------' 10 10];
            str = this.cc.toString(str);
            str = [str sprintf(' Minimum number of satellite per epoch:            %d\n', this.min_n_sat)];
            str = [str sprintf(' Cut-off [degrees]:                                %d\n', this.cut_off)];
            str = [str sprintf(' Signal-to-noise ratio threshold [dB]:             %d\n', this.snr_thr)];
            str = [str sprintf(' Minimum number of epoch in an arc of observations %d\n\n', this.min_arc)];

            str = [str '---- ADV DATA SELECTION --------------------------------------------------' 10 10];
            str = [str sprintf(' Enable Outlier detection                          %d\n', this.flag_outlier)];
            str = [str sprintf(' Threshold on code LS estimation error [m]:        %g\n', this.pp_spp_thr)];
            str = [str sprintf(' Threshold on maximum residual of code obs [m]:    %g\n', this.pp_max_code_err_thr)];
            str = [str sprintf(' Threshold on maximum residual of phase obs [m]:   %g\n\n', this.pp_max_phase_err_thr)];

            str = [str '---- PROCESSING PARAMETERS -----------------------------------------------' 10 10];
            str = [str sprintf(' Using %s\n\n', this.W_SMODE{this.w_mode+1})];
            str = [str sprintf(' Enable solide earth tides corrections:            %d\n', this.flag_solid_earth)];
            str = [str sprintf(' Enable pole tide corrections:                     %d\n', this.flag_pole_tide)];
            str = [str sprintf(' Enable phase wind up corrections:                 %d\n', this.flag_phase_wind)];
            str = [str sprintf(' Enable shapiro delay corrections:                 %d\n', this.flag_shapiro)];
            str = [str sprintf(' Enable ocean loading corrections:                 %d\n', this.flag_ocean_load)];
            str = [str sprintf(' Enable atmospheric loading corrections:           %d\n', this.flag_atm_load)];
            str = [str sprintf(' Enable high order ionosphere and bending:         %d\n', this.flag_hoi)];
            str = [str sprintf(' Enable Receiver pcv/pco corrections:              %d\n\n', this.flag_rec_pcv)];

            str = [str '---- ATMOSPHERE ----------------------------------------------------------' 10 10];
            str = [str sprintf(' Estimate tropospheric delay                       %d\n', this.flag_tropo)];
            str = [str sprintf(' Estimate tropospheric delay gradient              %d\n\n', this.flag_tropo_gradient)];
            str = [str sprintf(' Ionospheric model                                 %s\n', this.IONO_SMODE{this.iono_model+1})];
            str = [str sprintf(' Iono Mangement                                    %s\n', this.IE_SMODE{this.iono_management})];
            str = [str sprintf(' A-priori Zenith model                             %s\n', this.ZD_SMODE{this.zd_model})];
            str = [str sprintf(' Tropospheric model                                %s\n', this.MF_SMODE{this.mapping_function})];
            str = [str sprintf(' Meteo data model                                  %s\n\n', this.MD_SMODE{this.meteo_data})];
            
            str = [str '---- ADV ATMOSPHERE ------------------------------------------------------' 10 10];
            %str = [str sprintf(' STD of a priori tropospheric delay:               %g\n', this.sigma0_tropo)];
            str = [str sprintf(' STD of tropospheric delay:                        %g\n', this.std_tropo)];
            %str = [str sprintf(' STD of a priori tropospheric gradient:            %g\n', this.sigma0_tropo_gradient)];
            str = [str sprintf(' STD of tropospheric gradient:                     %g\n\n', this.std_tropo_gradient)];
            str = this.toString@Command_Settings(str);
        end

        function str_cell = exportIO_project(this, str_cell)
            % Export IO project parameters as ini file syntax
            %
            % SYNTAX
            %   str_cell = exportIO_project(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            fnp = File_Name_Processor;

            % PROJECT
            str_cell = Ini_Manager.toIniStringSection('PROJECT', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of the project', str_cell);
            str_cell = Ini_Manager.toIniString('prj_name', this.prj_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Home of the project', str_cell);
            str_cell = Ini_Manager.toIniString('prj_home', fnp.getRelDirPath(this.prj_home, pwd), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('NOTES ON THE NAMING CONVENTIONS:', str_cell);
            str_cell = File_Name_Processor.toIniStringComment(str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end

        function str_cell = exportIO_session(this, str_cell)
            % Export IO session parameters as ini file syntax
            %
            % SYNTAX
            %   str_cell = exportIO_session(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            % SESSION
            str_cell = Ini_Manager.toIniStringSection('SESSION', str_cell);
            str_cell = Ini_Manager.toIniStringComment('"sss_" parameters define the session of observation, they are used to substitute special keywords in file names', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Working session - first data of observation to consider (yyyy-mm-dd <HH:MM:SS>)', str_cell);
            str_cell = Ini_Manager.toIniStringComment('mainly used to detect the name of the file to process', str_cell);
            str_cell = Ini_Manager.toIniString('sss_date_start', this.sss_date_start.toString('yyyy-mm-dd HH:MM:SS'), str_cell);
            %str_cell = Ini_Manager.toIniString('sss_date_start', this.sss_date_start.toString('yyyy-mm-dd'), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Working session - last data of observation to consider (yyyy-mm-dd <HH:MM:SS>)', str_cell);
            str_cell = Ini_Manager.toIniString('sss_date_stop', this.sss_date_stop.toString('yyyy-mm-dd HH:MM:SS'), str_cell);
            %str_cell = Ini_Manager.toIniString('sss_date_stop', this.sss_date_stop.toString('yyyy-mm-dd'), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Id character sequence to be use for the session $(S) special keyword (e.g. "01233456789ABCabc")', str_cell);
            str_cell = Ini_Manager.toIniString('sss_id_list', this.sss_id_list, str_cell);
            str_cell = Ini_Manager.toIniStringComment('First session id (char of sss_id_list)', str_cell);
            str_cell = Ini_Manager.toIniString('sss_id_start', this.sss_id_start, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Last session id (char of sss_id_list)', str_cell);
            str_cell = Ini_Manager.toIniString('sss_id_stop', this.sss_id_stop, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to base the sessions on the RINEX files', str_cell);
            str_cell = Ini_Manager.toIniString('sss_file_based',this.sss_file_based,str_cell);
            str_cell = Ini_Manager.toIniStringComment('Session duration in seconds', str_cell);
            str_cell = Ini_Manager.toIniString('sss_duration',this.sss_duration  ,str_cell);
            str_cell = Ini_Manager.toIniStringComment('Session buffer in second [left right]', str_cell);
            str_cell = Ini_Manager.toIniString('sss_buffer',this.sss_buffer, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag DEBUG (0/1) to keep in memory all the processed receiver', str_cell);
            str_cell = Ini_Manager.toIniStringComment('WARNING: When the sessions are long do not use this feature', str_cell);
            str_cell = Ini_Manager.toIniStringComment('         this flag could cause memory problems', str_cell);
            str_cell = Ini_Manager.toIniString('flag_keep_rec_list', this.flag_keep_rec_list, str_cell);
            str_cell = Ini_Manager.toIniString('flag_smooth_tropo_out', this.flag_smooth_tropo_out, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end

        function str_cell = exportIO_station(this, str_cell)
            % Export IO stations parameters as ini file syntax
            %
            % SYNTAX
            %   str_cell = exportIO_station(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            fnp = File_Name_Processor;

            % STATION
            str_cell = Ini_Manager.toIniStringSection('STATION', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory containing the data (static)', str_cell);
            str_cell = Ini_Manager.toIniString('obs_dir', fnp.getRelDirPath(this.obs_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('File name of the receivers (can contain special keywords)', str_cell);
            str_cell = Ini_Manager.toIniString('obs_name', this.obs_name, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            str_cell = Ini_Manager.toIniStringComment('Directory of coordinates files', str_cell);
            str_cell = Ini_Manager.toIniString('crd_dir', fnp.getRelDirPath(this.crd_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of coordinates (CRD) file', str_cell);
            str_cell = Ini_Manager.toIniString('crd_name', this.crd_name, str_cell);
            
            str_cell = Ini_Manager.toIniStringComment('Set the a-priori information on the motion of the receiver', str_cell);
            str_cell = Ini_Manager.toIniString('rec_dyn_mode', this.rec_dyn_mode, str_cell);
            for i = 1 : numel(this.DYN_MODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.DYN_MODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of meteorological data', str_cell);
            str_cell = Ini_Manager.toIniString('met_dir', fnp.getRelDirPath(this.met_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Meteorological file (when found it will be used)', str_cell);
            str_cell = Ini_Manager.toIniString('met_name', this.met_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of ocean loading files', str_cell);
            str_cell = Ini_Manager.toIniString('ocean_dir', fnp.getRelDirPath(this.ocean_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of ocean loading file (when found it will be used)', str_cell);
            str_cell = Ini_Manager.toIniString('ocean_name', this.ocean_name, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end

        function str_cell = exportIO_reference(this, str_cell)
            % Export IO reference parameters as ini file syntax
            %
            % SYNTAX
            %   str_cell = exportIO_reference(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            fnp = File_Name_Processor;

            % REFERENCE
            str_cell = Ini_Manager.toIniStringSection('REFERENCE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Earth rotation/orientation parameters (ERP) files', str_cell);
            str_cell = Ini_Manager.toIniString('erp_dir', fnp.getRelDirPath(this.erp_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of International Geomagnetic Reference Frame (IGRF) files', str_cell);
            str_cell = Ini_Manager.toIniString('igrf_dir', fnp.getRelDirPath(this.igrf_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of IGRF file', str_cell);
            str_cell = Ini_Manager.toIniString('igrf_name', this.igrf_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Geoid files', str_cell);
            str_cell = Ini_Manager.toIniString('geoid_dir', fnp.getRelDirPath(this.geoid_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Filename in Geoid dir containing the map of ondulation of the geoid', str_cell);
            str_cell = Ini_Manager.toIniString('geoid_name', this.geoid_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Ionospheric Models files', str_cell);
            str_cell = Ini_Manager.toIniString('iono_dir', fnp.getRelDirPath(this.iono_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Atmospheric Loading Models files', str_cell);
            str_cell = Ini_Manager.toIniString('atm_load_dir', fnp.getRelDirPath(this.atm_load_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of VMF Coeficents', str_cell);
            str_cell = Ini_Manager.toIniString('vmf_dir', fnp.getRelDirPath(this.vmf_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end

        function str_cell = exportIO_computation_center(this, str_cell)
            % Export IO computation center parameters as ini file syntax
            %
            % SYNTAX
            %   str_cell = exportIO_computation_center(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            fnp = File_Name_Processor;

            % COMPUTATION CENTER
            str_cell = Ini_Manager.toIniStringSection('COMPUTATION_CENTER', str_cell);
            str_cell = Ini_Manager.toIniStringComment('List of the computeation center to be used for ephemeris retrival', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Every product is searched locally, when not found is downloaded', str_cell);
            str_cell = Ini_Manager.toIniStringComment('When the file is not found, the system fall back on the next available', str_cell);
            str_cell = Ini_Manager.toIniStringComment('The config file "remote_resource.ini" of the products is stored in:', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' => "%s"', fnp.getRelDirPath(this.remote_res_conf_dir, this.prj_home)), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred ephemeris type, valid only for source "igs",', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_EPH)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_eph', this.preferred_eph, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred ionospheric type,', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_IONO)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_iono', this.preferred_iono, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred center type (e.g. default, igs_glo, igs_gps, code, code_mgex, gfz, jaxa', str_cell);
            str_cell = Ini_Manager.toIniString('preferred_center', this.preferred_center, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end

        function str_cell = exportIO_satellite(this, str_cell)
            % Export IO satellite parameters as ini file syntax
            %
            % SYNTAX
            %   str_cell = exportIO_satellite(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            fnp = File_Name_Processor;

            % SATELLITES
            str_cell = Ini_Manager.toIniStringSection('SATELLITE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Ephemeris files', str_cell);
            str_cell = Ini_Manager.toIniString('eph_dir', fnp.getRelDirPath(this.eph_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of clock offset files', str_cell);
            str_cell = Ini_Manager.toIniString('clk_dir', fnp.getRelDirPath(this.clk_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of CRX files (containing satellite problems)', str_cell);
            str_cell = Ini_Manager.toIniString('crx_dir', fnp.getRelDirPath(this.crx_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DCB files (Differential Code Biases)', str_cell);
            str_cell = Ini_Manager.toIniString('dcb_dir', fnp.getRelDirPath(this.dcb_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of EMS files (EGNOS Message Server).', str_cell);
            str_cell = Ini_Manager.toIniString('ems_dir', fnp.getRelDirPath(this.ems_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end

        function str_cell = exportIO_antenna(this, str_cell)
            % Export IO antenna parameters as ini file syntax
            %
            % SYNTAX
            %   str_cell = exportIO_antenna(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            fnp = File_Name_Processor;

            % ANTENNA
            str_cell = Ini_Manager.toIniStringSection('ANTENNA', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of PCO - PCV antex (ATX) files', str_cell);
            str_cell = Ini_Manager.toIniString('atx_dir', fnp.getRelDirPath(this.atx_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('PCO - PCV antex (ATX) file', str_cell);
            str_cell = Ini_Manager.toIniString('atx_name', this.atx_name, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end

        function str_cell = exportIO_output(this, str_cell)
            % Export IO output parameters as ini file syntax
            %
            % SYNTAX
            %   str_cell = exportIO_output(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            fnp = File_Name_Processor;

            % OUTPUT
            str_cell = Ini_Manager.toIniStringSection('OUTPUT', str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Base dir that is going to store the ouput data files',str_cell);
            str_cell = Ini_Manager.toIniString('out_dir', fnp.getRelDirPath(this.out_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Prefix ("name") to add to the output (can contain special keywords / subfolders)',str_cell);
            str_cell = Ini_Manager.toIniString('out_prefix', this.out_prefix, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Current run number, when empty it will be automatically updated to avoid overwrite', str_cell);
            str_cell = Ini_Manager.toIniStringComment('the run_counter value is added as a 3 digit number to the output file name (after the prefix)', str_cell);
            str_cell = Ini_Manager.toIniStringComment('WARNING: when set it will be used, and can cause overwrites', str_cell);
            str_cell = Ini_Manager.toIniString('run_counter', iif(this.run_counter_is_set, this.run_counter, []), str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end

        function str_cell = exportIO(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the obj
            %
            % SYNTAX
            %   str_cell = exportIO(this, str_cell)
            if (nargin == 1)
                str_cell = {};
            end

            str_cell = this.exportIO_project(str_cell);
            str_cell = this.exportIO_session(str_cell);
            str_cell = this.exportIO_station(str_cell);
            str_cell = this.exportIO_reference(str_cell);
            str_cell = this.exportIO_computation_center(str_cell);
            str_cell = this.exportIO_satellite(str_cell);
            str_cell = this.exportIO_antenna(str_cell);
            str_cell = this.exportIO_output(str_cell);
        end
        
        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the this
            %
            % SYNTAX
            %   s_obj.export(str_cell)

            if (nargin == 1)
                str_cell = {};
            end

            this.check(); % check before export
            str_cell = Ini_Manager.toIniStringSection('SOFTWARE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('goGPS config file', str_cell);
            str_cell = Ini_Manager.toIniString('version', Core.GO_GPS_VERSION, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            str_cell = this.exportIO(str_cell);
            
            % RECEIVER DEFAULT PARAMETERS
            str_cell = Ini_Manager.toIniStringSection('ADV RECEIVERS', str_cell);
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
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            % DATA SELECTION
            str_cell = Ini_Manager.toIniStringSection('DATA_SELECTION', str_cell);
            str_cell = this.cc.export(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Minimum number of satellite per epoch', str_cell);
            str_cell = Ini_Manager.toIniString('min_n_sat', this.min_n_sat, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Cut-off [degrees]', str_cell);
            str_cell = Ini_Manager.toIniString('cut_off', this.cut_off, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Signal-to-noise ratio threshold [dB]', str_cell);
            str_cell = Ini_Manager.toIniString('snr_thr', this.snr_thr, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Minimum length an arc (a satellite to be used must be seen for a number of consecutive epochs equal or greater than this value)', str_cell);
            str_cell = Ini_Manager.toIniString('min_arc', this.min_arc, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

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
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Processing using weighting mode:', str_cell);
            str_cell = Ini_Manager.toIniString('w_mode', this.w_mode, str_cell);
            for i = 1 : numel(this.W_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %d: %s', i - 1, this.W_SMODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            str_cell = Ini_Manager.toIniStringComment('Enable corrections', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable solid earth tide corrections', str_cell);
            str_cell = Ini_Manager.toIniString('flag_solid_earth', this.flag_solid_earth, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable pole tide corrections', str_cell);
            str_cell = Ini_Manager.toIniString('flag_pole_tide', this.flag_pole_tide, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable phase wind up corrections', str_cell);
            str_cell = Ini_Manager.toIniString('flag_phase_wind', this.flag_phase_wind, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable Shapiro delay corrections', str_cell);
            str_cell = Ini_Manager.toIniString('flag_shapiro', this.flag_shapiro, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable ocean loading corrections', str_cell);
            str_cell = Ini_Manager.toIniString('flag_ocean_load', this.flag_ocean_load, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable atmospheric loading corrections', str_cell);
            str_cell = Ini_Manager.toIniString('flag_atm_load', this.flag_atm_load, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable high order ionospheric and bending corrections', str_cell);
            str_cell = Ini_Manager.toIniString('flag_hoi', this.flag_hoi, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Enable receiver pcv corrections', str_cell);
            str_cell = Ini_Manager.toIniString('flag_rec_pcv', this.flag_rec_pcv, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            % ATMOSPHERE
            str_cell = Ini_Manager.toIniStringSection('ATMOSPHERE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Management of ionosphere', str_cell);
            str_cell = Ini_Manager.toIniString('iono_management', this.iono_management, str_cell);
             for i = 1 : numel(this.IE_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.IE_SMODE{i}), str_cell);
             end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Ionospheric model', str_cell);
            str_cell = Ini_Manager.toIniString('iono_model', this.iono_model, str_cell);
             for i = 1 : numel(this.IONO_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.IONO_SMODE{i}), str_cell);
            end
            
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);           
            str_cell = Ini_Manager.toIniStringComment('Compute tropospheric indicators (e.g. ZTD):', str_cell);
            str_cell = Ini_Manager.toIniString('flag_tropo', this.flag_tropo, str_cell);
            str_cell = Ini_Manager.toIniString('flag_tropo_gradient', this.flag_tropo_gradient, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('A-priori zenith delay model', str_cell);
            str_cell = Ini_Manager.toIniString('zd_model', this.zd_model, str_cell);
            for i = 1 : numel(this.ZD_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.ZD_SMODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Mapping function', str_cell);
            str_cell = Ini_Manager.toIniString('mapping_function', this.mapping_function, str_cell);
            for i = 1 : numel(this.MF_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.MF_SMODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Meteo data', str_cell);
            str_cell = Ini_Manager.toIniString('meteo_data', this.meteo_data, str_cell);
            for i = 1 : numel(this.MD_SMODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.MD_SMODE{i}), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % ATMOSPHERE
            str_cell = Ini_Manager.toIniStringSection('ADV ATMOSPHERE', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(sprintf('Standard deviation of a priori tropospheric delay (default = %.3f)', this.SIGMA0_TROPO), str_cell);
            %str_cell = Ini_Manager.toIniString('sigma0_tropo', this.sigma0_tropo, str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('Standard deviation of tropospheric delay [m/h] (default = %.4f)', this.STD_TROPO), str_cell);
            str_cell = Ini_Manager.toIniString('std_tropo', this.std_tropo, str_cell);
            %str_cell = Ini_Manager.toIniStringComment(sprintf('Standard deviation of a priori tropospheric gradient (default = %.3f)', this.SIGMA0_TROPO_GRADIENT), str_cell);
            %str_cell = Ini_Manager.toIniString('sigma0_tropo_gradient', this.sigma0_tropo_gradient, str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('Standard deviation of tropospheric gradient [m/h] (default = %.4f)', this.STD_TROPO_GRADIENT), str_cell);
            str_cell = Ini_Manager.toIniString('std_tropo_gradient', this.std_tropo_gradient, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = this.export@Command_Settings(str_cell);
        end
    end

    % =========================================================================
    %%  ADDITIONAL PROTECTED
    % =========================================================================
    methods (Access = 'protected')
        function postImportInit(this)
            % Operations to run after the import of new parameters
            %
            % SYNTAX
            %   this.postImportInit
            this.check(); % check after import
        end
    end

    % =========================================================================
    %%  ADDITIONAL PUBLIC METHODS
    % =========================================================================
    methods (Access = 'public')
        function ini = save(this, file_path)
            % Save to a file (in INI fomat) the content of the settings object
            %
            % SYNTAX
            %   <ini> = this.save(<file_path>);
            %
            % when file_path is not specified settings are saved on the
            % current settings file stored whose location is stored into the
            % property "cur_ini" defined in the superclass Main_Settings
            % return optionally the ini manager object used by the save function

            if (nargin == 1)
                file_path = this.cur_ini;
            end
            [dir_path, ~, ~] = fileparts(file_path);
            if ~isempty(dir_path) && ~exist(dir_path, 'dir')
                mkdir(dir_path);
            end
            this.check(); % check before saving
            this.setFilePath(file_path);
            ini = this.save@Settings_Interface(file_path);
        end

        function importIniFile(this, file_path)
            % Import from an INI file the content of the settings object
            %
            % SYNTAX
            %   this.importIniFile(<file_path>);
            % when file_path is not specified settings are saved on the
            % current settings file stored whose location is stored into the
            % property "cur_ini" defined in the superclass Main_Settings            if (nargin == 1)
            if (nargin == 1)
                file_path = this.cur_ini;
            end
            this.setFilePath(file_path);
            if (exist(file_path, 'file') == 2)
                this.importIniFile@Settings_Interface(file_path);
                this.postImportInit();
                this.log.addStatusOk(sprintf('File "%s" found, settings imported!', file_path));
            else
                this.log.addWarning(sprintf('File "%s" not found, settings not imported!', file_path));
            end
        end

        function importLegacyFile(this, file_path)
            % Import from an old mat settings file the content of the Settings object
            %
            % SYNTAX
            %   this.importLegacyFile(file_path);
            try
                load(file_path, 'state');
                this.legacyImport(state);
            catch ex
                this.log.addError(sprintf('Failed to load state variable from legacy ".mat" file - %s', ex.message))
            end
            this.updateExternals();
            this.postImportInit();
        end
    end

    % =========================================================================
    %%  CHECKING FUNCTIONS
    % =========================================================================
    methods (Access = 'protected')
        function checkLogicalField(this, field_name)
            % Check if a logical field of the object is a valid logical number
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % SYNTAX
            %   this.checkLogicalField(string_field_name);
            this.(field_name) = this.checkLogical(field_name, this.(field_name), this.(upper(field_name)));
        end

        function checkCellStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % SYNTAX
            %   this.Field(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            switch nargin
                case 2, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid);
                case 4, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid, check_existence);
                otherwise, error('Settings checkCellStringField called with the wrong number of parameters');
            end
        end

        function is_existing = checkStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % check_existence 0: do not check
            %                 1: check and do nothing
            %                 2: check and try to correct
            %
            % SYNTAX
            %
            %   this.checkStringField(string_field_name, <empty_is_valid == false>, <check_existence == false>);            
            switch nargin
                case 2, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), this.(upper(field_name)));
                case 3, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid);
                case 4, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid, check_existence);
                otherwise, error('Settings checkStringField called with the wrong number of parameters');
            end
        end

        function is_existing = checkPathField(this, field_name, empty_is_valid, check_existence)
            % Check if a string path field of the object is a valid path
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % SYNTAX
            %   this.checkPathField(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            fnp = File_Name_Processor();
            this.(field_name) = fnp.getFullDirPath(this.(field_name), this.prj_home, [], fnp.getFullDirPath(this.(upper(field_name))));
            switch nargin
                case 2, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), fnp.getFullDirPath(this.(upper(field_name)), this.prj_home));
                case 3, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), fnp.getFullDirPath(this.(upper(field_name)), this.prj_home), empty_is_valid);
                case 4, [this.(field_name), is_existing] = this.checkString(field_name, this.(field_name), fnp.getFullDirPath(this.(upper(field_name)), this.prj_home), empty_is_valid, check_existence);
                otherwise, error('Settings checkStringField called with the wrong number of parameters');
            end
        end

        function checkNumericField(this, field_name, limits, valid_val)
            % Check if a numeric field of the object is valid
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            %
            % SYNTAX
            %   this.checkNumericField(string_field_name, <limits>, <valid_values>);
            switch nargin
                case 2, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits);
                case 4, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits, valid_val);
                otherwise, error('Settings checkNumericField called with the wrong number of parameters');
            end
        end
        
        % File type specific ----------------------------------------------

        function file_path = checkCrdPath(this, file_path)
            % Check if the crd file exists, if not try to look for it into the default dirs
            %
            % SYNTAX
            %   file_path = this.checkCrdPath(<file_path>)
            fnp = File_Name_Processor();
            file_path = fnp.checkPath(file_path);
            if ~isempty(file_path) && ~exist(file_path, 'file')
                [~, name, ext] = fileparts(file_path);
                % check for existence in the local project folder standard location
                tmp_path = fnp.checkPath([this.prj_home this.CRD_DIR(length(Main_Settings.DEFAULT_DIR_IN)+1:end) filesep name ext]);
                if exist(tmp_path, 'file')
                    file_path = tmp_path;
                else
                    % check for existence in the data folder standard location
                    tmp_path = fnp.checkPath([this.crd_dir filesep name ext]);
                    if exist(tmp_path, 'file')
                        file_path = tmp_path;
                    else
                        % the file cannot be found
                    end
                end
            end
            if strcmp(file_path, filesep)
                file_path = '';
            end
        end

        function file_path = checkAtxPath(this, file_path)
            % Check if the atx file exists, if not try to look for it into the default dirs
            %
            % SYNTAX
            %   file_path = this.checkAtxPath(<file_path>)
            fnp = File_Name_Processor();
            file_path = fnp.checkPath(file_path);

            if ~isempty(file_path) && ~exist(file_path, 'file')
                [~, name, ext] = fileparts(file_path);
                % check for existence in the local project folder standard location
                tmp_path = fnp.checkPath([this.prj_home this.ATX_DIR(length(Main_Settings.DEFAULT_DIR_IN)+1:end) filesep name ext]);
                if exist(tmp_path, 'file')
                    file_path = tmp_path;
                else
                    % check for existence in the data folder standard location
                    tmp_path = fnp.checkPath([this.atx_dir filesep name ext]);
                    if exist(tmp_path, 'file')
                        file_path = tmp_path;
                    else
                        % the file cannot be found
                    end
                end
            end
            if strcmp(file_path, filesep)
                file_path = '';
            end
        end

        function file_path = checkOceanPath(this, file_path)
            % Check if the atx file exists, if not try to look for it into the default dirs
            %
            % SYNTAX
            %   file_path = this.checkAtxPath(<file_path>)
            fnp = File_Name_Processor();
            file_path = fnp.checkPath(file_path);

            if ~isempty(file_path) && ~exist(file_path, 'file')
                fnp = File_Name_Processor();
                [~, name, ext] = fileparts(file_path);
                % check for existence in the local project folder standard location
                tmp_path = fnp.checkPath([this.prj_home this.OCEAN_DIR(length(Main_Settings.DEFAULT_DIR_IN)+1:end) filesep name ext]);
                if exist(tmp_path, 'file')
                    file_path = tmp_path;
                else
                    % check for existence in the data folder standard location
                    tmp_path = fnp.checkPath([this.ocean_dir filesep name ext]);
                    if exist(tmp_path, 'file')
                        file_path = tmp_path;
                    else
                        % the file cannot be found
                    end
                end
            end
            if strcmp(file_path, filesep)
                file_path = '';
            end
        end

        function file_path = checkMetPath(this, file_path)
            % Check if the met file exists, if not try to look for it into the default dirs
            %
            % SYNTAX
            %   file_path = this.checkMetPath(<file_path>)
            fnp = File_Name_Processor();
            file_path = fnp.checkPath(file_path);

            if ~isempty(file_path) && ~exist(file_path, 'file')
                fnp = File_Name_Processor();
                [~, name, ext] = fileparts(file_path);
                % check for existence in the local project folder standard location
                tmp_path = File_Name_Processor.checkPath([this.prj_home this.MET_DIR(length(Main_Settings.DEFAULT_DIR_IN)+1:end) filesep name ext]);
                if exist(tmp_path, 'file')
                    file_path = tmp_path;
                else
                    % check for existence in the data folder standard location
                    tmp_path = fnp.checkPath([this.met_dir filesep name ext]);
                    if exist(tmp_path, 'file')
                        file_path = tmp_path;
                    else
                        % the file cannot be found
                    end
                end
            end
            if strcmp(file_path, filesep)
                file_path = '';
            end
        end
    end

    % =========================================================================
    %%  CHECKING FUNCTIONS IO
    % =========================================================================
    methods (Access = 'public')
        function checkIO(this)
            % Check the validity of the fields
            %
            % SYNTAX
            %   this.check();
            EMPTY_IS_VALID = true;
            EMPTY_IS_NOT_VALID = false;
            CHECK_EXISTENCE = 1;

            this.checkStringField('prj_name', EMPTY_IS_NOT_VALID);
            is_existing = this.checkStringField('prj_home', EMPTY_IS_NOT_VALID, CHECK_EXISTENCE);
            CHECK_EXISTENCE = iif(is_existing, 2, 1);
            this.checkStringField('cur_ini', EMPTY_IS_NOT_VALID);

            this.checkCellStringField('preferred_eph', EMPTY_IS_NOT_VALID);
            this.checkCellStringField('preferred_iono', EMPTY_IS_NOT_VALID);
            this.checkCellStringField('preferred_center', EMPTY_IS_NOT_VALID);

            this.checkStringField('sss_id_list', EMPTY_IS_NOT_VALID);
            this.checkStringField('sss_id_start', EMPTY_IS_NOT_VALID);
            this.checkStringField('sss_id_stop', EMPTY_IS_NOT_VALID);

            this.checkLogicalField('sss_file_based');
            this.checkNumericField('sss_duration', [0 365*86400]);
            this.checkNumericField('sss_buffer', [0 86400*10]);
        
            this.checkPathField('obs_dir', EMPTY_IS_NOT_VALID, CHECK_EXISTENCE);
            this.checkCellStringField('obs_name', EMPTY_IS_NOT_VALID);

            this.checkPathField('atx_dir', EMPTY_IS_NOT_VALID);
            this.checkStringField('atx_name', EMPTY_IS_NOT_VALID);

            this.checkPathField('eph_dir', EMPTY_IS_NOT_VALID);
            % When the ephemeris file inserted here is not found -> the automatic downloader will dowload the proper file
            this.checkStringField('eph_name', EMPTY_IS_VALID);
            this.checkPathField('clk_dir', EMPTY_IS_NOT_VALID);
            this.checkStringField('clk_name', EMPTY_IS_VALID);
            this.checkPathField('erp_dir', EMPTY_IS_NOT_VALID);
            this.checkStringField('erp_name', EMPTY_IS_VALID);
            this.checkPathField('crx_dir', EMPTY_IS_NOT_VALID);
            this.checkPathField('dcb_dir', EMPTY_IS_NOT_VALID);
            this.checkPathField('ems_dir', EMPTY_IS_VALID);

            this.checkNumericField('rec_dyn_mode', [0 numel(this.DYN_MODE)-1]);
            if numel(this.rec_dyn_mode) < this.getRecCount()
                if numel(this.rec_dyn_mode) == 0
                    this.rec_dyn_mode(end : this.getRecCount()) = 0; % set all static
                else
                    this.rec_dyn_mode(end : this.getRecCount()) = this.rec_dyn_mode(1);
                end
            elseif numel(this.rec_dyn_mode) > this.getRecCount() % cut if longer than receiver number
                this.rec_dyn_mode = this.rec_dyn_mode(1 : this.getRecCount());
            end
            
            this.checkLogicalField('flag_keep_rec_list');
            this.checkLogicalField('flag_smooth_tropo_out');

            this.checkPathField('crd_dir', EMPTY_IS_NOT_VALID);
            this.checkPathField('met_dir', EMPTY_IS_NOT_VALID);
            this.checkStringField('ocean_dir', EMPTY_IS_NOT_VALID);

            this.checkPathField('igrf_dir', EMPTY_IS_NOT_VALID);
            this.checkStringField('igrf_name', EMPTY_IS_NOT_VALID);
            this.checkPathField('geoid_dir', EMPTY_IS_NOT_VALID, CHECK_EXISTENCE);
            this.checkStringField('geoid_name', EMPTY_IS_NOT_VALID);

            this.checkStringField('out_prefix', EMPTY_IS_VALID);

            if (this.run_counter_is_set) || ~(isempty(this.run_counter))
                this.checkNumericField('run_counter',[0 1e6]);
            end
        end

        function status = checkRecFiles(this, go_verbose)
            % check the availability of all the rinex files
            %
            % SYNTAX
            %   status = this.checkReceiverFiles(obs_type, go_verbose)
            % status is an array containing the file status for each receiver
            %
            %   0: all file are present
            %  -1: no file are present
            %   1: at least one file is present but not all

            if nargin == 2
                go_verbose = false;
            end

            n_rec = this.getRecCount();
            file_name_all = this.getRecPath();

            fnp = File_Name_Processor();

            % If no receiver have been found
            if n_rec == 0
                status = -1;
            else
                status = 0;
                for r = 1 : n_rec

                    if go_verbose
                        this.log.addMessage(sprintf('Checking files for receiver %d', r));
                    end
                    file_name = file_name_all{r};
                    file_count = 0;
                    for f = 1 : numel(file_name)
                        full_path = fnp.checkPath(file_name{f});
                        file_ok = exist(full_path, 'file') == 2;
                        file_count = file_count + uint16(logical(file_ok));
                        if go_verbose
                            if file_ok
                                this.log.addStatusOk(sprintf('%s is present', full_path));
                            else
                                this.log.addError(sprintf('%s does not exist', full_path));
                            end
                        end
                    end
                    if (file_count == 0)
                        status = -1;
                    end
                end
            end
        end

        function eph_ok = checkNavEphFiles(this)
            % check whether or not all the ephemeris files are available
            %
            % SYNTAX
            %   eph_ok = checkNavEphFiles(this)
            eph_ok = true;

            state = Global_Configuration.getCurrentSettings();
            file_name = this.getFullNavEphPath();
            file_name_rel = File_Name_Processor.getRelDirPath(file_name, state.getHomeDir());

            if isempty(file_name)
                eph_ok = false;
            elseif isempty(file_name{1})
                eph_ok = false;
            else
                this.log.addMarkedMessage(sprintf('Checking navigational files from %s', state.getHomeDir()));
                this.log.newLine();
                i = 0;
                while (i < numel(file_name) && eph_ok)
                    i = i + 1;
                    eph_ok = exist(file_name{i}, 'file') == 2;
                    if eph_ok
                        this.log.addStatusOk(sprintf('%s', file_name_rel{i}));
                    else
                        if ~(exist(file_name{i}, 'file') == 7) % if it's not a folder
                            this.log.addWarning(sprintf('%s does not exist', file_name{i}));
                        else
                            this.log.addWarning(sprintf('%s it''s a folder, no file name have been declared', file_name{i}));
                        end
                    end
                end
                this.log.newLine();
            end
        end

        function clk_ok = checkNavClkFiles(this)
            % check whether or not all the navigational clock files are available
            %
            % SYNTAX
            %   clk_ok = checkNavClkFiles(this)

            clk_ok = true;
            state = Global_Configuration.getCurrentSettings();
            file_name = this.getFullNavClkPath();
            file_name_rel = File_Name_Processor.getRelDirPath(file_name, state.getHomeDir());

            if isempty(file_name)
                clk_ok = true;
            elseif isempty(file_name{1})
                clk_ok = true;
            else
                this.log.addMarkedMessage(sprintf('Checking clock offsets files from %s', state.getHomeDir()));
                this.log.newLine();
                i = 0;
                while (i < numel(file_name) && clk_ok)
                    i = i + 1;
                    clk_ok = exist(file_name{i}, 'file') == 2;
                    if clk_ok
                        this.log.addStatusOk(sprintf('%s', file_name_rel{i}));
                    else
                        if ~(exist(file_name{i}, 'file') == 7) % if it's not a folder
                            this.log.addWarning(sprintf('%s does not exist', file_name{i}));
                        else
                            this.log.addWarning(sprintf('%s it''s a folder, no file name have been declared', file_name{i}));
                        end
                    end
                end
                this.log.newLine();
            end
        end

        function erp_ok = checkErpFiles(this)
            % check whether or not all the ERP files are available
            %
            % SYNTAX
            %   erp_ok = checkErpFiles(this)

            erp_ok = true;
            state = Global_Configuration.getCurrentSettings();
            file_name = this.getFullErpPath();
            file_name_rel = File_Name_Processor.getRelDirPath(file_name, state.getHomeDir());

            if isempty(file_name)
                erp_ok = true;
            elseif isempty(file_name{1})
                erp_ok = true;
            else
                this.log.addMarkedMessage(sprintf('Checking Earth rotation parameters files from %s', state.getHomeDir()));
                this.log.newLine();
                i = 0;
                while (i < numel(file_name) && erp_ok)
                    i = i + 1;
                    erp_ok = exist(file_name{i}, 'file') == 2;
                    if erp_ok
                        this.log.addStatusOk(sprintf('%s', file_name_rel{i}));
                    else
                        if ~(exist(file_name{i}, 'file') == 7) % if it's not a folder
                            this.log.addWarning(sprintf('%s does not exist', file_name{i}));
                        else
                            this.log.addWarning(sprintf('%s it''s a folder, no file name have been declared', file_name{i}));
                        end
                    end
                end
                this.log.newLine();
            end
        end
    end
    
    % =========================================================================
    %%  TEST PARAMETERS VALIDITY
    % =========================================================================

    methods (Access = 'public')
        function check(this)
            % Check the validity of the fields
            %
            % SYNTAX
            %   this.check();

            this.checkIO();
            
            % ADV RECEIVER DEFAULT PARAMETERS
            this.checkNumericField('std_code',[0 1e50]);
            this.checkNumericField('std_phase',[0 1e50]);
            this.checkNumericField('std_phase_if',[0 1e50]);
            this.checkNumericField('sigma0_clock',[0 1e50]);
            this.checkNumericField('sigma0_r_clock',[0 1e50]);

            % DATA SELECTION
            this.checkNumericField('min_n_sat',[1 300]);
            this.checkNumericField('cut_off',[0 90]);
            this.checkNumericField('snr_thr',[0 70]);
            this.checkNumericField('min_arc',[1 1800]);

            % ADV DATA SELECTION
            this.checkLogicalField('flag_outlier');
            this.checkNumericField('pp_spp_thr',[0.001 1e50]);
            this.checkNumericField('pp_max_code_err_thr',[0.001 1e50]);
            this.checkNumericField('pp_max_phase_err_thr',[0.001 1e50]);

            % PROCESSING PARAMETERS
            this.checkNumericField('w_mode',[1 numel(this.W_SMODE)]);
            
            this.checkLogicalField('flag_smooth_tropo_out');
            if this.isRinexSession && this.isSmoothTropoOut
                this.setSmoothTropoOut(false)
                this.log.addWarning('Smoothing of tropposphere is not possible when RINEX based sessions are requested');
            end
                
            this.checkNumericField('iono_management');
            this.checkLogicalField('flag_solid_earth');
            this.checkLogicalField('flag_pole_tide');
            this.checkLogicalField('flag_phase_wind');
            this.checkLogicalField('flag_shapiro');
            this.checkLogicalField('flag_ocean_load');
            this.checkLogicalField('flag_atm_load');
            this.checkLogicalField('flag_hoi');
            this.checkLogicalField('flag_rec_pcv');

            % ATMOSPHERE
            this.checkNumericField('iono_model',[1 numel(this.IONO_SMODE)]);
            this.checkNumericField('zd_model',[1 numel(this.ZD_SMODE)]);
            this.checkNumericField('mapping_function',[1 numel(this.MF_SMODE)]);
            this.checkNumericField('meteo_data',[1 numel(this.MD_SMODE)]);
            this.checkLogicalField('flag_tropo');
            this.checkLogicalField('flag_tropo_gradient');
            if this.flag_tropo_gradient && ~this.flag_tropo
                this.flag_tropo_gradient = false;
                this.log.addWarning('Gradients estimation appears to be requested\nbut troposphere estimation is disabled\n=> disabling gradients estimation');
            end
            
            

            % ADV ATMOSPHERE
            % this.checkNumericField('sigma0_tropo',[1e-11 10]);
            this.checkNumericField('std_tropo',[1e-12 1e50]);
            % this.checkNumericField('sigma0_tropo_gradient',[1e-11 10]);
            this.checkNumericField('std_tropo_gradient',[1e-12 1e50]);
        end
        
        function n_missing = checkDir(this, field_name, field_text, flag_verbose)
            % Check the validity of the fields
            %
            % SYNTAX
            %   n_missing = this.checkDir(field_name, field_text, flag_verbose);
            
            if nargin < 4
                flag_verbose = true;
            end
            n_missing = checkPath(this, field_name, field_text, flag_verbose, false);
        end
                
        function n_missing = checkDirErr(this, field_name, field_text, flag_verbose)
            % Check the validity of the fields
            %
            % SYNTAX
            %   n_missing = this.checkDir(field_name, field_text, flag_verbose);
            
            if nargin < 4
                flag_verbose = true;
            end
            n_missing = checkPath(this, field_name, field_text, flag_verbose, false, true);
        end
        
        function n_missing = checkFile(this, field_name, field_text, flag_verbose)
            % Check the validity of the fields
            %
            % SYNTAX
            %   n_missing = this.checkFile(field_name, field_text, flag_verbose);
            %   n_missing = this.checkFile({field_dir, field_name}, field_text, flag_verbose);
            
            if nargin < 4
                flag_verbose = true;
            end
            n_missing = checkPath(this, field_name, field_text, flag_verbose, true);
        end
        
        function n_missing = checkPath(this, field_name, field_text, flag_verbose, is_file, flag_error)
            % Check the validity of the fields
            %
            % SYNTAX 
            %   n_missing = this.checkPath(field_name, field_text, flag_verbose);
            %   n_missing = this.checkPath({field_dir, field_name}, field_text, flag_verbose);
            
            if nargin < 4
                flag_verbose = true;
            end
            if nargin < 6
                flag_error = false;
            end
            if ~iscell(field_name)
                field_name = {field_name};                
            end
            
            fnp = File_Name_Processor();
            date_start = this.getSessionsStartExt;
            date_stop = this.getSessionsStartExt;
            if numel(field_name) == 2
                file_name = this.(field_name{end});
                if ~iscell(file_name)
                    file_name = {file_name};
                end
                
                file_path = {};
                for i = 1 : numel(file_name)
                    file_path(i) = fnp.dateKeyRepBatch(fnp.checkPath(strcat(File_Name_Processor.getFullDirPath(this.(field_name{1}), this.getHomeDir), filesep, file_name{i})), date_start,  date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop); %#ok<AGROW>
                end
            else
                file_path = File_Name_Processor.getFullDirPath(this.(field_name{end}), this.getHomeDir);
                file_path = fnp.dateKeyRepBatch(file_path, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
            end
                        
            if ~iscell(file_path)
                file_path = {file_path};
            end
            
            n_missing = -numel(file_path);
            for d = 1 : numel(file_path)
                if exist(file_path{d}, 'file') == iif(is_file, 2, 7)
                    n_missing = n_missing + 1;
                end
            end
            
            if n_missing == 0
                if flag_verbose
                    this.log.addStatusOk([field_text ' is present']);
                end
            else
                if flag_verbose
                    if flag_error
                        this.log.addError(sprintf('%s is missing (%d)',field_text, -n_missing));
                    else
                        this.log.addWarning(sprintf('%s is missing (%d)',field_text, -n_missing));
                    end
                end
            end
        end
    end

    % =========================================================================
    %%  GETTERS IO
    % =========================================================================
    methods        
        function file_dir = getHomeDir(this)
            % Get the base directory containing the project
            %
            % SYNTAX
            %   file_dir = getHomeDir(this)
            file_dir = this.prj_home;
        end

        function remote_source_file = getRemoteSourceFile(this)
            % get remote source files
            %
            % SYNTAX
            %   remote_source_file = getRemoteSourceFile(this)
            fnp = File_Name_Processor();
            remote_ini_path = [this.remote_res_conf_dir iif(isempty(this.remote_res_conf_dir), '', filesep) 'remote_resource.ini'];
            remote_ini_path = fnp.getFullDirPath(remote_ini_path, this.getHomeDir);
            remote_source_file = fnp.checkPath(remote_ini_path);
            if ~exist(remote_source_file, 'file')
                remote_source_file = which('remote_resource.ini');
            end
        end       
        
        function remote_center = getRemoteCenter(this)
            %get remote center
            %
            % SYNTAX
            %   remote_center = getRemoteCenter(this)
            remote_center = this.preferred_center;
            if ~iscell(remote_center)
                remote_center = {remote_center};
            end
        end

        function preferred_eph = getPreferredEph(this)
            % get preferred ephemeris
            %
            % SYNTAX
            %   preferred_eph = getPreferredEph(this)
            preferred_eph = this.preferred_eph;
            if ~iscell(preferred_eph)
                preferred_eph = {preferred_eph};
            end
        end

        function dir = getFileDir(this, filename)
            % get file dir
            %
            % SYNTAX
            %   dir = getFileDir(this, filename)
            if length(filename) < 1
                dir = '';
                return
            end
            [~, name,ext] = fileparts(filename);
            if strcmpi(ext,'.sp3') || strcmp(ext,'.eph') || strcmp(ext,'.pre')
                dir = this.getNavEphDir();
            elseif strcmpi(ext,'.erp')
                dir = this.getErpDir();
            elseif instr(lower(ext),'.clk')
                dir = this.getNavClkDir();
            elseif strcmpi(ext,'.CRX')

            elseif ~isempty(regexp(ext,'\.\d\di', 'once')) || strcmpi(ext,'.${YY}i')
                dir = this.getIonoDir();
            elseif strcmpi(ext,'.DCB') || (strcmpi(ext,'.BSX')) || (strcmpi(ext,'.BIA'))
                dir = this.getDcbDir();
            elseif strcmpi(ext,'.${YY}p') || strcmpi(ext,'.${YY}n') || strcmpi(ext,'.${YY}l') || ~isempty(regexpi(ext,'\.\d\dp')) || ~isempty(regexpi(ext,'\.\d\dn')) || ~isempty(regexp(ext,'\.\d\dl', 'once'))
                dir = this.getNavEphDir();
            elseif strcmpi(ext,'.apl')
                dir = this.getAtmLoadDir();
            elseif instr(name,'VMFG') && instr(ext,'.H')
                 dir = this.getVMFDir();
            end

        end

        function base_rinex_dir = getRinexBaseDir(this)
            % Get the base directory containing RINEX files
            %
            % SYNTAX
            %   base_rinex_dir = getRinexBaseDir(this)
            base_rinex_dir = this.obs_dir();
        end

        function num_session = getSessionCount(this)
            % Get the number of sessions
            %
            % SYNTAX
            %   num_receiver = getSessionCount(this)
            if this.isRinexSession()
                file_name = this.getRecPath();
                num_session = numel(file_name{1});
            else
                num_session = ceil((this.getSessionsStop - this.getSessionsStart) / this.sss_duration);
            end
        end

        function num_stations = getRecCount(this)
            % Get the number of the receivers
            %
            % SYNTAX
            %   num_stations = getRecCount(this)
            num_stations = numel(this.getRecPath());
        end

        function file_name = getRecPath(this, rec_num, session)
            % Get the file list of receivers files
            %
            % SYNTAX
            %   file_name = this.getRecPath()
            %
            % A cell for each receiver containing the list of names as cell
            if isempty(this.obs_full_name) || isempty(this.obs_full_name{1})
                this.updateObsFileName();
            end
            file_name = this.obs_full_name;
            if ~iscell(file_name)
                file_name = {{file_name}};
            end
            if ~iscell(file_name{1})
                file_name{1} = {file_name};
            end
            if nargin >= 2
                file_name = file_name{rec_num};
                if nargin == 3
                    file_name = file_name{session};
                end
            end
        end
        
        function dyn_mode = getDynMode(this, rec_num)
            % Get the a-priori information on the motion of the receiver
            %
            % SYNTAX
            %   dyn_mode = getDynMode(this, rec_num)
            dyn_mode = this.rec_dyn_mode(rec_num);
        end

        function out = getNavEphType(this)
            % Get the order of preference of orbits files to search for
            %
            % SYNTAX
            %   out = getNavEphType(this)
            out = this.preferred_eph;
        end

        function out = getNavErpType(this)
            % Get the order of preference of erp files to search for
            %
            % SYNTAX
            %   out = getNavErpType(this)
            out = this.preferred_erp;
        end

        function file_name = getFullNavEphPath(this, id)
            % Get the file list of ephemeris files
            %
            % SYNTAX
            %   file_name = this.getFullNavEphPath(id)
            if isempty(this.eph_full_name)
                this.updateNavFileName();
            end
            file_name = this.eph_full_name;
            if (nargin == 2)
                file_name = file_name{id};
            end
        end

        function file_name = getFullNavClkPath(this, id)
            % Get the file list of clock files
            %
            % SYNTAX
            %   file_name = this.getFullNavClkPath(id)
            if isempty(this.clk_full_name)
                this.updateNavFileName();
            end
            file_name = this.clk_full_name;
            if (nargin == 2)
                file_name = file_name{id};
            end
        end

        function file_name = getFullErpPath(this, id)
            % Get the file list of ephemeris files
            %
            % SYNTAX
            %   file_name = this.getFullErpPath(id)
            if isempty(this.erp_full_name)
                this.updateErpFileName();
            end
            file_name = this.erp_full_name;
            if (nargin == 2)
                file_name = file_name{id};
            end
        end

        function file_name = getFullIonoPath(this, id)
            % Get the file list of ephemeris files
            %
            % SYNTAX
            %   file_name = this.getFullErpPath(id)
            if isempty(this.erp_full_name)
                this.updateErpFileName();
            end
            file_name = this.erp_full_name;
            if (nargin == 2)
                file_name = file_name{id};
            end
        end

        function out = getNavEphFile(this)
            % Get the file name of the navigational files
            %
            % SYNTAX
            %   nav_path = this.getNavPath()
            out = this.eph_name;
        end

        function clk_file = getNavClkFile(this)
            % Get the file name of the clock files
            %
            % SYNTAX
            %   clk_path = this.getClkPath()
            clk_file = this.clk_name;
        end

        function erp_file = getErpFile(this)
            % Get the file name of the ERP files
            %
            % SYNTAX
            %    erp_path = this.getErpPath()
            erp_file = this.erp_name;
        end

        function dcb_file = getDcbFile(this)
            % Get the file name of the ERP files
            %
            % SYNTAX
            %   erp_path = this.getErpPath()
            dcb_file = this.dcb_name;
        end

        function out = getNavEphDir(this)
            % Get the path to the navigational files
            %
            % SYNTAX
            %   nav_path = this.getNavEphDir()
            out = this.eph_dir;
        end

        function out = getNavClkDir(this)
            % Get the path to the clock files
            %
            % SYNTAX
            %    nav_path = this.getClkPath()
            out = this.clk_dir;
        end

        function out = getErpDir(this)
            % Get the path to the ERP files
            %
            % SYNTAX
            %   erp_path = this.getErpPath()
            out = this.erp_dir;
        end

        function out = getDcbDir(this)
            % Get the path to the DCB files
            %
            % SYNTAX
            %   dcb_path = this.getDcbPath()
            out = this.dcb_dir;
        end

        function out = getIonoDir(this)
            % Get the path to the DCB files
            %
            % SYNTAX
            %   dcb_path = this.getDcbPath()
            out = this.iono_dir;
        end
        
        function out = getAtmLoadDir(this)
            % Get the path to the DCB files
            %
            % SYNTAX
            %   dcb_path = this.getDcbPath()
            out = this.atm_load_dir;
        end
        
        function out = getVMFDir(this)
            % Get the path to the DCB files
            %
            % SYNTAX
            %   dcb_path = this.getDcbPath()
            out = this.vmf_dir;
        end

        function out = getNavEphPath(this)
            % Get the path to the navigational files
            %
            % SYNTAX
            %   nav_path = this.getNavEphPath()
            out = File_Name_Processor.checkPath(strcat(this.eph_dir, filesep, this.eph_name));
        end

        function out = getNavClkPath(this)
            % Get the path to the clock files
            %
            % SYNTAX
            %   nav_path = this.getNavClkPath()
            out = File_Name_Processor.checkPath(strcat(this.clk_dir, filesep, this.clk_name));
        end

        function out = getErpPath(this)
            % Get the path to the ERP files
            %
            % SYNTAX
            %   erp_path = this.getErpPath()
            out = File_Name_Processor.checkPath(strcat(this.erp_dir, filesep, this.erp_name));
        end

        function out = getCrdFile(this)
            % Get the path of the stations coordinates file
            %
            % SYNTAX
            %   file_path = this.getCrdFile()
            if (isempty(this.crd_name))
                out = '';
            else
                crd_dir = File_Name_Processor.getFullDirPath(this.crd_dir, this.getHomeDir);
                out = File_Name_Processor.checkPath(strcat(crd_dir, filesep, this.crd_name));
            end
        end

        function out = getAtxFile(this)
            % Get the path of the antex file
            %
            % SYNTAX
            %   file_path = this.getAtxFile()
            if (isempty(this.atx_name))
                out = '';
            else
                out = this.checkAtxPath(strcat(this.atx_dir, filesep, this.atx_name));
            end
        end

        function out = getOceanFile(this)
            % Get the path of the ocean loading file
            %
            % SYNTAX
            %   file_path = this.getOceanFile()
            if (isempty(this.ocean_name))
                out = '';
            else
                out = this.checkOceanPath(strcat(this.ocean_dir, filesep, this.ocean_name));
            end
        end

        function out = getMetDir(this)
            % Get the path of the meteorological file dir
            %
            % SYNTAX
            %   file_path = this.getMetDir()
            if (isempty(this.met_dir))
                out = '';
            else
                out = this.checkMetPath(strcat(this.met_dir));
            end
        end

        function out = getMetFile(this, id)
            % Get the path of the meteorological file
            %
            % SYNTAX
            %   file_path = this.getMetFile()
            if (isempty(this.met_name))
                out = '';
            else
                if isempty(this.met_full_name)
                    this.updateMetFileName();
                end
                out = this.met_full_name;
                if iscell(out{1})
                    out = out{1};
                end
                if (nargin == 2)
                    if (id > length(out))
                        out = out{end};
                        this.log.addWarning(sprintf('The session "%d" is non-existent, using %s', id, out));
                    else
                        out = out{id};
                    end
                end
            end
        end

        function file_path = getIgrfFile(this)
            % Get the file name of the Mg
            %
            % SYNTAX
            %   file_path = this.getIgrfFile()
            file_path = File_Name_Processor.checkPath(strcat(this.igrf_dir, filesep, this.igrf_name));
        end

        function out = getGeoidDir(this)
            % Get the dir of the geoid file
            %
            % SYNTAX
            %   out = this.getGeoidDir()
            out = File_Name_Processor.checkPath(strcat(this.geoid_dir, filesep));
        end
        
        function file_path = getGeoidFile(this)
            % Get the path of the geoid file
            %
            % SYNTAX
            %   file_path = this.getGeoidFile()
            file_path = File_Name_Processor.checkPath(strcat(this.geoid_dir, filesep, this.geoid_name));
        end

        function out_dir = getOutDir(this)
            % Get the path of the out folder
            %
            % SYNTAX
            %   out_dir = this.getOutDir()
            out_dir = File_Name_Processor.checkPath(this.out_dir);
        end

        function out_prefix = getOutPrefix(this)
            % Get the path of the out_prefix
            %
            % SYNTAX
            %   out_prefix = this.getOutPrefix()
            fnp = File_Name_Processor;
            out_prefix = fnp.checkPath(this.out_prefix);
        end

        function updateOutPath(this, date, session)
            % Update the full prefix of the putput files (replacing special keywords)
            %
            % SYNTAX
            %   this.updateOutPath(date, session);
            %
            % NOTE: when no date is specified special keywords are substituted considering a date 0 (0000/00/00 00:00:00)
            %       when no session is specified special keywords are substituted considering a session "0" (char)

            fnp = File_Name_Processor;

            % get the output prefix
            narginchk(1,3);

            if (nargin < 2)
                date = GPS_Time(0);
            end
            if (nargin < 3)
                session = '0';
            end

            this.out_full_path = fnp.dateKeyRep(fnp.checkPath([this.out_dir filesep this.out_prefix]), date, session);

            if ~(this.run_counter_is_set)
                % make sure to have the name of the file and the name of the
                % output folder
                [~, out_prefix, ~] = fileparts(this.out_full_path); %#ok<PROPLC>
                % list all the files that match the prefix
                file_list = dir([this.out_full_path '*']);
                % if there are no files in the putput folder
                if isempty(file_list)
                    this.run_counter = 0; % set the counter of the output == 0
                else
                    % put the cell of the file in a single string
                    file_list = fnp.checkPath(strCell2Str({file_list(:).name},''));
                    % parse with regexp for output numbers -> get the maximum
                    this.run_counter = max(str2double(unique(regexp(file_list, [ '(?<=' out_prefix '_)[0-9]*(?=_)'], 'match')))) + 1; %#ok<PROPLC>
                    this.run_counter = iif(isempty(this.run_counter), this.RUN_COUNTER, this.run_counter);
                end
            end
            this.out_full_path = [ this.out_full_path '_' sprintf('%03d', this.run_counter)];
        end

        function out = getFullOutPath(this)
            % Get the path of the out folder composed with the prefix and the count number
            % update the run counter if necessary
            % SYNTAX: out_prefix = this.getOutPath()

            if isempty(this.out_full_path)
                this.log.addWarning('Output prefix has not yet been computed! It should have been done before.');
                this.updateOutPath();
            end
            out = this.out_full_path;
        end

        function counter = getRunCounter(this)
            % Get the currGPS_Time(0)ent run counter
            %
            % SYNTAX
            %   counter = getRunCounter(this)
            counter = this.run_counter;
        end
        
        function keep = isKeepRecList(this)
            % Get the flag to keep in memory all the receivers (CEBUG)
            %
            % SYNTAX
            %   keep = this.isKeepRecList();
            keep = this.flag_keep_rec_list;
        end
        
        function date = getSessionsStart(this)
            % Get the beginning of all the sessions
            % not considering buffer
            %
            % SYNTAX
            %   date = getSessionsStart(this)
            date = this.sss_date_start.getCopy;
        end
        
        function date = getSessionsStartExt(this)
            % Get the beginning of all the sessions
            % considering buffer
            %
            % SYNTAX
            %   date = getSessionsStartExt(this)
            date = this.sss_date_start.getCopy;
            buf_lft = this.getBuffer();
            date.addSeconds(-buf_lft); % left buffer
        end
        
        function response = isRinexSession(this)
            % is the program using session based on rinex files?
            %
            % SYNTAX
            %  response = this.isRinexSession()
            response = this.sss_file_based;
        end
       
        function date = getSessionsStop(this)
            % Get the end of all the sessions
            % not considering buffer
            %
            % SYNTAX
            %   date = getSessionsStop(this)
            date = this.sss_date_stop.getCopy;
            if this.sss_date_stop == this.sss_date_start
                date.addSeconds(86400 - 1e-4); % If session is empty set it to be one day long
            end            
        end
        
        function date = getSessionsStopExt(this)
            % Get the end of all the sessions
            % considering buffer
            %
            % SYNTAX
            %   date = getSessionsStopExt(this)
            date = this.sss_date_stop.getCopy;
            if this.sss_date_stop == this.sss_date_start
                date.addSeconds(86400 - 1e-4); % If session is empty set it to be one day long
            end
            [~,buf_rgt] = this.getBuffer();
            date.addSeconds(buf_rgt); % right buffer
        end

        function [sss_ext_lim, sss_lim] = getSessionLimits(this, n)
            % Get the time limits for a specific session 
            %
            % OUTPUT
            %   sss_ext_lim     limits of the data to be used for the computation
            %   sss_lim         limits of the session results (tipically smaller of "buffer" seconds)
            %
            % SYNTAX
            %   [sss_ext_lim, sss_lim] = getSessionLimits(this, <n>)
            if nargin < 2
                n = 1;
            end
            [buf_lft, but_rgt] = this.getBuffer();
            sss_lim = this.getSessionsStart;
            sss_lim.addSeconds((n-1) * this.sss_duration);
            sss_ext_lim = sss_lim.getCopy();
            sss_ext_lim.addSeconds(-buf_lft); % left buffer
            time_stop = sss_lim.getCopy();
            time_stop.addSeconds(this.sss_duration - 1e-4);
            if time_stop > this.getSessionsStop
                time_stop = this.getSessionsStop;
            end
            time_ext_stop = time_stop.getCopy();
            time_ext_stop.addSeconds(but_rgt); % right buffer
            
            sss_ext_lim.append(time_ext_stop);
            sss_lim.append(time_stop);
        end
        
        function [buf_lft, buf_rgt] = getBuffer(this)
            % get the session buffer
            %
            % SYNTAX
            %  [buf_lft, but_rgt] = this.getBuffer()
            buf_lft = this.sss_buffer(1)*(~this.isRinexSession());
            buf_rgt = this.sss_buffer(end)*(~this.isRinexSession());
        end

        function iono_model = getIonoModel(this)
            % SYNTAX
            %   iono_model = this.getIonoModel()
            iono_model = this.iono_model;
        end
        
        function iono_management = getIonoManagement(this)
            % SYNTAX
            %   date = getIonoManagement(this)
            iono_management = this.iono_management;
        end

        function eph_full_name = getEphFileName(this, date_start, date_stop)
            % Get the full name of the ephemerides files (replacing special keywords)
            %
            % SYNTAX
            %   eph_full_name = getEphFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.eph_dir, filesep, this.eph_name));
            step_sec = fnp.getStepSec(file_name);

            if (~isempty(strfind(file_name, fnp.GPS_WD)) || ~isempty(strfind(file_name, fnp.GPS_WEEK)))
                date_start = date_start.getCopy; date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
                date_stop = date_stop.getCopy; date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin
            end
            eph_full_name = fnp.dateKeyRepBatch(file_name, date_start,  date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end
        
         function met_full_name = getMetFileName(this, date_start, date_stop)
            % Get the full name of the ephemerides files (replacing special keywords)
            %
            % SYNTAX
            %   met_full_name = getEphFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            if ~iscell(this.met_name)
                met_name = {this.met_name};
            else
                met_name = this.met_name;
            end
            step_sec = fnp.getStepSec(met_name);

            if (~isempty(strfind(this.met_dir, fnp.GPS_WD)) || ~isempty(strfind(this.met_dir, fnp.GPS_WEEK)))
                date_start = date_start.getCopy; date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
                date_stop = date_stop.getCopy; date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin
            end
            met_full_name = {};
            for i = 1 : numel(met_name)
                met_full_name{i} = fnp.dateKeyRepBatch(fnp.checkPath(strcat(this.met_dir, filesep, met_name{i})), date_start,  date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop); %#ok<AGROW>
            end
        end

        function clk_full_name = getClkFileName(this, date_start, date_stop)
            % Get the full name of the clock offset files (replacing special keywords)
            %
            % SYNTAX
            %   clk_full_name = getClkFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.clk_dir, filesep, this.clk_name));
            step_sec = fnp.getStepSec(file_name);

            if (~isempty(strfind(file_name, fnp.GPS_WD)) || ~isempty(strfind(file_name, fnp.GPS_WEEK)))
                date_start = date_start.getCopy; date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
                date_stop = date_stop.getCopy; date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin
            end
            clk_full_name = fnp.dateKeyRepBatch(file_name, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end

        function erp_full_name = getErpFileName(this, date_start, date_stop)
            % Get the full name of the ERP files (replacing special keywords)
            %
            % SYNTAX
            %   erp_full_name = getErpFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.erp_dir, filesep, this.erp_name));

            if (~isempty(strfind(file_name, fnp.GPS_WD)) || ~isempty(strfind(file_name, fnp.GPS_WEEK)))
                date_start = date_start.getCopy;
                date_stop = date_stop.getCopy;
            end
            erp_full_name = fnp.dateKeyRepBatch(file_name, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end

        function iono_full_name = getIonoFileName(this, date_start, date_stop)
            % Get the full name of the ERP files (replacing special keywords)
            %
            % SYNTAX
            %   erp_full_name = getErpFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.iono_dir, filesep, this.iono_name));

            if (~isempty(strfind(file_name, fnp.GPS_WD)) || ~isempty(strfind(file_name, fnp.GPS_WEEK)))
                date_start = date_start.getCopy;
                date_stop = date_stop.getCopy;
            end
            iono_full_name = fnp.dateKeyRepBatch(file_name, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end
        
        function stm_load_full_name = getNTAtmLoadFileName(this, date_start, date_stop)
            % Get the full name of the ERP files (replacing special keywords)
            %
            % SYNTAX
            %   erp_full_name = getErpFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.atm_load_dir, filesep, this.atm_load_name_nt));

            if (~isempty(strfind(file_name, fnp.GPS_WD)) || ~isempty(strfind(file_name, fnp.GPS_WEEK)))
                date_start = date_start.getCopy;
                date_stop = date_stop.getCopy;
            end
            stm_load_full_name = fnp.dateKeyRepBatch(file_name, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end
        
        function stm_load_full_name = getTAtmLoadFileName(this)
            % Get the full name of the ERP files (replacing special keywords)
            %
            % SYNTAX
            %   erp_full_name = getErpFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            stm_load_full_name = fnp.checkPath(strcat(strrep(this.atm_load_dir, '${YYYY}',''), filesep, this.atm_load_name_t));
        end
        
        function vmf_full_name = getVMFFileName(this, date_start, date_stop)
            % Get the full name of the ERP files (replacing special keywords)
            %
            % SYNTAX
            %   erp_full_name = getErpFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.vmf_dir, filesep, this.vmf_name));

            date_start = date_start.getCopy;
            date_stop = date_stop.getCopy;
            fnp = File_Name_Processor();
            step_s = fnp.getStepSec(file_name);
            date_start.addSeconds(-step_s);
            date_stop.addSeconds(step_s);
            vmf_full_name = fnp.dateKeyRepBatch(file_name, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end
        
        function vmf_height_name = getVMFHeightFileName(this)
            % Get the full name of the ERP files (replacing special keywords)
            % SYNTAX: erp_full_name = getErpFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            vmf_height_name = fnp.checkPath(strcat(strrep(this.vmf_dir, '${YYYY}',''), filesep, 'orography_ell'));

        end

        function crx_full_name = getCrxFileName(this, date_start, date_stop)
            % Get the full name of the ERP files (replacing special keywords)
            %
            % SYNTAX
            %   erp_full_name = getErpFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.crx_dir, filesep, this.crx_name));

            if (~isempty(strfind(file_name, fnp.GPS_WD)) || ~isempty(strfind(file_name, fnp.GPS_WEEK)))
                date_start = date_start.getCopy;
                date_stop = date_stop.getCopy;
            end
            crx_full_name = fnp.dateKeyRepBatch(file_name, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end
    end
    
    % =========================================================================
    %%  SETTERS IO
    % =========================================================================
    methods        
        function setFile(this, filename)
            % set file
            %
            % SYNTAX
            %   setFile(this, filename)
            [~, fname, ext] = fileparts(filename);
            if strcmpi(ext,'.sp3') || strcmpi(ext,'.eph')  || strcmpi(ext,'.pre')
                this.setNavEphFile(filename);
            elseif strcmpi(ext,'.erp')
                this.setErpFile(filename);
            elseif instr(lower(ext),'.clk')
                this.setNavClkFile(filename);
            elseif strcmpi(ext,'.CRX')
            elseif strcmpi(ext,'.apl')
                this.setAtmLoadFile(filename);
            elseif ~isempty(regexp(ext,'\.\d\di', 'once')) || strcmpi(ext,'.${YY}i')
                this.setIonoFile(filename);
            elseif strcmpi(ext,'.DCB') || (strcmpi(ext,'.SNX') && strcmpi(name(1:3),'DCB'))
                this.setDcbFile(filename);
            elseif instr(fname,'VMFG_')
                this.setVMFFile(filename);
            end
        end

        function setNavPath(this, nav_dir)
            % Set the path to the navigational files
            %
            % SYNTAX
            %   this.getNavPath(nav_path)
            this.eph_dir = nav_dir;
        end

        function setMetDir(this, met_dir)
            % Set the path of the meteorological file dir
            %
            % SYNTAX
            %   this.setMetDir(met_dir);
            this.met_dir = met_dir;
        end

        function setNavEphFile(this, nav_name)
            % Set the file name of the navigational files
            %
            % SYNTAX
            %   this.setNavEphFile(nav_name)
            this.eph_name = nav_name;
        end

        function setNavClkPath(this, clk_dir)
            % Set the path to the clock files
            % 
            % SYNTAX
            %   this.getClkPath(nav_path)
            this.clk_dir = clk_dir;
        end

        function setNavClkFile(this, clk_name)
            % Set the file name of the clock files
            % 
            % SYNTAX
            %   this.getClkFile(nav_name)
            this.clk_name = clk_name;
        end

        function setErpPath(this, erp_dir)
            % Set the path to the clock files
            % 
            % SYNTAX
            %   this.getClkPath(nav_path)
            this.erp_dir = erp_dir;
        end

        function setErpFile(this, erp_name)
            % Set the file name of the clock files
            %
            % SYNTAX
            %   this.getClkFile(erp_name)
            this.erp_name = erp_name;
        end

        function setDcbFile(this, dcb_name)
            % Set the file name of the clock files
            % SYNTAX
            %   this.getClkFile(erp_name)
            this.dcb_name = dcb_name;
        end

        function setIonoFile(this, iono_file)
            % Set the file name of the clock files
            %
            % SYNTAX
            %   this.getClkFile(erp_name)
            this.iono_name = iono_file;
        end
        
        function setAtmLoadFile(this, atm_load_file)
            % Set the file name of the clock files
            %
            % SYNTAX
            %   this.getClkFile(erp_name)
            this.atm_load_name_nt = atm_load_file;
        end
        
        function setVMFFile(this, vmf_file)
            % Set the file name of the clock files
            %
            % SYNTAX
            %   this.getClkFile(erp_name)
            this.vmf_name = vmf_file;
        end

        function setIGRFFile(this, igrf_name)
            % Set the file name of the clock files
            %
            % SYNTAX
            %   this.getClkFile(erp_name)
            this.igrf_name = igrf_name;
        end

        function setOutPrefix(this, out_prefix)
            % Set the path of the out_prefix
            %
            % SYNTAX
            %   out_prefix = this.setOutPrefix(out_prefix)
            this.out_prefix = File_Name_Processor.checkPath(out_prefix);
        end

        function updateObsFileName(this)
            % Update the full name of the observations files (replacing special keywords)
            %
            % SYNTAX
            %   this.updateObsFileName();
            this.obs_full_name = {};
            fnp = File_Name_Processor();
            if ~iscell(this.obs_name)
                this.obs_name = {this.obs_name};
            end
            for i = 1 : numel(this.obs_name)
                this.obs_full_name{i} = fnp.dateKeyRepBatch(fnp.checkPath(strcat(this.obs_dir, filesep, this.obs_name{i})), this.getSessionsStartExt,  this.getSessionsStopExt, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
            end
        end

        function updateMetFileName(this)
            % Update the full name of the observations files (replacing special keywords)
            %
            % SYNTAX
            %   this.updateObsFileName();
            this.met_full_name = {};
            fnp = File_Name_Processor();
            if ~iscell(this.met_name)
                this.met_name = {this.met_name};
            end
            this.met_full_name = {};
            for i = 1 : numel(this.met_name)
                this.met_full_name = [this.met_full_name; fnp.dateKeyRepBatch(fnp.checkPath(strcat(this.met_dir, filesep, this.met_name{i})), this.getSessionsStartExt,  this.getSessionsStopExt, this.sss_id_list, this.sss_id_start, this.sss_id_stop)];
            end
            %this.met_full_name = this.met_full_name(~isempty(this.met_full_name));
        end

        function updateNavFileName(this)
            % Update the full name of the navigational files (replacing special keywords)
            %
            % SYNTAX
            %   this.updateNavFileName();
            this.updateEphFileName();
            this.updateClkFileName();
        end

        function updateEphFileName(this)
            % Update the full name of the ephemerides files (replacing special keywords)
            %
            % SYNTAX
            %   this.updateEphFileName();
            this.eph_full_name = this.getEphFileName(this.getSessionsStartExt, this.getSessionsStopExt);
        end

        function updateClkFileName(this)
            % Update the full name of the clock offset files (replacing special keywords)
            %  
            % SYNTAX
            %   this.updateClkFileName();
            this.clk_full_name = this.getClkFileName(this.getSessionsStartExt, this.getSessionsStopExt);
        end

        function updateErpFileName(this)
            % Update the full name of the ERP files (replacing special keywords)
            %
            % SYNTAX
            %   this.updateClkFileName();
            this.erp_full_name = this.getErpFileName(this.getSessionsStartExt, this.getSessionsStopExt);
        end

        function date = setSessionStart(this, date)
            % set session start
            %
            % SYNTAX
            %   setSessionStart(this, date)
            this.sss_date_start = date.getCopy;
        end

        function date = setSessionStop(this, date)
            % set session stop
            %
            % SYNTAX
            %   setSessionStop(this, date)
            this.sss_date_stop = date.getCopy();
        end
        
        function keep = setKeepRecList(this, keep)
            % Get the flag to keep in memory all the receivers (CEBUG)
            %
            % SYNTAX
            %   keep = this.isKeepRecList();
            this.flag_keep_rec_list = keep;
        end
        
        function setSmoothTropoOut(this, is_smt)
            % Should the troposphere paramteres be smoothed
            %
            % SYNTAX
            %   this.isSmoothTropoOut(is_smt)
            
            this.flag_smooth_tropo_out = is_smt;
        end
        
        function setRemoteSourceDir(this, dir_path)
            this.remote_res_conf_dir = fnp.getFullDirPath(dir_path, this.getHomeDir);
        end
        
        function setRemCheck(this, flag)
            % Set the Remote Check flag
            %
            % SYNTAX
            %   this.setRemoteCheck(flag)
            this.flag_check_remote = flag;
        end
 
        function setPrjHome(this, prj_home)
            % Set home folder of the project
            %
            % SYNTAX
            %   dir = this.setPrjHome(prj_home)
            this.prj_home = prj_home;
        end

        function setFilePath(this, file_path)
            % Set the file name of the current settings
            %   
            % SYNTAX
            %   this.setFilePath(file_path)
            [path_str, name, ~] = fileparts(file_path);
            this.cur_ini = [path_str filesep name '.ini'];
        end

        function updatePrj(this, file_path)
            % Set the file project name / home / file_path from file_path
            %
            % SYNTAX
            %   this.autoUpdatePrj(this, file_path)
            fnp = File_Name_Processor();

            [path_str, name, ~] = fileparts(fnp.checkPath(file_path));
            this.cur_ini = [path_str filesep name '.ini'];
            path_parts = strsplit(path_str,filesep);
            if numel(path_parts) > 3
                this.prj_home = fnp.checkPath(fullfile(path_parts{1:end-1}, filesep));
                this.prj_name = path_parts{end-1};
                this.log.addMessage('Trying to guess project name / home / ini');
                this.log.addMessage(sprintf(' name: %s', this.prj_name));
                this.log.addMessage(sprintf(' home: %s', this.prj_home));
                this.log.addMessage(sprintf(' ini:  %s', this.cur_ini));
            end
        end
    end
    
    
    % =========================================================================
    %%  GETTERS
    % =========================================================================
    methods (Access = 'public')
        function name = getPrjName(this)
            % get the project name
            %
            % SYNTAX
            %   name = this.getPrjName();

            name = this.prj_name;
        end

        function cc = getConstellationCollector(this)
            % Get the constellation collector object
            %
            % SYNTAX
            %   cc = this.getConstellationCollector()
            cc = handle(this.cc);
        end

        function cc = getCC(this)
            % Get the constellation collector object (short name)
            %
            % SYNTAX
            %   cc = this.getCC();
            cc = handle(this.cc);
        end

        function cut_off = getCutOff(this)
            % Get the cut off
            %
            % SYNTAX
            %   cut_off = this.getCutOff();
            cut_off = this.cut_off;
        end
        
        function snr_thr = getSnrThr(this)
            % Get the snr threshold
            %
            % SYNTAX
            %   snr_thr = this.getSnrThr();
            snr_thr = this.snr_thr;
        end

        function min_arc = getMinArc(this)
            % Get the minimum arc legnth to be kept
            %
            % SYNTAX
            %   min_arc = this.getMinArc()
            min_arc = this.min_arc;
        end

        function err_thr = getMaxCodeErrThr(this)
            % Get the maximum error acceptable on code observations
            %
            % SYNTAX
            %   err_thr = this.getMaxCodeErrThr()
            err_thr = this.pp_max_code_err_thr;
        end

        function err_thr = getMaxPhaseErrThr(this)
            % Get the maximum error acceptable on phase observations
            %
            % SYNTAX
            %   err_thr = this.getMaxPhaseErrThr()
            err_thr = this.pp_max_phase_err_thr;
        end
        
        function is_rem_check = isRemCheck(this)
            % Get the Remote Check flag
            %
            % SYNTAX
            %   is_rem_check = this.isRemCheck()
            is_rem_check = this.flag_check_remote;
        end
        
        function is_vmf = isVMF(this)
            % Get the VMF flag
            %
            % SYNTAX
            %   is_vmf = this.isVMF()
            is_vmf = this.mapping_function == 2 || this.zd_model == 2;
        end
        
        function is_seamless = isSeamlessKF(this)
            % Get the Seamless Rate flag
            %
            % SYNTAX
            %   is_seamless = this.isSeamlessKF();
            is_seamless = this.flag_seamless_proc;
        end

        function is_fb = isForwardBackwardKF(this)
            % Get the Forward Backward flag
            %
            % SYNTAX
            %   is_fb = this.isForwardBackwardKF();
            is_fb = this.flag_kf_fb ~= 0;
        end

        function kf_fb = getForwardBackwardKF(this)
            % Get the Forward Backward flag
            %
            % SYNTAX
            %   kf_fb = this.getForwardBackwardKF();
            kf_fb = this.flag_kf_fb;
        end

        function is_static = isStaticKF(this)
            % Check wether the current KF mode is static (PP)
            %
            % SYNTAX
            %   is_static = this.isStaticKF();
            is_static = (~this.isModeMonitor() && this.kf_mode == 0);
        end

        function is_variable = isVariableKF(this)
            % Check wether the current KF mode is variable
            %
            % SYNTAX
            %   is_variable = this.isVariableKF();
            is_variable = (this.isModeMonitor() && this.kf_mode == 1) || (~this.isModeMonitor() && this.kf_mode == 3);
        end

        function is_tropo = isTropoOn(this)
            % Check whether the tropospheric delay estimation is enabled
            %
            % SYNTAX
            %   is_tropo = this.isTropoOn();
            is_tropo = this.flag_tropo;
        end
        
        function is_iono_free = isIonoFree(this)
            % Check whether the iono free combination is enabled
            %
            % SYNTAX
            %   is_iono_free = isIonoFree(this)
            is_iono_free = this.iono_management == 1;
        end
        
        function is_iono_ext = isIonoExtModel(this)
            % Check whether the iono external model is enabled
            %
            % SYNTAX
            %   is_iono_ext = isIonoExtModel(this)
            is_iono_ext = this.iono_management == 3;
        end
        function is_solid_earth = isSolidEarth(this)
            % Check whether the iono free combination is enabled
            %
            % SYNTAX
            %   is_solid_earth = isSolidEarth(this)
            is_solid_earth = this.flag_solid_earth;
        end
        
        function is_pole_tide = isPoleTide(this)
            % Check whether the iono free combination is enabled
            %
            % SYNTAX
            %   is_pole_tide = isPoleTide(this)
            is_pole_tide = this.flag_pole_tide;
        end
        
        function is_phase_wind = isPhaseWind(this)
            % Check whether the iono free combination is enabled
            %
            % SYNTAX
            %   is_phase_wind = isPhaseWind(this)
            is_phase_wind = this.flag_phase_wind;
        end
        
        function is_shapiro = isShapiro(this)
            % Check whether the iono free combination is enabled
            %
            % SYNTAX
            %   is_shapiro = isShapiro(this)
            is_shapiro = this.flag_shapiro;
        end
        
        function is_met = isMet(this)
            % Check whether the meteorological usage is active
            %
            % SYNTAX
            %   is_met = isMet(this)
            is_met = this.meteo_data == 3;
        end             
        
        function is_ocean_load = isOceanLoading(this)
            % Check whether the iono free combination is enabled
            %
            % SYNTAX
            %   is_ocean_load = isOceanLoading(this)
            is_ocean_load = this.flag_ocean_load;
        end
        
        function is_atm_load = isAtmLoading(this)
            % Check whether the iono free combination is enabled
            %
            % SYNTAX
            %   is_atm_load = isAtmLoading(this)
            is_atm_load = this.flag_atm_load;
        end
        
        function is_hoi = isHOI(this)
            % Check whether high order ionspheric delays are enabled
            %
            % SYNTAX
            %   is_hoi = isHOI(this)
            is_hoi = this.flag_hoi;
        end
        
        function is_smt = isSmoothTropoOut(this)
            % Should the troposphere paramteres be smoothed
            %
            % SYNTAX
            %   is_smt = this.isSmoothTropoOut()
            
            is_smt = this.flag_smooth_tropo_out;            
        end

        function need_iono = needIonoMap(this)
            % Check if ionospheric map are needed
            % SYNTAX: need_iono = needIonoMap(this)
            need_iono = this.isHOI || (this.iono_model == 3 && this.iono_management == 3) || true;
        end
        
        function is_rec_pcv = isRecPCV(this)
            % Check whether the iono free combination is enabled
            %
            % SYNTAX
            %   is_rec_pcv= isRecPCV(this)
            is_rec_pcv = this.flag_rec_pcv;
        end       

        function is_tropo_gradient = isTropoGradientEnabled(this)
            % Check whether the tropospheric delay gradient estimation is enabled
            %
            % SYNTAX
            %   is_tropo_gradient = this.isTropoGradientEnabled();
            is_tropo_gradient = this.flag_tropo_gradient;
        end

        function s_rate = getSolutionRate(this)
            % Get the solution rate to be exported (high rate) in go Block
            %
            % SYNTAX
            %   s_rate = this.isTropoGradientEnabled();
            s_rate = this.s_rate;
        end

        function flag = isPreCleaningOn(this)
            % Try to correct cycle slips / discontinuities in the observations and increase spike variance
            %
            % SYNTAX
            %   flag = this.isPreCleaningOn()
            flag = this.block_pre_cleaning;
        end

        function n_loop = getBlockPostCleaningLoops(this)
            % Try to correct cycle slips / discontinuities in the observations and increase spike variance
            %
            % SYNTAX
            %   flag = this.useBlockPreCleaning()
            n_loop = this.block_post_cleaning_loops;
        end

        function flag = isSeamlessHR(this)
            % Compute ambiguities and the high rate solution as a unique system (true) / compute independent goBlock high rate solution (false)
            %
            % SYNTAX
            %   flag = this.useBlockSeamlessHR()
            flag = this.block_seamless_hr;
        end

        function flag = getFullSlipSplit(this)
            % Compute ambiguities and the high rate solution as a unique system (true) / compute independent goBlock high rate solution (false)
            %
            % SYNTAX
            %   flag = this.useBlockSeamlessHR()
            flag = this.block_full_slip_split;
        end

        function flag = isOutlierRejectionOn(this)
            % Get the status of outlier rejection
            %
            % SYNTAX
            %   flag = this.isOutlierRejectionOn()
            flag = this.flag_outlier;
        end

        function flag = isBlockForceStabilizationOn(this)
            % Get the status of outlier rejection OF UNSTABLE ARCS
            % SYNTAX
            %   flag = this.isBlockForceStabilizationOn()
            flag = this.block_force_stabilization;
        end

        function flag = isBlockOneArc(this)
            % Get the status of "one arc" approach
            %
            % SYNTAX
            %   flag = this.isBlockOneArc()
            flag = this.block_one_arc;
        end
    end
    
    % =========================================================================
    %%  TEST
    % =========================================================================
    methods (Static, Access = 'public')
        function test()
            % Test the class
            %
            % SYNTAX
            %   state.test()
            s = Main_Settings();
            s.testInterfaceRoutines();
        end
    end
end
