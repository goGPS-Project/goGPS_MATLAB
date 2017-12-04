%   CLASS IO_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the file_paths and folders
%
% EXAMPLE
%   settings = IO_Settings();
%
% FOR A LIST OF CONSTANTS and METHODS use doc IO_Settings

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef IO_Settings < Settings_Interface

    properties (Constant, Access = 'protected')
        % id to string of out modes
        DEFAULT_DIR_IN = ['..' filesep '..' filesep ];
        DEFAULT_DIR_OUT = ['..' filesep '..' filesep];
    end
    
    % Default values for each field - useful to restore corrupted field
    properties (Constant, Access = 'public')

        % PROJECT
        PRJ_NAME = 'Defauld DD';  % Name of the project
        PRJ_HOME = [fileparts(which('goGPS.m')) filesep '..' filesep 'data' filesep 'project' filesep 'default_DD' filesep]; % Location of the project <relative path from goGPS folder>
        CUR_INI = [IO_Settings.PRJ_HOME 'Config' filesep 'settings.ini']; % Location of the current ini file

        % DEPRECATE
        INPUT_FILE_INI_PATH = ''; % deprecate INI - it contains some additional setting (yet not imported in the new settings system)

        % RECEIVERS
        SSS_DATE_START = GPS_Time(); % Start of the processing session
        SSS_DATE_STOP = GPS_Time();  % End of the processing session
        SSS_ID_LIST = '0';   % id character sequence to be use for the session $(S) special keyword
        SSS_ID_START = '0';  % first session id (char of sss_id_list)
        SSS_ID_STOP = '0';   % last session id (char of sss_id_list)

        OBS_DIR = 'RINEX';
        OBS_NAME = {'yamatogawa_master.obs' 'yamatogawa_rover.obs'};
        REC_TARGET = 0;
        REC_MASTER = 1;
        REC_REFERENCE = 2;
        OBS_TYPE_LIST = {'Target', ...
                         'Master', ...
                         'SEID reference'};
        OBS_TYPE = [1 0];

        ATX_DIR = [IO_Settings.DEFAULT_DIR_IN 'antenna' filesep 'ATX' filesep]; % Location of the antex files
        ATX_NAME = 'igs14_1941.atx';    % Name antex file

        XYZ_ANT = zeros(3, 1);
        XYZ_EV_POINT = zeros(3, 1);

        % COMPUTATION CENTERS
        % With official products for orbits and clocks
        PREFERRED_ARCHIVE = {'cddis', 'custom'}
        PREFERRED_GPS = {'igs', 'cod', 'jpl', 'esa', 'emx', 'gbm'}
        PREFERRED_GLO = {'igs', 'cod', 'esa', 'emx', 'gbm'}
        PREFERRED_MXD = {'gbm'}
        PREFERRED_CLK = {'clk_05s', 'clk_30s', 'clk'}
        PREFERRED_ERP = {'final', 'rapid', 'ultra'}
        PREFERRED_EPH = {'final', 'rapid', 'ultra', 'broadcast'}

        CUSTOM_ADDR = 'cddis.gsfc.nasa.gov/'
        CUSTOM_PORT = '21'
        CUSTOM_PATH = 'pub/gps/products/'
        CUSTOM_NAME_EPH = '${WWWW}/igs${WWWWD}.sp3';
        CUSTOM_NAME_CLK = '${WWWW}/igs${WWWWD}.clk_30s';
        CUSTOM_NAME_ERP = '${WWWW}/igs${WWWW}7.erp';

        % SATELLITES
        EPH_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'EPH' filesep]; % Path to Ephemeris files folder
        EPH_NAME = ''; % Name for Ephemeris files
        CLK_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'CLK' filesep]; % Path to Clock Offset files folder
        CLK_NAME = ''; % Name of Clock Offset files
        CRX_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'CRX' filesep]; % Path to CRX folder containing files of Satellites problems
        DCB_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'DCB' filesep]; % Path to DCB folder containing files of Differential Code Biases
        EMS_DIR = [IO_Settings.DEFAULT_DIR_IN 'satellite' filesep 'SBAS' filesep 'EMS' filesep]; % Path to EMS folder containing files of EGNOS Message Server.
        % STATIONS
        CRD_DIR = [IO_Settings.DEFAULT_DIR_IN 'station' filesep 'CRD' filesep]; % Path to Ephemeris files folder
        CRD_NAME = '';    % Location of the stations coordinate file
        MET_DIR = [IO_Settings.DEFAULT_DIR_IN 'station' filesep 'MET' filesep]; % Path to Clock Offset files folder
        MET_NAME = '';    % Location of the meteorological file
        OCEAN_DIR = [IO_Settings.DEFAULT_DIR_IN 'station' filesep 'ocean' filesep]; % Path to CRX folder containing files of Satellites problems
        OCEAN_NAME = '';  % Location of the ocean loading file
        % REFERENCE
        REF_GRAPH_FILE = [IO_Settings.DEFAULT_DIR_IN 'reference' filesep 'ref_path' filesep 'ref_path.mat']; % Reference path constraints
        ERP_DIR = [IO_Settings.DEFAULT_DIR_IN 'reference' filesep 'ERP' filesep]; % Earth Rotation Parameters
        ERP_NAME = ''; % Name of ERP files
        GEOID_DIR = [IO_Settings.DEFAULT_DIR_IN 'reference' filesep 'geoid' filesep]; % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        GEOID_NAME = 'geoid_EGM2008_05.mat'; % File name of the Geoid containing the geoid to be used for the computation of hortometric heighs

        % DTM (SET PATH AND LOAD PARAMETER FILES)
        DTM_DIR = [IO_Settings.DEFAULT_DIR_IN 'reference' filesep 'DTM' filesep]; % Path to DTM folder containing DTM files
        % UI IMAGES
        IMG_DIR = [IO_Settings.DEFAULT_DIR_IN 'img' filesep];  % Path to images used by the interface
        
        % OUT PATH
        OUT_DIR = [IO_Settings.DEFAULT_DIR_OUT  'project' filesep 'default_DD' filesep 'out' filesep]; % Directory containing the output of the project
        OUT_PREFIX = 'out';  % Every time a solution is computed a folder with prefix followed by the run number is created
        RUN_COUNTER = [];     % This parameter store the current run number

        % OUT FLAGS
        FLAG_OUT_POSITION = true;           % Flag to export position files
        FLAG_OUT_SETTINGS = true;           % Flag to export the used setting file
        FLAG_OUT_PDF_REPORT = false;        % Flag to export the PDF report
        FLAG_OUT_PDF_TROPO = false;         % Flag to export PWV retrival
        FLAG_OUT_PDF_DT = false;            % Flag to export the PDF of the estimated dt
        FLAG_OUT_PDF_CODE_RES = false;      % Flag to export the PDF code residuals
        FLAG_OUT_PDF_PH_RES = false;        % Flag to export the PDF phase residuals
        FLAG_OUT_BLOCK_OBJ = false;         % Flag to export the Core_Block object used for the processing
        FLAG_OUT_KML = false;               % Flag to export the KML file
        FLAG_OUT_NMEA = false;              % Flag to export the NMEA file
        
        % EXTERNAL INFO as imported from the input ini file does not have default values
    end

    properties (Constant, Access = 'public')
        % Location of the latest project (the ini contains just a reference to the default project ini file - that is actually a settings file
        LAST_SETTINGS = 'last_settings.ini';
    end


    properties (SetAccess = protected, GetAccess = protected)
        % Location of the current ini file
        cur_ini = IO_Settings.CUR_INI;
    end

    properties (SetAccess = public, GetAccess = public)
        %------------------------------------------------------------------
        % PROJECT
        %------------------------------------------------------------------

        % Name of the project
        prj_name = IO_Settings.PRJ_NAME;

        % Location of the project <relative path from goGPS folder>
        prj_home = IO_Settings.PRJ_HOME;

        %------------------------------------------------------------------
        % COMPUTATION CENTERS
        %------------------------------------------------------------------
        % Centers for computation of orbits and other related parameters

        preferred_archive = IO_Settings.PREFERRED_ARCHIVE ; % order of the ftpo serv3r to use
        preferred_gps = IO_Settings.PREFERRED_GPS;          % order of products to use (the first found is used) for GPS only solutions
        preferred_glo = IO_Settings.PREFERRED_GLO;          % order of products to use (the first found is used) for GLONASS-GPS only solutions
        preferred_mxd = IO_Settings.PREFERRED_MXD;          % order of products to use (the first found is used) for MultiConstellation only solutions
        preferred_eph = IO_Settings.PREFERRED_EPH;          % kind of orbits to prefer
        preferred_erp = IO_Settings.PREFERRED_ERP;          % kind of erp to prefer
        preferred_clk = IO_Settings.PREFERRED_CLK;          % kind of clocks to prefer

        % Custom entry for a server

        custom_addr = IO_Settings.CUSTOM_ADDR; % ftp address for a custom server
        custom_port = IO_Settings.CUSTOM_PORT; % ftp port for a custom server
        custom_path = IO_Settings.CUSTOM_PATH; % remote path
        custom_name_eph = IO_Settings.CUSTOM_NAME_EPH; % name of the ephemeris file in the remote dir
        custom_name_clk = IO_Settings.CUSTOM_NAME_CLK; % name of the clock file in the remote dir
        custom_name_erp = IO_Settings.CUSTOM_NAME_ERP; % name of the Earth Rotation Parameters file in the remote dir

        %------------------------------------------------------------------
        % DEPRECATE
        %------------------------------------------------------------------
        % deprecate INI - it contains some additional setting (yet not imported in the new settings system)
        % the ini file is instantiated in the ext_ini property of this class
        input_file_ini_path = IO_Settings.INPUT_FILE_INI_PATH;

        %------------------------------------------------------------------
        % RECEIVERS
        %------------------------------------------------------------------

        % Observation files of the Receivers

        % session

        sss_date_start = IO_Settings.SSS_DATE_START;    % start of the processing session
        sss_date_stop =  IO_Settings.SSS_DATE_STOP;     % end of the processing session
        sss_id_list =    IO_Settings.SSS_ID_LIST;       % id character sequence to be use for the session $(S) special keyworc
        sss_id_start =   IO_Settings.SSS_ID_START;      % first session id (char of sss_id_list)
        sss_id_stop =    IO_Settings.SSS_ID_STOP;       % last session id (char of sss_id_list)

        % reference receivers (e.g. master, SEID reference)

        obs_dir = IO_Settings.OBS_DIR;    % Directory containing the data (static)
        obs_name = IO_Settings.OBS_NAME;  % File name of the receivers (can contain special keywords)
        obs_full_name;                    % Full name of the observations generated during runtime from the provided parameters
        obs_type = IO_Settings.OBS_TYPE;  % Array of observations type (target / master / reference)

        atx_dir = IO_Settings.ATX_DIR;  % Location of the antex file
        atx_name = IO_Settings.ATX_NAME;  % Location of the antex file

        %------------------------------------------------------------------
        % GEOMETRY
        %------------------------------------------------------------------

        xyz_ant = IO_Settings.XYZ_ANT;            % Relative position of each recever antenna
        xyz_ev_point = IO_Settings.XYZ_EV_POINT;  % Position of the evaluation point

        %------------------------------------------------------------------
        % DATA CENTER SERVERS
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % SATELLITES
        %------------------------------------------------------------------


        eph_dir = IO_Settings.EPH_DIR;    % Path to Ephemeris files folder
        eph_name = IO_Settings.EPH_NAME;  % File name of ephemeris
        eph_full_name;                    % Full name of the ephemeris generated during runtime from the provided parameters

        clk_dir = IO_Settings.CLK_DIR;    % Path to Clock Offset files folder
        clk_name = IO_Settings.CLK_NAME;  % File name of clock offsets
        clk_full_name;                    % Full name of the clock offsets generated during runtime from the provided parameters

        erp_dir = IO_Settings.ERP_DIR;    % Path to ERP files folder
        erp_name = IO_Settings.ERP_NAME;  % File name of ERP
        erp_full_name;                    % Full name of ERPs generated during runtime from the provided parameters

        % Path to CRX folder containing files of Satellites problems
        crx_dir = IO_Settings.CRX_DIR;
        % Path to DCB folder containing files of Differential Code Biases
        dcb_dir = IO_Settings.DCB_DIR;
        % Path to EMS folder containing files of EGNOS Message Server.
        ems_dir = IO_Settings.EMS_DIR;

        %------------------------------------------------------------------
        % STATIONS
        %------------------------------------------------------------------


        % Path to stations coordinates files
        crd_dir = IO_Settings.CRD_DIR;
        % Location of the stations coordinate file
        crd_name = IO_Settings.CRD_NAME;

        % Path to stations meteorological files
        met_dir = IO_Settings.MET_DIR;
        % Location of the meteorological file
        met_name =  IO_Settings.MET_NAME;
        met_full_name; % Full name of the met file generated during runtime from the provided parameters

        % Path to stations ocean loading files
        ocean_dir = IO_Settings.OCEAN_DIR;

        %------------------------------------------------------------------
        % REFERENCE
        %------------------------------------------------------------------

        % Path to file containing the reference path
        ref_graph_file = IO_Settings.REF_GRAPH_FILE;

        % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        geoid_dir = IO_Settings.GEOID_DIR;
        % Name of the Geoid file containing the geoid to be used for the computation of hortometric heighs
        geoid_name = IO_Settings.GEOID_NAME;

        % Location of the ocean loading file
        ocean_name =  IO_Settings.OCEAN_NAME;

        %------------------------------------------------------------------
        % DTM (SET PATH AND LOAD PARAMETER FILES)
        %------------------------------------------------------------------

        % Path to DTM folder containing DTM files
        dtm_dir = IO_Settings.DTM_DIR;

        %------------------------------------------------------------------
        % UI IMAGES
        %------------------------------------------------------------------

        % Path to images used by the interface
        img_dir = IO_Settings.IMG_DIR;

        %------------------------------------------------------------------
        % OUTPUT
        %------------------------------------------------------------------

        out_dir = IO_Settings.OUT_DIR;        % Directory containing the output of the project
        out_prefix = IO_Settings.OUT_PREFIX;  % Every time a solution is computed a folder with prefix followed by the run number is created
        out_full_path;                        % Full prefix of the putput files generated during runtime from the provided parameters

        % This parameter store the current run number
        run_counter = IO_Settings.RUN_COUNTER;
        run_counter_is_set = false; % When importing the run counter, check if is set -> when set overwrite output

        flag_out_position = IO_Settings.FLAG_OUT_POSITION;          % Flag to export position files
        flag_out_settings = IO_Settings.FLAG_OUT_SETTINGS;          % Flag to export the used    setting file
        flag_out_pdf_report = IO_Settings.FLAG_OUT_PDF_REPORT;      % Flag to export the PDF report
        flag_out_pdf_tropo = IO_Settings.FLAG_OUT_PDF_TROPO;        % Flag to export PWV retrival
        flag_out_pdf_code_res = IO_Settings.FLAG_OUT_PDF_CODE_RES;  % Flag to export the PDF code residuals
        flag_out_pdf_ph_res = IO_Settings.FLAG_OUT_PDF_PH_RES;      % Flag to export the PDF phase residuals
        flag_out_pdf_dt = IO_Settings.FLAG_OUT_PDF_DT;              % Flag to export the PDF of the estimated dt
        flag_out_block_obj = IO_Settings.FLAG_OUT_BLOCK_OBJ;        % Flag to export the Core_Block object used for the processing
        flag_out_kml = IO_Settings.FLAG_OUT_KML;                    % Flag to export the KML file
        flag_out_nmea = IO_Settings.FLAG_OUT_NMEA;                  % Flag to export the NMEA file

        %------------------------------------------------------------------
        % EXTERNAL INFO as imported from INPUT FILE INI
        %------------------------------------------------------------------

        % Ini file object containing info about external data
        ext_ini;
    end

    % =========================================================================
    %  INIT
    % =========================================================================
    methods
        function this = IO_Settings()
            % Creator of IO_settings - verbosity level (true/false) can be set or ini file
            this.initLogger();
        end
    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods
        function import(this, settings)
            % This function import IO (only) settings from another setting object
            this.importIO(settings);
        end

        function importIO(this, settings)
            % This function import IO (only) settings from another setting object
            fnp = File_Name_Processor();

            if isa(settings, 'Ini_Manager')
                % PROJECT
                this.prj_name   = fnp.checkPath(settings.getData('prj_name'));
                this.prj_home   = fnp.getFullDirPath(settings.getData('prj_home'), pwd);
                if ~exist(this.prj_home, 'dir')
                    this.log.addWarning(sprintf('Project home "%s" does not exist\nusing prj_home = "%s"', this.prj_home, pwd));
                    this.prj_home = pwd;
                end
                % COMPUTATION CENTERS
                this.preferred_archive = fnp.checkPath(settings.getData('preferred_archive'));
                this.preferred_gps = fnp.checkPath(settings.getData('preferred_gps'));
                this.preferred_glo = fnp.checkPath(settings.getData('preferred_glo'));
                this.preferred_mxd = fnp.checkPath(settings.getData('preferred_mxd'));
                this.preferred_eph = fnp.checkPath(settings.getData('preferred_eph'));
                this.preferred_erp = fnp.checkPath(settings.getData('preferred_erp'));
                this.preferred_clk = fnp.checkPath(settings.getData('preferred_clk'));
                % Custom entry for a server
                this.custom_addr = fnp.checkPath(settings.getData('custom_addr'));
                this.custom_port = fnp.checkPath(settings.getData('custom_port'));
                this.custom_path = fnp.checkPath(settings.getData('custom_path'));
                this.custom_name_eph = fnp.checkPath(settings.getData('custom_name_eph'));
                this.custom_name_clk = fnp.checkPath(settings.getData('custom_name_clk'));
                this.custom_name_erp = fnp.checkPath(settings.getData('custom_name_erp'));
                % DEPRECATE
                this.input_file_ini_path = fnp.checkPath(settings.getData('input_file_ini_path'));
                % RECEIVERS
                if ~(isempty(settings.getData('sss_date_start')))
                    this.sss_date_start = GPS_Time(datenum(settings.getData('sss_date_start')));
                end
                if ~(isempty(settings.getData('sss_date_stop')))
                    this.sss_date_stop = GPS_Time(datenum(settings.getData('sss_date_stop')));
                end
                this.sss_id_list = settings.getData('sss_id_list');
                this.sss_id_start = settings.getData('sss_id_start');
                this.sss_id_stop = settings.getData('sss_id_stop');
                this.obs_dir = fnp.getFullDirPath(settings.getData('obs_dir'), this.prj_home, pwd);
                this.obs_name = fnp.checkPath(settings.getData('obs_name'));
                this.obs_full_name = {};
                this.obs_type = settings.getData('obs_type');
                this.atx_dir    = fnp.getFullDirPath(settings.getData('atx_dir'), this.prj_home, pwd);
                this.atx_name   = fnp.checkPath(settings.getData('atx_name'));
                % GEOMETRY
                tmp_xyz_ant = zeros(3, this.getTargetCount());
                for r = 1 : this.getTargetCount()
                    tmp = settings.getData(sprintf('xyz_ant_%02d', r));
                    if (size(tmp,1) == 3) && (size(tmp,2) == 1)
                        tmp_xyz_ant(:,r) = tmp;
                    end
                end
                this.xyz_ant = tmp_xyz_ant;
                clear tmp tmp_xyz_ant;
                this.xyz_ev_point = settings.getData('xyz_ev_point');
                % SATELLITES
                this.eph_dir    = fnp.getFullDirPath(settings.getData('eph_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('eph_dir')), this.prj_home));
                this.eph_name   = fnp.checkPath(settings.getData('eph_name'));
                this.clk_dir    = fnp.getFullDirPath(settings.getData('clk_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('clk_dir')), this.prj_home));
                this.clk_name   = fnp.checkPath(settings.getData('clk_name'));
                this.crx_dir    = fnp.getFullDirPath(settings.getData('crx_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('crx_dir')), this.prj_home));
                this.dcb_dir    = fnp.getFullDirPath(settings.getData('dcb_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('dcb_dir')), this.prj_home));
                this.ems_dir    = fnp.getFullDirPath(settings.getData('ems_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('ems_dir')), this.prj_home));
                % STATIONS
                this.crd_dir    = fnp.getFullDirPath(settings.getData('crd_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('crd_dir')), this.prj_home));
                this.crd_name    = fnp.checkPath(settings.getData('crd_name'));
                this.met_dir    = fnp.getFullDirPath(settings.getData('met_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('met_dir')), this.prj_home));
                this.met_name    = fnp.checkPath(settings.getData('met_name'));
                this.met_full_name = {};
                this.ocean_dir  = fnp.getFullDirPath(settings.getData('ocean_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('ocean_dir')), this.prj_home));
                this.ocean_name  = fnp.checkPath(settings.getData('ocean_name'));
                % REFERENCE
                this.ref_graph_file  = fnp.checkPath(settings.getData('ref_graph_file'));
                this.erp_dir    = fnp.getFullDirPath(settings.getData('erp_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('erp_dir')), this.prj_home));
                this.erp_name   = fnp.checkPath(settings.getData('erp_name'));
                this.geoid_dir  = fnp.getFullDirPath(settings.getData('geoid_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('geoid_dir')), this.prj_home));
                this.geoid_name = fnp.checkPath(settings.getData('geoid_name'));
                this.dtm_dir    = fnp.getFullDirPath(settings.getData('dtm_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('dtm_dir')), this.prj_home));
                % UI
                this.img_dir    = fnp.getFullDirPath(settings.getData('img_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('img_dir')), this.prj_home));
                % OUTPUT
                this.out_dir = fnp.getFullDirPath(settings.getData('out_dir'), this.prj_home, [], fnp.getFullDirPath(this.(upper('out_dir')), this.prj_home));
                if ~exist(this.out_dir, 'dir')
                    % fallback of fallback
                    this.out_dir = fnp.getFullDirPath(settings.getData('out_dir'), this.prj_home);
                end
                this.out_prefix = fnp.checkPath(settings.getData('out_prefix'));
                this.run_counter = settings.getData('run_counter');
                this.run_counter_is_set = ~isempty(this.run_counter);
                this.flag_out_position = settings.getData('flag_out_position');
                this.flag_out_settings = settings.getData('flag_out_settings');
                this.flag_out_pdf_report = settings.getData('flag_out_pdf_report');
                this.flag_out_pdf_tropo = settings.getData('flag_out_pdf_tropo');
                this.flag_out_pdf_dt = settings.getData('flag_out_pdf_dt');
                this.flag_out_pdf_code_res = settings.getData('flag_out_pdf_code_res');
                this.flag_out_pdf_ph_res = settings.getData('flag_out_pdf_ph_res');
                this.flag_out_block_obj = settings.getData('flag_out_block_obj');
                this.flag_out_kml = settings.getData('flag_out_kml');
                this.flag_out_nmea = settings.getData('flag_out_nmea');
            else
                % PROJECT
                this.prj_name   = settings.prj_name;
                this.prj_home   = settings.prj_home;
                %this.cur_ini   = settings.cur_ini;

                % COMPUTATION CENTERS
                this.preferred_archive = settings.preferred_archive;
                this.preferred_gps = settings.preferred_gps;
                this.preferred_glo = settings.preferred_glo;
                this.preferred_mxd = settings.preferred_mxd;
                this.preferred_eph = settings.preferred_eph;
                this.preferred_erp = settings.preferred_erp;
                this.preferred_clk = settings.preferred_clk;
                % Custom entry for a server
                this.custom_addr = settings.custom_addr;
                this.custom_port = settings.custom_port;
                this.custom_path = settings.custom_path;
                this.custom_name_eph = settings.custom_name_eph;
                this.custom_name_clk = settings.custom_name_clk;
                this.custom_name_erp = settings.custom_name_erp;
                % DEPRECATE
                this.input_file_ini_path = settings.input_file_ini_path;
                % RECEIVERS
                this.sss_date_start = settings.sss_date_start;
                this.sss_date_stop = settings.sss_date_stop;
                this.sss_id_list = settings.sss_id_list;
                this.sss_id_start = settings.sss_id_start;
                this.sss_id_stop = settings.sss_id_stop;
                this.obs_dir = settings.obs_dir;
                this.obs_name = settings.obs_name;
                this.obs_full_name = {};
                this.obs_type = settings.obs_type;
                this.atx_dir     = settings.atx_dir;
                this.atx_name    = settings.atx_name;
                this.xyz_ant = settings.xyz_ant;
                this.xyz_ev_point = settings.xyz_ev_point;
                % SATELLITES
                this.eph_dir     = settings.eph_dir;
                this.eph_name    = settings.eph_name;
                this.clk_dir     = settings.clk_dir;
                this.clk_name    = settings.clk_name;
                this.crx_dir     = settings.crx_dir;
                this.dcb_dir     = settings.dcb_dir;
                this.ems_dir     = settings.ems_dir;
                % STATIONS
                this.crd_dir     = settings.crd_dir;
                this.crd_name    = settings.crd_name;
                this.met_dir     = settings.met_dir;
                this.met_name    = settings.met_name;
                this.met_full_name = {};
                this.ocean_dir   = settings.ocean_dir;
                this.ocean_name  = settings.ocean_name;
                % REFERENCE
                this.ref_graph_file = settings.ref_graph_file;
                this.geoid_dir  = settings.geoid_dir;
                this.geoid_name = settings.geoid_name;
                this.dtm_dir    = settings.dtm_dir;
                % UI
                this.img_dir    = settings.img_dir;
                % OUTPUT
                this.out_dir = settings.out_dir;
                this.out_prefix = settings.out_prefix;
                this.run_counter = settings.run_counter;
                this.run_counter_is_set = ~isempty(this.run_counter);

                this.flag_out_position = settings.flag_out_position;
                this.flag_out_settings = settings.flag_out_settings;
                this.flag_out_pdf_report = settings.flag_out_pdf_report;
                this.flag_out_pdf_tropo = settings.flag_out_pdf_tropo;
                this.flag_out_pdf_dt = settings.flag_out_pdf_dt;
                this.flag_out_pdf_code_res = settings.flag_out_pdf_code_res;
                this.flag_out_pdf_ph_res = settings.flag_out_pdf_ph_res;
                this.flag_out_block_obj = settings.flag_out_block_obj;
                this.flag_out_kml = settings.flag_out_kml;
                this.flag_out_nmea = settings.flag_out_nmea;
            end
            this.check();
            this.updateExternals();
            this.eph_full_name = '';
            this.clk_full_name = '';
            this.erp_full_name = '';
            this.updateObsFileName();
        end

        function str = toString(this, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end
            fnp = File_Name_Processor;

            str = [str '---- PROJECT --------------------------------------------------------------' 10 10];
            str = [str sprintf(' Project name:                                     %s\n', this.prj_name)];
            str = [str sprintf(' Project home:                                     %s\n', fnp.getRelDirPath(this.prj_home, pwd))];
            str = [str sprintf(' Path to the current project ini file:             %s\n\n', this.cur_ini)];
            str = [str '---- DEPRECATE ------------------------------------------------------------' 10 10];
            str = [str sprintf(' Deprecate ini (of additional parameters):         %s\n\n', this.input_file_ini_path)];
            str = [str '---- GNSS RECEIVER FILES  -------------------------------------------------' 10 10];
            str = [str sprintf(' Definition of the file names to be parsed\n')];
            if ~(this.sss_date_start.isempty)
                str = [str sprintf(' Session start at                              %s \n', this.sss_date_start.toString())];
                str = [str sprintf(' Session end at                                %s \n', this.sss_date_start.toString())];
            end
            str = [str sprintf(' Character sequence to be used for the sessions    %s \n', this.sss_id_list)];
            str = [str sprintf(' First session char                                %c \n', this.sss_id_start)];
            str = [str sprintf(' Last session char                                 %c \n\n', this.sss_id_start)];
            str = [str sprintf(' Directory of the observation files                %s \n', fnp.getRelDirPath(this.obs_dir, this.prj_home))];
            str = [str sprintf(' Name of the observation files                     %s \n', strCell2Str(this.obs_name, ', '))];
            str = [str sprintf(' Array of observation types                        %s \n\n', strCell2Str(this.OBS_TYPE_LIST(this.obs_type + 1), ', '))];
            str = [str sprintf(' Directory of antennas (atx) files                 %s \n', this.atx_dir)];
            str = [str sprintf(' Antenna antex (ATX) file                          %s \n\n', this.atx_name)];
            if this.getTargetCount() > 1
                str = [str '---- GEOMETRY  ------------------------------------------------------------' 10 10];
                str = [str sprintf(' When multiple antennas are in use the positions on the structure are here defined:\n')];
                for r = 1 : size(this.xyz_ant, 2)
                    str = [str sprintf(' Position of the %02d receiver                       %d, %d, %d\n', r, this.xyz_ant(1, r), this.xyz_ant(2, r), this.xyz_ant(3, r))]; %#ok<AGROW>
                end
                str = [str sprintf(' Evaluation Point                                  %s \n\n', this.xyz_ev_point)];
            end
            str = [str '---- COMPUTATION CENTER  --------------------------------------------------' 10 10];
            str = [str sprintf(' List of server to be used for downloading ephemeris\n')];
            str = [str sprintf(' Preferred order of servers:                       %s\n', strCell2Str(this.preferred_archive))];
            str = [str sprintf(' Preferred order of GPS products:                  %s\n', strCell2Str(this.preferred_gps))];
            str = [str sprintf(' Preferred order of GPS/GLONASS products:          %s\n', strCell2Str(this.preferred_glo))];
            str = [str sprintf(' Preferred order of Multiconstellation products:   %s\n', strCell2Str(this.preferred_mxd))];
            str = [str sprintf(' Preferred order for orbits products:              %s\n', strCell2Str(this.preferred_eph))];
            str = [str sprintf(' Preferred order for erp products:                 %s\n', strCell2Str(this.preferred_erp))];
            str = [str sprintf(' Preferred order for clock products:               %s\n\n', strCell2Str(this.preferred_clk))];
            str = [str sprintf(' Custom server parameters:\n')];
            str = [str sprintf('  address:                                         %s\n', this.custom_addr)];
            str = [str sprintf('  port:                                            %s\n', this.custom_port)];
            str = [str sprintf('  path:                                            %s\n', this.custom_path)];
            str = [str sprintf('  eph name:                                        %s\n', this.custom_name_eph)];
            str = [str sprintf('  clk name:                                        %s\n', this.custom_name_clk)];
            str = [str sprintf('  erp name:                                        %s\n\n', this.custom_name_erp)];
            str = [str '---- INPUT FOLDERS: SATELLITE ---------------------------------------------' 10 10];
            str = [str sprintf(' Directory of Ephemeris files:                     %s\n', fnp.getRelDirPath(this.eph_dir, this.prj_home))];
            str = [str sprintf(' Name of Ephemeris files:                          %s\n', this.eph_name)];
            str = [str sprintf(' Directory of Satellite clock offsets:             %s\n', fnp.getRelDirPath(this.clk_dir, this.prj_home))];
            str = [str sprintf(' Name of Satellite clock offsets:                  %s\n', this.clk_name)];
            str = [str sprintf(' Directory of CRX (satellite problems):            %s\n', fnp.getRelDirPath(this.crx_dir, this.prj_home))];
            str = [str sprintf(' Directory of DCB (Differential Code Biases):      %s\n', fnp.getRelDirPath(this.dcb_dir, this.prj_home))];
            str = [str sprintf(' Directory of EMS (EGNOS Message Server):          %s\n\n', fnp.getRelDirPath(this.ems_dir, this.prj_home))];
            str = [str '---- INPUT FOLDERS: STATIONS ----------------------------------------------' 10 10];
            str = [str sprintf(' Directory of coordinates file:                    %s\n', fnp.getRelDirPath(this.crd_dir, this.prj_home))];
            str = [str sprintf(' Name of coordinate (CRD) file:                    %s\n', this.crd_name)];
            str = [str sprintf(' Directory of meteorological data:                 %s\n', fnp.getRelDirPath(this.met_dir, this.prj_home))];
            str = [str sprintf(' Name of meteorological (met) files:               %s\n', strCell2Str(this.met_name))];
            str = [str sprintf(' Directory of ocean loading files:                 %s\n', fnp.getRelDirPath(this.ocean_dir, this.prj_home))];
            str = [str sprintf(' Name of ocean loading file:                       %s\n\n', this.ocean_name)];
            str = [str '---- INPUT FOLDERS: REFERENCE ---------------------------------------------' 10 10];
            str = [str sprintf(' File contraining the reference graph:             %s\n', this.ref_graph_file)];
            str = [str sprintf(' Directory of ERP files:                           %s\n', fnp.getRelDirPath(this.erp_dir, this.prj_home))];
            str = [str sprintf(' Name of ERP files:                                %s\n', this.erp_name)];
            str = [str sprintf(' Directory of Geoid models:                        %s\n', fnp.getRelDirPath(this.geoid_dir, this.prj_home))];
            str = [str sprintf(' Name of the Geoid map file:                       %s\n', this.geoid_name)];
            str = [str sprintf(' Directory of DTM data:                            %s\n\n', fnp.getRelDirPath(this.dtm_dir, this.prj_home))];
            str = [str '---- INPUT FOLDERS: UI ----------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of images for UI:                       %s\n\n', fnp.getRelDirPath(this.img_dir, this.prj_home))];
            str = [str '---- OUTPUT SETTINGS ------------------------------------------------------' 10 10];
            str = [str sprintf(' Directory containing the output of the project:   %s\n', fnp.getRelDirPath(this.out_dir, this.prj_home))];
            str = [str sprintf(' Prefix of each run:                               %s\n', this.out_prefix)];
            str = [str sprintf(' Run counter:                                      %d\n', this.run_counter)];
            if (this.run_counter_is_set)
                str = [str sprintf(' Run counter has been set manually => overwriting output\n\n')];
            else
                str = [str sprintf(' Run counter has not been previously set \n => it will be set automatically to avoid overwriting of the oputputs\n\n')];
            end
            str = [str sprintf(' Export positions:                                 %d\n', this.flag_out_position)];
            str = [str sprintf(' Export settings:                                  %d\n', this.flag_out_settings)];
            str = [str sprintf(' Export PDF report:                                %d\n', this.flag_out_pdf_report)];
            str = [str sprintf(' Export PDF troposphere:                           %d\n', this.flag_out_pdf_tropo)];
            str = [str sprintf(' Export PDF clock dt:                              %d\n', this.flag_out_pdf_dt)];
            str = [str sprintf(' Export PDF code residuals:                        %d\n', this.flag_out_pdf_code_res)];
            str = [str sprintf(' Export PDF phase residuals:                       %d\n', this.flag_out_pdf_ph_res)];
            str = [str sprintf(' Export goBlock OBJECT:                            %d\n', this.flag_out_block_obj)];
            str = [str sprintf(' Export KML:                                       %d\n', this.flag_out_kml)];
            str = [str sprintf(' Export NMEA:                                      %d\n\n', this.flag_out_nmea)];
        end

        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the obj
            str_cell = this.exportIO(str_cell);
        end

        function str_cell = exportIO(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the obj
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
            % DEPRECATE
            str_cell = Ini_Manager.toIniStringSection('DEPRECATE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Deprecate ini - path - it contains some legacy setting (now imported in the new settings system)', str_cell);
            str_cell = Ini_Manager.toIniStringComment('[ WW ] if not empty the system will import its parameter overriding the one written in here', str_cell);
            str_cell = Ini_Manager.toIniString('input_file_ini_path', this.input_file_ini_path, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % RECEIVERxS
            str_cell = Ini_Manager.toIniStringSection('RECEIVER_FILES', str_cell);
            str_cell = Ini_Manager.toIniStringComment('"sss_" parameters define the session of observation, they are used to substitute special keywords in file names', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Working session - first data of observation to consider (yyyy-mm-dd <HH:MM:SS>)', str_cell);
            str_cell = Ini_Manager.toIniStringComment('mainly used to detect the name of the file to process', str_cell);
            str_cell = Ini_Manager.toIniString('sss_date_start', this.sss_date_start.toString('yyyy-mm-dd HH:MM:SS'), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Working session - last data of observation to consider (yyyy-mm-dd <HH:MM:SS>)', str_cell);
            str_cell = Ini_Manager.toIniString('sss_date_stop', this.sss_date_stop.toString('yyyy-mm-dd HH:MM:SS'), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Id character sequence to be use for the session $(S) special keyword (e.g. "01233456789ABCabc")', str_cell);
            str_cell = Ini_Manager.toIniString('sss_id_list', this.sss_id_list, str_cell);
            str_cell = Ini_Manager.toIniStringComment('First session id (char of sss_id_list)', str_cell);
            str_cell = Ini_Manager.toIniString('sss_id_start', this.sss_id_start, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Last session id (char of sss_id_list)', str_cell);
            str_cell = Ini_Manager.toIniString('sss_id_stop', this.sss_id_stop, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = File_Name_Processor.toIniStringComment(str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory containing the data (static)', str_cell);
            str_cell = Ini_Manager.toIniString('obs_dir', fnp.getRelDirPath(this.obs_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('File name of the receivers (can contain special keywords)', str_cell);
            str_cell = Ini_Manager.toIniString('obs_name', this.obs_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('Array of observations type (%s)', strCell2EnumStr(this.OBS_TYPE_LIST,', ')), str_cell);
            str_cell = Ini_Manager.toIniString('obs_type', this.obs_type, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of PCO - PCV antex (ATX) files', str_cell);
            str_cell = Ini_Manager.toIniString('atx_dir', fnp.getRelDirPath(this.atx_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('PCO - PCV antex (ATX) file', str_cell);
            str_cell = Ini_Manager.toIniString('atx_name', this.atx_name, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            % GEOMETRY
            str_cell = Ini_Manager.toIniStringSection('GEOMETRY', str_cell);
            str_cell = Ini_Manager.toIniStringComment('When using multiple receivers it is possible to set the relative position of each', str_cell);
            str_cell = Ini_Manager.toIniStringComment('receiver w.r.t. a local reference frame', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Positions are expressed in a XYZ RF as an array "xyz_ant_${NN}" where $NN is', str_cell);
            str_cell = Ini_Manager.toIniStringComment('the (target) receiver number (as ordered in obs_name) starting from 01', str_cell);
            for r = 1 : this.getTargetCount()
                str_cell = Ini_Manager.toIniString(sprintf('xyz_ant_%02d', r), this.xyz_ant(:, r), str_cell);
            end
            str_cell = Ini_Manager.toIniStringComment('Indicate the "baricenter" of the structure -> evaluation point', str_cell);
            str_cell = Ini_Manager.toIniString('xyz_ev_point', this.xyz_ev_point, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);

            % COMPUTATION CENTERS
            str_cell = Ini_Manager.toIniStringSection('COMPUTATION_CENTER', str_cell);
            str_cell = Ini_Manager.toIniStringComment('List of the computeation center to be used for ephemeris retrival', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Every product is searched locally, when not found is downloaded', str_cell);
            str_cell = Ini_Manager.toIniStringComment('When the file is not found, the system fall back on the next available', str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred list of web archives,', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_ARCHIVE)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_archive', this.preferred_archive, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred source for GPS only solution, ', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_GPS)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_gps', this.preferred_gps, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred source for GPS/GLONASS only solution,', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_GLO)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_glo', this.preferred_glo, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred source for Multiconstellation solution,', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_MXD)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_mxd', this.preferred_mxd, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred ephemeris type, valid only for source "igs",', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_EPH)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_eph', this.preferred_eph, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred Earth rotation parameters type, rapid are not always available,', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_ERP)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_erp', this.preferred_erp, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Preferred clock types, valid but for "igs" glonass,', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('accepted values: %s', Ini_Manager.strCell2Str(this.PREFERRED_CLK)), str_cell);
            str_cell = Ini_Manager.toIniString('preferred_clk', this.preferred_clk, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('A custom center/product can also be used', str_cell);
            str_cell = Ini_Manager.toIniStringComment('(all the fields are strings)', str_cell);
            str_cell = Ini_Manager.toIniString('custom_addr', this.custom_addr, str_cell);
            str_cell = Ini_Manager.toIniString('custom_port', this.custom_port, str_cell);
            str_cell = Ini_Manager.toIniString('custom_path', this.custom_path, str_cell);
            str_cell = Ini_Manager.toIniStringComment('The "variable" part of the path should be in the name, e.g. ${WWWW}/igs${WWWWD}.sp3', str_cell);
            str_cell = Ini_Manager.toIniString('custom_name_eph', this.custom_name_eph, str_cell);
            str_cell = Ini_Manager.toIniString('custom_name_clk', this.custom_name_clk, str_cell);
            str_cell = Ini_Manager.toIniString('custom_name_erp', this.custom_name_erp, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % SATELLITES
            str_cell = Ini_Manager.toIniStringSection('INPUT_SATELLITE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Ephemeris files', str_cell);
            str_cell = Ini_Manager.toIniString('eph_dir', fnp.getRelDirPath(this.eph_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of Ephemeris files - special keywords can be used', str_cell);
            str_cell = Ini_Manager.toIniStringComment('If not found, goGPS will try to download them following COMPUTATION_CENTER section', str_cell);
            str_cell = Ini_Manager.toIniString('eph_name', this.eph_name, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of clock offset files', str_cell);
            str_cell = Ini_Manager.toIniString('clk_dir', fnp.getRelDirPath(this.clk_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('If not found, goGPS will try to download them following COMPUTATION_CENTER section', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of clock offset files - special keywords can be used', str_cell);
            str_cell = Ini_Manager.toIniString('clk_name', this.clk_name, str_cell);
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
            % STATIONS
            str_cell = Ini_Manager.toIniStringSection('INPUT_STATIONS', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of coordinates files', str_cell);
            str_cell = Ini_Manager.toIniString('crd_dir', fnp.getRelDirPath(this.crd_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of coordinates (CRD) file', str_cell);
            str_cell = Ini_Manager.toIniString('crd_name', this.crd_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of meteorological data', str_cell);
            str_cell = Ini_Manager.toIniString('met_dir', fnp.getRelDirPath(this.met_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Meteorological file', str_cell);
            str_cell = Ini_Manager.toIniString('met_name', this.met_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of ocean loading files', str_cell);
            str_cell = Ini_Manager.toIniString('ocean_dir', fnp.getRelDirPath(this.ocean_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of ocean loading file', str_cell);
            str_cell = Ini_Manager.toIniString('ocean_name', this.ocean_name, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % REFERENCE
            str_cell = Ini_Manager.toIniStringSection('INPUT_REFERENCE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Earth rotation/orientation parameters (ERP) files', str_cell);
            str_cell = Ini_Manager.toIniString('erp_dir', fnp.getRelDirPath(this.erp_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('If not found, goGPS will try to download them following COMPUTATION_CENTER section', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of ERP files - special keywords can be used', str_cell);
            str_cell = Ini_Manager.toIniString('erp_name', this.erp_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Geoid files', str_cell);
            str_cell = Ini_Manager.toIniString('geoid_dir', fnp.getRelDirPath(this.geoid_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Filename in Geoid dir containing the map of ondulation of the geoid', str_cell);
            str_cell = Ini_Manager.toIniString('geoid_name', this.geoid_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DTM data', str_cell);
            str_cell = Ini_Manager.toIniString('dtm_dir', fnp.getRelDirPath(this.dtm_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('File containing a graph of path constraints', str_cell);
            str_cell = Ini_Manager.toIniString('ref_graph_file', this.ref_graph_file, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % UI
            str_cell = Ini_Manager.toIniStringSection('INPUT_UI', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of images for UI', str_cell);
            str_cell = Ini_Manager.toIniString('img_dir', fnp.getRelDirPath(this.img_dir, this.prj_home), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Path to the image of the logo 64x64 px', str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
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
            str_cell = Ini_Manager.toIniStringComment('Flag to export positions results', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_position', this.flag_out_position, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to export the used settings as ini file', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_settings', this.flag_out_settings, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to write the PDF report', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_pdf_report', this.flag_out_pdf_report, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to export PDF troposphere results (when available)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_pdf_tropo', this.flag_out_pdf_tropo, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to write the PDF of the code residuals', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_pdf_code_res', this.flag_out_pdf_code_res, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to write the PDF of the phase residuals', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_pdf_ph_res', this.flag_out_pdf_ph_res, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to write the PDF of the estimated dt', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_pdf_dt', this.flag_out_pdf_dt, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to export the object that coputed the block solution (if computed)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_block_obj', this.flag_out_block_obj, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to export the KML file (when expected)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_kml', this.flag_out_kml, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Flag to export the NMEA file (when expected)', str_cell);
            str_cell = Ini_Manager.toIniString('flag_out_nmea', this.flag_out_nmea, str_cell);
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

            % RECEIVER DEFAULT PARAMETERS ---------------------------------
            try
                if isfield(state, 'gogps_data_output')
                    this.out_dir = state.gogps_data_output;
                end
                if isfield(state, 'gogps_data_output_prefix')
                    this.out_prefix = state.gogps_data_output_prefix;
                end
                if isfield(state, 'INIsettings')
                    this.input_file_ini_path = state.INIsettings;
                end
            catch ex
                this.log.addWarning(['Legacy import "IO file / folders" failed - ', ex.message])
            end
            this.check();
            this.updateExternals();
        end
    end

    % =========================================================================
    %  GETTERS
    % =========================================================================
    methods
        function file_dir = getHomeDir(this)
            % Get the base directory containing the project
            file_dir = this.prj_home;
        end

        function base_rinex_dir = getRinexBaseDir(this)
            % Get the base directory containing RINEX files
            base_rinex_dir = this.obs_dir();
        end

        function num_receiver = getTargetCount(this)
            % Get the number of Target receivers
            num_receiver = sum(this.obs_type == this.REC_TARGET);
        end

        function num_receiver = getMasterCount(this)
            % Get the number of Master receivers
            num_receiver = sum(this.obs_type == this.REC_MASTER);
        end

        function num_receiver = getReferenceCount(this)
            % Get the number of Reference receivers
            num_receiver = sum(this.obs_type == this.REC_REFERENCE);
        end

        function file_name = getTargetPath(this, id)
            % Get the file list of target receivers files
            % SYNTAX: file_name = this.getTargetPath()
            % A cell for each receiver containing the list of names as cell
            if isempty(this.obs_full_name)
                this.updateObsFileName();
            end
            file_name = this.obs_full_name(this.obs_type == this.REC_TARGET);
            if (nargin == 2)
                file_name = file_name{id};
            end
        end

        function file_name = getMasterPath(this)
            % Get the file list of master receivers files
            % SYNTAX: file_name = this.getMasterPath()
            if isempty(this.obs_full_name)
                this.updateObsFileName();
            end
            file_name = this.obs_full_name( this.obs_type == this.REC_MASTER);
        end

        function file_name = getReferencePath(this)
            % Get the file list of the reference receivers files
            % SYNTAX: file_name = this.getReferencePath()
            if isempty(this.obs_full_name)
                this.updateObsFileName();
            end
            file_name = this.obs_full_name(this.obs_type == this.REC_REFERENCE);
        end

        function [geometry, ev_point] = getGeometry(this)
            % Get the receiver coordinates in instrumental RF
            % SYNTAX: [geometry, ev_point] = this.getGeometry()

            geometry = this.xyz_ant;
            ev_point = this.xyz_ev_point;
        end

        function out = getImgDir(this)
            % Get the directory of UI images
            % SYNTAX: dir = this.getImgDir()
            out = this.img_dir;
        end

        function out = getNavArchive(this)
            % Get the list of archives to use for the search of navigational files
            out = this.preferred_archive;
        end

        function out = getNavGpsProvider(this)
            % Get the list of gps provider to use during the search for valid navigational files of GPS satellites
            out = this.preferred_gps;
        end

        function out = getNavGloProvider(this)
            % Get the list of gps provider to use during the search for valid navigational files of GLONASS satellites
            out = this.preferred_glo;
        end

        function out = getNavMixedProvider(this)
            % Get the list of gps provider to use during the search for valid navigational files of GPS satellites
            out = this.preferred_mxd;
        end

        function out = getNavClkType(this)
            % Get the order of preference of clock files to search for
            out = this.preferred_clk;
        end

        function out = getNavEphType(this)
            % Get the order of preference of orbits files to search for
            out = this.preferred_eph;
        end

        function out = getNavErpType(this)
            % Get the order of preference of erp files to search for
            out = this.preferred_erp;
        end
        
        function [addr, port, path, eph_name, clk_name, erp_name] = getCustomArchive(this)
            % Get the custom navigational provider parameters
            addr = this.custom_addr;
            port = this.custom_port;
            path = this.custom_path;
            eph_name = this.custom_name_eph;
            clk_name = this.custom_name_clk;
            erp_name = this.custom_name_erp;
        end

        function file_name = getFullNavEphPath(this, id)
            % Get the file list of ephemeris files
            % SYNTAX: file_name = this.getFullNavEphPath(id)
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
            % SYNTAX: file_name = this.getFullNavClkPath(id)
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
            % SYNTAX: file_name = this.getFullErpPath(id)
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
            % SYNTAX: nav_path = this.getNavPath()
            out = this.eph_name;
        end

        function clk_file = getNavClkFile(this)
            % Get the file name of the clock files
            % SYNTAX: clk_path = this.getClkPath()
            clk_file = this.clk_name;
        end
        
        function erp_file = getErpFile(this)
            % Get the file name of the ERP files
            % SYNTAX: erp_path = this.getErpPath()
            erp_file = this.erp_name;
        end

        function out = getNavEphDir(this)
            % Get the path to the navigational files
            % SYNTAX: nav_path = this.getNavEphDir()
            out = this.eph_dir;
        end

        function out = getNavClkDir(this)
            % Get the path to the clock files
            % SYNTAX: nav_path = this.getClkPath()
            out = this.clk_dir;
        end
        
        function out = getErpDir(this)
            % Get the path to the ERP files
            % SYNTAX: erp_path = this.getErpPath()
            out = this.erp_dir;
        end
        
        function out = getDcbDir(this)
            % Get the path to the DCB files
            % SYNTAX: dcb_path = this.getDcbPath()
            out = this.dcb_dir;
        end

        function out = getNavEphPath(this)
            % Get the path to the navigational files
            % SYNTAX: nav_path = this.getNavEphPath()
            out = File_Name_Processor.checkPath(strcat(this.eph_dir, filesep, this.eph_name));
        end

        function out = getNavClkPath(this)
            % Get the path to the clock files
            % SYNTAX: nav_path = this.getNavClkPath()
            out = File_Name_Processor.checkPath(strcat(this.clk_dir, filesep, this.clk_name));
        end
        
        function out = getErpPath(this)
            % Get the path to the ERP files
            % SYNTAX: erp_path = this.getErpPath()
            out = File_Name_Processor.checkPath(strcat(this.erp_dir, filesep, this.erp_name));
        end

        function out = getCrdFile(this)
            % Get the path of the stations coordinates file
            % SYNTAX: file_path = this.getCrdFile()
            if (isempty(this.crd_name))
                out = '';
            else
                out = this.checkCrdPath(strcat(this.crd_dir, filesep, this.crd_name));
            end
        end

        function out = getAtxFile(this)
            % Get the path of the antex file
            % SYNTAX: file_path = this.getAtxFile()
            if (isempty(this.atx_name))
                out = '';
            else
                out = this.checkAtxPath(strcat(this.atx_dir, filesep, this.atx_name));
            end
        end

        function out = getOceanFile(this)
            % Get the path of the ocean loading file
            % SYNTAX: file_path = this.getOceanFile()
            if (isempty(this.ocean_name))
                out = '';
            else
                out = this.checkOceanPath(strcat(this.ocean_dir, filesep, this.ocean_name));
            end
        end

        function out = getMetDir(this)
            % Get the path of the meteorological file dir
            % SYNTAX: file_path = this.getMetDir()
            if (isempty(this.met_dir))
                out = '';
            else
                out = this.checkMetPath(strcat(this.met_dir));
            end
        end

        function out = getMetFile(this, id)
            % Get the path of the meteorological file
            % SYNTAX: file_path = this.getMetFile()
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

        function out = getDtmPath(this)
            % Get the path of the ocean loading file
            % SYNTAX: file_path = this.getDtmPath()
            out = this.dtm_dir;
        end

        function out = getRefFile(this)
            % Get the path of the reference path file
            % SYNTAX: file_path = this.getRefPath()
            out = File_Name_Processor.checkPath(this.ref_graph_file);
        end

        function file_path = getGeoidFile(this)
            % Get the path of the geoid file
            % SYNTAX: file_path = this.getGeoidFile()
            file_path = File_Name_Processor.checkPath(strcat(this.geoid_dir, filesep, this.geoid_name));
        end


        function out_dir = getOutDir(this)
            % Get the path of the out folder
            % SYNTAX: out_dir = this.getOutDir()
            out_dir = File_Name_Processor.checkPath(this.out_dir);
        end

        function out_prefix = getOutPrefix(this)
            % Get the path of the out_prefix
            % SYNTAX: out_prefix = this.getOutPrefix()
            fnp = File_Name_Processor;
            out_prefix = fnp.checkPath(this.out_prefix);
        end

        function updateOutPath(this, date, session)
            % Update the full prefix of the putput files (replacing special keywords)
            % SYNTAX: this.updateOutPath(date, session);
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
            % SYNTAX: counter = getRunCounter(this)
            counter = this.run_counter;
        end
        
        function flag = isOutPosition(this)
            % Get the export status for the positions
            % SYNTAX: flag = this.isOutPosition()
            flag = this.flag_out_position;
        end
        function flag = isOutReportPDF(this)
            % Get the export status for the PDF report
            % SYNTAX: flag = this.isOutReportPDF()
            flag = this.flag_out_pdf_report;
        end

        function flag = isOutSettings(this)
            % Get the export status for the ini settings
            % SYNTAX: flag = this.isOutSettings()
            flag = this.flag_out_settings;
        end
        
        function flag = isOutCodeResPDF(this)
            % Get the export status for the code residuals plots
            % SYNTAX: flag = this.isOutCodeResPDF()
            flag = this.flag_out_pdf_code_res;
        end
        
        function flag = isOutPhaseResPDF(this)
            % Get the export status for the phase residuals plots
            % SYNTAX: flag = this.isOutPhaseResPDF()
            flag = this.flag_out_pdf_ph_res;
        end
        
        function flag = isOutDtPDF(this)
            % Get the export status for the Dt
            % SYNTAX: flag = this.isOutDtPDF()
            flag = this.flag_out_pdf_dt;
        end
        
        function flag = isOutTropoPDF(this)
            % Get the export status for the Troposphere results
            % SYNTAX: flag = this.isOutTropoPDF()
            flag = this.flag_out_pdf_tropo;
        end
        
        function flag = isOutBlockObj(this)
            % Get the export the obbject used for the block solution
            % SYNTAX: flag = this.isOutBlockObj()
            flag = this.flag_out_block_obj;
        end
        
        function flag = isOutKML(this)
            % Get the export status for the KML file
            % SYNTAX: flag = this.isOutKML()
            flag = this.flag_out_kml;
        end
        
        function flag = isOutNMEA(this)
            % Get the export status for the NMEA file
            % SYNTAX: flag = this.isOutNMEA()
            flag = this.flag_out_nmea;
        end
        
    end

    % =========================================================================
    %  SETTERS
    % =========================================================================
    methods
        function setNavPath(this, nav_dir)
            % Set the path to the navigational files
            % SYNTAX: this.getNavPath(nav_path)
            this.eph_dir = nav_dir;
        end

        function setNavEphFile(this, nav_name)
            % Set the file name of the navigational files
            % SYNTAX: this.setNavEphFile(nav_name)
            this.eph_name = nav_name;
        end

        function setNavClkPath(this, clk_dir)
            % Set the path to the clock files
            % SYNTAX: this.getClkPath(nav_path)
            this.clk_dir = clk_dir;
        end

        function setNavClkFile(this, clk_name)
            % Set the file name of the clock files
            % SYNTAX: this.getClkFile(nav_name)
            this.clk_name = clk_name;
        end
        
        function setErpPath(this, erp_dir)
            % Set the path to the clock files
            % SYNTAX: this.getClkPath(nav_path)
            this.erp_dir = erp_dir;
        end

        function setErpFile(this, erp_name)
            % Set the file name of the clock files
            % SYNTAX: this.getClkFile(erp_name)
            this.erp_name = erp_name;
        end

        function setOutPrefix(this, out_prefix)
            % Set the path of the out_prefix
            % SYNTAX: out_prefix = this.setOutPrefix(out_prefix)
            this.out_prefix = File_Name_Processor.checkPath(out_prefix);
        end

        function updateObsFileName(this)
            % Update the full name of the observations files (replacing special keywords)
            % SYNTAX: this.updateObsFileName();
            this.obs_full_name = {};
            fnp = File_Name_Processor();
            if ~iscell(this.obs_name)
                this.obs_name = {this.obs_name};
            end
            for i = 1 : numel(this.obs_name)
                this.obs_full_name{i} = fnp.dateKeyRepBatch(fnp.checkPath(strcat(this.obs_dir, filesep, this.obs_name{i})), this.sss_date_start,  this.sss_date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
            end
        end

        function updateMetFileName(this)
            % Update the full name of the observations files (replacing special keywords)
            % SYNTAX: this.updateObsFileName();
            this.met_full_name = {};
            fnp = File_Name_Processor();
            if ~iscell(this.met_name)
                this.met_name = {this.met_name};
            end
            this.met_full_name = {};
            for i = 1 : numel(this.met_name)
                this.met_full_name = [this.met_full_name; fnp.dateKeyRepBatch(fnp.checkPath(strcat(this.met_dir, filesep, this.met_name{i})), this.sss_date_start,  this.sss_date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop)];
            end
            %this.met_full_name = this.met_full_name(~isempty(this.met_full_name));
        end

        function updateNavFileName(this)
            % Update the full name of the navigational files (replacing special keywords)
            % SYNTAX: this.updateNavFileName();
            this.updateEphFileName();
            this.updateClkFileName();
        end

        function updateEphFileName(this)
            % Update the full name of the ephemerides files (replacing special keywords)
            % SYNTAX: this.updateEphFileName();
            this.eph_full_name = this.getEphFileName(this.sss_date_start, this.sss_date_stop);
        end

        function updateClkFileName(this)
            % Update the full name of the clock offset files (replacing special keywords)
            % SYNTAX: this.updateClkFileName();
            this.clk_full_name = this.getClkFileName(this.sss_date_start, this.sss_date_stop);
        end
        
        function updateErpFileName(this)
            % Update the full name of the ERP files (replacing special keywords)
            % SYNTAX: this.updateClkFileName();
            this.erp_full_name = this.getErpFileName(this.sss_date_start, this.sss_date_stop);
        end

        function date = getSessionStart(this)
            % SYNTAX: date = getSessionStart(this)
            date = this.sss_date_start;
        end

        function date = getSessionStop(this)
            % SYNTAX: date = getSessionStop(this)
            date = this.sss_date_stop;
        end

        function date = getSessionLimits(this)
            % SYNTAX: date = getSessionLimits(this)
            date = this.sss_date_start.getCopy;
            date.append(this.sss_date_stop);
        end

        function eph_full_name = getEphFileName(this, date_start, date_stop)
            % Get the full name of the ephemerides files (replacing special keywords)
            % SYNTAX: eph_full_name = getEphFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.eph_dir, filesep, this.eph_name));
            step_sec = fnp.getStepSec(file_name);

            if (~isempty(strfind(file_name, fnp.GPS_WD)) || ~isempty(strfind(file_name, fnp.GPS_WEEK)))
                date_start = date_start.getCopy; date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
                date_stop = date_stop.getCopy; date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin
            end
            eph_full_name = fnp.dateKeyRepBatch(file_name, date_start,  date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end

        function clk_full_name = getClkFileName(this, date_start, date_stop)
            % Get the full name of the clock offset files (replacing special keywords)
            % SYNTAX: clk_full_name = getClkFileName(this, date_start, date_stop)
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
            % SYNTAX: erp_full_name = getErpFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(strcat(this.erp_dir, filesep, this.erp_name));

            if (~isempty(strfind(file_name, fnp.GPS_WD)) || ~isempty(strfind(file_name, fnp.GPS_WEEK)))
                date_start = date_start.getCopy;
                date_stop = date_stop.getCopy;
            end
            erp_full_name = fnp.dateKeyRepBatch(file_name, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end

        function updateExternals(this)
            % Import the value of the external input files (stored in inputFile.ini)
            % SYNTAX: this.updateExternals();
            if ~isempty(this.input_file_ini_path)
                this.log.addWarning('Legacy importing input ini files');
                this.ext_ini = Ini_Manager(this.input_file_ini_path);
                if ~this.ext_ini.readFile()
                    this.log.addError(sprintf('Legacy import failed - "%s" can not be read!!!', this.ext_ini.getFileName()));
                    this.input_file_ini_path = '';
                else

                    fnp = File_Name_Processor();

                    % Import Rinex from old ini file
                    dir_receiver = this.ext_ini.getData('Receivers', 'data_path');
                    if ~isempty(dir_receiver)
                        this.obs_dir = dir_receiver;
                    end
                    dir_master = this.ext_ini.getData('Master', 'data_path');
                    if ~(isempty(dir_master)) && ~strcmp(dir_master, this.obs_dir)
                        this.log.addWarning('Importing legacy input file Master data_path seems different from Rover data_path - fix settings file manually');
                    end

                    % import Receivers/SEID names
                    name_receiver = this.ext_ini.getData('Receivers', 'file_name');
                    this.obs_name = {};
                    this.obs_full_name = {};
                    if ~isempty(name_receiver)
                        if iscell(name_receiver)
                            this.obs_name = name_receiver;
                            if this.isModeSEID()
                                this.obs_type = this.REC_REFERENCE * ones(1, numel(name_receiver));
                            else
                                this.obs_type = this.REC_TARGET  * ones(1, numel(name_receiver));
                            end
                        else
                            this.obs_name = {name_receiver};
                            if this.isModeSEID()
                                this.obs_type = this.REC_REFERENCE;
                            else
                                this.obs_type = this.REC_TARGET;
                            end
                        end
                    end
                    name_master = this.ext_ini.getData('Master', 'file_name');
                    if ~isempty(name_master)
                        if iscell(name_master)
                            this.obs_name = [this.obs_name name_master];
                            if this.isModeSEID()
                                this.obs_type = [this.obs_type this.REC_TARGET * ones(1, numel(name_master))];
                            else
                                this.obs_type = [this.obs_type this.REC_MASTER  * ones(1, numel(name_master))];
                            end
                        else
                            this.obs_name = [this.obs_name {name_master}];
                            if this.isModeSEID()
                                this.obs_type = [this.obs_type this.REC_TARGET];
                            else
                                this.obs_type = [this.obs_type this.REC_MASTER];
                            end
                        end
                    end


                    % get the output prefix
                    full_out_prefix = fnp.checkPath([this.out_dir filesep this.out_prefix]);
                    % make sure to have the name of the file and the name of the
                    % output folder
                    [~, out_prefix, ~] = fileparts(full_out_prefix); %#ok<PROP>
                    % list all the files that match the prefix
                    file_list = dir([full_out_prefix '*']);
                    % if there are no files in the putput folder
                    if isempty(file_list)
                        this.run_counter = this.RUN_COUNTER; % set the counter of the output == 0
                    else
                        % put the cell of the file in a single string
                        file_list = fnp.checkPath(strCell2Str({file_list(:).name},''));
                        % parse with regexp for output numbers -> get the maximum
                        this.run_counter = max(str2double(unique(regexp(file_list, [ '(?<=' out_prefix '_)[0-9]*(?=_)'], 'match')))) + 1; %#ok<PROP>
                        this.run_counter = iif(isempty(this.run_counter), this.RUN_COUNTER, this.run_counter);
                    end

                    % import receivers position
                    tmp_xyz_ant = zeros(3,this.getTargetCount());
                    ispresent = false;
                    for r = 1:this.getTargetCount()
                        tmp = this.ext_ini.getData('Antennas RF',['XYZ_ant' num2str(r)]);
                        if ~isempty(tmp)
                            ispresent = true;
                            tmp_xyz_ant(:, r) = tmp;
                        end
                    end
                    if ispresent
                        this.xyz_ant = tmp_xyz_ant;
                    end
                    tmp_ev_point = this.ext_ini.getData('Antennas RF','XYZ_ev_point');
                    if ~isempty(tmp_ev_point)
                        this.ev_point = tmp_ev_point;
                    end

                    % import file location and check for default folders
                    dir_path = this.ext_ini.getData('STATIONS_file', 'data_path');
                    file_path = this.ext_ini.getData('STATIONS_file', 'file_name');
                    file_name = File_Name_Processor.getFullPath(dir_path, file_path);
                    if ~isempty(file_name)
                        [file_dir, name, ext] = fileparts(file_name);
                        this.crd_dir = file_dir;
                        this.crd_name = strcat(name, ext);
                    end

                    % import file location and check for default folders
                    dir_path = this.ext_ini.getData('PCO_PCV_file', 'data_path');
                    file_path = this.ext_ini.getData('PCO_PCV_file', 'file_name');
                    file_name = File_Name_Processor.getFullPath(dir_path, file_path);
                    if ~isempty(file_name)
                        [file_dir, name, ext] = fileparts(file_name);
                        this.atx_dir = fnp.checkPath(strcat(file_dir, filesep));
                        this.atx_name = fnp.checkPath(strcat(name, ext));
                    end

                    % import file location and check for default folders
                    dir_path = this.ext_ini.getData('OCEAN_LOADING_file', 'data_path');
                    file_path = this.ext_ini.getData('OCEAN_LOADING_file', 'file_name');
                    file_name = File_Name_Processor.getFullPath(dir_path, file_path);
                    if ~isempty(file_name)
                        [file_dir, name, ext] = fileparts(file_name);
                        this.ocean_dir = file_dir;
                        this.ocean_name = strcat(name, ext);
                    end

                    % import file location and check for default folders
                    dir_path = this.ext_ini.getData('METEOROLOGICAL_file', 'data_path');
                    file_path = this.ext_ini.getData('METEOROLOGICAL_file', 'file_name');
                    file_name = File_Name_Processor.getFullPath(dir_path, file_path);
                    if ~isempty(file_name)
                        [file_dir, name, ext] = fileparts(this.checkMetPath(file_name));
                        this.met_dir = file_dir;
                        this.met_name = strcat(name, ext);
                        this.met_full_name = {};
                    end

                    % import DTM folders
                    dir_path = fnp.checkPath(this.ext_ini.getData('DTM','data_path'));
                    if ~isempty(dir_path)
                        this.dtm_dir = dir_path;
                    end

                    % import reference folder
                    dir_path = this.ext_ini.getData('RefPath', 'data_path');
                    file_path = this.ext_ini.getData('RefPath', 'file_name');
                    file_name = File_Name_Processor.getFullPath(dir_path, file_path);
                    if ~isempty(dir_path)
                        this.ref_graph_file = fnp.checkPath(file_name);
                    end

                    % import file location for navigational files
                    nav_path = this.ext_ini.getData('Navigational', 'data_path');
                    if ~isempty(nav_path)
                        this.eph_dir = nav_path;
                        this.clk_dir = nav_path;
                    end
                    nav_name = this.ext_ini.getData('Navigational', 'file_name');
                    if ~isempty(nav_name)
                        this.eph_name = nav_name;
                        this.clk_name = nav_name;
                    end
                    this.input_file_ini_path = '';
                    this.eph_full_name = '';
                    this.clk_full_name = '';
                    this.check();
                    this.updateObsFileName();
                end
            end
        end

        function setProcessingTime(this, first_epoch, last_epoch, update_iif_smaller)
            % Set the first/last epoch of processing
            % SYNTAX: dir = this.setProcessingTime(first_epoch, last_epoch, <update_iif_smaller == false>)
            if nargin == 3
                update_iif_smaller = false;
            end

            if (~update_iif_smaller) || (this.sss_date_start.isempty()) || (this.sss_date_start.getMatlabTime() < first_epoch.getMatlabTime())
                this.sss_date_start = first_epoch.getCopy();
            end
            if (~update_iif_smaller) || (this.sss_date_stop.isempty()) || (this.sss_date_stop.getMatlabTime() < last_epoch.getMatlabTime())
                this.sss_date_stop = last_epoch.getCopy();
            end
        end

        function setPrjHome(this, prj_home)
            % Set home folder of the project
            % SYNTAX: dir = this.setPrjHome(prj_home)
            this.prj_home = prj_home;
        end

        function setDeprecateIniPath(this, new_path)
            % Set the legacy deprecate ini of the input files
            % SYNTAX: dir = this.setDeprecateIniPath(new_path)
            this.input_file_ini_path = new_path;
        end

        function setFilePath(this, file_path)
            % Set the file name of the current settings
            % SYNTAX: this.setFilePath(file_path)
            [path_str, name, ~] = fileparts(file_path);
            this.cur_ini = [path_str filesep name '.ini'];
        end

        function updatePrj(this, file_path)
            % Set the file project name / home / file_path from file_path
            % SYNTAX: this.autoUpdatePrj(this, file_path)
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
                    this.log.addMessage(sprintf(' - DTM tile header in %s have been read', [dtm_dir filesep 'tiles' filesep 'tile_header']));
                    load([dtm_dir filesep 'tiles' filesep 'tile_georef'], 'tile_georef');
                    this.tile_georef = tile_georef;
                    this.log.addMessage(sprintf(' - DTM tile georef in %s have been read', [dtm_dir filesep 'tiles' filesep 'tile_georef']));
                catch
                    this.log.addWarning(sprintf('Failed to read DTM stored in %s', [dtm_dir '/tiles/']));
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

        function checkCellStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.Field(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            switch nargin
                case 2, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid);
                case 4, this.(field_name) = this.checkCellString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid, check_existence);
                otherwise, error('Settings checkCellStringField called with the wrong number of parameters');
            end
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

        function checkPathField(this, field_name, empty_is_valid, check_existence)
            % Check if a string path field of the object is a valid path
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.checkPathField(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            fnp = File_Name_Processor();
            this.(field_name) = fnp.getFullDirPath(this.(field_name), this.prj_home, [], fnp.getFullDirPath(this.(upper(field_name))));
            switch nargin
                case 2, this.(field_name) = this.checkString(field_name, this.(field_name), fnp.getFullDirPath(this.(upper(field_name)), this.prj_home));
                case 3, this.(field_name) = this.checkString(field_name, this.(field_name), fnp.getFullDirPath(this.(upper(field_name)), this.prj_home), empty_is_valid);
                case 4, this.(field_name) = this.checkString(field_name, this.(field_name), fnp.getFullDirPath(this.(upper(field_name)), this.prj_home), empty_is_valid, check_existence);
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

        % File type specific ----------------------------------------------

        function file_path = checkCrdPath(this, file_path)
            % Check if the crd file exists, if not try to look for it into the default dirs
            % SYNTAX: file_path = this.checkCrdPath(<file_path>)
            fnp = File_Name_Processor();
            file_path = fnp.checkPath(file_path);
            if ~isempty(file_path) && ~exist(file_path, 'file')
                [~, name, ext] = fileparts(file_path);
                % check for existence in the local project folder standard location
                tmp_path = fnp.checkPath([this.prj_home this.CRD_DIR(length(IO_Settings.DEFAULT_DIR_IN)+1:end) filesep name ext]);
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
            % SYNTAX: file_path = this.checkAtxPath(<file_path>)
            fnp = File_Name_Processor();
            file_path = fnp.checkPath(file_path);

            if ~isempty(file_path) && ~exist(file_path, 'file')
                [~, name, ext] = fileparts(file_path);
                % check for existence in the local project folder standard location
                tmp_path = fnp.checkPath([this.prj_home this.ATX_DIR(length(IO_Settings.DEFAULT_DIR_IN)+1:end) filesep name ext]);
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
            % SYNTAX: file_path = this.checkAtxPath(<file_path>)
            fnp = File_Name_Processor();
            file_path = fnp.checkPath(file_path);

            if ~isempty(file_path) && ~exist(file_path, 'file')
                fnp = File_Name_Processor();
                [~, name, ext] = fileparts(file_path);
                % check for existence in the local project folder standard location
                tmp_path = fnp.checkPath([this.prj_home this.OCEAN_DIR(length(IO_Settings.DEFAULT_DIR_IN)+1:end) filesep name ext]);
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
            % SYNTAX: file_path = this.checkMetPath(<file_path>)
            fnp = File_Name_Processor();
            file_path = fnp.checkPath(file_path);

            if ~isempty(file_path) && ~exist(file_path, 'file')
                fnp = File_Name_Processor();
                [~, name, ext] = fileparts(file_path);
                % check for existence in the local project folder standard location
                tmp_path = File_Name_Processor.checkPath([this.prj_home this.MET_DIR(length(IO_Settings.DEFAULT_DIR_IN)+1:end) filesep name ext]);
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
    %  TEST PARAMETERS VALIDITY
    % =========================================================================
    methods (Access = 'public')
        function check(this)
            % Check the validity of the fields
            % SYNTAX: this.check();

            this.checkStringField('prj_name', false);
            this.checkStringField('prj_home', false, true);
            this.checkStringField('cur_ini', false);

            this.checkCellStringField('preferred_archive', false);
            this.checkCellStringField('preferred_gps', false);
            this.checkCellStringField('preferred_glo', false);
            this.checkCellStringField('preferred_mxd', false);
            this.checkCellStringField('preferred_eph', false);
            this.checkCellStringField('preferred_erp', false);
            this.checkCellStringField('preferred_clk', false);

            this.checkStringField('custom_addr', false);
            this.checkStringField('custom_port', false);
            this.checkStringField('custom_path', false);
            this.checkStringField('custom_name_eph', false);
            this.checkStringField('custom_name_clk', false);
            this.checkStringField('custom_name_erp', false);

            this.checkStringField('sss_id_list', false);
            this.checkStringField('sss_id_start', false);
            this.checkStringField('sss_id_stop', false);

            this.checkPathField('obs_dir', false, true);
            this.checkCellStringField('obs_name', false);

            if numel(this.obs_name) > numel(this.obs_type)
                this.log.addError('I read more obs file names than obs types -> fixing setting all the other types to zeros\nplease review the applied solution!!!');
                tmp = this.obs_type;
                this.obs_type = zeros(numel(this.obs_name),1);
                this.obs_type(1:numel(tmp)) = tmp(:);
            end
            if numel(this.obs_name) < numel(this.obs_type)
                this.log.addError('I read less obs file names than obs types -> fixing cutting type array\nplease review the applied solution!!!');
                this.obs_type = this.obs_type(1:numel(this.obs_name));
            end

            this.checkPathField('atx_dir', false);
            this.checkStringField('atx_name', false);

            this.checkStringField('input_file_ini_path', true);

            this.checkPathField('eph_dir', false, true);
            % When the ephemeris file inserted here is not found -> the automatic downloader will dowload the proper file
            this.checkStringField('eph_name', true);
            this.checkPathField('clk_dir', false, true);
            this.checkStringField('clk_name', true);
            this.checkPathField('erp_dir', false, true);
            this.checkStringField('erp_name', true);
            this.checkPathField('crx_dir', false);
            this.checkPathField('dcb_dir', false);
            this.checkPathField('ems_dir', true);

            this.checkPathField('crd_dir', false);
            this.checkPathField('met_dir', false);
            this.checkStringField('ocean_dir', false);

            this.checkStringField('ref_graph_file', true);
            this.checkPathField('geoid_dir', false);
            this.checkStringField('geoid_name', false);
            this.checkPathField('dtm_dir', true);

            this.checkPathField('img_dir', false, true);

            this.checkStringField('out_prefix', true);

            if (this.run_counter_is_set) || ~(isempty(this.run_counter))
                this.checkNumericField('run_counter',[0 1e6]);
            end
            
            this.checkLogicalField('flag_out_position');
            this.checkLogicalField('flag_out_pdf_report');
            this.checkLogicalField('flag_out_settings');
            this.checkLogicalField('flag_out_pdf_code_res');
            this.checkLogicalField('flag_out_pdf_ph_res');
            this.checkLogicalField('flag_out_pdf_dt');
            this.checkLogicalField('flag_out_pdf_tropo');
            this.checkLogicalField('flag_out_block_obj');
            this.checkLogicalField('flag_out_kml');
            this.checkLogicalField('flag_out_nmea');


            % Check size of xyz antenna
            if (size(this.xyz_ant,2) ~= this.getTargetCount())
                tmp_xyz_ant = zeros(3,this.getTargetCount());
                % check num of coordinates
                if (size(this.xyz_ant, 1) == 3)
                    tmp_xyz_ant(:, 1 : min(this.getTargetCount(), size(this.xyz_ant,2))) = this.xyz_ant(:, 1 : min(this.getTargetCount(), size(this.xyz_ant,2)));
                    this.xyz_ant = tmp_xyz_ant;
                else
                    this.xyz_ant = tmp_xyz_ant;
                end
            end
            if (size(this.xyz_ev_point,1) ~= 3)
                this.xyz_ev_point = this.XYZ_EV_POINT;
            end
        end

        function status = checkTargetFiles(this, go_verbose)
            % check the availability of all the rinex files
            % SYNTAX: status = this.checkTargetFiles(go_verbose)
            % status is an array containing the file status for each receiver
            %   0: all file are present
            %  -1: no file are present
            %   1: at least one file is present but not all
            if nargin == 1
                go_verbose = false;
            end
            status = checkReceiverFiles(this, this.REC_TARGET, go_verbose);
        end

        function status = checkMasterFiles(this, go_verbose)
            % check the availability of all the rinex files
            % SYNTAX: status = this.checkMasterFiles(go_verbose)
            % status is an array containing the file status for each receiver
            %   0: all file are present
            %  -1: no file are present
            %   1: at least one file is present but not all
            if nargin == 1
                go_verbose = false;
            end
            status = checkReceiverFiles(this, this.REC_MASTER, go_verbose);
        end

        function status = checkReferenceFiles(this, go_verbose)
            % check the availability of all the rinex files
            % SYNTAX: status = this.checkReferenceFiles(go_verbose)
            % status is an array containing the file status for each receiver
            %   0: all file are present
            %  -1: no file are present
            %   1: at least one file is present but not all
            if nargin == 1
                go_verbose = false;
            end
            status = checkReceiverFiles(this, this.REC_REFERENCE, go_verbose);
        end

        function [status_trg, status_ref, status_mst] = checkAllReceiversFiles(this, go_verbose)
            % check the availability of all the rinex files
            % SYNTAX: status = this.checkAllReceiversFiles(go_verbose)
            % status is an array containing the file status for each receiver
            %   0: all file are present
            %  -1: no file are present
            %   1: at least one file is present but not all
            if nargin == 1
                go_verbose = false;
            end
            status_trg = checkReceiverFiles(this, this.REC_TARGET, go_verbose);
            status_ref = checkReceiverFiles(this, this.REC_REFERENCE, go_verbose);
            status_mst = checkReceiverFiles(this, this.REC_MASTER, go_verbose);
        end

        function status = checkReceiverFiles(this, obs_type, go_verbose)
            % check the availability of all the rinex files
            % SYNTAX: status = this.checkReceiverFiles(obs_type, go_verbose)
            % status is an array containing the file status for each receiver
            %   0: all file are present
            %  -1: no file are present
            %   1: at least one file is present but not all

            if nargin == 2
                go_verbose = false;
            end

            switch obs_type
                case this.REC_TARGET
                    obs_type_name = 'target';
                    n_rec = this.getTargetCount();
                    file_name_all = this.getTargetPath();
                case this.REC_MASTER
                    obs_type_name = 'master';
                    n_rec = this.getMasterCount();
                    file_name_all = this.getMasterPath();
                case this.REC_REFERENCE
                    obs_type_name = 'reference';
                    n_rec = this.getReferenceCount();
                    file_name_all = this.getReferencePath();
            end

            fnp = File_Name_Processor();

            % If no receiver have been found
            if n_rec == 0
                status = -1;
            else
                status = 0;
                for r = 1 : n_rec

                    if go_verbose
                        this.log.addMessage(sprintf('Checking files for %s receiver %d', obs_type_name, r));
                    end
                    file_name = file_name_all{r};
                    file_count = 0;
                    for f = 1 : numel(file_name)
                        full_path = fnp.checkPath(file_name{f});
                        file_ok = exist(full_path, 'file');
                        file_count = file_count + uint16(logical(file_ok));
                        if go_verbose
                            if file_ok
                                this.log.addStatusOk(sprintf('%s is present', full_path));
                                status = 1;
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
            eph_ok = true;

            state = Go_State.getCurrentSettings();
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

            clk_ok = true;
            state = Go_State.getCurrentSettings();
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

            erp_ok = true;
            state = Go_State.getCurrentSettings();
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
    %  TEST
    % =========================================================================
    methods (Static, Access = 'public')
        function test()
            % test the class
            s = IO_Settings();
            s.testInterfaceRoutines();
        end
    end
end
